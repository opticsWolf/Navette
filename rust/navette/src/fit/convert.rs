// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::convert — instrument data into the solver's conventions
//! (decision 3, plan §3.5.5).
//!
//! The engine works in the solver's conventions only: wavelengths in nm,
//! angles and Psi/Delta in radians, Delta in `(-pi, pi]` with the solver's
//! sign, reflectances and transmittances as fractions. Data an instrument
//! wrote enters through [`to_solver`], which maps its
//! [`InstrumentConvention`] onto those, propagates σ through the map, and
//! records on each [`SpectrumData`] how it was read (`source`, carried
//! into the fit report). [`from_solver`] is the inverse.
//!
//! Presets are data (`conventions.json`, [`InstrumentConvention::preset`]);
//! the explicit form takes any combination. No file formats here: reading
//! a vendor file is the Python door's job, or a later reader module.
//!
//! **`tan Ψ, cos Δ` loses Δ's sign**: it comes back in `[0, π]`. Fit it
//! only where the model's Δ has that sign, or measure N, C, S.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::smatrix::angle::Radians;

use super::data::{SpectrumData, wrap_angle};
use super::model::DEFAULT_ENV;

/// `h·c` in nm·eV.
const HC_NM_EV: f64 = 1239.841984;

const PRESETS: &str = include_str!("conventions.json");

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum WavelengthUnit {
    Nm,
    Um,
    Angstrom,
    /// Photon energy.
    Ev,
    /// cm⁻¹.
    Wavenumber,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AngleUnit {
    Deg,
    Rad,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ValueScale {
    Fraction,
    Percent,
}

/// The range Δ is written in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DeltaRange {
    /// `(-180°, 180°]`.
    Symmetric,
    /// `[0°, 360°)`.
    Positive,
}

/// Δ's sign against the solver's (the time convention of the instrument's
/// software, `exp(+iωt)` vs `exp(-iωt)`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DeltaSign {
    Solver,
    Opposite,
}

/// What an ellipsometric record's columns hold.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EllForm {
    /// Ψ, Δ.
    PsiDelta,
    /// tan Ψ, cos Δ.
    TanPsiCosDelta,
    /// N = cos 2Ψ, C = sin 2Ψ cos Δ, S = sin 2Ψ sin Δ.
    Ncs,
}

/// Where σ comes from (in the instrument's units; converted with the
/// values).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SigmaSource {
    /// The raw record's `sigma` columns.
    Column,
    Constant(f64),
    /// `fraction · |value|` (with a floor at 1e-12 of the unit).
    Relative(f64),
    /// σ = 1 in solver units: unit weight.
    None,
}

/// How an instrument wrote its data.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InstrumentConvention {
    pub wavelength_unit: WavelengthUnit,
    /// The incidence angle's unit, and Ψ's and Δ's.
    pub angle_unit: AngleUnit,
    /// R, T, A.
    pub value_scale: ValueScale,
    pub delta_range: DeltaRange,
    pub delta_sign: DeltaSign,
    pub ell_form: EllForm,
    #[serde(default = "sigma_none")]
    pub sigma: SigmaSource,
}

fn default_env() -> String {
    DEFAULT_ENV.to_string()
}

fn sigma_none() -> SigmaSource {
    SigmaSource::None
}

impl InstrumentConvention {
    /// A named preset from `conventions.json`.
    pub fn preset(name: &str) -> Result<Self, String> {
        let all: serde_json::Map<String, serde_json::Value> =
            serde_json::from_str(PRESETS).expect("conventions.json parses");
        let v = all.get(name).ok_or_else(|| {
            format!(
                "no convention preset {name:?}; presets: {:?}",
                Self::preset_names()
            )
        })?;
        serde_json::from_value(v.clone()).map_err(|e| format!("preset {name:?}: {e}"))
    }

    /// Every preset name, in file order.
    pub fn preset_names() -> Vec<String> {
        let all: serde_json::Map<String, serde_json::Value> =
            serde_json::from_str(PRESETS).expect("conventions.json parses");
        all.keys().cloned().collect()
    }

    pub fn with_sigma(mut self, sigma: SigmaSource) -> Self {
        self.sigma = sigma;
        self
    }

    /// One line, for reports.
    pub fn describe(&self) -> String {
        serde_json::to_string(self).expect("plain data")
    }

    fn angle_to_rad(&self) -> f64 {
        match self.angle_unit {
            AngleUnit::Deg => PI / 180.0,
            AngleUnit::Rad => 1.0,
        }
    }

    fn sign(&self) -> f64 {
        match self.delta_sign {
            DeltaSign::Solver => 1.0,
            DeltaSign::Opposite => -1.0,
        }
    }

    fn nm_of(&self, x: f64) -> f64 {
        match self.wavelength_unit {
            WavelengthUnit::Nm => x,
            WavelengthUnit::Um => 1000.0 * x,
            WavelengthUnit::Angstrom => 0.1 * x,
            WavelengthUnit::Ev => HC_NM_EV / x,
            WavelengthUnit::Wavenumber => 1e7 / x,
        }
    }

    fn unit_of(&self, nm: f64) -> f64 {
        match self.wavelength_unit {
            WavelengthUnit::Nm => nm,
            WavelengthUnit::Um => nm / 1000.0,
            WavelengthUnit::Angstrom => 10.0 * nm,
            WavelengthUnit::Ev => HC_NM_EV / nm,
            WavelengthUnit::Wavenumber => 1e7 / nm,
        }
    }
}

/// Which side an ellipsometric record measures.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Side {
    R,
    T,
}

/// What a raw record holds.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RawQuantity {
    /// One column of a solver channel (`Rs`, `T_avg`, ...).
    Spectrum(String),
    /// Two (Ψ, Δ or tan Ψ, cos Δ) or three (N, C, S) columns.
    Ellipsometry(Side),
}

/// One record as an instrument wrote it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RawDataset {
    pub id: String,
    /// The sample; [`DEFAULT_ENV`] when a record does not name one.
    #[serde(default = "default_env")]
    pub env: String,
    pub quantity: RawQuantity,
    /// In the convention's angle unit.
    pub angle: f64,
    /// In the convention's wavelength unit.
    pub abscissa: Vec<f64>,
    pub columns: Vec<Vec<f64>>,
    /// Same shape as `columns` for [`SigmaSource::Column`], else empty.
    #[serde(default)]
    pub sigma: Vec<Vec<f64>>,
}

impl RawDataset {
    pub fn new(
        id: impl Into<String>,
        quantity: RawQuantity,
        angle: f64,
        abscissa: Vec<f64>,
        columns: Vec<Vec<f64>>,
    ) -> Self {
        Self {
            id: id.into(),
            env: DEFAULT_ENV.to_string(),
            quantity,
            angle,
            abscissa,
            columns,
            sigma: Vec::new(),
        }
    }
}

/// Convert `raw` into solver-convention datasets: one for a spectrum,
/// `<id>.psi` and `<id>.delta` for an ellipsometric record.
pub fn to_solver(
    raw: &RawDataset,
    conv: &InstrumentConvention,
) -> Result<Vec<SpectrumData>, String> {
    let bad = |why: String| format!("raw dataset {:?}: {why}", raw.id);
    let n = raw.abscissa.len();
    let want = match (&raw.quantity, conv.ell_form) {
        (RawQuantity::Spectrum(_), _) => 1,
        (RawQuantity::Ellipsometry(_), EllForm::Ncs) => 3,
        (RawQuantity::Ellipsometry(_), _) => 2,
    };
    if raw.columns.len() != want {
        return Err(bad(format!(
            "{} columns; {:?} in this convention needs {want}",
            raw.columns.len(),
            raw.quantity
        )));
    }
    if let Some(c) = raw.columns.iter().find(|c| c.len() != n) {
        return Err(bad(format!(
            "a column of {} values for {n} abscissae",
            c.len()
        )));
    }
    if let Some(x) = raw.abscissa.iter().find(|x| !(**x > 0.0 && x.is_finite())) {
        return Err(bad(format!("abscissa {x} is not a positive finite number")));
    }
    // σ in the instrument's units, per column and point.
    let sigma: Vec<Vec<f64>> = match conv.sigma {
        SigmaSource::Column => {
            if raw.sigma.len() != want || raw.sigma.iter().any(|s| s.len() != n) {
                return Err(bad(
                    "the convention reads sigma columns; give one per value column".into(),
                ));
            }
            raw.sigma.clone()
        }
        SigmaSource::Constant(s) => vec![vec![s; n]; want],
        SigmaSource::Relative(f) => raw
            .columns
            .iter()
            .map(|c| c.iter().map(|v| (f * v.abs()).max(1e-12)).collect())
            .collect(),
        SigmaSource::None => Vec::new(),
    };
    let wl: Vec<f64> = raw.abscissa.iter().map(|&x| conv.nm_of(x)).collect();
    let a = conv.angle_to_rad();
    let angle = Radians::new(raw.angle * a).map_err(bad)?;
    let source = Some(conv.describe());
    let make = |id: String, channel: String, values: Vec<f64>, sig: Option<Vec<f64>>| {
        let mut d = SpectrumData::new(id, channel, angle, wl.clone(), values, 1.0);
        d.env = raw.env.clone();
        if let Some(s) = sig {
            d.sigma = s;
        }
        d.source = source.clone();
        d
    };

    match &raw.quantity {
        RawQuantity::Spectrum(channel) => {
            let k = match conv.value_scale {
                ValueScale::Fraction => 1.0,
                ValueScale::Percent => 0.01,
            };
            let values = raw.columns[0].iter().map(|v| v * k).collect();
            let sig = sigma.first().map(|s| s.iter().map(|v| v * k).collect());
            Ok(vec![make(raw.id.clone(), channel.clone(), values, sig)])
        }
        RawQuantity::Ellipsometry(side) => {
            let s = conv.sign();
            let has_sigma = !sigma.is_empty();
            let mut psi = Vec::with_capacity(n);
            let mut delta = Vec::with_capacity(n);
            let mut s_psi = Vec::with_capacity(n);
            let mut s_delta = Vec::with_capacity(n);
            for i in 0..n {
                let col = |c: usize| raw.columns[c][i];
                let sg = |c: usize| if has_sigma { sigma[c][i] } else { 0.0 };
                let (p, d, sp, sd) = match conv.ell_form {
                    EllForm::PsiDelta => {
                        (col(0) * a, wrap_angle(s * col(1) * a), sg(0) * a, sg(1) * a)
                    }
                    EllForm::TanPsiCosDelta => {
                        let (t, c) = (col(0), col(1));
                        if !(t >= 0.0) || c.abs() > 1.0 + 1e-12 {
                            return Err(bad(format!(
                                "point {i}: tan Psi {t} / cos Delta {c} out of range"
                            )));
                        }
                        let c = c.clamp(-1.0, 1.0);
                        (
                            t.atan(),
                            c.acos(),
                            sg(0) / (1.0 + t * t),
                            sg(1) / (1.0 - c * c).sqrt().max(1e-6),
                        )
                    }
                    EllForm::Ncs => {
                        let (nn, cc, ss) = (col(0), col(1), col(2));
                        let rho2 = cc * cc + ss * ss;
                        let rho = rho2.sqrt();
                        let r2 = nn * nn + rho2;
                        if !(r2 > 0.0) {
                            return Err(bad(format!("point {i}: N = C = S = 0")));
                        }
                        // Ψ = ½·atan2(ρ, N), Δ = atan2(S, C); σ through the
                        // Jacobian of that map, columns independent.
                        let (sn, sc, sss) = (sg(0), sg(1), sg(2));
                        let dpsi_dn = -0.5 * rho / r2;
                        let dpsi_drho = 0.5 * nn / r2;
                        let (drho_dc, drho_ds) = if rho > 0.0 {
                            (cc / rho, ss / rho)
                        } else {
                            (0.0, 0.0)
                        };
                        let sp = ((dpsi_dn * sn).powi(2)
                            + (dpsi_drho * drho_dc * sc).powi(2)
                            + (dpsi_drho * drho_ds * sss).powi(2))
                        .sqrt();
                        let sd = if rho2 > 0.0 {
                            ((ss / rho2 * sc).powi(2) + (cc / rho2 * sss).powi(2)).sqrt()
                        } else {
                            f64::INFINITY
                        };
                        (0.5 * rho.atan2(nn), wrap_angle(s * ss.atan2(cc)), sp, sd)
                    }
                };
                psi.push(p);
                delta.push(d);
                s_psi.push(sp);
                s_delta.push(sd);
            }
            let suffix = match side {
                Side::R => "R",
                Side::T => "T",
            };
            let (sig_p, sig_d) = if has_sigma {
                if let Some(s) = s_psi
                    .iter()
                    .chain(&s_delta)
                    .find(|s| !(**s > 0.0 && s.is_finite()))
                {
                    return Err(bad(format!(
                        "a propagated sigma is {s}; the record is degenerate there"
                    )));
                }
                (Some(s_psi), Some(s_delta))
            } else {
                (None, None)
            };
            Ok(vec![
                make(
                    format!("{}.psi", raw.id),
                    format!("Psi_{suffix}"),
                    psi,
                    sig_p,
                ),
                make(
                    format!("{}.delta", raw.id),
                    format!("Delta_{suffix}"),
                    delta,
                    sig_d,
                ),
            ])
        }
    }
}

/// The inverse of [`to_solver`]: solver-convention values written as
/// `conv` would. `data` is one spectrum, or a `Psi_*` and a `Delta_*` on
/// the same wavelengths and angle. σ is not written back.
pub fn from_solver(
    id: &str,
    data: &[&SpectrumData],
    conv: &InstrumentConvention,
) -> Result<RawDataset, String> {
    let first = data.first().ok_or("from_solver: no data")?;
    let a = conv.angle_to_rad();
    let abscissa: Vec<f64> = first.wavelengths.iter().map(|&w| conv.unit_of(w)).collect();
    let angle = first.angle.get() / a;
    let (quantity, columns) = match data {
        [d] if !d.channel.starts_with("Psi_") && !d.channel.starts_with("Delta_") => {
            let k = match conv.value_scale {
                ValueScale::Fraction => 1.0,
                ValueScale::Percent => 100.0,
            };
            (
                RawQuantity::Spectrum(d.channel.clone()),
                vec![d.values.iter().map(|v| v * k).collect()],
            )
        }
        [p, d] if p.channel.starts_with("Psi_") && d.channel.starts_with("Delta_") => {
            let side = match (&p.channel[4..], &d.channel[6..]) {
                ("R", "R") => Side::R,
                ("T", "T") => Side::T,
                _ => return Err("from_solver: Psi and Delta of different sides".into()),
            };
            if p.wavelengths != d.wavelengths || p.angle != d.angle {
                return Err("from_solver: Psi and Delta on different points".into());
            }
            let s = conv.sign();
            let inst_delta: Vec<f64> = d
                .values
                .iter()
                .map(|&x| {
                    let w = wrap_angle(s * x);
                    match conv.delta_range {
                        DeltaRange::Symmetric => w,
                        DeltaRange::Positive if w < 0.0 => w + 2.0 * PI,
                        DeltaRange::Positive => w,
                    }
                })
                .collect();
            let columns = match conv.ell_form {
                EllForm::PsiDelta => vec![
                    p.values.iter().map(|v| v / a).collect(),
                    inst_delta.iter().map(|v| v / a).collect(),
                ],
                EllForm::TanPsiCosDelta => vec![
                    p.values.iter().map(|v| v.tan()).collect(),
                    inst_delta.iter().map(|v| v.cos()).collect(),
                ],
                EllForm::Ncs => vec![
                    p.values.iter().map(|v| (2.0 * v).cos()).collect(),
                    p.values
                        .iter()
                        .zip(&inst_delta)
                        .map(|(v, d)| (2.0 * v).sin() * d.cos())
                        .collect(),
                    p.values
                        .iter()
                        .zip(&inst_delta)
                        .map(|(v, d)| (2.0 * v).sin() * d.sin())
                        .collect(),
                ],
            };
            (RawQuantity::Ellipsometry(side), columns)
        }
        _ => {
            return Err(
                "from_solver: give one spectrum, or a Psi_* and a Delta_* dataset in that order"
                    .into(),
            );
        }
    };
    Ok(RawDataset {
        id: id.to_string(),
        env: first.env.clone(),
        quantity,
        angle,
        abscissa,
        columns,
        sigma: Vec::new(),
    })
}
