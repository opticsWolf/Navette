// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::data — measured spectra and the residual rows they produce.
//!
//! Values are in the solver's conventions (decision 3): reflectances and
//! transmittances as fractions, Psi, Delta and phases in radians with the
//! solver's sign. Anything an instrument wrote is converted before it gets
//! here (`fit::convert`, phase 3), never inside a residual.
//!
//! Every row is `sqrt(weight) * (sim - value) / sigma`. Angular channels
//! (Delta, phases, differential phase) wrap `sim - value` into `(-pi, pi]`
//! first, so a Delta measured as `-179 deg` against a model at `+179 deg`
//! costs 2 deg, not 358.
//!
//! [`NkData`] carries measured optical constants of one material (n/k,
//! ε₁/ε₂, or one of them) for dispersion fits; its rows follow the same
//! convention, column after column.

use std::f64::consts::PI;

use num_complex::Complex64;

use crate::smatrix::angle::Radians;
use crate::smatrix::core_engine::*;

/// Channel name (a `Solution::f64maps` key) and the request bit it needs.
const CHANNELS: [(&str, u64); 37] = [
    ("Rs", REQ_RS),
    ("Rp", REQ_RP),
    ("Ts", REQ_TS),
    ("Tp", REQ_TP),
    ("R_avg", REQ_R_AVG),
    ("T_avg", REQ_T_AVG),
    ("A_s", REQ_A_S),
    ("A_p", REQ_A_P),
    ("A_avg", REQ_A_AVG),
    ("Psi_R", REQ_PSI_R),
    ("Psi_T", REQ_PSI_T),
    ("Delta_R", REQ_DELTA_R),
    ("Delta_T", REQ_DELTA_T),
    ("DOP_R", REQ_DOP_R),
    ("DOP_T", REQ_DOP_T),
    ("Diattenuation_R", REQ_DIATT_R),
    ("Diattenuation_T", REQ_DIATT_T),
    ("S0_R", REQ_S0_R),
    ("S1_R", REQ_S1_R),
    ("S2_R", REQ_S2_R),
    ("S3_R", REQ_S3_R),
    ("S0_T", REQ_S0_T),
    ("S1_T", REQ_S1_T),
    ("S2_T", REQ_S2_T),
    ("S3_T", REQ_S3_T),
    ("Retardance_R", REQ_RETARD_R),
    ("Retardance_T", REQ_RETARD_T),
    ("phi_rs", REQ_PHI_RS),
    ("phi_rp", REQ_PHI_RP),
    ("phi_ts", REQ_PHI_TS),
    ("phi_tp", REQ_PHI_TP),
    ("phi_rbs", REQ_PHI_RBS),
    ("phi_rbp", REQ_PHI_RBP),
    ("phi_tbs", REQ_PHI_TBS),
    ("phi_tbp", REQ_PHI_TBP),
    ("PDts", REQ_PD_TS),
    ("PDtp", REQ_PD_TP),
];

/// Request bit for a channel name; unknown names list the known ones.
pub fn channel_request(channel: &str) -> Result<u64, String> {
    CHANNELS
        .iter()
        .find(|(name, _)| *name == channel)
        .map(|(_, bit)| *bit)
        .ok_or_else(|| {
            format!(
                "unknown channel {channel:?}; known: {}",
                CHANNELS.map(|(n, _)| n).join(", ")
            )
        })
}

/// Whether a channel's residual wraps modulo 2 pi.
pub(crate) fn is_angular(channel: &str) -> bool {
    channel.starts_with("Delta_") || channel.starts_with("phi_") || channel.starts_with("PD")
}

/// `d` wrapped into `(-pi, pi]`.
pub(crate) fn wrap_angle(d: f64) -> f64 {
    let w = d - 2.0 * PI * (d / (2.0 * PI)).round();
    if w <= -PI { w + 2.0 * PI } else { w }
}

/// One measured curve: a channel at one angle over wavelengths.
#[derive(Debug, Clone, PartialEq)]
pub struct SpectrumData {
    pub id: String,
    /// A `Solution` map key: `Rs`, `Tp`, `R_avg`, `Psi_R`, `Delta_R`, `S1_R`, ...
    pub channel: String,
    pub angle: Radians,
    /// nm.
    pub wavelengths: Vec<f64>,
    pub values: Vec<f64>,
    /// One value for every point, or one per point.
    pub sigma: Vec<f64>,
    pub weight: f64,
    /// Keep only points with `lo <= wavelength <= hi` (nm).
    pub window: Option<(f64, f64)>,
}

impl SpectrumData {
    /// A curve with one sigma for every point, weight 1, no window.
    pub fn new(
        id: impl Into<String>,
        channel: impl Into<String>,
        angle: Radians,
        wavelengths: Vec<f64>,
        values: Vec<f64>,
        sigma: f64,
    ) -> Self {
        Self {
            id: id.into(),
            channel: channel.into(),
            angle,
            wavelengths,
            values,
            sigma: vec![sigma],
            weight: 1.0,
            window: None,
        }
    }

    pub fn sigma_at(&self, i: usize) -> f64 {
        if self.sigma.len() == 1 {
            self.sigma[0]
        } else {
            self.sigma[i]
        }
    }

    /// Indices of the points inside the window.
    pub(crate) fn kept(&self) -> Vec<usize> {
        kept(&self.wavelengths, self.window)
    }

    /// Refuse, naming the dataset, anything that cannot make residuals.
    pub fn check(&self) -> Result<(), String> {
        let bad = |why: String| format!("dataset {:?}: {why}", self.id);
        channel_request(&self.channel).map_err(bad)?;
        let n = self.wavelengths.len();
        if n == 0 {
            return Err(bad("no points".into()));
        }
        if self.values.len() != n {
            return Err(bad(format!(
                "{} values for {n} wavelengths",
                self.values.len()
            )));
        }
        if self.sigma.len() != 1 && self.sigma.len() != n {
            return Err(bad(format!(
                "{} sigmas for {n} points; give one, or one per point",
                self.sigma.len()
            )));
        }
        if let Some(s) = self.sigma.iter().find(|s| !(**s > 0.0 && s.is_finite())) {
            return Err(bad(format!("sigma {s} is not a positive finite number")));
        }
        if let Some(w) = self
            .wavelengths
            .iter()
            .find(|w| !(**w > 0.0 && w.is_finite()))
        {
            return Err(bad(format!(
                "wavelength {w} is not a positive finite number"
            )));
        }
        if let Some(v) = self.values.iter().find(|v| !v.is_finite()) {
            return Err(bad(format!("value {v} is not finite")));
        }
        if !(self.weight > 0.0 && self.weight.is_finite()) {
            return Err(bad(format!(
                "weight {} is not a positive finite number",
                self.weight
            )));
        }
        if let Some((lo, hi)) = self.window
            && !(lo <= hi)
        {
            return Err(bad(format!("window ({lo}, {hi}) is empty")));
        }
        if self.kept().is_empty() {
            return Err(bad("no point inside the window".into()));
        }
        Ok(())
    }
}

fn kept(wavelengths: &[f64], window: Option<(f64, f64)>) -> Vec<usize> {
    (0..wavelengths.len())
        .filter(|&i| match window {
            Some((lo, hi)) => lo <= wavelengths[i] && wavelengths[i] <= hi,
            None => true,
        })
        .collect()
}

/// What one column of an [`NkData`] holds, from `n + ik` (`ε = (n + ik)²`,
/// so `ε₂ = 2nk >= 0` for an absorbing material).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NkQuantity {
    N,
    K,
    Eps1,
    Eps2,
}

impl NkQuantity {
    pub fn of(self, nk: Complex64) -> f64 {
        match self {
            NkQuantity::N => nk.re,
            NkQuantity::K => nk.im,
            NkQuantity::Eps1 => (nk * nk).re,
            NkQuantity::Eps2 => (nk * nk).im,
        }
    }
}

/// One measured optical-constant curve of an [`NkData`].
#[derive(Debug, Clone, PartialEq)]
pub struct NkColumn {
    pub quantity: NkQuantity,
    pub values: Vec<f64>,
    /// One value for every point, or one per point.
    pub sigma: Vec<f64>,
}

impl NkColumn {
    pub fn new(quantity: NkQuantity, values: Vec<f64>, sigma: f64) -> Self {
        Self {
            quantity,
            values,
            sigma: vec![sigma],
        }
    }

    pub fn sigma_at(&self, i: usize) -> f64 {
        if self.sigma.len() == 1 {
            self.sigma[0]
        } else {
            self.sigma[i]
        }
    }
}

/// Measured optical constants of one material of the model: n and k, ε₁
/// and ε₂, or any one of them, on their own wavelengths. Rows are
/// `sqrt(weight) * (sim - value) / sigma`, column after column.
#[derive(Debug, Clone, PartialEq)]
pub struct NkData {
    pub id: String,
    /// The model's material these values describe.
    pub material: String,
    /// nm.
    pub wavelengths: Vec<f64>,
    pub columns: Vec<NkColumn>,
    pub weight: f64,
    /// Keep only points with `lo <= wavelength <= hi` (nm).
    pub window: Option<(f64, f64)>,
}

impl NkData {
    /// Any columns, weight 1, no window.
    pub fn new(
        id: impl Into<String>,
        material: impl Into<String>,
        wavelengths: Vec<f64>,
        columns: Vec<NkColumn>,
    ) -> Self {
        Self {
            id: id.into(),
            material: material.into(),
            wavelengths,
            columns,
            weight: 1.0,
            window: None,
        }
    }

    /// n and k, one sigma each.
    pub fn nk(
        id: impl Into<String>,
        material: impl Into<String>,
        wavelengths: Vec<f64>,
        n: Vec<f64>,
        k: Vec<f64>,
        sigma_n: f64,
        sigma_k: f64,
    ) -> Self {
        let columns = vec![
            NkColumn::new(NkQuantity::N, n, sigma_n),
            NkColumn::new(NkQuantity::K, k, sigma_k),
        ];
        Self::new(id, material, wavelengths, columns)
    }

    /// ε₁ and ε₂, one sigma each.
    pub fn eps(
        id: impl Into<String>,
        material: impl Into<String>,
        wavelengths: Vec<f64>,
        eps1: Vec<f64>,
        eps2: Vec<f64>,
        sigma1: f64,
        sigma2: f64,
    ) -> Self {
        let columns = vec![
            NkColumn::new(NkQuantity::Eps1, eps1, sigma1),
            NkColumn::new(NkQuantity::Eps2, eps2, sigma2),
        ];
        Self::new(id, material, wavelengths, columns)
    }

    /// Indices of the points inside the window.
    pub(crate) fn kept(&self) -> Vec<usize> {
        kept(&self.wavelengths, self.window)
    }

    /// Refuse, naming the dataset, anything that cannot make residuals.
    pub fn check(&self) -> Result<(), String> {
        let bad = |why: String| format!("nk dataset {:?}: {why}", self.id);
        let n = self.wavelengths.len();
        if n == 0 {
            return Err(bad("no points".into()));
        }
        if self.columns.is_empty() {
            return Err(bad("no columns".into()));
        }
        for (c, col) in self.columns.iter().enumerate() {
            if self.columns[..c].iter().any(|o| o.quantity == col.quantity) {
                return Err(bad(format!("{:?} is given twice", col.quantity)));
            }
            if col.values.len() != n {
                return Err(bad(format!(
                    "{} {:?} values for {n} wavelengths",
                    col.values.len(),
                    col.quantity
                )));
            }
            if col.sigma.len() != 1 && col.sigma.len() != n {
                return Err(bad(format!(
                    "{} {:?} sigmas for {n} points; give one, or one per point",
                    col.sigma.len(),
                    col.quantity
                )));
            }
            if let Some(s) = col.sigma.iter().find(|s| !(**s > 0.0 && s.is_finite())) {
                return Err(bad(format!("sigma {s} is not a positive finite number")));
            }
            if let Some(v) = col.values.iter().find(|v| !v.is_finite()) {
                return Err(bad(format!("value {v} is not finite")));
            }
        }
        if let Some(w) = self
            .wavelengths
            .iter()
            .find(|w| !(**w > 0.0 && w.is_finite()))
        {
            return Err(bad(format!(
                "wavelength {w} is not a positive finite number"
            )));
        }
        if !(self.weight > 0.0 && self.weight.is_finite()) {
            return Err(bad(format!(
                "weight {} is not a positive finite number",
                self.weight
            )));
        }
        if let Some((lo, hi)) = self.window
            && !(lo <= hi)
        {
            return Err(bad(format!("window ({lo}, {hi}) is empty")));
        }
        if self.kept().is_empty() {
            return Err(bad("no point inside the window".into()));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrap_lands_in_half_open_interval() {
        assert_eq!(wrap_angle(PI), PI);
        assert_eq!(wrap_angle(-PI), PI);
        assert!((wrap_angle(2.0 * PI - 0.1) + 0.1).abs() < 1e-12);
        let d = (179.0f64 - (-179.0)).to_radians();
        assert!((wrap_angle(d) + 2.0f64.to_radians()).abs() < 1e-12);
    }

    #[test]
    fn every_channel_names_one_bit() {
        let mut seen = 0u64;
        for (name, bit) in CHANNELS {
            assert_eq!(bit.count_ones(), 1, "{name}");
            assert_eq!(seen & bit, 0, "{name} reuses a bit");
            seen |= bit;
        }
        assert!(channel_request("Delta_R").is_ok());
        assert!(channel_request("R").unwrap_err().contains("Rs"));
    }

    #[test]
    fn bad_data_names_the_dataset() {
        let ok = SpectrumData::new(
            "R8",
            "Rs",
            Radians::ZERO,
            vec![500.0, 600.0],
            vec![0.1, 0.2],
            0.01,
        );
        ok.check().unwrap();
        let mut d = ok.clone();
        d.sigma = vec![0.0];
        assert!(d.check().unwrap_err().contains("\"R8\""));
        let mut d = ok.clone();
        d.values.pop();
        assert!(d.check().unwrap_err().contains("1 values for 2"));
        let mut d = ok.clone();
        d.window = Some((700.0, 800.0));
        assert!(d.check().unwrap_err().contains("window"));
        let mut d = ok;
        d.channel = "Reflectance".into();
        assert!(d.check().unwrap_err().contains("unknown channel"));
    }
}
