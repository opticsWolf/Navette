// SPDX-License-Identifier: LGPL-3.0-or-later
//! synthesis::weights — per-point weight profiles (synthesis targets plan W4).
//!
//! A target's point weights multiply its scalar weight point by point:
//! a pointwise residual becomes `r_i·√(weight·w_i/count)`, an integral
//! target constrains the weighted mean `Σw_i d_i / Σw_i` over
//! `Σw_i tol_i / Σw_i`. Tolerances keep their meaning (scale and band);
//! weights carry emphasis.
//!
//! A target gives either explicit `weights` (one per point) or a
//! [`WeightProfile`], inline or by name from the set's `weight_profiles`.
//! The profile is evaluated on the target's own axis: wavelength in nm
//! for a spectral target, angle in degrees for an angular one.

use serde::{Deserialize, Serialize};

use crate::interpolate::UniInterpolator;

/// A weight as a function of the target's axis (nm, or degrees).
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum WeightProfile {
    /// The same weight everywhere.
    Constant { value: f64 },
    /// Linear from `w_lo` at `lo` to `w_hi` at `hi`, held flat outside.
    Ramp {
        lo: f64,
        hi: f64,
        w_lo: f64,
        w_hi: f64,
    },
    /// `base + (peak − base)·exp(−4 ln2 (x − center)²/fwhm²)`.
    Gaussian {
        center: f64,
        fwhm: f64,
        peak: f64,
        base: f64,
    },
    /// Weight `w` on each closed `[lo, hi]`, `outside` elsewhere. Bands
    /// may touch but not overlap; at a shared edge the lower band wins.
    Bands {
        bands: Vec<(f64, f64, f64)>,
        outside: f64,
    },
    /// Tabulated weights, interpolated (`pchip`, which does not overshoot,
    /// or `linear`) and held flat past the ends.
    Table {
        x: Vec<f64>,
        w: Vec<f64>,
        #[serde(default = "d_interp")]
        interp: String,
    },
    /// Pointwise product, e.g. an emphasis times a window.
    Product { factors: Vec<WeightProfile> },
}

fn d_interp() -> String {
    "pchip".to_string()
}

fn finite(name: &str, v: f64) -> Result<(), String> {
    if v.is_finite() {
        Ok(())
    } else {
        Err(format!("{name} must be finite (got {v})"))
    }
}

fn weight(name: &str, v: f64) -> Result<(), String> {
    if v.is_finite() && v >= 0.0 {
        Ok(())
    } else {
        Err(format!("{name} must be finite and >= 0 (got {v})"))
    }
}

impl WeightProfile {
    /// Check the parameters alone (no axis needed).
    pub fn validate(&self) -> Result<(), String> {
        match self {
            WeightProfile::Constant { value } => weight("constant: value", *value),
            WeightProfile::Ramp { lo, hi, w_lo, w_hi } => {
                finite("ramp: lo", *lo)?;
                finite("ramp: hi", *hi)?;
                if lo >= hi {
                    return Err(format!("ramp: lo < hi required (got {lo}, {hi})"));
                }
                weight("ramp: w_lo", *w_lo)?;
                weight("ramp: w_hi", *w_hi)
            }
            WeightProfile::Gaussian {
                center,
                fwhm,
                peak,
                base,
            } => {
                finite("gaussian: center", *center)?;
                if !(fwhm.is_finite() && *fwhm > 0.0) {
                    return Err(format!(
                        "gaussian: fwhm must be finite and > 0 (got {fwhm})"
                    ));
                }
                weight("gaussian: peak", *peak)?;
                weight("gaussian: base", *base)
            }
            WeightProfile::Bands { bands, outside } => {
                weight("bands: outside", *outside)?;
                for (i, &(lo, hi, w)) in bands.iter().enumerate() {
                    finite(&format!("bands[{i}]: lo"), lo)?;
                    finite(&format!("bands[{i}]: hi"), hi)?;
                    if lo >= hi {
                        return Err(format!("bands[{i}]: lo < hi required (got {lo}, {hi})"));
                    }
                    weight(&format!("bands[{i}]: w"), w)?;
                }
                let mut sorted: Vec<&(f64, f64, f64)> = bands.iter().collect();
                sorted.sort_by(|a, b| a.0.total_cmp(&b.0));
                for pair in sorted.windows(2) {
                    if pair[1].0 < pair[0].1 {
                        return Err(format!(
                            "bands: [{}, {}] and [{}, {}] overlap",
                            pair[0].0, pair[0].1, pair[1].0, pair[1].1
                        ));
                    }
                }
                Ok(())
            }
            WeightProfile::Table { x, w, interp } => {
                if x.len() < 2 || x.len() != w.len() {
                    return Err(format!(
                        "table: x and w need the same length, at least 2 (got {} and {})",
                        x.len(),
                        w.len()
                    ));
                }
                for (i, &v) in x.iter().enumerate() {
                    finite(&format!("table: x[{i}]"), v)?;
                }
                if x.windows(2).any(|p| p[1] <= p[0]) {
                    return Err("table: x must be strictly increasing".into());
                }
                for (i, &v) in w.iter().enumerate() {
                    weight(&format!("table: w[{i}]"), v)?;
                }
                if !matches!(interp.as_str(), "pchip" | "linear") {
                    return Err(format!(
                        "table: interp must be \"pchip\" or \"linear\" (got {interp:?})"
                    ));
                }
                Ok(())
            }
            WeightProfile::Product { factors } => {
                if factors.is_empty() {
                    return Err("product: needs at least one factor".into());
                }
                factors.iter().try_for_each(WeightProfile::validate)
            }
        }
    }

    /// The weight at each `x` (nm for spectral targets, degrees for
    /// angular ones). Validates first; every value is finite and >= 0.
    pub fn evaluate(&self, x: &[f64]) -> Result<Vec<f64>, String> {
        self.validate()?;
        let w = self.evaluate_unchecked(x)?;
        if let Some((i, v)) = w
            .iter()
            .enumerate()
            .find(|(_, v)| !(v.is_finite() && **v >= 0.0))
        {
            return Err(format!("weight profile gave {v} at x = {}", x[i]));
        }
        Ok(w)
    }

    fn evaluate_unchecked(&self, x: &[f64]) -> Result<Vec<f64>, String> {
        Ok(match self {
            WeightProfile::Constant { value } => vec![*value; x.len()],
            WeightProfile::Ramp { lo, hi, w_lo, w_hi } => x
                .iter()
                .map(|&v| {
                    let f = ((v - lo) / (hi - lo)).clamp(0.0, 1.0);
                    w_lo + (w_hi - w_lo) * f
                })
                .collect(),
            WeightProfile::Gaussian {
                center,
                fwhm,
                peak,
                base,
            } => {
                let a = 4.0 * std::f64::consts::LN_2 / (fwhm * fwhm);
                x.iter()
                    .map(|&v| base + (peak - base) * (-a * (v - center).powi(2)).exp())
                    .collect()
            }
            WeightProfile::Bands { bands, outside } => {
                let mut sorted: Vec<&(f64, f64, f64)> = bands.iter().collect();
                sorted.sort_by(|a, b| a.0.total_cmp(&b.0));
                x.iter()
                    .map(|&v| {
                        sorted
                            .iter()
                            .find(|(lo, hi, _)| *lo <= v && v <= *hi)
                            .map_or(*outside, |b| b.2)
                    })
                    .collect()
            }
            WeightProfile::Table { x: tx, w, interp } => {
                let y = ndarray017::Array2::from_shape_vec((1, w.len()), w.clone())
                    .map_err(|e| e.to_string())?;
                let it = UniInterpolator::new(
                    ndarray017::Array1::from(tx.clone()),
                    y,
                    false,
                    interp,
                    false,
                    0,
                    "clamp",
                )?;
                it.evaluate(x, 0, None)?.row(0).to_vec()
            }
            WeightProfile::Product { factors } => {
                let mut out = vec![1.0; x.len()];
                for f in factors {
                    for (o, v) in out.iter_mut().zip(f.evaluate_unchecked(x)?) {
                        *o *= v;
                    }
                }
                out
            }
        })
    }
}

/// A target's `weight_profile`: inline, or a name in the set's
/// `weight_profiles`.
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
#[serde(untagged)]
pub enum ProfileRef {
    Name(String),
    Inline(WeightProfile),
}

/// Check resolved point weights: one per point, finite, >= 0, not all
/// zero (a target is muted with `weight = 0`, not with zero weights).
pub(crate) fn check_point_weights(w: &[f64], n: usize, label: &str) -> Result<(), String> {
    if w.len() != n {
        return Err(format!(
            "{label}: weights shape mismatch: {} != {n}.",
            w.len()
        ));
    }
    if let Some(v) = w.iter().find(|v| !(v.is_finite() && **v >= 0.0)) {
        return Err(format!(
            "{label}: weights must be finite and >= 0 (got {v})."
        ));
    }
    if w.iter().all(|v| *v == 0.0) {
        return Err(format!(
            "{label}: weights are all zero; mute a target with weight=0 instead."
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn eval(json: &str, x: &[f64]) -> Result<Vec<f64>, String> {
        serde_json::from_str::<WeightProfile>(json)
            .map_err(|e| e.to_string())?
            .evaluate(x)
    }

    #[test]
    fn each_profile_evaluates() {
        let x = [400.0, 450.0, 500.0, 550.0, 600.0];
        assert_eq!(
            eval(r#"{"kind":"constant","value":2}"#, &x).unwrap(),
            vec![2.0; 5]
        );
        assert_eq!(
            eval(r#"{"kind":"ramp","lo":450,"hi":550,"w_lo":1,"w_hi":3}"#, &x).unwrap(),
            vec![1.0, 1.0, 2.0, 3.0, 3.0]
        );
        let g = eval(
            r#"{"kind":"gaussian","center":500,"fwhm":100,"peak":5,"base":1}"#,
            &x,
        )
        .unwrap();
        assert_eq!(g[2], 5.0);
        assert!((g[1] - 3.0).abs() < 1e-12 && (g[3] - 3.0).abs() < 1e-12);
        assert_eq!(
            eval(
                r#"{"kind":"bands","bands":[[450,500,4],[500,550,2]],"outside":1}"#,
                &x
            )
            .unwrap(),
            vec![1.0, 4.0, 4.0, 2.0, 1.0]
        );
        let lin = eval(
            r#"{"kind":"table","x":[400,600],"w":[1,5],"interp":"linear"}"#,
            &[300.0, 500.0, 700.0],
        )
        .unwrap();
        assert_eq!(lin, vec![1.0, 3.0, 5.0]);
        let p = eval(
            r#"{"kind":"product","factors":[{"kind":"constant","value":2},
                {"kind":"ramp","lo":450,"hi":550,"w_lo":0,"w_hi":1}]}"#,
            &x,
        )
        .unwrap();
        assert_eq!(p, vec![0.0, 0.0, 1.0, 2.0, 2.0]);
    }

    #[test]
    fn pchip_tables_stay_inside_the_data() {
        let w = eval(
            r#"{"kind":"table","x":[0,1,2,3],"w":[0,0,4,4]}"#,
            &(0..=30).map(|i| i as f64 / 10.0).collect::<Vec<_>>(),
        )
        .unwrap();
        assert!(w.iter().all(|v| (0.0..=4.0).contains(v)));
    }

    #[test]
    fn bad_profiles_refuse_by_name() {
        let x = [500.0];
        for (json, want) in [
            (
                r#"{"kind":"ramp","lo":5,"hi":5,"w_lo":1,"w_hi":1}"#,
                "ramp: lo < hi",
            ),
            (
                r#"{"kind":"gaussian","center":0,"fwhm":0,"peak":1,"base":0}"#,
                "fwhm",
            ),
            (r#"{"kind":"constant","value":-1}"#, "constant: value"),
            (
                r#"{"kind":"bands","bands":[[0,10,1],[5,20,1]],"outside":0}"#,
                "overlap",
            ),
            (
                r#"{"kind":"table","x":[0,0],"w":[1,1]}"#,
                "strictly increasing",
            ),
            (
                r#"{"kind":"table","x":[0,1],"w":[1,1],"interp":"cubic"}"#,
                "interp",
            ),
            (r#"{"kind":"product","factors":[]}"#, "at least one"),
        ] {
            let e = eval(json, &x).unwrap_err();
            assert!(e.contains(want), "{json}: {e}");
        }
        assert!(eval(r#"{"kind":"ramp","lo":0,"hi":1,"w_lo":1}"#, &x).is_err());
        assert!(eval(r#"{"kind":"constant","value":1,"extra":2}"#, &x).is_err());
    }

    #[test]
    fn point_weights_refuse_bad_shapes() {
        assert!(check_point_weights(&[1.0, 2.0], 2, "t").is_ok());
        assert!(
            check_point_weights(&[1.0], 2, "t")
                .unwrap_err()
                .contains("shape")
        );
        assert!(
            check_point_weights(&[1.0, -1.0], 2, "t")
                .unwrap_err()
                .contains(">= 0")
        );
        assert!(
            check_point_weights(&[0.0, 0.0], 2, "t")
                .unwrap_err()
                .contains("all zero")
        );
    }
}
