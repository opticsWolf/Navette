// SPDX-License-Identifier: LGPL-3.0-or-later
//! Scaling/error policy shared by a set of layers.
//!
//! Mirrors `navette.structure.models.Group` exactly (factors, summands,
//! masks, error laws + params, floors, state keys). Deliberate differences:
//! - Randomness is full-Rust (§9.2): draws take `&mut dyn RngCore`
//!   (`StdRng::seed_from_u64` = reproducible, `rand::rng()` = thread).
//!   Streams differ from NumPy by algorithm (ChaCha vs PCG64) —
//!   acceptance is statistical agreement + per-side determinism.
//! - Unknown error laws are unrepresentable (typed `ErrorType`); Python
//!   falls through to unperturbed values. Fail-closed by construction.
//! - `validate` returns typed issues; `set_properties` returns warnings.
//! - States require complete param maps (Python replaces whole dicts and
//!   fails later at draw time; refusing at load is fail-closed).

use std::collections::BTreeMap;
use std::fmt;

use num_complex::Complex64;
use rand::RngCore;
use rand_distr::{Distribution, Normal, Uniform};
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::structure::enums::ErrorType;
use crate::structure::validation::ValidationIssue;
use crate::structure::version::{SCHEMA_VERSION, check_schema_version};

fn one() -> f64 {
    1.0
}

/// Eight-parameter error vocabulary, shared by every channel.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ErrorParams {
    pub abs_mean_delta_g: f64,
    pub abs_std_dev: f64,
    pub rel_mean_delta_g: f64,
    pub rel_std_dev: f64,
    pub abs_mean_delta_h: f64,
    pub abs_variance: f64,
    pub rel_mean_delta_h: f64,
    pub rel_variance: f64,
}

impl ErrorParams {
    /// Default law params (all channels except roughness).
    ///
    /// The relative spreads are unit-free *fractions* of the value, not
    /// percentages: `rel_std_dev: 0.01` is a 1% one-sigma relative scatter,
    /// matching the absolute channel's 0.01 nm. They read 1.0 from the first
    /// Python upload through 0.7.16 -- a 100% one-sigma spread, one sigma
    /// covering the entire nominal thickness, which is not a fabrication
    /// tolerance anyone has. It stayed invisible because `error_mask`
    /// defaults to all-zero, so no channel draws at all until it is switched
    /// on, at which point the first draw was wild.
    pub fn standard() -> Self {
        Self {
            abs_mean_delta_g: 0.0,
            abs_std_dev: 0.01,
            rel_mean_delta_g: 0.0,
            rel_std_dev: 0.01,
            abs_mean_delta_h: 0.0,
            abs_variance: 0.01,
            rel_mean_delta_h: 0.0,
            rel_variance: 0.01,
        }
    }

    /// Roughness-channel defaults: `abs_*` x0.1 so the physical magnitude is
    /// unchanged by the Å→nm switch (0.01 Å == 0.001 nm); `rel_*` untouched,
    /// since a fraction has no unit to convert.
    pub fn roughness() -> Self {
        Self {
            abs_mean_delta_g: 0.0,
            abs_std_dev: 0.001,
            rel_mean_delta_g: 0.0,
            rel_std_dev: 0.01,
            abs_mean_delta_h: 0.0,
            abs_variance: 0.001,
            rel_mean_delta_h: 0.0,
            rel_variance: 0.01,
        }
    }

    /// Extinction-channel defaults: `abs_*` x0.01 relative to `standard()`.
    ///
    /// The absolute spreads carry the unit of the quantity they perturb, and
    /// `standard()`'s 0.01 is sized for a thickness in nanometres. Applied to
    /// `k`, which runs from about 1e-4 to 1e-2 in the visible, a 0.01
    /// absolute scatter is one to two orders of magnitude *larger than the
    /// value*, so roughly half of all draws land at `k < 0`. That is optical
    /// gain, which the solver door refuses outright
    /// (`GAIN_MEDIUM_EXPLANATION`) -- so enabling the k error channel used to
    /// abort a tolerance run rather than perturb it, and the refusal named
    /// gain rather than the tolerance that caused it.
    ///
    /// 0.0001 keeps the absolute scatter below a typical `k`. The relative
    /// channel needed nothing: `k * (1 + g)` scales with the value and was
    /// already well behaved.
    pub fn extinction() -> Self {
        Self {
            abs_mean_delta_g: 0.0,
            abs_std_dev: 0.0001,
            rel_mean_delta_g: 0.0,
            rel_std_dev: 0.01,
            abs_mean_delta_h: 0.0,
            abs_variance: 0.0001,
            rel_mean_delta_h: 0.0,
            rel_variance: 0.01,
        }
    }
}

/// Wrap a group in a shared handle (structure storage).
pub fn shared_group(group: Group) -> crate::structure::SharedGroup {
    std::rc::Rc::new(std::cell::RefCell::new(group))
}

/// Scaling/error policy shared by a set of layers (material-keyed).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Group {
    pub group_name: String,
    pub thick_factor: f64,
    pub thick_summand: f64,
    pub n_factor: f64,
    pub k_factor: f64,
    pub inh_delta_summand: f64,
    pub roughness_summand: f64,
    pub interface_summand: f64,
    pub error_mask: [i32; 6],
    pub optimization_mask: [i32; 7],
    pub thickness_error_type: ErrorType,
    pub n_error_type: ErrorType,
    pub k_error_type: ErrorType,
    pub inh_delta_error_type: ErrorType,
    pub roughness_error_type: ErrorType,
    pub interface_error_type: ErrorType,
    pub thickness_error_params: ErrorParams,
    pub inh_delta_error_params: ErrorParams,
    pub roughness_error_params: ErrorParams,
    pub interface_error_params: ErrorParams,
    pub n_error_params: ErrorParams,
    pub k_error_params: ErrorParams,
}

impl Group {
    /// Named group at identity (factors (1,1), zero summands, GAUSSIAN laws).
    pub fn new(group_name: impl Into<String>) -> Self {
        Self {
            group_name: group_name.into(),
            thick_factor: 1.0,
            thick_summand: 0.0,
            n_factor: 1.0,
            k_factor: 1.0,
            inh_delta_summand: 0.0,
            roughness_summand: 0.0,
            interface_summand: 0.0,
            error_mask: [0; 6],
            optimization_mask: [1; 7],
            thickness_error_type: ErrorType::Gaussian,
            n_error_type: ErrorType::Gaussian,
            k_error_type: ErrorType::Gaussian,
            inh_delta_error_type: ErrorType::Gaussian,
            roughness_error_type: ErrorType::Gaussian,
            interface_error_type: ErrorType::Gaussian,
            thickness_error_params: ErrorParams::standard(),
            inh_delta_error_params: ErrorParams::standard(),
            roughness_error_params: ErrorParams::roughness(),
            interface_error_params: ErrorParams::standard(),
            n_error_params: ErrorParams::standard(),
            k_error_params: ErrorParams::extinction(),
        }
    }

    /// `(n_factor, k_factor)` as a complex multiplier (Python `nk_factor`).
    pub fn nk_factor(&self) -> Complex64 {
        Complex64::new(self.n_factor, self.k_factor)
    }

    /// Factor domains: identity (1,1); negatives unphysical; NaN invalid;
    /// optimization mask must be 7 binary entries.
    pub fn validate(&self) -> Vec<ValidationIssue> {
        let mut issues = Vec::new();
        if self.thick_factor < 0.0 || self.thick_factor.is_nan() {
            issues.push(ValidationIssue::error(format!(
                "Group '{}': thick_factor {} < 0 (NaN counts as invalid).",
                self.group_name, self.thick_factor
            )));
        }
        if self.n_factor < 0.0 || self.n_factor.is_nan() {
            issues.push(ValidationIssue::error(format!(
                "Group '{}': n_factor {} < 0 (no negative-index media).",
                self.group_name, self.n_factor
            )));
        }
        if self.k_factor < 0.0 || self.k_factor.is_nan() {
            issues.push(ValidationIssue::error(format!(
                "Group '{}': k_factor {} < 0 (no gain media).",
                self.group_name, self.k_factor
            )));
        }
        if self.optimization_mask.iter().any(|v| *v != 0 && *v != 1) {
            issues.push(ValidationIssue::error(format!(
                "Group '{}': optimization_mask must be 7 binary entries (see OptMask).",
                self.group_name
            )));
        }
        issues
    }

    /// Draw one perturbation (Python `_apply_error`, order-preserved).
    /// Non-positive spreads contribute their mean deterministically (NumPy
    /// draws exactly at zero spread; no RNG consumed either way that matters
    /// — streams are per-side by §9.2).
    ///
    /// Both laws are `(centre, spread)`: Gaussian is
    /// `N(*_mean_delta_g, *_std_dev)`, uniform is
    /// `U(*_mean_delta_h ± *_variance)`. The uniform centre was dropped by
    /// every version before 0.7.15 — see `unif_draw`.
    ///
    /// `Cascaded` draws the same four numbers as `Combined`, in the same
    /// order, and composes the two relative laws as factors instead of
    /// summing them — see [`ErrorType`] for the table and the reasoning.
    pub fn apply_error(
        value: f64,
        error_type: ErrorType,
        params: &ErrorParams,
        rng: &mut dyn RngCore,
    ) -> f64 {
        let g_abs = gauss_draw(params.abs_mean_delta_g, params.abs_std_dev, rng);
        match error_type {
            ErrorType::Gaussian => {
                value + g_abs + gauss_draw(params.rel_mean_delta_g, params.rel_std_dev, rng) * value
            }
            ErrorType::Uniform => {
                value
                    + unif_draw(params.abs_mean_delta_h, params.abs_variance, rng)
                    + unif_draw(params.rel_mean_delta_h, params.rel_variance, rng) * value
            }
            ErrorType::Combined => {
                value
                    + g_abs
                    + gauss_draw(params.rel_mean_delta_g, params.rel_std_dev, rng) * value
                    + unif_draw(params.abs_mean_delta_h, params.abs_variance, rng)
                    + unif_draw(params.rel_mean_delta_h, params.rel_variance, rng) * value
            }
            ErrorType::Cascaded => {
                // Same four draws as Combined, in the same order, so the two
                // laws stay stream-comparable; only the composition differs.
                let g_rel = gauss_draw(params.rel_mean_delta_g, params.rel_std_dev, rng);
                let u_abs = unif_draw(params.abs_mean_delta_h, params.abs_variance, rng);
                let u_rel = unif_draw(params.rel_mean_delta_h, params.rel_variance, rng);
                value * (1.0 + g_rel) * (1.0 + u_rel) + g_abs + u_abs
            }
        }
    }

    /// Perturbed thickness (floored at 0).
    pub fn thickness_error(&self, value: f64, rng: &mut dyn RngCore) -> f64 {
        Self::apply_error(
            value,
            self.thickness_error_type,
            &self.thickness_error_params,
            rng,
        )
        .max(0.0)
    }

    /// Perturbed grading strength (unfloored, like Python).
    pub fn inh_delta_error(&self, value: f64, rng: &mut dyn RngCore) -> f64 {
        Self::apply_error(
            value,
            self.inh_delta_error_type,
            &self.inh_delta_error_params,
            rng,
        )
    }

    /// Perturbed surface roughness [nm] (floored at 0).
    pub fn sr_roughness_error(&self, value: f64, rng: &mut dyn RngCore) -> f64 {
        Self::apply_error(
            value,
            self.roughness_error_type,
            &self.roughness_error_params,
            rng,
        )
        .max(0.0)
    }

    /// Perturbed interface width [nm] (floored at 0).
    pub fn interface_error(&self, value: f64, rng: &mut dyn RngCore) -> f64 {
        Self::apply_error(
            value,
            self.interface_error_type,
            &self.interface_error_params,
            rng,
        )
        .max(0.0)
    }

    /// Perturbed index with n floored at 0 (k untouched by the floor).
    pub fn nk_error(&self, nk_value: Complex64, rng: &mut dyn RngCore) -> Complex64 {
        let n =
            Self::apply_error(nk_value.re, self.n_error_type, &self.n_error_params, rng).max(0.0);
        let k = Self::apply_error(nk_value.im, self.k_error_type, &self.k_error_params, rng);
        Complex64::new(n, k)
    }

    /// Serialize all slots (Python `get_state`: schema_version + slots).
    pub fn to_state(&self) -> Value {
        let mut v = serde_json::to_value(self).expect("Group serialization is infallible");
        v.as_object_mut()
            .expect("Group serializes as a map")
            .insert("schema_version".to_string(), Value::from(SCHEMA_VERSION));
        v
    }

    /// Rebuild from state (Python `from_state`): version-checked, unknown
    /// keys ignored, deep-copied by value. Missing `group_name` → "default".
    pub fn from_state(value: &Value) -> Result<Self, String> {
        let found = value
            .get("schema_version")
            .and_then(|v| v.as_u64())
            .map(|v| v as u32);
        check_schema_version(found, "Group")?;
        #[derive(Deserialize)]
        struct Raw {
            #[serde(default = "default_group_name")]
            group_name: String,
            #[serde(default = "one")]
            thick_factor: f64,
            #[serde(default)]
            thick_summand: f64,
            #[serde(default = "one")]
            n_factor: f64,
            #[serde(default = "one")]
            k_factor: f64,
            #[serde(default)]
            inh_delta_summand: f64,
            #[serde(default)]
            roughness_summand: f64,
            #[serde(default)]
            interface_summand: f64,
            #[serde(default)]
            error_mask: [i32; 6],
            #[serde(default = "default_opt_mask")]
            optimization_mask: [i32; 7],
            #[serde(default = "gaussian")]
            thickness_error_type: ErrorType,
            #[serde(default = "gaussian")]
            n_error_type: ErrorType,
            #[serde(default = "gaussian")]
            k_error_type: ErrorType,
            #[serde(default = "gaussian")]
            inh_delta_error_type: ErrorType,
            #[serde(default = "gaussian")]
            roughness_error_type: ErrorType,
            #[serde(default = "gaussian")]
            interface_error_type: ErrorType,
            #[serde(default = "ErrorParams::standard")]
            thickness_error_params: ErrorParams,
            #[serde(default = "ErrorParams::standard")]
            inh_delta_error_params: ErrorParams,
            #[serde(default = "ErrorParams::roughness")]
            roughness_error_params: ErrorParams,
            #[serde(default = "ErrorParams::standard")]
            interface_error_params: ErrorParams,
            #[serde(default = "ErrorParams::standard")]
            n_error_params: ErrorParams,
            #[serde(default = "ErrorParams::standard")]
            k_error_params: ErrorParams,
        }
        fn default_group_name() -> String {
            "default".to_string()
        }
        fn default_opt_mask() -> [i32; 7] {
            [1; 7]
        }
        fn gaussian() -> ErrorType {
            ErrorType::Gaussian
        }
        let raw: Raw = serde_json::from_value(value.clone())
            .map_err(|e| format!("Group: malformed state ({e})."))?;
        Ok(Self {
            group_name: raw.group_name,
            thick_factor: raw.thick_factor,
            thick_summand: raw.thick_summand,
            n_factor: raw.n_factor,
            k_factor: raw.k_factor,
            inh_delta_summand: raw.inh_delta_summand,
            roughness_summand: raw.roughness_summand,
            interface_summand: raw.interface_summand,
            error_mask: raw.error_mask,
            optimization_mask: raw.optimization_mask,
            thickness_error_type: raw.thickness_error_type,
            n_error_type: raw.n_error_type,
            k_error_type: raw.k_error_type,
            inh_delta_error_type: raw.inh_delta_error_type,
            roughness_error_type: raw.roughness_error_type,
            interface_error_type: raw.interface_error_type,
            thickness_error_params: raw.thickness_error_params,
            inh_delta_error_params: raw.inh_delta_error_params,
            roughness_error_params: raw.roughness_error_params,
            interface_error_params: raw.interface_error_params,
            n_error_params: raw.n_error_params,
            k_error_params: raw.k_error_params,
        })
    }

    /// Bulk-set known properties; unknown keys become returned warnings.
    pub fn set_properties(&mut self, props: &BTreeMap<String, Value>) -> Vec<ValidationIssue> {
        let mut warnings = Vec::new();
        let bad = |msg: String| ValidationIssue::warning(format!("Group.set_properties: {msg}"));
        for (key, value) in props {
            match key.as_str() {
                "group_name" => match value.as_str() {
                    Some(s) => self.group_name = s.to_string(),
                    None => warnings.push(bad("ignoring non-string 'group_name'.".to_string())),
                },
                "thick_factor" | "thick_summand" | "n_factor" | "k_factor"
                | "inh_delta_summand" | "roughness_summand" | "interface_summand" => {
                    match value.as_f64() {
                        Some(v) => match key.as_str() {
                            "thick_factor" => self.thick_factor = v,
                            "thick_summand" => self.thick_summand = v,
                            "n_factor" => self.n_factor = v,
                            "k_factor" => self.k_factor = v,
                            "inh_delta_summand" => self.inh_delta_summand = v,
                            "roughness_summand" => self.roughness_summand = v,
                            _ => self.interface_summand = v,
                        },
                        None => warnings.push(bad(format!("ignoring non-numeric '{key}'."))),
                    }
                }
                "error_mask" => match parse_i32_array::<6>(value) {
                    Some(m) => self.error_mask = m,
                    None => warnings.push(bad("ignoring malformed 'error_mask'.".to_string())),
                },
                "optimization_mask" => match parse_i32_array::<7>(value) {
                    Some(m) => self.optimization_mask = m,
                    None => {
                        warnings.push(bad("ignoring malformed 'optimization_mask'.".to_string()))
                    }
                },
                t @ ("thickness_error_type"
                | "n_error_type"
                | "k_error_type"
                | "inh_delta_error_type"
                | "roughness_error_type"
                | "interface_error_type") => {
                    match value
                        .as_i64()
                        .map(|v| v as i32)
                        .map(ErrorType::try_from_i32)
                    {
                        Some(Ok(e)) => match t {
                            "thickness_error_type" => self.thickness_error_type = e,
                            "n_error_type" => self.n_error_type = e,
                            "k_error_type" => self.k_error_type = e,
                            "inh_delta_error_type" => self.inh_delta_error_type = e,
                            "roughness_error_type" => self.roughness_error_type = e,
                            _ => self.roughness_error_type = e,
                        },
                        _ => warnings.push(bad(format!("unknown {t}; ignoring."))),
                    }
                }
                p @ ("thickness_error_params"
                | "inh_delta_error_params"
                | "roughness_error_params"
                | "interface_error_params"
                | "n_error_params"
                | "k_error_params") => match serde_json::from_value::<ErrorParams>(value.clone()) {
                    Ok(ep) => match p {
                        "thickness_error_params" => self.thickness_error_params = ep,
                        "inh_delta_error_params" => self.inh_delta_error_params = ep,
                        "roughness_error_params" => self.roughness_error_params = ep,
                        "interface_error_params" => self.interface_error_params = ep,
                        "n_error_params" => self.n_error_params = ep,
                        _ => self.k_error_params = ep,
                    },
                    Err(_) => warnings.push(bad(format!("ignoring malformed '{p}'."))),
                },
                other => warnings.push(bad(format!("ignoring unknown attribute '{other}'."))),
            }
        }
        warnings
    }
}

pub(crate) fn gauss_draw<R: RngCore + ?Sized>(mean: f64, std: f64, rng: &mut R) -> f64 {
    if std <= 0.0 {
        mean
    } else {
        Normal::new(mean, std)
            .map(|d| d.sample(rng))
            .unwrap_or(mean)
    }
}

/// One uniform draw on `[mean - half_width, mean + half_width]`.
///
/// `mean` is the law's systematic offset, the uniform counterpart of
/// `gauss_draw`'s. It was declared as `abs_mean_delta_h` / `rel_mean_delta_h`
/// from the first Python upload and never read by any version, Python or
/// Rust: the draw was hard-centred on zero, so a configured bias was
/// accepted, stored, serialized and then silently dropped. Both defaults are
/// `0.0`, where `U(-w, w)` and `U(mean-w, mean+w)` are the same distribution,
/// which is why nothing ever noticed.
///
/// Zero or negative width contributes the mean deterministically and consumes
/// no RNG, exactly as `gauss_draw` does at zero spread.
pub(crate) fn unif_draw<R: RngCore + ?Sized>(mean: f64, half_width: f64, rng: &mut R) -> f64 {
    if half_width <= 0.0 {
        mean
    } else {
        Uniform::new(mean - half_width, mean + half_width)
            .map(|d| d.sample(rng))
            .unwrap_or(mean)
    }
}

fn parse_i32_array<const N: usize>(value: &Value) -> Option<[i32; N]> {
    let arr = value.as_array()?;
    if arr.len() != N {
        return None;
    }
    let mut out = [0i32; N];
    for (i, v) in arr.iter().enumerate() {
        out[i] = v.as_i64()? as i32;
    }
    Some(out)
}

impl fmt::Display for Group {
    /// Python `__repr__` format, verbatim.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Group(name='{}', thick_factor={:.3})",
            self.group_name, self.thick_factor
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use serde_json::json;

    fn rng() -> StdRng {
        StdRng::seed_from_u64(7)
    }

    #[test]
    fn defaults_match_python_ctor() {
        let g = Group::new("TiO2");
        assert_eq!((g.thick_factor, g.thick_summand), (1.0, 0.0));
        assert_eq!((g.n_factor, g.k_factor), (1.0, 1.0));
        assert_eq!(g.nk_factor(), Complex64::new(1.0, 1.0));
        assert_eq!(g.error_mask, [0; 6]);
        assert_eq!(g.optimization_mask, [1; 7]);
        assert_eq!(g.thickness_error_type, ErrorType::Gaussian);
        assert_eq!(g.thickness_error_params.abs_std_dev, 0.01);
        // Roughness abs defaults x0.1 (Å→nm magnitude preservation).
        assert_eq!(g.roughness_error_params.abs_std_dev, 0.001);
        assert_eq!(g.roughness_error_params.abs_variance, 0.001);
        assert_eq!(g.roughness_error_params.rel_std_dev, 0.01);
        assert_eq!(g.to_string(), "Group(name='TiO2', thick_factor=1.000)");
    }

    #[test]
    fn validate_domains_match_python() {
        assert!(Group::new("g").validate().is_empty());
        let mut g = Group::new("g");
        g.thick_factor = -1.0;
        g.n_factor = -0.5;
        g.k_factor = f64::NAN;
        g.optimization_mask[3] = 2;
        let issues = g.validate();
        assert_eq!(issues.len(), 4);
        assert!(issues.iter().all(|i| i.is_error()));
    }

    #[test]
    fn error_laws_agree_statistically() {
        // GAUSSIAN: abs N(1,2) only → mean ≈ 11, std ≈ 2 on value 10.
        let p = ErrorParams {
            abs_mean_delta_g: 1.0,
            abs_std_dev: 2.0,
            rel_mean_delta_g: 0.0,
            rel_std_dev: 0.0,
            ..ErrorParams::standard()
        };
        // Zero rel spread contributes its mean (0.0) deterministically.
        let mut r = rng();
        let n = 20_000;
        let mut sum = 0.0;
        let mut sum2 = 0.0;
        for _ in 0..n {
            let d = Group::apply_error(10.0, ErrorType::Gaussian, &p, &mut r);
            sum += d;
            sum2 += d * d;
        }
        let mean = sum / n as f64;
        let std = ((sum2 / n as f64) - mean * mean).sqrt();
        assert!((mean - 11.0).abs() < 0.1, "mean {mean}");
        assert!((std - 2.0).abs() < 0.1, "std {std}");
        // UNIFORM: abs U(-3,3) → bounded.
        let pu = ErrorParams {
            abs_variance: 3.0,
            rel_variance: 0.0,
            ..ErrorParams::standard()
        };
        let mut r = rng();
        for _ in 0..1000 {
            let d = Group::apply_error(10.0, ErrorType::Uniform, &pu, &mut r);
            assert!((7.0..=13.0).contains(&d), "out of range {d}");
        }
        // COMBINED runs finite.
        let mut r = rng();
        for _ in 0..100 {
            let d = Group::apply_error(10.0, ErrorType::Combined, &ErrorParams::standard(), &mut r);
            assert!(d.is_finite());
        }
    }

    /// The uniform law is `U(mean ± width)`, not `U(±width)`.
    ///
    /// `abs_mean_delta_h` was declared from the first Python upload and read
    /// by no version until 0.7.15: a configured bias was accepted, stored and
    /// dropped. Both bounds are asserted, so a regression to the old
    /// zero-centred draw fails on the lower one rather than merely widening.
    #[test]
    fn the_uniform_law_is_centred_on_its_mean() {
        let p = ErrorParams {
            abs_mean_delta_h: 5.0,
            abs_variance: 1.0,
            rel_mean_delta_h: 0.0,
            rel_variance: 0.0,
            ..ErrorParams::standard()
        };
        let mut r = rng();
        let n = 20_000;
        let mut sum = 0.0;
        for _ in 0..n {
            let d = Group::apply_error(10.0, ErrorType::Uniform, &p, &mut r);
            // Old behaviour drew U(-1, 1) and landed in [9, 11].
            assert!((14.0..=16.0).contains(&d), "out of the shifted band: {d}");
            sum += d;
        }
        let mean = sum / n as f64;
        assert!((mean - 15.0).abs() < 0.05, "mean {mean}");
    }

    /// The relative channel carries its own centre, scaled by the value.
    #[test]
    fn the_relative_uniform_law_is_centred_too() {
        let p = ErrorParams {
            abs_mean_delta_h: 0.0,
            abs_variance: 0.0,
            rel_mean_delta_h: 0.1,
            rel_variance: 0.01,
            ..ErrorParams::standard()
        };
        let mut r = rng();
        for _ in 0..2_000 {
            let d = Group::apply_error(10.0, ErrorType::Uniform, &p, &mut r);
            // 10 + U(0.1 ± 0.01)*10 = 10 + [0.9, 1.1] = [10.9, 11.1].
            assert!((10.9..=11.1).contains(&d), "out of the shifted band: {d}");
        }
    }

    /// Zero width contributes the mean deterministically, as `gauss_draw`
    /// does at zero spread. Before 0.7.15 it contributed 0.0 and the bias
    /// vanished entirely.
    #[test]
    fn a_zero_width_uniform_contributes_its_mean() {
        let p = ErrorParams {
            abs_mean_delta_h: 3.0,
            abs_variance: 0.0,
            rel_mean_delta_h: 0.0,
            rel_variance: 0.0,
            ..ErrorParams::standard()
        };
        let mut r = rng();
        for _ in 0..16 {
            assert_eq!(
                Group::apply_error(10.0, ErrorType::Uniform, &p, &mut r),
                13.0
            );
        }
    }

    /// The fix is inert at the shipped defaults, which is why it can land
    /// without moving a single pinned number: both centres are 0.0, and
    /// `U(0 ± w)` is `U(±w)`.
    #[test]
    fn zero_centres_leave_the_draw_where_it_was() {
        let p = ErrorParams::standard();
        assert_eq!(p.abs_mean_delta_h, 0.0);
        assert_eq!(p.rel_mean_delta_h, 0.0);
        let mut a = rng();
        let mut b = rng();
        for _ in 0..64 {
            let got = Group::apply_error(10.0, ErrorType::Uniform, &p, &mut a);
            // Reproduce the pre-fix expression from the same stream. The
            // absolute Gaussian is drawn unconditionally, before the match,
            // so the mirror has to consume it too or the streams desync.
            let _ = gauss_draw(p.abs_mean_delta_g, p.abs_std_dev, &mut b);
            let want = 10.0
                + unif_draw(0.0, p.abs_variance, &mut b)
                + unif_draw(0.0, p.rel_variance, &mut b) * 10.0;
            assert_eq!(got, want);
        }
    }

    /// `Cascaded` is `Combined` plus the cross term `v * G_rel * U_rel`,
    /// drawing the same four numbers in the same order. Replaying the two
    /// laws on parallel streams pins that difference. The identity is exact
    /// in real arithmetic but not bitwise: the product associates the
    /// multiplications differently from the sum, so it holds to a rounding.
    #[test]
    fn cascaded_is_combined_plus_the_cross_term() {
        let p = ErrorParams::standard();
        let mut a = rng();
        let mut b = rng();
        for _ in 0..256 {
            let casc = Group::apply_error(10.0, ErrorType::Cascaded, &p, &mut a);
            // The same stream, decomposed: apply_error's order is the
            // unconditional absolute Gaussian, then g_rel, u_abs, u_rel.
            let g_abs = gauss_draw(p.abs_mean_delta_g, p.abs_std_dev, &mut b);
            let g_rel = gauss_draw(p.rel_mean_delta_g, p.rel_std_dev, &mut b);
            let u_abs = unif_draw(p.abs_mean_delta_h, p.abs_variance, &mut b);
            let u_rel = unif_draw(p.rel_mean_delta_h, p.rel_variance, &mut b);
            let combined = 10.0 + g_abs + g_rel * 10.0 + u_abs + u_rel * 10.0;
            let want = combined + 10.0 * g_rel * u_rel;
            assert!(
                (casc - want).abs() <= 1e-14 * want.abs().max(1.0),
                "{casc} vs {want}"
            );
        }
    }

    /// One relative channel off makes its factor exactly 1, so the product
    /// collapses and `Cascaded` becomes `Combined` -- to a rounding, since
    /// `v * (1 + g)` and `v + g * v` are the same number but not the same
    /// sequence of operations. This is what makes the new law a superset
    /// rather than a rival: turning a channel off cannot change the answer.
    #[test]
    fn cascaded_collapses_to_combined_when_a_factor_is_one() {
        for p in [
            // Gaussian relative off: (1 + 0) * (1 + U_rel).
            ErrorParams {
                rel_mean_delta_g: 0.0,
                rel_std_dev: 0.0,
                ..ErrorParams::standard()
            },
            // Uniform relative off: (1 + G_rel) * (1 + 0).
            ErrorParams {
                rel_mean_delta_h: 0.0,
                rel_variance: 0.0,
                ..ErrorParams::standard()
            },
        ] {
            let mut a = rng();
            let mut b = rng();
            for _ in 0..128 {
                let casc = Group::apply_error(10.0, ErrorType::Cascaded, &p, &mut a);
                let comb = Group::apply_error(10.0, ErrorType::Combined, &p, &mut b);
                assert!(
                    (casc - comb).abs() <= 1e-14 * comb.abs().max(1.0),
                    "{casc} vs {comb}"
                );
            }
        }
    }

    /// The composition is what it says: two multiplicative stages in series.
    /// Systematic centres with no spread make the draw deterministic, so the
    /// product is checkable against arithmetic rather than a statistic.
    #[test]
    fn cascaded_multiplies_its_systematic_centres() {
        let p = ErrorParams {
            abs_mean_delta_g: 0.0,
            abs_std_dev: 0.0,
            rel_mean_delta_g: 0.1,
            rel_std_dev: 0.0,
            abs_mean_delta_h: 0.0,
            abs_variance: 0.0,
            rel_mean_delta_h: 0.2,
            rel_variance: 0.0,
        };
        let mut r = rng();
        for _ in 0..16 {
            // 100 * 1.1 * 1.2 = 132, where Combined would give 100 * 1.3 = 130.
            let d = Group::apply_error(100.0, ErrorType::Cascaded, &p, &mut r);
            assert!((d - 132.0).abs() < 1e-12, "{d}");
        }
    }

    /// The expansion path builds the same law out of a scalar `(abs, rel)`
    /// pair, since `(1 + g)(1 + u) = 1 + (g + u + g*u)`. Both call sites must
    /// therefore agree on the value, whatever order they draw in.
    #[test]
    fn cascaded_agrees_across_both_draw_paths() {
        let p = ErrorParams::standard();
        let mut a = rng();
        let mut b = rng();
        for _ in 0..128 {
            // apply_error's order: g_abs, g_rel, u_abs, u_rel.
            let direct = Group::apply_error(10.0, ErrorType::Cascaded, &p, &mut a);
            let g_abs = gauss_draw(p.abs_mean_delta_g, p.abs_std_dev, &mut b);
            let g_rel = gauss_draw(p.rel_mean_delta_g, p.rel_std_dev, &mut b);
            let u_abs = unif_draw(p.abs_mean_delta_h, p.abs_variance, &mut b);
            let u_rel = unif_draw(p.rel_mean_delta_h, p.rel_variance, &mut b);
            // The expansion form: one scalar rel, applied as v + abs + rel*v.
            let rel = g_rel + u_rel + g_rel * u_rel;
            let via_pair = 10.0 + (g_abs + u_abs) + rel * 10.0;
            assert!((direct - via_pair).abs() <= 1e-12 * direct.abs().max(1.0));
        }
    }

    #[test]
    fn draws_are_deterministic_per_seed() {
        let p = ErrorParams::standard();
        let mut a = rng();
        let mut b = StdRng::seed_from_u64(7);
        for _ in 0..50 {
            let x = Group::apply_error(10.0, ErrorType::Combined, &p, &mut a);
            let y = Group::apply_error(10.0, ErrorType::Combined, &p, &mut b);
            assert_eq!(x, y);
        }
    }

    #[test]
    fn floors_match_python() {
        // Extreme negative systematic offset forces every floor.
        let mut g = Group::new("g");
        for params in [
            &mut g.thickness_error_params,
            &mut g.roughness_error_params,
            &mut g.interface_error_params,
            &mut g.n_error_params,
        ] {
            params.abs_mean_delta_g = -100.0;
            params.abs_std_dev = 0.0;
        }
        let mut r = rng();
        assert_eq!(g.thickness_error(50.0, &mut r), 0.0);
        assert_eq!(g.sr_roughness_error(1.0, &mut r), 0.0);
        assert_eq!(g.interface_error(2.0, &mut r), 0.0);
        let nk = g.nk_error(Complex64::new(2.0, 5.0), &mut r);
        assert_eq!(nk.re, 0.0);
        // k untouched by the floor (default k law: zero-mean, unit-scale —
        // finite and unclamped is the assertion).
        assert!(nk.im.is_finite());
        // inh_delta unfloored.
        g.inh_delta_error_params.abs_mean_delta_g = -100.0;
        g.inh_delta_error_params.abs_std_dev = 0.0;
        assert!(g.inh_delta_error(0.2, &mut r) < 0.0);
    }

    #[test]
    fn state_round_trip_and_fingerprint() {
        let g = Group::new("TiO2");
        let v = g.to_state();
        let mut keys: Vec<_> = v.as_object().unwrap().keys().cloned().collect();
        keys.sort();
        assert_eq!(
            keys,
            [
                "error_mask",
                "group_name",
                "inh_delta_error_params",
                "inh_delta_error_type",
                "inh_delta_summand",
                "interface_error_params",
                "interface_error_type",
                "interface_summand",
                "k_error_params",
                "k_error_type",
                "k_factor",
                "n_error_params",
                "n_error_type",
                "n_factor",
                "optimization_mask",
                "roughness_error_params",
                "roughness_error_type",
                "roughness_summand",
                "schema_version",
                "thick_factor",
                "thick_summand",
                "thickness_error_params",
                "thickness_error_type"
            ]
        );
        let back = Group::from_state(&v).unwrap();
        assert_eq!(back, g);
        // Unknown keys ignored; version enforced; group_name defaults.
        let mut with_bogus = v.clone();
        with_bogus["bogus"] = json!(1);
        assert_eq!(Group::from_state(&with_bogus).unwrap(), g);
        let mut nover = v.clone();
        nover.as_object_mut().unwrap().remove("schema_version");
        assert!(Group::from_state(&nover).is_err());
        let minimal = json!({"schema_version": 1});
        let d = Group::from_state(&minimal).unwrap();
        assert_eq!(d.group_name, "default");
        assert_eq!((d.thick_factor, d.n_factor, d.k_factor), (1.0, 1.0, 1.0));
        // Bad enum discriminant refused (stricter than Python — fail-closed).
        let mut bad_enum = v.clone();
        bad_enum["n_error_type"] = json!(9);
        assert!(Group::from_state(&bad_enum).is_err());
    }
}
