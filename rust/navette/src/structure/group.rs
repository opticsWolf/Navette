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
    /// The eight field names, in declaration order.
    pub const FIELDS: [&'static str; 8] = [
        "abs_mean_delta_g",
        "abs_std_dev",
        "rel_mean_delta_g",
        "rel_std_dev",
        "abs_mean_delta_h",
        "abs_variance",
        "rel_mean_delta_h",
        "rel_variance",
    ];

    /// The fields a JSON object does not carry, in declaration order.
    ///
    /// Both live doors -- `Group::set_properties` and the `set_error_params`
    /// binding -- deserialize a whole `ErrorParams`, so an incomplete block
    /// is refused. Until 0.7.26 they refused it as "malformed", which is the
    /// wrong word for a well-formed subset and sent the reader looking for a
    /// syntax error that was not there; serde names only the first missing
    /// field, so even a careful reader fixed them one round trip at a time.
    /// This names all of them at once.
    ///
    /// Empty for a non-object, where the block really is malformed and
    /// serde's own message is the better one.
    pub fn missing_fields(value: &Value) -> Vec<&'static str> {
        match value.as_object() {
            Some(map) => Self::FIELDS
                .iter()
                .copied()
                .filter(|f| !map.contains_key(*f))
                .collect(),
            None => Vec::new(),
        }
    }

    /// The base law params, and the reference all six per-channel
    /// constructors below derive from.
    ///
    /// Derive, literally: each spells out only the fields it changes and
    /// takes the rest with `..Self::standard()`. Three of them wrote all
    /// eight as literals until 0.7.24, which meant a change to the relative
    /// spreads here would have reached three channels and silently skipped
    /// the other three. The relative defaults are the shared part, so they
    /// belong in exactly one place.
    ///
    /// Every channel has its own constructor as of 0.7.20, because the
    /// *absolute* spreads carry the unit of the quantity they perturb and the
    /// six channels do not share one: thickness, roughness and interface
    /// width are lengths in nm; `n`, `k` and the grading amplitude are
    /// dimensionless. A single shared default necessarily had five of them
    /// wrong. Only the relative spreads are genuinely common, being unit-free
    /// fractions in every channel.
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
            abs_std_dev: 0.001,
            abs_variance: 0.001,
            ..Self::standard()
        }
    }

    /// Thickness-channel defaults: `abs_*` 0.5 nm.
    ///
    /// The absolute spread here is a physical thickness in nm, and 0.5 nm is
    /// a one-sigma figure that matches what deposition control actually
    /// achieves -- the process sits somewhere around 0.1 to 1 nm absolute,
    /// depending on monitoring. It read 0.01 nm from the first Python upload
    /// through 0.7.19, fifty times tighter than any real chamber, which made
    /// the default an optimistic answer rather than a neutral one.
    ///
    /// This is the only channel whose default was widened rather than
    /// narrowed in the 0.7.18-0.7.21 pass. The other five were mis-scaled
    /// because they inherited a length's spread for a quantity that is not a
    /// length; this one was in the right unit all along and simply carried
    /// an unrealistic value. It also means the change is not inert for
    /// anyone who had enabled the thickness error channel and relied on the
    /// default: tolerance spreads widen by 50x in the absolute term, which
    /// is the point.
    pub fn thickness() -> Self {
        Self {
            abs_std_dev: 0.5,
            abs_variance: 0.5,
            ..Self::standard()
        }
    }

    /// Interface-channel defaults: `abs_*` x0.1 relative to `standard()`.
    ///
    /// `interface_thickness` is a width in nm, the same physical quantity
    /// class as `roughness`, and the two are gated by the same rules in
    /// `Layer::validate`. They now share a scale as well: an interface is a
    /// sub-nanometre-to-few-nanometre feature, so it takes roughness's
    /// absolute spread rather than a full layer thickness's.
    pub fn interface() -> Self {
        Self {
            abs_std_dev: 0.001,
            abs_variance: 0.001,
            ..Self::standard()
        }
    }

    /// Grading-amplitude defaults: `abs_*` x0.1 relative to `standard()`.
    ///
    /// This channel does not perturb the authored `inh_delta`. Expansion
    /// computes `current_delta = (delta_layer + inh_delta_summand) * 0.5`,
    /// clamps it, and perturbs *that* -- the ramp half-amplitude, which runs
    /// around 0.05 to 0.1 for a typical authored delta of 0.1 to 0.2. It is
    /// dimensionless: writing `c` for that half-amplitude, the profile scales
    /// the complex index by `1 - c ..= 1 + c`. So `standard()`'s 0.01 was a
    /// 10-20% absolute scatter in nanometre units on a quantity that has no
    /// unit; 0.001 is a percent-level modulation, in line with `index()` on
    /// the index it modulates.
    ///
    /// The letter matters here. `Layer::validate` writes the authored
    /// `inh_delta` as `d` and its window as `1 - d/2 ..= 1 + d/2`; this is
    /// the same physical ramp seen after the halving, so `c = d/2` and the
    /// two windows are one window. Both were written `d` until 0.7.24, which
    /// put a factor of two between two identically-named quantities in one
    /// module.
    ///
    /// Unlike the length channels this one is signed at the point of the
    /// draw -- the cap is `clamp(-cap, cap)` and a negative amplitude simply
    /// reverses the ramp -- which is why `inh_delta_error` has no floor.
    pub fn inh_delta() -> Self {
        Self {
            abs_std_dev: 0.001,
            abs_variance: 0.001,
            ..Self::standard()
        }
    }

    /// Index-channel defaults: `abs_*` x0.1 relative to `standard()`.
    ///
    /// Same unit argument as `extinction()`, one order of magnitude milder.
    /// `standard()`'s absolute 0.01 is sized for a thickness in nanometres;
    /// against a refractive index of 1.5 to 2.4 it is a scatter of under 1%,
    /// which is survivable but still a nanometre spread wearing an index's
    /// clothes. 0.001 is a defensible index tolerance -- the fourth decimal
    /// is where dispersion data itself usually stops being trustworthy --
    /// and it stays sane for a low-contrast film near `n = 1`, where 0.01
    /// starts to matter and `nk_error`'s floor at 0 is the only backstop.
    pub fn index() -> Self {
        Self {
            abs_std_dev: 0.001,
            abs_variance: 0.001,
            ..Self::standard()
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
            abs_std_dev: 0.0001,
            abs_variance: 0.0001,
            ..Self::standard()
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
            thickness_error_params: ErrorParams::thickness(),
            inh_delta_error_params: ErrorParams::inh_delta(),
            roughness_error_params: ErrorParams::roughness(),
            interface_error_params: ErrorParams::interface(),
            n_error_params: ErrorParams::index(),
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

    /// The scalar field one `set_properties` key names.
    ///
    /// These three resolvers exist so the dispatch cannot drift: each spells
    /// out every name its caller's guard pattern admits, and the catch-all is
    /// `unreachable!` instead of a field. The previous form dispatched all but
    /// the last alternative explicitly and let a catch-all stand for the
    /// remainder, which is only correct while the guard pattern's final
    /// alternative and the catch-all's target agree by hand. In
    /// `*_error_type` they did not: `interface_error_type` fell through to
    /// `roughness_error_type`. Adding a channel, or reordering the
    /// alternatives, now fails to compile rather than silently writing a
    /// neighbouring field.
    fn scalar_slot(&mut self, key: &str) -> &mut f64 {
        match key {
            "thick_factor" => &mut self.thick_factor,
            "thick_summand" => &mut self.thick_summand,
            "n_factor" => &mut self.n_factor,
            "k_factor" => &mut self.k_factor,
            "inh_delta_summand" => &mut self.inh_delta_summand,
            "roughness_summand" => &mut self.roughness_summand,
            "interface_summand" => &mut self.interface_summand,
            other => unreachable!("scalar_slot on unguarded key '{other}'"),
        }
    }

    /// The error-law field one `set_properties` key names. See [`Self::scalar_slot`].
    fn error_type_slot(&mut self, key: &str) -> &mut ErrorType {
        match key {
            "thickness_error_type" => &mut self.thickness_error_type,
            "n_error_type" => &mut self.n_error_type,
            "k_error_type" => &mut self.k_error_type,
            "inh_delta_error_type" => &mut self.inh_delta_error_type,
            "roughness_error_type" => &mut self.roughness_error_type,
            "interface_error_type" => &mut self.interface_error_type,
            other => unreachable!("error_type_slot on unguarded key '{other}'"),
        }
    }

    /// The error-params block one `set_properties` key names. See [`Self::scalar_slot`].
    fn error_params_slot(&mut self, key: &str) -> &mut ErrorParams {
        match key {
            "thickness_error_params" => &mut self.thickness_error_params,
            "inh_delta_error_params" => &mut self.inh_delta_error_params,
            "roughness_error_params" => &mut self.roughness_error_params,
            "interface_error_params" => &mut self.interface_error_params,
            "n_error_params" => &mut self.n_error_params,
            "k_error_params" => &mut self.k_error_params,
            other => unreachable!("error_params_slot on unguarded key '{other}'"),
        }
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
                s @ ("thick_factor" | "thick_summand" | "n_factor" | "k_factor"
                | "inh_delta_summand" | "roughness_summand" | "interface_summand") => {
                    match value.as_f64() {
                        Some(v) => *self.scalar_slot(s) = v,
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
                        Some(Ok(e)) => *self.error_type_slot(t) = e,
                        _ => warnings.push(bad(format!("unknown {t}; ignoring."))),
                    }
                }
                p @ ("thickness_error_params"
                | "inh_delta_error_params"
                | "roughness_error_params"
                | "interface_error_params"
                | "n_error_params"
                | "k_error_params") => match serde_json::from_value::<ErrorParams>(value.clone()) {
                    Ok(ep) => *self.error_params_slot(p) = ep,
                    Err(e) => {
                        let missing = ErrorParams::missing_fields(value);
                        // Only absence gets the incompleteness message. A bad
                        // *value* on a field that IS present keeps serde's own
                        // diagnosis, or naming the absent fields would bury the
                        // one thing actually wrong with what was written.
                        warnings.push(bad(
                            if missing.is_empty() || !e.to_string().starts_with("missing field") {
                                format!("ignoring malformed '{p}': {e}.")
                            } else {
                                format!(
                                    "ignoring incomplete '{p}': missing {}. A params block \
                                 replaces the channel whole, so all eight fields are \
                                 required; to change a few, edit a copy of the current \
                                 block rather than naming only what moves.",
                                    missing.join(", ")
                                )
                            },
                        ));
                    }
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

    /// Every `Group::new` default, pinned.
    ///
    /// It began as a Python-parity oracle and was still called
    /// `defaults_match_python_ctor` at 0.7.23, by which point its own doc
    /// comment had to concede that the error params deliberately no longer
    /// match Python's. The factors, summands, masks and law choice still do,
    /// and those assertions are unchanged; the error params are pinned
    /// per channel instead, so a change has to be deliberate rather than
    /// merely passing.
    #[test]
    fn constructor_defaults_are_pinned() {
        let g = Group::new("TiO2");
        assert_eq!((g.thick_factor, g.thick_summand), (1.0, 0.0));
        assert_eq!((g.n_factor, g.k_factor), (1.0, 1.0));
        assert_eq!(g.nk_factor(), Complex64::new(1.0, 1.0));
        assert_eq!(g.error_mask, [0; 6]);
        assert_eq!(g.optimization_mask, [1; 7]);
        assert_eq!(g.thickness_error_type, ErrorType::Gaussian);

        // The error-param defaults are no longer Python's: the absolute
        // spreads carry the unit of the quantity they perturb, and one shared
        // value was wrong for five of the six channels (0.7.18-0.7.21).
        // Pinned here per channel so a change has to be deliberate.
        for (channel, abs, want) in [
            ("thickness", g.thickness_error_params.abs_std_dev, 0.5),
            ("roughness", g.roughness_error_params.abs_std_dev, 0.001),
            ("interface", g.interface_error_params.abs_std_dev, 0.001),
            ("inh_delta", g.inh_delta_error_params.abs_std_dev, 0.001),
            ("n", g.n_error_params.abs_std_dev, 0.001),
            ("k", g.k_error_params.abs_std_dev, 0.0001),
        ] {
            assert_eq!(abs, want, "{channel} abs_std_dev");
        }
        // The uniform half-width takes the same number as the Gaussian sigma
        // in every channel, so the two laws stay within a small factor of one
        // another and switching `ErrorType` cannot change the scatter by an
        // order of magnitude.
        //
        // It does not make them equal, and an earlier version of this comment
        // claimed it did. A uniform of half-width `w` has standard deviation
        // `w/sqrt(3)`, so matching the numbers leaves the uniform law 42%
        // NARROWER in sigma, on bounded rather than unbounded support; and
        // `Combined`/`Cascaded` add the two in quadrature, so they are wider
        // than either. Measured on a 100 nm layer at the thickness defaults
        // (abs 0.5 nm, rel 0.01):
        //
        //   Gaussian sigma 1.105   Uniform sigma 0.647  (-41%)
        //   Combined sigma 1.283   Cascaded sigma 1.283  (+16% vs Gaussian)
        //
        // Equal *numbers* is the property worth pinning -- it is what keeps a
        // channel's two spreads from drifting apart by a factor of ten -- so
        // the assertion stands; only the reason given for it was wrong.
        for (channel, params) in [
            ("thickness", &g.thickness_error_params),
            ("roughness", &g.roughness_error_params),
            ("interface", &g.interface_error_params),
            ("inh_delta", &g.inh_delta_error_params),
            ("n", &g.n_error_params),
            ("k", &g.k_error_params),
        ] {
            assert_eq!(
                params.abs_variance, params.abs_std_dev,
                "{channel} abs_variance vs abs_std_dev"
            );
            // The relative spreads are unit-free fractions, common to all.
            assert_eq!(params.rel_std_dev, 0.01, "{channel} rel_std_dev");
            assert_eq!(params.rel_variance, 0.01, "{channel} rel_variance");
        }
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

    /// The factored law reduces to one scalar `(abs, rel)` pair, since
    /// `(1 + g)(1 + u) = 1 + (g + u + g*u)`. That is what lets
    /// `expansion::channel_draws` carry `Cascaded` without changing its
    /// signature or its call sites.
    ///
    /// This checks the *algebra*, from `apply_error`'s draw order on both
    /// sides. It was called `cascaded_agrees_across_both_draw_paths` until
    /// 0.7.24, which overstated it twice over: it never calls
    /// `channel_draws`, and the two paths do not in fact produce the same
    /// number from the same stream state. They draw in deliberately
    /// different orders -- `channel_draws` takes the abs pair then the rel
    /// pair, `apply_error` takes `g_abs, g_rel, u_abs, u_rel` -- so the same
    /// seed feeds the four draws to different slots. Nothing is broken by
    /// that: the two serve different channels (`channel_draws` does n and k
    /// at expansion, `apply_error` the other four and the `nk_error` probe)
    /// on independent per-side streams, and it is equally true of
    /// `Combined`, which predates this law. But no test should be read as
    /// promising they agree numerically, because they do not.
    #[test]
    fn cascaded_reduces_to_one_scalar_rel_pair() {
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

    /// Every `*_error_type` key lands in its own channel and moves no other.
    ///
    /// `set_properties` dispatched five of the six names explicitly and let a
    /// catch-all take the sixth, which only works while the guard pattern's
    /// last alternative and the catch-all's target agree. They did not:
    /// `interface_error_type` fell through to `roughness_error_type`, so
    /// asking for interface silently wrote roughness and left interface
    /// untouched -- no error, no warning, the wrong field mutated.
    #[test]
    fn set_properties_lands_each_error_type_in_its_own_channel() {
        let read: [(&str, fn(&Group) -> ErrorType); 6] = [
            ("thickness_error_type", |g| g.thickness_error_type),
            ("n_error_type", |g| g.n_error_type),
            ("k_error_type", |g| g.k_error_type),
            ("inh_delta_error_type", |g| g.inh_delta_error_type),
            ("roughness_error_type", |g| g.roughness_error_type),
            ("interface_error_type", |g| g.interface_error_type),
        ];
        for (key, _) in read {
            let mut g = Group::new("g");
            let mut props = BTreeMap::new();
            props.insert(key.to_string(), Value::from(ErrorType::Cascaded.as_i32()));
            let warnings = g.set_properties(&props);
            assert!(warnings.is_empty(), "{key}: {warnings:?}");
            for (other, get) in read {
                let want = if other == key {
                    ErrorType::Cascaded
                } else {
                    ErrorType::Gaussian
                };
                assert_eq!(get(&g), want, "setting {key} moved {other}");
            }
        }
    }

    /// An incomplete params block is named as incomplete, with every missing
    /// field listed, not as "malformed" with one of them.
    #[test]
    fn an_incomplete_params_block_names_every_missing_field() {
        let mut g = Group::new("g");
        let mut props = BTreeMap::new();
        props.insert(
            "k_error_params".to_string(),
            serde_json::json!({"abs_std_dev": 0.5, "rel_std_dev": 0.02}),
        );
        let warnings = g.set_properties(&props);
        assert_eq!(warnings.len(), 1);
        let msg = warnings[0].message.clone();
        assert!(msg.contains("incomplete"), "{msg}");
        assert!(!msg.contains("malformed"), "{msg}");
        // The six absent fields, all of them, and neither of the two present.
        for f in [
            "abs_mean_delta_g",
            "rel_mean_delta_g",
            "abs_mean_delta_h",
            "abs_variance",
            "rel_mean_delta_h",
            "rel_variance",
        ] {
            assert!(msg.contains(f), "missing field {f} not named: {msg}");
        }
        // Refusal is unchanged: the block is still not applied.
        assert_eq!(g.k_error_params, ErrorParams::extinction());

        // A block that is genuinely malformed keeps serde's own diagnosis.
        let mut props = BTreeMap::new();
        props.insert("k_error_params".to_string(), serde_json::json!("nonsense"));
        let warnings = g.set_properties(&props);
        assert!(warnings[0].message.contains("malformed"), "{warnings:?}");

        // A bad VALUE on a present field keeps serde's diagnosis too, even
        // though the block is also incomplete -- otherwise listing the absent
        // fields would bury the one thing actually wrong with what was typed.
        let mut props = BTreeMap::new();
        props.insert(
            "k_error_params".to_string(),
            serde_json::json!({"abs_std_dev": "not a number"}),
        );
        let warnings = g.set_properties(&props);
        let msg = warnings[0].message.clone();
        assert!(msg.contains("malformed"), "{msg}");
        assert!(msg.contains("invalid type"), "{msg}");
        assert!(!msg.contains("incomplete"), "{msg}");
    }

    #[test]
    fn missing_fields_lists_declaration_order_and_nothing_for_a_complete_block() {
        let complete = serde_json::to_value(ErrorParams::standard()).unwrap();
        assert!(ErrorParams::missing_fields(&complete).is_empty());
        assert_eq!(
            ErrorParams::missing_fields(&serde_json::json!({})),
            ErrorParams::FIELDS.to_vec()
        );
        // Not an object: serde's message is the better one, so nothing here.
        assert!(ErrorParams::missing_fields(&serde_json::json!(7)).is_empty());
    }

    /// The same trap in the sibling arms of `set_properties`: a catch-all
    /// standing for the guard pattern's last alternative. Both were correct,
    /// and both were one reordering away from the defect above.
    #[test]
    fn set_properties_lands_each_summand_and_params_block_in_its_own_field() {
        let mut g = Group::new("g");
        let mut props = BTreeMap::new();
        props.insert("interface_summand".to_string(), Value::from(3.5));
        assert!(g.set_properties(&props).is_empty());
        assert_eq!(g.interface_summand, 3.5);
        assert_eq!(g.roughness_summand, 0.0);
        assert_eq!(g.inh_delta_summand, 0.0);

        let mut g = Group::new("g");
        let mut props = BTreeMap::new();
        let ep = ErrorParams {
            abs_std_dev: 7.0,
            ..ErrorParams::standard()
        };
        props.insert(
            "k_error_params".to_string(),
            serde_json::to_value(&ep).unwrap(),
        );
        assert!(g.set_properties(&props).is_empty());
        assert_eq!(g.k_error_params.abs_std_dev, 7.0);
        assert_eq!(g.n_error_params.abs_std_dev, 0.001);
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
