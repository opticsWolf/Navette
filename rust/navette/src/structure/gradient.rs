// SPDX-License-Identifier: LGPL-3.0-or-later
//! Gradient mixture profiles (F1.1): spec types, the sublayer-resolution
//! rule, validation, and the EMA application.
//!
//! A `GradientSpec` turns one design layer into a span of sublayers whose
//! complex indices interpolate between two provider-resolved endpoint
//! spectra, `material_a` (host, `f = 0`) and `material_b` (inclusion,
//! `f = 1`), through one of the [`MixRule`] kernels.
//!
//! HOST/INCLUSION ASYMMETRY (§D9.4, load-bearing). `material_a` is ALWAYS
//! the host and `f` is the volume fraction of B in A: `f = 0.1` means 10 %
//! B inclusions in 90 % A. Maxwell-Garnett and Mori-Tanaka are
//! host/inclusion-asymmetric by construction, so swapping `material_a`
//! and `material_b` is NOT the same physics, and a silent swap corrupts
//! fits. The kernels here are therefore called as
//! `ema(nk_b, nk_a, f)` - inclusion first, mirroring `MaterialSpec`'s
//! inclusion/host argument order - so the vocabulary users already know
//! from `EffectiveMaterial` (host / inclusion / fraction) transfers
//! directly.
//!
//! F1.1 ships `GradientMode::FixedSpan` and [`ProfileShape::Linear`]
//! only; `RateCapped` (F1.2/F1.3) and further shapes arrive as new enum
//! variants with no schema break.

use num_complex::Complex64;
use serde::{Deserialize, Serialize};

use crate::materials::MixRule;
use crate::materials::ema;
use crate::structure::validation::ValidationIssue;

/// Intra-span profile shape. V1 ships Linear only (§D9.2, locked); curved
/// profiles (S-curve, exponential, ...) arrive as new variants - no schema
/// break, `mode` and `shape` are already enums, and each new variant gets
/// its own twin batch.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ProfileShape {
    /// Linear in physical depth: `f(z)` interpolates uniformly between the
    /// mode's endpoints.
    Linear,
}

/// Profile kind of a gradient span. F1.1 ships the thickness-independent
/// variant; `RateCapped` (slope in absolute depth, clamped) is F1.2.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum GradientMode {
    /// (a) `f_start` at the deposition face, `f_end` at the far face, over
    /// any thickness: `f(z) = f_start + (f_end - f_start) * z / thickness`.
    /// Endpoints depend only on `f_start`/`f_end`, never on thickness
    /// (only the sublayer count changes) - which is what makes the mode
    /// scale freely at F1.6.
    FixedSpan { f_start: f64, f_end: f64 },
    /// (b) `f(z) = f_start + rate * (z / ref_thickness)`, clamped to
    /// `[f_min, f_max]`. The slope is thickness-RELATIVE (e.g. 0.3 per
    /// 100 nm); a thick film saturates into a flat pure-material tail
    /// where the cap binds - which is both the physics and the merge
    /// test (N2: the tail rows are nk-identical, so `merge_adjacent`
    /// collapses them). NOT scale-free (F1.6's table: the fraction
    /// depends on absolute `z`) - F1.7's refresh territory.
    RateCapped {
        f_start: f64,
        rate: f64,
        ref_thickness: f64,
        f_min: f64,
        f_max: f64,
    },
}

/// Mixture gradient specification for one layer.
///
/// The docstring obligation (§D9.4) lives on the module: `material_a` is
/// always the HOST, `f` is the volume fraction of B in A, and the
/// asymmetric kernels (Maxwell-Garnett, Mori-Tanaka) change physics
/// under an A/B swap.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GradientSpec {
    /// Host material (provider key); the `f = 0` endpoint. Names the same
    /// provider the layer's own material resolves through.
    pub material_a: String,
    /// Inclusion material (provider key); the `f = 1` endpoint.
    pub material_b: String,
    /// Mixing kernel (A6: the whole `MixRule` enum rides along, all six
    /// variants work on day one, parameters live in the variants).
    /// Default Bruggeman (§D9.1, locked).
    pub ema: MixRule,
    /// Profile mode (F1.1: `FixedSpan` only).
    pub mode: GradientMode,
    /// Profile shape (F1.1: `Linear` only).
    pub shape: ProfileShape,
    /// Sublayer-count override. `Some(n)` is clamped to `[2, 256]` at
    /// expansion; outside that window it is an advisory finding, not an
    /// error (see [`Self::issues`]).
    pub sublayers: Option<u32>,
}

impl GradientSpec {
    /// The FixedSpan spec shorthand used by tests and callers that mean
    /// the default kernel and shape.
    pub fn fixed_span(
        material_a: impl Into<String>,
        material_b: impl Into<String>,
        f_start: f64,
        f_end: f64,
    ) -> Self {
        Self {
            material_a: material_a.into(),
            material_b: material_b.into(),
            ema: MixRule::Bruggeman {
                max_iter: 100,
                tol: 1e-9,
            },
            mode: GradientMode::FixedSpan { f_start, f_end },
            shape: ProfileShape::Linear,
            sublayers: None,
        }
    }

    /// The RateCapped constructor shorthand (F1.2): slope per
    /// `ref_thickness` with the cap window; defaults per D0(b).
    pub fn rate_capped(
        material_a: impl Into<String>,
        material_b: impl Into<String>,
        f_start: f64,
        rate: f64,
        ref_thickness: f64,
        f_min: f64,
        f_max: f64,
    ) -> Self {
        Self {
            material_a: material_a.into(),
            material_b: material_b.into(),
            ema: MixRule::Bruggeman {
                max_iter: 100,
                tol: 1e-9,
            },
            mode: GradientMode::RateCapped {
                f_start,
                rate,
                ref_thickness,
                f_min,
                f_max,
            },
            shape: ProfileShape::Linear,
            sublayers: None,
        }
    }

    /// The endpoints of the span, in traversal order (deposition face
    /// first). A `match` so a new variant cannot be forgotten here.
    pub fn endpoints(&self) -> (f64, f64) {
        match self.mode {
            GradientMode::FixedSpan { f_start, f_end } => (f_start, f_end),
            // RateCapped's "endpoints" for the advisory surface: the
            // nominal start and the upper cap (f_min unused here).
            GradientMode::RateCapped { f_start, f_max, .. } => (f_start, f_max),
        }
    }

    /// The profile value at ABSOLUTE depth `z_nm` of a layer
    /// `thickness_nm` thick (deposition face at 0), clamped for
    /// `RateCapped`. One `match` per mode.
    pub fn f_at(&self, z_nm: f64, thickness_nm: f64) -> f64 {
        match self.mode {
            GradientMode::FixedSpan { f_start, f_end } => {
                let x = if thickness_nm > 0.0 {
                    z_nm / thickness_nm
                } else {
                    0.0
                };
                f_start + (f_end - f_start) * x
            }
            GradientMode::RateCapped {
                f_start,
                rate,
                ref_thickness,
                f_min,
                f_max,
            } => (f_start + rate * (z_nm / ref_thickness)).clamp(f_min, f_max),
        }
    }

    /// Mid-profile value used by the design path's homogenize (§D4.3):
    /// the mean of the two endpoints for `FixedSpan`; the mean of the
    /// CLAMPED profile values at both faces for `RateCapped`.
    pub fn f_mid(&self, thickness_nm: f64) -> f64 {
        match self.mode {
            GradientMode::FixedSpan { f_start, f_end } => (f_start + f_end) * 0.5,
            GradientMode::RateCapped { .. } => {
                (self.f_at(0.0, thickness_nm) + self.f_at(thickness_nm, thickness_nm)) * 0.5
            }
        }
    }

    /// The profile value of sublayer `i` of `n` - INCLUSIVE endpoints:
    /// the first sublayer IS `f_start`, the last IS `f_end`. This is the
    /// legacy graded branch's own convention (`factors[i]` runs
    /// `i / (sub-1)`), and it is what makes the endpoint rows bitwise
    /// thickness-independent (the F1.1 gate) and the mode scale freely
    /// at fixed row count (F1.6's table: "normalized depth between two
    /// endpoints"). §D3's "midpoint z_i" reading is deliberately NOT
    /// implemented: with midpoint sampling no row sits at an endpoint,
    /// and the gate's bitwise claim would be unimplementable.
    pub fn f_sublayer(&self, i: u32, n: u32, thickness_nm: f64) -> f64 {
        match self.mode {
            GradientMode::FixedSpan { f_start, f_end } => {
                if n <= 1 || i == 0 {
                    // f_start + 0.0 is bitwise f_start anyway; stating it
                    // makes the contract explicit.
                    f_start
                } else if i + 1 == n {
                    // The ratio form would give f_start + (f_end -
                    // f_start) * 1.0, which is f_end only to one ulp.
                    // The last sublayer IS the f_end endpoint.
                    f_end
                } else {
                    f_start + (f_end - f_start) * f64::from(i) / f64::from(n - 1)
                }
            }
            // RateCapped samples the SAME inclusive depth grid (the row
            // itself is the profile value at z_i = t * i / (n - 1), the
            // clamp binding where it binds - saturated tail rows come out
            // bitwise pure material). One subtlety: unlike FixedSpan,
            // f(0) is clamp(f_start), which is NOT f_start when f_start
            // sits outside the caps - the caps govern, per the mode
            // formula.
            GradientMode::RateCapped { .. } => {
                let z_i = if n <= 1 {
                    0.0
                } else {
                    thickness_nm * f64::from(i) / f64::from(n - 1)
                };
                self.f_at(z_i, thickness_nm)
            }
        }
    }

    /// The fatal profile checks, as one message (the first hit). Run at
    /// expansion time (D2: validation "at parse + at expand") so a spec
    /// assembled around the gates still cannot emit. `inhomogen` is the
    /// carrying layer's flag - the two profile engines are mutually
    /// exclusive (§D2) and this is the last line of that defence.
    pub fn expansion_error(&self, inhomogen: bool) -> Option<String> {
        if inhomogen {
            return Some(
                "gradient and inhomogen are both set: two profile engines on \
                 one layer. Single-material drift = inhomogen alone; a mixture \
                 profile = gradient alone - set one."
                    .to_string(),
            );
        }
        if self.material_a == self.material_b {
            return Some(format!(
                "gradient materials are identical ('{}'): a gradient of a \
                 material with itself is a spec bug; single-material drift \
                 is inhomogen.",
                self.material_a
            ));
        }
        match self.mode {
            GradientMode::FixedSpan { f_start, f_end } => {
                for (name, v) in [("f_start", f_start), ("f_end", f_end)] {
                    if !(0.0..=1.0).contains(&v) {
                        return Some(format!("gradient {name} {v} is outside [0, 1]."));
                    }
                }
                if f_start == f_end {
                    return Some(format!(
                        "gradient f_start == f_end ({f_start}): a constant profile - \
                         the span is degenerate and merges away; use a plain layer."
                    ));
                }
            }
            GradientMode::RateCapped {
                f_start,
                rate,
                ref_thickness,
                f_min,
                f_max,
            } => {
                if !rate.is_finite() {
                    return Some(format!("gradient rate {rate} is not finite."));
                }
                if !(ref_thickness > 0.0) {
                    return Some(format!(
                        "gradient ref_thickness {ref_thickness} must be > 0."
                    ));
                }
                if !(0.0..=1.0).contains(&f_start) {
                    return Some(format!("gradient f_start {f_start} is outside [0, 1]."));
                }
                for (name, v) in [("f_min", f_min), ("f_max", f_max)] {
                    if !(0.0..=1.0).contains(&v) {
                        return Some(format!("gradient {name} {v} is outside [0, 1]."));
                    }
                }
                if f_min > f_max {
                    return Some(format!(
                        "gradient f_min {f_min} > f_max {f_max}: the cap window \
                         is empty."
                    ));
                }
                if rate == 0.0 {
                    return Some(
                        "gradient rate 0: a constant profile - the span is \
                         degenerate and merges away; use a plain layer."
                            .to_string(),
                    );
                }
            }
        }
        None
    }

    /// The self-contained findings (N5: the one rule surface for the
    /// numbers). Unprefixed; callers prefix with their own label. The
    /// provider-existence check is deliberately NOT here - it needs the
    /// provider and lives at the two provider doors, sharing
    /// [`endpoint_error`] for its wording.
    pub fn issues(&self) -> Vec<ValidationIssue> {
        let mut issues = Vec::new();
        let bad = |m: String| ValidationIssue::error(m);
        let note = |m: String| ValidationIssue::warning(m);
        if let Some(msg) = self.expansion_error(false) {
            issues.push(bad(msg));
        }
        if let Some(n) = self.sublayers
            && !(2..=256).contains(&n)
        {
            issues.push(note(format!(
                "gradient sublayers {n} outside [2, 256]; clamped at expansion."
            )));
        }
        issues
    }
}

/// The provider-door refusal helper (N5): ONE wording, used by the
/// structure path's `expand` and the design path's `from_design`, naming
/// the layer, the span's endpoints, and the provider's own complaint
/// (which names the absent material).
pub(crate) fn endpoint_error(layer_label: &str, grad: &GradientSpec, err: String) -> String {
    format!(
        "layer '{}': gradient {} -> {}: {err}",
        layer_label, grad.material_a, grad.material_b
    )
}

/// Default sublayer step for gradient spans [nm]: never coarser than a
/// tenth of the shortest optical wavelength in the coarsest material on
/// the grid, capped so thick films cannot explode the row count.
pub const MAX_GRADIENT_STEP_NM: f64 = 20.0;

/// Sublayer count of a gradient span (D3 step 2, formula fixed at F1.1):
///
/// ```text
/// n_sub     = clamp(ceil(thickness_nm / max_step_nm), 3, 64)
/// max_step  = min(20.0, lambda_min / (10 * n_max_re))
/// ```
///
/// The step must be small against the local optical wavelength or the
/// staircase itself becomes the physics; `lambda_min / 10n` is the
/// conventional floor and the 20 nm cap keeps thick films bounded.
/// `sublayers: Some(n)` overrides, clamped to `[2, 256]` (the advisory
/// for an out-of-range override is [`GradientSpec::issues`]'s job).
///
/// Pinned by a differential test over randomized thicknesses and grids
/// (`gradient_count_matches_formula`): integer boundaries must agree
/// bit-exactly.
pub fn gradient_sub_layer_count(
    thickness_nm: f64,
    wavelengths: &[f64],
    nk_a: &[Complex64],
    nk_b: &[Complex64],
    sublayers: Option<u32>,
) -> u32 {
    if let Some(n) = sublayers {
        return n.clamp(2, 256);
    }
    let lambda_min = wavelengths.iter().copied().fold(f64::INFINITY, f64::min);
    let n_max_re = nk_a
        .iter()
        .chain(nk_b.iter())
        .map(|z| z.re)
        .fold(0.0f64, f64::max);
    let step = if n_max_re > 0.0 && lambda_min.is_finite() {
        MAX_GRADIENT_STEP_NM.min(lambda_min / (10.0 * n_max_re))
    } else {
        MAX_GRADIENT_STEP_NM
    };
    let raw = (thickness_nm / step).ceil();
    // Non-finite or negative leftovers clamp into the window; the doors
    // refuse non-finite thickness before this runs, and zero-thickness
    // layers never reach the gradient branch (expansion falls back to
    // the uniform row for them).
    raw.clamp(3.0, 64.0) as u32
}

/// One sublayer's complex index: `EMA(nk_b, nk_a, f)` over the whole
/// grid, inclusion B first (§D9.4's argument order, mirrored from
/// `MaterialSpec`), then the sqrt the Python `EffectiveMaterial` applies.
/// The F1.1 oracle twins call THIS against the expansion branch - same
/// kernel, same inputs, bitwise.
pub fn mix_row(rule: MixRule, nk_b: &[Complex64], nk_a: &[Complex64], f: f64) -> Vec<Complex64> {
    use ndarray::Array1;
    let b = Array1::from_vec(nk_b.to_vec());
    let a = Array1::from_vec(nk_a.to_vec());
    let eps = match rule {
        MixRule::Bruggeman { max_iter, tol } => {
            ema::bruggeman(b.view(), a.view(), f, max_iter, tol)
        }
        MixRule::MaxwellGarnett => ema::maxwell_garnett(b.view(), a.view(), f),
        MixRule::Looyenga => ema::looyenga(b.view(), a.view(), f),
        MixRule::Lichtenecker => ema::lichtenecker(b.view(), a.view(), f),
        MixRule::MoriTanaka { l } => ema::mori_tanaka(b.view(), a.view(), f, l),
        MixRule::PowerLaw { alpha } => ema::general_power_law(b.view(), a.view(), f, alpha),
    };
    ema::eps_to_nk(eps.view()).to_vec()
}

/// A5 (§8.8): the design path's gradient transport.
///
/// The synthesis design path has no materials library - every nk is
/// film-supplied - so a gradient design film's inclusion spectrum rides
/// the film dict itself (`nk_b`, evaluated Python-side in `_film_dicts`),
/// and `assemble_stack` injects `material_b -> nk_b` into the provider
/// entries. The host is expected to name the film itself (whose nk is
/// already registered); any other `material_a` fails the provider-door
/// refusal rather than falling back to a half-mixture.
#[derive(Clone, Debug, PartialEq)]
pub struct GradientJson {
    /// Host provider key; expected to name the film itself (A5).
    pub material_a: String,
    /// Inclusion provider key.
    pub material_b: String,
    /// The inclusion spectrum on the simulation grid (Python-evaluated).
    pub nk_b: Vec<Complex64>,
    /// Mixing kernel.
    pub ema: MixRule,
    /// The profile mode (F1.2: FixedSpan or RateCapped, tagged).
    pub mode: GradientMode,
    /// Sublayer-count override (clamped at expansion).
    pub sublayers: Option<u32>,
}

impl GradientJson {
    /// The spec the expansion branch consumes.
    pub fn to_spec(&self) -> GradientSpec {
        GradientSpec {
            material_a: self.material_a.clone(),
            material_b: self.material_b.clone(),
            ema: self.ema,
            mode: self.mode,
            shape: ProfileShape::Linear,
            sublayers: self.sublayers,
        }
    }
}

/// Application mode for the legacy scaling drift (physics frozen).
///
/// `Fixed` is today's behavior, byte for byte: the grading strength is
/// the authored `inh_delta`. `RateCapped` (F1.3) makes the strength grow
/// with thickness - D2's `delta_layer = min(rate * thickness /
/// ref_thickness, cap)` - a deposition drift that scales with the film
/// instead of being a fixed property. The clamp that binds the
/// COMBINED nominal is two-sided ([-cap, cap]): a negative rate
/// inverts the drift direction (the frozen ramp `1 - delta ..= 1 +
/// delta` runs inverted), and the magnitude still saturates at `cap`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum InhMode {
    /// Current behavior, bit-identical: `delta = inh_delta`.
    #[default]
    Fixed,
    /// `delta_layer = min(rate * thickness / ref_thickness, cap)`
    /// (D2's formula: the positive side saturates at `cap` in the
    /// delta itself; the nominal clamp below bounds both sides).
    RateCapped {
        rate: f64,
        ref_thickness: f64,
        cap: f64,
    },
}

impl InhMode {
    /// Whether the mode needs the nominal clamp (and its bound); the
    /// delta itself is [`Layer::delta_layer`] (it reads the layer's
    /// thickness and, for `Fixed`, the authored `inh_delta`).
    pub fn nominal_cap(&self) -> Option<f64> {
        match *self {
            InhMode::Fixed => None,
            InhMode::RateCapped { cap, .. } => Some(cap),
        }
    }

    /// The self-contained findings (the rate/ref/cap numbers). The
    /// provider-independent half of the validation.
    pub fn issues(&self) -> Vec<ValidationIssue> {
        let mut issues = Vec::new();
        let bad = |m: String| ValidationIssue::error(m);
        if let InhMode::RateCapped {
            rate,
            ref_thickness,
            cap,
        } = self
        {
            if !rate.is_finite() {
                issues.push(bad(format!("inh rate {rate} is not finite.")));
            }
            if !(*ref_thickness > 0.0) {
                issues.push(bad(format!(
                    "inh ref_thickness {ref_thickness} must be > 0."
                )));
            }
            if !(0.0..=1.0).contains(cap) {
                issues.push(bad(format!("inh cap {cap} is outside [0, 1].")));
            }
        }
        issues
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f_at_and_f_mid_are_pure_functions_of_the_fraction() {
        let g = GradientSpec::fixed_span("A", "B", 0.2, 0.8);
        assert_eq!(g.f_at(0.0, 100.0), 0.2);
        assert_eq!(g.f_at(100.0, 100.0), 0.8);
        assert_eq!(g.f_at(50.0, 100.0), 0.5);
        assert_eq!(g.f_mid(100.0), 0.5);
        // Bitwise: f_at(midpoint of the endpoints) == f_mid.
        assert_eq!(g.f_at(50.0, 100.0), g.f_mid(100.0));
        let (a, b) = g.endpoints();
        assert_eq!((a, b), (0.2, 0.8));
    }

    /// F1.2 - RateCapped: the D0(b) formula, the clamp binding where it
    /// binds, and the saturation gate (the tail rows are bitwise pure
    /// material). 200 nm film, f_start 0.1, rate 0.3 per 100 nm: the far
    /// face sits at 0.7; the same spec at 400 nm saturates to 1.0 from
    /// z = 300 nm.
    #[test]
    fn rate_capped_formula_and_saturation() {
        let g = GradientSpec::rate_capped("A", "B", 0.1, 0.3, 100.0, 0.0, 1.0);
        // The 200 nm readout: f(0) = 0.1, f(200) = 0.1 + 0.3 * 2 = 0.7.
        assert_eq!(g.f_at(0.0, 200.0), 0.1);
        assert_eq!(g.f_at(200.0, 200.0), 0.1 + 0.3 * 2.0);
        assert_eq!(g.f_at(100.0, 200.0), 0.1 + 0.3);
        assert_eq!(g.f_mid(200.0), (0.1 + (0.1 + 0.3 * 2.0)) * 0.5);
        // The 400 nm film: the clamp binds where the raw value exceeds
        // the cap. NOTE the knee arithmetic: 0.1 + 0.3 * 3.0 is
        // 0.9999999999999999 in IEEE (one ulp short of 1.0), so at
        // z = 300 the profile value is the RAW value, not the cap - the
        // plan's "saturates at 1.0 from z = 300" is idealized. The tail
        // (raw clearly >= cap) is bitwise the cap, and those are the
        // rows the saturation gate calls pure material_b.
        assert_eq!(g.f_at(300.0, 400.0), 0.1 + 0.3 * 3.0);
        assert_eq!(g.f_at(400.0, 400.0), 1.0, "raw 1.1 clamps to exactly 1.0");
        assert_eq!(g.f_at(1000.0, 400.0), 1.0, "deep tail is bitwise the cap");
        // Sub-caps bind too: a negative rate saturates to f_min.
        let gneg = GradientSpec::rate_capped("A", "B", 0.9, -0.3, 100.0, 0.2, 1.0);
        assert_eq!(gneg.f_at(0.0, 400.0), 0.9);
        assert_eq!(gneg.f_at(300.0, 400.0), 0.2);
        assert_eq!(gneg.f_at(400.0, 400.0), 0.2);
        // f(0) respects the clamp when f_start sits outside the caps.
        let gb = GradientSpec::rate_capped("A", "B", 0.1, 0.3, 100.0, 0.3, 1.0);
        assert_eq!(gb.f_at(0.0, 400.0), 0.3, "the cap governs, not f_start");
    }

    /// RateCapped sublayer sampling: the same inclusive depth grid as
    /// FixedSpan (z_i = t * i / (n - 1)), the clamp binding in the tail.
    /// rate 0.25 keeps the knee arithmetic exact (0.25 * k is binary
    /// exact), so the saturated tail rows are unambiguously bitwise
    /// pure material_b.
    #[test]
    fn rate_capped_sublayer_sampling_hits_the_knee() {
        let g = GradientSpec::rate_capped("A", "B", 0.1, 0.25, 100.0, 0.0, 1.0);
        // t = 400, n = 4: z = 0, 133.33, 266.67, 400 -> f = 0.1, 0.4333..,
        // 0.7666.., 1.0 (the raw value at z = 400 is 1.1, clamped).
        let z1 = 400.0 * 1.0 / 3.0;
        let z2 = 400.0 * 2.0 / 3.0;
        assert_eq!(g.f_sublayer(0, 4, 400.0), 0.1);
        assert_eq!(g.f_sublayer(1, 4, 400.0), 0.1 + 0.25 * (z1 / 100.0));
        assert_eq!(g.f_sublayer(2, 4, 400.0), 0.1 + 0.25 * (z2 / 100.0));
        assert_eq!(g.f_sublayer(3, 4, 400.0), 1.0, "tail row is bitwise pure");
    }

    /// The RateCapped refusals (D2's validation list).
    #[test]
    fn rate_capped_refusals() {
        let refuses = |g: GradientSpec| g.expansion_error(false).unwrap();
        let ok = GradientSpec::rate_capped("A", "B", 0.1, 0.3, 100.0, 0.0, 1.0);
        assert!(ok.expansion_error(false).is_none());
        let m = refuses(GradientSpec::rate_capped(
            "A",
            "B",
            0.1,
            f64::NAN,
            100.0,
            0.0,
            1.0,
        ));
        assert!(m.contains("rate") && m.contains("finite"), "{m}");
        let m = refuses(GradientSpec::rate_capped("A", "B", 0.1, 0.3, 0.0, 0.0, 1.0));
        assert!(m.contains("ref_thickness"), "{m}");
        let m = refuses(GradientSpec::rate_capped(
            "A", "B", 1.2, 0.3, 100.0, 0.0, 1.0,
        ));
        assert!(m.contains("f_start"), "{m}");
        let m = refuses(GradientSpec::rate_capped(
            "A", "B", 0.1, 0.3, 100.0, 0.8, 0.2,
        ));
        assert!(m.contains("f_min") && m.contains("f_max"), "{m}");
        let m = refuses(GradientSpec::rate_capped(
            "A", "B", 0.1, 0.3, 100.0, 0.0, 1.5,
        ));
        assert!(m.contains("f_max") && m.contains("[0, 1]"), "{m}");
        let m = refuses(GradientSpec::rate_capped(
            "A", "B", 0.1, 0.0, 100.0, 0.0, 1.0,
        ));
        assert!(m.contains("rate 0") && m.contains("degenerate"), "{m}");
    }

    /// The RateCapped serde form (one-key tagged map, like MixRule).
    #[test]
    fn rate_capped_serde_representation() {
        use serde_json::json;
        let g = GradientSpec::rate_capped("A", "B", 0.1, 0.3, 100.0, 0.0, 1.0);
        let v = serde_json::to_value(g.mode).unwrap();
        assert_eq!(
            v,
            json!({"RateCapped": {"f_start": 0.1, "rate": 0.3,
                                    "ref_thickness": 100.0, "f_min": 0.0, "f_max": 1.0}})
        );
        let back: GradientMode = serde_json::from_value(v).unwrap();
        assert_eq!(back, g.mode);
    }

    /// The inclusive sublayer convention: first sublayer IS f_start,
    /// last IS f_end (the legacy graded convention, bitwise across any
    /// sublayer count).
    #[test]
    fn f_sublayer_hits_the_endpoints_exactly() {
        let g = GradientSpec::fixed_span("A", "B", 0.1, 0.9);
        for n in [2u32, 3, 7, 64, 256] {
            assert_eq!(g.f_sublayer(0, n, 123.0), 0.1);
            assert_eq!(g.f_sublayer(n - 1, n, 123.0), 0.9);
            assert_eq!(g.f_sublayer(1, n, 123.0), 0.1 + 0.8 / f64::from(n - 1));
        }
        assert_eq!(g.f_sublayer(0, 1, 123.0), 0.1); // degenerate guard
    }

    /// The sublayer-count differential pin (same technique as the legacy
    /// `sub_layer_count` pin): integer boundaries must agree bit-exactly
    /// with the formula over randomized thicknesses and grids.
    #[test]
    fn gradient_count_matches_formula() {
        let oracle = |t: f64, wl: &[f64], na: &[Complex64], nb: &[Complex64]| {
            let lambda_min = wl.iter().copied().fold(f64::INFINITY, f64::min);
            let n_max = na.iter().chain(nb).map(|z| z.re).fold(0.0f64, f64::max);
            let step = MAX_GRADIENT_STEP_NM.min(lambda_min / (10.0 * n_max));
            (t / step).ceil().clamp(3.0, 64.0) as u32
        };
        // A fixed oracle grid: nk_a = 1.5, nk_b = 2.35 (n_max = 2.35).
        let wl: Vec<f64> = (0..9).map(|i| 400.0 + i as f64 * 100.0).collect();
        let na = vec![Complex64::new(1.5, 0.0); wl.len()];
        let nb = vec![Complex64::new(2.35, 0.0); wl.len()];
        // Sweep thicknesses across many boundaries, including exact
        // multiples of the step and values just off them.
        let step = MAX_GRADIENT_STEP_NM.min(400.0 / (10.0 * 2.35));
        let mut t = step * 0.25;
        let mut cases = 0;
        while t < 1400.0 {
            assert_eq!(
                gradient_sub_layer_count(t, &wl, &na, &nb, None),
                oracle(t, &wl, &na, &nb),
                "boundary divergence at t = {t}"
            );
            cases += 1;
            t += step * 0.125;
        }
        assert!(cases > 50, "sweep covered only {cases} cases");
        // Both extremes: the clamps bind.
        assert_eq!(gradient_sub_layer_count(1.0, &wl, &na, &nb, None), 3);
        assert_eq!(gradient_sub_layer_count(1e6, &wl, &na, &nb, None), 64);
        // A wider grid does not change the count (lambda_min governs).
        let wl2: Vec<f64> = wl.iter().chain([2000.0].iter()).copied().collect();
        assert_eq!(
            gradient_sub_layer_count(300.0, &wl2, &na, &nb, None),
            gradient_sub_layer_count(300.0, &wl, &na, &nb, None)
        );
    }

    #[test]
    fn override_clamps_into_two_to_256() {
        let wl = [550.0];
        let na = [Complex64::new(1.5, 0.0)];
        let nb = [Complex64::new(2.35, 0.0)];
        assert_eq!(gradient_sub_layer_count(300.0, &wl, &na, &nb, Some(1)), 2);
        assert_eq!(gradient_sub_layer_count(300.0, &wl, &na, &nb, Some(5)), 5);
        assert_eq!(
            gradient_sub_layer_count(300.0, &wl, &na, &nb, Some(1000)),
            256
        );
    }

    /// The expansion refusals, each naming both sides (F1.1 gate list).
    #[test]
    fn expansion_refusals_name_both_sides() {
        let ok = GradientSpec::fixed_span("A", "B", 0.0, 1.0);
        assert!(ok.expansion_error(false).is_none());
        let two_engines = ok.expansion_error(true).unwrap();
        assert!(two_engines.contains("gradient") && two_engines.contains("inhomogen"));
        let range = GradientSpec::fixed_span("A", "B", -0.1, 1.0);
        assert!(range.expansion_error(false).unwrap().contains("f_start"));
        assert!(range.expansion_error(false).unwrap().contains("[0, 1]"));
        let hi = GradientSpec::fixed_span("A", "B", 0.0, 1.5);
        assert!(hi.expansion_error(false).unwrap().contains("f_end"));
        let nan = GradientSpec::fixed_span("A", "B", f64::NAN, 1.0);
        assert!(nan.expansion_error(false).unwrap().contains("f_start"));
        let same = GradientSpec::fixed_span("TiO2", "TiO2", 0.0, 1.0);
        let m = same.expansion_error(false).unwrap();
        assert!(m.contains("TiO2") && m.contains("inhomogen"));
        let degenerate = GradientSpec::fixed_span("A", "B", 0.5, 0.5);
        let m = degenerate.expansion_error(false).unwrap();
        assert!(m.contains("f_start == f_end") && m.contains("0.5"));
        // The override advisory is separate from the fatal set.
        assert_eq!(ok.issues().len(), 0);
        let mut wide = ok.clone();
        wide.sublayers = Some(1000);
        let iss = wide.issues();
        assert_eq!(iss.len(), 1);
        assert!(!iss[0].is_error());
        assert!(iss[0].message.contains("clamped at expansion"));
    }

    /// MixRule serde: the schema-visible surface, named deliberately
    /// (A6). Unit variants ride as bare strings; parameterized variants
    /// as one-key maps (serde's externally-tagged representation).
    #[test]
    fn mix_rule_serde_representation_is_pinned() {
        use serde_json::json;
        assert_eq!(
            serde_json::to_value(MixRule::Looyenga).unwrap(),
            json!("Looyenga")
        );
        assert_eq!(
            serde_json::to_value(MixRule::Bruggeman {
                max_iter: 100,
                tol: 1e-9
            })
            .unwrap(),
            json!({"Bruggeman": {"max_iter": 100, "tol": 1e-9}})
        );
        assert_eq!(
            serde_json::to_value(MixRule::MoriTanaka { l: 1.0 / 3.0 }).unwrap(),
            json!({"MoriTanaka": {"l": 1.0 / 3.0}})
        );
        assert_eq!(
            serde_json::to_value(MixRule::PowerLaw { alpha: 0.5 }).unwrap(),
            json!({"PowerLaw": {"alpha": 0.5}})
        );
        let back: MixRule = serde_json::from_value(json!("Looyenga")).unwrap();
        assert_eq!(back, MixRule::Looyenga);
        let back: MixRule = serde_json::from_value(json!({"PowerLaw": {"alpha": 0.5}})).unwrap();
        assert_eq!(back, MixRule::PowerLaw { alpha: 0.5 });
    }

    /// The per-row oracle itself is checked against the raw kernels
    /// bitwise in `expansion.rs` (same file as the branch); here we pin
    /// the dispatch: every variant returns the same value its kernel
    /// does, and the default is Bruggeman.
    #[test]
    fn mix_row_dispatches_every_variant_bitwise() {
        let wl_len = 4;
        let nk_a = vec![Complex64::new(1.5, 0.01); wl_len];
        let nk_b = vec![Complex64::new(2.35, 0.05); wl_len];
        use ndarray::Array1;
        let b = Array1::from_vec(nk_b.clone());
        let a = Array1::from_vec(nk_a.clone());
        let cases: Vec<MixRule> = vec![
            MixRule::Bruggeman {
                max_iter: 100,
                tol: 1e-9,
            },
            MixRule::MaxwellGarnett,
            MixRule::Looyenga,
            MixRule::Lichtenecker,
            MixRule::MoriTanaka { l: 0.25 },
            MixRule::PowerLaw { alpha: 0.5 },
        ];
        for rule in cases {
            let want = match rule {
                MixRule::Bruggeman { max_iter, tol } => {
                    ema::eps_to_nk(ema::bruggeman(b.view(), a.view(), 0.3, max_iter, tol).view())
                        .to_vec()
                }
                MixRule::MaxwellGarnett => {
                    ema::eps_to_nk(ema::maxwell_garnett(b.view(), a.view(), 0.3).view()).to_vec()
                }
                MixRule::Looyenga => {
                    ema::eps_to_nk(ema::looyenga(b.view(), a.view(), 0.3).view()).to_vec()
                }
                MixRule::Lichtenecker => {
                    ema::eps_to_nk(ema::lichtenecker(b.view(), a.view(), 0.3).view()).to_vec()
                }
                MixRule::MoriTanaka { l } => {
                    ema::eps_to_nk(ema::mori_tanaka(b.view(), a.view(), 0.3, l).view()).to_vec()
                }
                MixRule::PowerLaw { alpha } => {
                    ema::eps_to_nk(ema::general_power_law(b.view(), a.view(), 0.3, alpha).view())
                        .to_vec()
                }
            };
            assert_eq!(mix_row(rule, &nk_b, &nk_a, 0.3), want, "{rule:?}");
        }
    }
}
