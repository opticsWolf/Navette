// SPDX-License-Identifier: LGPL-3.0-or-later
//! Two-phase stack expansion: design layers → solver rows.
//!
//! Transliteration of `navette.structure.expander._LayerExpander`:
//! phase 1 resolves every entry's bulk (group scaling, no RNG); phase 2
//! emits rows in traversal order with error draws. The mirror rule
//! (inverted-run planes take the traversal predecessor's flags; incident
//! edges clean), donor-side carve with buffered-row rescale, owner+carrier
//! Looyenga mix at f=0.5, roughness-follows-plane, and the draw order
//! (thick, nk, rough, iface, inhg) are all preserved exactly.
//!
//! Randomness is full-Rust: `seed = Some` → reproducible `StdRng`,
//! `None` → thread RNG. Forward output must equal Python bit-for-bit
//! (pinned below); Monte-Carlo paths agree statistically (§9.2).

use std::collections::HashMap;

use num_complex::Complex64;
use rand::RngCore;
use rand::SeedableRng;
use rand::rngs::{StdRng, ThreadRng};

use crate::structure::enums::{ErrorType, RoughnessType};
use crate::structure::group::{Group, gauss_draw, unif_draw};
use crate::structure::layer::Layer;
use crate::structure::providers::MaterialProvider;

/// One emitted solver row-block: rows `[start, end)` belong to logical
/// entry `logical` (interface slices resolve to their carrier).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    pub start: usize,
    pub end: usize,
    pub logical: usize,
    /// Row `start` is an interface slice (derived row: never a free
    /// parameter, never a needle host). Set at emission time.
    pub slice: bool,
    /// First bulk row (B9, F0.1). Rows `[bulk_start, end)` carry the
    /// layer's own (group-scaled, possibly graded) index; row `start` —
    /// when `slice` — is the interface slice, carved from the carrier's
    /// authored thickness (expansion carves, it does not append). Carried
    /// as a field so later items read one source instead of re-deriving
    /// `start + slice as usize`; `is_singleton_bulk` asserts the equality
    /// so the two cannot silently drift apart.
    pub bulk_start: usize,
}

impl Span {
    /// N1: a span is *singleton-bulk* iff it holds exactly one non-slice
    /// row. A plain film (1 row) and an interface-carrying plain film
    /// (2 rows, one of them the slice) both qualify; graded and gradient
    /// spans do not. Say "singleton-bulk", never "singleton": the row
    /// count is not the predicate.
    pub fn is_singleton_bulk(&self) -> bool {
        debug_assert_eq!(
            self.bulk_start,
            self.start + usize::from(self.slice),
            "Span::bulk_start disagrees with the arithmetic predicate"
        );
        self.end - self.bulk_start == 1
    }
}

/// Engine-ready arrays: row-major `indices` (`n_rows × n_wavelengths`).
#[derive(Debug, Clone, PartialEq)]
pub struct SolverArrays {
    pub thicknesses: Vec<f64>,
    pub indices: Vec<Complex64>,
    pub n_wavelengths: usize,
    pub incoherent: Vec<bool>,
    pub rough_types: Vec<i32>,
    pub rough_vals: Vec<f64>,
}

impl SolverArrays {
    pub fn n_rows(&self) -> usize {
        self.thicknesses.len()
    }

    pub fn row(&self, i: usize) -> &[Complex64] {
        let s = i * self.n_wavelengths;
        &self.indices[s..s + self.n_wavelengths]
    }
}

/// Expansion options: error draws on/off + draw stream.
#[derive(Debug, Clone, Copy)]
pub struct ExpandOptions {
    pub apply_errors: bool,
    pub seed: Option<u64>,
}

impl ExpandOptions {
    pub fn deterministic() -> Self {
        Self {
            apply_errors: false,
            seed: None,
        }
    }
}

enum AnyRng {
    Seeded(Box<StdRng>),
    Thread(ThreadRng),
}

impl AnyRng {
    fn rng(&mut self) -> &mut dyn RngCore {
        match self {
            AnyRng::Seeded(r) => r,
            AnyRng::Thread(r) => r,
        }
    }
}

/// Expand `(layer, inverted)` traversal entries to solver rows + spans.
///
/// `wavelengths` is the simulation grid (provider arrays resolve on it).
/// Empty sequences are refused. Deterministic output (errors off) must
/// equal the Python expander bit-for-bit.
pub fn expand(
    seq: &[(Layer, bool)],
    provider: &dyn MaterialProvider,
    wavelengths: &[f64],
    groups: &HashMap<String, Group>,
    opts: ExpandOptions,
) -> Result<(SolverArrays, Vec<Span>), String> {
    if seq.is_empty() {
        return Err(
            "_LayerExpander.expand: No layers to expand. Empty layer sequence provided."
                .to_string(),
        );
    }
    let default_group = Group::new("_default_");
    let group_of = |material: &str| groups.get(material).unwrap_or(&default_group);

    // ---- Phase 1: deterministic bulk resolution (no RNG). ----
    let mut bulk_nk: Vec<Vec<Complex64>> = Vec::with_capacity(seq.len());
    let mut bulk_t: Vec<f64> = Vec::with_capacity(seq.len());
    let mut owner_of: Vec<Option<usize>> = Vec::with_capacity(seq.len());
    for (k, (layer, inv)) in seq.iter().enumerate() {
        let group = group_of(layer.material.as_str());
        let base = provider.nk(layer.material.as_str(), wavelengths)?;
        let scaled = if group.n_factor != 1.0 || group.k_factor != 1.0 {
            base.iter()
                .map(|z| Complex64::new(z.re * group.n_factor, z.im * group.k_factor))
                .collect()
        } else {
            base
        };
        bulk_nk.push(scaled);
        bulk_t.push(layer.thickness * group.thick_factor + group.thick_summand);
        owner_of.push(if *inv && k > 0 && seq[k - 1].1 {
            Some(k - 1)
        } else if !inv && k > 0 {
            Some(k)
        } else {
            None
        });
    }

    // ---- Phase 2: emission in traversal order (RNG draws here). ----
    //
    // `seed = None` is deliberately asymmetric, and the asymmetry is the honest
    // reading of the request rather than an oversight:
    //   * errors ON  -> the thread RNG. The caller asked for randomness and did
    //     not pin it, so the run is genuinely not reproducible and does not
    //     pretend to be by silently seeding itself with a constant.
    //   * errors OFF -> a seeded(0) RNG that nothing ever draws from. Cheaper
    //     than making the field an Option, and it cannot reach the output.
    let mut rng = match opts.seed {
        Some(s) => AnyRng::Seeded(Box::new(StdRng::seed_from_u64(s))),
        None if opts.apply_errors => AnyRng::Thread(rand::rng()),
        None => AnyRng::Seeded(Box::new(StdRng::seed_from_u64(0))), // unused; errors off
    };
    let mut em = Emission::with_capacity(seq.len());
    let mut prev_eff_nk: Option<Vec<Complex64>> = None;

    for k in 0..seq.len() {
        emit_entry(
            &mut em,
            seq,
            k,
            groups,
            &bulk_nk[k],
            bulk_t[k],
            owner_of[k],
            opts,
            &mut rng,
            &mut prev_eff_nk,
            provider,
            wavelengths,
        )?;
    }

    let n_rows = em.col_thick.len();
    let sa = SolverArrays {
        thicknesses: em.col_thick,
        indices: em.col_nk,
        n_wavelengths: wavelengths.len(),
        incoherent: em.col_coh.iter().map(|c| !c).collect(),
        rough_types: em.col_r_type,
        rough_vals: em.col_r_val,
    };
    debug_assert_eq!(sa.indices.len(), n_rows * wavelengths.len());
    Ok((sa, em.spans))
}

/// Mutable per-`expand` emission state: exactly the accumulators one
/// iteration of the emission loop reads and writes (the F0.1 extraction
/// contract, open decision 15). Owning the vectors keeps `expand`'s
/// locals and `emit_entry`'s body byte-for-byte the old code.
struct Emission {
    col_thick: Vec<f64>,
    col_nk: Vec<Complex64>,
    col_coh: Vec<bool>,
    col_r_val: Vec<f64>,
    col_r_type: Vec<i32>,
    spans: Vec<Span>,
    bulk_spans: Vec<(usize, usize)>,
    err_nk: Vec<Vec<Complex64>>,
    err_t: Vec<f64>,
}

impl Emission {
    fn with_capacity(seq_len: usize) -> Self {
        Self {
            col_thick: Vec::new(),
            col_nk: Vec::new(),
            col_coh: Vec::new(),
            col_r_val: Vec::new(),
            col_r_type: Vec::new(),
            spans: Vec::new(),
            bulk_spans: Vec::new(),
            err_nk: Vec::with_capacity(seq_len),
            err_t: Vec::with_capacity(seq_len),
        }
    }
}

/// Emit entry `k`'s rows into the running column set.
///
/// Extracted verbatim from `expand`'s emission loop at F0.1 so that F1.7's
/// `refresh_profiles` and construction are the same code by construction
/// (open decision 15: "refresh equals construction" as a property, not a
/// test result). Pure code motion: the parameter list is exactly what the
/// loop body reads, nothing else; no push reordered, no draw order
/// changed.
///
/// F1.1 adds one branch beside the legacy graded one (the gradient span)
/// and, with it, the first way this function can fail (a gradient spec
/// that the gates did not catch, or an endpoint the provider does not
/// carry). The `Result` return is that failure path only: every legacy
/// branch is unchanged and returns `Ok(())` - no push reordered, no draw
/// order changed, still true by inspection.
#[allow(clippy::too_many_arguments)]
fn emit_entry(
    em: &mut Emission,
    seq: &[(Layer, bool)],
    k: usize,
    groups: &HashMap<String, Group>,
    bulk_nk_k: &[Complex64],
    bulk_t_k: f64,
    owner_k: Option<usize>,
    opts: ExpandOptions,
    rng: &mut AnyRng,
    prev_eff_nk: &mut Option<Vec<Complex64>>,
    provider: &dyn MaterialProvider,
    wavelengths: &[f64],
) -> Result<(), String> {
    let default_group = Group::new("_default_");
    let group_of = |material: &str| groups.get(material).unwrap_or(&default_group);
    let (layer, inv) = (&seq[k].0, &seq[k].1);
    let Emission {
        col_thick,
        col_nk,
        col_coh,
        col_r_val,
        col_r_type,
        spans,
        bulk_spans,
        err_nk,
        err_t,
    } = em;

    let group = group_of(layer.material.as_str());
    let mut layer_nk = bulk_nk_k.to_vec();
    let mut layer_thickness = bulk_t_k;
    if opts.apply_errors {
        if group.error_mask[crate::structure::enums::ErrorMask::Thickness as usize] != 0 {
            layer_thickness = group.thickness_error(layer_thickness, rng.rng());
        }
        let me = crate::structure::enums::ErrorMask::NReal as usize;
        let ke = crate::structure::enums::ErrorMask::NImag as usize;
        if group.error_mask[me] != 0 || group.error_mask[ke] != 0 {
            perturb_nk(&mut layer_nk, group, rng.rng());
        }
    }
    layer_thickness = layer_thickness.max(0.0);
    err_nk.push(layer_nk.clone());
    err_t.push(layer_thickness);
    let o = owner_k;

    // Bulk roughness (owner-group draws; RNG order thick, nk, rough, ...).
    let (current_roughness, rtype) = if *inv {
        match o {
            Some(oi) => {
                let (olayer, _) = &seq[oi];
                let ogroup = group_of(olayer.material.as_str());
                let mut r = (olayer.roughness + ogroup.roughness_summand).max(0.0);
                if opts.apply_errors
                    && ogroup.error_mask[crate::structure::enums::ErrorMask::Roughness as usize]
                        != 0
                {
                    r = ogroup.sr_roughness_error(r, rng.rng());
                }
                (r, olayer.rough_type as i32)
            }
            None => (0.0, RoughnessType::None as i32),
        }
    } else {
        let mut r = (layer.roughness + group.roughness_summand).max(0.0);
        if opts.apply_errors
            && group.error_mask[crate::structure::enums::ErrorMask::Roughness as usize] != 0
        {
            r = group.sr_roughness_error(r, rng.rng());
        }
        (r, layer.rough_type as i32)
    };

    let start = col_thick.len();
    let mut emitted_slice = false;

    // Plane slice (flag owner's group governs summand + draws).
    if let Some(oi) = o {
        let (olayer, _) = &seq[oi];
        let ogroup = group_of(olayer.material.as_str());
        if olayer.interface {
            let mut t_interface = olayer.interface_thickness + ogroup.interface_summand;
            if opts.apply_errors
                && ogroup.error_mask[crate::structure::enums::ErrorMask::Interface as usize] != 0
            {
                t_interface = ogroup.interface_error(t_interface, rng.rng());
            }
            let carve_total = if oi == k { layer_thickness } else { err_t[oi] };
            t_interface = t_interface.min(carve_total);
            let mix = if oi == k {
                layer_thickness -= t_interface;
                let prev = prev_eff_nk.as_deref().unwrap_or(&err_nk[k]);
                looyenga_mix(&layer_nk, prev)
            } else {
                let (start_o, end_o) = bulk_spans[oi];
                if carve_total > 0.0 && end_o > start_o {
                    let scale = (carve_total - t_interface) / carve_total;
                    for row in &mut col_thick[start_o..end_o] {
                        *row *= scale;
                    }
                }
                looyenga_mix(&err_nk[oi], &err_nk[k])
            };
            push_row(
                col_thick,
                col_nk,
                col_coh,
                col_r_val,
                col_r_type,
                t_interface,
                mix,
                true,
                0.0,
                RoughnessType::None as i32,
            );
            emitted_slice = true;
        }
    }

    let bulk_start = col_thick.len();
    let sub = layer.sub_layer_count();
    if let Some(grad) = &layer.gradient
        && layer_thickness > 0.0
    {
        // F1.1: the gradient span. Beside the legacy graded branch, not
        // inside it (N7: the two branches deliberately use different
        // sublayer-thickness conventions - the legacy one divides
        // uniformly, so sum(d) == thickness only to float precision;
        // this one has the LAST sublayer absorb the remainder, so
        // sum(d) == thickness exactly).
        //
        // Endpoint resolution (D3 step 1): each endpoint resolves
        // through the provider like every other lookup - same grid
        // assertion, and the carrier's GROUP policy applies through each
        // endpoint's OWN material (scaling factors, and the nk error
        // channels when errors are on - a gradient film in an error run
        // shifts as one body, the way the fabrication offset applies to
        // any other material). When material_a names the carrier, its
        // endpoint IS this layer's phase-1 nk, already scaled and drawn.
        // The delta channels are inert here: mixtures have no inh_delta
        // (D2) - the group's inh policy simply never runs on this branch.
        if let Some(msg) = grad.expansion_error(layer.inhomogen) {
            return Err(format!("layer '{}': {msg}", layer.material));
        }
        let resolve = |name: &str, rng: &mut AnyRng| -> Result<Vec<Complex64>, String> {
            let g = group_of(name);
            let base = provider.nk(name, wavelengths).map_err(|e| {
                crate::structure::gradient::endpoint_error(&layer.material, grad, e)
            })?;
            let mut nk = if g.n_factor != 1.0 || g.k_factor != 1.0 {
                base.iter()
                    .map(|z| Complex64::new(z.re * g.n_factor, z.im * g.k_factor))
                    .collect()
            } else {
                base
            };
            if opts.apply_errors {
                let me = crate::structure::enums::ErrorMask::NReal as usize;
                let ke = crate::structure::enums::ErrorMask::NImag as usize;
                if g.error_mask[me] != 0 || g.error_mask[ke] != 0 {
                    perturb_nk(&mut nk, g, rng.rng());
                }
            }
            Ok(nk)
        };
        // Draw order thick -> nk -> rough -> iface (legacy order) gains
        // one nk draw per endpoint, in A-then-B order; only gradient
        // layers draw here, so no legacy stream is affected.
        let nk_a = if grad.material_a == layer.material {
            layer_nk.clone()
        } else {
            resolve(&grad.material_a, rng)?
        };
        let nk_b = if grad.material_b == layer.material {
            layer_nk.clone()
        } else {
            resolve(&grad.material_b, rng)?
        };
        let n_sub = crate::structure::gradient::gradient_sub_layer_count(
            layer_thickness,
            wavelengths,
            &nk_a,
            &nk_b,
            grad.sublayers,
        );
        let step = layer_thickness / f64::from(n_sub);
        // Deposition-order rows (f_start face first); `inv` reverses the
        // ROW order afterwards - the physical flip, same as legacy - so
        // the roughness-bearing interface row stays the one adjacent to
        // what was emitted before this span ("first sublayer", legacy
        // convention, keeps interface physics identical both ways).
        // Profile sampling is INCLUSIVE (f_sublayer): the first row IS
        // f_start, the last IS f_end - the legacy factors' convention,
        // and what the endpoint-row bitwise gates sample.
        let mut rows: Vec<(f64, Vec<Complex64>)> = Vec::with_capacity(n_sub as usize);
        for i in 0..n_sub {
            let d_i = if i + 1 == n_sub {
                layer_thickness - step * f64::from(n_sub - 1) // remainder
            } else {
                step
            };
            let f_i = grad.f_sublayer(i, n_sub, layer_thickness);
            rows.push((
                d_i,
                crate::structure::gradient::mix_row(grad.ema, &nk_b, &nk_a, f_i),
            ));
        }
        if *inv {
            rows.reverse();
        }
        for (ix, (d_i, row_nk)) in rows.into_iter().enumerate() {
            push_row(
                col_thick,
                col_nk,
                col_coh,
                col_r_val,
                col_r_type,
                d_i,
                row_nk,
                layer.coherent,
                if ix == 0 { current_roughness } else { 0.0 },
                if ix == 0 {
                    rtype
                } else {
                    RoughnessType::None as i32
                },
            );
        }
    } else if layer.inhomogen && sub > 1 {
        // F1.3: the grading strength is the mode's `delta_layer` (the
        // ONE source every reader consumes - B6), flowed through the
        // FROZEN combination order. `Fixed` reads the authored
        // `inh_delta` (byte-identical path); `RateCapped` computes
        // D2's min-formula and then clamps the COMBINED nominal to
        // [-cap, cap] BEFORE the stochastic draw - the cap binds the
        // nominal, the noise is allowed to exceed it (clamping draws
        // would bias Monte-Carlo statistics).
        let delta_layer = layer.delta_layer();
        let mut current_delta = (delta_layer + group.inh_delta_summand) * 0.5;
        if let Some(cap) = layer.inh_mode.nominal_cap() {
            current_delta = current_delta.clamp(-cap, cap);
        }
        if opts.apply_errors
            && group.error_mask[crate::structure::enums::ErrorMask::InhDelta as usize] != 0
        {
            current_delta = group.inh_delta_error(current_delta, rng.rng());
        }
        let mut factors: Vec<f64> = (0..sub)
            .map(|i| 1.0 - current_delta + 2.0 * current_delta * f64::from(i) / f64::from(sub - 1))
            .collect();
        if *inv {
            factors.reverse();
        }
        let step_t = layer_thickness / f64::from(sub);
        for (ix, f) in factors.iter().enumerate() {
            let row_nk: Vec<Complex64> = layer_nk.iter().map(|z| z * f).collect();
            push_row(
                col_thick,
                col_nk,
                col_coh,
                col_r_val,
                col_r_type,
                step_t,
                row_nk,
                layer.coherent,
                if ix == 0 { current_roughness } else { 0.0 },
                if ix == 0 {
                    rtype
                } else {
                    RoughnessType::None as i32
                },
            );
        }
    } else {
        push_row(
            col_thick,
            col_nk,
            col_coh,
            col_r_val,
            col_r_type,
            layer_thickness,
            layer_nk.clone(),
            layer.coherent,
            current_roughness,
            rtype,
        );
    }
    spans.push(Span {
        start,
        end: col_thick.len(),
        logical: k,
        slice: emitted_slice,
        bulk_start,
    });
    bulk_spans.push((bulk_start, col_thick.len()));

    *prev_eff_nk = Some(layer_nk);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn push_row(
    col_thick: &mut Vec<f64>,
    col_nk: &mut Vec<Complex64>,
    col_coh: &mut Vec<bool>,
    col_r_val: &mut Vec<f64>,
    col_r_type: &mut Vec<i32>,
    thickness: f64,
    nk: Vec<Complex64>,
    coherent: bool,
    rough_val: f64,
    rough_type: i32,
) {
    col_thick.push(thickness);
    col_nk.extend(nk);
    col_coh.push(coherent);
    col_r_val.push(rough_val);
    col_r_type.push(rough_type);
}

/// Owner+carrier Looyenga mix at f=0.5 (native kernel; eps→n at insertion).
fn looyenga_mix(owner_nk: &[Complex64], other_nk: &[Complex64]) -> Vec<Complex64> {
    use ndarray::Array1;
    let ni = Array1::from_vec(owner_nk.to_vec());
    let nh = Array1::from_vec(other_nk.to_vec());
    let eps = crate::materials::ema::looyenga(ni.view(), nh.view(), 0.5);
    crate::materials::ema::eps_to_nk(eps.view()).to_vec()
}

/// Array-wide nk perturbation with scalar draws (systematic fabrication
/// offset, not per-wavelength noise — mirrors `Group._apply_error` on
/// arrays: n floored at 0, k untouched).
fn perturb_nk(nk: &mut [Complex64], group: &Group, rng: &mut dyn RngCore) {
    use crate::structure::enums::ErrorMask;
    let me = ErrorMask::NReal as usize;
    let ke = ErrorMask::NImag as usize;
    if group.error_mask[me] != 0 {
        let (dn_abs, dn_rel) = channel_draws(group.n_error_type, &group.n_error_params, rng);
        for z in nk.iter_mut() {
            z.re = (z.re + dn_abs + dn_rel * z.re).max(0.0);
        }
    }
    if group.error_mask[ke] != 0 {
        let (dk_abs, dk_rel) = channel_draws(group.k_error_type, &group.k_error_params, rng);
        for z in nk.iter_mut() {
            z.im += dk_abs + dk_rel * z.im;
        }
    }
}

/// Scalar (abs, rel) draws for one channel in legacy order.
fn channel_draws(
    error_type: ErrorType,
    params: &crate::structure::group::ErrorParams,
    rng: &mut dyn RngCore,
) -> (f64, f64) {
    match error_type {
        ErrorType::Gaussian => (
            gauss_draw(params.abs_mean_delta_g, params.abs_std_dev, rng),
            gauss_draw(params.rel_mean_delta_g, params.rel_std_dev, rng),
        ),
        ErrorType::Uniform => (
            unif_draw(params.abs_variance, rng),
            unif_draw(params.rel_variance, rng),
        ),
        ErrorType::Combined => (
            gauss_draw(params.abs_mean_delta_g, params.abs_std_dev, rng)
                + unif_draw(params.abs_variance, rng),
            gauss_draw(params.rel_mean_delta_g, params.rel_std_dev, rng)
                + unif_draw(params.rel_variance, rng),
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::layer::Layer;
    use crate::structure::providers::{DictProvider, Entry};
    use std::collections::{HashMap, HashSet};

    const WL: [f64; 2] = [1000.0, 1500.0];

    fn mats() -> DictProvider {
        let mut entries = HashMap::new();
        entries.insert(
            "glass".to_string(),
            Entry::Array(vec![Complex64::new(1.52, 0.0), Complex64::new(1.51, 0.0)]),
        );
        entries.insert(
            "TiO2".to_string(),
            Entry::Array(vec![
                Complex64::new(2.35, 0.01),
                Complex64::new(2.33, 0.008),
            ]),
        );
        DictProvider::with_grid(entries, WL.to_vec()).unwrap()
    }

    fn flat_seq() -> Vec<(Layer, bool)> {
        vec![
            (Layer::film(0.0, "glass"), false),
            (Layer::film(50.0, "TiO2"), false),
            (Layer::film(0.0, "glass"), false),
        ]
    }

    fn close(a: &[f64], b: &[f64]) {
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b) {
            assert!((x - y).abs() < 1e-12, "{x} vs {y}");
        }
    }

    fn close_c(a: &[Complex64], re: &[f64], im: &[f64]) {
        assert_eq!(a.len(), re.len());
        for (z, (r, i)) in a.iter().zip(re.iter().zip(im)) {
            assert!((z.re - r).abs() < 1e-12, "{} vs {r}", z.re);
            assert!((z.im - i).abs() < 1e-12, "{} vs {i}", z.im);
        }
    }

    /// Oracle twin: flat forward (Python FLAT).
    #[test]
    fn flat_forward_matches_python() {
        let (sa, spans) = expand(
            &flat_seq(),
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        assert_eq!(sa.n_rows(), 3);
        close(&sa.thicknesses, &[0.0, 50.0, 0.0]);
        close_c(sa.row(1), &[2.35, 2.33], &[0.01, 0.008]);
        assert_eq!(sa.incoherent, vec![false; 3]);
        assert_eq!(sa.rough_types, vec![0; 3]);
        assert_eq!(
            spans
                .iter()
                .map(|s| (s.start, s.end, s.logical, s.bulk_start))
                .collect::<Vec<_>>(),
            vec![(0, 1, 0, 0), (1, 2, 1, 1), (2, 3, 2, 2)]
        );
        assert!(spans.iter().all(|s| !s.slice));
        assert!(spans.iter().all(|s| s.is_singleton_bulk()));
    }

    /// F0.1: the predicate's four cases, no pipeline involved. The second
    /// is the N1 case that would otherwise have been discovered by a
    /// moved fingerprint: an interface-carrying plain film is a 2-row
    /// span and STILL singleton-bulk.
    #[test]
    fn is_singleton_bulk_four_cases() {
        let cases: Vec<(Vec<(Layer, bool)>, bool)> = vec![
            // plain film
            (vec![(Layer::film(50.0, "glass"), false)], true),
            // plain film with interface (slice + bulk)
            (
                vec![(
                    Layer {
                        interface: true,
                        interface_thickness: 5.0,
                        ..Layer::film(50.0, "glass")
                    },
                    false,
                )],
                true,
            ),
            // graded film
            (
                vec![(
                    Layer {
                        inhomogen: true,
                        inh_delta: 0.2,
                        ..Layer::film(50.0, "glass")
                    },
                    false,
                )],
                false,
            ),
            // graded film with interface (slice + graded bulk)
            (
                vec![(
                    Layer {
                        inhomogen: true,
                        inh_delta: 0.2,
                        interface: true,
                        interface_thickness: 5.0,
                        ..Layer::film(50.0, "glass")
                    },
                    false,
                )],
                false,
            ),
        ];
        for (seq, want) in cases {
            let (_, spans) = expand(
                &seq,
                &mats(),
                &WL,
                &HashMap::new(),
                ExpandOptions::deterministic(),
            )
            .unwrap();
            assert_eq!(spans.len(), 1, "{seq:?}");
            assert_eq!(
                spans[0].is_singleton_bulk(),
                want,
                "span {:?} (case: {:?})",
                spans[0],
                want
            );
        }
    }

    /// Oracle twin: scaling + slice + grading + roughness (Python FULL).
    #[test]
    fn full_stack_matches_python() {
        let mut layers = vec![(Layer::film(0.0, "glass"), false)];
        let mut l = Layer::film(50.0, "TiO2");
        l.roughness = 2.0;
        l.rough_type = RoughnessType::Gaussian;
        l.interface = true;
        l.interface_thickness = 6.0;
        l.inhomogen = true;
        l.inh_delta = 0.2;
        layers.push((l, false));
        let mut groups = HashMap::new();
        let mut g = Group::new("TiO2");
        g.thick_factor = 1.1;
        g.thick_summand = 1.0;
        g.n_factor = 1.05;
        g.k_factor = 0.5;
        g.roughness_summand = 1.0;
        g.interface_summand = 2.0;
        groups.insert("TiO2".to_string(), g);
        let (sa, spans) = expand(
            &layers,
            &mats(),
            &WL,
            &groups,
            ExpandOptions::deterministic(),
        )
        .unwrap();
        assert_eq!(sa.n_rows(), 13);
        assert!(spans[1].slice && spans[1].start == 1);
        // B9: the bulk range is carried on the span, and the slice leads it.
        assert_eq!(spans[1].bulk_start, 2);
        assert!(
            spans
                .iter()
                .all(|s| s.is_singleton_bulk() == (s.end - s.start - usize::from(s.slice) == 1))
        );
        assert!(!spans[1].is_singleton_bulk());
        let mut want_t = vec![0.0, 8.0];
        want_t.extend(vec![48.0 / 11.0; 11]);
        close(&sa.thicknesses, &want_t);
        // Slice mix + first/last graded factors (2.4675/2.4465 x 0.9 … x 1.1).
        close_c(
            sa.row(1),
            &[1.974736420464, 1.959531653889],
            &[0.002321082511, 0.00185737228],
        );
        close_c(sa.row(2), &[2.22075, 2.20185], &[0.0045, 0.0036]);
        close_c(sa.row(12), &[2.71425, 2.69115], &[0.0055, 0.0044]);
        assert_eq!(sa.rough_types[2], 4);
        assert!((sa.rough_vals[2] - 3.0).abs() < 1e-12);
        assert!(sa.rough_vals[3..].iter().all(|v| *v == 0.0));
    }

    /// Oracle twin: whole-chain mirror appends cleanly (Python INV).
    #[test]
    fn inverted_mirror_matches_python() {
        let mut seq = flat_seq();
        let mut mirror: Vec<(Layer, bool)> = flat_seq()
            .into_iter()
            .rev()
            .map(|(l, _)| (l, true))
            .collect();
        seq.append(&mut mirror);
        let (sa, _) = expand(
            &seq,
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        assert_eq!(sa.n_rows(), 6);
        close(&sa.thicknesses, &[0.0, 50.0, 0.0, 0.0, 50.0, 0.0]);
        close_c(sa.row(4), &[2.35, 2.33], &[0.01, 0.008]);
    }

    #[test]
    fn empty_sequence_refused() {
        assert!(
            expand(
                &[],
                &mats(),
                &WL,
                &HashMap::new(),
                ExpandOptions::deterministic()
            )
            .is_err()
        );
    }

    #[test]
    fn error_paths_deterministic_per_seed() {
        let mut groups = HashMap::new();
        let mut g = Group::new("TiO2");
        g.error_mask = [1, 1, 1, 1, 1, 1];
        groups.insert("TiO2".to_string(), g);
        let opts = ExpandOptions {
            apply_errors: true,
            seed: Some(11),
        };
        let (a, _) = expand(&flat_seq(), &mats(), &WL, &groups, opts).unwrap();
        let (b, _) = expand(&flat_seq(), &mats(), &WL, &groups, opts).unwrap();
        assert_eq!(a, b);
        // Draws perturbed something (thickness floored draws or nk offsets).
        let plain = mats().nk("TiO2", &WL).unwrap();
        let mut expected = Vec::new();
        for _ in 0..3 {
            expected.extend(plain.clone());
        }
        assert!(a.thicknesses != vec![0.0, 50.0, 0.0] || a.indices != expected);
    }

    // ------------------------------------------------------------------
    // F1.1 - gradient spans
    // ------------------------------------------------------------------

    use crate::materials::MixRule;
    use crate::structure::gradient::{
        GradientMode, GradientSpec, InhMode, ProfileShape, gradient_sub_layer_count, mix_row,
    };

    fn grad_seq(t: f64, f_start: f64, f_end: f64, inv: bool) -> Vec<(Layer, bool)> {
        let mut l = Layer::film(t, "TiO2");
        l.gradient = Some(GradientSpec::fixed_span("glass", "TiO2", f_start, f_end));
        vec![
            (Layer::film(0.0, "glass"), false),
            (l, inv),
            (Layer::film(0.0, "glass"), false),
        ]
    }

    /// Oracle twin: every emitted sublayer row is bitwise the direct
    /// `mix_row` call at its formula `f_i` (same kernel, same inputs),
    /// the `f_i` recompute is from the mode formula, sum(d) == thickness
    /// EXACTLY (N7: the last sublayer absorbs the remainder), and the
    /// roughness lands on the first emitted sublayer only.
    #[test]
    fn gradient_rows_are_bitwise_the_ema_oracle() {
        let t = 137.0; // deliberately not a multiple of any step
        let (sa, spans) = expand(
            &grad_seq(t, 0.1, 0.9, false),
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        let nk_a = mats().nk("glass", &WL).unwrap();
        let nk_b = mats().nk("TiO2", &WL).unwrap();
        let n = gradient_sub_layer_count(t, &WL, &nk_a, &nk_b, None) as usize;
        assert_eq!(spans.len(), 3);
        let sp = &spans[1];
        assert_eq!(sp.end - sp.start, n);
        assert!(!sp.is_singleton_bulk());
        assert_eq!(sa.n_rows(), n + 2);
        // sum(d) == t exactly, in f64 addition order.
        let mut acc = 0.0f64;
        for r in sp.start..sp.end {
            acc += sa.thicknesses[r];
        }
        assert_eq!(acc, t, "sublayer thicknesses must sum to the layer exactly");
        for (ix, r) in (sp.start..sp.end).enumerate() {
            let _d_i = if ix + 1 == n {
                t - (t / n as f64) * (n - 1) as f64
            } else {
                t / n as f64
            };
            // Inclusive endpoint sampling (f_sublayer): f_i depends only
            // on the row index, never on the thickness; the last row is
            // exactly f_end.
            let f_i = if ix == 0 {
                0.1
            } else if ix + 1 == n {
                0.9
            } else {
                0.1 + 0.8 * ix as f64 / (n - 1) as f64
            };
            let want = mix_row(
                MixRule::Bruggeman {
                    max_iter: 100,
                    tol: 1e-9,
                },
                &nk_b,
                &nk_a,
                f_i,
            );
            assert_eq!(
                sa.row(r),
                &want[..],
                "row {ix} must be bitwise the EMA at f={f_i}"
            );
        }
        // The profile direction: first emitted row's f < last's.
        let first = sa.row(sp.start)[0].re;
        let last = sa.row(sp.end - 1)[0].re;
        assert!(
            first < last,
            "profile must run f_start -> f_end in emission order"
        );
        // Roughness on the FIRST emitted sublayer only (legacy convention).
        let mut l = Layer::film(t, "TiO2");
        l.roughness = 2.0;
        l.rough_type = RoughnessType::Gaussian;
        l.gradient = Some(GradientSpec::fixed_span("glass", "TiO2", 0.0, 1.0));
        let (sa2, sp2) = expand(
            &[(l, false)],
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        assert_eq!(
            sa2.rough_types[sp2[0].start],
            RoughnessType::Gaussian as i32
        );
        assert!(sa2.rough_vals[sp2[0].start] > 0.0);
        assert!(
            sa2.rough_vals[sp2[0].start + 1..sp2[0].end]
                .iter()
                .all(|v| *v == 0.0)
        );
        assert_eq!(
            sa2.rough_types[sp2[0].start + 1],
            RoughnessType::None as i32
        );
    }

    /// Thickness-independence (D6 G-thickness, type a): the same
    /// FixedSpan spec at two thicknesses - endpoint rows bitwise equal,
    /// counts differ.
    #[test]
    fn gradient_fixed_span_is_thickness_independent() {
        let (sa50, sp50) = expand(
            &grad_seq(50.0, 0.0, 1.0, false),
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        let (sa500, sp500) = expand(
            &grad_seq(500.0, 0.0, 1.0, false),
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        let a = sa50.row(sp50[1].start);
        let b = sa500.row(sp500[1].start);
        assert_eq!(
            a, b,
            "first sublayer (f = f_start) must not depend on thickness"
        );
        let a = sa50.row(sp50[1].end - 1);
        let b = sa500.row(sp500[1].end - 1);
        assert_eq!(
            a, b,
            "last sublayer (f = f_end) must not depend on thickness"
        );
        assert_ne!(sp50[1].end - sp50[1].start, sp500[1].end - sp500[1].start);
    }

    /// Inversion reverses the ROW order (the physical flip): the
    /// inverted span's rows are the forward span's rows reversed, so the
    /// profile direction flips and the interface row stays adjacent to
    /// what precedes the span.
    #[test]
    fn gradient_inversion_reverses_row_order() {
        let (fwd, _) = expand(
            &grad_seq(200.0, 0.0, 1.0, false),
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        let (inv, spinv) = expand(
            &grad_seq(200.0, 0.0, 1.0, true),
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        assert_eq!(fwd.n_rows(), inv.n_rows());
        let n = spinv[1].end - spinv[1].start;
        for i in 0..n {
            assert_eq!(
                fwd.row(1 + i),
                inv.row(1 + (n - 1 - i)),
                "inverted row {i} must be the forward run's row {n}/{i}"
            );
        }
    }

    /// Zero thickness never enters the gradient branch (no division by
    /// the layer thickness): the uniform fallback emits the carrier row.
    #[test]
    fn gradient_zero_thickness_falls_back_to_uniform() {
        let mut l = Layer::film(0.0, "TiO2");
        l.gradient = Some(GradientSpec::fixed_span("glass", "TiO2", 0.0, 1.0));
        let (sa, spans) = expand(
            &[(l, false)],
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        assert_eq!(sa.n_rows(), 1);
        assert!(spans[0].is_singleton_bulk());
        let plain = mats().nk("TiO2", &WL).unwrap();
        assert_eq!(sa.row(0), &plain[..]);
    }

    /// The refusals, each naming both sides (F1.1 gate list).
    #[test]
    fn gradient_refusals_at_expand() {
        let refuses = |l: Layer| {
            expand(
                &[(l, false)],
                &mats(),
                &WL,
                &HashMap::new(),
                ExpandOptions::deterministic(),
            )
            .expect_err("expected a refusal")
        };
        let mut l = Layer::film(50.0, "TiO2");
        l.gradient = Some(GradientSpec::fixed_span("glass", "TiO2", 0.0, 1.0));
        l.inhomogen = true;
        let m = refuses(l);
        assert!(m.contains("gradient") && m.contains("inhomogen"), "{m}");
        let mut l = Layer::film(50.0, "TiO2");
        l.gradient = Some(GradientSpec::fixed_span("TiO2", "TiO2", 0.0, 1.0));
        let m = refuses(l);
        assert!(m.contains("TiO2") && m.contains("identical"), "{m}");
        let mut l = Layer::film(50.0, "TiO2");
        l.gradient = Some(GradientSpec::fixed_span("glass", "TiO2", 0.4, 0.4));
        let m = refuses(l);
        assert!(m.contains("f_start == f_end"), "{m}");
        let mut l = Layer::film(50.0, "TiO2");
        l.gradient = Some(GradientSpec::fixed_span("glass", "TiO2", -0.2, 1.0));
        let m = refuses(l);
        assert!(m.contains("[0, 1]"), "{m}");
        // Provider door: the absent endpoint is named, with the layer.
        let mut l = Layer::film(50.0, "TiO2");
        l.gradient = Some(GradientSpec::fixed_span("unobtainium", "TiO2", 0.0, 1.0));
        let m = refuses(l);
        assert!(m.contains("layer 'TiO2'"), "{m}");
        assert!(m.contains("unobtainium"), "{m}");
    }

    /// Error runs: the endpoints draw under their OWN groups' nk
    /// channels (A-then-B), so a gradient film in an error run shifts as
    /// one body; output is deterministic per seed.
    #[test]
    fn gradient_error_draws_follow_endpoint_groups() {
        let mut groups = HashMap::new();
        let mut g = Group::new("TiO2");
        g.error_mask = [0, 0, 1, 0, 0, 0]; // NReal only
        groups.insert("TiO2".to_string(), g);
        let opts = ExpandOptions {
            apply_errors: true,
            seed: Some(7),
        };
        let (a, _) = expand(
            &grad_seq(50.0, 0.0, 1.0, false),
            &mats(),
            &WL,
            &groups,
            opts,
        )
        .unwrap();
        let (b, _) = expand(
            &grad_seq(50.0, 0.0, 1.0, false),
            &mats(),
            &WL,
            &HashMap::new(),
            opts,
        )
        .unwrap();
        // The material_b group draws: rows differ from the plain spectra.
        let plain_b = mats().nk("TiO2", &WL).unwrap();
        let shifted = (0..a.n_rows()).any(|r| a.row(r)[0] != plain_b[0]);
        assert!(shifted, "the nk draw must reach the emitted rows");
        // Deterministic per seed (same config twice).
        let (c, _) = expand(
            &grad_seq(50.0, 0.0, 1.0, false),
            &mats(),
            &WL,
            &groups,
            opts,
        )
        .unwrap();
        assert_eq!(a, c);
        let _ = b; // the second run exercises the no-group path
    }

    /// F1.2 (D6 G-thickness, type b): the saturating RateCapped film's
    /// tail rows are BITWISE pure `material_b` (the saturation proof -
    /// EMA at exactly the cap), the head runs the formula, and the
    /// thickness that doubles the film moves the knee, not the tail.
    #[test]
    fn rate_capped_saturation_tail_is_bitwise_pure() {
        let mk = |t: f64| {
            let mut l = Layer::film(t, "TiO2");
            l.gradient = Some(GradientSpec::rate_capped(
                "glass", "TiO2", 0.1, 0.25, 100.0, 0.0, 1.0,
            ));
            l
        };
        let nk_a = mats().nk("glass", &WL).unwrap();
        let nk_b = mats().nk("TiO2", &WL).unwrap();
        let pure_b = mix_row(
            MixRule::Bruggeman {
                max_iter: 100,
                tol: 1e-9,
            },
            &nk_b,
            &nk_a,
            1.0,
        );
        for t in [400.0, 800.0] {
            let (sa, spans) = expand(
                &[(mk(t), false)],
                &mats(),
                &WL,
                &HashMap::new(),
                ExpandOptions::deterministic(),
            )
            .unwrap();
            let sp = &spans[0];
            let n = (sp.end - sp.start) as u32;
            // Tail rows (raw value past the cap) are bitwise pure_b.
            let tail_start = (0..n)
                .find(|i| 0.1 + 0.25 * (t * f64::from(*i) / f64::from(n - 1)) / 100.0 >= 1.0)
                .unwrap();
            for i in tail_start..n {
                assert_eq!(
                    sa.row(sp.start + i as usize),
                    &pure_b[..],
                    "tail row {i} must be bitwise pure material_b"
                );
            }
            // The head rows follow the formula.
            for i in 0..tail_start {
                let f_i = 0.1 + 0.25 * (t * f64::from(i) / f64::from(n - 1)) / 100.0;
                let want = mix_row(
                    MixRule::Bruggeman {
                        max_iter: 100,
                        tol: 1e-9,
                    },
                    &nk_b,
                    &nk_a,
                    f_i,
                );
                assert_eq!(sa.row(sp.start + i as usize), &want[..], "head row {i}");
            }
        }
        // Negative rate saturates to pure material_a at the low end.
        let mut l = Layer::film(400.0, "TiO2");
        l.gradient = Some(GradientSpec::rate_capped(
            "glass", "TiO2", 0.9, -0.25, 100.0, 0.2, 1.0,
        ));
        let (sa, spans) = expand(
            &[(l, false)],
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        let pure_a = mix_row(
            MixRule::Bruggeman {
                max_iter: 100,
                tol: 1e-9,
            },
            &nk_b,
            &nk_a,
            0.2,
        );
        let sp = &spans[0];
        let n = (sp.end - sp.start) as u32;
        let tail_start = (0..n)
            .find(|i| (0.9 - 0.25 * (400.0 * f64::from(*i) / f64::from(n - 1)) / 100.0) <= 0.2)
            .unwrap();
        for i in tail_start..n {
            assert_eq!(
                sa.row(sp.start + i as usize),
                &pure_a[..],
                "low tail row {i}"
            );
        }
    }

    // ------------------------------------------------------------------
    // F1.3 - InhMode::RateCapped: the B6 reader agreement + the draw
    // ------------------------------------------------------------------

    /// The frozen legacy arithmetic parameterized by the mode's delta:
    /// the emitted ramp is `1 - delta_nom ..= 1 + delta_nom` with
    /// `delta_nom = (delta_layer + summand) * 0.5` (the RateCapped
    /// nominal clamped to [-cap, cap]), and the emitted row count is
    /// exactly `sub_layer_count()` (B6 reader 1: emission is the truth).
    #[test]
    fn rate_capped_emission_matches_frozen_arithmetic() {
        let mut l = Layer::film(400.0, "TiO2");
        l.inhomogen = true;
        l.inh_mode = InhMode::RateCapped {
            rate: 0.05,
            ref_thickness: 100.0,
            cap: 0.3,
        };
        // delta_layer = min(0.05*400/100, 0.3) = 0.2; nominal = 0.1.
        let (sa, spans) = expand(
            &[(l.clone(), false)],
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        let sp = &spans[0];
        assert_eq!(
            sp.end - sp.start,
            l.sub_layer_count() as usize,
            "B6: emission == prediction"
        );
        let base = mats().nk("TiO2", &WL).unwrap();
        let sub = sp.end - sp.start;
        for i in 0..sub {
            let f = 1.0 - 0.1 + 2.0 * 0.1 * i as f64 / (sub - 1) as f64;
            let want: Vec<Complex64> = base.iter().map(|z| z * f).collect();
            assert_eq!(sa.row(sp.start + i), &want[..], "row {i} (f={f})");
        }
    }

    /// B6: the statistical twin - the cap binds the NOMINAL, the draw
    /// is unclamped, so at saturation the drawn spread exceeds the
    /// clamped nominal and its mean sits AT the nominal (clamping the
    /// draw would bias the mean below it).
    #[test]
    fn rate_capped_error_draw_survives_the_clamp_boundary() {
        let mut groups = HashMap::new();
        let mut g = Group::new("TiO2");
        g.error_mask = [0, 0, 0, 0, 1, 0]; // InhDelta only
        g.inh_delta_error_type = crate::structure::enums::ErrorType::Gaussian;
        g.inh_delta_error_params.abs_mean_delta_g = 0.0;
        g.inh_delta_error_params.abs_std_dev = 0.05;
        groups.insert("TiO2".to_string(), g);
        let mut l = Layer::film(1200.0, "TiO2");
        l.inhomogen = true;
        // Saturated: delta_layer = cap = 0.3 -> nominal 0.15.
        l.inh_mode = InhMode::RateCapped {
            rate: 0.05,
            ref_thickness: 100.0,
            cap: 0.3,
        };
        let sub = l.sub_layer_count();
        let mut deltas: Vec<f64> = Vec::new();
        for seed in 0..200u64 {
            let (sa, spans) = expand(
                &[(l.clone(), false)],
                &mats(),
                &WL,
                &groups,
                ExpandOptions {
                    apply_errors: true,
                    seed: Some(seed + 1),
                },
            )
            .unwrap();
            let sp = &spans[0];
            assert_eq!(sp.end - sp.start, sub as usize);
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for r in sp.start..sp.end {
                let f = sa.row(r)[0].re / 2.35; // the factor
                lo = lo.min(f);
                hi = hi.max(f);
            }
            deltas.push((hi - lo) / 2.0); // = |current_delta| drawn
        }
        let mean: f64 = deltas.iter().sum::<f64>() / deltas.len() as f64;
        // The nominal is 0.15; a clamped draw would pull the mean below
        // it, an unclamped zero-mean draw leaves it there.
        assert!(
            mean > 0.14,
            "mean drawn delta {mean} collapsed below the nominal"
        );
        assert!(
            deltas.iter().any(|d| *d > 0.16),
            "no draw exceeded the nominal cap - the draw is being clamped"
        );
    }

    // ------------------------------------------------------------------
    // F1.6 - the scale-invariance twin, written before the feature
    // ------------------------------------------------------------------

    /// THE PHYSICS CLAIM of F1.6 (U2), asserted before the optimizer
    /// touches any of it: a profile whose fraction depends only on the
    /// sublayer's FRACTIONAL position scales exactly - stretch the layer
    /// and every sublayer's thickness scales while every sublayer's nk
    /// stays bitwise put.
    ///
    /// Engine: `GradientMode::FixedSpan` with the `sublayers` override
    /// pinning the count (the plan's "row count pinned to the same
    /// value"). Thicknesses 128/256 nm are binary-exact so the row
    /// thicknesses double BITWISE, including the remainder row
    /// (128 = 4 x 32: no rounding anywhere).
    #[test]
    fn f16_scale_invariance_fixed_span_bitwise() {
        let mk = |t: f64| {
            let mut l = Layer::film(t, "TiO2");
            l.gradient = Some(GradientSpec {
                material_a: "glass".to_string(),
                material_b: "TiO2".to_string(),
                ema: MixRule::Bruggeman {
                    max_iter: 100,
                    tol: 1e-12,
                },
                mode: GradientMode::FixedSpan {
                    f_start: 0.3,
                    f_end: 0.7,
                },
                shape: ProfileShape::Linear,
                sublayers: Some(4),
            });
            l
        };
        let (sa1, _) = expand(
            &[(mk(128.0), false)],
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        let (sa2, _) = expand(
            &[(mk(256.0), false)],
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        assert_eq!(sa1.n_rows(), 4);
        assert_eq!(sa2.n_rows(), 4, "the pinned count");
        for r in 0..4 {
            assert_eq!(sa1.row(r), sa2.row(r), "row {r}: nk bitwise equal");
        }
        for r in 0..4 {
            assert_eq!(
                sa2.thicknesses[r],
                2.0 * sa1.thicknesses[r],
                "row {r}: thickness exactly doubled"
            );
            assert_eq!(sa1.thicknesses[r], 32.0, "binary-exact split");
        }
    }

    /// The legacy engine's half of the same claim. Its row count has no
    /// override, and no thickness doubling preserves the count
    /// (t^0.4 grows by 2^0.5 per doubling, so t and 2t always ceil
    /// differently) - the plan's 100/200 example was idealized. Same
    /// count pair instead: 100 nm and 120 nm both give
    /// ceil(t^0.4) = 7, i.e. the SAME 7-row ramp. The nk rows are then
    /// bitwise equal (the ramp depends on i/(sub-1) only) and every row
    /// thickness is bitwise the direct division t/sub, which is what
    /// makes the profile scale-free at fixed count.
    #[test]
    fn f16_scale_invariance_inhomogen_same_count() {
        let mk = |t: f64| {
            let mut l = Layer::film(t, "TiO2");
            l.inhomogen = true;
            l.inh_delta = 0.2;
            l
        };
        assert_eq!(mk(100.0).sub_layer_count(), mk(120.0).sub_layer_count());
        let (sa1, _) = expand(
            &[(mk(100.0), false)],
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        let (sa2, _) = expand(
            &[(mk(120.0), false)],
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        assert_eq!(sa1.n_rows(), sa2.n_rows());
        for r in 0..sa1.n_rows() {
            assert_eq!(sa1.row(r), sa2.row(r), "row {r}: nk bitwise equal");
        }
        let sub = mk(100.0).sub_layer_count() as f64;
        for r in 0..sa1.n_rows() {
            assert_eq!(sa1.thicknesses[r], 100.0 / sub);
            assert_eq!(sa2.thicknesses[r], 120.0 / sub);
        }
    }

    #[test]
    fn spans_cover_rows_exactly_once() {
        let (sa, spans) = expand(
            &flat_seq(),
            &mats(),
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap();
        let mut seen = HashSet::new();
        for s in &spans {
            for r in s.start..s.end {
                assert!(seen.insert(r), "row {r} covered twice");
            }
        }
        assert_eq!(seen.len(), sa.n_rows());
    }
}
