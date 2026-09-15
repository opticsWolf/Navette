// SPDX-License-Identifier: LGPL-3.0-or-later
//! Navette -- Rust Rewrite of Numba-optimized thin-film optical solver
//!
//! synthesis::cycle — one needle pass: repeated analytic insertions.
//!
//! Port of NeedleSynthesizer.run() with compute_p_function replaced by the
//! analytic sweep (needle_pass). The Python convergence test measured the
//! MF drop from inserting a 1 nm TEST needle:
//!     improvement_py = MF_now − MF(test needle, 1 nm) ≈ −P(z*)·1 nm
//! (first-order, since P = dF/dδ). We generalize to the seed thickness:
//!     predicted_improvement = −P_best · δ_seed
//! and stop when it falls below `convergence_threshold`. Thresholds are in
//! merit-units-per-nm-of-seed — recalibrate when porting configs.

use std::collections::HashMap;
use std::sync::Arc;

use num_complex::Complex64;

use crate::smatrix::synthesis::context::DesignContext;
use crate::smatrix::synthesis::needle_pass::{
    NeedlePassInput, build_needle_targets_env, build_scan_sites, run_needle_pass,
};
use crate::smatrix::synthesis::pipeline::SpectralInputs;
use crate::smatrix::synthesis::structure::{DesignStack, LayerSpec};

/// Knobs for the inner insertion loop (Python `NeedleConfig` subset).
#[derive(Clone, Debug)]
pub struct NeedleCycleConfig {
    /// Max insertions per pass (pipeline forwards `needles_per_cycle`).
    pub max_needles: usize,
    /// Stop when predicted improvement drops below this (merit units).
    pub convergence_threshold: f64,
    /// Seed layer thickness (nm).
    pub needle_seed_thickness_nm: f64,
    /// Scan resolution (nm).
    pub scan_step_nm: f64,
    /// Re-fold needle targets against the live sim each cycle (op-point
    /// masking for one-sided/banded/integral kinds + live PD gain shift).
    /// `false` keeps the static conservative fold (A/B benchmark switch).
    pub refold_per_cycle: bool,
}

impl Default for NeedleCycleConfig {
    fn default() -> Self {
        NeedleCycleConfig {
            max_needles: 10,
            convergence_threshold: 1e-4,
            needle_seed_thickness_nm: 5.0,
            scan_step_nm: 2.0,
            refold_per_cycle: true,
        }
    }
}

/// Record of one insertion cycle — mirrors Python `NeedleCycleResult`
/// with the scan-MF replaced by its analytic analog.
#[derive(Clone, Debug)]
pub struct NeedleCycleResult {
    pub cycle: usize,
    pub merit_before: f64,
    pub merit_after: f64,
    /// Most-negative P value at the chosen site (dF/dδ there).
    pub best_p: Option<f64>,
    /// −best_p · seed_thickness — the recalibrated convergence metric.
    pub predicted_improvement: Option<f64>,
    pub layer_count: usize,
    pub insertion: Option<Insertion>,
}

#[derive(Clone, Debug)]
pub struct Insertion {
    pub film_idx: usize,
    pub depth_into_layer_nm: f64,
    pub material: Arc<str>,
}

/// Contrast-material table: film material name → needle LayerSpec template.
/// (nk arrays must be evaluated on the simulation grid; thickness/flags of
/// the template are ignored — the seed is built fresh.)
pub type ContrastMap = HashMap<Arc<str>, LayerSpec>;

/// Run one full needle pass on `stack`.
///
/// Mirrors `synth.run(max_needles = cfg.max_needles)`:
/// initial optimization → [re-fold → scan → select → insert → optimize]×N.
///
/// The fold refreshes against the live sim each cycle (when
/// `cfg.refold_per_cycle` and the context simulates): one-sided/banded
/// demands mask at the operating point and PD gain shifts go live, so
/// insertion tracks the merit the optimizer actually sees. Any simulate
/// or fold failure keeps the previous fold (mock contexts without a
/// solver run the static conservative fold throughout).
/// F2.3: one candidate locus, agreed on by every environment.
///
/// A locus is a *shared design parameter* plus a step into it, never a row
/// — row 7 of one assembly and row 3 of another are the same physical
/// place when both answer to the same slot.
#[derive(Clone, Debug)]
struct JointCandidate {
    /// Index into `CompiledEnvironments::slots`.
    slot: usize,
    depth_into_layer_nm: f64,
    /// Seed material: the contrast partner whose sweep won.
    material: Arc<str>,
    /// P summed over environments.
    p: f64,
}

/// F2.3: the joint analytic sweep — K scans, one candidate.
///
/// Each environment scans ITS OWN assembly with ITS OWN fold (an operating
/// point masks one-sided and banded kinds, and environment 1's operating
/// point is not environment 0's), then every site that lands in the shared
/// design is routed **by design parameter** into a shared bucket and summed.
/// The sum is the right combination because the joint merit is a sum of the
/// per-environment merits, so dF/dδ at a shared locus is the sum of the
/// per-environment slopes — the same identity F2.2's residual concatenation
/// rests on, one derivative up.
///
/// Surroundings never deposit: their rows answer `None` to `slot_of_row`,
/// which is the second lock on the door the `needle` flag already holds.
///
/// Returns the winner (most negative ΣP, `None` when nothing improves), the
/// gain shift summed over the K folds for the convergence test, and the
/// number of admissible design sites — zero means "stack too thin".
fn joint_needle_sweep<C: DesignContext + ?Sized>(
    ctx: &C,
    stack: &DesignStack,
    spectral: &SpectralInputs,
    contrast: &ContrastMap,
    cfg: &NeedleCycleConfig,
) -> Result<(Option<JointCandidate>, f64, usize), String> {
    let envs = ctx
        .environments()
        .ok_or_else(|| "joint_needle_sweep: the context has no environments".to_string())?;
    let k = envs.n_envs();
    let stacks = envs.expand(stack)?;

    // Seed materials come from the DESIGN, not from the assemblies. A cover
    // material that happens to carry a contrast entry is not a candidate
    // seed, and sweeping it in the one environment that has it would open a
    // bucket no other environment can fill.
    let mats: Vec<&LayerSpec> = {
        let mut seen: Vec<Arc<str>> = Vec::new();
        let mut out = Vec::new();
        for (row, l) in stacks[0].films().iter().enumerate() {
            if envs.slot_of_row(0, row).is_some()
                && contrast.contains_key(&l.material)
                && !seen.contains(&l.material)
            {
                seen.push(l.material.clone());
                out.push(contrast.get(&l.material).unwrap());
            }
        }
        out
    };

    let mut gtot = 0.0f64;
    // (seed material, slot, step into the host) -> (depth, ΣP, environments seen)
    let mut buckets: HashMap<(Arc<str>, usize, usize), (f64, f64, usize)> = HashMap::new();
    let mut n_sites = 0usize;

    for (e, st) in stacks.iter().enumerate() {
        let mut fold = spectral.folds[e].clone();
        if cfg.refold_per_cycle
            && let Ok(sim) = ctx.simulate(st)
            && let Ok(f) = build_needle_targets_env(
                &spectral.spec,
                &spectral.angles_deg,
                &spectral.wavls,
                Some(&sim),
                e as u32,
            )
        {
            fold = f;
        }
        gtot += fold.phi_gain_shift.iter().sum::<f64>();

        let sa = st.solver_arrays();
        for mat in &mats {
            let input = NeedlePassInput {
                n_stack_cache: &sa.n_stack_cache,
                thicknesses: &sa.thicknesses,
                rough_types: &sa.rough_types,
                rough_vals: &sa.rough_vals,
                n_layers: sa.n_layers as usize,
                wavls: spectral.wavls.as_slice(),
                sin_theta: spectral.sin_theta.as_slice(),
                fold: &fold,
                needle_n_per_wav: &mat.nk,
                start_idx: 0,
                end_idx: (sa.n_layers - 1) as usize,
                calc_s: true,
                calc_p: true,
            };
            let res = run_needle_pass(&input, st.films(), st.spans(), cfg.scan_step_nm)?;
            // `build_scan_sites` walks rows in order and emits steps
            // 1..n_steps inside each admissible row, so the running count
            // per row IS the step index — and it is the same index in
            // every environment, because a design row has the same
            // thickness everywhere by construction (`expand` copies it).
            let mut step: HashMap<usize, usize> = HashMap::new();
            for (site, &p) in res.sites.iter().zip(res.p_profile.iter()) {
                let n = step.entry(site.film_idx).or_insert(0);
                let ordinal = *n;
                *n += 1;
                let Some(slot) = envs.slot_of_row(e, site.film_idx) else {
                    continue; // surroundings: inadmissible, and never shared
                };
                if e == 0 {
                    n_sites += 1;
                }
                let entry = buckets
                    .entry((mat.material.clone(), slot, ordinal))
                    .or_insert((site.depth_into_layer_nm, 0.0, 0));
                entry.1 += p;
                entry.2 += 1;
            }
        }
    }

    // §4.4's alignment assert, in its span-aware form: the same design
    // object must produce the same scan grid in every environment, so every
    // bucket is filled exactly K times. A short bucket means the assemblies
    // drifted apart, and summing it would sum over different loci.
    if let Some(((m, slot, ord), (_, _, n))) = buckets.iter().find(|(_, v)| v.2 != k) {
        return Err(format!(
            "joint needle sweep: the locus (slot {slot}, step {ord}, seed {m:?}) \
             appears in {n} of {k} environments - the assemblies no longer \
             share the design's scan grid"
        ));
    }

    // Deterministic order before the min: a hash map's iteration order is
    // not stable, and two loci can tie on P exactly (identical
    // environments, mirrored stacks).
    let mut items: Vec<((Arc<str>, usize, usize), (f64, f64, usize))> =
        buckets.into_iter().collect();
    items.sort_by(|a, b| a.0.cmp(&b.0));

    let mut best: Option<JointCandidate> = None;
    for ((material, slot, _), (depth, p, _)) in items {
        if p >= 0.0 {
            continue; // not an improvement anywhere
        }
        let better = match &best {
            None => true,
            Some(b) => p < b.p,
        };
        if better {
            best = Some(JointCandidate {
                slot,
                depth_into_layer_nm: depth,
                material,
                p,
            });
        }
    }
    Ok((best, gtot, n_sites))
}

pub fn run_needle_cycles<C: DesignContext + ?Sized>(
    ctx: &mut C,
    stack: &mut DesignStack,
    spectral: &SpectralInputs,
    contrast: &ContrastMap,
    cfg: &NeedleCycleConfig,
) -> Result<Vec<NeedleCycleResult>, String> {
    use crate::smatrix::synthesis::needle_pass::build_needle_targets;

    let wavls = &spectral.wavls;
    let sin_theta = &spectral.sin_theta;
    let mut fold = spectral.folds[0].clone();
    let mut history = Vec::new();

    // F2.3: §4.6's one branch, read once, outside every loop. `stack` is
    // the shared design object and IS environment 0's assembly; the other
    // K-1 are re-expressed from it inside the sweep.
    let multi = ctx.environments().is_some();

    // Initial optimization.
    ctx.optimize_thicknesses(stack)?;

    for cycle in 0..cfg.max_needles {
        if multi {
            if joint_cycle(ctx, stack, spectral, contrast, cfg, cycle, &mut history)? {
                break;
            }
            continue;
        }
        // Refresh the fold against the live operating point.
        if cfg.refold_per_cycle
            && let Ok(sim) = ctx.simulate(stack)
            && let Ok(f) = build_needle_targets(
                &spectral.spec,
                &spectral.angles_deg,
                &spectral.wavls,
                Some(&sim),
            )
        {
            fold = f;
        }
        // 1. Build candidate sites restricted to films whose material has a
        //    contrast entry (Python skips layers without a mapping) AND
        //    that are needle hosts. The flag check is load-bearing:
        //    interface slices and pinned graded rows share their carrier's
        //    material (hence a contrast entry) but must never host seeds.
        let sites = build_scan_sites(stack.films(), stack.spans(), cfg.scan_step_nm);
        let sites: Vec<_> = sites
            .into_iter()
            .filter(|s| {
                stack
                    .films()
                    .get(s.film_idx)
                    .map(|l| l.needle && contrast.contains_key(&l.material))
                    .unwrap_or(false)
            })
            .collect();
        if sites.is_empty() {
            break; // "Stack too thin for needle insertion."
        }

        // 2. Analytic sweep (both polarizations summed) per DISTINCT contrast
        //    material among admissible hosts; global best kept across sweeps.
        let sa = stack.solver_arrays();

        // Distinct contrast materials among admissible hosts.
        let mats: Vec<&LayerSpec> = {
            let mut seen: Vec<Arc<str>> = Vec::new();
            let mut out = Vec::new();
            for l in stack.films() {
                if contrast.contains_key(&l.material) && !seen.contains(&l.material) {
                    seen.push(l.material.clone());
                    out.push(contrast.get(&l.material).unwrap());
                }
            }
            out
        };

        let mut best: Option<(Insertion, f64)> = None;
        for mat in &mats {
            let input = NeedlePassInput {
                n_stack_cache: &sa.n_stack_cache,
                thicknesses: &sa.thicknesses,
                rough_types: &sa.rough_types,
                rough_vals: &sa.rough_vals,
                n_layers: sa.n_layers as usize,
                wavls: wavls.as_slice(),
                sin_theta: sin_theta.as_slice(),
                fold: &fold,
                needle_n_per_wav: &mat.nk,
                start_idx: 0,
                end_idx: (sa.n_layers - 1) as usize,
                calc_s: true,
                calc_p: true,
            };
            let res = run_needle_pass(&input, stack.films(), stack.spans(), cfg.scan_step_nm)?;
            if let Some((site, p)) = res.best() {
                let better = match &best {
                    None => true,
                    Some((_, bp)) => p < *bp,
                };
                if better {
                    best = Some((
                        Insertion {
                            film_idx: site.film_idx,
                            depth_into_layer_nm: site.depth_into_layer_nm,
                            material: mat.material.clone(),
                        },
                        p,
                    ));
                }
            }
        }

        // 3. Convergence check on the predicted improvement.
        let current_mf = ctx.evaluate_merit(stack)?;
        let result = match best {
            None => {
                // No improving site anywhere.
                history.push(NeedleCycleResult {
                    cycle: cycle + 1,
                    merit_before: current_mf,
                    merit_after: current_mf,
                    best_p: None,
                    predicted_improvement: None,
                    layer_count: stack.films().len(),
                    insertion: None,
                });
                break;
            }
            Some((ins, p)) => {
                // Inserting a seed of thickness δ also grows the
                // equivalent-medium reference by δ, so the true slope is
                // P(z*) + Σ gain_shift (0 for absolute-phase demand sets;
                // the site itself is unaffected, only this bookkeeping).
                let gtot: f64 = fold.phi_gain_shift.iter().sum();
                let predicted = -(p + gtot) * cfg.needle_seed_thickness_nm;
                if predicted < cfg.convergence_threshold {
                    history.push(NeedleCycleResult {
                        cycle: cycle + 1,
                        merit_before: current_mf,
                        merit_after: current_mf,
                        best_p: Some(p),
                        predicted_improvement: Some(predicted),
                        layer_count: stack.films().len(),
                        insertion: None,
                    });
                    break; // "Convergence reached — stopping."
                }

                // 4. Insert the seed (split host, Python `_insert_needle`).
                let seed_nk: Arc<[Complex64]> = contrast
                    .get(&stack.films()[ins.film_idx].material)
                    .map(|m| m.nk.clone())
                    .unwrap_or_else(|| vec![Complex64::new(1.0, 0.0); wavls.len()].into());
                let seed = LayerSpec {
                    material: ins.material.clone(),
                    nk: seed_nk,
                    d_nm: cfg.needle_seed_thickness_nm,
                    coherent: true,
                    rough_type: 0,
                    rough_val: 0.0,
                    optimize: true,
                    needle: true,
                };
                stack.insert_needle_seed(ins.film_idx, ins.depth_into_layer_nm, seed)?;

                // 5. Re-optimize.
                let new_mf = ctx.optimize_thicknesses(stack)?;

                NeedleCycleResult {
                    cycle: cycle + 1,
                    merit_before: current_mf,
                    merit_after: new_mf,
                    best_p: Some(p),
                    predicted_improvement: Some(predicted),
                    layer_count: stack.films().len(),
                    insertion: Some(ins),
                }
            }
        };
        history.push(result);
    }

    Ok(history)
}

/// F2.3: one cycle of the joint arm. Returns `true` when the pass is over.
///
/// Same five steps as the K == 1 body — sweep, convergence test, insert,
/// re-optimize, record — with two differences and no others. The sweep is
/// [`joint_needle_sweep`]'s (K scans summed by design parameter), and the
/// insertion happens TWICE: once on the shared design object and once on
/// the compile, so every environment's template splits at its own row for
/// the same parameter. Nothing else about the insertion changes: the same
/// seed, the same depth, the same `insert_needle_seed`, the same
/// re-optimization through the context.
fn joint_cycle<C: DesignContext + ?Sized>(
    ctx: &mut C,
    stack: &mut DesignStack,
    spectral: &SpectralInputs,
    contrast: &ContrastMap,
    cfg: &NeedleCycleConfig,
    cycle: usize,
    history: &mut Vec<NeedleCycleResult>,
) -> Result<bool, String> {
    let (best, gtot, n_sites) = joint_needle_sweep(&*ctx, stack, spectral, contrast, cfg)?;
    if n_sites == 0 {
        return Ok(true); // "Stack too thin for needle insertion."
    }
    let current_mf = ctx.evaluate_merit(stack)?;
    let stalled = |history: &mut Vec<NeedleCycleResult>, p, predicted| {
        history.push(NeedleCycleResult {
            cycle: cycle + 1,
            merit_before: current_mf,
            merit_after: current_mf,
            best_p: p,
            predicted_improvement: predicted,
            layer_count: stack.films().len(),
            insertion: None,
        });
    };

    let Some(cand) = best else {
        stalled(history, None, None);
        return Ok(true); // no improving locus in any environment
    };
    // The gain shift is summed over the K folds for the same reason P is:
    // growing the seed grows every environment's equivalent-medium
    // reference, and the joint merit charges all of them.
    let predicted = -(cand.p + gtot) * cfg.needle_seed_thickness_nm;
    if predicted < cfg.convergence_threshold {
        stalled(history, Some(cand.p), Some(predicted));
        return Ok(true); // "Convergence reached - stopping."
    }

    // Locus -> the shared design's own row. `host_row` answers the BULK
    // row, so an interface-carrying design film is still a legal host (N1).
    let row0 = ctx
        .environments()
        .and_then(|e| e.host_row(0, cand.slot))
        .ok_or_else(|| {
            format!(
                "joint needle: design parameter {} has no host row in \
                 environment 0",
                cand.slot
            )
        })?;
    let seed_nk: Arc<[Complex64]> = contrast
        .get(&stack.films()[row0].material)
        .map(|m| m.nk.clone())
        .unwrap_or_else(|| vec![Complex64::new(1.0, 0.0); spectral.wavls.len()].into());
    let seed = LayerSpec {
        material: cand.material.clone(),
        nk: seed_nk,
        d_nm: cfg.needle_seed_thickness_nm,
        coherent: true,
        rough_type: 0,
        rough_val: 0.0,
        optimize: true,
        needle: true,
    };
    stack.insert_needle_seed(row0, cand.depth_into_layer_nm, seed.clone())?;
    ctx.environments_mut()
        .ok_or_else(|| "joint needle: the context lost its environments".to_string())?
        .insert_seed(cand.slot, cand.depth_into_layer_nm, &seed)?;

    let new_mf = ctx.optimize_thicknesses(stack)?;
    history.push(NeedleCycleResult {
        cycle: cycle + 1,
        merit_before: current_mf,
        merit_after: new_mf,
        best_p: Some(cand.p),
        predicted_improvement: Some(predicted),
        layer_count: stack.films().len(),
        insertion: Some(Insertion {
            film_idx: row0,
            depth_into_layer_nm: cand.depth_into_layer_nm,
            material: cand.material,
        }),
    });
    Ok(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    // Only the tests build NeedleTargets by hand; the pass returns them.
    use crate::smatrix::synthesis::merit::{MeritSpec, SimCurves};
    use crate::smatrix::synthesis::needle_pass::NeedleTargets;
    use std::sync::Arc;

    /// Frozen context: merit constant, optimization a no-op (keeps the
    /// analytic AR condition intact so the scan sees exact zeros).
    struct StillCtx;
    impl DesignContext for StillCtx {
        fn evaluate_merit(&self, _s: &DesignStack) -> Result<f64, String> {
            Ok(0.0)
        }
        fn simulate(&self, _s: &DesignStack) -> Result<SimCurves, String> {
            Err("mock context has no simulator".into())
        }
        fn optimize_thicknesses(&mut self, s: &mut DesignStack) -> Result<f64, String> {
            self.evaluate_merit(s)
        }
    }

    /// air | G(200 nm, n = 1.52, host) | glass: all matched, so R is
    /// exactly the single-interface value ((n−1)/(n+1))² ≈ 0.04258 and
    /// T = 1 − R (lossless) at normal incidence.
    /// (A lossless R = 0 point would NOT discriminate: T = 1 − R
    /// identically there, so every intensity gradient vanishes with R's.)
    fn glass_stack() -> (DesignStack, f64) {
        let nw = 1;
        let r0 = ((1.52_f64 - 1.0) / (1.52_f64 + 1.0)).powi(2);
        let mut ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
        ambient.optimize = false;
        ambient.needle = false;
        let mut substrate = LayerSpec::constant("sub", 1.52, 0.0, 0.0, nw);
        substrate.optimize = false;
        substrate.needle = false;
        let stack = DesignStack::with_films(
            ambient,
            substrate,
            vec![LayerSpec::constant("G", 1.52, 0.0, 200.0, nw)],
        )
        .unwrap();
        (stack, r0)
    }

    fn fold_with(bucket: &str, target: f64) -> NeedleTargets {
        let zero = || (vec![0.0f64; 1], vec![0.0f64; 1]);
        let demand = || (vec![target; 1], vec![1.0f64; 1]);
        let (r, t) = match bucket {
            "r" => (demand(), zero()),
            "t" => (zero(), demand()),
            _ => unreachable!(),
        };
        NeedleTargets {
            r,
            t,
            a: zero(),
            rb: zero(),
            tb: zero(),
            ab: zero(),
            phi: [zero(), zero(), zero(), zero()],
            phi_gain_shift: [0.0; 4],
            grad_r: vec![0.0f64; 1],
            grad_t: vec![0.0f64; 1],
        }
    }

    fn contrast_h() -> ContrastMap {
        let mut m = ContrastMap::new();
        m.insert(Arc::from("G"), LayerSpec::constant("H", 2.35, 0.0, 0.0, 1));
        m
    }

    fn cycle_cfg() -> NeedleCycleConfig {
        NeedleCycleConfig {
            max_needles: 2,
            convergence_threshold: 1e-4,
            needle_seed_thickness_nm: 5.0,
            scan_step_nm: 2.0,
            refold_per_cycle: true,
        }
    }

    /// SpectralInputs wrapping a hand-built fold (empty spec — the
    /// StillCtx simulator errors, so the static fold stands).
    fn spectral_of(fold: NeedleTargets) -> SpectralInputs {
        SpectralInputs {
            wavls: vec![1000.0],
            sin_theta: vec![0.0],
            folds: vec![fold],
            spec: MeritSpec::new(),
            angles_deg: vec![0.0],
        }
    }

    /// Non-host rows (needle=false) are never insertion sites even when
    /// their material has a contrast entry (interface slices, pinned
    /// graded spans). The run breaks with the stack untouched.
    #[test]
    fn non_host_rows_never_seed() {
        let nw = 1;
        let mut ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
        ambient.optimize = false;
        ambient.needle = false;
        let mut substrate = LayerSpec::constant("sub", 1.52, 0.0, 0.0, nw);
        substrate.optimize = false;
        substrate.needle = false;
        let mut host = LayerSpec::constant("G", 1.52, 0.0, 200.0, nw);
        host.needle = false; // contrast entry exists, host flag refuses
        let mut stack = DesignStack::with_films(ambient, substrate, vec![host]).unwrap();
        let mut ctx = StillCtx;
        let spectral = spectral_of(fold_with("t", 2.0));
        let hist = run_needle_cycles(&mut ctx, &mut stack, &spectral, &contrast_h(), &cycle_cfg())
            .unwrap();
        assert!(hist.is_empty());
        assert_eq!(stack.films().len(), 1);
        assert!((stack.films()[0].d_nm - 200.0).abs() < 1e-12);
    }

    #[test]
    fn satisfied_r_demand_inserts_nothing() {
        // R demanded at its exact value → residual ~1e-16 → dust profile,
        // convergence gate stops: no insertion, film count untouched.
        let (mut stack, r0) = glass_stack();
        let mut ctx = StillCtx;
        let spectral = spectral_of(fold_with("r", r0));
        let hist = run_needle_cycles(&mut ctx, &mut stack, &spectral, &contrast_h(), &cycle_cfg())
            .unwrap();
        assert_eq!(hist.len(), 1);
        assert!(hist[0].insertion.is_none());
        assert!(hist[0].best_p.map(|p| p > -1e-9).unwrap_or(true));
        assert_eq!(stack.films().len(), 1);
    }

    #[test]
    fn live_refold_masks_satisfied_above_demand() {
        // Glass R ≈ 0.0426 against Above 0.01: satisfied (sim ≥ target)
        // at the operating point. The static conservative fold would
        // insert (Above folds active without a sim); the per-cycle live
        // re-fold masks it, so nothing is inserted. Real solver context.
        use crate::smatrix::synthesis::evaluator::SmatrixContext;
        use crate::smatrix::synthesis::merit::{
            ConstraintKind, MeritKey, MeritTarget, SimTransform,
        };
        use crate::smatrix::synthesis::thick_opt::LmConfig;

        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: crate::smatrix::synthesis::merit::CurveId::Rs,
        });
        spec.add_target(MeritTarget {
            key_idx: k as u32,
            env_idx: 0,
            wavelengths: vec![1000.0].into(),
            kind: ConstraintKind::Above,
            transform: SimTransform::Linear,
            norm_factor: 1.0,
            normalized_targets: vec![0.01].into(),
            tolerances: vec![0.05].into(),
            band: vec![].into(),
            phase: false,
            differential_passes: None,
            integral: false,
            weight: 1.0,
            count_norm: None,
        })
        .unwrap();
        let mut ctx = SmatrixContext {
            wavls: vec![1000.0],
            sin_theta: vec![0.0],
            spec: spec.clone(),
            clamp_min_nm: 2.0,
            clamp_max_nm: 1000.0,
            lm: LmConfig::default(),
            clamp_accumulator: crate::smatrix::synthesis::structure::ClampReport::default(),
            envs: None,
            thin_layer_policy: crate::smatrix::synthesis::config::ThinLayerPolicy::Remove,
        };
        let spectral = SpectralInputs::from_spec(&spec, &[0.0], &[1000.0]).unwrap();
        let (mut stack, _) = glass_stack();
        let hist = run_needle_cycles(&mut ctx, &mut stack, &spectral, &contrast_h(), &cycle_cfg())
            .unwrap();
        assert!(hist.iter().all(|h| h.insertion.is_none()));
        assert_eq!(stack.films().len(), 1);
    }

    #[test]
    fn violated_t_demand_drives_insertion() {
        // Same stack, T ≈ 0.957 against a T = 0 demand: inserting H
        // disrupts the matching and lowers T, so the scan must find a
        // negative-P site and the host film must split.
        let (mut stack, _) = glass_stack();
        let mut ctx = StillCtx;
        let spectral = spectral_of(fold_with("t", 0.0));
        let hist = run_needle_cycles(&mut ctx, &mut stack, &spectral, &contrast_h(), &cycle_cfg())
            .unwrap();
        // max_needles = 2 and T stays violated → two insertions
        // (1 → 3 → 5 films), both H into the original host lineage.
        assert_eq!(hist.len(), 2);
        for h in &hist {
            let ins = h.insertion.as_ref().expect("expected an insertion");
            assert!(h.best_p.unwrap() < 0.0);
            assert_eq!(ins.material.as_ref(), "H");
        }
        assert_eq!(stack.films().len(), 5);
    }
}
