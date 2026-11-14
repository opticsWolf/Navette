// SPDX-License-Identifier: LGPL-3.0-or-later
//! synthesis::staged — A1: the needle pipeline as callable stages.
//!
//! The macro-loop in `cycle` is convenient for a full run but too coarse for
//! a memetic operator: it inserts, re-optimizes, and repeats internally. This
//! module exposes the same operations one at a time:
//!
//! * [`scan`] — analytic P sweep over the current stack, returning every
//!   improving `(site, material)` candidate, best first.
//! * [`insert`] — put one candidate into the stack (split a host or add a
//!   boundary film), without re-optimizing.
//! * [`polish`] — bounded LM finish over a chosen free-set subset.
//! * [`cleanup`] — merge/prune/re-optimize, with the clamp report attached.
//!
//! `run_needle_cycles` remains the macro-loop; this module is the additive
//! stage surface the GA bridge calls. It does not change the loop's
//! behaviour, so the existing benches and tests stay valid.

use std::sync::Arc;

use crate::smatrix::synthesis::cleanup::{CleanupResult, cleanup_design};
use crate::smatrix::synthesis::context::DesignContext;
use crate::smatrix::synthesis::cycle::{
    Insertion, MaterialPool, NeedleCycleConfig, new_needle, rough_boundary_start_nm,
};
use crate::smatrix::synthesis::needle_pass::{
    NeedlePassInput, ScanSite, build_boundary_sites, build_needle_targets, build_scan_sites,
    needle_pass_scan,
};
use crate::smatrix::synthesis::pipeline::SpectralInputs;
use crate::smatrix::synthesis::structure::{ClampReport, DesignStack};

/// One scored insertion candidate.
///
/// `p` is the analytic P value. It is `None` for a rough-boundary candidate,
/// which is scored by evaluating the merit with the needle in place: P is
/// first order in the seed thickness and cannot see the roughness jump.
#[derive(Clone, Debug)]
pub struct NeedleCandidate {
    pub site: ScanSite,
    pub material_index: usize,
    pub material: Arc<str>,
    /// Analytic P at the site; negative means insertion lowers the merit.
    pub p: Option<f64>,
    /// Seed thickness the candidate was scored at.
    pub seed_nm: f64,
    /// Predicted merit change at `seed_nm` (negative = improvement).
    pub predicted_delta_merit: f64,
}

impl NeedleCandidate {
    pub fn is_improving(&self) -> bool {
        self.predicted_delta_merit < 0.0
    }
}

/// The result of one stage-1 scan.
#[derive(Clone, Debug, Default)]
pub struct NeedleScan {
    /// Improving candidates, best (most negative predicted delta) first.
    pub candidates: Vec<NeedleCandidate>,
    /// Merit of the stack before any insertion.
    pub current_merit: f64,
    /// Number of `(site, material)` pairs actually swept.
    pub sites_scanned: usize,
}

impl NeedleScan {
    pub fn best(&self) -> Option<&NeedleCandidate> {
        self.candidates.first()
    }

    pub fn is_empty(&self) -> bool {
        self.candidates.is_empty()
    }
}

/// Which films a [`polish`] call may move.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PolishScope {
    /// Every film currently flagged `optimize`.
    All,
    /// The last `k` currently-free films, substrate side first.
    LastK(usize),
    /// Currently-free films whose material name or label matches one of these.
    Named(Vec<String>),
}

/// What one polish did.
#[derive(Clone, Debug)]
pub struct PolishReport {
    pub scope: PolishScope,
    pub merit_before: f64,
    pub merit_after: f64,
    pub free_layers: usize,
}

/// Cleanup knobs; `Default` is the cleanup pass used by the pipeline.
#[derive(Clone, Copy, Debug, Default)]
pub struct CleanupOptions {
    pub min_thickness_nm: Option<f64>,
    pub max_removals: Option<usize>,
    pub reoptimize: bool,
}

/// Cleanup result plus whatever clamp work it caused.
#[derive(Clone, Debug)]
pub struct CleanupReport {
    pub result: CleanupResult,
    pub clamp: Option<ClampReport>,
}

/// A1 scan: analytic P over the current stack, every improving candidate.
///
/// The flat (one-environment) path only. A joint spec is refused here rather
/// than folded against environment 0; the joint macro-loop remains the place
/// where K environments are handled together.
pub fn scan<C: DesignContext + ?Sized>(
    ctx: &C,
    stack: &DesignStack,
    spectral: &SpectralInputs,
    pool: &MaterialPool,
    cfg: &NeedleCycleConfig,
) -> Result<NeedleScan, String> {
    if spectral.spec.n_envs() != 1 {
        return Err(format!(
            "staged::scan: the spec carries {} environments; the staged arm is \
             flat-only (use run_needle_cycles for a joint run)",
            spectral.spec.n_envs()
        ));
    }
    let current_merit = ctx.evaluate_merit(stack)?;
    if stack.films().is_empty() {
        return Ok(NeedleScan {
            current_merit,
            ..NeedleScan::default()
        });
    }

    let mut fold = spectral
        .folds
        .first()
        .cloned()
        .ok_or_else(|| "staged::scan: spectral.folds is empty".to_string())?;
    if cfg.refold_per_cycle && ctx.has_simulator() {
        let sim = ctx.simulate(stack)?;
        fold = build_needle_targets(
            &spectral.spec,
            &spectral.angles,
            &spectral.wavls,
            Some(&sim),
        )?;
    }
    let phi_gain_shift = fold.phi_gain_shift();

    let sa = stack.solver_arrays();
    let mut sites = build_scan_sites(stack.films(), stack.spans(), cfg.scan_step_nm);
    let mut rough_sites = Vec::new();
    if cfg.boundary_needles {
        for site in build_boundary_sites(stack.films(), stack.spans()) {
            if rough_boundary_start_nm(&sa, site.film_idx, 0.0, 0.0).is_some() {
                rough_sites.push(site);
            } else {
                sites.push(site);
            }
        }
    }

    let mut candidates = Vec::new();
    let mut sites_scanned = 0usize;

    for (material_index, material) in pool.iter().enumerate() {
        let material_sites = pool.sites_for(material_index, stack, &sites)?;
        if !material_sites.is_empty() {
            let input = NeedlePassInput {
                n_stack_cache: &sa.n_stack_cache,
                thicknesses: &sa.thicknesses,
                rough_types: &sa.rough_types,
                rough_vals: &sa.rough_vals,
                n_layers: sa.n_layers as usize,
                wavls: spectral.wavls.as_slice(),
                sin_theta: spectral.sin_theta.as_slice(),
                fold: &fold,
                needle_n_per_wav: material.nk.as_ref(),
                start_idx: 0,
                end_idx: (sa.n_layers - 1) as usize,
                incoherent_flags: &sa.incoherent_flags,
                calc_s: true,
                calc_p: true,
            };
            let result = needle_pass_scan(&input, &material_sites)?;
            sites_scanned += material_sites.len();
            for (site, &p) in result.sites.iter().zip(&result.p_profile) {
                if p >= 0.0 {
                    continue;
                }
                candidates.push(NeedleCandidate {
                    site: *site,
                    material_index,
                    material: material.material.clone(),
                    p: Some(p),
                    seed_nm: cfg.needle_seed_thickness_nm,
                    predicted_delta_merit: -(p + phi_gain_shift) * cfg.needle_seed_thickness_nm,
                });
            }
        }

        // Rough boundary sites are trial-scored: P cannot see the roughness
        // jump, so the candidate is the merit difference at its start
        // thickness. The trial is a stack clone; nothing is mutated.
        let rough_material_sites = pool.sites_for(material_index, stack, &rough_sites)?;
        for site in rough_material_sites {
            let d0 = rough_boundary_start_nm(
                &sa,
                site.film_idx,
                cfg.needle_seed_thickness_nm,
                stack.roughness_cap(),
            )
            .ok_or_else(|| "staged::scan: rough site without a rough boundary".to_string())?;
            let mut trial = stack.clone();
            trial
                .insert_film_at_boundary(site.film_idx, new_needle(stack.films(), material, d0))?;
            let gain = current_merit - ctx.evaluate_merit(&trial)?;
            sites_scanned += 1;
            if gain > 0.0 {
                candidates.push(NeedleCandidate {
                    site,
                    material_index,
                    material: material.material.clone(),
                    p: None,
                    seed_nm: d0,
                    predicted_delta_merit: -gain,
                });
            }
        }
    }

    candidates.sort_by(|a, b| {
        a.predicted_delta_merit
            .total_cmp(&b.predicted_delta_merit)
            .then_with(|| a.material.cmp(&b.material))
    });

    Ok(NeedleScan {
        candidates,
        current_merit,
        sites_scanned,
    })
}

/// A1 insert: put one scanned candidate into the stack, no re-optimization.
pub fn insert(
    stack: &mut DesignStack,
    pool: &MaterialPool,
    candidate: &NeedleCandidate,
) -> Result<Insertion, String> {
    let material = pool.get(candidate.material_index).ok_or_else(|| {
        format!(
            "staged::insert: material index {} is outside the pool",
            candidate.material_index
        )
    })?;
    let seed = new_needle(stack.films(), material, candidate.seed_nm);
    if candidate.site.boundary {
        stack.insert_film_at_boundary(candidate.site.film_idx, seed)?;
    } else {
        stack.insert_needle_seed(
            candidate.site.film_idx,
            candidate.site.depth_into_layer_nm,
            seed,
        )?;
    }
    Ok(Insertion {
        film_idx: candidate.site.film_idx,
        depth_into_layer_nm: candidate.site.depth_into_layer_nm,
        material: material.material.clone(),
        boundary: candidate.site.boundary,
        seed_nm: candidate.seed_nm,
    })
}

/// A1 polish: bounded LM over one free-set subset.
///
/// For `All` the original stack is refined in place. For a subset, a clone is
/// refined with the other films temporarily frozen; thicknesses are copied
/// back only when the topology is unchanged. A polish that would remove a
/// layer is refused - that is cleanup's job, not a step of the memetic
/// operator.
pub fn polish<C: DesignContext + ?Sized>(
    ctx: &mut C,
    stack: &mut DesignStack,
    scope: PolishScope,
) -> Result<PolishReport, String> {
    let n = stack.films().len();
    if n == 0 {
        return Err("staged::polish: the stack has no films".into());
    }
    let eligible: Vec<bool> = stack.films().iter().map(|film| film.optimize).collect();
    let selected = select_scope(&scope, stack, &eligible)?;
    let free_layers = selected.iter().filter(|&&s| s).count();
    if free_layers == 0 {
        return Err("staged::polish: the scope selects no free films".into());
    }

    let merit_before = ctx.evaluate_merit(stack)?;
    let merit_after = if matches!(scope, PolishScope::All) {
        ctx.refine_thicknesses(stack)?
    } else {
        let original: Vec<bool> = stack.films().iter().map(|film| film.optimize).collect();
        let mut trial = stack.clone();
        for (row, &free) in selected.iter().enumerate() {
            trial.set_row_optimize(row, free);
        }
        let trial_merit = ctx.refine_thicknesses(&mut trial)?;
        if trial.films().len() != n {
            return Err(format!(
                "staged::polish: the LM removed {} layer(s); use cleanup instead",
                n.saturating_sub(trial.films().len())
            ));
        }
        for (row, &optimize) in original.iter().enumerate() {
            trial.set_row_optimize(row, optimize);
        }
        for row in 0..n {
            stack.set_thickness(row, trial.films()[row].d_nm)?;
        }
        stack.refresh_profiles()?;
        trial_merit
    };

    Ok(PolishReport {
        scope,
        merit_before,
        merit_after,
        free_layers,
    })
}

fn select_scope(
    scope: &PolishScope,
    stack: &DesignStack,
    eligible: &[bool],
) -> Result<Vec<bool>, String> {
    let mut selected = vec![false; eligible.len()];
    match scope {
        PolishScope::All => selected.copy_from_slice(eligible),
        PolishScope::LastK(k) => {
            if *k == 0 {
                return Err("staged::polish: LastK(0) selects no films".into());
            }
            let mut taken = 0usize;
            for row in (0..eligible.len()).rev() {
                if eligible[row] {
                    selected[row] = true;
                    taken += 1;
                    if taken == *k {
                        break;
                    }
                }
            }
        }
        PolishScope::Named(names) => {
            if names.is_empty() {
                return Err("staged::polish: Named([]) selects no films".into());
            }
            for (row, film) in stack.films().iter().enumerate() {
                let label_matches = film
                    .label
                    .as_deref()
                    .is_some_and(|label| names.iter().any(|name| name == label));
                let name_matches = names.iter().any(|name| name == film.material.as_ref());
                if eligible[row] && (label_matches || name_matches) {
                    selected[row] = true;
                }
            }
        }
    }
    Ok(selected)
}

/// A1 cleanup: merge, prune, optionally re-optimize, and drain the clamp
/// report.
pub fn cleanup<C: DesignContext + ?Sized>(
    ctx: &mut C,
    stack: &mut DesignStack,
    options: CleanupOptions,
) -> Result<CleanupReport, String> {
    let result = cleanup_design(
        ctx,
        stack,
        options.min_thickness_nm,
        options.max_removals,
        options.reoptimize,
    )?;
    let clamp = ctx.take_clamp_report();
    Ok(CleanupReport { result, clamp })
}

// -- DesignStack conveniences ------------------------------------------------

impl DesignStack {
    /// [`scan`] on this stack.
    pub fn needle_scan<C: DesignContext + ?Sized>(
        &self,
        ctx: &C,
        spectral: &SpectralInputs,
        pool: &MaterialPool,
        cfg: &NeedleCycleConfig,
    ) -> Result<NeedleScan, String> {
        scan(ctx, self, spectral, pool, cfg)
    }

    /// [`insert`] one candidate into this stack.
    pub fn insert_needle_candidate(
        &mut self,
        pool: &MaterialPool,
        candidate: &NeedleCandidate,
    ) -> Result<Insertion, String> {
        insert(self, pool, candidate)
    }

    /// [`polish`] this stack.
    pub fn polish_needles<C: DesignContext + ?Sized>(
        &mut self,
        ctx: &mut C,
        scope: PolishScope,
    ) -> Result<PolishReport, String> {
        polish(ctx, self, scope)
    }

    /// [`cleanup`] this stack.
    pub fn cleanup_needles<C: DesignContext + ?Sized>(
        &mut self,
        ctx: &mut C,
        options: CleanupOptions,
    ) -> Result<CleanupReport, String> {
        cleanup(ctx, self, options)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lsq::lm::LmConfig;
    use crate::smatrix::angle::Radians;
    use crate::smatrix::synthesis::config::ThinLayerPolicy;
    use crate::smatrix::synthesis::evaluator::SmatrixContext;
    use crate::smatrix::synthesis::merit::{
        ConstraintKind, CurveId, MeritKey, MeritSpec, MeritTarget, SimTransform,
    };
    use crate::smatrix::synthesis::structure::LayerSpec;

    const WL: f64 = 550.0;

    fn constant(name: &str, n: f64, d_nm: f64) -> LayerSpec {
        LayerSpec::constant(name, n, 0.0, d_nm, 1)
    }

    fn stack(films: Vec<LayerSpec>) -> DesignStack {
        let mut air = constant("air", 1.0, 0.0);
        air.optimize = false;
        air.needle = false;
        let mut substrate = constant("substrate", 1.52, 0.0);
        substrate.optimize = false;
        substrate.needle = false;
        DesignStack::with_films(air, substrate, films).unwrap()
    }

    fn zero_spec() -> MeritSpec {
        let mut spec = MeritSpec::new();
        let key = spec.add_key(MeritKey {
            angle: Radians::ZERO,
            curve: CurveId::Rs,
        });
        spec.add_target(MeritTarget {
            key_idx: key as u32,
            env_idx: 0,
            wavelengths: vec![WL].into(),
            kind: ConstraintKind::Exact,
            transform: SimTransform::Linear,
            norm_factor: 1.0,
            normalized_targets: vec![0.0].into(),
            tolerances: vec![1e-3].into(),
            band: Vec::new().into(),
            phase: false,
            differential_passes: None,
            retardance: false,
            ratio: None,
            dispersion: None,
            stokes: None,
            point_weights: Vec::new().into(),
            integral: false,
            weight: 1.0,
            count_norm: None,
        })
        .unwrap();
        spec
    }

    fn ctx(spec: &MeritSpec) -> SmatrixContext {
        SmatrixContext {
            wavls: vec![WL],
            sin_theta: vec![0.0],
            spec: spec.clone(),
            clamp_min_nm: 2.0,
            clamp_max_nm: 1000.0,
            thin_layer_policy: ThinLayerPolicy::Remove,
            lm: LmConfig::default(),
            clamp_accumulator: ClampReport::default(),
            envs: None,
        }
    }

    fn spectral(spec: &MeritSpec) -> SpectralInputs {
        SpectralInputs::from_spec(spec, &[Radians::ZERO], &[WL]).unwrap()
    }

    fn cfg() -> NeedleCycleConfig {
        NeedleCycleConfig {
            scan_step_nm: 10.0,
            max_needles: 1,
            ..NeedleCycleConfig::default()
        }
    }

    fn pool() -> MaterialPool {
        MaterialPool::from(vec![constant("H", 2.2, 0.0), constant("L", 1.46, 0.0)])
    }

    #[test]
    fn scan_reports_sorted_candidates_and_sweeps_sites() {
        let spec = zero_spec();
        let c = ctx(&spec);
        let sp = spectral(&spec);
        let st = stack(vec![constant("H", 2.2, 80.0), constant("L", 1.46, 120.0)]);
        let scan = scan(&c, &st, &sp, &pool(), &cfg()).unwrap();

        assert!(scan.sites_scanned > 0, "{scan:?}");
        assert!(
            scan.candidates
                .windows(2)
                .all(|w| { w[0].predicted_delta_merit <= w[1].predicted_delta_merit })
        );
        assert!(scan.current_merit.is_finite());
    }

    #[test]
    fn insert_splits_the_host_without_reoptimizing() {
        let spec = zero_spec();
        let c = ctx(&spec);
        let sp = spectral(&spec);
        let st = stack(vec![constant("H", 2.2, 80.0), constant("L", 1.46, 120.0)]);
        let scan = scan(&c, &st, &sp, &pool(), &cfg()).unwrap();
        if let Some(best) = scan.best() {
            let mut trial = st.clone();
            let before = trial.films().len();
            let insertion = insert(&mut trial, &pool(), best).unwrap();
            let expected = before + if insertion.boundary { 1 } else { 2 };
            assert_eq!(trial.films().len(), expected);
        } else {
            let mut trial = st.clone();
            let candidate = NeedleCandidate {
                site: ScanSite {
                    film_idx: 0,
                    depth_into_layer_nm: 40.0,
                    z_nm: 40.0,
                    boundary: false,
                },
                material_index: 0,
                material: Arc::from("H"),
                p: Some(-1.0),
                seed_nm: 5.0,
                predicted_delta_merit: -5.0,
            };
            insert(&mut trial, &pool(), &candidate).unwrap();
            assert_eq!(trial.films().len(), st.films().len() + 2);
        }
    }

    #[test]
    fn polish_all_recovers_a_detuned_quarter_wave() {
        let spec = zero_spec();
        let mut c = ctx(&spec);
        let n = 1.52_f64.sqrt();
        let d = WL / (4.0 * n);
        let mut st = stack(vec![constant("L", n, 0.8 * d)]);
        let report = polish(&mut c, &mut st, PolishScope::All).unwrap();

        assert!(report.merit_after < report.merit_before);
        assert_eq!(report.free_layers, 1);
        assert!((st.films()[0].d_nm - d).abs() < 2.0);
    }

    #[test]
    fn polish_last_k_freezes_the_other_films() {
        let spec = zero_spec();
        let mut c = ctx(&spec);
        let mut st = stack(vec![constant("H", 2.2, 50.0), constant("L", 1.46, 80.0)]);
        let first_before = st.films()[0].d_nm;
        let report = polish(&mut c, &mut st, PolishScope::LastK(1)).unwrap();

        assert_eq!(report.free_layers, 1);
        assert_eq!(st.films()[0].d_nm, first_before);
    }

    #[test]
    fn cleanup_removes_a_thin_layer() {
        let spec = zero_spec();
        let mut c = ctx(&spec);
        let mut st = stack(vec![constant("H", 2.2, 0.4), constant("L", 1.46, 100.0)]);
        let before = st.films().len();
        let report = cleanup(
            &mut c,
            &mut st,
            CleanupOptions {
                min_thickness_nm: Some(2.0),
                max_removals: Some(1),
                reoptimize: false,
            },
        )
        .unwrap();

        assert_eq!(report.result.layers_before, before);
        assert!(report.result.layers_after < before, "{report:?}");
    }

    #[test]
    fn scan_of_an_empty_stack_is_empty() {
        let spec = zero_spec();
        let c = ctx(&spec);
        let sp = spectral(&spec);
        let st = stack(Vec::new());
        let scan = scan(&c, &st, &sp, &pool(), &cfg()).unwrap();
        assert!(scan.is_empty());
        assert_eq!(scan.sites_scanned, 0);
    }
}
