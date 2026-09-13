// SPDX-License-Identifier: LGPL-3.0-or-later
//! Navette -- Rust Rewrite of Numba-optimized thin-film optical solver
//!
//! synthesis::pipeline — NeedlePipeline macro-loop.
//!
//! Verbatim port of needle_pipeline.py section 6. Each macro-cycle:
//!   1. Needle pass   (cycle.rs: analytic insertions + optimize)
//!   2. Cleanup       (optional; threshold defaults to clamp_min;
//!                     post-cleanup clamp sweep)
//!   3. Inflate       (optional; clamps BEFORE and after re-optimize —
//!                     the ClampedNeedleSynthesizer overrides)
//!   4. Stagnation    record mf_end → check (divergence → oscillation → plateau)
//!
//! Budget checks run pre-flight, after the needle phase, and post-cycle.
//! The loop ends with a final optimize + clamp sweep + evaluation.
//! User abort = callback returning Err (the KeyboardInterrupt analog).

// clippy::doc_overindented_list_items: the list above is a hand-aligned
// table (name / description / parenthetical). Clippy's two-space
// continuation rule would break the columns, and the columns are what
// make the block readable in source -- which is where it is read.
#![allow(clippy::doc_overindented_list_items)]

use crate::smatrix::synthesis::cleanup::{CleanupResult, cleanup_design};
use crate::smatrix::synthesis::config::{PipelineConfig, TerminationReason, ThinLayerPolicy};
use crate::smatrix::synthesis::context::DesignContext;
use crate::smatrix::synthesis::cycle::{
    ContrastMap, NeedleCycleConfig, NeedleCycleResult, run_needle_cycles,
};
use crate::smatrix::synthesis::inflate::{InflateResult, inflate_design};
use crate::smatrix::synthesis::merit::MeritSpec;
use crate::smatrix::synthesis::needle_pass::{NeedleTargets, build_needle_targets};
use crate::smatrix::synthesis::stagnation::StagnationDetector;
use crate::smatrix::synthesis::structure::{ClampReport, DesignStack};

/// Record of one macro-cycle — mirrors Python `PipelinePhaseResult`.
#[derive(Clone, Debug)]
pub struct PipelinePhaseResult {
    pub macro_cycle: usize,
    pub mf_after_needle: f64,
    pub mf_after_cleanup: Option<f64>,
    pub mf_after_inflate: Option<f64>,
    pub mf_end: f64,
    pub layer_count: usize,
    pub total_thickness_nm: f64,
    pub needle_results: Vec<NeedleCycleResult>,
    pub cleanup_result: Option<CleanupResult>,
    pub inflate_result: Option<InflateResult>,
    /// F0.2: every clamp pass this phase ran - the pipeline's own
    /// post-cleanup / post-inflate sweeps plus the clamps inside each
    /// `optimize_thicknesses`, drained from the context at record time.
    /// `None` when nothing was removed and nothing was capped (B4), so a
    /// no-span run's phase dicts are byte-identical.
    pub clamp_report: Option<ClampReport>,
}

/// Full output of a pipeline run — mirrors Python `PipelineResult`.
#[derive(Clone, Debug)]
pub struct PipelineResult {
    pub phases: Vec<PipelinePhaseResult>,
    pub termination: TerminationReason,
    pub final_mf: f64,
    pub final_layer_count: usize,
    pub final_total_thickness_nm: f64,
    pub stagnation_detail: Option<String>,
    /// F0.2: the final clamp sweep (after the last optimization). `None`
    /// when it removed and capped nothing (B4).
    pub final_clamp_report: Option<ClampReport>,
}

/// Continuous iterative needle synthesis pipeline.
pub struct NeedlePipeline {
    pub stack: DesignStack,
    pub cfg: PipelineConfig,
    pub needle_cfg: NeedleCycleConfig,
    pub contrast: ContrastMap,
    /// Fixed spectral problem definition for the needle pass.
    pub spectral: SpectralInputs,
    detector: StagnationDetector,
}

impl NeedlePipeline {
    pub fn new(
        stack: DesignStack,
        spectral: SpectralInputs,
        cfg: PipelineConfig,
        needle_cfg: NeedleCycleConfig,
        contrast: ContrastMap,
    ) -> Result<Self, String> {
        let cfg = cfg.validated()?;
        // F0.2 (licence item 2) narrowed by F1.6: a PINNED profile above
        // the manufacturing ceiling is an authoring error - the user
        // asked for a film the machine cannot make, and silently
        // squeezing it in answers a question nobody asked. A SCALABLE
        // span (every bulk row optimize-flagged - F1.6's optimized
        // profiled carriers) instead gets the ceiling as an LM bound.
        for sp in stack.spans() {
            if sp.end - sp.start <= 1 || stack.span_is_scalable(sp) {
                continue;
            }
            let d: f64 = stack.films()[sp.start..sp.end].iter().map(|l| l.d_nm).sum();
            if d > cfg.clamp_max_nm {
                return Err(format!(
                    "NeedlePipeline: span '{}' is {:.1} nm thick, above the                      clamp_max_nm ceiling of {:.1} nm - refusing rather than                      rescaling a pinned profile",
                    stack.films()[sp.start].material,
                    d,
                    cfg.clamp_max_nm
                ));
            }
        }
        // F0.3: `ClampUpAlways` silently disables layer elimination, so
        // needle runs only ever grow - a bad seed parks at the floor
        // permanently instead of being rejected by the optimizer
        // shrinking it through the floor. Refuse the combination, name
        // both settings, and point at the variant that works.
        if cfg.thin_layer_policy == ThinLayerPolicy::ClampUpAlways && cfg.needles_per_cycle > 0 {
            return Err(format!(
                "NeedlePipeline: thin_layer_policy 'clamp_up_always' conflicts with \
                 needles_per_cycle = {} - with the floor as a hard LM bound nothing \
                 can ever be eliminated, so a rejected seed parks at the floor and \
                 the layer count only ever grows. Use 'clamp_up_final' (the \
                 documented default for needle runs) or set needles_per_cycle = 0.",
                cfg.needles_per_cycle
            ));
        }
        let detector = StagnationDetector::new(
            cfg.stagnation_window,
            cfg.stagnation_gradient_tol,
            cfg.stagnation_oscillation_ratio,
            cfg.stagnation_divergence_count,
        );
        Ok(NeedlePipeline {
            stack,
            spectral,
            cfg,
            needle_cfg,
            contrast,
            detector,
        })
    }

    fn check_budgets<C: DesignContext + ?Sized>(
        &self,
        ctx: &C,
    ) -> Result<Option<TerminationReason>, String> {
        // F0.2 (U5): the budget is a manufacturability limit - it means
        // physical layers, and one graded film is one physical layer, so
        // it counts spans, not solver rows.
        if self.stack.spans().len() >= self.cfg.max_film_layers {
            return Ok(Some(TerminationReason::LayerBudgetReached));
        }
        let total: f64 = self.stack.films().iter().map(|l| l.d_nm).sum();
        if total >= self.cfg.max_total_thickness_nm {
            return Ok(Some(TerminationReason::ThicknessBudgetReached));
        }
        if self.cfg.merit_target > 0.0 {
            let mf = ctx.evaluate_merit(&self.stack)?;
            if mf <= self.cfg.merit_target {
                return Ok(Some(TerminationReason::MeritTargetReached));
            }
        }
        Ok(None)
    }

    /// Execute the pipeline.
    ///
    /// `callback(macro_cycle, &phase, &detector)` runs after each completed
    /// cycle; an `Err` return aborts the run as [`TerminationReason::UserAbort`].
    pub fn run<C: DesignContext + ?Sized>(
        &mut self,
        ctx: &mut C,
        mut callback: impl FnMut(usize, &PipelinePhaseResult, &StagnationDetector) -> Result<(), String>,
    ) -> Result<PipelineResult, String> {
        let mut phases: Vec<PipelinePhaseResult> = Vec::new();
        self.detector.reset();

        // Per-cycle needle budget.
        let mut needle_cfg = self.needle_cfg.clone();
        needle_cfg.max_needles = self.cfg.needles_per_cycle;

        let mut termination = TerminationReason::MaxIterationsReached;
        let mut stag_detail: Option<String> = None;
        let mut user_abort = false;

        'main: for cycle_i in 1..=self.cfg.max_macro_cycles {
            // -- Pre-flight budget check --
            if let Some(reason) = self.check_budgets(ctx)? {
                termination = reason;
                break;
            }

            // F0.2: this phase's clamp aggregate - the context's
            // accumulator (clamps inside every optimize_thicknesses) plus
            // the pipeline's own sweeps below, merged at record time.
            let mut clamp_report = ctx.take_clamp_report().unwrap_or_default();

            // ── Phase 1: Needle pass ──
            let needle_results = run_needle_cycles(
                ctx,
                &mut self.stack,
                &self.spectral,
                &self.contrast,
                &needle_cfg,
            )?;
            let mf_needle = ctx.evaluate_merit(&self.stack)?;

            // Budget check after needle.
            if let Some(reason) = self.check_budgets(ctx)? {
                let total: f64 = self.stack.films().iter().map(|l| l.d_nm).sum();
                let phase = PipelinePhaseResult {
                    macro_cycle: cycle_i,
                    mf_after_needle: mf_needle,
                    mf_after_cleanup: None,
                    mf_after_inflate: None,
                    mf_end: mf_needle,
                    layer_count: self.stack.spans().len(),
                    total_thickness_nm: total,
                    needle_results,
                    cleanup_result: None,
                    inflate_result: None,
                    clamp_report: ctx.take_clamp_report(),
                };
                phases.push(phase);
                let phase = phases.last().unwrap();
                self.detector.record(phase.mf_end);
                // KeyboardInterrupt analog: callback failure aborts the loop.
                if callback(cycle_i, phase, &self.detector).is_err() {
                    user_abort = true;
                }
                termination = reason;
                break 'main;
            }

            // ── Phase 2: Cleanup (optional) ──
            let cleanup_result = if self.cfg.enable_cleanup {
                let r = cleanup_design(
                    ctx,
                    &mut self.stack,
                    self.cfg.cleanup_min_nm,
                    self.cfg.cleanup_max_removals,
                    true,
                )?;
                // Post-cleanup clamp (Clamped override semantics). F0.2:
                // the report joins the phase's aggregate. F0.3: during
                // the run the floor clamps up only under ClampUpAlways.
                // F1.6: scalable spans follow the same pass structure.
                let rep = self.stack.clamp_all_policy(
                    self.cfg.clamp_min_nm,
                    self.cfg.clamp_max_nm,
                    self.cfg.thin_layer_policy,
                    false,
                )?;
                if !rep.is_empty() {
                    clamp_report.merge(rep);
                }
                Some(r)
            } else {
                None
            };
            let mf_cleanup = cleanup_result.as_ref().map(|r| r.merit_after);

            // ── Phase 3: Inflate (optional) ──
            let inflate_result = if self.cfg.enable_inflate {
                let r = inflate_design(
                    ctx,
                    &mut self.stack,
                    &self.spectral.wavls,
                    self.cfg.inflate_addon_qwot,
                    self.cfg.inflate_reference_wl,
                    self.cfg.inflate_max_layers,
                    true,
                )?;
                // Clamp BEFORE re-optimize happened in Python before the call
                // ordering above; enforce the AFTER clamp here too. F0.2:
                // the report joins the phase's aggregate. F0.3: during
                // the run the floor clamps up only under ClampUpAlways.
                let rep = self.stack.clamp_all_policy(
                    self.cfg.clamp_min_nm,
                    self.cfg.clamp_max_nm,
                    self.cfg.thin_layer_policy,
                    false,
                )?;
                if !rep.is_empty() {
                    clamp_report.merge(rep);
                }
                Some(r)
            } else {
                None
            };
            let mf_inflate = inflate_result.as_ref().map(|r| r.merit_after);

            let mf_cleanup_s = mf_cleanup;
            let mf_inflate_s = mf_inflate;

            // ── Record phase ──
            let total: f64 = self.stack.films().iter().map(|l| l.d_nm).sum();
            let mf_end = mf_inflate.or(mf_cleanup).unwrap_or(mf_needle);
            clamp_report.merge(ctx.take_clamp_report().unwrap_or_default());
            phases.push(PipelinePhaseResult {
                macro_cycle: cycle_i,
                mf_after_needle: mf_needle,
                mf_after_cleanup: mf_cleanup_s,
                mf_after_inflate: mf_inflate_s,
                mf_end,
                layer_count: self.stack.spans().len(),
                total_thickness_nm: total,
                needle_results,
                cleanup_result,
                inflate_result,
                clamp_report: if clamp_report.is_empty() {
                    None
                } else {
                    Some(clamp_report)
                },
            });

            // ── Stagnation check ──
            self.detector.record(mf_end);
            {
                let phase = phases.last().unwrap();
                if callback(cycle_i, phase, &self.detector).is_err() {
                    // KeyboardInterrupt analog: skip straight to finalization.
                    user_abort = true;
                    break 'main;
                }
            }

            if let Some(stag) = self.detector.check() {
                termination = stag;
                stag_detail = Some(self.detector.summary());
                break 'main;
            }

            // Post-cycle budget check.
            if let Some(reason) = self.check_budgets(ctx)? {
                termination = reason;
                break 'main;
            }
        }

        // Exception-style abort overrides any other reason (Python catches
        // KeyboardInterrupt around the whole loop).
        if user_abort {
            termination = TerminationReason::UserAbort;
            stag_detail = None;
        }

        // ── Final optimisation + clamp sweep ──
        ctx.optimize_thicknesses(&mut self.stack)?;
        // F0.3: the final pass substitutes clamping for removal under
        // both clamp-up policies; the reported merit is evaluated AFTER
        // this sweep (the B7 ordering pin asserts it). F1.6: scalable
        // spans ride the same pass structure - scaled up/capped instead
        // of removed/refused.
        let rep_final = self.stack.clamp_all_policy(
            self.cfg.clamp_min_nm,
            self.cfg.clamp_max_nm,
            self.cfg.thin_layer_policy,
            true,
        )?;
        let final_mf = ctx.evaluate_merit(&self.stack)?;

        // F0.2: the final sweep's own report plus any clamps the final
        // optimization ran, surfaced on the run result (B4: absent when
        // empty).
        let mut final_clamp_report = ctx.take_clamp_report().unwrap_or_default();
        final_clamp_report.merge(rep_final);
        let final_clamp_report = if final_clamp_report.is_empty() {
            None
        } else {
            Some(final_clamp_report)
        };

        Ok(PipelineResult {
            phases,
            termination,
            final_mf,
            final_layer_count: self.stack.spans().len(),
            final_total_thickness_nm: self.stack.films().iter().map(|l| l.d_nm).sum(),
            stagnation_detail: stag_detail,
            final_clamp_report,
        })
    }
}

/// Spectral problem definition needed by the needle pass: the fixed grid
/// plus the full folded demand set (all quantities — R/T/A, back-incidence
/// siblings, per-channel phase pairs). Folded once at construction
/// (conservative form: no operating-point sim, so one-sided/banded kinds
/// take their conservative arm and masking resolves inside the scan's
/// merit evaluations, exactly as the standalone fold documents).
#[derive(Clone, Debug)]
pub struct SpectralInputs {
    pub wavls: Vec<f64>,
    /// Sines of incidence angles.
    pub sin_theta: Vec<f64>,
    /// Conservative fold of `spec` on (`angles_deg`, `wavls`) — the
    /// starting fold; re-folded against the live sim each needle cycle
    /// (see `run_needle_cycles`).
    pub fold: NeedleTargets,
    /// Kept for per-cycle re-folds (spec + degree-convention angles).
    pub spec: MeritSpec,
    pub angles_deg: Vec<f64>,
}

impl SpectralInputs {
    /// Build from a merit spec: `angles_deg` matches the spec's key
    /// convention (degrees, as produced by the Python converter).
    pub fn from_spec(spec: &MeritSpec, angles_deg: &[f64], wavls: &[f64]) -> Result<Self, String> {
        let fold = build_needle_targets(spec, angles_deg, wavls, None)?;
        Ok(SpectralInputs {
            wavls: wavls.to_vec(),
            sin_theta: angles_deg
                .iter()
                .map(|a| (a * std::f64::consts::PI / 180.0).sin())
                .collect(),
            fold,
            spec: spec.clone(),
            angles_deg: angles_deg.to_vec(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smatrix::synthesis::merit::SimCurves;
    use crate::smatrix::synthesis::structure::LayerSpec;

    const NW: usize = 3;

    fn dummy_spectral() -> SpectralInputs {
        let zero = || (vec![0.0f64; NW], vec![0.0f64; NW]);
        SpectralInputs {
            wavls: vec![500.0; NW],
            sin_theta: vec![0.0],
            fold: NeedleTargets {
                r: (vec![0.0; NW], vec![1.0; NW]),
                t: zero(),
                a: zero(),
                rb: zero(),
                tb: zero(),
                ab: zero(),
                phi: [zero(), zero(), zero(), zero()],
                phi_gain_shift: [0.0; 4],
                grad_r: vec![0.0f64; NW],
                grad_t: vec![0.0f64; NW],
            },
            spec: MeritSpec::new(),
            angles_deg: vec![0.0],
        }
    }

    fn air() -> LayerSpec {
        let mut l = LayerSpec::constant("air", 1.0, 0.0, 0.0, NW);
        l.optimize = false;
        l
    }

    fn sub() -> LayerSpec {
        let mut l = LayerSpec::constant("sub", 1.52, 0.0, 0.0, NW);
        l.optimize = false;
        l
    }

    /// No-op context: MF constant (2.5), optimize does nothing.
    struct FlatCtx;
    impl DesignContext for FlatCtx {
        fn evaluate_merit(&self, _s: &DesignStack) -> Result<f64, String> {
            Ok(2.5)
        }
        fn simulate(&self, _s: &DesignStack) -> Result<SimCurves, String> {
            Err("mock context has no simulator".into())
        }
        fn optimize_thicknesses(&mut self, s: &mut DesignStack) -> Result<f64, String> {
            self.evaluate_merit(s)
        }
    }

    use crate::smatrix::synthesis::merit::MeritSpec;
    use crate::smatrix::synthesis::needle_pass::NeedleTargets;

    fn pipeline(cfg_over: impl FnOnce(&mut PipelineConfig)) -> NeedlePipeline {
        let mut cfg = PipelineConfig::default();
        cfg_over(&mut cfg);
        let stack = DesignStack::with_films(
            air(),
            sub(),
            vec![LayerSpec::constant("H", 2.35, 0.0, 100.0, NW)],
        )
        .unwrap();
        NeedlePipeline::new(
            stack,
            dummy_spectral(),
            cfg,
            NeedleCycleConfig::default(),
            ContrastMap::new(), // no contrast → needle pass inserts nothing
        )
        .unwrap()
    }

    #[test]
    fn max_macro_cycles_is_default_termination() {
        let mut p = pipeline(|c| c.max_macro_cycles = 2);
        let res = p.run(&mut FlatCtx, |_, _, _| Ok(())).unwrap();
        assert_eq!(res.termination, TerminationReason::MaxIterationsReached);
        assert_eq!(res.phases.len(), 2);
        assert!((res.final_mf - 2.5).abs() < 1e-12);
    }

    #[test]
    fn layer_budget_fires_preflight() {
        let mut p = pipeline(|c| {
            c.max_macro_cycles = 5;
            c.max_film_layers = 1; // 1 film present ≥ 1
        });
        let res = p.run(&mut FlatCtx, |_, _, _| Ok(())).unwrap();
        assert_eq!(res.termination, TerminationReason::LayerBudgetReached);
        assert!(res.phases.is_empty()); // fired BEFORE cycle 1
    }

    #[test]
    fn merit_target_fires_preflight() {
        let mut p = pipeline(|c| {
            c.max_macro_cycles = 5;
            c.merit_target = 3.0; // initial MF 2.5 ≤ 3 → stop immediately
        });
        let res = p.run(&mut FlatCtx, |_, _, _| Ok(())).unwrap();
        assert_eq!(res.termination, TerminationReason::MeritTargetReached);
        assert!(res.phases.is_empty());
    }

    #[test]
    fn thickness_budget_fires() {
        let mut p = pipeline(|c| {
            c.max_total_thickness_nm = 50.0; // stack is 100 nm thick already
        });
        let res = p.run(&mut FlatCtx, |_, _, _| Ok(())).unwrap();
        assert_eq!(res.termination, TerminationReason::ThicknessBudgetReached);
    }

    #[test]
    fn callback_error_aborts_as_user_abort() {
        let mut p = pipeline(|c| {
            // Need one completed phase before callback fires post-cycle:
            // disable cleanup/inflate so phases complete; but pre-flight
            // budgets must pass: defaults ok.
            c.enable_cleanup = false;
            c.stagnation_window = usize::MAX; // never fire stagnation first
        });
        let res = p.run(&mut FlatCtx, |_, _, _| Err("stop!".into())).unwrap();
        assert_eq!(res.termination, TerminationReason::UserAbort);
    }

    #[test]
    fn plateau_terminates_flat_trajectory() {
        let mut p = pipeline(|c| {
            c.enable_cleanup = false;
            c.stagnation_window = 2;
            c.max_macro_cycles = 10;
        });
        let res = p.run(&mut FlatCtx, |_, _, _| Ok(())).unwrap();
        assert_eq!(res.termination, TerminationReason::StagnationPlateau);
        assert!(res.stagnation_detail.is_some());
        // Window of 2 samples recorded then plateau detected at cycle 2.
        assert_eq!(res.phases.len(), 2);
    }

    #[test]
    fn final_optimize_and_clamp_run_after_loop() {
        // Context whose optimizer caps thicknesses via recorded calls.
        struct CountingCtx {
            opt_calls: usize,
        }
        impl DesignContext for CountingCtx {
            fn evaluate_merit(&self, _: &DesignStack) -> Result<f64, String> {
                Ok(1.0)
            }
            fn simulate(&self, _: &DesignStack) -> Result<SimCurves, String> {
                Err("mock context has no simulator".into())
            }
            fn optimize_thicknesses(&mut self, s: &mut DesignStack) -> Result<f64, String> {
                self.opt_calls += 1;
                s.set_thickness(0, 5000.0)?;
                self.evaluate_merit(s)
            }
        }
        // clamp_max small → final clamp sweep must cap the 5000 nm film.
        let mut p = pipeline(|c| {
            c.clamp_max_nm = 800.0;
            c.max_macro_cycles = 1;
            c.enable_cleanup = false;
            c.stagnation_window = usize::MAX;
        });
        let mut ctx = CountingCtx { opt_calls: 0 };
        let res = p.run(&mut ctx, |_, _, _| Ok(())).unwrap();
        // needle pass optimize (initial) + per-cycle + FINAL = counted
        assert!(ctx.opt_calls >= 2);
        assert!(res.final_total_thickness_nm <= 800.0);
    }

    // ------------------------------------------------------------------
    // F0.2 - the layer budget counts spans; the ceiling refuses at the door
    // ------------------------------------------------------------------

    /// Licence item 3, measured: one 1000 nm graded film at delta = 0.5
    /// expands to 57 rows. Before F0.2, `max_film_layers = 40` terminated
    /// the run on the pre-flight of cycle 1 before any work; the budget is
    /// a manufacturability limit - one graded film is ONE physical layer.
    #[test]
    fn f02_budget_counts_spans_not_rows() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = std::collections::HashMap::new();
        nk.insert(
            std::sync::Arc::from("TiO2"),
            vec![num_complex::Complex64::new(2.35, 0.0); NW],
        );
        let graded = vec![crate::structure::Layer {
            inhomogen: true,
            inh_delta: 0.5,
            ..crate::structure::Layer::film(1000.0, "TiO2")
        }];
        let bg: std::collections::HashSet<String> = ["TiO2".to_string()].into_iter().collect();
        let (stack, _) = DesignStack::from_design(
            air(),
            sub(),
            &graded,
            &nk,
            &std::collections::HashMap::new(),
            &wl,
            &bg,
        )
        .unwrap();
        assert_eq!(stack.films().len(), 57, "the plan's measured row count");
        assert_eq!(stack.spans().len(), 1);
        let cfg = PipelineConfig {
            max_macro_cycles: 1,
            max_film_layers: 40,  // < 57 rows, > 1 span: the old trip
            clamp_max_nm: 1500.0, // the span total must pass the ceiling
            max_total_thickness_nm: 10_000.0,
            stagnation_window: usize::MAX,
            ..Default::default()
        };
        let mut p = NeedlePipeline::new(
            stack,
            dummy_spectral(),
            cfg,
            NeedleCycleConfig::default(),
            ContrastMap::new(),
        )
        .unwrap();
        let res = p.run(&mut FlatCtx, |_, _, _| Ok(())).unwrap();
        // The run PROCEEDS (the old binary terminated before cycle 1).
        assert_eq!(res.termination, TerminationReason::MaxIterationsReached);
        assert_eq!(res.phases.len(), 1);
        // Licence item 4: layer_count reads the physical-layer count.
        assert_eq!(res.phases[0].layer_count, 1);
        assert_eq!(res.final_layer_count, 1);
    }

    /// Licence item 2 at the door: a pinned profile above the ceiling is
    /// an authoring error - refused at `NeedlePipeline::new`, naming the
    /// span material and both numbers. Same span, larger ceiling: builds.
    #[test]
    fn f02_new_refuses_pinned_span_above_the_ceiling() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = std::collections::HashMap::new();
        nk.insert(
            std::sync::Arc::from("TiO2"),
            vec![num_complex::Complex64::new(2.35, 0.0); NW],
        );
        let graded = vec![crate::structure::Layer {
            inhomogen: true,
            inh_delta: 0.5,
            ..crate::structure::Layer::film(1000.0, "TiO2")
        }];
        let bg: std::collections::HashSet<String> = ["TiO2".to_string()].into_iter().collect();
        let (stack, _) = DesignStack::from_design(
            air(),
            sub(),
            &graded,
            &nk,
            &std::collections::HashMap::new(),
            &wl,
            &bg,
        )
        .unwrap();
        let cfg = PipelineConfig {
            clamp_max_nm: 300.0,
            ..Default::default()
        };
        let err = match NeedlePipeline::new(
            stack,
            dummy_spectral(),
            cfg,
            NeedleCycleConfig::default(),
            ContrastMap::new(),
        ) {
            Err(e) => e,
            Ok(_) => panic!("expected the ceiling refusal"),
        };
        assert!(err.contains("'TiO2'"), "{err}");
        assert!(err.contains("1000.0"), "{err}");
        assert!(err.contains("300.0"), "{err}");
        // Control: the same span under a 1500 nm ceiling constructs.
        let (stack2, _) = DesignStack::from_design(
            air(),
            sub(),
            &graded,
            &nk,
            &std::collections::HashMap::new(),
            &wl,
            &bg,
        )
        .unwrap();
        let cfg2 = PipelineConfig {
            clamp_max_nm: 1500.0,
            max_total_thickness_nm: 10_000.0,
            ..Default::default()
        };
        assert!(
            NeedlePipeline::new(
                stack2,
                dummy_spectral(),
                cfg2,
                NeedleCycleConfig::default(),
                ContrastMap::new(),
            )
            .is_ok()
        );
    }

    /// F1.6: the same span, marked optimize=true (NOT background - bg
    /// empty), is a SCALABLE span: it no longer refuses at the door, it
    /// gets the ceiling as an LM bound. The pinned twin above is
    /// unchanged.
    #[test]
    fn f16_new_builds_scalable_span_above_the_ceiling() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = std::collections::HashMap::new();
        nk.insert(
            std::sync::Arc::from("TiO2"),
            vec![num_complex::Complex64::new(2.35, 0.0); NW],
        );
        let graded = vec![crate::structure::Layer {
            inhomogen: true,
            inh_delta: 0.5,
            ..crate::structure::Layer::film(1000.0, "TiO2")
        }];
        let (stack, warns) = DesignStack::from_design(
            air(),
            sub(),
            &graded,
            &nk,
            &std::collections::HashMap::new(),
            &wl,
            &std::collections::HashSet::new(), // NOT background
        )
        .unwrap();
        assert!(
            warns.is_empty(),
            "an optimized profiled carrier must not homogenize: {warns:?}"
        );
        assert_eq!(stack.films().len(), 57, "the profile kept all 57 rows");
        assert!(stack.span_is_scalable(stack.spans().first().unwrap()));
        let cfg = PipelineConfig {
            clamp_max_nm: 300.0,
            ..Default::default()
        };
        assert!(
            NeedlePipeline::new(
                stack,
                dummy_spectral(),
                cfg,
                NeedleCycleConfig::default(),
                ContrastMap::new(),
            )
            .is_ok(),
            "a scalable span gets the ceiling as a bound, not a refusal"
        );
    }

    /// Licence item 6 (B4): the phase's clamp report is present when the
    /// floor removed something, absent otherwise - a no-span run's phase
    /// dict is byte-identical. The film is PINNED (optimize = false):
    /// cleanup's flag-guarded removal skips it, so the clamp is what
    /// takes it, and the report is what says so.
    #[test]
    fn f02_clamp_report_present_only_when_something_happened() {
        // Floor above the film: one-row span removed and NAMED.
        let mut pinned = LayerSpec::constant("H", 2.35, 0.0, 100.0, NW);
        pinned.optimize = false;
        pinned.needle = false;
        let cfg = PipelineConfig {
            max_macro_cycles: 1,
            enable_inflate: false,
            clamp_min_nm: 200.0,
            stagnation_window: usize::MAX,
            ..Default::default()
        };
        let mut p = NeedlePipeline::new(
            DesignStack::with_films(air(), sub(), vec![pinned]).unwrap(),
            dummy_spectral(),
            cfg,
            NeedleCycleConfig::default(),
            ContrastMap::new(),
        )
        .unwrap();
        let res = p.run(&mut FlatCtx, |_, _, _| Ok(())).unwrap();
        let rep = res.phases[0].clamp_report.as_ref().unwrap();
        assert_eq!(rep.spans_removed, vec!["H (100.0 nm)".to_string()]);
        assert_eq!(rep.rows_removed, 1);
        // After the phase removed the only film, the final sweep does nothing.
        assert!(res.final_clamp_report.is_none());

        // Default floor: nothing removed, nothing capped - report absent.
        let mut p = pipeline(|c| {
            c.max_macro_cycles = 1;
            c.enable_inflate = false;
            c.stagnation_window = usize::MAX;
        });
        let res = p.run(&mut FlatCtx, |_, _, _| Ok(())).unwrap();
        assert!(res.phases[0].clamp_report.is_none());
        assert!(res.final_clamp_report.is_none());
    }

    // ------------------------------------------------------------------
    // F0.3 - ThinLayerPolicy
    // ------------------------------------------------------------------

    /// `ClampUpFinal` twin: the search runs exactly as today (elimination
    /// and all); only the FINAL clamp pass sets a surviving sub-minimum
    /// film to `clamp_min_nm` instead of removing it. The mock optimizer
    /// drives the film to 0.8 nm in the final solve - under `Remove` the
    /// returned stack has no film; under `ClampUpFinal` it has one at
    /// exactly the floor, and the reported merit describes the stack the
    /// user receives (B7: the pin on the ordering the tree already keeps).
    #[test]
    fn f03_clamp_up_final_lands_the_film_on_the_floor() {
        struct ThinningCtx;
        impl DesignContext for ThinningCtx {
            fn evaluate_merit(&self, _s: &DesignStack) -> Result<f64, String> {
                Ok(1.0)
            }
            fn simulate(&self, _s: &DesignStack) -> Result<SimCurves, String> {
                Err("mock context has no simulator".into())
            }
            fn optimize_thicknesses(&mut self, s: &mut DesignStack) -> Result<f64, String> {
                s.set_thickness(0, 0.8)?;
                self.evaluate_merit(s)
            }
        }
        for (policy, want_count, want_d) in [
            ("remove", 0usize, None),
            ("clamp_up_final", 1usize, Some(2.0)),
        ] {
            let cfg = PipelineConfig {
                max_macro_cycles: 1,
                enable_cleanup: false,
                needles_per_cycle: 0,
                stagnation_window: usize::MAX,
                thin_layer_policy: ThinLayerPolicy::parse(policy).unwrap(),
                ..Default::default()
            };
            let mut p = NeedlePipeline::new(
                DesignStack::with_films(
                    air(),
                    sub(),
                    vec![LayerSpec::constant("H", 2.35, 0.0, 100.0, NW)],
                )
                .unwrap(),
                dummy_spectral(),
                cfg,
                NeedleCycleConfig::default(),
                ContrastMap::new(),
            )
            .unwrap();
            let mut ctx = ThinningCtx;
            let res = p.run(&mut ctx, |_, _, _| Ok(())).unwrap();
            assert_eq!(res.final_layer_count, want_count, "policy {policy}");
            if let Some(d) = want_d {
                assert!(
                    (p.stack.films()[0].d_nm - d).abs() < 1e-12,
                    "policy {policy}: film at {}, wanted {d}",
                    p.stack.films()[0].d_nm
                );
            }
            // B7 pin: the reported merit describes the RETURNED stack.
            if want_count == 1 {
                let fresh = match ctx.evaluate_merit(&p.stack) {
                    Ok(v) => v,
                    Err(_) => panic!("mock evaluation failed"),
                };
                assert_eq!(res.final_mf, fresh, "final_mf describes the returned stack");
            }
        }
    }

    /// Conflict-refusal twin: `ClampUpAlways` + needles = layer count only
    /// ever grows. The refusal names both settings and points at
    /// `ClampUpFinal`.
    #[test]
    fn f03_clamp_up_always_is_refused_alongside_needles() {
        let cfg = PipelineConfig {
            thin_layer_policy: ThinLayerPolicy::ClampUpAlways,
            needles_per_cycle: 3,
            ..Default::default()
        };
        let err = match NeedlePipeline::new(
            DesignStack::with_films(
                air(),
                sub(),
                vec![LayerSpec::constant("H", 2.35, 0.0, 100.0, NW)],
            )
            .unwrap(),
            dummy_spectral(),
            cfg,
            NeedleCycleConfig::default(),
            ContrastMap::new(),
        ) {
            Err(e) => e,
            Ok(_) => panic!("expected the conflict refusal"),
        };
        assert!(err.contains("clamp_up_always"), "{err}");
        assert!(err.contains("needles_per_cycle = 3"), "{err}");
        assert!(err.contains("clamp_up_final"), "{err}");
        // The same policy with needles off constructs.
        let cfg = PipelineConfig {
            thin_layer_policy: ThinLayerPolicy::ClampUpAlways,
            needles_per_cycle: 0,
            ..Default::default()
        };
        assert!(
            NeedlePipeline::new(
                DesignStack::with_films(
                    air(),
                    sub(),
                    vec![LayerSpec::constant("H", 2.35, 0.0, 100.0, NW)],
                )
                .unwrap(),
                dummy_spectral(),
                cfg,
                NeedleCycleConfig::default(),
                ContrastMap::new(),
            )
            .map(|_| ())
            .is_ok()
        );
    }

    /// Span deferral twin: an under-thickness PINNED graded span (the
    /// carrier is background: optimize forced false) is removed whole
    /// with the F0.2 report under EVERY policy. F1.6 landed the scale-up
    /// for SCALABLE spans (f16_under_thickness_scalable_span_is_scaled_to
    /// _the_floor); a pinned span's profile cannot be re-parameterized by
    /// the optimizer, so removal is still the honest repair here.
    #[test]
    fn f03_thin_graded_span_is_removed_whole_under_clamp_up_too() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = std::collections::HashMap::new();
        nk.insert(
            std::sync::Arc::from("TiO2"),
            vec![num_complex::Complex64::new(2.35, 0.0); NW],
        );
        let graded = vec![crate::structure::Layer {
            inhomogen: true,
            inh_delta: 0.1,
            ..crate::structure::Layer::film(5.0, "TiO2")
        }];
        let bg: std::collections::HashSet<String> = ["TiO2".to_string()].into_iter().collect();
        let (stack, _) = DesignStack::from_design(
            air(),
            sub(),
            &graded,
            &nk,
            &std::collections::HashMap::new(),
            &wl,
            &bg,
        )
        .unwrap();
        assert_eq!(stack.spans().len(), 1);
        // The final-pass clamp_up branch is the strongest form of the
        // deferral: even told to clamp up, a multi-row span goes whole.
        let mut stack = stack;
        let rep = stack.clamp_all(6.0, 1000.0, true).unwrap();
        assert_eq!(rep.spans_removed.len(), 1);
        assert_eq!(rep.rows_removed, 4);
        assert!(stack.films().is_empty());
    }
}
