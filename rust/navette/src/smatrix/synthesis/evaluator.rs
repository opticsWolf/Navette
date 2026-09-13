// SPDX-License-Identifier: LGPL-3.0-or-later
//! Navette -- Rust Rewrite of Numba-optimized thin-film optical solver
//!
//! synthesis::evaluator — SmatrixContext: the REAL DesignContext.
//!
//! Wires together, all in-process (no Python round trip):
//!   * coherent_block::solve_coherent_block_fields_dual  — one dual-pol
//!     interface sweep per spectral point (rayon over angle×λ grid)
//!   * MeritSpec::merit / ::residuals                    — targets engine
//!   * thick_opt::levenberg_marquardt                    — bounded LM
//!
//! optimize_thicknesses mirrors ClampedNeedleSynthesizer.optimize_thicknesses:
//! LM over optimize-flagged films with bounds [0, clamp_max], then a
//! clamp_all(clamp_min, clamp_max) sweep that REMOVES sub-min layers.

use std::sync::Arc;

use num_complex::Complex64;
use rayon::prelude::*;

use crate::smatrix::coherent_block::solve_coherent_block_fields_dual;
use crate::smatrix::needle_operator::{
    block_flux_factors, build_stack_fields_range, needle_slopes4_ddz,
};
use crate::smatrix::synthesis::config::ThinLayerPolicy;
use crate::smatrix::synthesis::context::DesignContext;
use crate::smatrix::synthesis::jacobian::CurveDeposits;
use crate::smatrix::synthesis::jacobian::assemble_jacobian_mapped;
use crate::smatrix::synthesis::merit::{CurveId, MeritSpec, SimCurves};
use crate::smatrix::synthesis::optimizer::{OptimizerResult, run_optimizer};
use crate::smatrix::synthesis::structure::{ClampReport, DesignStack, SolverArrays};
use crate::smatrix::synthesis::thick_opt::{JacobianMode, JacobianSource, LmConfig, NoJacobian};

/// Solver + merit context for one synthesis problem.
#[derive(Clone)]
pub struct SmatrixContext {
    pub wavls: Vec<f64>,
    /// Sines of incidence angles (solver convention).
    pub sin_theta: Vec<f64>,
    pub spec: MeritSpec,
    pub clamp_min_nm: f64,
    pub clamp_max_nm: f64,
    /// F0.3: the floor policy - `ClampUpAlways` also moves the LM lower
    /// bound to `clamp_min_nm` (see the limit-cycle argument on the enum).
    pub thin_layer_policy: ThinLayerPolicy,
    pub lm: LmConfig,
    /// F0.2: clamp reports from every `optimize_thicknesses` sweep since
    /// the last drain. The pipeline drains it into the phase result and
    /// the accumulator resets; the residual closures' clones carry a
    /// snapshot that is never read back.
    pub clamp_accumulator: ClampReport,
}

impl SmatrixContext {
    /// Simulate the stack on the fixed grid → SimCurves.
    ///
    /// Fully-coherent path: block = [0, nl−1) with the substrate as
    /// half-space. Rs/Rp/Ts/Tp assembled row-major (k = a·nw + w).
    pub fn simulate(&self, stack: &DesignStack) -> Result<SimCurves, String> {
        self.simulate_inner(stack, &[]).map(|(sim, _)| sim)
    }

    /// Simulate, and on the same sweep read off ∂(curve)/∂(thickness) for the
    /// films in `par_films` (indices into `stack.films()`, in the optimizer's
    /// parameter order) — the right half of the analytic Jacobian (R4.5).
    ///
    /// The deposits cost one extra `StackFields` decomposition per point and
    /// polarization, plus O(1) per film: the derivative of a layer's
    /// thickness is the needle operator with the needle material set to the
    /// host's own index, where ρ̂ vanishes and τ̂ is the bare propagation
    /// slope iβ_j. So the whole Jacobian is a constant multiple of one
    /// simulate, not the 2n of the central-difference path.
    pub fn simulate_with_deposits(
        &self,
        stack: &DesignStack,
        par_films: &[usize],
    ) -> Result<(SimCurves, CurveDeposits), String> {
        if par_films.is_empty() {
            return Err("simulate_with_deposits: no parameter films".into());
        }
        let n_films = stack.films().len();
        if let Some(&bad) = par_films.iter().find(|&&i| i >= n_films) {
            return Err(format!("simulate_with_deposits: film {bad} out of range"));
        }
        let (sim, dep) = self.simulate_inner(stack, par_films)?;
        dep.ok_or_else(|| "simulate_with_deposits: deposits not produced".to_string())
            .map(|d| (sim, d))
    }

    fn simulate_inner(
        &self,
        stack: &DesignStack,
        par_films: &[usize],
    ) -> Result<(SimCurves, Option<CurveDeposits>), String> {
        let sa: SolverArrays = stack.solver_arrays();
        let nl = sa.n_layers as usize;
        if nl < 2 {
            return Err("degenerate stack".into());
        }
        let nw = self.wavls.len();
        let na = self.sin_theta.len();
        let start = 0usize;
        let end = nl - 1;
        if end <= start {
            return Err("stack needs at least ambient + substrate".into());
        }

        // Per-point intensities + forward-transmission amplitudes via the
        // dual solver. BlockResult order is (rf, tb, tf, rb, R, Tb, T, Rb):
        // `.2` is the complex forward t (front incidence) — `.1`/`.5` are
        // the BACKWARD (reciprocal) quantities, equal in intensity by
        // reciprocity but NOT in phase. PDts/PDtp need `.2`.
        struct Pt {
            rs: f64,
            rp: f64,
            ts: f64,
            tp: f64,
            tfs: Complex64,
            tfp: Complex64,
        }
        let n_par = par_films.len();
        let swept: Vec<(Pt, Vec<f64>)> = (0..na * nw)
            .into_par_iter()
            .map(|k| {
                let a = k / nw;
                let w = k % nw;
                let lam = self.wavls[w];
                let base = w * nl * 2;
                let n_slice: Vec<Complex64> = (0..nl)
                    .map(|l| {
                        Complex64::new(
                            sa.n_stack_cache[base + l * 2],
                            sa.n_stack_cache[base + l * 2 + 1],
                        )
                    })
                    .collect();
                let inv_n_slice: Vec<Complex64> = n_slice.iter().map(|&n| 1.0 / n).collect();
                let nsin_fi = n_slice[0] * Complex64::new(self.sin_theta[a], 0.0);
                let (s_res, p_res) = solve_coherent_block_fields_dual(
                    start,
                    end,
                    &n_slice,
                    &inv_n_slice,
                    &sa.thicknesses,
                    &sa.rough_vals,
                    &sa.rough_types,
                    lam,
                    nsin_fi,
                );
                let pt = Pt {
                    rs: s_res.4,
                    rp: p_res.4,
                    ts: s_res.5,
                    tp: p_res.5,
                    tfs: s_res.2,
                    tfp: p_res.2,
                };
                // Thickness deposits, when asked for. Nothing above this line
                // is touched by the request: the values are the same solver
                // call either way, and the fingerprint says so.
                let dep = if n_par == 0 {
                    Vec::new()
                } else {
                    deposit_row(&n_slice, &sa, start, end, lam, nsin_fi, par_films, n_par)
                };
                (pt, dep)
            })
            .collect();
        let deposits = if n_par == 0 {
            None
        } else {
            let rows: Vec<Vec<f64>> = swept.iter().map(|(_, d)| d.clone()).collect();
            Some(CurveDeposits::from_point_rows(rows, n_par)?)
        };
        let pts: Vec<Pt> = swept.into_iter().map(|(p, _)| p).collect();

        let mk = |id: CurveId| -> Arc<[f64]> {
            let idx = id.index();
            let pick = |p: &Pt| match id {
                CurveId::Rs => p.rs,
                CurveId::Rp => p.rp,
                CurveId::Ts => p.ts,
                CurveId::Tp => p.tp,
                _ => unreachable!("only R/T s/p curves produced"),
            };
            let _ = idx;
            let v: Vec<f64> = pts.iter().map(pick).collect();
            v.into()
        };
        let mut curves = [None, None, None, None, None, None, None, None, None];
        curves[CurveId::Rs.index()] = Some(mk(CurveId::Rs));
        curves[CurveId::Rp.index()] = Some(mk(CurveId::Rp));
        curves[CurveId::Ts.index()] = Some(mk(CurveId::Ts));
        curves[CurveId::Tp.index()] = Some(mk(CurveId::Tp));
        // Complex forward-t rows for (differential-)phase demands — ONLY
        // when the spec asks (two allocations + O(grid) copies saved per
        // merit call for intensity-only optimizations; the amplitudes
        // themselves come from the dual solver regardless).
        let mut cplx: [Option<Arc<[Complex64]>>; 6] = [None, None, None, None, None, None];
        if self.spec.uses_phase() {
            cplx[CurveId::Ts.index()] = Some(pts.iter().map(|p| p.tfs).collect::<Vec<_>>().into());
            cplx[CurveId::Tp.index()] = Some(pts.iter().map(|p| p.tfp).collect::<Vec<_>>().into());
        }
        // Stack metadata for the PD reference: ambient/substrate thickness
        // entries are zero, so the plain sum is the coating thickness D;
        // ambient index at the centre wavelength (dispersive ambients are
        // pathological — the scalar is a documented approximation).
        // Gated the same way (defaults zero the reference anyway).
        let (total_d, n_front_re, n_back_re) = if self.spec.uses_differential() {
            let total_d: f64 = sa.thicknesses.iter().sum();
            let n_front_re = sa.n_stack_cache[(nw / 2) * nl * 2];
            // Substrate exit index (back-phase reference; unused front-only).
            let n_back_re = sa.n_stack_cache[(nw / 2) * nl * 2 + (nl - 1) * 2];
            (total_d, n_front_re, n_back_re)
        } else {
            (0.0, 1.0, 1.0)
        };

        Ok((
            SimCurves {
                angles: self.sin_theta.clone().into(),
                wavelengths: self.wavls.clone().into(),
                total_d,
                n_front_re,
                n_back_re,
                curves,
                cplx,
                ..Default::default()
            },
            deposits,
        ))
    }
}

/// ∂(Rs, Rp, Ts, Tp)/∂(thickness) at one grid point, laid out
/// `[channel · n_par + par]` for [`CurveDeposits::from_point_rows`].
///
/// Increasing film `j`'s thickness is the needle operator with the needle
/// material set to the host's own index: r₁₂ = 0, so ρ̂ vanishes and
/// τ̂ = iβ_j is the bare propagation slope. The composition
/// `U ⊗ N_dual ⊗ L` then differentiates the whole block — every
/// multiple-reflection path included — through the same Redheffer star
/// product the forward solver uses, which is why this cannot drift from it.
///
/// The intensities follow the solver's own `finalize`: R = |r_f|² and the
/// stored T is |t_back|²·f_back, with `f_back` the boundary-admittance ratio.
/// That factor depends only on the block's two half-spaces, so it is
/// constant under a film thickness and rides outside the derivative.
#[allow(clippy::too_many_arguments)]
fn deposit_row(
    n_slice: &[Complex64],
    sa: &SolverArrays,
    start: usize,
    end: usize,
    lam: f64,
    nsin_fi: Complex64,
    par_films: &[usize],
    n_par: usize,
) -> Vec<f64> {
    let mut row = vec![0.0f64; 4 * n_par];
    for pol in [0i32, 1i32] {
        let pi = pol as usize;
        let fields = build_stack_fields_range(
            start,
            end,
            n_slice,
            &sa.thicknesses,
            &sa.rough_vals,
            &sa.rough_types,
            lam,
            nsin_fi,
            pol,
        );
        let m = fields.s_left[fields.end];
        let f_back = block_flux_factors(&fields, pol)[1];
        let r_conj = m.0.conj();
        let t_conj = m.1.conj();
        for (kp, &fi) in par_films.iter().enumerate() {
            // Film `fi` is solver slot `fi + 1` (ambient occupies slot 0).
            let j = fi + 1;
            let s4 = needle_slopes4_ddz(&fields, nsin_fi, j, 0.0, fields.n[j], pol, lam);
            // d|z|²/dd = 2·Re(conj(z)·dz/dd).
            row[pi * n_par + kp] = 2.0 * (r_conj * s4[0]).re;
            row[(2 + pi) * n_par + kp] = 2.0 * f_back * (t_conj * s4[1]).re;
        }
    }
    row
}

impl DesignContext for SmatrixContext {
    fn evaluate_merit(&self, stack: &DesignStack) -> Result<f64, String> {
        let sim = self.simulate(stack)?;
        Ok(self.spec.merit(&sim, 1e6))
    }

    fn simulate(&self, stack: &DesignStack) -> Result<SimCurves, String> {
        SmatrixContext::simulate(self, stack)
    }

    fn optimize_thicknesses(&mut self, stack: &mut DesignStack) -> Result<f64, String> {
        self.optimize_thicknesses_report(stack).map(|(mf, _)| mf)
    }

    fn take_clamp_report(&mut self) -> Option<ClampReport> {
        if self.clamp_accumulator.is_empty() {
            return None;
        }
        Some(std::mem::take(&mut self.clamp_accumulator))
    }
}

impl SmatrixContext {
    /// `optimize_thicknesses`, plus the solver's own account of the run.
    ///
    /// The report is `None` when there was nothing to optimize (no
    /// optimize-flagged films), which is not a solve and has no diagnostics.
    /// Otherwise it is the solver's `OptimizerResult` — `backend` says which
    /// solver ran, `analytic_jacobians` which Jacobian path it took, and
    /// `gain_ratio` / `termination` how the run ended. The trait method drops
    /// it; callers that want to know keep it.
    pub fn optimize_thicknesses_report(
        &mut self,
        stack: &mut DesignStack,
    ) -> Result<(f64, Option<OptimizerResult>), String> {
        // F1.6 (U2/B5): the parameter list stops being a row list. A
        // profiled carrier's span is ONE physical layer, so its total
        // thickness is one parameter, distributed over the span's bulk
        // rows by the frozen fractions. Everything else is today's row
        // parameter, and the list is the same list it is today.
        let params = build_params(stack);
        if params.is_empty() {
            return self.evaluate_merit(stack).map(|mf| (mf, None));
        }

        // Owned copies for the (Send+Sync) residual closure.
        let base_stack = stack.clone();
        let spec = self.spec.clone();
        let ctx_self = self.clone();
        let params_owned = params.clone();

        let residuals = move |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            let mut st = base_stack.clone();
            apply_params(&mut st, &params_owned, x)?;
            let sim = ctx_self.simulate(&st)?;
            spec.residuals(&sim, out)
                .map_err(|c| format!("missing curve {c:?}"))
        };

        let x0: Vec<f64> = params
            .iter()
            .map(|p| p.value0(stack).clamp(0.0, self.clamp_max_nm))
            .collect();
        // F0.3 (U1): the lower bound moves with the policy. Under
        // `ClampUpAlways` the floor is a hard bound - without it, LM drives
        // a film to 0.5 nm, the clamp puts it back to 2.0, and the pair
        // oscillates until the stagnation detector terminates the run with
        // a true report of a false condition. With it, a design whose
        // optimum wants a 0.5 nm film converges to the floor in ONE
        // optimization. `Remove`/`ClampUpFinal` keep today's `0.0` - the
        // search runs exactly as before, elimination and all.
        let lb_floor = if self.thin_layer_policy == ThinLayerPolicy::ClampUpAlways {
            self.clamp_min_nm
        } else {
            0.0
        };
        let lb = vec![lb_floor; x0.len()];
        let ub = vec![self.clamp_max_nm; x0.len()];

        // The analytic Jacobian, when the spec is fully covered by it. The
        // source decides per call: a spec with a phase target or a color
        // demand declines and the driver differences that iteration, so the
        // mode is a preference, not a promise.
        //
        // Which *solver* consumes them is `self.lm.backend`'s business, and
        // `run_optimizer`'s: the built-in LM is the default and the only one
        // with bounds semantics of its own (R4.4c).
        let res = match self.lm.jacobian {
            JacobianMode::Analytic => {
                let src = DepositJacobian {
                    ctx: self.clone(),
                    base_stack: stack.clone(),
                    params: params.clone(),
                    n_wav: self.wavls.len(),
                };
                run_optimizer(&residuals, Some(&src), &x0, &lb, &ub, &self.lm)?
            }
            JacobianMode::Fd => {
                run_optimizer(&residuals, None::<&NoJacobian>, &x0, &lb, &ub, &self.lm)?
            }
        };

        // Write back, then clamp sweep (removes sub-min, caps above-max).
        // F0.2: the sweep's report accumulates into the context and the
        // pipeline drains it into the phase result at record time - a
        // per-call message here would be noise (this fires after every
        // thickness optimization, dozens of times per cycle).
        apply_params(stack, &params, &res.x)?;
        // During the run: F0.3's pass structure (only ClampUpAlways
        // clamps up mid-run; the pipeline's final pass carries
        // ClampUpFinal). F1.6: scalable spans follow the same structure.
        let rep = stack.clamp_all_policy(
            self.clamp_min_nm,
            self.clamp_max_nm,
            self.thin_layer_policy,
            false,
        )?;
        self.clamp_accumulator.merge(rep);

        self.evaluate_merit(stack).map(|mf| (mf, Some(res)))
    }
}

/// F1.6 (U2/B5): one LM parameter.
///
/// - `Row(i)` is film i's own thickness — today's parameter, and the only
///   kind a stack with no scalable span produces, in the same order.
/// - `Span` is one profiled carrier's TOTAL thickness `D`, distributed over
///   the span's bulk rows by fractions frozen at build time
///   (`phi_r = d_r / D0`). The interface slice is not in `rows` and is not
///   scaled: it is an interface property, not part of the layer's
///   thickness. Scalable = the span has at least two bulk rows and every
///   one of them is optimize-flagged — the shape `from_design` gives a
///   profiled (`inhomogen` or `gradient`) carrier with `optimize = true`
///   whose mode scales exactly (`InhMode::Fixed`, `GradientMode::FixedSpan`);
///   the RateCapped modes depend on absolute depth and stay homogenized
///   until F1.7.
#[derive(Clone, Debug)]
pub(crate) enum Param {
    Row(usize),
    Span {
        rows: Vec<usize>,
        fractions: Vec<f64>,
    },
}

impl Param {
    /// The parameter's starting value: the row's thickness, or the span
    /// total D0 the fractions were frozen at.
    fn value0(&self, stack: &DesignStack) -> f64 {
        match self {
            Param::Row(i) => stack.films()[*i].d_nm,
            Param::Span { rows, .. } => rows.iter().map(|&r| stack.films()[r].d_nm).sum(),
        }
    }
}

/// The parameter list for one `optimize_thicknesses` call.
///
/// Span parameters first-class: each scalable span contributes ONE
/// parameter (its bulk rows in row order, fractions frozen from the
/// current row thicknesses). Every other optimize-flagged row is a `Row`
/// parameter. A stack with no scalable span produces exactly today's
/// `opt_indices` list, as `Row` params.
pub(crate) fn build_params(stack: &DesignStack) -> Vec<Param> {
    let films = stack.films();
    let mut params = Vec::new();
    for sp in stack.spans() {
        if stack.span_is_scalable(sp) {
            let rows: Vec<usize> = (sp.bulk_start..sp.end).collect();
            let d0: f64 = rows.iter().map(|&r| films[r].d_nm).sum();
            let fractions: Vec<f64> = rows.iter().map(|&r| films[r].d_nm / d0).collect();
            params.push(Param::Span { rows, fractions });
        } else {
            for r in sp.bulk_start..sp.end {
                if films[r].optimize {
                    params.push(Param::Row(r));
                }
            }
        }
    }
    params
}

/// Write one parameter vector into a stack. `Row` writes its row; a span
/// writes `phi_r * D` per bulk row (frozen fractions, moving total). At
/// `x[j] == D0` the rebuild agrees with the original stack to the
/// multiply-divide rounding (~1 ulp per row) — the one place this path is
/// not bitwise, and no fingerprint touches a span stack.
pub(crate) fn apply_params(
    stack: &mut DesignStack,
    params: &[Param],
    x: &[f64],
) -> Result<(), String> {
    for (j, p) in params.iter().enumerate() {
        match p {
            Param::Row(i) => stack.set_thickness(*i, x[j])?,
            Param::Span { rows, fractions } => {
                for (&r, &f) in rows.iter().zip(fractions.iter()) {
                    stack.set_thickness(r, f * x[j])?;
                }
            }
        }
    }
    Ok(())
}

/// The analytic Jacobian of one `optimize_thicknesses` problem: rebuild the
/// stack at `x`, take one sweep with deposits, and multiply the two halves.
///
/// Declines (`Ok(None)`) when the merit spec has rows whose dependence the
/// sensitivity pass does not carry — phase targets and color demands. That is
/// a property of the spec, not of `x`, so a declining source declines every
/// iteration and the run is simply the finite-difference run.
struct DepositJacobian {
    ctx: SmatrixContext,
    base_stack: DesignStack,
    params: Vec<Param>,
    n_wav: usize,
}

impl JacobianSource for DepositJacobian {
    fn fill(&self, x: &[f64], jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
        let mut st = self.base_stack.clone();
        apply_params(&mut st, &self.params, x)?;
        // The deposits are per FILM ROW (the flat list: each parameter's
        // rows in parameter order), then the mapped assembly contracts
        // row columns into parameter columns with the frozen weights.
        let mut flat_rows = Vec::new();
        for p in &self.params {
            match p {
                Param::Row(i) => flat_rows.push(*i),
                Param::Span { rows, .. } => flat_rows.extend_from_slice(rows),
            }
        }
        let (sim, dep) = self.ctx.simulate_with_deposits(&st, &flat_rows)?;
        let sens = self
            .ctx
            .spec
            .curve_sensitivity(&sim)
            .map_err(|c| format!("missing curve {c:?}"))?;
        if !sens.is_complete() {
            return Ok(None);
        }
        assemble_jacobian_mapped(&sens, &dep, self.n_wav, &self.params, jac)?;
        Ok(Some(sens.rows.len()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smatrix::synthesis::merit::{ConstraintKind, MeritKey, MeritTarget, SimTransform};
    use crate::smatrix::synthesis::structure::LayerSpec;

    fn ar_spec(angle: f64, wl: f64) -> MeritSpec {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle,
            curve: CurveId::Rs,
        });
        spec.add_target(MeritTarget {
            key_idx: k as u32,
            wavelengths: vec![wl].into(),
            kind: ConstraintKind::Exact,
            transform: SimTransform::Linear,
            norm_factor: 1.0,
            normalized_targets: vec![0.0].into(), // R = 0 demanded
            tolerances: vec![0.01].into(),
            band: vec![].into(),
            phase: false,
            differential_passes: None,
            integral: false,
            weight: 1.0,
            count_norm: None,
        })
        .unwrap();
        spec
    }

    /// air | L(d, n=1.2329) | sub(n=1.52) — ideal single-layer AR.
    fn ar_stack(d: f64) -> DesignStack {
        let nw = 3;
        let n_l = 1.52_f64.sqrt(); // perfect AR index at normal incidence
        let mut ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
        ambient.optimize = false;
        ambient.needle = false;
        let mut substrate = LayerSpec::constant("sub", 1.52, 0.0, 0.0, nw);
        substrate.optimize = false;
        substrate.needle = false;
        DesignStack::with_films(
            ambient,
            substrate,
            vec![LayerSpec::constant("L", n_l, 0.0, d, nw)],
        )
        .unwrap()
    }

    fn ar_ctx(clamp_max: f64) -> SmatrixContext {
        SmatrixContext {
            wavls: vec![900.0, 1000.0, 1100.0],
            sin_theta: vec![0.0],
            spec: ar_spec(0.0, 1000.0),
            clamp_min_nm: 2.0,
            clamp_max_nm: clamp_max,
            lm: LmConfig::default(),
            clamp_accumulator: ClampReport::default(),
            thin_layer_policy: crate::smatrix::synthesis::config::ThinLayerPolicy::Remove,
        }
    }

    /// Context whose spec demands PDts (Ts, phase, 1 pass): `simulate()`
    /// must then fill complex-t rows + reference metadata (gated assembly).
    fn ar_ctx_pd() -> SmatrixContext {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ts,
        });
        spec.add_target(MeritTarget {
            key_idx: k as u32,
            wavelengths: vec![1000.0].into(),
            kind: ConstraintKind::Exact,
            transform: SimTransform::Phase,
            norm_factor: 1.0,
            normalized_targets: vec![0.0].into(),
            tolerances: vec![0.05].into(),
            band: vec![].into(),
            phase: true,
            differential_passes: Some(1.0),
            integral: false,
            weight: 1.0,
            count_norm: None,
        })
        .unwrap();
        SmatrixContext {
            wavls: vec![900.0, 1000.0, 1100.0],
            sin_theta: vec![0.0],
            spec,
            clamp_min_nm: 2.0,
            clamp_max_nm: 1000.0,
            lm: LmConfig::default(),
            clamp_accumulator: ClampReport::default(),
            thin_layer_policy: crate::smatrix::synthesis::config::ThinLayerPolicy::Remove,
        }
    }

    #[test]
    fn energy_conservation_lossless_stack() {
        // R + T == 1 for a lossless coherent stack at normal incidence.
        let ctx = ar_ctx(1000.0);
        let stack = ar_stack(200.0);
        let sim = ctx.simulate(&stack).unwrap();

        let rs = sim.curve(CurveId::Rs).unwrap();
        let ts = sim.curve(CurveId::Ts).unwrap();
        for k in 0..rs.len() {
            assert!(
                (rs[k] + ts[k] - 1.0).abs() < 1e-10,
                "k={k} R={} T={}",
                rs[k],
                ts[k]
            );
        }
        // All values physical.
        for arr in [
            sim.curve(CurveId::Rs).unwrap(),
            sim.curve(CurveId::Rp).unwrap(),
        ] {
            for &v in arr.iter() {
                assert!((0.0..=1.0).contains(&v));
            }
        }
    }

    #[test]
    fn quarter_wave_ar_is_zero_reflection() {
        // Exact analytic anchor: d = λ/(4n) gives R = 0 for n₁ = √(n₀n_s).
        let n_l = 1.52_f64.sqrt();
        let d_qw = 1000.0 / (4.0 * n_l);
        let ctx = ar_ctx(1000.0);
        let stack = ar_stack(d_qw);
        let mf = ctx.evaluate_merit(&stack).unwrap();
        assert!(mf < 1e-16, "mf={mf}");

        // Off-quarter-wave has nonzero reflection.
        let stack_off = ar_stack(d_qw * 0.7);
        assert!(ctx.evaluate_merit(&stack_off).unwrap() > 1e-4);
    }

    #[test]
    fn optimizer_recovers_quarter_wave_from_bad_start() {
        // Start far off (350 nm); LM must find the quarter-wave thickness.
        let n_l = 1.52_f64.sqrt();
        let d_qw = 1000.0 / (4.0 * n_l);
        let mut ctx = ar_ctx(1000.0);
        let mut stack = ar_stack(350.0);

        let mf = ctx.optimize_thicknesses(&mut stack).unwrap();
        assert!(mf < 1e-6, "mf={mf}");
        let d_final = stack.films()[0].d_nm;
        // FD-Jacobian LM lands within a fraction of a nm here (smooth 1-D).
        assert!(
            (d_final - d_qw).abs() < 1.0,
            "d_final={d_final}, expected {d_qw}"
        );
    }

    /// Independent 2×2 characteristic-matrix oracle (s-pol, normal
    /// incidence): M = [[cosδ, i·sinδ/n],[i·n·sinδ, cosδ]],
    /// t = 2·n₀/(n₀·M₁₁ + n₀·n_s·M₁₂ + M₂₁ + n_s·M₂₂).
    fn oracle_tf(n0: f64, n1: f64, ns: f64, d: f64, lam: f64) -> Complex64 {
        use num_complex::Complex64 as C;
        let delta = 2.0 * std::f64::consts::PI * n1 * d / lam;
        let (s, c) = delta.sin_cos();
        let m11 = C::new(c, 0.0);
        let m12 = C::new(0.0, s / n1);
        let m21 = C::new(0.0, s * n1);
        let m22 = C::new(c, 0.0);
        C::new(2.0 * n0, 0.0)
            / (C::new(n0, 0.0) * m11 + C::new(n0 * ns, 0.0) * m12 + m21 + C::new(ns, 0.0) * m22)
    }

    #[test]
    fn solver_propagation_sign_matches_reference() {
        // All-matched stack (n = 1 everywhere, film D = 500): tf is pure
        // propagation — its arg IS the solver's propagation-phase sign,
        // and `reference_phase` must reproduce it (convention lock for
        // every differential-phase demand in the crate).
        let nw = 3;
        let ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
        let substrate = LayerSpec::constant("sub", 1.0, 0.0, 0.0, nw);
        let slab = DesignStack::with_films(
            ambient,
            substrate,
            vec![LayerSpec::constant("F", 1.0, 0.0, 500.0, nw)],
        )
        .unwrap();
        // NOTE: empty spec → gated assembly skips complex rows, so this
        // calibration context carries a (value-irrelevant) phase demand.
        let mut pd_spec = MeritSpec::new();
        let pk = pd_spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ts,
        });
        pd_spec
            .add_target(MeritTarget {
                key_idx: pk as u32,
                wavelengths: vec![400.0].into(),
                kind: ConstraintKind::Exact,
                transform: SimTransform::Phase,
                norm_factor: 1.0,
                normalized_targets: vec![0.0].into(),
                tolerances: vec![0.05].into(),
                band: vec![].into(),
                phase: true,
                differential_passes: None,
                integral: false,
                weight: 1.0,
                count_norm: None,
            })
            .unwrap();
        let ctx = SmatrixContext {
            wavls: vec![400.0],
            sin_theta: vec![0.0],
            spec: pd_spec,
            clamp_min_nm: 2.0,
            clamp_max_nm: 1000.0,
            lm: LmConfig::default(),
            clamp_accumulator: ClampReport::default(),
            thin_layer_policy: crate::smatrix::synthesis::config::ThinLayerPolicy::Remove,
        };
        let sim = ctx.simulate(&slab).unwrap();
        let tf = sim.cplx[CurveId::Ts.index()].as_ref().unwrap()[0];
        // kD = 2π·500/400 = 2.5π → +π/2 vs −π/2, unambiguous.
        assert!(
            (tf.arg() - std::f64::consts::PI / 2.0).abs() < 1e-9,
            "tf={tf}"
        );
        // `reference_phase` is unwrapped (2.5π here) while `arg()` wraps:
        // compare in wrapped space, exactly as the merit kernel does.
        let r = crate::smatrix::optics_core::reference_phase(400.0, 1.0, 0.0, 500.0, 1.0);
        let rw = r - std::f64::consts::TAU * (r / std::f64::consts::TAU).round();
        assert!((rw - tf.arg()).abs() < 1e-9, "ref={r} arg={}", tf.arg());
    }

    #[test]
    fn gated_assembly_skips_unrequested_rows() {
        // Intensity-only spec: no complex rows, default metadata (the
        // premise — unrequested paths stay dark in the hot LM loop).
        let ctx = ar_ctx(1000.0);
        assert!(!ctx.spec.uses_phase());
        assert!(!ctx.spec.uses_differential());
        let sim = ctx.simulate(&ar_stack(200.0)).unwrap();
        assert!(sim.cplx.iter().all(|c| c.is_none()));
        assert_eq!(sim.total_d, 0.0);
        assert_eq!((sim.n_front_re, sim.n_back_re), (1.0, 1.0));
        // Phase-demanding spec flips both gates.
        let ctx_pd = ar_ctx_pd();
        assert!(ctx_pd.spec.uses_phase());
        assert!(ctx_pd.spec.uses_differential());
    }

    #[test]
    fn simulate_fills_pd_metadata_and_complex_t() {
        let n_l = 1.52_f64.sqrt();
        let d = 200.0;
        let ctx = ar_ctx_pd();
        let stack = ar_stack(d);
        let sim = ctx.simulate(&stack).unwrap();
        // Metadata: coating thickness + media.
        assert!((sim.total_d - d).abs() < 1e-12, "total_d={}", sim.total_d);
        assert!((sim.n_front_re - 1.0).abs() < 1e-12);
        assert!((sim.n_back_re - 1.52).abs() < 1e-12);
        // Complex-t rows present for Ts/Tp, consistent with intensities:
        // |tf|² × flux (n_s/n_0 at normal incidence) == Ts row.
        let ts = sim.curve(CurveId::Ts).unwrap();
        let tfs = &sim.cplx[CurveId::Ts.index()].as_ref().unwrap();
        assert_eq!(tfs.len(), ts.len());
        for k in 0..ts.len() {
            assert!((tfs[k].norm_sqr() * 1.52 - ts[k]).abs() < 1e-10, "k={k}");
        }
        // Independent oracle per wavelength (pins tf-vs-tb: the backward
        // amplitude has a different phase in asymmetric stacks).
        // `.conj()`: the solver's forward-propagation convention is the
        // conjugate of Macleod textbooks (see `reference_phase` docs) —
        // magnitudes/physics identical, phase sign flipped crate-wide.
        for (k, &lam) in [900.0, 1000.0, 1100.0].iter().enumerate() {
            let t_oracle = oracle_tf(1.0, n_l, 1.52, d, lam).conj();
            let diff = (tfs[k] - t_oracle).norm();
            assert!(diff < 1e-9, "k={k} lam={lam} diff={diff}");
        }
    }

    #[test]
    fn pd_merit_matches_hand_delta_phi() {
        // Δφ = arg(t_oracle) − 2π·D/λ (air, normal) demanded exactly → ~0.
        let n_l = 1.52_f64.sqrt();
        let d = 200.0;
        let mut ctx = ar_ctx_pd();
        let stack = ar_stack(d);
        let sim = ctx.simulate(&stack).unwrap();
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ts,
        });
        let wl = vec![900.0, 1000.0, 1100.0];
        let tgt: Vec<f64> = wl
            .iter()
            .map(|&lam| {
                // `.conj()`: solver convention (see oracle test above).
                oracle_tf(1.0, n_l, 1.52, d, lam).conj().arg()
                    - 2.0 * std::f64::consts::PI * d / lam
            })
            .collect();
        spec.add_target(MeritTarget {
            key_idx: k as u32,
            wavelengths: wl.into(),
            kind: ConstraintKind::Exact,
            transform: SimTransform::Phase,
            norm_factor: 1.0,
            normalized_targets: tgt.into(),
            tolerances: vec![0.05, 0.05, 0.05].into(),
            band: vec![].into(),
            phase: true,
            differential_passes: Some(1.0),
            integral: false,
            weight: 1.0,
            count_norm: None,
        })
        .unwrap();
        assert!(spec.merit(&sim, 1e6) < 1e-20, "m={}", spec.merit(&sim, 1e6));
        ctx.spec = spec;
        assert!(ctx.evaluate_merit(&stack).unwrap() < 1e-20);
    }

    #[test]
    fn optimizer_recovers_thickness_from_pd_target() {
        // End-to-end PD loop: demand the QW design's Δφ, start at 350 nm.
        // Convergence is on MERIT (phase wraps admit 2π-branch solutions).
        let n_l = 1.52_f64.sqrt();
        let d_qw = 1000.0 / (4.0 * n_l);
        let qw_stack = ar_stack(d_qw);
        let probe = ar_ctx_pd();
        let qw_sim = probe.simulate(&qw_stack).unwrap();
        let qw_phase = qw_sim.cplx[CurveId::Ts.index()].as_ref().unwrap()[1].arg();
        let ref_qw = 2.0 * std::f64::consts::PI * d_qw / 1000.0;
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Ts,
        });
        spec.add_target(MeritTarget {
            key_idx: k as u32,
            wavelengths: vec![1000.0].into(),
            kind: ConstraintKind::Exact,
            transform: SimTransform::Phase,
            norm_factor: 1.0,
            normalized_targets: vec![qw_phase - ref_qw].into(),
            tolerances: vec![0.05].into(),
            band: vec![].into(),
            phase: true,
            differential_passes: Some(1.0),
            integral: false,
            weight: 1.0,
            count_norm: None,
        })
        .unwrap();
        let mut ctx = ar_ctx(1000.0);
        ctx.spec = spec;
        let mut stack = ar_stack(350.0);
        let mf = ctx.optimize_thicknesses(&mut stack).unwrap();
        assert!(mf < 1e-6, "mf={mf}");
    }

    // -- R4.5 increment ii: thickness deposits ---------------------------------
    mod deposits {
        use super::*;
        use crate::smatrix::synthesis::jacobian::{
            CurveDeposits, DEPOSIT_CHANNELS, assemble_jacobian, deposit_channel,
        };

        const WLS: [f64; 5] = [480.0, 520.0, 560.0, 600.0, 640.0];
        const SINS: [f64; 2] = [0.0, 0.45];

        /// air | H(n=2.30) | L(n=1.46) | H(n=2.30, k=0.004) | glass(1.52).
        ///
        /// Three optimizable films of different index, one of them weakly
        /// absorbing (so A = 1 - R - T is not identically zero and the two
        /// polarizations genuinely differ off normal).
        fn stack3(d0: f64, d1: f64, d2: f64) -> DesignStack {
            let nw = WLS.len();
            let mut ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
            ambient.optimize = false;
            ambient.needle = false;
            let mut substrate = LayerSpec::constant("glass", 1.52, 0.0, 0.0, nw);
            substrate.optimize = false;
            substrate.needle = false;
            DesignStack::with_films(
                ambient,
                substrate,
                vec![
                    LayerSpec::constant("H", 2.30, 0.0, d0, nw),
                    LayerSpec::constant("L", 1.46, 0.0, d1, nw),
                    LayerSpec::constant("Habs", 2.30, 0.004, d2, nw),
                ],
            )
            .unwrap()
        }

        fn ctx3(spec: MeritSpec) -> SmatrixContext {
            SmatrixContext {
                wavls: WLS.to_vec(),
                sin_theta: SINS.to_vec(),
                spec,
                clamp_min_nm: 2.0,
                clamp_max_nm: 1000.0,
                lm: LmConfig::default(),
                clamp_accumulator: ClampReport::default(),
                thin_layer_policy: crate::smatrix::synthesis::config::ThinLayerPolicy::Remove,
            }
        }

        /// A pointwise Exact demand on `curve` over the whole grid at `angle`.
        fn demand(spec: &mut MeritSpec, angle: f64, curve: CurveId, target: f64) {
            demand_at(spec, angle, curve, vec![target; WLS.len()]);
        }

        /// As `demand`, with a per-wavelength target vector.
        fn demand_at(spec: &mut MeritSpec, angle: f64, curve: CurveId, targets: Vec<f64>) {
            let k = spec.add_key(MeritKey { angle, curve });
            spec.add_target(MeritTarget {
                key_idx: k as u32,
                wavelengths: WLS.to_vec().into(),
                kind: ConstraintKind::Exact,
                transform: SimTransform::Linear,
                norm_factor: 1.0,
                normalized_targets: targets.into(),
                tolerances: vec![0.02; WLS.len()].into(),
                band: vec![].into(),
                phase: false,
                differential_passes: None,
                integral: false,
                weight: 1.0,
                count_norm: None,
            })
            .unwrap();
        }

        #[test]
        fn thickness_deposits_are_a_finite_difference_of_the_simulated_curves() {
            // The right half of J, against its own oracle: `simulate()` centrally
            // differenced in THICKNESS space. No merit, no optimizer -- if this
            // and the merit-side finite difference both hold, the product holds.
            let ctx = ctx3(ar_spec(0.0, 560.0));
            let d = [118.0, 203.0, 64.0];
            let stack = stack3(d[0], d[1], d[2]);
            let par = [0usize, 1, 2];
            let (_, dep) = ctx.simulate_with_deposits(&stack, &par).unwrap();
            assert_eq!(dep.n_points(), SINS.len() * WLS.len());
            assert_eq!(dep.n_par(), par.len());

            let h = 1e-4;
            let mut worst = 0.0f64;
            let mut live = 0usize;
            for (kp, &fi) in par.iter().enumerate() {
                let mut sp = stack.clone();
                sp.set_thickness(fi, d[fi] + h).unwrap();
                let mut sm = stack.clone();
                sm.set_thickness(fi, d[fi] - h).unwrap();
                let a = ctx.simulate(&sp).unwrap();
                let b = ctx.simulate(&sm).unwrap();
                for (ch, &id) in DEPOSIT_CHANNELS.iter().enumerate() {
                    let ca = a.curve(id).unwrap();
                    let cb = b.curve(id).unwrap();
                    let fd: Vec<f64> = ca
                        .iter()
                        .zip(cb.iter())
                        .map(|(p, m)| (p - m) / (2.0 * h))
                        .collect();
                    let scale = fd.iter().fold(1e-6f64, |acc, v| acc.max(v.abs()));
                    for (pt, &f) in fd.iter().enumerate() {
                        let an = dep.get(ch, pt, kp);
                        let dev = (an - f).abs() / scale;
                        worst = worst.max(dev);
                        assert!(dev < 1e-6, "film {fi}, {id:?}, point {pt}: {an} vs {f}");
                        if f.abs() > 1e-6 {
                            live += 1;
                        }
                    }
                }
            }
            assert!(
                live > 40,
                "only {live} live entries -- the sweep barely moved"
            );
            println!("  deposits: {live} live entries, worst rel dev {worst:.3e}");
        }

        #[test]
        fn the_assembled_jacobian_matches_a_finite_difference_of_the_residuals() {
            // Both halves together, against the oracle the LM would otherwise use:
            // central differences of the whole residual vector in thickness space.
            // R, T and A demands at two angles, both polarizations.
            let mut spec = MeritSpec::new();
            demand(&mut spec, SINS[0], CurveId::Rs, 0.01);
            demand(&mut spec, SINS[1], CurveId::Rp, 0.02);
            demand(&mut spec, SINS[1], CurveId::Ts, 0.95);
            demand(&mut spec, SINS[0], CurveId::As, 0.0);
            let ctx = ctx3(spec);

            let d = [118.0, 203.0, 64.0];
            let stack = stack3(d[0], d[1], d[2]);
            let par = [0usize, 1, 2];
            let (sim, dep) = ctx.simulate_with_deposits(&stack, &par).unwrap();
            let sens = ctx.spec.curve_sensitivity(&sim).unwrap();
            assert!(sens.is_complete());
            let mut jac = Vec::new();
            assemble_jacobian(&sens, &dep, WLS.len(), &mut jac).unwrap();

            let m = sens.rows.len();
            assert_eq!(m, 4 * WLS.len());
            assert_eq!(jac.len(), m * par.len());

            let h = 1e-4;
            let mut worst = 0.0f64;
            let mut live = 0usize;
            for (kp, &fi) in par.iter().enumerate() {
                let mut sp = stack.clone();
                sp.set_thickness(fi, d[fi] + h).unwrap();
                let mut sm = stack.clone();
                sm.set_thickness(fi, d[fi] - h).unwrap();
                let mut rp = Vec::new();
                let mut rm = Vec::new();
                ctx.spec
                    .residuals(&ctx.simulate(&sp).unwrap(), &mut rp)
                    .unwrap();
                ctx.spec
                    .residuals(&ctx.simulate(&sm).unwrap(), &mut rm)
                    .unwrap();
                let fd: Vec<f64> = rp
                    .iter()
                    .zip(&rm)
                    .map(|(p, q)| (p - q) / (2.0 * h))
                    .collect();
                let scale = fd.iter().fold(1e-6f64, |acc, v| acc.max(v.abs()));
                for i in 0..m {
                    let an = jac[i * par.len() + kp];
                    let dev = (an - fd[i]).abs() / scale;
                    worst = worst.max(dev);
                    assert!(dev < 1e-6, "row {i}, film {fi}: {an} vs {}", fd[i]);
                    if fd[i].abs() > 1e-6 {
                        live += 1;
                    }
                }
            }
            assert!(live > 30, "only {live} live entries");
            println!("  analytic J: {live} live entries, worst rel dev {worst:.3e}");
        }

        #[test]
        fn absorption_rows_pick_up_both_companions_through_the_deposits() {
            // A = 1 - R - T is the one demand whose row reads two channels. Its
            // Jacobian row must equal -(dR/dd + dT/dd)/tol, a different number
            // from either channel alone -- a single-channel bug would look
            // plausible and be wrong.
            let mut spec = MeritSpec::new();
            demand(&mut spec, SINS[0], CurveId::As, 0.0);
            let ctx = ctx3(spec);
            let stack = stack3(118.0, 203.0, 64.0);
            let par = [2usize];
            let (sim, dep) = ctx.simulate_with_deposits(&stack, &par).unwrap();
            let sens = ctx.spec.curve_sensitivity(&sim).unwrap();
            let mut jac = Vec::new();
            assemble_jacobian(&sens, &dep, WLS.len(), &mut jac).unwrap();

            let cr = deposit_channel(CurveId::Rs).unwrap();
            let ct = deposit_channel(CurveId::Ts).unwrap();
            for w in 0..WLS.len() {
                let expect = -(dep.get(cr, w, 0) + dep.get(ct, w, 0)) / 0.02;
                assert!(
                    (jac[w] - expect).abs() <= 1e-12 * expect.abs().max(1.0),
                    "row {w}: {} vs {expect}",
                    jac[w]
                );
            }
        }

        #[test]
        fn a_phase_demand_refuses_the_analytic_path_instead_of_zeroing_it() {
            // The fallback contract: an uncovered row must stop the assembly, not
            // come back as a row of zeros that the LM would read as "flat".
            let ctx = ar_ctx_pd();
            let stack = stack3(118.0, 203.0, 64.0);
            let par = [0usize];
            let (sim, dep) = ctx.simulate_with_deposits(&stack, &par).unwrap();
            let sens = ctx.spec.curve_sensitivity(&sim).unwrap();
            assert!(!sens.is_complete());
            let mut jac = Vec::new();
            let err = assemble_jacobian(&sens, &dep, WLS.len(), &mut jac).unwrap_err();
            assert!(err.contains("no curve sensitivity"), "{err}");
        }

        #[test]
        fn the_deposits_do_not_disturb_the_simulated_curves() {
            // `simulate_with_deposits` shares one sweep with `simulate`; the values
            // must come back bit-identical, or the fingerprint would depend on
            // whether a Jacobian was asked for.
            let ctx = ctx3(ar_spec(0.0, 560.0));
            let stack = stack3(118.0, 203.0, 64.0);
            let plain = ctx.simulate(&stack).unwrap();
            let (with, _) = ctx.simulate_with_deposits(&stack, &[0, 1, 2]).unwrap();
            for id in DEPOSIT_CHANNELS {
                let a = plain.curve(id).unwrap();
                let b = with.curve(id).unwrap();
                assert_eq!(a.len(), b.len());
                for (x, y) in a.iter().zip(b.iter()) {
                    assert_eq!(x.to_bits(), y.to_bits(), "{id:?}: {x} vs {y}");
                }
            }
        }

        /// F0.3 (U1): the limit cycle and its absence, on the real
        /// machinery. The merit wants the bare substrate's reflectance -
        /// reachable only as the film thins toward zero - so the optimizer
        /// pulls below the floor and the two behaviours part ways:
        /// `ClampUpAlways` moves `lb` to `clamp_min_nm`, so ONE bounded
        /// solve converges to the floor and the merit history is monotone
        /// (the stagnation detector never sees an oscillation); with
        /// today's `Remove` the same stack eliminates the film (the
        /// search keeps its rejection mechanism).
        #[test]
        fn clamp_up_always_moves_the_lm_bound_and_never_oscillates() {
            let r_bare: f64 = {
                let t: f64 = (1.0 - 1.52) / (1.0 + 1.52);
                t * t
            };
            let mut spec = MeritSpec::new();
            let k = spec.add_key(MeritKey {
                angle: 0.0,
                curve: CurveId::Rs,
            });
            spec.add_target(MeritTarget {
                key_idx: k as u32,
                wavelengths: vec![900.0].into(),
                kind: ConstraintKind::Exact,
                transform: SimTransform::Linear,
                norm_factor: 1.0,
                normalized_targets: vec![r_bare].into(),
                tolerances: vec![0.01].into(),
                band: vec![].into(),
                phase: false,
                differential_passes: None,
                integral: false,
                weight: 1.0,
                count_norm: None,
            })
            .map_err(|e| panic!("{e}"))
            .unwrap();

            // ClampUpAlways: lb = clamp_min_nm -> converges to the floor
            // in one solve; five sweeps stay monotone.
            let mut ctx = SmatrixContext {
                wavls: vec![900.0],
                sin_theta: vec![0.0],
                spec: spec.clone(),
                clamp_min_nm: 2.0,
                clamp_max_nm: 1000.0,
                thin_layer_policy: ThinLayerPolicy::ClampUpAlways,
                lm: LmConfig::default(),
                clamp_accumulator: ClampReport::default(),
            };
            // Start THIN: on the thin side of the AR dip the nearest
            // exact-R_bare solution is the film's own disappearance, so
            // the pull is toward the floor (from 200 nm the bounded solve
            // finds the thick-side crossing instead).
            let mut stack = ar_stack(10.0);
            let mut merits: Vec<f64> = Vec::new();
            for _ in 0..5 {
                let mf = ctx.optimize_thicknesses(&mut stack).unwrap();
                merits.push(mf);
            }
            for w in merits.windows(2) {
                assert!(
                    w[1] <= w[0] + 1e-12,
                    "merit history must be monotone (limit cycle's absence): {merits:?}"
                );
            }
            // One optimization reaches the floor; the clamp-up sweep lands
            // the film EXACTLY on it and keeps it there.
            assert!(
                (stack.films()[0].d_nm - 2.0).abs() < 1e-9,
                "film parked at {}, not the floor",
                stack.films()[0].d_nm
            );

            // Remove (default): the same pull eliminates the film - the
            // search keeps its rejection mechanism under the default.
            let mut ctx = SmatrixContext {
                wavls: vec![900.0],
                sin_theta: vec![0.0],
                spec,
                clamp_min_nm: 2.0,
                clamp_max_nm: 1000.0,
                thin_layer_policy: ThinLayerPolicy::Remove,
                lm: LmConfig::default(),
                clamp_accumulator: ClampReport::default(),
            };
            let mut stack = ar_stack(10.0);
            ctx.optimize_thicknesses(&mut stack).unwrap();
            assert!(
                stack.films().is_empty(),
                "the default must keep elimination"
            );
        }

        /// C2 (V3-corrected wording): the interface variant of the bound
        /// twin. The LM floor binds the BULK row (the parameter); the
        /// interface slice is fixed and not a parameter. So a thin
        /// interface-carrying film parks at `clamp_min_nm` BULK - total
        /// `clamp_min_nm + t_slice` - not at the floor: the two floors are
        /// defined on different quantities and the bound is conservative
        /// (C2/V3, stated in `ThinLayerPolicy`'s doc comment). The film is
        /// never deleted and the merit history stays monotone, which is
        /// what the bound exists for.
        #[test]
        fn clamp_up_always_interface_film_binds_the_bulk_row_not_the_total() {
            use std::collections::{HashMap, HashSet};
            let r_bare: f64 = {
                let t: f64 = (1.0 - 1.52) / (1.0 + 1.52);
                t * t
            };
            let mut spec = MeritSpec::new();
            let k = spec.add_key(MeritKey {
                angle: 0.0,
                curve: CurveId::Rs,
            });
            spec.add_target(MeritTarget {
                key_idx: k as u32,
                wavelengths: vec![900.0].into(),
                kind: ConstraintKind::Exact,
                transform: SimTransform::Linear,
                norm_factor: 1.0,
                normalized_targets: vec![r_bare].into(),
                tolerances: vec![0.01].into(),
                band: vec![].into(),
                phase: false,
                differential_passes: None,
                integral: false,
                weight: 1.0,
                count_norm: None,
            })
            .map_err(|e| panic!("{e}"))
            .unwrap();

            let mut ctx = SmatrixContext {
                wavls: vec![900.0],
                sin_theta: vec![0.0],
                spec: spec.clone(),
                clamp_min_nm: 2.0,
                clamp_max_nm: 1000.0,
                thin_layer_policy: ThinLayerPolicy::ClampUpAlways,
                lm: LmConfig::default(),
                clamp_accumulator: ClampReport::default(),
            };

            // air | G 20 (pinned lead) | slice 0.5 | bulk L | sub - the
            // interface slice needs a lead entry (the flag owner resolves
            // only from the second entry on), and the lead must be pinned
            // (not a second parameter). The lead's nk is EXACTLY the
            // ambient's (1.0 + 0j): optically invisible apart from a
            // global phase, and the slice shares the carrier's nk, so the
            // two same-nk rows reduce (transfer matrices multiply) to ONE
            // L layer of 0.5 + bulk - the plain AR structure of the twin
            // above, with the demand's zero now at bulk = -0.5. The pull
            // toward disappearance therefore binds at the floor, and the
            // bound is on the BULK parameter only (V3).
            let nw = 3usize;
            let n_l = 1.52_f64.sqrt();
            let mut ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
            ambient.optimize = false;
            ambient.needle = false;
            let mut substrate = LayerSpec::constant("sub", 1.52, 0.0, 0.0, nw);
            substrate.optimize = false;
            substrate.needle = false;
            let mut nk = HashMap::new();
            nk.insert(
                std::sync::Arc::<str>::from("G"),
                vec![Complex64::new(1.0, 0.0); nw],
            );
            nk.insert(
                std::sync::Arc::<str>::from("L"),
                vec![Complex64::new(n_l, 0.0); nw],
            );
            let mut lead = crate::structure::Layer::film(20.0, "G");
            lead.optimize = false;
            lead.needle = false;
            let mut carrier = crate::structure::Layer::film(10.0, "L");
            carrier.interface = true;
            carrier.interface_thickness = 0.5;
            let wl: Vec<f64> = vec![900.0, 1000.0, 1100.0];
            let (mut stack, warns) = DesignStack::from_design(
                ambient,
                substrate,
                &[lead, carrier],
                &nk,
                &HashMap::new(),
                &wl,
                &HashSet::new(),
            )
            .unwrap();
            assert!(warns.is_empty());
            let sp = stack.spans()[1];
            assert!(sp.slice && sp.is_singleton_bulk());
            let slice_d = stack.films()[sp.start].d_nm;
            assert_eq!(slice_d, 0.5);

            let mut merits: Vec<f64> = Vec::new();
            for _ in 0..5 {
                let mf = ctx.optimize_thicknesses(&mut stack).unwrap();
                merits.push(mf);
            }
            for w in merits.windows(2) {
                assert!(
                    w[1] <= w[0] + 1e-12,
                    "merit history must be monotone (limit cycle's absence): {merits:?}"
                );
            }
            // The BULK parks at the bound (the lb is on the parameter);
            // the total is the floor PLUS the untouched slice. Exactly the
            // V3 asymmetry the clamp-up branch's widening is careful not
            // to promise away.
            let bulk = stack.films()[sp.bulk_start].d_nm;
            assert!(
                (bulk - 2.0).abs() < 1e-9,
                "bulk parked at {bulk}, not the bound"
            );
            assert_eq!(stack.films()[sp.start].d_nm, slice_d, "slice untouched");
            assert_eq!(
                stack.films().len(),
                3,
                "never deleted (lead + slice + bulk)"
            );

            // Control: under Remove the carrier is eliminated whole
            // (slice included - F0.2's span rule); the lead survives.
            let mut ctx2 = SmatrixContext {
                wavls: vec![900.0],
                sin_theta: vec![0.0],
                spec,
                clamp_min_nm: 2.0,
                clamp_max_nm: 1000.0,
                thin_layer_policy: ThinLayerPolicy::Remove,
                lm: LmConfig::default(),
                clamp_accumulator: ClampReport::default(),
            };
            let nw = 3usize;
            let mut ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
            ambient.optimize = false;
            ambient.needle = false;
            let mut substrate = LayerSpec::constant("sub", 1.52, 0.0, 0.0, nw);
            substrate.optimize = false;
            substrate.needle = false;
            let mut lead = crate::structure::Layer::film(20.0, "G");
            lead.optimize = false;
            lead.needle = false;
            let mut carrier = crate::structure::Layer::film(10.0, "L");
            carrier.interface = true;
            carrier.interface_thickness = 0.5;
            let (mut stack, _) = DesignStack::from_design(
                ambient,
                substrate,
                &[lead, carrier],
                &nk,
                &HashMap::new(),
                &wl,
                &HashSet::new(),
            )
            .unwrap();
            ctx2.optimize_thicknesses(&mut stack).unwrap();
            assert_eq!(
                stack.films().len(),
                1,
                "the lead survives; the interface carrier is gone whole"
            );
            assert_eq!(stack.films()[0].material.as_ref(), "G");
        }

        /// F0.2: the accumulator drains once and resets. With no
        /// optimize-flagged films there is no solve and no sweep, so no
        /// report; a report that arrives is handed to the pipeline
        /// exactly once.
        #[test]
        fn the_clamp_accumulator_drains_and_resets() {
            let mut ctx = ar_ctx(1000.0);
            let mut stack = ar_stack(200.0);
            stack.set_row_optimize_for_test(0, false);
            ctx.optimize_thicknesses(&mut stack).unwrap();
            assert!(ctx.take_clamp_report().is_none());

            ctx.clamp_accumulator
                .merge(crate::smatrix::synthesis::structure::ClampReport {
                    spans_removed: vec!["H (1.0 nm)".to_string()],
                    spans_capped: 0,
                    rows_removed: 1,
                });
            let rep = ctx.take_clamp_report().unwrap();
            assert_eq!(rep.rows_removed, 1);
            assert!(ctx.take_clamp_report().is_none(), "the drain resets");
        }

        #[test]
        fn deposits_refuse_a_film_index_that_is_not_there() {
            let ctx = ctx3(ar_spec(0.0, 560.0));
            let stack = stack3(118.0, 203.0, 64.0);
            assert!(ctx.simulate_with_deposits(&stack, &[]).is_err());
            assert!(ctx.simulate_with_deposits(&stack, &[0, 3]).is_err());
        }

        #[test]
        fn the_deposit_source_declines_exactly_the_specs_it_cannot_cover() {
            // The coverage gate, at the seam the optimizer actually uses. A spec
            // whose channels are all deposited answers with a Jacobian; one with a
            // phase target answers `None`, which is what makes the run fall back
            // to differences instead of optimizing against zeros.
            let stack = stack3(118.0, 203.0, 64.0);
            let x = [118.0, 203.0, 64.0];

            let mut spec = MeritSpec::new();
            demand(&mut spec, SINS[0], CurveId::Rs, 0.0);
            demand(&mut spec, SINS[1], CurveId::As, 0.0);
            let covered = DepositJacobian {
                ctx: ctx3(spec),
                base_stack: stack.clone(),
                params: vec![Param::Row(0), Param::Row(1), Param::Row(2)],
                n_wav: WLS.len(),
            };
            let mut jac = Vec::new();
            let m = covered
                .fill(&x, &mut jac)
                .unwrap()
                .expect("covered spec declined");
            assert_eq!(m, 2 * WLS.len());
            assert_eq!(jac.len(), m * 3);
            assert!(
                jac.iter().any(|v| v.abs() > 1e-9),
                "a covered spec gave a flat J"
            );

            let phase = DepositJacobian {
                ctx: ar_ctx_pd(),
                base_stack: stack,
                params: vec![Param::Row(0)],
                n_wav: WLS.len(),
            };
            assert!(
                phase.fill(&[118.0], &mut jac).unwrap().is_none(),
                "a phase spec must decline, not answer"
            );
        }

        #[test]
        fn the_two_jacobian_modes_optimize_to_the_same_stack() {
            // §8's question for this change: does the exact Jacobian move a pinned
            // optimum? Same problem, same start, same bounds — only the Jacobian
            // differs, and the two solvers must agree to far better than the
            // difference noise they are being compared across.
            // Targets taken from a reference stack, so the optimum is interior,
            // unique and sitting at merit ~ 0: no bound activity, no layer the
            // clamp sweep wants to remove, one basin. An "R = 0 everywhere"
            // demand is none of those things — it is multimodal, and it drives a
            // film to zero.
            let reference = stack3(118.0, 203.0, 64.0);
            let probe = ctx3(ar_spec(0.0, 560.0));
            let sim0 = probe.simulate(&reference).unwrap();
            let rs0 = sim0.curve(CurveId::Rs).unwrap().to_vec();
            let rp0 = sim0.curve(CurveId::Rp).unwrap().to_vec();
            let nw = WLS.len();
            let mut spec = MeritSpec::new();
            demand_at(&mut spec, SINS[0], CurveId::Rs, rs0[..nw].to_vec());
            demand_at(&mut spec, SINS[1], CurveId::Rp, rp0[nw..].to_vec());

            let mut analytic = ctx3(spec.clone());
            analytic.lm.jacobian = JacobianMode::Analytic;
            let mut differenced = ctx3(spec);
            differenced.lm.jacobian = JacobianMode::Fd;
            // The post-LM clamp sweep REMOVES sub-minimum films, which changes the
            // parameter count and makes the two answers incomparable as vectors.
            // Compare the optimizers, not the cleanup: floor the removal out.
            analytic.clamp_min_nm = 1e-9;
            differenced.clamp_min_nm = 1e-9;

            // Start a few nm off that optimum, so both runs are unambiguously in
            // its basin. Reflectance against thickness is oscillatory: started far
            // away, two local solvers legitimately land in different minima, and
            // comparing across basins measures the landscape rather than the
            // Jacobian. (With an "R = 0" demand and a distant start they do
            // diverge here — and the analytic run finds the better minimum.)
            let mut sa = stack3(122.0, 198.0, 67.0);
            let mut sb = sa.clone();
            let ma = analytic.optimize_thicknesses(&mut sa).unwrap();
            let mb = differenced.optimize_thicknesses(&mut sb).unwrap();

            assert_eq!(sa.films().len(), sb.films().len(), "one run dropped a film");
            assert!(
                (ma - mb).abs() <= 1e-6 * mb.abs().max(1e-12),
                "merit {ma} vs {mb}"
            );
            for (i, (fa, fb)) in sa.films().iter().zip(sb.films()).enumerate() {
                assert!(
                    (fa.d_nm - fb.d_nm).abs() < 1e-4,
                    "film {i}: {} vs {} nm",
                    fa.d_nm,
                    fb.d_nm
                );
            }
            println!("  modes agree: merit {ma:.6e} vs {mb:.6e}");
        }

        #[test]
        fn the_deposit_rows_have_to_be_the_shape_the_channels_expect() {
            let bad = vec![vec![0.0; 5], vec![0.0; 5]];
            assert!(CurveDeposits::from_point_rows(bad, 2).is_err());
            let good = vec![vec![1.0; 8], vec![2.0; 8]];
            let d = CurveDeposits::from_point_rows(good, 2).unwrap();
            assert_eq!(d.n_points(), 2);
            assert_eq!(d.n_par(), 2);
            assert_eq!(d.get(3, 1, 1), 2.0);
        }

        // ------------------------------------------------------------
        // F1.6 - one thickness parameter per graded span
        // ------------------------------------------------------------

        /// A span stack: one graded carrier (optimize=true, the F1.6
        /// posture) expanded by `from_design` on the WLS grid, flanked by
        /// pinned ambient/substrate.
        use std::collections::{HashMap, HashSet};

        /// The F1.6 test grid: ar_ctx's wavelengths (ar_spec's 1000 nm
        /// demand is on-grid here).
        const GWL: [f64; 3] = [900.0, 1000.0, 1100.0];

        /// A span stack on the GWL grid: one graded carrier
        /// (optimize=true, the F1.6 posture) expanded by `from_design`,
        /// flanked by pinned ambient/substrate. `delta = None` builds a
        /// PLAIN film (singleton span, Row parameter) for the no-span
        /// control.
        fn span_stack(delta: Option<f64>, d: f64) -> DesignStack {
            let nw = GWL.len();
            let mut nk = HashMap::new();
            nk.insert(
                std::sync::Arc::<str>::from("H"),
                vec![Complex64::new(2.35, 0.01); nw],
            );
            let mut carrier = crate::structure::Layer::film(d, "H");
            if let Some(delta) = delta {
                carrier.inhomogen = true;
                carrier.inh_delta = delta;
            }
            let mut ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
            ambient.optimize = false;
            ambient.needle = false;
            let mut substrate = LayerSpec::constant("sub", 1.52, 0.0, 0.0, nw);
            substrate.optimize = false;
            substrate.needle = false;
            let (stack, warns) = DesignStack::from_design(
                ambient,
                substrate,
                std::slice::from_ref(&carrier),
                &nk,
                &HashMap::new(),
                &GWL,
                &HashSet::new(),
            )
            .unwrap();
            assert!(
                warns.is_empty(),
                "the F1.6 posture must not warn: {warns:?}"
            );
            stack
        }

        /// The span param's (rows, fractions, D0) for a stack with
        /// exactly one span parameter.
        fn span_param(stack: &DesignStack) -> (Vec<usize>, Vec<f64>, f64) {
            let params = build_params(stack);
            assert_eq!(params.len(), 1, "one span, one parameter");
            match &params[0] {
                Param::Span { rows, fractions } => {
                    let d0: f64 = rows.iter().map(|&r| stack.films()[r].d_nm).sum();
                    (rows.clone(), fractions.clone(), d0)
                }
                Param::Row(_) => panic!("expected a span parameter"),
            }
        }

        /// J3 - the no-span bit-exactness gate, asserted not assumed: a
        /// stack with only Row parameters must assemble the SAME Jacobian
        /// through the parameter-map path as the pre-refactor row-level
        /// assembly, to the last bit.
        #[test]
        fn f16_no_span_jacobian_is_bitwise_the_row_assembly() {
            let stack = span_stack(None, 118.0);
            let params = build_params(&stack);
            assert_eq!(params.len(), 1);
            assert!(matches!(params[0], Param::Row(0)));
            let ctx = ar_ctx(1000.0);
            let mut flat_rows = Vec::new();
            for p in &params {
                if let Param::Row(i) = p {
                    flat_rows.push(*i);
                }
            }
            let (sim, dep) = ctx.simulate_with_deposits(&stack, &flat_rows).unwrap();
            let sens = ctx.spec.curve_sensitivity(&sim).unwrap();
            assert!(sens.is_complete());
            let mut old_style = Vec::new();
            assemble_jacobian(&sens, &dep, GWL.len(), &mut old_style).unwrap();
            let mut mapped = Vec::new();
            assemble_jacobian_mapped(&sens, &dep, GWL.len(), &params, &mut mapped).unwrap();
            assert_eq!(old_style, mapped, "bitwise, not close");
        }

        /// J2 - the weight twin: the span's analytic column is
        /// `sum_r phi_r * (the analytic row columns)`, accumulated in the
        /// span's row order, BITWISE. This is the twin that catches a
        /// wrong weight.
        #[test]
        fn f16_span_column_is_the_frozen_weighted_sum_bitwise() {
            let stack = span_stack(Some(0.2), 200.0);
            let (rows, fractions, _) = span_param(&stack);
            assert!(rows.len() >= 2);
            let ctx = ar_ctx(1000.0);
            let params = build_params(&stack);
            let (sim, dep) = ctx.simulate_with_deposits(&stack, &rows).unwrap();
            let sens = ctx.spec.curve_sensitivity(&sim).unwrap();
            assert!(sens.is_complete());
            // The row-level columns (the pre-refactor contract).
            let mut row_jac = Vec::new();
            assemble_jacobian(&sens, &dep, GWL.len(), &mut row_jac).unwrap();
            // The mapped columns.
            let mut mapped = Vec::new();
            assemble_jacobian_mapped(&sens, &dep, GWL.len(), &params, &mut mapped).unwrap();
            let m = sens.rows.len();
            for i in 0..m {
                let mut acc = 0.0f64;
                for (&r, &w) in rows.iter().zip(&fractions) {
                    acc += w * row_jac[i * rows.len() + r];
                }
                assert_eq!(
                    mapped[i], acc,
                    "residual row {i}: span column bitwise the weighted sum"
                );
            }
        }

        /// J1 - the FD twin: the analytic span column against a central
        /// difference on D, to 1e-6 relative. The deposits make the
        /// analytic column O(1) per row; the difference would cost two
        /// full simulates - exact AND cheaper, which is why analytic is
        /// the right answer here.
        #[test]
        fn f16_span_column_matches_central_difference() {
            let stack = span_stack(Some(0.2), 200.0);
            let (rows, fractions, d0) = span_param(&stack);
            let ctx = ar_ctx(1000.0);
            let params = build_params(&stack);
            let (sim, dep) = ctx.simulate_with_deposits(&stack, &rows).unwrap();
            let sens = ctx.spec.curve_sensitivity(&sim).unwrap();
            let mut mapped = Vec::new();
            assemble_jacobian_mapped(&sens, &dep, GWL.len(), &params, &mut mapped).unwrap();
            // Central difference on D through the SAME parameterization
            // the residual closure uses.
            let spec = ctx.spec.clone();
            let residuals = |dd: f64| -> Vec<f64> {
                let mut st = stack.clone();
                for (&r, &f) in rows.iter().zip(&fractions) {
                    st.set_thickness(r, f * dd).unwrap();
                }
                let sim = ctx.simulate(&st).unwrap();
                let mut out = Vec::new();
                spec.residuals(&sim, &mut out).unwrap();
                out
            };
            let h = 1e-4 * d0.max(1.0);
            let plus = residuals(d0 + h);
            let minus = residuals(d0 - h);
            let m = sens.rows.len();
            for i in 0..m {
                let fd = (plus[i] - minus[i]) / (2.0 * h);
                let an = mapped[i];
                assert!(
                    (fd - an).abs() <= 1e-6 * an.abs().max(1e-12),
                    "residual row {i}: analytic {an} vs fd {fd}"
                );
            }
        }

        /// The parameter builder's no-span behavior: today's list, as
        /// Row params, in today's order - including mixed-flag spans.
        #[test]
        fn f16_build_params_no_span_is_todays_list() {
            let stack = ar_stack(118.0);
            let ps = build_params(&stack);
            assert_eq!(ps.len(), 1);
            assert!(matches!(ps[0], Param::Row(0)));
            // A pinned span (background carrier: no optimize rows)
            // contributes NOTHING, and neighboring free rows still do.
            let nw = WLS.len();
            let mut ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
            ambient.optimize = false;
            ambient.needle = false;
            let mut substrate = LayerSpec::constant("sub", 1.52, 0.0, 0.0, nw);
            substrate.optimize = false;
            substrate.needle = false;
            let mut a = LayerSpec::constant("A", 2.0, 0.0, 50.0, nw);
            a.optimize = true;
            let mut b = LayerSpec::constant("B", 2.0, 0.0, 50.0, nw);
            b.optimize = false;
            let mut c = LayerSpec::constant("C", 2.0, 0.0, 50.0, nw);
            c.optimize = true;
            let stack = DesignStack::with_films(ambient, substrate, vec![a, b, c]).unwrap();
            let ps = build_params(&stack);
            assert_eq!(ps.len(), 2);
            assert!(matches!(ps[0], Param::Row(0)));
            assert!(matches!(ps[1], Param::Row(2)));
        }

        /// P - the profile-preservation twin: after an LM solve that
        /// moves D substantially, every row is bitwise the frozen-fraction
        /// writeback (`phi_r * D`) and every nk row is untouched. The
        /// plan's literal "d_r / D unchanged to the last bit" is not
        /// well-posed - D is re-summed from the rows, so `d_r/D` picks up
        /// the summation's rounding; the writeback contract is the
        /// bitwise oracle, and the re-summed fractions agree to 1e-12.
        #[test]
        fn f16_lm_solve_preserves_the_profile() {
            let mut ctx = ar_ctx(1000.0);
            let mut stack = span_stack(Some(0.2), 200.0);
            let (rows, fractions, d0) = span_param(&stack);
            let nk_before: Vec<Vec<Complex64>> =
                rows.iter().map(|&r| stack.films()[r].nk.to_vec()).collect();

            // Part 1 - the 40% move through the writeback contract
            // itself: the plan's "moves D by 40%" exercised on the SAME
            // operation the residual closure and the post-solve writeback
            // run. The LM's own move size is design-dependent (here the
            // Rs=0 attractor is local), so the fixed-ratio move is
            // asserted through apply_params; part 2 checks a real solve.
            let d_moved = 1.4 * d0;
            {
                let mut st = stack.clone();
                apply_params(&mut st, &build_params(&stack), &[d_moved]).unwrap();
                for (&r, &f) in rows.iter().zip(&fractions) {
                    assert_eq!(
                        st.films()[r].d_nm,
                        f * d_moved,
                        "row {r}: bitwise the frozen-fraction writeback"
                    );
                    assert_eq!(
                        st.films()[r].nk.to_vec(),
                        nk_before[rows.iter().position(|&q| q == r).unwrap()],
                        "row {r}: nk untouched by the scale"
                    );
                }
            }

            // Part 2 - a real LM solve: whatever move it makes, every
            // row is bitwise the frozen-fraction writeback of the
            // parameter it returned, and the profile is undeformed.
            let (mf, res) = ctx.optimize_thicknesses_report(&mut stack).unwrap();
            let res = res.expect("there is a span to optimize");
            let d_new = res.x[0];
            assert!(
                (d_new - d0).abs() > 1e-9,
                "the solve moved D at all (d0 {d0}, d_new {d_new}, mf {mf})"
            );
            for (&r, &f) in rows.iter().zip(&fractions) {
                assert_eq!(
                    stack.films()[r].d_nm,
                    f * d_new,
                    "row {r}: bitwise the frozen-fraction writeback"
                );
                assert_eq!(
                    stack.films()[r].nk.to_vec(),
                    nk_before[rows.iter().position(|&q| q == r).unwrap()],
                    "row {r}: nk untouched"
                );
            }
            // The re-summed fractions agree with the frozen ones to
            // 1e-12 relative (the summation's own rounding).
            let total: f64 = rows.iter().map(|&r| stack.films()[r].d_nm).sum();
            for (&r, &f) in rows.iter().zip(&fractions) {
                let f_now = stack.films()[r].d_nm / total;
                assert!(
                    (f_now - f).abs() <= 1e-12 * f.abs().max(1e-12),
                    "row {r}: fraction drifted from {f} to {f_now}"
                );
            }
        }

        /// The bound twin: the F0.2 cap case (1000 nm graded, ceiling
        /// 300) now marked optimize=true constructs INSTEAD of refusing
        /// (no homogenize warning, no NeedlePipeline-style refusal), and
        /// the optimizer never returns a D above the ceiling.
        #[test]
        fn f16_cap_case_is_a_bound_not_a_refusal() {
            let mut ctx = ar_ctx(300.0);
            let mut stack = span_stack(Some(0.5), 1000.0);
            assert_eq!(stack.films().len(), 57, "the plan's measured row count");
            // The LM runs, respects the ceiling, and the clamp caps the
            // span (fractions preserved) rather than refusing.
            let (_mf, res) = ctx.optimize_thicknesses_report(&mut stack).unwrap();
            let res = res.expect("there is a span to optimize");
            assert!(res.x[0] <= 300.0, "D above the ceiling: {}", res.x[0]);
            let total: f64 = (0..stack.films().len())
                .map(|r| stack.films()[r].d_nm)
                .sum();
            assert!(
                total <= 300.0 + 1e-9,
                "span total {total} above the ceiling"
            );
            // Every row scaled by the same factor: the profile survived.
            let factors: Vec<f64> = (0..stack.films().len())
                .map(|r| stack.films()[r].d_nm)
                .collect();
            let f0 = factors[0];
            for f in &factors {
                assert!((*f / f0 - 1.0).abs() < 1e-9, "uniform scaling");
            }
        }
    }

    // ------------------------------------------------------------------
    // F1.7 - profile refresh for the rate modes (U3/U4)
    // ------------------------------------------------------------------

    use std::collections::{HashMap, HashSet};

    /// The three-point grid `ar_ctx` runs on (mirrors the `deposits`
    /// module's GWL; a local copy so these twins stand alone).
    const F17WL: [f64; 3] = [900.0, 1000.0, 1100.0];

    /// A gradient RateCapped carrier stack on the real GWL AR context:
    /// the carrier keeps its profile (optimize=true, the F1.6/F1.7
    /// posture) as ONE scalable span whose profile is rate-type.
    fn rate_span_stack(d: f64, rate: f64) -> DesignStack {
        let nw = F17WL.len();
        let mut nk = HashMap::new();
        nk.insert(
            std::sync::Arc::<str>::from("H"),
            vec![Complex64::new(2.35, 0.01); nw],
        );
        nk.insert(
            std::sync::Arc::<str>::from("L"),
            vec![Complex64::new(1.46, 0.0); nw],
        );
        let mut carrier = crate::structure::Layer::film(d, "H");
        carrier.gradient = Some(crate::structure::GradientSpec::rate_capped(
            "H", "L", 0.0, rate, 100.0, 0.0, 1.0,
        ));
        let mut ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
        ambient.optimize = false;
        ambient.needle = false;
        let mut substrate = LayerSpec::constant("sub", 1.52, 0.0, 0.0, nw);
        substrate.optimize = false;
        substrate.needle = false;
        let (stack, warns) = DesignStack::from_design(
            ambient,
            substrate,
            std::slice::from_ref(&carrier),
            &nk,
            &HashMap::new(),
            &F17WL,
            &HashSet::new(),
        )
        .unwrap();
        assert!(
            warns.is_empty(),
            "the rate posture must not warn: {warns:?}"
        );
        stack
    }

    /// Row-count-frozen twin: across one real `optimize_thicknesses`
    /// call the span's row count is FROZEN, even though the solve moves
    /// the span's total across the count rule's ceil() boundaries (the
    /// U4 requirement: the residual must not change length mid-solve,
    /// and a live profile would put step discontinuities at the
    /// boundary). The scenario: 12 rows at 221 nm (max_step =
    /// min(20, 900/(10*2.35)) = 19.149); the AR demand's attractor pulls
    /// the solve UP to 260 nm - across the 12/13 boundary - and the
    /// count stays 12 through the whole solve; the refresh AFTER it
    /// re-derives 14. The frozen-then-healed pair is the whole U4
    /// contract in one assertion.
    #[test]
    fn f17_row_count_frozen_across_one_solve_then_healed_by_refresh() {
        // At 221 nm: max_step = min(20, 900/(10*2.35)) = 19.1489361...
        // raw = ceil(221/19.1489...) = ceil(11.54) = 12 rows.
        let mut stack = rate_span_stack(221.0, 0.3);
        let sp0 = stack.spans()[0];
        assert_eq!(sp0.end - sp0.start, 12, "count(221 nm) = 12 rows");
        let d0: f64 = stack.films().iter().map(|l| l.d_nm).sum();
        let mut ctx = ar_ctx(1000.0);
        let (mf, res) = ctx.optimize_thicknesses_report(&mut stack).unwrap();
        let res = res.expect("there is a span to optimize");
        let d_new = res.x[0];
        assert!(
            (d_new - d0).abs() > 1e-9,
            "the solve moved the span at all (d0 {d0}, d_new {d_new}, mf {mf})"
        );
        // Frozen: the row count did not change during the solve, even
        // though the parameter crossed the ceil boundary.
        let sp1 = stack.spans()[0];
        assert_eq!(
            sp1.end - sp1.start,
            12,
            "the row count is FROZEN for the duration of one solve (U4)"
        );
        // The parameter did cross the boundary the count rule would
        // now re-derive at (the run was deliberately parameterized to).
        let would_now = crate::structure::gradient::gradient_sub_layer_count(
            d_new,
            &F17WL,
            &[num_complex::Complex64::new(2.35, 0.01); F17WL.len()],
            &[num_complex::Complex64::new(1.46, 0.0); F17WL.len()],
            None,
        );
        assert_eq!(
            usize::try_from(would_now).unwrap(),
            14,
            "the post-solve extent sits two ceil buckets up ({d_new})"
        );
        // Healed: the construction-point refresh re-derives it.
        stack.refresh_profiles().unwrap();
        let sp2 = stack.spans()[0];
        assert_eq!(
            sp2.end - sp2.start,
            14,
            "the refresh re-derived the count at the moved extent"
        );
        stack.assert_spans_partition();
    }

    /// Drift-bound twin: the staleness a frozen profile carries within
    /// one macro cycle is bounded and small - the plan's number is the
    /// mixing-fraction error `rate * dD / ref` (0.015 at the deepest
    /// sublayer for a 5 nm excursion at rate 0.3, ref 100). Measured on
    /// this scenario (GWL AR demand, 200 nm carrier, rate 0.3): the 5 nm
    /// excursion moves the merit by 26.9 and the staleness correction is
    /// 6.9 - second-order, under half the excursion's own response and
    /// under 3% of the merit. The twin pins BOTH readings (the
    /// plan's "documented tolerance" is the pair: correction < half the
    /// excursion response, and < 5% relative); if a future scenario
    /// breaks it, the cycle is too long for that rate and the docs say
    /// so rather than the code hiding it.
    #[test]
    fn f17_drift_bound_within_one_cycle() {
        let mut stack = rate_span_stack(200.0, 0.3);
        let ctx = ar_ctx(1000.0);
        let m0 = ctx.evaluate_merit(&stack).unwrap();
        // Simulate the LM's within-cycle excursion: a 5 nm move
        // (dD/ref = 0.05, the plan's worked number).
        let sp = stack.spans()[0];
        let d0: f64 = (sp.start..sp.end).map(|r| stack.films()[r].d_nm).sum();
        let f = (d0 + 5.0) / d0;
        for r in sp.start..sp.end {
            let d = stack.films()[r].d_nm * f;
            stack.set_thickness(r, d).unwrap();
        }
        let m_stale = ctx.evaluate_merit(&stack).unwrap();
        // The construction-point refresh reconciles the profile.
        stack.refresh_profiles().unwrap();
        let m_fresh = ctx.evaluate_merit(&stack).unwrap();
        let drift = (m_fresh - m_stale).abs();
        let excursion = (m_stale - m0).abs();
        // The excursion itself was real: the merit moved when the span
        // moved (the scenario is not a no-op).
        assert!(
            excursion > 1e-12,
            "the 5 nm excursion changed the merit at all"
        );
        // The staleness correction is second-order: smaller than half
        // the excursion's own merit response...
        assert!(
            drift < 0.5 * excursion,
            "staleness {drift} < half the excursion response {excursion} (m0 {m0},              stale {m_stale}, fresh {m_fresh})"
        );
        // ...and under 5% of the merit it rides on.
        assert!(
            drift / m_stale.abs().max(1e-12) < 0.05,
            "the within-cycle staleness stays under the documented 5% (drift {drift},              stale {m_stale}, fresh {m_fresh})"
        );
        let _ = m_fresh;
    }
}
