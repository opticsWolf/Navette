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
use crate::smatrix::synthesis::context::DesignContext;
use crate::smatrix::synthesis::jacobian::CurveDeposits;
use crate::smatrix::synthesis::merit::{CurveId, MeritSpec, SimCurves};
use crate::smatrix::synthesis::structure::{DesignStack, SolverArrays};
use crate::smatrix::synthesis::jacobian::assemble_jacobian;
use crate::smatrix::synthesis::thick_opt::{
    levenberg_marquardt_with, JacobianMode, JacobianSource, LmConfig, LmResult, NoJacobian,
};

/// Solver + merit context for one synthesis problem.
#[derive(Clone)]
pub struct SmatrixContext {
    pub wavls: Vec<f64>,
    /// Sines of incidence angles (solver convention).
    pub sin_theta: Vec<f64>,
    pub spec: MeritSpec,
    pub clamp_min_nm: f64,
    pub clamp_max_nm: f64,
    pub lm: LmConfig,
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
                    .map(|l| Complex64::new(sa.n_stack_cache[base + l * 2], sa.n_stack_cache[base + l * 2 + 1]))
                    .collect();
                let inv_n_slice: Vec<Complex64> =
                    n_slice.iter().map(|&n| 1.0 / n).collect();
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
                let pt = Pt { rs: s_res.4, rp: p_res.4, ts: s_res.5, tp: p_res.5,
                              tfs: s_res.2, tfp: p_res.2 };
                // Thickness deposits, when asked for. Nothing above this line
                // is touched by the request: the values are the same solver
                // call either way, and the fingerprint says so.
                let dep = if n_par == 0 {
                    Vec::new()
                } else {
                    deposit_row(
                        &n_slice, &sa, start, end, lam, nsin_fi, par_films, n_par,
                    )
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
            cplx[CurveId::Ts.index()] =
                Some(pts.iter().map(|p| p.tfs).collect::<Vec<_>>().into());
            cplx[CurveId::Tp.index()] =
                Some(pts.iter().map(|p| p.tfp).collect::<Vec<_>>().into());
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

    fn optimize_thicknesses(
        &mut self,
        stack: &mut DesignStack,
    ) -> Result<f64, String> {
        self.optimize_thicknesses_report(stack).map(|(mf, _)| mf)
    }
}

impl SmatrixContext {
    /// `optimize_thicknesses`, plus the solver's own account of the run.
    ///
    /// The report is `None` when there was nothing to optimize (no
    /// optimize-flagged films), which is not a solve and has no diagnostics.
    /// Otherwise it is the LM's `LmResult` — `analytic_jacobians` says which
    /// Jacobian path actually ran, `gain_ratio` and `termination` say how the
    /// run ended. The trait method drops it; callers that want to know keep
    /// it.
    pub fn optimize_thicknesses_report(
        &mut self,
        stack: &mut DesignStack,
    ) -> Result<(f64, Option<LmResult>), String> {
        // Collect optimize-flagged film indices and their starting values.
        let opt_indices: Vec<usize> = stack
            .films()
            .iter()
            .enumerate()
            .filter(|(_, l)| l.optimize)
            .map(|(i, _)| i)
            .collect();
        if opt_indices.is_empty() {
            return self.evaluate_merit(stack).map(|mf| (mf, None));
        }

        // Owned copies for the (Send+Sync) residual closure.
        let base_stack = stack.clone();
        let spec = self.spec.clone();
        let ctx_self = self.clone();
        let indices = opt_indices.clone();

        let residuals = move |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            let mut st = base_stack.clone();
            for (j, &i) in indices.iter().enumerate() {
                st.set_thickness(i, x[j])?;
            }
            let sim = ctx_self.simulate(&st)?;
            spec.residuals(&sim, out).map_err(|c| format!("missing curve {c:?}"))
        };

        let x0: Vec<f64> = opt_indices
            .iter()
            .map(|&i| stack.films()[i].d_nm.clamp(0.0, self.clamp_max_nm))
            .collect();
        let lb = vec![0.0f64; x0.len()];
        let ub = vec![self.clamp_max_nm; x0.len()];

        // The analytic Jacobian, when the spec is fully covered by it. The
        // source decides per call: a spec with a phase target or a color
        // demand declines and the driver differences that iteration, so the
        // mode is a preference, not a promise.
        let res = match self.lm.jacobian {
            JacobianMode::Analytic => {
                let src = DepositJacobian {
                    ctx: self.clone(),
                    base_stack: stack.clone(),
                    indices: opt_indices.clone(),
                    n_wav: self.wavls.len(),
                };
                levenberg_marquardt_with(&residuals, Some(&src), &x0, &lb, &ub, &self.lm)?
            },
            JacobianMode::Fd => {
                levenberg_marquardt_with(&residuals, None::<&NoJacobian>, &x0, &lb, &ub, &self.lm)?
            },
        };

        // Write back, then clamp sweep (removes sub-min, caps above-max).
        for (j, &i) in opt_indices.iter().enumerate() {
            if i < stack.films().len() {
                stack.set_thickness(i, res.x[j])?;
            }
        }
        stack.clamp_all(self.clamp_min_nm, self.clamp_max_nm);

        self.evaluate_merit(stack).map(|mf| (mf, Some(res)))
    }
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
    indices: Vec<usize>,
    n_wav: usize,
}

impl JacobianSource for DepositJacobian {
    fn fill(&self, x: &[f64], jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
        let mut st = self.base_stack.clone();
        for (j, &i) in self.indices.iter().enumerate() {
            st.set_thickness(i, x[j])?;
        }
        let (sim, dep) = self.ctx.simulate_with_deposits(&st, &self.indices)?;
        let sens = self
            .ctx
            .spec
            .curve_sensitivity(&sim)
            .map_err(|c| format!("missing curve {c:?}"))?;
        if !sens.is_complete() {
            return Ok(None);
        }
        assemble_jacobian(&sens, &dep, self.n_wav, jac)?;
        Ok(Some(sens.rows.len()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smatrix::synthesis::structure::LayerSpec;
    use crate::smatrix::synthesis::merit::{ConstraintKind, MeritKey, MeritTarget, SimTransform};

    fn ar_spec(angle: f64, wl: f64) -> MeritSpec {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey { angle, curve: CurveId::Rs });
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
        }
    }

    /// Context whose spec demands PDts (Ts, phase, 1 pass): `simulate()`
    /// must then fill complex-t rows + reference metadata (gated assembly).
    fn ar_ctx_pd() -> SmatrixContext {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey { angle: 0.0, curve: CurveId::Ts });
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
            assert!((rs[k] + ts[k] - 1.0).abs() < 1e-10, "k={k} R={} T={}", rs[k], ts[k]);
        }
        // All values physical.
        for arr in [sim.curve(CurveId::Rs).unwrap(), sim.curve(CurveId::Rp).unwrap()] {
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
            ambient, substrate, vec![LayerSpec::constant("F", 1.0, 0.0, 500.0, nw)],
        ).unwrap();
        // NOTE: empty spec → gated assembly skips complex rows, so this
        // calibration context carries a (value-irrelevant) phase demand.
        let mut pd_spec = MeritSpec::new();
        let pk = pd_spec.add_key(MeritKey { angle: 0.0, curve: CurveId::Ts });
        pd_spec.add_target(MeritTarget {
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
        };
        let sim = ctx.simulate(&slab).unwrap();
        let tf = sim.cplx[CurveId::Ts.index()].as_ref().unwrap()[0];
        // kD = 2π·500/400 = 2.5π → +π/2 vs −π/2, unambiguous.
        assert!((tf.arg() - std::f64::consts::PI / 2.0).abs() < 1e-9, "tf={tf}");
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
        let k = spec.add_key(MeritKey { angle: 0.0, curve: CurveId::Ts });
        let wl = vec![900.0, 1000.0, 1100.0];
        let tgt: Vec<f64> = wl.iter().map(|&lam| {
            // `.conj()`: solver convention (see oracle test above).
            oracle_tf(1.0, n_l, 1.52, d, lam).conj().arg() - 2.0 * std::f64::consts::PI * d / lam
        }).collect();
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
        let k = spec.add_key(MeritKey { angle: 0.0, curve: CurveId::Ts });
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
      assemble_jacobian, deposit_channel, CurveDeposits, DEPOSIT_CHANNELS,
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
      }
    }

    /// A pointwise Exact demand on `curve` over the whole grid at `angle`.
    fn demand(spec: &mut MeritSpec, angle: f64, curve: CurveId, target: f64) {
      demand_at(spec, angle, curve, vec![target; WLS.len()]);
    }

    /// As `demand`, with a per-wavelength target vector.
    fn demand_at(spec: &mut MeritSpec, angle: f64, curve: CurveId, targets: Vec<f64>) {
      let k = spec.add_key(MeritKey { angle, curve });
      spec
        .add_target(MeritTarget {
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
          let fd: Vec<f64> =
            ca.iter().zip(cb.iter()).map(|(p, m)| (p - m) / (2.0 * h)).collect();
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
      assert!(live > 40, "only {live} live entries -- the sweep barely moved");
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
        ctx.spec.residuals(&ctx.simulate(&sp).unwrap(), &mut rp).unwrap();
        ctx.spec.residuals(&ctx.simulate(&sm).unwrap(), &mut rm).unwrap();
        let fd: Vec<f64> =
          rp.iter().zip(&rm).map(|(p, q)| (p - q) / (2.0 * h)).collect();
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
        assert!((jac[w] - expect).abs() <= 1e-12 * expect.abs().max(1.0),
                "row {w}: {} vs {expect}", jac[w]);
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
        indices: vec![0, 1, 2],
        n_wav: WLS.len(),
      };
      let mut jac = Vec::new();
      let m = covered.fill(&x, &mut jac).unwrap().expect("covered spec declined");
      assert_eq!(m, 2 * WLS.len());
      assert_eq!(jac.len(), m * 3);
      assert!(jac.iter().any(|v| v.abs() > 1e-9), "a covered spec gave a flat J");

      let phase = DepositJacobian {
        ctx: ar_ctx_pd(),
        base_stack: stack,
        indices: vec![0],
        n_wav: WLS.len(),
      };
      assert!(phase.fill(&[118.0], &mut jac).unwrap().is_none(),
              "a phase spec must decline, not answer");
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
      assert!((ma - mb).abs() <= 1e-6 * mb.abs().max(1e-12),
              "merit {ma} vs {mb}");
      for (i, (fa, fb)) in sa.films().iter().zip(sb.films()).enumerate() {
        assert!((fa.d_nm - fb.d_nm).abs() < 1e-4,
                "film {i}: {} vs {} nm", fa.d_nm, fb.d_nm);
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
  }
}
