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
use std::sync::atomic::{AtomicBool, Ordering};

use num_complex::Complex64;
use rayon::prelude::*;

use crate::smatrix::coherent_block::{BlockResult, solve_coherent_block_fields_dual};
use crate::smatrix::needle_operator::{
    block_flux_factors, build_stack_fields_range_dual, partition_blocks, spacer_tau, star_real,
    thickness_slopes_ddz, thickness_slopes4_ddz,
};
use crate::smatrix::synthesis::config::ThinLayerPolicy;
use crate::smatrix::synthesis::context::DesignContext;
use crate::smatrix::synthesis::environments::CompiledEnvironments;
use crate::smatrix::synthesis::jacobian::assemble_jacobian_mapped;
use crate::smatrix::synthesis::jacobian::{ALL_SLOTS, BACK_SLOT, CurveDeposits};
use crate::smatrix::synthesis::merit::{CurveId, CurveSet, MeritSpec, Radians, SimCurves};
use crate::smatrix::synthesis::optimizer::{OptimizerResult, run_optimizer};
use crate::smatrix::synthesis::structure::{ClampReport, DesignStack, SolverArrays};
use crate::smatrix::synthesis::thick_opt::{JacobianMode, JacobianSource, LmConfig, NoJacobian};

/// Largest exit-medium extinction coefficient a back-incidence demand
/// accepts. Glass data with residual k ≈ 1e-7…1e-5 passes; coloured glass
/// and every metal or semiconductor in its absorbing band is refused.
pub const BACK_EXIT_MAX_K: f64 = 1e-3;

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
    /// F2.2: the other environments, or `None` for the flat path.
    ///
    /// This field IS multi_environment_plan.md §4.6's single branch. It is
    /// read once per eval, outside every loop, and `None` takes the
    /// pre-F2.2 call sequence op for op — same functions, same signatures,
    /// no slice, no `Vec`.
    ///
    /// It lives on the context rather than on `NeedlePipeline` because
    /// every phase of the pipeline already reaches the solver through
    /// `DesignContext`: putting K here makes the budget check, cleanup and
    /// inflate joint without any of them learning what an environment is,
    /// where a pipeline field would have had to be threaded through all
    /// three (and through `cleanup_design`/`inflate_design`, which take a
    /// `&mut DesignStack` and would then need K of them).
    ///
    /// `Arc` because the LM residual closure is `Send + Sync` and clones
    /// the context per call.
    pub envs: Option<Arc<CompiledEnvironments>>,
}

impl SmatrixContext {
    /// Rows `simulate_inner` can fill: every front and back intensity
    /// (s, p, unpolarized; absorption derives from the R/T pairs) and the
    /// s/p complex amplitudes of R, T, RB and TB. Unpolarized phase has no
    /// meaning and stays out.
    pub fn provides() -> CurveSet {
        CurveSet::default()
            .with_intensity(&[
                CurveId::Rs,
                CurveId::Rp,
                CurveId::Ru,
                CurveId::Ts,
                CurveId::Tp,
                CurveId::Tu,
                CurveId::RBs,
                CurveId::RBp,
                CurveId::RBu,
                CurveId::TBs,
                CurveId::TBp,
                CurveId::TBu,
            ])
            .with_complex(&[
                CurveId::Rs,
                CurveId::Rp,
                CurveId::Ts,
                CurveId::Tp,
                CurveId::RBs,
                CurveId::RBp,
                CurveId::TBs,
                CurveId::TBp,
            ])
    }

    /// Refuse back-incidence demands when the exit medium absorbs over
    /// their wavelengths. Incidence from inside a semi-infinite absorber
    /// has no well-defined incident flux: `T_b = |t_b|²·Re(y_0)/Re(y_exit)`
    /// divides by a vanishing `Re(y_exit)` (T_b = 84 from bare Ag at
    /// 400 nm), and `R_b` and the back phases are referenced to a wave
    /// that cannot reach the stack.
    fn check_back_exit(&self, sa: &SolverArrays, nl: usize) -> Result<(), String> {
        let Some((lo, hi)) = self.spec.back_demand_span() else {
            return Ok(());
        };
        let exit = nl - 1;
        for (w, &lam) in self.wavls.iter().enumerate() {
            if lam < lo || lam > hi {
                continue;
            }
            let k = sa.n_stack_cache[w * nl * 2 + exit * 2 + 1].abs();
            if k > BACK_EXIT_MAX_K {
                return Err(format!(
                    "back-incidence demand: the exit medium absorbs at {lam} nm \
                     (k = {k:.3e} > {BACK_EXIT_MAX_K:e}); incidence from inside a \
                     semi-infinite absorber is undefined, so RB/TB/AB demands \
                     need a transparent exit medium"
                ));
            }
        }
        Ok(())
    }

    /// W3: refuse a spec this evaluator cannot score, before any solve.
    /// Every door that builds a `SmatrixContext` from user input calls it.
    pub fn check_spec(spec: &MeritSpec) -> Result<(), String> {
        spec.check_supported(&Self::provides())
    }

    /// Simulate the stack on the fixed grid → SimCurves.
    ///
    /// Fully-coherent path: block = [0, nl−1) with the substrate as
    /// half-space. Rs/Rp/Ts/Tp assembled row-major (k = a·nw + w).
    pub fn simulate(&self, stack: &DesignStack) -> Result<SimCurves, String> {
        self.simulate_inner(stack, &[], None).map(|(sim, _)| sim)
    }

    /// `simulate`, reading only environment `env`'s demands for which
    /// rows to build (D4). The run door has dropped the phase demands of
    /// every environment with an incoherent plate, so that environment
    /// solves through the cascade while the others keep their phase rows.
    pub fn simulate_env(&self, stack: &DesignStack, env: usize) -> Result<SimCurves, String> {
        self.simulate_inner(stack, &[], Some(env))
            .map(|(sim, _)| sim)
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
        let (sim, dep) = self.simulate_inner(stack, par_films, None)?;
        dep.ok_or_else(|| "simulate_with_deposits: deposits not produced".to_string())
            .map(|d| (sim, d))
    }

    fn simulate_inner(
        &self,
        stack: &DesignStack,
        par_films: &[usize],
        env: Option<usize>,
    ) -> Result<(SimCurves, Option<CurveDeposits>), String> {
        // Complex rows are built for the demands this solve serves: the
        // whole spec on the flat path, one environment's in a joint run.
        let uses_phase = match env {
            None => self.spec.uses_phase(),
            Some(e) => self.spec.uses_phase_in(e),
        };
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
        // `.2` is the complex forward t (front incidence) and `.6` its
        // intensity |t_fwd|²·Re(y_last)/Re(y_first), the engine's Ts. The
        // BACKWARD `.1`/`.5` agree with them only when both half-spaces are
        // transparent: into an absorbing substrate `.5` divides by its
        // near-zero Re(y_last) (T = 84 on bare Ag at 400 nm). Back demands
        // are therefore refused on an absorbing exit medium, below.
        //
        // Every point keeps both polarizations' whole BlockResult; which
        // rows are assembled from them is decided once, from what the spec
        // reads (`required_rows`), so an intensity-only spec builds exactly
        // the four rows it always did.
        let need = self.spec.required_rows();
        if need.any_back() {
            self.check_back_exit(&sa, nl)?;
        }
        // W6: incoherent spacers split the stack into coherent blocks
        // joined by the engine's real intensity cascade. Only intensities
        // cross a spacer, so phase demands and analytic deposits are
        // refused here, by name, rather than read off the first block.
        let incoherent = sa.incoherent_flags.contains(&1);
        let partition = if incoherent {
            if uses_phase {
                return Err("a phase demand on a stack with an incoherent film: phase \
                     is not defined across an incoherent layer, so only intensity, \
                     absorption and colour demands can be synthesized"
                    .into());
            }
            if !par_films.is_empty() {
                return Err("simulate_with_deposits: thickness deposits need a fully \
                     coherent stack (an incoherent film differences instead)"
                    .into());
            }
            Some(partition_blocks(&sa.incoherent_flags))
        } else {
            None
        };
        struct Pt {
            s: BlockResult,
            p: BlockResult,
        }
        let n_par = par_films.len();
        // Amplitude deposits for phase demands, only when the spec has any:
        // an intensity-only sweep collects the same four channels as before.
        let with_phase = n_par > 0 && uses_phase;
        // Back channels (and with them every front amplitude slot) only
        // when the spec reads a back curve.
        let with_back = n_par > 0 && need.any_back();
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
                let (s_res, p_res) = match &partition {
                    None => solve_coherent_block_fields_dual(
                        start,
                        end,
                        &n_slice,
                        &inv_n_slice,
                        &sa.thicknesses,
                        &sa.rough_vals,
                        &sa.rough_types,
                        lam,
                        nsin_fi,
                    ),
                    Some((blocks, spacers)) => solve_cascade_dual(
                        blocks,
                        spacers,
                        &n_slice,
                        &inv_n_slice,
                        &sa,
                        lam,
                        nsin_fi,
                    ),
                };
                let pt = Pt { s: s_res, p: p_res };
                // Thickness deposits, when asked for. Nothing above this line
                // is touched by the request: the values are the same solver
                // call either way, and the fingerprint says so.
                let dep = if n_par == 0 {
                    Vec::new()
                } else {
                    deposit_row(
                        &n_slice, &sa, start, end, lam, nsin_fi, par_films, n_par, with_phase,
                        with_back,
                    )
                };
                (pt, dep)
            })
            .collect();
        let (pts, rows): (Vec<Pt>, Vec<Vec<f64>>) = swept.into_iter().unzip();
        let deposits = if n_par == 0 {
            None
        } else {
            Some(CurveDeposits::from_point_rows(rows, n_par)?)
        };

        // Intensity of one curve at one point. Front: `.4` = |r_f|², `.6` =
        // |t_f|²·f_fwd (the engine's Rs/Ts). Back incidence (from the exit
        // medium): `.7` = |r_b|², `.5` = |t_b|²·f_back — the conventions of
        // the needle operators' RB/TB terms. Unpolarized is the engine's
        // `0.5 * (s + p)` (R_avg/T_avg).
        fn intensity(id: CurveId, p: &Pt) -> f64 {
            match id {
                CurveId::Rs => p.s.4,
                CurveId::Rp => p.p.4,
                CurveId::Ts => p.s.6,
                CurveId::Tp => p.p.6,
                CurveId::RBs => p.s.7,
                CurveId::RBp => p.p.7,
                CurveId::TBs => p.s.5,
                CurveId::TBp => p.p.5,
                _ => {
                    let (a, b) = id
                        .unpolarized_parts()
                        .expect("absorption curves are derived, never stored");
                    0.5 * (intensity(a, p) + intensity(b, p))
                }
            }
        }
        // Complex amplitude: r_f `.0`, t_f `.2`, r_b `.3`, t_b `.1` (the
        // engine's rs_c…tp_c and rbs_c…tbp_c, first coherent block).
        fn amplitude(id: CurveId, p: &Pt) -> Complex64 {
            let r = match id {
                CurveId::Rs | CurveId::Ts | CurveId::RBs | CurveId::TBs => &p.s,
                _ => &p.p,
            };
            match id {
                CurveId::Rs | CurveId::Rp => r.0,
                CurveId::Ts | CurveId::Tp => r.2,
                CurveId::RBs | CurveId::RBp => r.3,
                CurveId::TBs | CurveId::TBp => r.1,
                _ => unreachable!("complex rows exist for s/p R/T/RB/TB only"),
            }
        }
        let row = |id: CurveId| -> Arc<[f64]> {
            pts.iter()
                .map(|p| intensity(id, p))
                .collect::<Vec<f64>>()
                .into()
        };
        let crow = |id: CurveId| -> Arc<[Complex64]> {
            pts.iter()
                .map(|p| amplitude(id, p))
                .collect::<Vec<Complex64>>()
                .into()
        };

        // The four front s/p rows always (every caller reads them); the
        // rest only when the spec reads them.
        let mut curves: [Option<Arc<[f64]>>; 9] = Default::default();
        let mut back: [Option<Arc<[f64]>>; 6] = Default::default();
        for id in CurveId::ALL {
            if id.is_absorption() {
                continue;
            }
            let always = matches!(id, CurveId::Rs | CurveId::Rp | CurveId::Ts | CurveId::Tp);
            if !(always || need.has_intensity(id)) {
                continue;
            }
            match id.back_index() {
                Some(i) => back[i] = Some(row(id)),
                None => curves[id.index()] = Some(row(id)),
            }
        }
        // Complex rows for (differential-)phase demands. Forward t keeps its
        // pre-W2 gate (both rows whenever the spec has any phase demand);
        // r_f and the back amplitudes are built only when read.
        let mut cplx: [Option<Arc<[Complex64]>>; 6] = Default::default();
        let mut cplx_back: [Option<Arc<[Complex64]>>; 4] = Default::default();
        if uses_phase {
            cplx[CurveId::Ts.index()] = Some(crow(CurveId::Ts));
            cplx[CurveId::Tp.index()] = Some(crow(CurveId::Tp));
            for id in [CurveId::Rs, CurveId::Rp] {
                if need.has_complex(id) {
                    cplx[id.index()] = Some(crow(id));
                }
            }
            for (i, id) in [CurveId::RBs, CurveId::RBp, CurveId::TBs, CurveId::TBp]
                .into_iter()
                .enumerate()
            {
                if need.has_complex(id) {
                    cplx_back[i] = Some(crow(id));
                }
            }
        }
        // Stack metadata for the PD reference: the coating thickness D is
        // the sum of the INTERIOR thicknesses (`[1..nl-1]`) - the same
        // expression the engine's PD keys use (solver.rs), so the merit
        // op point and compute(PDts) agree even on the degenerate
        // direct-DesignStack door that lets a non-zero half-space
        // thickness through (review PD G2; bitwise-identical whenever
        // the half-spaces are zero, which every real stack is);
        // incidence/exit indices as PER-WAVELENGTH columns (PD1) - a
        // dispersive medium now carries its real index at every λ instead
        // of the centre-λ scalar, which made Δφ wrong by
        // err(λ) = 2π·D·cosθ·[n(λ) − n(λ_centre)]/λ. Gated the same way
        // (defaults zero the reference anyway).
        let (total_d, n_front_re, n_back_re) = if self.spec.uses_differential() {
            let total_d: f64 = sa.thicknesses[1..nl - 1].iter().sum();
            // Layer 0's real index per wavelength (front reference)...
            let n_front_re: Arc<[f64]> = (0..nw)
                .map(|w| sa.n_stack_cache[w * nl * 2])
                .collect::<Vec<f64>>()
                .into();
            // ...and the substrate's (exit reference; unused front-only).
            let n_back_re: Arc<[f64]> = (0..nw)
                .map(|w| sa.n_stack_cache[w * nl * 2 + (nl - 1) * 2])
                .collect::<Vec<f64>>()
                .into();
            (total_d, n_front_re, n_back_re)
        } else {
            (0.0, Arc::from([1.0_f64]), Arc::from([1.0_f64]))
        };

        Ok((
            SimCurves {
                // The angle axis, rebuilt from the solver's sines.
                angles: self
                    .sin_theta
                    .iter()
                    .map(|&s| Radians::from_checked_sin(s))
                    .collect::<Vec<Radians>>()
                    .into(),
                wavelengths: self.wavls.clone().into(),
                total_d,
                n_front_re,
                n_back_re,
                curves,
                back,
                cplx,
                cplx_back,
            },
            deposits,
        ))
    }
}

/// W6: one point of a stack with incoherent spacers, both polarizations,
/// in the dual solver's `BlockResult` layout. The intensities `.4..=.7`
/// are the cascade totals (R, T_b, T, R_b), walked exactly as the engine's
/// `solve_point` walks them: each coherent block's flux-normalized
/// intensities, then the spacer's attenuation, joined by the real star
/// product from the identity. The amplitudes `.0..=.3` are the first
/// block's, as the engine reports them; `simulate_inner` refuses every
/// demand that would read them.
fn solve_cascade_dual(
    blocks: &[(usize, usize)],
    spacers: &[Option<usize>],
    n_slice: &[Complex64],
    inv_n_slice: &[Complex64],
    sa: &SolverArrays,
    lam: f64,
    nsin_fi: Complex64,
) -> (BlockResult, BlockResult) {
    let mut ig = [[0.0, 1.0, 1.0, 0.0]; 2];
    let mut first: Option<(BlockResult, BlockResult)> = None;
    for (&(bs, be), sp) in blocks.iter().zip(spacers) {
        let (s, p) = solve_coherent_block_fields_dual(
            bs,
            be,
            n_slice,
            inv_n_slice,
            &sa.thicknesses,
            &sa.rough_vals,
            &sa.rough_types,
            lam,
            nsin_fi,
        );
        ig[0] = star_real(ig[0], [s.4, s.5, s.6, s.7]);
        ig[1] = star_real(ig[1], [p.4, p.5, p.6, p.7]);
        first.get_or_insert((s, p));
        if let Some(sp) = *sp {
            let tau = spacer_tau(n_slice[sp], sa.thicknesses[sp], lam, nsin_fi);
            for g in &mut ig {
                *g = star_real(*g, [0.0, tau, tau, 0.0]);
            }
        }
    }
    let (s, p) = first.expect("a stack has at least one coherent block");
    let with = |b: BlockResult, v: [f64; 4]| (b.0, b.1, b.2, b.3, v[0], v[1], v[2], v[3]);
    (with(s, ig[0]), with(p, ig[1]))
}

/// ∂(Rs, Rp, Ts, Tp)/∂(thickness) at one grid point, laid out
/// `[channel · n_par + par]` for [`CurveDeposits::from_point_rows`]; with
/// `with_phase`, then ∂(Re, Im) of the complex forward `t_s`, `t_p`, then
/// of the front reflection `r_s`, `r_p` (see `PHASE_DEPOSIT_CHANNELS`);
/// with `with_back`, all of those and then the back channels
/// (`BACK_DEPOSIT_CHANNELS`): ∂(RBs, RBp, TBs, TBp), then ∂(Re, Im) of
/// `r_b` and `t_b` per polarization.
///
/// Increasing film `j`'s thickness is the needle operator with the needle
/// material set to the host's own index: r₁₂ = 0, so ρ̂ vanishes and
/// τ̂ = iβ_j is the bare propagation slope. The composition
/// `U ⊗ N_dual ⊗ L` then differentiates the whole block — every
/// multiple-reflection path included — through the same Redheffer star
/// product the forward solver uses, which is why this cannot drift from it.
///
/// The intensities follow the solver's own `finalize`: R = |r_f|² and the
/// stored T is |t_fwd|²·f_fwd, with `f_fwd` the boundary-admittance ratio;
/// RB = |r_b|² and TB = |t_b|²·f_back. Those factors depend only on the
/// block's two half-spaces, so they are constant under a film thickness
/// and ride outside the derivative.
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
    with_phase: bool,
    with_back: bool,
) -> Vec<f64> {
    let with_phase = with_phase || with_back;
    let n_ch = if with_back {
        ALL_SLOTS
    } else if with_phase {
        BACK_SLOT
    } else {
        4
    };
    let mut row = vec![0.0f64; n_ch * n_par];
    let (fields_s, fields_p) = build_stack_fields_range_dual(
        start,
        end,
        n_slice,
        &sa.thicknesses,
        &sa.rough_vals,
        &sa.rough_types,
        lam,
        nsin_fi,
    );
    for (pol, fields) in [(0i32, &fields_s), (1i32, &fields_p)] {
        let pi = pol as usize;
        let m = fields.s_left[fields.end];
        let flux = block_flux_factors(fields, pol);
        let f_fwd = flux[2];
        let r_conj = m.0.conj();
        let t_conj = m.2.conj();
        for (kp, &fi) in par_films.iter().enumerate() {
            // Film `fi` is solver slot `fi + 1` (ambient occupies slot 0).
            let j = fi + 1;
            let [dr, dt] = if with_back {
                let [dr, dtb, dt, drb] = thickness_slopes4_ddz(fields, j, lam);
                let b = BACK_SLOT;
                row[(b + pi) * n_par + kp] = 2.0 * (m.3.conj() * drb).re;
                row[(b + 2 + pi) * n_par + kp] = 2.0 * flux[1] * (m.1.conj() * dtb).re;
                row[(b + 4 + 2 * pi) * n_par + kp] = drb.re;
                row[(b + 5 + 2 * pi) * n_par + kp] = drb.im;
                row[(b + 8 + 2 * pi) * n_par + kp] = dtb.re;
                row[(b + 9 + 2 * pi) * n_par + kp] = dtb.im;
                [dr, dt]
            } else {
                thickness_slopes_ddz(fields, j, lam)
            };
            // d|z|²/dd = 2·Re(conj(z)·dz/dd).
            row[pi * n_par + kp] = 2.0 * (r_conj * dr).re;
            row[(2 + pi) * n_par + kp] = 2.0 * f_fwd * (t_conj * dt).re;
            if with_phase {
                row[(4 + 2 * pi) * n_par + kp] = dt.re;
                row[(5 + 2 * pi) * n_par + kp] = dt.im;
                row[(8 + 2 * pi) * n_par + kp] = dr.re;
                row[(9 + 2 * pi) * n_par + kp] = dr.im;
            }
        }
    }
    row
}

impl DesignContext for SmatrixContext {
    fn evaluate_merit(&self, stack: &DesignStack) -> Result<f64, String> {
        // §4.6's one branch. K = 1 is the pre-F2.2 pair of calls, reached
        // without allocating a slice to hold one simulation.
        match &self.envs {
            None => {
                let sim = self.simulate(stack)?;
                Ok(self.spec.merit(&sim, 1e6))
            }
            Some(_) => {
                let sims = self.simulate_all(stack)?;
                Ok(self.spec.merit_multi(&sims, 1e6))
            }
        }
    }

    fn simulate(&self, stack: &DesignStack) -> Result<SimCurves, String> {
        SmatrixContext::simulate(self, stack)
    }

    fn simulate_env(&self, stack: &DesignStack, env: usize) -> Result<SimCurves, String> {
        SmatrixContext::simulate_env(self, stack, env)
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

    fn environments(&self) -> Option<&CompiledEnvironments> {
        self.envs.as_deref()
    }

    fn environments_mut(&mut self) -> Option<&mut CompiledEnvironments> {
        // `make_mut` clones only when a residual closure is holding the
        // other end of the `Arc`, and insertions are rare - once per
        // needle cycle against dozens of LM iterations.
        self.envs.as_mut().map(Arc::make_mut)
    }
}

impl SmatrixContext {
    /// F2.2: one solve per environment, in roster order.
    ///
    /// `stack` is the shared design — environment 0's stack. The other
    /// environments are re-expanded from it (`CompiledEnvironments::expand`)
    /// and solved in turn. §6 limitation 6: the ×K cost is inherent, and
    /// this is where it is paid, once per eval.
    ///
    /// On the flat path this is one `simulate` in a one-element `Vec`; the
    /// callers that care about the extra allocation branch on `envs`
    /// themselves rather than going through here.
    pub fn simulate_all(&self, stack: &DesignStack) -> Result<Vec<SimCurves>, String> {
        match &self.envs {
            None => Ok(vec![self.simulate(stack)?]),
            // W6: each assembly solves on its own, so an environment with
            // an incoherent plate takes the cascade and the others do not.
            Some(envs) => envs
                .expand(stack)?
                .iter()
                .enumerate()
                .map(|(e, st)| self.simulate_env(st, e))
                .collect(),
        }
    }

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

        // F2.2: the same one branch, hoisted OUT of the closure - the
        // discriminant cannot change between iterations, so the flat path
        // pays one bool per call and keeps its exact call sequence.
        let multi = self.envs.is_some();
        let residuals = move |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            let mut st = base_stack.clone();
            apply_params(&mut st, &params_owned, x)?;
            if multi {
                // Every environment sees the same parameter vector: the
                // design object is shared, and `expand` is what writes one
                // thickness into K places. Residuals come back env-major.
                let sims = ctx_self.simulate_all(&st)?;
                spec.residuals_multi(&sims, out)
                    .map_err(|c| format!("missing curve {c:?}"))
            } else {
                let sim = ctx_self.simulate(&st)?;
                spec.residuals(&sim, out)
                    .map_err(|c| format!("missing curve {c:?}"))
            }
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
        //
        // F2.3: a joint run takes the same posture, for a different
        // reason. Elimination is the INVERSE of `insert_seed`, and the
        // compile has only the forward direction: a span that vanishes
        // from the shared design here would leave K templates still
        // carrying it, and the next `expand` would route thicknesses into
        // spans that no longer line up. So while K > 1 the floor is a hard
        // bound and the sweep below clamps up instead of removing. The
        // cost is F0.3's: a seed the optimizer wants to reject parks at
        // the floor instead of disappearing.
        let joint = self.envs.is_some();
        let lb_floor = if joint || self.thin_layer_policy == ThinLayerPolicy::ClampUpAlways {
            self.clamp_min_nm
        } else {
            0.0
        };
        let lb = vec![lb_floor; x0.len()];
        let ub = vec![self.clamp_max_nm; x0.len()];

        // The analytic Jacobian, when the spec is fully covered by it. The
        // source decides per call: a spec with a color demand declines and
        // the driver differences that iteration, so the mode is a
        // preference, not a promise.
        //
        // Which *solver* consumes them is `self.lm.backend`'s business, and
        // `run_optimizer`'s: the built-in LM is the default and the only one
        // with bounds semantics of its own (R4.4c).
        //
        // The pool is sized by the work one parallel call spreads: grid
        // points × solver layers, times the columns a differenced Jacobian
        // runs at once, times the environments a joint run simulates.
        let differenced = self.lm.jacobian == JacobianMode::Fd
            || multi
            || self.spec.sensitivity_incomplete()
            || stack.has_incoherent();
        let work = self.sin_theta.len()
            * self.wavls.len()
            * stack.solver_arrays().n_layers as usize
            * if differenced { params.len().max(1) } else { 1 }
            * self.envs.as_ref().map_or(1, |e| e.names().len());
        let res = in_lm_pool(work, || match self.lm.jacobian {
            JacobianMode::Analytic => {
                let src = DepositJacobian::new(
                    self.clone(),
                    stack.clone(),
                    params.clone(),
                    self.wavls.len(),
                );
                run_optimizer(&residuals, Some(&src), &x0, &lb, &ub, &self.lm)
            }
            JacobianMode::Fd => {
                run_optimizer(&residuals, None::<&NoJacobian>, &x0, &lb, &ub, &self.lm)
            }
        })?;

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
            if joint {
                ThinLayerPolicy::ClampUpAlways
            } else {
                self.thin_layer_policy
            },
            false,
        )?;
        self.clamp_accumulator.merge(rep);
        // The optimizer moved the primaries; the design's other copies of
        // the same films follow them.
        if let Some(envs) = &self.envs {
            envs.sync(stack)?;
        }

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
///   whose mode scales exactly (`InhMode::Fixed`, `GradientMode::FixedSpan`).
///   The RATE modes read absolute depth, so scaling one thickness moves
///   every row's `f` — those spans are scalable too since F1.7, and
///   `refresh_profiles` (at the construction points) rebuilds the profile
///   at the new depth before anything reads it: the rebuild is what makes
///   "one thickness per span" sound for every mode, not just the
///   scale-exact ones.
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

/// Threads for one optimizer run up to [`LM_BASE_WORK`]; the pool doubles at
/// every further 16× of work. Measured on 32 logical CPUs, ms per LM
/// iteration at 8 / 16 / 32 threads (work = points × layers × columns):
///
/// | case                         | work  | 8     | 16    | 32    |
/// |------------------------------|-------|-------|-------|-------|
/// | 16 films, 61 λ               | 1.1k  | 0.14  | 0.19  | 0.24  |
/// | 80 films, 61 λ               | 5.0k  | 1.36  | 1.46  | 1.66  |
/// | 16 graded films, 61 λ        | 8.4k  | 0.60  | 0.54  | 0.54  |
/// | 80 films, 401 λ              | 33k   | 6.34  | 5.81  | 6.10  |
/// | 32 graded films (1196), 61 λ | 73k   | 5.10  | 3.87  | 3.81  |
/// | 80 graded films (2990), 61 λ | 183k  | 17.3  | 22.4  | 13.8  |
/// | 16 films, 1001 λ × 30 angles | 541k  | 63.0  | 49.1  | 45.5  |
/// | 32 graded films, 201 λ × 30  | 7.3M  | 587   | 381   | 313   |
///
/// A small sweep loses to the cost of waking and coordinating threads that
/// get a few points each; a large one is throughput-bound and takes every
/// thread there is.
const LM_THREADS: usize = 8;

/// Work below which [`LM_THREADS`] is the pool: the 16-film case on 601
/// wavelengths (10.8k) sits at the crossover, a tie between 8 and 16.
const LM_BASE_WORK: usize = 8192;

/// The pool for `work`: [`LM_THREADS`], doubled per 16× of work above
/// [`LM_BASE_WORK`] (8k → 16, 131k → 32, 2.1M → 64), before the cap.
fn lm_threads_for(work: usize) -> usize {
    let mut threads = LM_THREADS;
    let mut bound = LM_BASE_WORK;
    while work >= bound && threads < 1 << 16 {
        threads *= 2;
        bound = bound.saturating_mul(16);
    }
    threads
}

/// Run `f` (one optimizer run) on a rayon worker of a small pool.
///
/// Two costs this avoids, both per parallel call and both paid hundreds of
/// times per run. Called from outside a pool, every sweep is injected into
/// the global queue while the caller sleeps on an OS latch; on a worker it
/// is pushed locally. And a big pool on a small sweep wakes threads that
/// cost more to coordinate than their share saves. The pool has
/// [`lm_threads_for`]`(work)` threads, or `NAVETTE_LM_THREADS` rounded up
/// to a power of two, capped at the global
/// pool's size (so `RAYON_NUM_THREADS` still bounds it); the pools are
/// built once and kept. Results do not depend on the thread count: every
/// sweep collects in index order.
///
/// Already on a worker (a caller that parallelizes over runs), `f` runs
/// where it is, so the caller's own parallelism is left alone.
fn in_lm_pool<R: Send>(work: usize, f: impl FnOnce() -> R + Send) -> R {
    use std::sync::OnceLock;
    const N_POOLS: usize = usize::BITS as usize;
    static POOLS: [OnceLock<Option<rayon::ThreadPool>>; N_POOLS] =
        [const { OnceLock::new() }; N_POOLS];

    if rayon::current_thread_index().is_some() {
        return f();
    }
    let global = rayon::current_num_threads();
    let want = std::env::var("NAVETTE_LM_THREADS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or_else(|| lm_threads_for(work))
        .max(1)
        .next_power_of_two();
    if want >= global {
        return rayon::scope(|_| f());
    }
    let pool = POOLS[want.trailing_zeros() as usize].get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(want)
            .thread_name(move |i| format!("navette-lm{want}-{i}"))
            .build()
            .ok()
    });
    match pool {
        Some(p) => p.install(f),
        None => rayon::scope(|_| f()),
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
/// sensitivity pass does not carry — color demands. That is
/// a property of the spec, not of `x`, so a declining source declines every
/// iteration and the run is simply the finite-difference run. The first
/// decline is remembered, so the later calls skip the deposits sweep whose
/// answer they already know.
struct DepositJacobian {
    ctx: SmatrixContext,
    base_stack: DesignStack,
    params: Vec<Param>,
    n_wav: usize,
    declined: AtomicBool,
}

impl DepositJacobian {
    fn new(ctx: SmatrixContext, base_stack: DesignStack, params: Vec<Param>, n_wav: usize) -> Self {
        Self {
            ctx,
            base_stack,
            params,
            n_wav,
            declined: AtomicBool::new(false),
        }
    }
}

impl JacobianSource for DepositJacobian {
    fn fill(&self, x: &[f64], jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
        // F2.2: a multi-environment spec is an uncovered spec for now.
        // `curve_sensitivity` walks the targets of ONE simulation and
        // knows nothing about `env_idx`; routing the deposits by design
        // -film name and summing them across environments is F2.3. Until
        // then the joint run is a finite-difference run, which is the
        // documented decline, not a silent single-environment Jacobian.
        // W6: thickness deposits are read off one coherent block; through
        // an incoherent cascade they would need the cascade's adjoint
        // weights too. A stack with an incoherent film differences instead.
        if self.ctx.envs.is_some()
            || self.declined.load(Ordering::Relaxed)
            || self.base_stack.has_incoherent()
        {
            return Ok(None);
        }
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
        // Uncovered rows (colour), or rows whose curves the sweep has no
        // deposits for: difference instead.
        if !sens.is_complete() || !dep.covers(&sens) {
            self.declined.store(true, Ordering::Relaxed);
            return Ok(None);
        }
        assemble_jacobian_mapped(&sens, &dep, self.n_wav, &self.params, jac)?;
        Ok(Some(sens.rows.len()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[allow(unused_imports)]
    use crate::smatrix::angle::{Radians, deg, degs};
    use crate::smatrix::synthesis::merit::{ConstraintKind, MeritKey, MeritTarget, SimTransform};
    use crate::smatrix::synthesis::structure::LayerSpec;

    #[test]
    fn lm_pool_doubles_per_sixteenfold_work() {
        assert_eq!(lm_threads_for(0), 8);
        assert_eq!(lm_threads_for(8191), 8);
        assert_eq!(lm_threads_for(8192), 16);
        assert_eq!(lm_threads_for(131_071), 16);
        assert_eq!(lm_threads_for(131_072), 32);
        assert_eq!(lm_threads_for(2_097_152), 64);
        assert!(lm_threads_for(usize::MAX) <= 1 << 16);
    }

    fn ar_spec(angle: f64, wl: f64) -> MeritSpec {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: deg(angle),
            curve: CurveId::Rs,
        });
        spec.add_target(MeritTarget {
            key_idx: k as u32,
            env_idx: 0,
            wavelengths: vec![wl].into(),
            kind: ConstraintKind::Exact,
            transform: SimTransform::Linear,
            norm_factor: 1.0,
            normalized_targets: vec![0.0].into(), // R = 0 demanded
            tolerances: vec![0.01].into(),
            band: vec![].into(),
            phase: false,
            differential_passes: None,
            retardance: false,
            point_weights: Vec::new().into(),
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
            envs: None,
            thin_layer_policy: crate::smatrix::synthesis::config::ThinLayerPolicy::Remove,
        }
    }

    /// Context whose spec demands PDts (Ts, phase, 1 pass): `simulate()`
    /// must then fill complex-t rows + reference metadata (gated assembly).
    fn ar_ctx_pd() -> SmatrixContext {
        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: Radians::ZERO,
            curve: CurveId::Ts,
        });
        spec.add_target(MeritTarget {
            key_idx: k as u32,
            env_idx: 0,
            wavelengths: vec![1000.0].into(),
            kind: ConstraintKind::Exact,
            transform: SimTransform::Phase,
            norm_factor: 1.0,
            normalized_targets: vec![0.0].into(),
            tolerances: vec![0.05].into(),
            band: vec![].into(),
            phase: true,
            differential_passes: Some(1.0),
            retardance: false,
            point_weights: Vec::new().into(),
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
            envs: None,
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
            angle: Radians::ZERO,
            curve: CurveId::Ts,
        });
        pd_spec
            .add_target(MeritTarget {
                key_idx: pk as u32,
                env_idx: 0,
                wavelengths: vec![400.0].into(),
                kind: ConstraintKind::Exact,
                transform: SimTransform::Phase,
                norm_factor: 1.0,
                normalized_targets: vec![0.0].into(),
                tolerances: vec![0.05].into(),
                band: vec![].into(),
                phase: true,
                differential_passes: None,
                retardance: false,
                point_weights: Vec::new().into(),
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
            envs: None,
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
        let r = crate::smatrix::optics_core::reference_phase(400.0, 1.0, Radians::ZERO, 500.0, 1.0);
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
        assert_eq!(sim.n_front_re.as_ref(), &[1.0_f64]);
        assert_eq!(sim.n_back_re.as_ref(), &[1.0_f64]);
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
        // PD1: per-λ columns, one entry per wavelength (air / 1.52 are
        // constant here, so every entry repeats).
        assert_eq!(sim.n_front_re.as_ref(), &[1.0_f64, 1.0, 1.0]);
        assert_eq!(sim.n_back_re.as_ref(), &[1.52_f64, 1.52, 1.52]);
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
            angle: Radians::ZERO,
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
            env_idx: 0,
            wavelengths: wl.into(),
            kind: ConstraintKind::Exact,
            transform: SimTransform::Phase,
            norm_factor: 1.0,
            normalized_targets: tgt.into(),
            tolerances: vec![0.05, 0.05, 0.05].into(),
            band: vec![].into(),
            phase: true,
            differential_passes: Some(1.0),
            retardance: false,
            point_weights: Vec::new().into(),
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
            angle: Radians::ZERO,
            curve: CurveId::Ts,
        });
        spec.add_target(MeritTarget {
            key_idx: k as u32,
            env_idx: 0,
            wavelengths: vec![1000.0].into(),
            kind: ConstraintKind::Exact,
            transform: SimTransform::Phase,
            norm_factor: 1.0,
            normalized_targets: vec![qw_phase - ref_qw].into(),
            tolerances: vec![0.05].into(),
            band: vec![].into(),
            phase: true,
            differential_passes: Some(1.0),
            retardance: false,
            point_weights: Vec::new().into(),
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
        #[allow(unused_imports)]
        use crate::smatrix::angle::{Radians, deg, degs};
        use crate::smatrix::synthesis::jacobian::{
            ALL_SLOTS, BACK_DEPOSIT_CHANNELS, BACK_SLOT, CurveDeposits, DEPOSIT_CHANNELS,
            assemble_jacobian, deposit_channel,
        };

        const WLS: [f64; 5] = [480.0, 520.0, 560.0, 600.0, 640.0];
        const SINS: [f64; 2] = [0.0, 0.45];

        /// The incidence angle whose sine the context holds.
        fn angle_of(sin_theta: f64) -> Radians {
            Radians::from_sin(sin_theta).unwrap()
        }

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
                envs: None,
                thin_layer_policy: crate::smatrix::synthesis::config::ThinLayerPolicy::Remove,
            }
        }

        /// A pointwise Exact demand on `curve` over the whole grid at `angle`.
        fn demand(spec: &mut MeritSpec, angle: Radians, curve: CurveId, target: f64) {
            demand_at(spec, angle, curve, vec![target; WLS.len()]);
        }

        /// A phase demand on `curve` at `angle` over `wls`, target 0 rad;
        /// differential with `passes`, or absolute.
        fn phase_demand(
            spec: &mut MeritSpec,
            angle: Radians,
            curve: CurveId,
            wls: Vec<f64>,
            passes: Option<f64>,
            integral: bool,
        ) {
            let n = wls.len();
            let k = spec.add_key(MeritKey { angle, curve });
            spec.add_target(MeritTarget {
                key_idx: k as u32,
                env_idx: 0,
                wavelengths: wls.into(),
                kind: ConstraintKind::Exact,
                transform: SimTransform::Phase,
                norm_factor: 1.0,
                normalized_targets: vec![0.0; n].into(),
                tolerances: vec![0.05; n].into(),
                band: vec![].into(),
                phase: true,
                differential_passes: passes,
                retardance: false,
                point_weights: Vec::new().into(),
                integral,
                weight: 1.0,
                count_norm: None,
            })
            .unwrap();
        }

        /// An Rs demand plus a color demand: the one kind of spec the
        /// sensitivity pass still cannot cover.
        fn color_spec() -> MeritSpec {
            use crate::smatrix::synthesis::color_merit::{
                ColorDemand, ColorDistance, ColorQuantity, ColorReference, E313_CX_D65_10,
                E313_CZ_D65_10,
            };
            let mut spec = MeritSpec::new();
            demand(&mut spec, angle_of(SINS[0]), CurveId::Rs, 0.0);
            let k = spec.add_key(MeritKey {
                angle: angle_of(SINS[0]),
                curve: CurveId::Rs,
            });
            let n = WLS.len();
            spec.add_color_demand(
                ColorDemand::new(
                    k as u32,
                    vec![[0.3, 0.4, 0.2]; n],
                    WLS.to_vec(),
                    vec![1.0; n],
                    WLS.to_vec(),
                    ColorQuantity::XyY,
                    ColorReference::Triple([0.33, 0.33, 0.1]),
                    ColorDistance::Channels,
                    1.0,
                    E313_CX_D65_10,
                    E313_CZ_D65_10,
                    None,
                )
                .unwrap(),
            )
            .unwrap();
            spec
        }

        /// As `demand`, with a per-wavelength target vector.
        fn demand_at(spec: &mut MeritSpec, angle: Radians, curve: CurveId, targets: Vec<f64>) {
            let k = spec.add_key(MeritKey { angle, curve });
            spec.add_target(MeritTarget {
                key_idx: k as u32,
                env_idx: 0,
                wavelengths: WLS.to_vec().into(),
                kind: ConstraintKind::Exact,
                transform: SimTransform::Linear,
                norm_factor: 1.0,
                normalized_targets: targets.into(),
                tolerances: vec![0.02; WLS.len()].into(),
                band: vec![].into(),
                phase: false,
                differential_passes: None,
                retardance: false,
                point_weights: Vec::new().into(),
                integral: false,
                weight: 1.0,
                count_norm: None,
            })
            .unwrap();
        }

        #[test]
        fn the_angle_axis_matches_the_keys_so_each_key_reads_its_own_row() {
            // An axis left in sines would send every key at a non-zero
            // angle to the last row.
            let ctx = ctx3(ar_spec(0.0, 560.0));
            let (sim, _) = ctx
                .simulate_with_deposits(&stack3(118.0, 203.0, 64.0), &[0])
                .unwrap();
            for (i, &s) in SINS.iter().enumerate() {
                assert!((sim.angles[i].get() - angle_of(s).get()).abs() < 1e-15);
                assert_eq!(sim.angle_row(angle_of(s)), i);
            }
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
            demand(&mut spec, angle_of(SINS[0]), CurveId::Rs, 0.01);
            demand(&mut spec, angle_of(SINS[1]), CurveId::Rp, 0.02);
            demand(&mut spec, angle_of(SINS[1]), CurveId::Ts, 0.95);
            demand(&mut spec, angle_of(SINS[0]), CurveId::As, 0.0);
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
            demand(&mut spec, angle_of(SINS[0]), CurveId::As, 0.0);
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
        fn phase_deposits_are_a_finite_difference_of_the_complex_rows() {
            // The amplitude half of a phase Jacobian: ∂(Re, Im) of the complex
            // forward t_s, t_p against `simulate()`'s complex rows differenced
            // in thickness. Only a phase spec collects them.
            let mut spec = MeritSpec::new();
            phase_demand(
                &mut spec,
                angle_of(SINS[0]),
                CurveId::Ts,
                WLS.to_vec(),
                None,
                false,
            );
            phase_demand(
                &mut spec,
                angle_of(SINS[1]),
                CurveId::Tp,
                WLS.to_vec(),
                None,
                false,
            );
            // Front reflection phase (W2): r_s, r_p ride in slots 8..11.
            for id in [CurveId::Rs, CurveId::Rp] {
                phase_demand(&mut spec, angle_of(SINS[1]), id, WLS.to_vec(), None, false);
            }
            let ctx = ctx3(spec);
            let d = [118.0, 203.0, 64.0];
            let stack = stack3(d[0], d[1], d[2]);
            let par = [0usize, 1, 2];
            let (_, dep) = ctx.simulate_with_deposits(&stack, &par).unwrap();
            assert_eq!(dep.n_channels(), 12);
            let (_, plain) = ctx3(ar_spec(0.0, 560.0))
                .simulate_with_deposits(&stack, &par)
                .unwrap();
            assert_eq!(
                plain.n_channels(),
                4,
                "an intensity spec grew phase deposits"
            );

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
                for (k, id) in [CurveId::Ts, CurveId::Tp, CurveId::Rs, CurveId::Rp]
                    .into_iter()
                    .enumerate()
                {
                    let ca = a.complex_curve(id).unwrap();
                    let cb = b.complex_curve(id).unwrap();
                    for pt in 0..ca.len() {
                        let fd = (ca[pt] - cb[pt]) / (2.0 * h);
                        for (ch, f) in [(4 + 2 * k, fd.re), (5 + 2 * k, fd.im)] {
                            let an = dep.get(ch, pt, kp);
                            let dev = (an - f).abs() / fd.norm().max(1e-6);
                            worst = worst.max(dev);
                            assert!(
                                dev < 1e-6,
                                "film {fi}, {id:?} ch {ch}, point {pt}: {an} vs {f}"
                            );
                            if f.abs() > 1e-6 {
                                live += 1;
                            }
                        }
                    }
                }
            }
            assert!(live > 80, "only {live} live entries");
            println!("  phase deposits: {live} live entries, worst rel dev {worst:.3e}");
        }

        #[test]
        fn the_assembled_phase_jacobian_matches_a_finite_difference_of_the_residuals() {
            // Phase rows end to end: absolute and differential (through D),
            // aligned and interpolated grids, pointwise and integral, both
            // polarizations, next to an intensity demand; plain and through
            // the optimizer's parameter map.
            let off = vec![500.0, 540.0, 580.0, 620.0];
            let mut spec = MeritSpec::new();
            demand(&mut spec, angle_of(SINS[1]), CurveId::Rs, 0.01);
            phase_demand(
                &mut spec,
                angle_of(SINS[0]),
                CurveId::Ts,
                WLS.to_vec(),
                Some(1.0),
                false,
            );
            phase_demand(
                &mut spec,
                angle_of(SINS[1]),
                CurveId::Tp,
                off.clone(),
                None,
                false,
            );
            phase_demand(
                &mut spec,
                angle_of(SINS[1]),
                CurveId::Ts,
                off,
                Some(2.0),
                true,
            );
            assert!(!spec.sensitivity_incomplete());
            let ctx = ctx3(spec);

            let d = [118.0, 203.0, 64.0];
            let stack = stack3(d[0], d[1], d[2]);
            let par = [0usize, 1, 2];
            let (sim, dep) = ctx.simulate_with_deposits(&stack, &par).unwrap();
            let sens = ctx.spec.curve_sensitivity(&sim).unwrap();
            assert!(sens.is_complete(), "uncovered {:?}", sens.uncovered);
            let mut jac = Vec::new();
            assemble_jacobian(&sens, &dep, WLS.len(), &mut jac).unwrap();
            let m = sens.rows.len();
            assert_eq!(m, WLS.len() + WLS.len() + 4 + 1);

            let src = DepositJacobian::new(
                ctx.clone(),
                stack.clone(),
                vec![Param::Row(0), Param::Row(1), Param::Row(2)],
                WLS.len(),
            );
            let mut mapped = Vec::new();
            assert_eq!(src.fill(&d, &mut mapped).unwrap(), Some(m));
            for (a, b) in jac.iter().zip(&mapped) {
                assert!((a - b).abs() <= 1e-12 * a.abs().max(1.0), "{a} vs {b}");
            }

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
            println!("  phase J: {live} live entries, worst rel dev {worst:.3e}");
        }

        /// W2: the rows an intensity-only spec never asked for stay unbuilt,
        /// and a spec that reads them gets them.
        #[test]
        fn new_rows_are_built_only_when_read() {
            let stack = stack3(118.0, 203.0, 64.0);
            let mut plain = MeritSpec::new();
            demand(&mut plain, angle_of(SINS[1]), CurveId::Rs, 0.0);
            let sim = ctx3(plain).simulate(&stack).unwrap();
            assert!(sim.curve(CurveId::Ru).is_none());
            assert!(sim.back.iter().all(Option::is_none));
            assert!(sim.cplx.iter().all(Option::is_none));
            assert!(sim.cplx_back.iter().all(Option::is_none));

            let mut spec = MeritSpec::new();
            demand(&mut spec, angle_of(SINS[1]), CurveId::Au, 0.0);
            demand(&mut spec, angle_of(SINS[1]), CurveId::ABp, 0.0);
            phase_demand(
                &mut spec,
                angle_of(SINS[1]),
                CurveId::RBs,
                WLS.to_vec(),
                None,
                false,
            );
            let sim = ctx3(spec).simulate(&stack).unwrap();
            for id in [CurveId::Ru, CurveId::Tu] {
                assert!(sim.curve(id).is_some(), "{id:?}");
            }
            for id in [CurveId::RBp, CurveId::TBp] {
                assert!(sim.back_curve(id).is_some(), "{id:?}");
            }
            assert!(sim.back_curve(CurveId::RBs).is_none());
            assert!(sim.complex_back_curve(CurveId::RBs).is_some());
            assert!(sim.complex_back_curve(CurveId::TBs).is_none());
            assert!(sim.complex_curve(CurveId::Rs).is_none());
        }

        /// W2 values: unpolarized rows are the engine's `0.5 * (s + p)`
        /// bit for bit; back transmission equals front transmission
        /// (reciprocity, any stack), back reflection differs from front
        /// reflection (the top film absorbs), and the back intensities are
        /// the moduli of the back amplitudes.
        #[test]
        fn new_rows_follow_averages_reciprocity_and_their_amplitudes() {
            let mut spec = MeritSpec::new();
            for curve in [
                CurveId::Ru,
                CurveId::Tu,
                CurveId::RBs,
                CurveId::RBp,
                CurveId::RBu,
                CurveId::TBs,
                CurveId::TBp,
                CurveId::TBu,
            ] {
                demand(&mut spec, angle_of(SINS[1]), curve, 0.0);
            }
            for curve in [CurveId::RBs, CurveId::RBp, CurveId::TBs, CurveId::TBp] {
                phase_demand(
                    &mut spec,
                    angle_of(SINS[1]),
                    curve,
                    WLS.to_vec(),
                    None,
                    false,
                );
            }
            let sim = ctx3(spec).simulate(&stack3(118.0, 203.0, 64.0)).unwrap();
            let f = |id: CurveId| -> Arc<[f64]> {
                sim.curve(id)
                    .or_else(|| sim.back_curve(id))
                    .unwrap()
                    .clone()
            };
            for (u, s, p) in [
                (CurveId::Ru, CurveId::Rs, CurveId::Rp),
                (CurveId::Tu, CurveId::Ts, CurveId::Tp),
                (CurveId::RBu, CurveId::RBs, CurveId::RBp),
                (CurveId::TBu, CurveId::TBs, CurveId::TBp),
            ] {
                let (u, s, p) = (f(u), f(s), f(p));
                for k in 0..u.len() {
                    assert_eq!(u[k].to_bits(), (0.5 * (s[k] + p[k])).to_bits());
                }
            }
            let mut r_gap = 0.0f64;
            for (fr, bk, rf, rb) in [
                (CurveId::Ts, CurveId::TBs, CurveId::Rs, CurveId::RBs),
                (CurveId::Tp, CurveId::TBp, CurveId::Rp, CurveId::RBp),
            ] {
                let (t, tb, r, rbk) = (f(fr), f(bk), f(rf), f(rb));
                for k in 0..t.len() {
                    assert!(
                        (t[k] - tb[k]).abs() < 1e-12,
                        "{fr:?}[{k}]: {} vs {}",
                        t[k],
                        tb[k]
                    );
                    r_gap = r_gap.max((r[k] - rbk[k]).abs());
                }
            }
            assert!(
                r_gap > 1e-4,
                "absorbing film should break R_b = R_f ({r_gap})"
            );
            let rb_c = sim.complex_back_curve(CurveId::RBs).unwrap();
            let rb = f(CurveId::RBs);
            for k in 0..rb.len() {
                assert_eq!(rb[k].to_bits(), rb_c[k].norm_sqr().to_bits());
            }
        }

        /// W2: incidence from inside an absorbing exit medium is refused at
        /// the first solve, over the back demands' wavelengths only.
        #[test]
        fn back_demands_refuse_an_absorbing_exit_medium() {
            let nw = WLS.len();
            let absorbing_sub = |k: f64| {
                let mut ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
                ambient.optimize = false;
                ambient.needle = false;
                let mut substrate = LayerSpec::constant("metal", 0.2, k, 0.0, nw);
                substrate.optimize = false;
                substrate.needle = false;
                DesignStack::with_films(
                    ambient,
                    substrate,
                    vec![LayerSpec::constant("L", 1.46, 0.0, 90.0, nw)],
                )
                .unwrap()
            };
            let mut back = MeritSpec::new();
            demand(&mut back, angle_of(SINS[0]), CurveId::TBs, 0.0);
            let err = ctx3(back).simulate(&absorbing_sub(3.0)).unwrap_err();
            assert!(err.contains("exit medium absorbs"), "{err}");

            // A front demand on the same stack is fine...
            let mut front = MeritSpec::new();
            demand(&mut front, angle_of(SINS[0]), CurveId::Rs, 0.0);
            ctx3(front).simulate(&absorbing_sub(3.0)).unwrap();
            // ...and so is a back demand on a nearly transparent exit.
            let mut back = MeritSpec::new();
            demand(&mut back, angle_of(SINS[0]), CurveId::RBp, 0.0);
            ctx3(back).simulate(&absorbing_sub(1e-6)).unwrap();
        }

        /// W2 analytic Jacobian: unpolarized rows split onto their s/p
        /// deposits, absorption on unpolarized companions, and front
        /// reflection phase on the new r_f slots — all against a finite
        /// difference of the residuals.
        #[test]
        fn the_jacobian_covers_unpolarized_rows_and_reflection_phase() {
            let mut spec = MeritSpec::new();
            demand(&mut spec, angle_of(SINS[1]), CurveId::Ru, 0.02);
            demand(&mut spec, angle_of(SINS[0]), CurveId::Tu, 0.9);
            demand(&mut spec, angle_of(SINS[1]), CurveId::Au, 0.01);
            phase_demand(
                &mut spec,
                angle_of(SINS[1]),
                CurveId::Rs,
                WLS.to_vec(),
                None,
                false,
            );
            phase_demand(
                &mut spec,
                angle_of(SINS[1]),
                CurveId::Rp,
                vec![500.0, 540.0, 580.0, 620.0],
                None,
                false,
            );
            assert!(!spec.sensitivity_incomplete());
            let ctx = ctx3(spec);
            let d = [118.0, 203.0, 64.0];
            let stack = stack3(d[0], d[1], d[2]);
            let par = [0usize, 1, 2];
            let (sim, dep) = ctx.simulate_with_deposits(&stack, &par).unwrap();
            let sens = ctx.spec.curve_sensitivity(&sim).unwrap();
            assert!(sens.is_complete() && dep.covers(&sens));
            let mut jac = Vec::new();
            assemble_jacobian(&sens, &dep, WLS.len(), &mut jac).unwrap();
            let m = sens.rows.len();

            let h = 1e-4;
            let mut live = 0usize;
            for (kp, &fi) in par.iter().enumerate() {
                let mut sp = stack.clone();
                sp.set_thickness(fi, d[fi] + h).unwrap();
                let mut sm = stack.clone();
                sm.set_thickness(fi, d[fi] - h).unwrap();
                let (mut rp, mut rm) = (Vec::new(), Vec::new());
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
                    assert!(
                        (an - fd[i]).abs() / scale < 1e-6,
                        "row {i}, film {fi}: {an} vs {}",
                        fd[i]
                    );
                    if fd[i].abs() > 1e-6 {
                        live += 1;
                    }
                }
            }
            assert!(live > 30, "only {live} live entries");
        }

        /// W2 step 4: the back channels' deposits — ∂RB, ∂TB and ∂(Re, Im)
        /// of r_b, t_b — against `simulate()`'s back rows differenced in
        /// thickness. The front slots of the same sweep stay bitwise those
        /// of a sweep without back channels.
        #[test]
        fn back_deposits_are_a_finite_difference_of_the_back_rows() {
            let mut spec = MeritSpec::new();
            for (i, curve) in BACK_DEPOSIT_CHANNELS.into_iter().enumerate() {
                let angle = angle_of(SINS[i % 2]);
                demand(&mut spec, angle, curve, 0.0);
                phase_demand(&mut spec, angle, curve, WLS.to_vec(), None, false);
            }
            demand(&mut spec, angle_of(SINS[1]), CurveId::Rs, 0.0);
            let ctx = ctx3(spec);
            let d = [118.0, 203.0, 64.0];
            let stack = stack3(d[0], d[1], d[2]);
            let par = [0usize, 1, 2];
            let (_, dep) = ctx.simulate_with_deposits(&stack, &par).unwrap();
            assert_eq!(dep.n_channels(), ALL_SLOTS);

            let mut front = MeritSpec::new();
            demand(&mut front, angle_of(SINS[1]), CurveId::Rs, 0.0);
            let (_, plain) = ctx3(front).simulate_with_deposits(&stack, &par).unwrap();
            for ch in 0..plain.n_channels() {
                for pt in 0..dep.n_points() {
                    for kp in 0..par.len() {
                        assert_eq!(
                            dep.get(ch, pt, kp).to_bits(),
                            plain.get(ch, pt, kp).to_bits(),
                            "front slot {ch} moved with the back channels on"
                        );
                    }
                }
            }

            let h = 1e-4;
            let mut worst = 0.0f64;
            let mut live = 0usize;
            let mut check = |an: f64, f: f64, scale: f64, what: &str| {
                let dev = (an - f).abs() / scale;
                worst = worst.max(dev);
                assert!(dev < 1e-6, "{what}: {an} vs {f}");
                if f.abs() > 1e-6 {
                    live += 1;
                }
            };
            for (kp, &fi) in par.iter().enumerate() {
                let mut sp = stack.clone();
                sp.set_thickness(fi, d[fi] + h).unwrap();
                let mut sm = stack.clone();
                sm.set_thickness(fi, d[fi] - h).unwrap();
                let a = ctx.simulate(&sp).unwrap();
                let b = ctx.simulate(&sm).unwrap();
                for (k, id) in BACK_DEPOSIT_CHANNELS.into_iter().enumerate() {
                    let (ia, ib) = (a.back_curve(id).unwrap(), b.back_curve(id).unwrap());
                    let fd: Vec<f64> = ia
                        .iter()
                        .zip(ib.iter())
                        .map(|(p, m)| (p - m) / (2.0 * h))
                        .collect();
                    let scale = fd.iter().fold(1e-6f64, |acc, v| acc.max(v.abs()));
                    let ch = deposit_channel(id).unwrap();
                    assert_eq!(ch, BACK_SLOT + k);
                    for (pt, &f) in fd.iter().enumerate() {
                        check(
                            dep.get(ch, pt, kp),
                            f,
                            scale,
                            &format!("film {fi}, {id:?}, point {pt}"),
                        );
                    }
                    let ca = a.complex_back_curve(id).unwrap();
                    let cb = b.complex_back_curve(id).unwrap();
                    for pt in 0..ca.len() {
                        let fd = (ca[pt] - cb[pt]) / (2.0 * h);
                        let ch = BACK_SLOT + 4 + 2 * k;
                        for (c, f) in [(ch, fd.re), (ch + 1, fd.im)] {
                            check(
                                dep.get(c, pt, kp),
                                f,
                                fd.norm().max(1e-6),
                                &format!("film {fi}, {id:?} slot {c}, point {pt}"),
                            );
                        }
                    }
                }
            }
            assert!(live > 150, "only {live} live entries");
            println!("  back deposits: {live} live entries, worst rel dev {worst:.3e}");
        }

        /// W2 step 4 end to end: back intensity (polarized, unpolarized,
        /// absorption) and back phase (absolute and differential, aligned
        /// and interpolated, pointwise and integral) rows of the analytic
        /// Jacobian against a finite difference of the residuals, and the
        /// optimizer's source answers instead of declining.
        #[test]
        fn the_jacobian_covers_back_demands() {
            let off = vec![500.0, 540.0, 580.0, 620.0];
            let mut spec = MeritSpec::new();
            demand(&mut spec, angle_of(SINS[1]), CurveId::RBu, 0.1);
            demand(&mut spec, angle_of(SINS[0]), CurveId::TBp, 0.9);
            demand(&mut spec, angle_of(SINS[1]), CurveId::ABs, 0.05);
            demand(&mut spec, angle_of(SINS[0]), CurveId::Rs, 0.01);
            phase_demand(
                &mut spec,
                angle_of(SINS[1]),
                CurveId::RBs,
                WLS.to_vec(),
                None,
                false,
            );
            phase_demand(
                &mut spec,
                angle_of(SINS[0]),
                CurveId::TBp,
                off.clone(),
                Some(1.0),
                false,
            );
            phase_demand(
                &mut spec,
                angle_of(SINS[1]),
                CurveId::TBs,
                off,
                Some(2.0),
                true,
            );
            assert!(!spec.sensitivity_incomplete());
            let ctx = ctx3(spec);
            let d = [118.0, 203.0, 64.0];
            let stack = stack3(d[0], d[1], d[2]);
            let par = [0usize, 1, 2];
            let (sim, dep) = ctx.simulate_with_deposits(&stack, &par).unwrap();
            let sens = ctx.spec.curve_sensitivity(&sim).unwrap();
            assert!(sens.is_complete(), "uncovered {:?}", sens.uncovered);
            assert!(dep.covers(&sens));
            let mut jac = Vec::new();
            assemble_jacobian(&sens, &dep, WLS.len(), &mut jac).unwrap();
            let m = sens.rows.len();
            assert_eq!(m, 4 * WLS.len() + WLS.len() + 4 + 1);

            let src = DepositJacobian::new(
                ctx.clone(),
                stack.clone(),
                vec![Param::Row(0), Param::Row(1), Param::Row(2)],
                WLS.len(),
            );
            let mut mapped = Vec::new();
            assert_eq!(src.fill(&d, &mut mapped).unwrap(), Some(m));
            for (a, b) in jac.iter().zip(&mapped) {
                assert_eq!(a.to_bits(), b.to_bits(), "{a} vs {b}");
            }

            let h = 1e-4;
            let mut worst = 0.0f64;
            let mut live = 0usize;
            for (kp, &fi) in par.iter().enumerate() {
                let mut sp = stack.clone();
                sp.set_thickness(fi, d[fi] + h).unwrap();
                let mut sm = stack.clone();
                sm.set_thickness(fi, d[fi] - h).unwrap();
                let (mut rp, mut rm) = (Vec::new(), Vec::new());
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
            assert!(live > 40, "only {live} live entries");
            println!("  back J: {live} live entries, worst rel dev {worst:.3e}");
        }

        /// W8a forward: the four retardance labels, compiled with target 0
        /// and tolerance 1, leave residuals equal to the engine's
        /// `Retardance_R/T` and `Delta_R/T` on the same stack, at 0°, 45°
        /// and 70°. At 0° Ret_R is identically 0 (s ≡ p).
        #[test]
        fn retardance_follows_the_engine() {
            use crate::smatrix::core_engine::{
                REQ_DELTA_R, REQ_DELTA_T, REQ_RETARD_R, REQ_RETARD_T,
            };
            use crate::smatrix::synthesis::targets::{TargetSet, compile_merit_spec};
            use num_complex::Complex64;
            let stack = stack3(118.0, 203.0, 64.0);
            let sa = stack.solver_arrays();
            let (nl, nw) = (sa.n_layers as usize, WLS.len());
            let indices: Vec<Complex64> = (0..nl * nw)
                .map(|i| {
                    let (l, w) = (i / nw, i % nw);
                    let b = (w * nl + l) * 2;
                    Complex64::new(sa.n_stack_cache[b], sa.n_stack_cache[b + 1])
                })
                .collect();
            let wrap = |d: f64| d - std::f64::consts::TAU * (d / std::f64::consts::TAU).round();
            for angle_deg in [0.0, 45.0, 70.0] {
                let (sol, _) = crate::smatrix::solver::solve_arrays(
                    &indices,
                    &sa.thicknesses,
                    &vec![false; nl],
                    &sa.rough_types,
                    &sa.rough_vals,
                    &WLS,
                    &[deg(angle_deg)],
                    REQ_RETARD_R | REQ_RETARD_T | REQ_DELTA_R | REQ_DELTA_T,
                    2,
                )
                .unwrap();
                let engine = |key: &str| -> Vec<f64> {
                    sol.f64maps
                        .iter()
                        .find(|(k, _)| k == key)
                        .unwrap()
                        .1
                        .clone()
                };
                for (label, key) in [
                    ("RetR", "Retardance_R"),
                    ("RetT", "Retardance_T"),
                    ("DeltaR", "Delta_R"),
                    ("DeltaT", "Delta_T"),
                ] {
                    let (zeros, ones) = (vec![0.0; nw], vec![1.0; nw]);
                    let set = TargetSet::from_value(serde_json::json!({
                        "spectral": [{
                            "wavelengths": WLS, "values": zeros, "tolerances": ones,
                            "angle_deg": angle_deg, "polarization": "u",
                            "spectral": label, "phase": true,
                        }],
                    }))
                    .unwrap();
                    let mut ctx = ctx3(compile_merit_spec(&set).unwrap());
                    ctx.sin_theta = vec![deg(angle_deg).get().sin()];
                    let mut r = Vec::new();
                    ctx.spec
                        .residuals(&ctx.simulate(&stack).unwrap(), &mut r)
                        .unwrap();
                    let want = engine(key);
                    assert_eq!(r.len(), want.len());
                    if angle_deg == 70.0 {
                        assert!(want.iter().any(|e| e.abs() > 0.05), "{key} flat: {want:?}");
                    }
                    for (i, (&got, &e)) in r.iter().zip(&want).enumerate() {
                        assert!(
                            wrap(got - e).abs() < 1e-10,
                            "{label} at {angle_deg} deg, point {i}: {got} vs engine {key} {e}"
                        );
                        if angle_deg == 0.0 && label == "RetR" {
                            assert!(got.abs() < 1e-12, "RetR at normal incidence: {got}");
                        }
                    }
                }
            }
        }

        /// A retardance demand on p key `curve` at `angle` over `wls`.
        fn retardance_demand(
            spec: &mut MeritSpec,
            angle: Radians,
            curve: CurveId,
            wls: Vec<f64>,
            target: f64,
            integral: bool,
        ) {
            let n = wls.len();
            let k = spec.add_key(MeritKey { angle, curve });
            spec.add_target(MeritTarget {
                key_idx: k as u32,
                env_idx: 0,
                wavelengths: wls.into(),
                kind: ConstraintKind::Exact,
                transform: SimTransform::Phase,
                norm_factor: 1.0,
                normalized_targets: vec![target; n].into(),
                tolerances: vec![0.05; n].into(),
                band: vec![].into(),
                phase: true,
                differential_passes: None,
                retardance: true,
                point_weights: Vec::new().into(),
                integral,
                weight: 1.0,
                count_norm: None,
            })
            .unwrap();
        }

        /// W8a LM: retardance rows of the analytic Jacobian (front R and T,
        /// back T; aligned and interpolated, pointwise and integral, one
        /// target across the wrap) against a finite difference of the
        /// residuals, and the optimizer's source answers.
        #[test]
        fn the_jacobian_covers_retardance_demands() {
            let off = vec![500.0, 540.0, 580.0, 620.0];
            let mut spec = MeritSpec::new();
            retardance_demand(
                &mut spec,
                angle_of(SINS[1]),
                CurveId::Rp,
                WLS.to_vec(),
                0.3,
                false,
            );
            retardance_demand(
                &mut spec,
                angle_of(SINS[1]),
                CurveId::Tp,
                off.clone(),
                3.1,
                false,
            );
            retardance_demand(&mut spec, angle_of(SINS[1]), CurveId::Rp, off, -0.2, true);
            retardance_demand(
                &mut spec,
                angle_of(SINS[1]),
                CurveId::TBp,
                WLS.to_vec(),
                0.1,
                false,
            );
            assert!(!spec.sensitivity_incomplete());
            let ctx = ctx3(spec);
            let d = [118.0, 203.0, 64.0];
            let stack = stack3(d[0], d[1], d[2]);
            let par = [0usize, 1, 2];
            let (sim, dep) = ctx.simulate_with_deposits(&stack, &par).unwrap();
            let sens = ctx.spec.curve_sensitivity(&sim).unwrap();
            assert!(sens.is_complete(), "uncovered {:?}", sens.uncovered);
            assert!(dep.covers(&sens));
            let mut jac = Vec::new();
            assemble_jacobian(&sens, &dep, WLS.len(), &mut jac).unwrap();
            let m = sens.rows.len();
            assert_eq!(m, 2 * WLS.len() + 4 + 1);
            let src = DepositJacobian::new(
                ctx.clone(),
                stack.clone(),
                vec![Param::Row(0), Param::Row(1), Param::Row(2)],
                WLS.len(),
            );
            let mut mapped = Vec::new();
            assert_eq!(src.fill(&d, &mut mapped).unwrap(), Some(m));

            let h = 1e-4;
            let mut worst = 0.0f64;
            let mut live = 0usize;
            for (kp, &fi) in par.iter().enumerate() {
                let mut sp = stack.clone();
                sp.set_thickness(fi, d[fi] + h).unwrap();
                let mut sm = stack.clone();
                sm.set_thickness(fi, d[fi] - h).unwrap();
                let (mut rp, mut rm) = (Vec::new(), Vec::new());
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
            assert!(live > 25, "only {live} live entries");
            println!("  retardance J: {live} live entries, worst rel dev {worst:.3e}");
        }

        #[test]
        fn a_color_demand_refuses_the_analytic_path_instead_of_zeroing_it() {
            // The fallback contract: an uncovered row must stop the assembly, not
            // come back as a row of zeros that the LM would read as "flat".
            let ctx = ctx3(color_spec());
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
                angle: Radians::ZERO,
                curve: CurveId::Rs,
            });
            spec.add_target(MeritTarget {
                key_idx: k as u32,
                env_idx: 0,
                wavelengths: vec![900.0].into(),
                kind: ConstraintKind::Exact,
                transform: SimTransform::Linear,
                norm_factor: 1.0,
                normalized_targets: vec![r_bare].into(),
                tolerances: vec![0.01].into(),
                band: vec![].into(),
                phase: false,
                differential_passes: None,
                retardance: false,
                point_weights: Vec::new().into(),
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
                envs: None,
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
                envs: None,
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
                angle: Radians::ZERO,
                curve: CurveId::Rs,
            });
            spec.add_target(MeritTarget {
                key_idx: k as u32,
                env_idx: 0,
                wavelengths: vec![900.0].into(),
                kind: ConstraintKind::Exact,
                transform: SimTransform::Linear,
                norm_factor: 1.0,
                normalized_targets: vec![r_bare].into(),
                tolerances: vec![0.01].into(),
                band: vec![].into(),
                phase: false,
                differential_passes: None,
                retardance: false,
                point_weights: Vec::new().into(),
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
                envs: None,
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
                envs: None,
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
            stack.set_row_optimize(0, false);
            ctx.optimize_thicknesses(&mut stack).unwrap();
            assert!(ctx.take_clamp_report().is_none());

            ctx.clamp_accumulator
                .merge(crate::smatrix::synthesis::structure::ClampReport {
                    spans_removed: vec!["H (1.0 nm)".to_string()],
                    spans_capped: 0,
                    rows_removed: 1,
                    edges_held: Vec::new(),
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
            // whose channels are all deposited answers with a Jacobian, a phase
            // spec included; one with a color demand answers `None`, which is
            // what makes the run fall back to differences instead of optimizing
            // against zeros.
            let stack = stack3(118.0, 203.0, 64.0);
            let x = [118.0, 203.0, 64.0];

            let mut spec = MeritSpec::new();
            demand(&mut spec, angle_of(SINS[0]), CurveId::Rs, 0.0);
            demand(&mut spec, angle_of(SINS[1]), CurveId::As, 0.0);
            let covered = DepositJacobian::new(
                ctx3(spec),
                stack.clone(),
                vec![Param::Row(0), Param::Row(1), Param::Row(2)],
                WLS.len(),
            );
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

            let mut pspec = MeritSpec::new();
            phase_demand(
                &mut pspec,
                angle_of(SINS[0]),
                CurveId::Ts,
                WLS.to_vec(),
                Some(1.0),
                false,
            );
            let phase =
                DepositJacobian::new(ctx3(pspec), stack.clone(), vec![Param::Row(0)], WLS.len());
            assert_eq!(phase.fill(&[118.0], &mut jac).unwrap(), Some(WLS.len()));
            assert!(
                jac.iter().any(|v| v.abs() > 1e-9),
                "a phase spec gave a flat J"
            );

            let color =
                DepositJacobian::new(ctx3(color_spec()), stack, vec![Param::Row(0)], WLS.len());
            assert!(
                color.fill(&[118.0], &mut jac).unwrap().is_none(),
                "a color spec must decline, not answer"
            );
            assert!(color.declined.load(Ordering::Relaxed));
            // Remembered: a later x declines too (without the sweep), and a
            // covered source never latches.
            assert!(color.fill(&[140.0], &mut jac).unwrap().is_none());
            assert!(!covered.declined.load(Ordering::Relaxed));
            assert!(!phase.declined.load(Ordering::Relaxed));
            // The up-front predicate the pool sizing reads agrees with all three.
            assert!(color.ctx.spec.sensitivity_incomplete());
            assert!(!phase.ctx.spec.sensitivity_incomplete());
            assert!(!covered.ctx.spec.sensitivity_incomplete());
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
            demand_at(
                &mut spec,
                angle_of(SINS[0]),
                CurveId::Rs,
                rs0[..nw].to_vec(),
            );
            demand_at(
                &mut spec,
                angle_of(SINS[1]),
                CurveId::Rp,
                rp0[nw..].to_vec(),
            );

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

        /// The same comparison on back demands (W2 step 4): RB and TB
        /// intensities taken from the reference stack, the analytic run on
        /// the back deposits, the other differenced.
        #[test]
        fn the_two_jacobian_modes_agree_on_back_demands() {
            let reference = stack3(118.0, 203.0, 64.0);
            let mut probe_spec = MeritSpec::new();
            demand(&mut probe_spec, angle_of(SINS[0]), CurveId::RBs, 0.0);
            demand(&mut probe_spec, angle_of(SINS[1]), CurveId::TBp, 0.0);
            let sim0 = ctx3(probe_spec).simulate(&reference).unwrap();
            let rb0 = sim0.back_curve(CurveId::RBs).unwrap().to_vec();
            let tb0 = sim0.back_curve(CurveId::TBp).unwrap().to_vec();
            let nw = WLS.len();
            let mut spec = MeritSpec::new();
            demand_at(
                &mut spec,
                angle_of(SINS[0]),
                CurveId::RBs,
                rb0[..nw].to_vec(),
            );
            demand_at(
                &mut spec,
                angle_of(SINS[1]),
                CurveId::TBp,
                tb0[nw..].to_vec(),
            );
            assert!(!spec.sensitivity_incomplete());

            let mut analytic = ctx3(spec.clone());
            analytic.lm.jacobian = JacobianMode::Analytic;
            let mut differenced = ctx3(spec);
            differenced.lm.jacobian = JacobianMode::Fd;
            analytic.clamp_min_nm = 1e-9;
            differenced.clamp_min_nm = 1e-9;

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
            println!("  back modes agree: merit {ma:.6e} vs {mb:.6e}");
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
            "staleness {drift} < half the excursion response {excursion} (m0 {m0}, \
             stale {m_stale}, fresh {m_fresh})"
        );
        // ...and under 5% of the merit it rides on.
        assert!(
            drift / m_stale.abs().max(1e-12) < 0.05,
            "the within-cycle staleness stays under the documented 5% (drift {drift}, \
             stale {m_stale}, fresh {m_fresh})"
        );
        let _ = m_fresh;
    }

    /// Where one LM iteration's time goes, part by part (a profile, not a
    /// check: `--ignored --nocapture` in a release build). The problem is
    /// `bench_lm_breakdown.py`'s: an `n_films` H/L AR coating on 61
    /// wavelengths at normal incidence, Rs and Rp demanded zero.
    /// `NAVETTE_PROFILE_FILMS` sets the film count (default 16).
    #[test]
    #[ignore]
    fn profile_lm_parts() {
        use std::time::Instant;
        let n_films: usize = std::env::var("NAVETTE_PROFILE_FILMS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(16);
        let nw = 61;
        let wl: Vec<f64> = (0..nw).map(|i| 400.0 + 5.0 * i as f64).collect();
        let mut spec = MeritSpec::new();
        for curve in [CurveId::Rs, CurveId::Rp] {
            let k = spec.add_key(MeritKey {
                angle: Radians::ZERO,
                curve,
            });
            spec.add_target(MeritTarget {
                key_idx: k as u32,
                env_idx: 0,
                wavelengths: wl.clone().into(),
                kind: ConstraintKind::Exact,
                transform: SimTransform::Linear,
                norm_factor: 1.0,
                normalized_targets: vec![0.0; nw].into(),
                tolerances: vec![0.002; nw].into(),
                band: vec![].into(),
                phase: false,
                differential_passes: None,
                retardance: false,
                point_weights: Vec::new().into(),
                integral: false,
                weight: 1.0,
                count_norm: None,
            })
            .unwrap();
        }
        let mut ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
        ambient.optimize = false;
        ambient.needle = false;
        let mut substrate = LayerSpec::constant("sub", 1.52, 0.0, 0.0, nw);
        substrate.optimize = false;
        substrate.needle = false;
        let films = (0..n_films)
            .map(|i| {
                let (name, n) = if i % 2 == 0 { ("H", 2.35) } else { ("L", 1.46) };
                let d = 550.0 / (4.0 * n) * (1.0 + 0.3 * (1.7 * i as f64 + 0.4).sin());
                LayerSpec::constant(name, n, 0.0, d, nw)
            })
            .collect();
        let stack = DesignStack::with_films(ambient, substrate, films).unwrap();
        let ctx = SmatrixContext {
            wavls: wl.clone(),
            sin_theta: vec![0.0],
            spec: spec.clone(),
            clamp_min_nm: 2.0,
            clamp_max_nm: 1000.0,
            lm: LmConfig::default(),
            clamp_accumulator: ClampReport::default(),
            envs: None,
            thin_layer_policy: ThinLayerPolicy::Remove,
        };
        let params = build_params(&stack);
        let x: Vec<f64> = params.iter().map(|p| p.value0(&stack)).collect();
        let rows: Vec<usize> = (0..n_films).collect();
        let src = DepositJacobian::new(ctx.clone(), stack.clone(), params.clone(), nw);

        fn time<T>(reps: usize, mut f: impl FnMut() -> T) -> f64 {
            for _ in 0..reps / 10 + 1 {
                std::hint::black_box(f());
            }
            let t0 = Instant::now();
            for _ in 0..reps {
                std::hint::black_box(f());
            }
            t0.elapsed().as_secs_f64() / reps as f64 * 1e6
        }
        let reps = 2000;
        let sim = ctx.simulate(&stack).unwrap();
        let (sim_d, dep) = ctx.simulate_with_deposits(&stack, &rows).unwrap();
        let sens = spec.curve_sensitivity(&sim_d).unwrap();
        let mut out = Vec::new();
        let mut jac = Vec::new();
        let parts = [
            (
                "clone + apply_params",
                time(reps, || {
                    let mut st = stack.clone();
                    apply_params(&mut st, &params, &x).unwrap();
                    st
                }),
            ),
            ("solver_arrays", time(reps, || stack.solver_arrays())),
            ("simulate", time(reps, || ctx.simulate(&stack).unwrap())),
            (
                "residuals",
                time(reps, || {
                    out.clear();
                    spec.residuals(&sim, &mut out).unwrap();
                }),
            ),
            (
                "simulate_with_deposits",
                time(reps, || ctx.simulate_with_deposits(&stack, &rows).unwrap()),
            ),
            (
                "curve_sensitivity",
                time(reps, || spec.curve_sensitivity(&sim_d).unwrap()),
            ),
            (
                "assemble_jacobian_mapped",
                time(reps, || {
                    assemble_jacobian_mapped(&sens, &dep, nw, &params, &mut jac).unwrap()
                }),
            ),
            (
                "DepositJacobian::fill",
                time(reps, || src.fill(&x, &mut jac).unwrap()),
            ),
        ];
        // One spectral point of the deposit sweep, piece by piece.
        let sa = stack.solver_arrays();
        let nl = sa.n_layers as usize;
        let lam = wl[nw / 2];
        let base = (nw / 2) * nl * 2;
        let n_slice: Vec<Complex64> = (0..nl)
            .map(|l| {
                Complex64::new(
                    sa.n_stack_cache[base + l * 2],
                    sa.n_stack_cache[base + l * 2 + 1],
                )
            })
            .collect();
        let inv_n: Vec<Complex64> = n_slice.iter().map(|&n| 1.0 / n).collect();
        let nsin = n_slice[0] * Complex64::new(0.0, 0.0);
        use crate::smatrix::needle_operator::build_stack_fields_range;
        let fields = build_stack_fields_range(
            0,
            nl - 1,
            &n_slice,
            &sa.thicknesses,
            &sa.rough_vals,
            &sa.rough_types,
            lam,
            nsin,
            0,
        );
        let point = [
            (
                "point: dual solve",
                time(reps * 10, || {
                    solve_coherent_block_fields_dual(
                        0,
                        nl - 1,
                        &n_slice,
                        &inv_n,
                        &sa.thicknesses,
                        &sa.rough_vals,
                        &sa.rough_types,
                        lam,
                        nsin,
                    )
                }),
            ),
            (
                "point: fields build (1 pol)",
                time(reps * 10, || {
                    build_stack_fields_range(
                        0,
                        nl - 1,
                        &n_slice,
                        &sa.thicknesses,
                        &sa.rough_vals,
                        &sa.rough_types,
                        lam,
                        nsin,
                        0,
                    )
                }),
            ),
            (
                "point: fields build (s + p)",
                time(reps * 10, || {
                    build_stack_fields_range_dual(
                        0,
                        nl - 1,
                        &n_slice,
                        &sa.thicknesses,
                        &sa.rough_vals,
                        &sa.rough_types,
                        lam,
                        nsin,
                    )
                }),
            ),
            (
                "point: slopes, all films (1 pol)",
                time(reps * 10, || {
                    (1..nl - 1)
                        .map(|j| thickness_slopes_ddz(&fields, j, lam)[0])
                        .sum::<Complex64>()
                }),
            ),
            (
                "point: deposit_row (2 pols)",
                time(reps * 10, || {
                    deposit_row(
                        &n_slice,
                        &sa,
                        0,
                        nl - 1,
                        lam,
                        nsin,
                        &rows,
                        n_films,
                        false,
                        false,
                    )
                }),
            ),
        ];
        let mut c = ctx.clone();
        let lm_us = time(5, || {
            let mut st = stack.clone();
            c.optimize_thicknesses_report(&mut st).unwrap()
        });
        let mut st = stack.clone();
        let (_, rep) = c.optimize_thicknesses_report(&mut st).unwrap();
        let rep = rep.unwrap();
        eprintln!(
            "--- {n_films} films, {nw} wavelengths, {} threads ---",
            rayon::current_num_threads()
        );
        for (name, us) in parts.into_iter().chain(point) {
            eprintln!("  {name:32} {us:10.3} us");
        }
        eprintln!(
            "  LM: {lm_us:.0} us, {} iterations, {} evals, {} jacobians -> {:.1} us/iteration",
            rep.iterations,
            rep.evals,
            rep.analytic_jacobians,
            lm_us / rep.iterations.max(1) as f64
        );
    }
}

/// W1's oracle: the needle P(z) the fold drives is the slope of the merit
/// the LM minimises, `dM/dδ` for a needle of thickness δ inserted at z,
/// for every demand type, at oblique incidence too. The reference inserts
/// a real seed and differences the merit (one-sided Richardson: the seed
/// cannot be thinner than zero).
#[cfg(test)]
mod needle_oracle {
    use super::*;
    use crate::smatrix::angle::{Radians, deg};
    use crate::smatrix::synthesis::color_merit::{
        ColorDemand, ColorDistance, ColorQuantity, ColorReference, E313_CX_D65_10, E313_CZ_D65_10,
    };
    use crate::smatrix::synthesis::merit::{ConstraintKind, MeritKey, MeritTarget, SimTransform};
    use crate::smatrix::synthesis::needle_pass::{
        NeedlePassInput, build_needle_targets, build_scan_sites, needle_pass_scan,
    };
    use crate::smatrix::synthesis::structure::LayerSpec;
    use num_complex::Complex64;

    const WLS: [f64; 5] = [480.0, 520.0, 560.0, 600.0, 640.0];
    const NEEDLE_N: f64 = 1.38;

    /// air | H(2.30) | L(1.46) | Habs(2.30 + 0.004i) | glass(1.52).
    fn stack() -> DesignStack {
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
                LayerSpec::constant("H", 2.30, 0.0, 118.0, nw),
                LayerSpec::constant("L", 1.46, 0.0, 203.0, nw),
                LayerSpec::constant("Habs", 2.30, 0.004, 64.0, nw),
            ],
        )
        .unwrap()
    }

    fn ctx(spec: MeritSpec, angle: Radians) -> SmatrixContext {
        SmatrixContext {
            wavls: WLS.to_vec(),
            sin_theta: vec![angle.get().sin()],
            spec,
            clamp_min_nm: 2.0,
            clamp_max_nm: 1000.0,
            lm: LmConfig::default(),
            clamp_accumulator: ClampReport::default(),
            envs: None,
            thin_layer_policy: crate::smatrix::synthesis::config::ThinLayerPolicy::Remove,
        }
    }

    /// One demand row: `curve` at `angle` over the grid.
    #[derive(Clone, Copy)]
    struct D {
        curve: CurveId,
        target: f64,
        kind: ConstraintKind,
        band: f64,
        phase: bool,
        passes: Option<f64>,
        integral: bool,
        retardance: bool,
        /// Non-uniform W4 point weights (0.5, 1.5, 2.5, ...) instead of none.
        weighted: bool,
        /// W5: a Log demand; `target` is then log10 of the level.
        log: bool,
    }

    fn d(curve: CurveId, target: f64) -> D {
        D {
            curve,
            target,
            kind: ConstraintKind::Exact,
            band: 0.0,
            phase: false,
            passes: None,
            integral: false,
            retardance: false,
            weighted: false,
            log: false,
        }
    }

    /// A Log demand on `curve` at the level `level` (log10 taken here).
    fn lg(curve: CurveId, level: f64) -> D {
        D {
            log: true,
            ..d(curve, level.log10())
        }
    }

    /// A retardance demand `arg(z_p · conj(z_s))` on p key `curve`.
    fn ret(curve: CurveId, target: f64) -> D {
        D {
            phase: true,
            retardance: true,
            ..d(curve, target)
        }
    }

    fn add(spec: &mut MeritSpec, angle: Radians, x: D) {
        let n = WLS.len();
        let k = spec.add_key(MeritKey {
            angle,
            curve: x.curve,
        });
        spec.add_target(MeritTarget {
            key_idx: k as u32,
            env_idx: 0,
            wavelengths: WLS.to_vec().into(),
            kind: x.kind,
            transform: if x.phase {
                SimTransform::Phase
            } else if x.log {
                SimTransform::Log
            } else {
                SimTransform::Linear
            },
            norm_factor: 1.0,
            normalized_targets: vec![x.target; n].into(),
            tolerances: vec![0.05; n].into(),
            band: if x.band > 0.0 {
                vec![x.band; n].into()
            } else {
                vec![].into()
            },
            phase: x.phase,
            differential_passes: x.passes,
            retardance: x.retardance,
            point_weights: if x.weighted {
                (0..n).map(|i| 0.5 + i as f64).collect::<Vec<_>>().into()
            } else {
                Vec::new().into()
            },
            integral: x.integral,
            weight: 1.0,
            count_norm: None,
        })
        .unwrap();
    }

    fn add_color(spec: &mut MeritSpec, angle: Radians, curve: CurveId) {
        let k = spec.add_key(MeritKey { angle, curve });
        let n = WLS.len();
        spec.add_color_demand(
            ColorDemand::new(
                k as u32,
                vec![[0.3, 0.4, 0.2]; n],
                WLS.to_vec(),
                vec![1.0; n],
                WLS.to_vec(),
                ColorQuantity::XyY,
                ColorReference::Triple([0.33, 0.33, 0.1]),
                ColorDistance::Channels,
                1.0,
                E313_CX_D65_10,
                E313_CZ_D65_10,
                None,
            )
            .unwrap(),
        )
        .unwrap();
    }

    /// Fold P and the merit slope at a few interior sites of every film.
    fn p_and_slope(spec: MeritSpec, angle: Radians) -> Vec<(f64, f64)> {
        p_and_slope_on(stack(), spec, angle)
    }

    fn p_and_slope_on(st: DesignStack, spec: MeritSpec, angle: Radians) -> Vec<(f64, f64)> {
        let c = ctx(spec, angle);
        let sim = c.simulate(&st).unwrap();
        let fold = build_needle_targets(&c.spec, &[angle], &WLS, Some(&sim)).unwrap();
        let sa = st.solver_arrays();
        let nk = vec![Complex64::new(NEEDLE_N, 0.0); WLS.len()];
        let sites: Vec<_> = build_scan_sites(st.films(), st.spans(), 37.0)
            .into_iter()
            .filter(|s| !s.boundary)
            .collect();
        let res = needle_pass_scan(
            &NeedlePassInput {
                n_stack_cache: &sa.n_stack_cache,
                thicknesses: &sa.thicknesses,
                rough_types: &sa.rough_types,
                rough_vals: &sa.rough_vals,
                n_layers: sa.n_layers as usize,
                wavls: &WLS,
                sin_theta: &c.sin_theta,
                fold: &fold,
                needle_n_per_wav: &nk,
                start_idx: 0,
                end_idx: (sa.n_layers - 1) as usize,
                incoherent_flags: &sa.incoherent_flags,
                calc_s: true,
                calc_p: true,
            },
            &sites,
        )
        .unwrap();
        let m0 = c.evaluate_merit(&st).unwrap();
        let merit_with = |site: &crate::smatrix::synthesis::needle_pass::ScanSite, delta: f64| {
            let mut s = st.clone();
            let seed = LayerSpec::constant("N", NEEDLE_N, 0.0, delta, WLS.len());
            s.insert_needle_seed(site.film_idx, site.depth_into_layer_nm, seed)
                .unwrap();
            c.evaluate_merit(&s).unwrap()
        };
        let h = 1e-4;
        sites
            .iter()
            .zip(&res.p_profile)
            .map(|(site, &p)| {
                let slope =
                    (4.0 * merit_with(site, h) - merit_with(site, 2.0 * h) - 3.0 * m0) / (2.0 * h);
                // The differential-phase shift is uniform in z (see
                // `PolFold::phi_gain_shift`); the cycle adds it the same way.
                (p + fold.phi_gain_shift(), slope)
            })
            .collect()
    }

    fn check(name: &str, spec: MeritSpec, angle: Radians) {
        check_on(name, stack(), spec, angle);
    }

    fn check_on(name: &str, st: DesignStack, spec: MeritSpec, angle: Radians) {
        let pairs = p_and_slope_on(st, spec, angle);
        let scale = pairs.iter().map(|&(_, s)| s.abs()).fold(0.0, f64::max);
        assert!(scale > 1e-6, "{name}: the merit does not move");
        for (i, &(p, s)) in pairs.iter().enumerate() {
            assert!(
                (p - s).abs() <= 1e-5 * scale,
                "{name} at {:.0} deg, site {i}: P {p} vs merit slope {s}",
                angle.get().to_degrees()
            );
        }
    }

    fn cases(angle: Radians) -> Vec<(&'static str, MeritSpec)> {
        use ConstraintKind::*;
        use CurveId::*;
        let one = |x: D| {
            let mut s = MeritSpec::new();
            add(&mut s, angle, x);
            s
        };
        let mut out = vec![
            ("Rs", one(d(Rs, 0.3))),
            ("Rp", one(d(Rp, 0.3))),
            ("Ru", one(d(Ru, 0.3))),
            ("Ts", one(d(Ts, 0.6))),
            ("Tp", one(d(Tp, 0.6))),
            ("Tu", one(d(Tu, 0.6))),
            ("As", one(d(As, 0.05))),
            ("Au", one(d(Au, 0.05))),
            ("RBs", one(d(RBs, 0.3))),
            ("TBp", one(d(TBp, 0.6))),
            ("RBu", one(d(RBu, 0.3))),
            ("TBu", one(d(TBu, 0.6))),
            ("ABu", one(d(ABu, 0.05))),
            (
                "Rs phase",
                one(D {
                    phase: true,
                    ..d(Rs, 0.5)
                }),
            ),
            (
                "Tp phase",
                one(D {
                    phase: true,
                    ..d(Tp, 0.5)
                }),
            ),
            (
                "PDts",
                one(D {
                    phase: true,
                    passes: Some(1.0),
                    ..d(Ts, 0.5)
                }),
            ),
            (
                "Ru above",
                one(D {
                    kind: Above,
                    ..d(Ru, 0.9)
                }),
            ),
            (
                "Ru below",
                one(D {
                    kind: Below,
                    ..d(Ru, 0.0)
                }),
            ),
            (
                "Ru range",
                one(D {
                    kind: Range,
                    band: 0.02,
                    ..d(Ru, 0.6)
                }),
            ),
            (
                "Tu centre band",
                one(D {
                    kind: CenterBand,
                    band: 0.5,
                    ..d(Tu, 0.6)
                }),
            ),
            (
                "Ru integral",
                one(D {
                    integral: true,
                    ..d(Ru, 0.3)
                }),
            ),
            (
                "Rp integral",
                one(D {
                    integral: true,
                    ..d(Rp, 0.3)
                }),
            ),
        ];
        for (name, curve) in [("colour Rs", Rs), ("colour Ru", Ru), ("colour Tp", Tp)] {
            let mut s = MeritSpec::new();
            add_color(&mut s, angle, curve);
            out.push((name, s));
        }
        // A polarizing beam splitter: the branches pull in opposite ways.
        let mut pbs = MeritSpec::new();
        add(&mut pbs, angle, d(Rs, 1.0));
        add(&mut pbs, angle, d(Rp, 0.0));
        out.push(("PBS", pbs));
        // Two phase demands on one element and branch, at least one of
        // them past the wrap wherever arg(r_p) is outside [-0.64, 0.64]:
        // they share one bucket, which must fold their wrapped residuals.
        let mut wrap_pair = MeritSpec::new();
        add(
            &mut wrap_pair,
            angle,
            D {
                phase: true,
                ..d(Rp, 2.5)
            },
        );
        add(
            &mut wrap_pair,
            angle,
            D {
                phase: true,
                ..d(Rp, -2.5)
            },
        );
        out.push(("Rp phase pair across the wrap", wrap_pair));
        // Integral phase (plan §5.9): arg(t_p) sits ~1.5 rad from 0 over the
        // band, past pi/N with N = 5, where the pair's N*G wrapped; arg(r_p)
        // crosses ±pi inside the band at 0, 30 and 60 deg, where the gap of
        // the raw-phase mean is 2*pi off the wrapped residuals' mean.
        out.push((
            "Tp phase integral, mean gap past pi/N",
            one(D {
                phase: true,
                integral: true,
                ..d(Tp, 0.0)
            }),
        ));
        out.push((
            "Rp phase integral, band across ±pi",
            one(D {
                phase: true,
                integral: true,
                ..d(Rp, 3.1)
            }),
        ));
        out.push((
            "Rp phase integral centre band, across ±pi",
            one(D {
                phase: true,
                integral: true,
                kind: CenterBand,
                band: 0.01,
                ..d(Rp, 3.1)
            }),
        ));
        let mut int_pair = MeritSpec::new();
        for target in [3.1, -1.0] {
            add(
                &mut int_pair,
                angle,
                D {
                    phase: true,
                    integral: true,
                    ..d(Rp, target)
                },
            );
        }
        out.push(("two integral Rp phases on one bucket", int_pair));
        // W4: non-uniform point weights, pointwise and integral, on an
        // intensity pair, an unpolarized demand, a kind and a phase.
        for (name, x) in [
            ("Rs weighted", d(Rs, 0.3)),
            ("Tu weighted", d(Tu, 0.6)),
            (
                "Ru range weighted",
                D {
                    kind: Range,
                    band: 0.02,
                    ..d(Ru, 0.6)
                },
            ),
            (
                "Tp phase weighted",
                D {
                    phase: true,
                    ..d(Tp, 0.5)
                },
            ),
            (
                "Ru integral weighted",
                D {
                    integral: true,
                    ..d(Ru, 0.3)
                },
            ),
            (
                "Rs above integral weighted",
                D {
                    kind: Above,
                    integral: true,
                    ..d(Rs, 0.9)
                },
            ),
        ] {
            out.push((
                name,
                one(D {
                    weighted: true,
                    ..x
                }),
            ));
        }
        // W5: Log demands fold their tangent at the operating point, on a
        // pair, an unpolarized demand, absorption, the back side, the kinds,
        // and the integral (plain and weighted). Levels sit near the sims'
        // values so the one-sided kinds are active somewhere and idle
        // elsewhere.
        for (name, x) in [
            ("Rs log", lg(Rs, 0.3)),
            ("Tu log", lg(Tu, 0.6)),
            ("Au log", lg(Au, 0.05)),
            ("RBp log", lg(RBp, 0.3)),
            (
                "Ru log below",
                D {
                    kind: Below,
                    ..lg(Ru, 0.01)
                },
            ),
            (
                "Tu log range",
                D {
                    kind: Range,
                    band: 0.02,
                    ..lg(Tu, 0.6)
                },
            ),
            (
                "Ru log centre band",
                D {
                    kind: CenterBand,
                    band: 0.1,
                    ..lg(Ru, 0.2)
                },
            ),
            (
                "Rs log weighted",
                D {
                    weighted: true,
                    ..lg(Rs, 0.3)
                },
            ),
            (
                "Ru log integral",
                D {
                    integral: true,
                    ..lg(Ru, 0.3)
                },
            ),
            (
                "Tp log integral weighted",
                D {
                    integral: true,
                    weighted: true,
                    ..lg(Tp, 0.9)
                },
            ),
            (
                "Rs log above integral",
                D {
                    kind: Above,
                    integral: true,
                    ..lg(Rs, 0.9)
                },
            ),
        ] {
            out.push((name, one(x)));
        }
        // Retardance (W8a) is identically zero at normal incidence, so the
        // merit cannot move there; `retardance_is_flat_at_normal_incidence`
        // checks that instead.
        if angle.get() != 0.0 {
            out.extend(retardance_cases(angle));
        }
        out
    }

    /// W8a: retardance demands on r_f, t_f and t_b, pointwise and
    /// integral, two kinds, and one target across the wrap (DeltaR is RetR
    /// with the target shifted by pi, so 3.0 stands for DeltaR = -0.14).
    fn retardance_cases(angle: Radians) -> Vec<(&'static str, MeritSpec)> {
        use ConstraintKind::*;
        use CurveId::*;
        let one = |x: D| {
            let mut s = MeritSpec::new();
            add(&mut s, angle, x);
            s
        };
        let mut both = MeritSpec::new();
        add(&mut both, angle, ret(Rp, 0.4));
        add(&mut both, angle, d(Rs, 0.3));
        add(
            &mut both,
            angle,
            D {
                phase: true,
                ..d(Rp, -0.5)
            },
        );
        vec![
            ("RetR", one(ret(Rp, 0.4))),
            ("RetT", one(ret(Tp, -0.3))),
            ("DeltaR (RetR across the wrap)", one(ret(Rp, 3.0))),
            ("RetTB", one(ret(TBp, 0.2))),
            (
                "RetR integral",
                one(D {
                    integral: true,
                    ..ret(Rp, 0.4)
                }),
            ),
            (
                "RetT above",
                one(D {
                    kind: Above,
                    ..ret(Tp, 0.5)
                }),
            ),
            (
                "RetR range",
                one(D {
                    kind: Range,
                    band: 0.01,
                    ..ret(Rp, 0.6)
                }),
            ),
            ("RetR + Rs + Rp phase", both),
            (
                "RetR weighted",
                one(D {
                    weighted: true,
                    ..ret(Rp, 0.4)
                }),
            ),
            (
                "RetT integral weighted",
                one(D {
                    integral: true,
                    weighted: true,
                    ..ret(Tp, -0.3)
                }),
            ),
            // Plan §5.9: a mean retardance gap past pi/N (N = 5).
            (
                "RetR integral, mean gap past pi/N",
                one(D {
                    integral: true,
                    ..ret(Rp, -1.5)
                }),
            ),
        ]
    }

    #[test]
    fn retardance_p_is_the_merit_slope_at_45_and_70_deg() {
        for angle in [deg(45.0), deg(70.0)] {
            for (name, spec) in retardance_cases(angle) {
                check(name, spec, angle);
            }
        }
    }

    /// At normal incidence s and p coincide, Ret_R is identically 0 and no
    /// seed moves it: the fold's P is zero to round-off, like the slope.
    #[test]
    fn retardance_is_flat_at_normal_incidence() {
        for (name, x) in [
            ("RetR", ret(CurveId::Rp, 0.4)),
            ("RetT", ret(CurveId::Tp, 0.4)),
        ] {
            let mut spec = MeritSpec::new();
            add(&mut spec, Radians::ZERO, x);
            for (i, (p, s)) in p_and_slope(spec, Radians::ZERO).into_iter().enumerate() {
                assert!(
                    p.abs() < 1e-9 && s.abs() < 1e-6,
                    "{name}, site {i}: P {p}, slope {s}"
                );
            }
        }
    }

    #[test]
    fn p_is_the_merit_slope_at_normal_incidence() {
        for (name, spec) in cases(Radians::ZERO) {
            check(name, spec, Radians::ZERO);
        }
    }

    #[test]
    fn p_is_the_merit_slope_at_30_deg() {
        for (name, spec) in cases(deg(30.0)) {
            check(name, spec, deg(30.0));
        }
    }

    #[test]
    fn p_is_the_merit_slope_at_60_deg() {
        for (name, spec) in cases(deg(60.0)) {
            check(name, spec, deg(60.0));
        }
    }

    // --- W6: a coated plate, both faces, in air ---------------------------

    /// Front films, the incoherent plate, back films, in air.
    fn plate_films(nw: usize) -> Vec<LayerSpec> {
        let mut plate = LayerSpec::constant("G", 1.52, 2e-6, 200_000.0, nw);
        plate.coherent = false;
        plate.optimize = false;
        plate.needle = false;
        vec![
            LayerSpec::constant("H", 2.30, 0.0, 118.0, nw),
            LayerSpec::constant("L", 1.46, 0.0, 203.0, nw),
            plate,
            LayerSpec::constant("Habs", 2.30, 0.004, 64.0, nw),
            LayerSpec::constant("L", 1.46, 0.0, 100.0, nw),
        ]
    }

    fn in_air(films: Vec<LayerSpec>) -> DesignStack {
        let nw = WLS.len();
        let half = || {
            let mut a = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
            a.optimize = false;
            a.needle = false;
            a
        };
        DesignStack::with_films(half(), half(), films).unwrap()
    }

    fn plate() -> DesignStack {
        in_air(plate_films(WLS.len()))
    }

    #[test]
    fn p_is_the_merit_slope_through_an_incoherent_plate() {
        for angle in [Radians::ZERO, deg(30.0)] {
            for (name, spec) in cases(angle) {
                if spec.uses_phase() {
                    continue;
                }
                check_on(name, plate(), spec, angle);
            }
        }
    }

    /// Front rows against the engine's own block walk; back rows against
    /// the front rows of the mirrored stack (air on both sides, so the
    /// flux factors and the angle carry over).
    #[test]
    fn incoherent_forward_matches_the_engine() {
        use crate::smatrix::core_engine::{MODE_B, solve_point};
        use CurveId::*;
        let angle = deg(30.0);
        let mut spec = MeritSpec::new();
        for id in [RBs, RBp, TBs, TBp] {
            add(&mut spec, angle, d(id, 0.3));
        }
        let c = ctx(spec, angle);
        let st = plate();
        let sim = c.simulate(&st).unwrap();
        let sa = st.solver_arrays();
        let nl = sa.n_layers as usize;
        for (w, &lam) in WLS.iter().enumerate() {
            let ns: Vec<Complex64> = (0..nl)
                .map(|l| {
                    let b = w * nl * 2 + l * 2;
                    Complex64::new(sa.n_stack_cache[b], sa.n_stack_cache[b + 1])
                })
                .collect();
            let inv: Vec<Complex64> = ns.iter().map(|&n| 1.0 / n).collect();
            let e = solve_point(
                nl - 1,
                lam,
                c.sin_theta[0],
                &ns,
                &inv,
                &sa.thicknesses,
                &sa.incoherent_flags,
                &sa.rough_types,
                &sa.rough_vals,
                MODE_B,
                true,
                true,
                false,
            );
            for (id, v) in [(Rs, e.rs), (Rp, e.rp), (Ts, e.ts), (Tp, e.tp)] {
                assert_eq!(sim.curve(id).unwrap()[w], v, "{id:?} at {lam}");
            }
        }
        let mut films = plate_films(WLS.len());
        films.reverse();
        let mirror = c.simulate(&in_air(films)).unwrap();
        for (back, front) in [(RBs, Rs), (RBp, Rp), (TBs, Ts), (TBp, Tp)] {
            for w in 0..WLS.len() {
                let (b, f) = (
                    sim.back_curve(back).unwrap()[w],
                    mirror.curve(front).unwrap()[w],
                );
                assert!(
                    (b - f).abs() < 1e-12,
                    "{back:?} {b} vs mirrored {front:?} {f}"
                );
            }
        }
        // The plate really is incoherent: a coherent one fringes.
        let mut films = plate_films(WLS.len());
        films[2].coherent = true;
        let coh = c.simulate(&in_air(films)).unwrap();
        let gap = (0..WLS.len())
            .map(|w| (coh.curve(Rs).unwrap()[w] - sim.curve(Rs).unwrap()[w]).abs())
            .fold(0.0, f64::max);
        assert!(gap > 1e-3, "coherent and incoherent plates agree ({gap})");
    }

    #[test]
    fn incoherent_stack_refuses_phase_and_deposits() {
        let angle = deg(30.0);
        let mut spec = MeritSpec::new();
        add(
            &mut spec,
            angle,
            D {
                phase: true,
                ..d(CurveId::Rs, 0.5)
            },
        );
        let err = ctx(spec, angle).simulate(&plate()).unwrap_err();
        assert!(err.contains("phase"), "{err}");
        let mut spec = MeritSpec::new();
        add(&mut spec, angle, d(CurveId::Rs, 0.3));
        let err = ctx(spec, angle)
            .simulate_with_deposits(&plate(), &[0])
            .unwrap_err();
        assert!(err.contains("coherent"), "{err}");
    }

    /// A boundary site on top of the plate is a needle at the bottom of
    /// the film above it.
    #[test]
    fn a_site_on_the_plate_is_the_bottom_of_the_film_above() {
        use crate::smatrix::synthesis::needle_pass::ScanSite;
        let angle = deg(30.0);
        let mut spec = MeritSpec::new();
        add(&mut spec, angle, d(CurveId::Ru, 0.3));
        add(&mut spec, angle, d(CurveId::TBs, 0.6));
        let c = ctx(spec, angle);
        let st = plate();
        let sim = c.simulate(&st).unwrap();
        let fold = build_needle_targets(&c.spec, &[angle], &WLS, Some(&sim)).unwrap();
        let sa = st.solver_arrays();
        let nk = vec![Complex64::new(NEEDLE_N, 0.0); WLS.len()];
        let z = 118.0 + 203.0;
        let sites = [
            ScanSite {
                film_idx: 2,
                depth_into_layer_nm: 0.0,
                z_nm: z,
                boundary: true,
            },
            ScanSite {
                film_idx: 1,
                depth_into_layer_nm: 203.0,
                z_nm: z,
                boundary: false,
            },
        ];
        let res = needle_pass_scan(
            &NeedlePassInput {
                n_stack_cache: &sa.n_stack_cache,
                thicknesses: &sa.thicknesses,
                rough_types: &sa.rough_types,
                rough_vals: &sa.rough_vals,
                n_layers: sa.n_layers as usize,
                wavls: &WLS,
                sin_theta: &c.sin_theta,
                fold: &fold,
                needle_n_per_wav: &nk,
                start_idx: 0,
                end_idx: (sa.n_layers - 1) as usize,
                incoherent_flags: &sa.incoherent_flags,
                calc_s: true,
                calc_p: true,
            },
            &sites,
        )
        .unwrap();
        assert!(res.p_profile[0] != 0.0);
        assert_eq!(res.p_profile[0], res.p_profile[1]);
    }

    /// LM through the plate: it differences, lowers the merit, and leaves
    /// the plate alone.
    #[test]
    fn lm_optimizes_the_coatings_around_a_fixed_plate() {
        let mut spec = MeritSpec::new();
        add(&mut spec, Radians::ZERO, d(CurveId::Ru, 0.02));
        let mut c = ctx(spec, Radians::ZERO);
        let mut st = plate();
        let m0 = c.evaluate_merit(&st).unwrap();
        let m1 = c.optimize_thicknesses(&mut st).unwrap();
        assert!(m1 < 0.5 * m0, "merit {m0} -> {m1}");
        let f = &st.films()[2];
        assert!(!f.coherent && f.d_nm == 200_000.0);
    }
}
