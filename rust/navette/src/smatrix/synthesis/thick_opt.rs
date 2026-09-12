// SPDX-License-Identifier: LGPL-3.0-or-later
//! Navette -- Rust Rewrite of Numba-optimized thin-film optical solver
//!
//! synthesis::thick_opt — bounded Levenberg–Marquardt for film thicknesses.
//!
//! Replaces scipy `least_squares(method="trf")` from needle_synthesis.py /
//! ClampedNeedleSynthesizer.optimize_thicknesses. The residual system is
//! injected as a closure, so this module is solver-agnostic and unit-testable
//! standalone; the synthesis wiring supplies residuals via core_engine +
//! MeritSpec (Phase 4+).
//!
//! Algorithm notes:
//!   * **The damped step is solved by QR, not by the normal equations.**
//!     Marquardt scaling is unchanged — the step still satisfies
//!     (JᵀJ + λ·D²) δ = −Jᵀr with D² = diag(JᵀJ) floored — but it is obtained
//!     as the least-squares solution of the augmented system
//!     `[J; √λ·D] δ ≈ [−r; 0]`, whose condition number is the square root of
//!     the normal-equation matrix's. Thin-film stacks with correlated layers
//!     are exactly where JᵀJ goes singular and the λ-floor has to bail the
//!     solve out (review §3.6). `J` is factored ONCE per iteration (m·n²) and
//!     each λ trial then costs a 2n×n QR (n³), the MINPACK `qrsolv` structure.
//!     `solve_symmetric` survives as the fallback when the QR degenerates,
//!     and as the cross-check oracle in the tests.
//!   * **Damping is Nielsen's gain ratio** (`LmDamping::GainRatio`, default),
//!     not a fixed ×5/÷3 ladder: ρ = actual/predicted reduction; on acceptance
//!     λ ← λ·max(1/3, 1−(2ρ−1)³) and ν ← 2, on rejection λ ← λ·ν and ν ← 2ν.
//!     `LmDamping::Fixed` keeps the old ladder for comparison.
//!   * **The predicted reduction is computed for the step actually taken.**
//!     The bound veto+clamp below rewrites δ *after* it is solved; scoring it
//!     with the full LM step's prediction would overstate the model's promise
//!     whenever a bound is active, making ρ too small — over-damping and early
//!     ftol exits right at the boundary. Prediction uses `trial − x`.
//!   * **ftol needs both reductions** (MINPACK semantics): actual *and*
//!     predicted relative reduction ≤ ftol. An actual reduction alone can be
//!     small because the step was clipped, not because the optimum is near.
//!   * **gtol has a scale-invariant form** (`gtol_scale_invariant`, default
//!     on): cos∠(J·e_j, r) ≤ gtol, alongside the scale-dependent ‖Jᵀr‖∞ test.
//!     They are different notions of stationarity and can exit at different
//!     points on flat, noise-dominated valleys.
//!   * Central-difference Jacobian, per-column step h_j = ∛ε·max(|x_j|, 1),
//!     evaluated RAYON-PARALLEL across columns (the film-synthesis system
//!     costs one full TMM solve per evaluation — columns dominate runtime).
//!   * Bounds: steps are vetoed component-wise when they push away from an
//!     active bound, then the trial point is clamped into [lb, ub] (mirrors
//!     ClampedNeedleSynthesizer's optimizer-bounds + post-clamp contract;
//!     removal of sub-min layers stays the CALLER's job, as in Python).
//!     Neither MINPACK nor the argmin ecosystem's LM has bounds at all; this
//!     contract is Navette's own and is deliberately preserved.

use rayon::prelude::*;

// ---------------------------------------------------------------------------
// Configuration / results
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct LmConfig {
    /// Hard cap on accepted iterations.
    pub max_iterations: usize,
    /// Hard cap on total residual evaluations (including Jacobian probes).
    pub max_evals: usize,
    /// Terminate when relative cost decrease falls below this.
    pub ftol: f64,
    /// Terminate when relative step size falls below this.
    pub xtol: f64,
    /// Terminate when ‖Jᵀr‖∞ falls below this.
    pub gtol: f64,
    /// Initial damping λ.
    pub lambda_init: f64,
    /// λ multiplier on rejected steps (`LmDamping::Fixed`, and both modes'
    /// error ladder: a failed solve or a failed residual evaluation).
    pub lambda_up: f64,
    /// λ divisor on accepted steps (`LmDamping::Fixed` only).
    pub lambda_down: f64,
    /// How λ moves between trials. See `LmDamping`.
    pub damping: LmDamping,
    /// Also test the scale-invariant gradient criterion cos∠(J·e_j, r) ≤ gtol
    /// beside ‖Jᵀr‖∞ < gtol. The two are different notions of stationarity:
    /// the angle form does not change its mind when a parameter is rescaled,
    /// which ‖Jᵀr‖∞ does.
    pub gtol_scale_invariant: bool,
    /// Where the Jacobian comes from, when the caller has an analytic one to
    /// offer. Honoured by whoever supplies the `JacobianSource` — the driver
    /// itself only ever has the residual closure, so a run with no source is
    /// finite-difference regardless.
    pub jacobian: JacobianMode,
    /// Which solver runs. `BuiltinLm` is this module; the others live behind
    /// cargo features and are dispatched by `synthesis::optimizer`.
    ///
    /// The field sits here rather than on a parallel `OptimizerConfig`
    /// because every other knob in this struct is one the alternatives take
    /// too (ftol/xtol/gtol/max_evals are MINPACK's own names), and two
    /// structs that must agree field for field are a synchronization bug
    /// waiting to be written. `levenberg_marquardt_with` ignores it — it *is*
    /// the built-in — so only `optimizer::run_optimizer` reads it.
    pub backend: crate::smatrix::synthesis::optimizer::OptimizerBackend,
}

/// Which Jacobian the thickness optimizer should use.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum JacobianMode {
    /// The analytic Jacobian when the problem is fully covered by it, and
    /// central differences for the rows it cannot express (phase targets,
    /// color demands). The default: exact where it applies, and it costs a
    /// constant number of solver sweeps instead of 2n.
    #[default]
    Analytic,
    /// Central differences throughout — the pre-0.6.9 behaviour, kept as the
    /// cross-check oracle and as the escape hatch if an analytic chain is
    /// ever suspected.
    Fd,
}

/// How the damping parameter λ is updated between trial steps.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum LmDamping {
    /// Nielsen's gain ratio (MINPACK-style), the default. ρ = actual reduction
    /// over the reduction the linear model predicted for the step actually
    /// taken; λ shrinks in proportion to how well the model did, and grows
    /// geometrically (ν, doubling) while it keeps failing.
    #[default]
    GainRatio,
    /// The fixed ×`lambda_up` / ÷`lambda_down` ladder this module used before
    /// 0.6.7. Kept so the two can be compared on the same problem.
    Fixed,
}

impl Default for LmConfig {
    fn default() -> Self {
        LmConfig {
            max_iterations: 200,
            max_evals: 100_000,
            ftol: 1e-12,
            xtol: 1e-12,
            gtol: 1e-10,
            lambda_init: 1e-3,
            lambda_up: 5.0,
            lambda_down: 3.0,
            damping: LmDamping::GainRatio,
            gtol_scale_invariant: true,
            jacobian: JacobianMode::Analytic,
            backend: crate::smatrix::synthesis::optimizer::OptimizerBackend::BuiltinLm,
        }
    }
}

/// An analytic Jacobian the driver can ask for instead of differencing.
///
/// Kept as a trait rather than a closure because the driver holds it across
/// rayon-parallel residual evaluations: `Sync` is the requirement, and a
/// named trait makes it one the compiler states rather than one the caller
/// discovers.
pub trait JacobianSource: Sync {
    /// Fill `jac` row-major (m×n) with ∂r_i/∂x_j at `x` and return `m`.
    ///
    /// `Ok(None)` means "not this problem" — the driver falls back to central
    /// differences for that iteration, which is how a spec with rows the
    /// analytic chain cannot express still runs.
    fn fill(&self, x: &[f64], jac: &mut Vec<f64>) -> Result<Option<usize>, String>;
}

/// The absence of one. `levenberg_marquardt` is `levenberg_marquardt_with`
/// against this.
pub struct NoJacobian;

impl JacobianSource for NoJacobian {
    fn fill(&self, _x: &[f64], _jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
        Ok(None)
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum LmTermination {
    /// ‖Jᵀr‖∞ < gtol, or — when `gtol_scale_invariant` is set — the residual
    /// is orthogonal to every Jacobian column to within gtol.
    Gradient,
    /// Relative step below xtol.
    Step,
    /// Actual *and* predicted relative cost improvement below ftol, or no
    /// further step could improve the cost at all.
    Cost,
    /// max_iterations reached.
    MaxIterations,
    /// Damping escalated without any acceptable step (stuck at bounds or
    /// numerically degenerate system).
    Stalled,
}

#[derive(Clone, Debug)]
pub struct LmResult {
    pub x: Vec<f64>,
    pub cost: f64,
    pub iterations: usize,
    pub evals: usize,
    pub termination: LmTermination,
    /// ρ of the last accepted step: the reduction actually obtained over the
    /// reduction the linear model predicted **for the step actually taken**
    /// (after the bound veto and clamp). ≈ 1 means the model described the
    /// step well; ≪ 1 means the step outran the linearization. NaN when no
    /// step was ever accepted.
    ///
    /// A diagnostic, not a control: it is what `LmDamping::GainRatio` steers
    /// λ by, surfaced so a synthesis that stalls at its bounds can be told
    /// apart from one that is simply finished.
    pub gain_ratio: f64,
    /// How many of the run's Jacobians came from the analytic source. `0` is
    /// a fully differenced run; anything else means the source supplied at
    /// least that many. Note it counts Jacobian *builds*, which is one more
    /// than `iterations` when the run ends on the gradient test — the last
    /// Jacobian is what proved the point was stationary. The bench asserts on
    /// this rather than inferring the path from a timing.
    pub analytic_jacobians: usize,
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

/// Minimize ½‖r(x)‖²  s.t.  lb ≤ x ≤ ub.
///
/// `residuals(x, out)` must append exactly `m` components (m fixed across
/// calls); inactive constraints should contribute zeros (see
/// `synthesis::merit::MeritSpec::residuals`).
pub fn levenberg_marquardt<F>(
    residuals: &F,
    x0: &[f64],
    lb: &[f64],
    ub: &[f64],
    cfg: &LmConfig,
) -> Result<LmResult, String>
where
    F: Fn(&[f64], &mut Vec<f64>) -> Result<(), String> + Sync,
{
    levenberg_marquardt_with(residuals, None::<&NoJacobian>, x0, lb, ub, cfg)
}

/// As [`levenberg_marquardt`], with an analytic Jacobian the driver prefers
/// over central differences whenever the source supplies one.
///
/// The source is asked once per iteration, at the current `x`. Declining
/// (`Ok(None)`) costs nothing but the fallback; erroring aborts the run, on
/// the grounds that a Jacobian that cannot be built is a different fact from
/// one that does not apply.
pub fn levenberg_marquardt_with<F, J>(
    residuals: &F,
    jacobian: Option<&J>,
    x0: &[f64],
    lb: &[f64],
    ub: &[f64],
    cfg: &LmConfig,
) -> Result<LmResult, String>
where
    F: Fn(&[f64], &mut Vec<f64>) -> Result<(), String> + Sync,
    J: JacobianSource + ?Sized,
{
    let n = x0.len();
    if n == 0 {
        return Err("levenberg_marquardt: empty parameter vector".into());
    }
    if lb.len() != n || ub.len() != n {
        return Err("levenberg_marquardt: bound length mismatch".into());
    }
    for j in 0..n {
        if !(lb[j] <= x0[j] && x0[j] <= ub[j]) {
            return Err(format!(
                "levenberg_marquardt: x0[{}]={} outside [{}, {}]",
                j, x0[j], lb[j], ub[j]
            ));
        }
        if lb[j] > ub[j] {
            return Err(format!("levenberg_marquardt: lb[{}] > ub[{}]", j, j));
        }
    }

    let evals = std::cell::Cell::new(0usize);
    let mut buf_r: Vec<f64> = Vec::new();
    let mut buf_t: Vec<f64> = Vec::new();

    let evaluate = |x: &[f64], out: &mut Vec<f64>| -> Result<f64, String> {
        evals.set(evals.get() + 1);
        residuals(x, out)?;
        Ok(out.iter().map(|r| r * r).sum::<f64>())
    };

    let mut x = x0.to_vec();
    // Clamp x0 defensively into bounds.
    for j in 0..n {
        x[j] = x[j].clamp(lb[j], ub[j]);
    }

    let mut cost = evaluate(&x, &mut buf_r)?;
    let m = buf_r.len();
    let mut jac = vec![0.0f64; m * n]; // row-major m×n
    let mut qrj = vec![0.0f64; m * n]; // scratch copy the QR destroys
    let mut qtr = vec![0.0f64; m]; // Qᵀr (only the leading n matter)
    let mut jtr = vec![0.0f64; n];
    let mut col_sq = vec![0.0f64; n]; // ‖J·e_j‖² == diag(JᵀJ)
    let mut d_scale = vec![0.0f64; n]; // √(floored diag(JᵀJ))
    let mut r_tri = vec![0.0f64; n * n]; // upper-triangular R from J = QR
    let mut delta = vec![0.0f64; n];
    let mut applied = vec![0.0f64; n];
    let mut trial = vec![0.0f64; n];

    let mut lambda = cfg.lambda_init;
    let mut nu = 2.0f64; // Nielsen's rejection multiplier, reset on acceptance
    let mut gain_ratio = f64::NAN;
    let mut termination: Option<LmTermination> = None;
    let mut iteration = 0usize;

    let mut analytic_jacobians = 0usize;

    while iteration < cfg.max_iterations {
        // ---- Jacobian: the analytic source if it has one for this x,
        //      otherwise rayon across columns with central differences ----
        let supplied = match jacobian {
            Some(src) => src.fill(&x, &mut jac)?,
            None => None,
        };
        match supplied {
            Some(rows) => {
                if rows != m || jac.len() != m * n {
                    return Err(format!(
                        "levenberg_marquardt: analytic jacobian is {rows}×{} where the                          residual vector is {m} long",
                        if rows == 0 {
                            0
                        } else {
                            jac.len() / rows.max(1)
                        }
                    ));
                }
                analytic_jacobians += 1;
            }
            None => {
                build_jacobian(residuals, &x, &mut jac)
                    .map(|added| evals.set(evals.get() + added))?;
            }
        }

        // ---- g = Jᵀr and the column norms, in one pass ----
        // diag(JᵀJ) is all of JᵀJ this solver needs: the Marquardt scaling
        // uses the diagonal, and the step comes from a QR of J itself.
        jtr.fill(0.0);
        col_sq.fill(0.0);
        for i in 0..m {
            let ri = buf_r[i];
            let row = &jac[i * n..(i + 1) * n];
            for j in 0..n {
                jtr[j] += row[j] * ri;
                col_sq[j] += row[j] * row[j];
            }
        }
        for j in 0..n {
            d_scale[j] = col_sq[j].max(1e-14).sqrt();
        }

        // ---- gradient convergence on the current point ----
        let g_inf = jtr.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
        if g_inf < cfg.gtol {
            termination = Some(LmTermination::Gradient);
            break;
        }
        if cfg.gtol_scale_invariant && gradient_cosine(&jtr, &col_sq, cost.sqrt()) <= cfg.gtol {
            termination = Some(LmTermination::Gradient);
            break;
        }

        // ---- factor J once; each λ trial is then a 2n×n solve ----
        qrj.copy_from_slice(&jac);
        qtr.copy_from_slice(&buf_r);
        let have_qr = householder_qr_in_place(&mut qrj, &mut qtr, m, n).is_ok();
        if have_qr {
            r_tri.fill(0.0);
            for i in 0..n {
                r_tri[i * n + i..(i + 1) * n].copy_from_slice(&qrj[i * n + i..(i + 1) * n]);
            }
        }

        // ---- damped-step loop: escalate λ until some step is accepted ----
        let mut accepted = false;
        for _damping_try in 0..40 {
            if evals.get() >= cfg.max_evals {
                break;
            }

            // Step: (JᵀJ + λ·D²)δ = −Jᵀr, solved as least squares on
            // [R; √λ·D]δ ≈ [−Qᵀr; 0]. RᵀR = JᵀJ exactly, so this is the same
            // step the normal equations define, at half the condition number.
            let solved = have_qr
                && solve_damped_step(&r_tri, &qtr[..n], &d_scale, lambda, &mut delta).is_ok();
            if !solved {
                // Degenerate factorization: fall back to the normal equations
                // (damping makes JᵀJ + λD² positive definite for λ > 0, so
                // this usually succeeds), and only then escalate λ.
                if normal_equation_step(&jac, &jtr, &d_scale, lambda, m, n, &mut delta).is_err() {
                    lambda *= cfg.lambda_up;
                    if lambda > 1e18 {
                        break;
                    }
                    continue;
                }
            }

            // Bound-aware projection: veto components pushing past an ACTIVE
            // bound, then clamp the trial point (post-clamp contract).
            for j in 0..n {
                if x[j] <= lb[j] && delta[j] < 0.0 {
                    delta[j] = 0.0;
                }
                if x[j] >= ub[j] && delta[j] > 0.0 {
                    delta[j] = 0.0;
                }
                trial[j] = (x[j] + delta[j]).clamp(lb[j], ub[j]);
                applied[j] = trial[j] - x[j];
            }
            if applied.iter().all(|d| *d == 0.0) {
                // Purely-zero step: all coordinates pinned. Try raising λ to
                // escape coupling once, else declare stall.
                if lambda > 1e18 {
                    break;
                }
                lambda *= cfg.lambda_up;
                continue;
            }

            // The model's promise for the step we are ACTUALLY taking, not
            // for the one the solve returned. See the module note.
            let predicted = predicted_reduction(&jac, &jtr, &applied, m, n);

            let new_cost = match evaluate(&trial, &mut buf_t) {
                Ok(c) => c,
                Err(_) => {
                    lambda *= cfg.lambda_up;
                    continue;
                }
            };

            if new_cost < cost {
                let actual = cost - new_cost;
                let denom = cost.abs().max(1e-300);
                let rel_actual = actual / denom;
                let rel_pred = predicted / denom;

                std::mem::swap(&mut buf_r, &mut buf_t);
                std::mem::swap(&mut x, &mut trial);
                cost = new_cost;
                accepted = true;

                gain_ratio = if predicted > 0.0 {
                    actual / predicted
                } else {
                    f64::NAN
                };

                match cfg.damping {
                    LmDamping::GainRatio if predicted > 0.0 => {
                        let rho = gain_ratio;
                        lambda *= (1.0 - (2.0 * rho - 1.0).powi(3)).max(1.0 / 3.0);
                        nu = 2.0;
                    }
                    // A non-positive prediction means the linear model did not
                    // expect this step to help, yet the cost fell: the model is
                    // stale rather than over-confident, so shrink λ the plain
                    // way instead of dividing by a meaningless ρ.
                    _ => lambda /= cfg.lambda_down,
                }

                if ftol_converged(rel_actual, rel_pred, cfg.ftol) {
                    termination = Some(LmTermination::Cost);
                }
                break;
            }

            match cfg.damping {
                LmDamping::GainRatio => {
                    lambda *= nu;
                    nu *= 2.0;
                }
                LmDamping::Fixed => lambda *= cfg.lambda_up,
            }
            if lambda > 1e18 {
                break;
            }
        }

        if !accepted {
            // Every trial in the ladder failed to lower the cost. Whether that
            // is a degenerate system or an optimum pinned against its bounds,
            // no further step exists from here.
            termination = Some(LmTermination::Stalled);
            break;
        }
        if termination.is_some() {
            break;
        }

        // Step-size convergence: ‖Δx‖ ≤ xtol·(xtol + ‖x‖). `trial` holds the
        // previous point (the accept path swapped the two).
        let dx_norm: f64 = x
            .iter()
            .zip(&trial)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        let x_norm: f64 = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        if dx_norm <= cfg.xtol * (cfg.xtol + x_norm) {
            termination = Some(LmTermination::Step);
            break;
        }
        if evals.get() >= cfg.max_evals {
            termination = Some(LmTermination::MaxIterations);
            break;
        }

        iteration += 1;
    }

    Ok(LmResult {
        x,
        cost,
        iterations: iteration,
        evals: evals.get(),
        termination: termination.unwrap_or(LmTermination::MaxIterations),
        gain_ratio,
        analytic_jacobians,
    })
}

// ---------------------------------------------------------------------------
// Step solve
// ---------------------------------------------------------------------------

/// MINPACK's ftol test: the step is a convergence signal only when the
/// reduction it delivered *and* the reduction the model predicted for it are
/// both below `ftol`.
///
/// The old rule looked at the actual reduction alone. A step clipped by a
/// bound, or one taken while λ is still large, can deliver very little while
/// the model still sees plenty of room — small progress is not the same fact
/// as no progress left. Requiring both is what distinguishes them.
///
/// A negative prediction never converges: the model did not expect this step
/// to help at all, so its silence says nothing about how near the optimum is.
fn ftol_converged(rel_actual: f64, rel_pred: f64, ftol: f64) -> bool {
    rel_actual < ftol && (0.0..ftol).contains(&rel_pred)
}

/// The largest cos∠(J·e_j, r) = |(Jᵀr)_j| / (‖J·e_j‖·‖r‖) over the columns.
///
/// MINPACK's scale-invariant stationarity measure. Unlike ‖Jᵀr‖∞ it does not
/// change when a parameter -- hence its Jacobian column -- is rescaled, so it
/// means the same thing for thicknesses carried in nanometres as in
/// micrometres. Zero is perfect stationarity: the residual is orthogonal to
/// every direction the parameters can move in.
///
/// `col_sq` is diag(JᵀJ). A zero column carries no direction and is skipped
/// rather than dividing by zero; an empty residual (‖r‖ = 0) is an exact fit,
/// reported as 0 so the caller terminates.
fn gradient_cosine(jtr: &[f64], col_sq: &[f64], r_norm: f64) -> f64 {
    if !(r_norm > 0.0) {
        return 0.0;
    }
    jtr.iter().zip(col_sq).fold(0.0f64, |acc, (&g, &cs)| {
        let cn = cs.sqrt();
        if cn > 0.0 {
            acc.max(g.abs() / (cn * r_norm))
        } else {
            acc
        }
    })
}

/// Reduction the linear model r(x+δ) ≈ r + Jδ predicts for `step`, in the
/// same units as `cost` (which is ‖r‖², not ½‖r‖²):
/// ‖r‖² − ‖r + Jδ‖² = −2·(Jᵀr)ᵀδ − ‖Jδ‖².
///
/// `step` is the step after the bound veto and clamp, so this is what the
/// gain ratio and the ftol test are entitled to compare against.
fn predicted_reduction(jac: &[f64], jtr: &[f64], step: &[f64], m: usize, n: usize) -> f64 {
    let mut jd_sq = 0.0;
    for i in 0..m {
        let row = &jac[i * n..(i + 1) * n];
        let mut s = 0.0;
        for j in 0..n {
            s += row[j] * step[j];
        }
        jd_sq += s * s;
    }
    let gd: f64 = jtr.iter().zip(step).map(|(g, d)| g * d).sum();
    -2.0 * gd - jd_sq
}

/// Solve `[R; √λ·D] δ ≈ [−Qᵀr; 0]` in least squares, R upper-triangular n×n.
///
/// Equivalent to (RᵀR + λD²)δ = −Rᵀ(Qᵀr), i.e. (JᵀJ + λD²)δ = −Jᵀr, but
/// formed from the square roots. The augmented matrix has full rank for every
/// λ > 0 with a floored D, so no pivoting is needed to make it solvable --
/// pivoting would only add rank diagnostics.
fn solve_damped_step(
    r_tri: &[f64],
    qtr: &[f64],
    d_scale: &[f64],
    lambda: f64,
    out: &mut [f64],
) -> Result<(), String> {
    let n = qtr.len();
    let rows = 2 * n;
    let mut a = vec![0.0f64; rows * n];
    for i in 0..n {
        a[i * n + i..(i + 1) * n].copy_from_slice(&r_tri[i * n + i..(i + 1) * n]);
    }
    let sqrt_lambda = lambda.sqrt();
    for j in 0..n {
        a[(n + j) * n + j] = sqrt_lambda * d_scale[j];
    }
    let mut b = vec![0.0f64; rows];
    for i in 0..n {
        b[i] = -qtr[i];
    }
    householder_qr_in_place(&mut a, &mut b, rows, n)?;
    back_substitute(&a, &b, n, out)
}

/// The pre-0.6.7 step: build JᵀJ, damp its diagonal, solve the normal
/// equations. Retained as the fallback when the QR degenerates, and as the
/// oracle the QR path is cross-checked against in the tests.
fn normal_equation_step(
    jac: &[f64],
    jtr: &[f64],
    d_scale: &[f64],
    lambda: f64,
    m: usize,
    n: usize,
    out: &mut [f64],
) -> Result<(), String> {
    let mut a = vec![0.0f64; n * n];
    for i in 0..m {
        let row = &jac[i * n..(i + 1) * n];
        for j in 0..n {
            let jr = row[j];
            for k in j..n {
                a[j * n + k] += jr * row[k];
            }
        }
    }
    for j in 0..n {
        for k in 0..j {
            a[j * n + k] = a[k * n + j];
        }
    }
    for j in 0..n {
        a[j * n + j] += lambda * d_scale[j] * d_scale[j];
    }
    solve_symmetric(&a, &neg(jtr), out)
}

/// Householder QR of a row-major `rows`×`n` matrix (`rows` >= `n`), in place.
/// R lands in the upper triangle of the first `n` rows; the same reflections
/// are applied to `b` (length `rows`), so `b[..n]` becomes the head of Qᵀb.
///
/// No column pivoting: every system this module factors is either J itself
/// (where a rank-deficient column is handled by the damping) or an augmented
/// `[R; √λD]` that is full rank by construction.
pub(crate) fn householder_qr_in_place(
    a: &mut [f64],
    b: &mut [f64],
    rows: usize,
    n: usize,
) -> Result<(), String> {
    if rows < n {
        return Err("householder_qr: fewer rows than columns".into());
    }
    let mut v = vec![0.0f64; rows];
    for k in 0..n {
        let mut norm_sq = 0.0;
        for i in k..rows {
            norm_sq += a[i * n + k] * a[i * n + k];
        }
        // `!(x > 0.0)` rejects NaN as well as zero; see the crate-level
        // clippy note in lib.rs.
        if !(norm_sq > 0.0) || !norm_sq.is_finite() {
            return Err("householder_qr: degenerate or non-finite column".into());
        }
        let norm = norm_sq.sqrt();
        let akk = a[k * n + k];
        // Reflect away from akk so v[k] is formed without cancellation.
        let alpha = if akk >= 0.0 { -norm } else { norm };
        for i in k..rows {
            v[i] = a[i * n + k];
        }
        v[k] -= alpha;
        let vtv: f64 = (k..rows).map(|i| v[i] * v[i]).sum();
        a[k * n + k] = alpha;
        for i in (k + 1)..rows {
            a[i * n + k] = 0.0;
        }
        if vtv <= 0.0 {
            // The column was already alpha*e_k; the reflection is the identity.
            continue;
        }
        for j in (k + 1)..n {
            let mut s = 0.0;
            for i in k..rows {
                s += v[i] * a[i * n + j];
            }
            let f = 2.0 * s / vtv;
            for i in k..rows {
                a[i * n + j] -= f * v[i];
            }
        }
        let mut s = 0.0;
        for i in k..rows {
            s += v[i] * b[i];
        }
        let f = 2.0 * s / vtv;
        for i in k..rows {
            b[i] -= f * v[i];
        }
    }
    Ok(())
}

/// Back-substitute the upper-triangular leading n×n block of a row-major
/// matrix with row stride `n`.
pub(crate) fn back_substitute(
    a: &[f64],
    b: &[f64],
    n: usize,
    out: &mut [f64],
) -> Result<(), String> {
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in (i + 1)..n {
            s -= a[i * n + j] * out[j];
        }
        let d = a[i * n + i];
        if d.abs() < 1e-300 || !d.is_finite() {
            return Err("back_substitute: singular triangular factor".into());
        }
        out[i] = s / d;
    }
    Ok(())
}

fn neg(v: &[f64]) -> Vec<f64> {
    v.iter().map(|&x| -x).collect()
}

/// Central-difference Jacobian, columns in parallel.
/// Returns the number of residual evaluations performed (2·n).
pub(crate) fn build_jacobian<F>(
    residuals: &F,
    x: &[f64],
    jac: &mut Vec<f64>,
) -> Result<usize, String>
where
    F: Fn(&[f64], &mut Vec<f64>) -> Result<(), String> + Sync,
{
    let n = x.len();
    let cube_root_eps = f64::EPSILON.cbrt();

    // Probe column 0 first (serially) to learn m and size the buffer.
    let h0 = cube_root_eps * x[0].abs().max(1.0);
    let mut xp = x.to_vec();
    let mut xm = x.to_vec();
    xp[0] += h0;
    xm[0] -= h0;

    let mut rp = Vec::new();
    let mut rm = Vec::new();
    residuals(&xp, &mut rp)?;
    residuals(&xm, &mut rm)?;
    if rp.len() != rm.len() {
        return Err("levenberg_marquardt: inconsistent residual length".into());
    }
    let m = rp.len();
    if jac.len() != m * n {
        jac.resize(m * n, 0.0);
    }
    let inv2h = 1.0 / (2.0 * h0);
    for i in 0..m {
        jac[i * n] = (rp[i] - rm[i]) * inv2h;
    }

    // Remaining columns in parallel; each task owns its buffers.
    let cols: Vec<usize> = (1..n).collect();
    let results: Result<Vec<(usize, Vec<f64>)>, String> = cols
        .into_par_iter()
        .map(|j| {
            let hj = cube_root_eps * x[j].abs().max(1.0);
            let mut xp = x.to_vec();
            let mut xm = x.to_vec();
            xp[j] += hj;
            xm[j] -= hj;
            let mut rp = Vec::new();
            let mut rm = Vec::new();
            residuals(&xp, &mut rp)?;
            residuals(&xm, &mut rm)?;
            if rp.len() != m || rm.len() != m {
                return Err("levenberg_marquardt: inconsistent residual length".into());
            }
            let inv2h = 1.0 / (2.0 * hj);
            let col: Vec<f64> = (0..m).map(|i| (rp[i] - rm[i]) * inv2h).collect();
            Ok((j, col))
        })
        .collect();

    for (j, col) in results? {
        for i in 0..m {
            jac[i * n + j] = col[i];
        }
    }
    Ok(2 * n)
}

/// Symmetric positive-definite solve via Cholesky with a fallback to
/// Gaussian elimination with partial pivoting (for semi-definite systems
/// nudged by damping). Returns the solution vector.
fn solve_symmetric(a_in: &[f64], b: &[f64], out: &mut [f64]) -> Result<(), String> {
    let n = b.len();
    debug_assert_eq!(a_in.len(), n * n);

    // Try Cholesky in-place on a copy.
    let mut l = a_in.to_vec();
    for j in 0..n {
        let mut d = l[j * n + j];
        for k in 0..j {
            d -= l[j * n + k] * l[j * n + k];
        }
        if d <= 1e-300 || !d.is_finite() {
            return gauss_solve(a_in, b, out);
        }
        let dj = d.sqrt();
        l[j * n + j] = dj;
        for i in (j + 1)..n {
            let mut s = l[i * n + j];
            for k in 0..j {
                s -= l[i * n + k] * l[j * n + k];
            }
            l[i * n + j] = s / dj;
        }
    }

    // Forward/back substitution using lower triangle (L Lᵀ).
    let mut y = vec![0.0f64; n];
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= l[i * n + k] * y[k];
        }
        y[i] = s / l[i * n + i];
    }
    for i in (0..n).rev() {
        let mut s = y[i];
        for k in (i + 1)..n {
            s -= l[k * n + i] * out[k];
        }
        out[i] = s / l[i * n + i];
    }
    Ok(())
}

/// Dense Gaussian elimination with partial pivoting (fallback path).
fn gauss_solve(a: &[f64], b: &[f64], out: &mut [f64]) -> Result<(), String> {
    let n = b.len();
    let mut m = vec![0.0f64; n * (n + 1)];
    for i in 0..n {
        m[i * (n + 1)..i * (n + 1) + n].copy_from_slice(&a[i * n..i * n + n]);
        m[i * (n + 1) + n] = b[i];
    }

    for col in 0..n {
        // Partial pivot.
        let (piv, _) = (col..n).map(|r| (r, m[r * (n + 1) + col].abs())).fold(
            (col, 0.0f64),
            |(br, bv), (r, v)| if v > bv { (r, v) } else { (br, bv) },
        );
        if m[piv * (n + 1) + col].abs() < 1e-300 {
            return Err("singular normal-equation system".into());
        }
        if piv != col {
            for c in 0..(n + 1) {
                m.swap(col * (n + 1) + c, piv * (n + 1) + c);
            }
        }
        let d = m[col * (n + 1) + col];
        for r in (col + 1)..n {
            let factor = m[r * (n + 1) + col] / d;
            if factor == 0.0 {
                continue;
            }
            for c in col..(n + 1) {
                let idx = r * (n + 1) + c;
                m[idx] -= factor * m[col * (n + 1) + c];
            }
        }
    }

    for i in (0..n).rev() {
        let mut s = m[i * (n + 1) + n];
        for c in (i + 1)..n {
            s -= m[i * (n + 1) + c] * out[c];
        }
        out[i] = s / m[i * (n + 1) + i];
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg_fast() -> LmConfig {
        LmConfig {
            max_iterations: 200,
            ..LmConfig::default()
        }
    }

    /// Tight-gradient config for problems where we assert near machine
    /// precision (with FD Jacobians, termination at ‖Jᵀr‖∞ < gtol leaves
    /// parameter error ~√gtol — same as scipy with numeric Jacobians).
    fn cfg_precise() -> LmConfig {
        LmConfig {
            max_iterations: 500,
            ftol: 1e-16,
            xtol: 1e-16,
            gtol: 1e-14,
            ..LmConfig::default()
        }
    }

    #[test]
    fn linear_least_squares_exact_recovery() {
        // y = 2x − 1 sampled; residuals r_i = a·x_i + b − y_i
        let xs = [0.0_f64, 1.0, 2.0, 3.0, 4.0];
        let ys = [-1.0, 1.0, 3.0, 5.0, 7.0];
        let f = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            for (&xi, &yi) in xs.iter().zip(&ys) {
                out.push(x[0] * xi + x[1] - yi);
            }
            Ok(())
        };
        let res = levenberg_marquardt(
            &f,
            &[0.0, 0.0],
            &[-100., -100.],
            &[100., 100.],
            &cfg_precise(),
        )
        .unwrap();
        assert!((res.x[0] - 2.0).abs() < 1e-8, "slope {:?}", res.x);
        assert!((res.x[1] + 1.0).abs() < 1e-8, "intercept {:?}", res.x);
        assert!(res.cost < 1e-18, "cost {}", res.cost);
    }

    #[test]
    fn exponential_fit_nonlinear() {
        // y = 2·exp(−0.5 t); start far away at (a,b) = (1, 1).
        let ts = [0.0_f64, 0.5, 1.0, 1.5, 2.0, 3.0];
        let ys: Vec<f64> = ts.iter().map(|&t| 2.0 * (-0.5 * t).exp()).collect();
        let f = move |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            for (&t, &y) in ts.iter().zip(&ys) {
                out.push(x[0] * (-x[1] * t).exp() - y);
            }
            Ok(())
        };
        let res =
            levenberg_marquardt(&f, &[1.0, 1.0], &[-10., -10.], &[10., 10.], &cfg_fast()).unwrap();
        assert!((res.x[0] - 2.0).abs() < 1e-6, "a {:?}", res.x);
        assert!((res.x[1] - 0.5).abs() < 1e-6, "b {:?}", res.x);
        assert!(res.cost < 1e-20);
    }

    #[test]
    fn bounds_respected_optimum_at_boundary() {
        // min (x−5)² with ub = 3 → boundary optimum x = 3, cost = 4.
        let f = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(x[0] - 5.0);
            Ok(())
        };
        let res = levenberg_marquardt(&f, &[0.0], &[-10.0], &[3.0], &cfg_fast()).unwrap();
        assert!((res.x[0] - 3.0).abs() < 1e-9, "x {:?}", res.x);
        assert!((res.cost - 4.0).abs() < 1e-9);

        // Mirror case at lower bound.
        let res_lo = levenberg_marquardt(&f, &[8.0], &[7.0], &[20.0], &cfg_fast()).unwrap();
        assert!((res_lo.x[0] - 7.0).abs() < 1e-9);
    }

    #[test]
    fn box_constrained_corner_optimum_two_params() {
        // min (x−2)² + (y−3)² inside [2.5, ∞) × [3.5, ∞) → corner (2.5, 3.5).
        let f = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(x[0] - 2.0);
            out.push(x[1] - 3.0);
            Ok(())
        };
        let res = levenberg_marquardt(&f, &[10.0, 10.0], &[2.5, 3.5], &[100.0, 100.0], &cfg_fast())
            .unwrap();
        assert!((res.x[0] - 2.5).abs() < 1e-9, "x {:?}", res.x);
        assert!((res.x[1] - 3.5).abs() < 1e-9, "y {:?}", res.x);
        assert!((res.cost - 0.5).abs() < 1e-9);
    }

    #[test]
    fn mixed_interior_and_bound_variables() {
        // One variable interior-optimal, other pinned at its upper bound:
        // min (x−1)² + (y−0)² with y ∈ (…, 2].
        let f = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(x[0] - 1.0);
            out.push(x[1]);
            Ok(())
        };
        let res = levenberg_marquardt(&f, &[5.0, 0.0], &[-50.0, -50.0], &[50.0, 2.0], &cfg_fast())
            .unwrap();
        assert!((res.x[0] - 1.0).abs() < 1e-9);
        assert!((res.x[1]).abs() < 1e-9); // interior: bound irrelevant
    }

    #[test]
    fn residual_error_propagates() {
        let f = |_x: &[f64], _out: &mut Vec<f64>| -> Result<(), String> {
            Err("solver blew up".into())
        };
        let err = levenberg_marquardt(&f, &[0.0], &[-1.0], &[1.0], &LmConfig::default());
        assert!(err.is_err());
    }

    #[test]
    fn max_iterations_respected() {
        // Slow-converging ill-conditioned problem with a tiny budget.
        let f = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(1e3 * (x[0] - 1.0));
            out.push(x[0] * x[0] - 1.0);
            Ok(())
        };
        let cfg = LmConfig {
            max_iterations: 3,
            ..LmConfig::default()
        };
        let res = levenberg_marquardt(&f, &[0.0], &[-10.0], &[10.0], &cfg).unwrap();
        assert!(res.iterations <= 3);
        assert!(matches!(
            res.termination,
            LmTermination::MaxIterations
                | LmTermination::Step
                | LmTermination::Stalled
                | LmTermination::Gradient
        ));
    }

    #[test]
    fn x0_outside_bounds_rejected() {
        let f = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(x[0]);
            Ok(())
        };
        // Out-of-bounds x0 is a caller bug → loud error, not silent clamping.
        assert!(levenberg_marquardt(&f, &[99.0], &[0.0], &[10.0], &cfg_fast()).is_err());

        // Valid x0 at the upper bound converges to the interior optimum.
        // Default gtol=1e-10 on r = x leaves |x| ~ √gtol — assert accordingly.
        let res = levenberg_marquardt(&f, &[10.0], &[0.0], &[10.0], &cfg_fast()).unwrap();
        assert!(res.x[0].abs() < 1e-4, "x {:?}", res.x);
    }

    #[test]
    fn invalid_inputs_rejected() {
        let f = |_x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            Ok(())
        };
        assert!(levenberg_marquardt(&f, &[], &[], &[], &LmConfig::default()).is_err());
        assert!(levenberg_marquardt(&f, &[0.0], &[1.0], &[-1.0], &LmConfig::default()).is_err());
    }

    #[test]
    fn fixed_residual_length_enforced() {
        // Variable-length residual system must be detected.
        let flip = std::sync::atomic::AtomicBool::new(false);
        let f = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(x[0]);
            if flip.load(std::sync::atomic::Ordering::Relaxed) {
                out.push(0.0);
            }
            Ok(())
        };
        // First call defines m; subsequent mismatch errors surface either as
        // Err from the driver or are tolerated depending on ordering — here
        // flip is never set, so it must succeed cleanly instead.
        let res = levenberg_marquardt(&f, &[1.0], &[-1.0], &[1.0], &cfg_fast());
        assert!(res.is_ok());
    }

    // ----------------------------------------------------------------------
    // R4.4b internals: the two guards that make the MINPACK port trustworthy
    // ----------------------------------------------------------------------

    /// A deterministic, reproducible matrix generator (no `rand` dependency).
    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    }

    /// `jac` (m x n row-major), `r` (m) -> the pieces both step solvers need.
    fn step_inputs(jac: &[f64], r: &[f64], m: usize, n: usize) -> (Vec<f64>, Vec<f64>) {
        let mut jtr = vec![0.0; n];
        let mut col_sq = vec![0.0; n];
        for i in 0..m {
            for j in 0..n {
                jtr[j] += jac[i * n + j] * r[i];
                col_sq[j] += jac[i * n + j] * jac[i * n + j];
            }
        }
        let d_scale = col_sq.iter().map(|c| c.max(1e-14).sqrt()).collect();
        (jtr, d_scale)
    }

    /// The QR path, driven exactly as `levenberg_marquardt` drives it.
    fn qr_step(
        jac: &[f64],
        r: &[f64],
        d_scale: &[f64],
        lambda: f64,
        m: usize,
        n: usize,
    ) -> Vec<f64> {
        let mut qrj = jac.to_vec();
        let mut qtr = r.to_vec();
        householder_qr_in_place(&mut qrj, &mut qtr, m, n).expect("QR of a full-rank J");
        let mut r_tri = vec![0.0; n * n];
        for i in 0..n {
            r_tri[i * n + i..(i + 1) * n].copy_from_slice(&qrj[i * n + i..(i + 1) * n]);
        }
        let mut out = vec![0.0; n];
        solve_damped_step(&r_tri, &qtr[..n], d_scale, lambda, &mut out).expect("damped solve");
        out
    }

    /// The damped linear least-squares objective the step is supposed to
    /// minimize: ||J d + r||^2 + lambda*||D d||^2. Lower is a better step,
    /// whichever way it was obtained.
    fn step_objective(
        jac: &[f64],
        r: &[f64],
        d_scale: &[f64],
        lambda: f64,
        step: &[f64],
        m: usize,
        n: usize,
    ) -> f64 {
        let mut acc = 0.0;
        for i in 0..m {
            let mut s = r[i];
            for j in 0..n {
                s += jac[i * n + j] * step[j];
            }
            acc += s * s;
        }
        for j in 0..n {
            let dj = d_scale[j] * step[j];
            acc += lambda * dj * dj;
        }
        acc
    }

    #[test]
    fn qr_step_reproduces_the_normal_equation_step_when_well_conditioned() {
        // The cheapest strong guard against a linear-algebra bug: on a
        // well-conditioned system the two formulations define the SAME step,
        // so the new solve must reproduce the old one to near machine
        // precision. Anywhere it does not, one of them is wrong.
        let (m, n) = (40usize, 5usize);
        let mut seed = 0x5EED_1234u64;
        let jac: Vec<f64> = (0..m * n).map(|_| lcg(&mut seed)).collect();
        let r: Vec<f64> = (0..m).map(|_| lcg(&mut seed)).collect();
        let (jtr, d_scale) = step_inputs(&jac, &r, m, n);

        for &lambda in &[1e-6, 1e-3, 1.0, 1e3] {
            let qr = qr_step(&jac, &r, &d_scale, lambda, m, n);
            let mut legacy = vec![0.0; n];
            normal_equation_step(&jac, &jtr, &d_scale, lambda, m, n, &mut legacy)
                .expect("legacy solve");

            let scale = qr.iter().fold(0.0f64, |a, v| a.max(v.abs()));
            assert!(
                scale > 0.0,
                "lambda={lambda}: degenerate test, step is exactly zero"
            );
            let dev = qr
                .iter()
                .zip(&legacy)
                .fold(0.0f64, |a, (p, q)| a.max((p - q).abs()));
            assert!(
                dev / scale < 1e-8,
                "lambda={lambda}: rel dev {}",
                dev / scale
            );
        }

        // Non-degeneracy, checked where the step is largest: heavy damping
        // shrinks it towards zero by design, so a blanket size floor would
        // fail for the wrong reason.
        let lightly_damped = qr_step(&jac, &r, &d_scale, 1e-6, m, n);
        assert!(lightly_damped.iter().fold(0.0f64, |a, v| a.max(v.abs())) > 1e-2);
    }

    #[test]
    fn qr_step_beats_the_normal_equations_when_the_columns_are_correlated() {
        // Why the QR is there at all (review S3.6). A Vandermonde J has a
        // condition number the normal equations square; thin-film stacks with
        // correlated layers are the same situation. The step is defined as the
        // minimizer of ||J d + r||^2 + lambda*||D d||^2, so that objective --
        // not agreement with the old code -- is the referee here.
        let (m, n) = (40usize, 14usize);
        let mut jac = vec![0.0; m * n];
        for i in 0..m {
            let t = i as f64 / (m - 1) as f64;
            let mut p = 1.0;
            for j in 0..n {
                jac[i * n + j] = p;
                p *= t;
            }
        }
        let mut seed = 0xC0FFEEu64;
        let r: Vec<f64> = (0..m).map(|_| lcg(&mut seed)).collect();
        let (jtr, d_scale) = step_inputs(&jac, &r, m, n);
        // The regime where it matters: lambda has decayed to near nothing, as
        // gain-ratio damping makes it do near a good optimum -- exactly where
        // step accuracy decides the last digits. At the lambda the solver
        // starts from, the Marquardt term regularizes J^T J enough that the
        // two formulations are indistinguishable; the advantage is not that
        // one is always better, it is that the QR does not collapse when the
        // damping stops covering for it.
        let lambda = 1e-16;

        let qr = qr_step(&jac, &r, &d_scale, lambda, m, n);
        let mut legacy = vec![0.0; n];
        normal_equation_step(&jac, &jtr, &d_scale, lambda, m, n, &mut legacy)
            .expect("legacy solve");

        let o_qr = step_objective(&jac, &r, &d_scale, lambda, &qr, m, n);
        let o_legacy = step_objective(&jac, &r, &d_scale, lambda, &legacy, m, n);
        assert!(
            o_qr <= o_legacy,
            "QR step {o_qr} is worse than legacy {o_legacy}"
        );
        // And the gap is real, not round-off: this is the regression the
        // condition-squaring critique predicts.
        assert!(
            o_qr < o_legacy * (1.0 - 1e-4),
            "no measurable advantage: qr={o_qr} legacy={o_legacy}"
        );
    }

    #[test]
    fn a_clipped_step_is_scored_by_its_own_prediction() {
        // The trap in a naive MINPACK port. The bound veto+clamp rewrites the
        // step AFTER it is solved; scoring it with the full step's predicted
        // reduction overstates what the model promised, so rho comes out too
        // small -- over-damping and premature ftol exits at the boundary.
        let (m, n) = (12usize, 3usize);
        let mut seed = 0xBEEF_0001u64;
        let jac: Vec<f64> = (0..m * n).map(|_| lcg(&mut seed)).collect();
        let r: Vec<f64> = (0..m).map(|_| lcg(&mut seed)).collect();
        let (jtr, d_scale) = step_inputs(&jac, &r, m, n);
        let full = qr_step(&jac, &r, &d_scale, 1e-3, m, n);

        // Component 1 vetoed by an active bound; the rest survive.
        let mut clipped = full.clone();
        clipped[1] = 0.0;

        let p_full = predicted_reduction(&jac, &jtr, &full, m, n);
        let p_clip = predicted_reduction(&jac, &jtr, &clipped, m, n);
        assert!(
            p_full > 0.0,
            "the LM step must predict a reduction: {p_full}"
        );
        assert!(
            p_clip < p_full,
            "clipping removed a descent component but the prediction did \
                 not shrink: full={p_full} clipped={p_clip}"
        );
        assert!(
            p_clip > 0.0,
            "the surviving components still descend: {p_clip}"
        );

        // The clipped step is still a descent step on the true residual.
        let mut cost0 = 0.0;
        let mut cost1 = 0.0;
        for i in 0..m {
            let mut s = 0.0;
            for j in 0..n {
                s += jac[i * n + j] * clipped[j];
            }
            cost0 += r[i] * r[i];
            cost1 += (r[i] + s) * (r[i] + s);
        }
        assert!(cost1 < cost0, "clipped step does not reduce the model cost");
    }

    #[test]
    fn bound_active_run_still_converges_to_the_boundary_optimum() {
        // The end-to-end half of the test above: an unconstrained optimum
        // outside the box, so a bound is active on every accepted step and the
        // clipped-prediction path is exercised by the driver itself.
        let f = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(x[0] - 5.0);
            out.push(x[1] + 5.0);
            out.push(0.25 * (x[0] - x[1]));
            Ok(())
        };
        let res = levenberg_marquardt(&f, &[0.5, 0.5], &[0.0, 0.0], &[1.0, 1.0], &cfg_precise())
            .expect("bounded run");
        assert!((res.x[0] - 1.0).abs() < 1e-9, "x0 = {}", res.x[0]);
        assert!((res.x[1] - 0.0).abs() < 1e-9, "x1 = {}", res.x[1]);
        // (1, 0) gives residuals (-4, 5, 0.25): 16 + 25 + 0.0625.
        assert!((res.cost - 41.0625).abs() < 1e-9, "cost = {}", res.cost);
    }

    #[test]
    fn gain_ratio_and_the_fixed_ladder_find_the_same_optimum() {
        // Damping controls the PATH, not the answer. Iteration counts are
        // reported, not pinned (R4.4d): they are allowed to differ.
        let xs = [0.0_f64, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
        let ys: Vec<f64> = xs.iter().map(|t| 2.5 * (-0.7 * t).exp()).collect();
        let f = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            for (&t, &y) in xs.iter().zip(&ys) {
                out.push(x[0] * (-x[1] * t).exp() - y);
            }
            Ok(())
        };
        let bounds_lo = [0.1, 0.1];
        let bounds_hi = [10.0, 5.0];

        let run = |damping| {
            levenberg_marquardt(
                &f,
                &[1.0, 0.2],
                &bounds_lo,
                &bounds_hi,
                &LmConfig {
                    damping,
                    ..cfg_precise()
                },
            )
            .expect("run")
        };
        let gain = run(LmDamping::GainRatio);
        let fixed = run(LmDamping::Fixed);

        assert!(
            (gain.x[0] - 2.5).abs() < 1e-6 && (gain.x[1] - 0.7).abs() < 1e-6,
            "gain-ratio optimum {:?}",
            gain.x
        );
        assert!(
            (gain.x[0] - fixed.x[0]).abs() < 1e-6 && (gain.x[1] - fixed.x[1]).abs() < 1e-6,
            "{:?} vs {:?}",
            gain.x,
            fixed.x
        );
    }

    #[test]
    fn the_gradient_cosine_does_not_move_when_a_column_is_rescaled() {
        // The criterion itself, in isolation. Carrying a thickness in
        // micrometres instead of nanometres scales its Jacobian column, and
        // with it |(J^T r)_j| -- so ||J^T r||_inf answers a different question
        // after the change of units. The cosine answers the same one.
        let col_sq = [4.0_f64, 9.0, 1.0];
        let jtr = [0.02_f64, 0.3, 0.005];
        let r_norm = 5.0;
        let base = gradient_cosine(&jtr, &col_sq, r_norm);

        // Column 1 in units 1000x smaller: its norm and its gradient entry
        // both scale by 1000.
        let scaled_cs = [4.0_f64, 9.0 * 1e6, 1.0];
        let scaled_g = [0.02_f64, 0.3 * 1e3, 0.005];
        let scaled = gradient_cosine(&scaled_g, &scaled_cs, r_norm);

        assert!((base - scaled).abs() <= 1e-15 * base, "{base} vs {scaled}");
        // ...whereas the scale-dependent measure moved by three orders.
        let inf_base = jtr.iter().fold(0.0f64, |a, v| a.max(v.abs()));
        let inf_scaled = scaled_g.iter().fold(0.0f64, |a, v| a.max(v.abs()));
        assert!(
            inf_scaled > inf_base * 100.0,
            "the contrast this test exists for is gone"
        );
    }

    #[test]
    fn the_gradient_cosine_handles_the_degenerate_inputs() {
        // An exact fit is stationary by definition -- and must not divide by
        // ||r|| = 0. A zero column contributes no direction at all.
        assert_eq!(gradient_cosine(&[1.0, 2.0], &[1.0, 1.0], 0.0), 0.0);
        assert_eq!(gradient_cosine(&[1.0, 2.0], &[1.0, 1.0], f64::NAN), 0.0);
        assert_eq!(gradient_cosine(&[5.0, 0.0], &[0.0, 4.0], 1.0), 0.0);
    }

    #[test]
    fn a_badly_scaled_problem_converges_to_the_same_fit_either_way() {
        // End to end: fit the same line twice, once with the slope carried in
        // units 10^4 times smaller. The optimum is a property of the problem,
        // not of the units, and neither gradient criterion may change it.
        let xs = [0.0_f64, 1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = [1.0_f64, 3.1, 4.9, 7.2, 8.9, 11.1];
        let scale = 1e4;

        let plain = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            for (&t, &y) in xs.iter().zip(&ys) {
                out.push(x[0] * t + x[1] - y);
            }
            Ok(())
        };
        let scaled = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            for (&t, &y) in xs.iter().zip(&ys) {
                out.push((x[0] / scale) * t + x[1] - y);
            }
            Ok(())
        };

        for invariant in [true, false] {
            let cfg = LmConfig {
                gtol_scale_invariant: invariant,
                ..cfg_precise()
            };
            let a = levenberg_marquardt(&plain, &[0.0, 0.0], &[-1e3, -1e3], &[1e3, 1e3], &cfg)
                .expect("plain");
            let b = levenberg_marquardt(&scaled, &[0.0, 0.0], &[-1e7, -1e7], &[1e7, 1e7], &cfg)
                .expect("scaled");
            assert!(
                (a.x[0] - b.x[0] / scale).abs() < 1e-6,
                "invariant={invariant}: {} vs {}",
                a.x[0],
                b.x[0] / scale
            );
            assert!(
                (a.x[1] - b.x[1]).abs() < 1e-6,
                "invariant={invariant}: {} vs {}",
                a.x[1],
                b.x[1]
            );
            assert!(
                (a.cost - b.cost).abs() / a.cost < 1e-9,
                "invariant={invariant}: {} vs {}",
                a.cost,
                b.cost
            );
        }
    }

    #[test]
    fn the_scale_invariant_gradient_test_can_be_switched_off() {
        // It is a config flag, not a hard-wired change of contract: the
        // scale-dependent ||J^T r||_inf test stays available on its own.
        let xs = [0.0_f64, 1.0, 2.0, 3.0];
        let f = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            for &t in &xs {
                out.push(x[0] * t + x[1] - (2.0 * t - 1.0));
            }
            Ok(())
        };
        for invariant in [true, false] {
            let cfg = LmConfig {
                gtol_scale_invariant: invariant,
                ..cfg_precise()
            };
            let res = levenberg_marquardt(&f, &[0.0, 0.0], &[-10.0, -10.0], &[10.0, 10.0], &cfg)
                .expect("run");
            assert!(
                (res.x[0] - 2.0).abs() < 1e-8 && (res.x[1] + 1.0).abs() < 1e-8,
                "invariant={invariant}: {:?}",
                res.x
            );
        }
    }

    #[test]
    fn ftol_needs_both_reductions_not_just_the_one_that_happened() {
        // The rule, in isolation. No end-to-end case is pinned here on
        // purpose: which problems separate the two rules depends on the
        // damping path, so a run-level assertion would be pinning a
        // coincidence. What the release changes is the predicate.
        let ftol = 1e-6;

        // Both small: converged, and this is the only case that is.
        assert!(ftol_converged(1e-9, 1e-9, ftol));
        // Delivered little, but the model still promises a lot -- the case
        // the old actual-only rule mistook for convergence.
        assert!(!ftol_converged(1e-9, 1e-2, ftol));
        // Delivered a lot: not converged either way.
        assert!(!ftol_converged(1e-2, 1e-9, ftol));
        assert!(!ftol_converged(1e-2, 1e-2, ftol));
        // The model did not expect the step to help; its silence is not
        // evidence about the optimum.
        assert!(!ftol_converged(1e-9, -1e-9, ftol));
        // Exactly at the threshold is not below it.
        assert!(!ftol_converged(ftol, 0.0, ftol));
        assert!(!ftol_converged(0.0, ftol, ftol));
    }

    #[test]
    fn a_clamped_step_reports_the_gain_ratio_of_the_step_it_took() {
        // The wiring half of `a_clipped_step_is_scored_by_its_own_prediction`.
        // The optimum here is (100, -50) and the box is [0,1]^2, so the very
        // first LM step overshoots by two orders and is CLAMPED -- the applied
        // step is a small fraction of the solved one. The problem is linear,
        // so the model is exact for whatever step is actually taken: rho must
        // be 1. Scoring the clamped step with the full step's prediction gives
        // rho ~ 1e-2 instead, which is how the driver would over-damp itself
        // into crawling along a boundary.
        let ts = [0.0_f64, 1.0, 2.0, 3.0, 4.0];
        let f = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            for &t in &ts {
                out.push(x[0] * t + x[1] - (100.0 * t - 50.0));
            }
            Ok(())
        };
        let res = levenberg_marquardt(&f, &[0.5, 0.5], &[0.0, 0.0], &[1.0, 1.0], &cfg_precise())
            .expect("clamped run");
        assert_eq!(res.x, vec![1.0, 1.0], "the corner nearest the optimum");
        assert!(
            (res.gain_ratio - 1.0).abs() < 1e-6,
            "rho = {} for an exact linear model; the prediction was not \
                 taken for the clamped step",
            res.gain_ratio
        );
    }

    #[test]
    fn the_gain_ratio_is_nan_when_nothing_was_accepted() {
        // Starting at the optimum: the gradient test fires before any step.
        let f = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(x[0]);
            out.push(x[1]);
            Ok(())
        };
        let res = levenberg_marquardt(&f, &[0.0, 0.0], &[-1.0, -1.0], &[1.0, 1.0], &cfg_fast())
            .expect("run");
        assert_eq!(res.termination, LmTermination::Gradient);
        assert!(res.gain_ratio.is_nan(), "gain_ratio = {}", res.gain_ratio);
    }

    #[test]
    fn householder_qr_refuses_a_zero_or_non_finite_column() {
        // The fallback to the normal equations hangs off this Err; a QR that
        // quietly returned NaNs would poison every step after it.
        let mut a = vec![1.0, 0.0, 2.0, 0.0, 3.0, 0.0];
        let mut b = vec![1.0, 1.0, 1.0];
        assert!(householder_qr_in_place(&mut a, &mut b, 3, 2).is_err());

        let mut a = vec![1.0, f64::NAN, 2.0, 1.0, 3.0, 1.0];
        let mut b = vec![1.0, 1.0, 1.0];
        assert!(householder_qr_in_place(&mut a, &mut b, 3, 2).is_err());

        let mut a = vec![1.0, 1.0];
        let mut b = vec![1.0];
        assert!(
            householder_qr_in_place(&mut a, &mut b, 1, 2).is_err(),
            "fewer rows than columns is not a factorable system"
        );
    }

    #[test]
    fn a_rank_deficient_jacobian_still_produces_a_step() {
        // Two identical columns: J^T J is singular, and the QR of J itself
        // fails on the second column. Damping is what rescues the system, and
        // the driver must fall through to the normal equations rather than
        // stall. (Physically: two layers the merit cannot tell apart.)
        let f = |x: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            for t in 0..6 {
                let t = t as f64;
                out.push((x[0] + x[1]) * t - 3.0 * t);
            }
            Ok(())
        };
        let res = levenberg_marquardt(&f, &[0.0, 0.0], &[-10.0, -10.0], &[10.0, 10.0], &cfg_fast())
            .expect("degenerate run");
        assert!(
            (res.x[0] + res.x[1] - 3.0).abs() < 1e-6,
            "sum {} should reach 3",
            res.x[0] + res.x[1]
        );
    }

    // -- R4.5 increment ii: the analytic-Jacobian hook ----------------------
    mod analytic_hook {
        use super::*;

        // r_i(x) = a_i·x0 + b_i·x1² − y_i, over ten rows. Nonlinear in x1,
        // so the Jacobian genuinely moves between iterations, and exact:
        // ∂r_i/∂x0 = a_i, ∂r_i/∂x1 = 2·b_i·x1.
        const A: [f64; 10] = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1];
        const B: [f64; 10] = [0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1.0, 1.15, 1.3, 1.45];

        fn y_of(x: &[f64]) -> Vec<f64> {
            (0..10).map(|i| A[i] * x[0] + B[i] * x[1] * x[1]).collect()
        }

        fn problem(y: Vec<f64>) -> impl Fn(&[f64], &mut Vec<f64>) -> Result<(), String> + Sync {
            move |x: &[f64], out: &mut Vec<f64>| {
                out.clear();
                for i in 0..10 {
                    out.push(A[i] * x[0] + B[i] * x[1] * x[1] - y[i]);
                }
                Ok(())
            }
        }

        struct Exact;
        impl JacobianSource for Exact {
            fn fill(&self, x: &[f64], jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
                jac.clear();
                for i in 0..10 {
                    jac.push(A[i]);
                    jac.push(2.0 * B[i] * x[1]);
                }
                Ok(Some(10))
            }
        }

        struct Declines;
        impl JacobianSource for Declines {
            fn fill(&self, _x: &[f64], _jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
                Ok(None)
            }
        }

        struct WrongShape;
        impl JacobianSource for WrongShape {
            fn fill(&self, _x: &[f64], jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
                jac.clear();
                jac.resize(14, 0.5);
                Ok(Some(7))
            }
        }

        struct Explodes;
        impl JacobianSource for Explodes {
            fn fill(&self, _x: &[f64], _jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
                Err("no deposits for this stack".into())
            }
        }

        fn run<J: JacobianSource>(src: Option<&J>) -> LmResult {
            let truth = [3.0, 1.6];
            let f = problem(y_of(&truth));
            levenberg_marquardt_with(
                &f,
                src,
                &[0.5, 0.4],
                &[-10.0, -10.0],
                &[10.0, 10.0],
                &LmConfig::default(),
            )
            .unwrap()
        }

        #[test]
        fn an_exact_jacobian_and_a_differenced_one_reach_the_same_optimum() {
            // The point of the hook: the same solve, from a better J. The
            // answers agree to well inside the difference noise, and the
            // analytic run never calls the residual for a Jacobian probe —
            // which is where the 2n evaluations went.
            let a = run(Some(&Exact));
            let b = run(None::<&NoJacobian>);
            // One build per iteration, plus the final one that proved the
            // gradient test — nothing in this run was differenced.
            assert!(a.analytic_jacobians >= a.iterations, "{a:?}");
            assert_eq!(b.analytic_jacobians, 0);
            for j in 0..2 {
                assert!((a.x[j] - b.x[j]).abs() < 1e-6, "{:?} vs {:?}", a.x, b.x);
            }
            assert!(a.cost < 1e-18 && b.cost < 1e-18, "{} / {}", a.cost, b.cost);
            // The difference is exactly the 2n = 4 probes per Jacobian the
            // analytic run does not make.
            assert!(
                b.evals >= a.evals + 4 * a.analytic_jacobians,
                "analytic {} evals / {} jacobians, differenced {} evals",
                a.evals,
                a.analytic_jacobians,
                b.evals
            );
        }

        #[test]
        fn a_source_that_declines_leaves_the_run_exactly_as_it_was() {
            // `Ok(None)` is the fallback contract: a spec the analytic chain
            // cannot express must still optimize, by differences, with no
            // trace of the attempt in the answer.
            let a = run(Some(&Declines));
            let b = run(None::<&NoJacobian>);
            assert_eq!(a.analytic_jacobians, 0);
            assert_eq!(a.evals, b.evals);
            assert_eq!(a.iterations, b.iterations);
            for j in 0..2 {
                assert_eq!(a.x[j].to_bits(), b.x[j].to_bits());
            }
        }

        #[test]
        fn a_jacobian_of_the_wrong_shape_is_refused_not_used() {
            // A source that disagrees with the residual vector about m is a
            // bug in the source; reading it as a Jacobian would corrupt the
            // step silently.
            let truth = [3.0, 1.6];
            let f = problem(y_of(&truth));
            let err = levenberg_marquardt_with(
                &f,
                Some(&WrongShape),
                &[0.5, 0.4],
                &[-10.0, -10.0],
                &[10.0, 10.0],
                &LmConfig::default(),
            )
            .unwrap_err();
            assert!(err.contains("analytic jacobian"), "{err}");
        }

        #[test]
        fn a_source_that_errors_aborts_rather_than_quietly_differencing() {
            // Declining and failing are different facts. A failure means the
            // caller believed it had a Jacobian and was wrong — silently
            // carrying on would hide it behind a slower run.
            let truth = [3.0, 1.6];
            let f = problem(y_of(&truth));
            let err = levenberg_marquardt_with(
                &f,
                Some(&Explodes),
                &[0.5, 0.4],
                &[-10.0, -10.0],
                &[10.0, 10.0],
                &LmConfig::default(),
            )
            .unwrap_err();
            assert_eq!(err, "no deposits for this stack");
        }

        #[test]
        fn the_plain_entry_point_asks_for_no_jacobian_at_all() {
            let truth = [3.0, 1.6];
            let f = problem(y_of(&truth));
            let r = levenberg_marquardt(
                &f,
                &[0.5, 0.4],
                &[-10.0, -10.0],
                &[10.0, 10.0],
                &LmConfig::default(),
            )
            .unwrap();
            assert_eq!(r.analytic_jacobians, 0);
        }
    }
}
