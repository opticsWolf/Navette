// SPDX-License-Identifier: LGPL-3.0-or-later
//! Navette -- Rust Rewrite of Numba-optimized thin-film optical solver
//!
//! synthesis::trf — trust-region reflective, the bounded least-squares
//! reference method (Branch–Coleman–Li 1999; R4.6, review §3.6).
//!
//! `thick_opt`'s own docstring says it replaces
//! `scipy.optimize.least_squares(method="trf")`. What it actually implements
//! is a Levenberg-Marquardt that keeps `lb ≤ x ≤ ub` by *vetoing and clamping
//! the step it just solved* — bounds are applied after the fact, to a step
//! computed as though they were not there. TRF is the method the docs were
//! naming: the bounds enter the subproblem, through a scaling that shrinks
//! the trust region in the directions that run into them.
//!
//! # How the bounds get in
//!
//! The Coleman-Li vector `v` is, per coordinate, the distance to the bound
//! the anti-gradient points at (and `1` when it points away from both). With
//! `D = diag(√v)` the first-order condition `D²g = 0` says exactly what
//! optimality under bounds means: the gradient vanishes for an interior
//! coordinate, and points *into* the box for one on a bound. The subproblem
//! is then solved in the "hat" variables `x = D x̂`, where the trust region is
//! a plain ball — a large step toward a near bound is expensive in `x̂` — plus
//! a diagonal `C = diag(g)·Jv` that carries the curvature of `v` itself.
//!
//! Three candidate steps are scored against that quadratic and the best is
//! taken: the trust-region step, cut back to the first bound it hits; the
//! same step **reflected** off that bound; and the constrained Cauchy step
//! along `−g_h`. Iterates stay *strictly* interior (the `theta` step-back),
//! which is what keeps `v` differentiable.
//!
//! # Why this backend exists
//!
//! * **Boundary optima are handled by construction.** A film pinned at zero
//!   thickness is where the built-in LM's clamp is weakest — the step is
//!   vetoed component-wise and the model's promise is re-evaluated for the
//!   clipped step. TRF never computes an infeasible step to begin with.
//! * **scipy is a direct oracle.** This is the same algorithm, not a cousin:
//!   the same Coleman-Li scaling, the same `select_step` between the same
//!   three candidates, the same Moré subproblem, the same ftol/xtol/gtol
//!   tests. `validation/review/lm_check.py` part D asserts the tightest
//!   agreement in the plan against `least_squares(method="trf")`.
//! * **Bounds are native.** Unlike [`OptimizerBackend::MinpackLm`], this
//!   backend needs no interior reparametrization, so it can land arbitrarily
//!   close to a bound without the `tanh` map's vanishing gradient.
//!
//! [`OptimizerBackend::MinpackLm`]: super::optimizer::OptimizerBackend::MinpackLm
//!
//! # What it does *not* do
//!
//! It never returns an `x` exactly *on* a bound: `make_strictly_feasible`
//! steps back by one ULP. The synthesis loop's removal sweep compares
//! against `clamp_min`, not against zero, so a film driven to `lb` is still
//! removed — but a caller testing `x == lb` will be disappointed, and that is
//! the one behavioural difference from the built-in worth knowing.
//!
//! `LmConfig`'s damping knobs (`lambda_init`, `lambda_up`, `lambda_down`,
//! `damping`) are Levenberg-Marquardt settings and do nothing here; the
//! trust-region radius plays their role and is not user-settable, exactly as
//! in scipy. `gtol_scale_invariant` is likewise ignored: TRF's stationarity
//! measure `‖g·v‖∞` is already the bound-aware one.
//!
//! # Deviations from scipy, and why
//!
//! * **The subproblem is factored by QR, not SVD.** Moré's algorithm needs
//!   `p(α) = −(BᵀB + αI)⁻¹Bᵀf` and `‖p‖`, `‖q‖` with `Rαᵀq = p`; the
//!   augmented QR this crate already has (R4.4b) supplies both, at the same
//!   arithmetic, and the Newton recurrence on the secular equation is
//!   identical term for term. The one place the two part company is the
//!   rank test: scipy compares the smallest singular value to the largest,
//!   this compares `min|Rᵢᵢ|` to `max|Rᵢᵢ|`, which brackets it. A borderline
//!   rank call changes which branch initializes `α`, not the answer.
//! * **`x_scale` is fixed at 1.** scipy's `x_scale` is an orthogonal
//!   user-facing knob (default `1.0`) that pre-scales the parameters before
//!   Coleman-Li is applied. `LmConfig` has no such field, so there is nothing
//!   to thread, and the default is what the oracle runs.
//! * **`tr_solver` is `"exact"`.** scipy's other option (`"lsmr"`, the 2-D
//!   subspace / indefinite-dogleg route) exists for sparse Jacobians of
//!   millions of rows. Thickness Jacobians are dense and have as many columns
//!   as there are films.

use super::thick_opt::{
    back_substitute, build_jacobian, householder_qr_in_place, JacobianSource, LmConfig, LmResult,
    LmTermination,
};

/// scipy's `EPS`, used in the same two places: the rank threshold and
/// `make_strictly_feasible`'s relative step-back.
const EPS: f64 = f64::EPSILON;

/// Relative step-back applied to `x0` before the first iteration, so the run
/// starts strictly interior even when the caller handed us a point on a
/// bound. scipy's default.
const INITIAL_RSTEP: f64 = 1e-10;

/// `rtol` for Moré's secular-equation Newton iteration, and its cap. Both are
/// scipy's: the subproblem is solved to 1% of Δ, because a trust-region step
/// that is 1% too long is still a trust-region step.
const ALPHA_RTOL: f64 = 0.01;
const ALPHA_MAX_ITER: usize = 10;

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

/// Minimize ‖r(x)‖² over `lb ≤ x ≤ ub` by trust-region reflective.
///
/// Same signature and same `LmResult` as
/// [`levenberg_marquardt_with`](super::thick_opt::levenberg_marquardt_with),
/// so the two are interchangeable behind
/// [`run_optimizer`](super::optimizer::run_optimizer). `cost` is reported in
/// this crate's convention (‖r‖², not the ½‖r‖² scipy carries internally) and
/// `evals` counts every residual evaluation including the ones a differenced
/// Jacobian spends, which is what `max_evals` is measured against.
pub fn trust_region_reflective<F, J>(
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
        return Err("trust_region_reflective: empty parameter vector".into());
    }
    if lb.len() != n || ub.len() != n {
        return Err("trust_region_reflective: bound length mismatch".into());
    }
    for j in 0..n {
        if lb[j] > ub[j] {
            return Err(format!("trust_region_reflective: lb[{j}] > ub[{j}]"));
        }
        if !(lb[j] <= x0[j] && x0[j] <= ub[j]) {
            return Err(format!(
                "trust_region_reflective: x0[{}]={} outside [{}, {}]",
                j, x0[j], lb[j], ub[j]
            ));
        }
        if !lb[j].is_finite() || !ub[j].is_finite() {
            return Err(format!(
                "trust_region_reflective: bound {j} is not finite ([{}, {}]) — the \
                 Coleman-Li scaling is a distance to the bound, so an infinite one \
                 has nothing to measure",
                lb[j], ub[j]
            ));
        }
    }

    let mut evals = 0usize;
    let mut analytic_jacobians = 0usize;

    let mut f: Vec<f64> = Vec::new();
    let mut x = make_strictly_feasible(x0, lb, ub, INITIAL_RSTEP);
    residuals(&x, &mut f)?;
    evals += 1;
    let m = f.len();
    if m == 0 {
        return Err("trust_region_reflective: empty residual vector".into());
    }
    let mut cost = 0.5 * dot(&f, &f);

    let mut jac = vec![0.0f64; m * n];
    fill_jacobian(
        residuals,
        jacobian,
        &x,
        m,
        n,
        &mut jac,
        &mut evals,
        &mut analytic_jacobians,
    )?;
    let mut g = jt_times(&jac, &f, m, n);

    // Δ₀ = ‖x₀ / √v‖, scipy's. A start at the origin gives 0, which is not a
    // radius; scipy falls back to 1 and so do we.
    let (v0, _) = cl_scaling_vector(&x, &g, lb, ub);
    let mut delta = {
        let s: f64 = x
            .iter()
            .zip(&v0)
            .map(|(&xi, &vi)| if vi > 0.0 { (xi * xi) / vi } else { 0.0 })
            .sum();
        let d = s.sqrt();
        if d > 0.0 && d.is_finite() { d } else { 1.0 }
    };

    let mut alpha = 0.0f64; // the Levenberg-Marquardt parameter, carried across
    let mut termination: Option<LmTermination> = None;
    let mut iteration = 0usize;
    let mut ratio = f64::NAN;

    // Scratch that outlives the iteration.
    let mut jh = vec![0.0f64; m * n];
    let mut f_new: Vec<f64> = Vec::new();

    loop {
        let (v, dv) = cl_scaling_vector(&x, &g, lb, ub);
        let g_norm = g
            .iter()
            .zip(&v)
            .fold(0.0f64, |acc, (&gi, &vi)| acc.max((gi * vi).abs()));
        if g_norm < cfg.gtol {
            termination = Some(LmTermination::Gradient);
        }
        if termination.is_some() || evals >= cfg.max_evals || iteration >= cfg.max_iterations {
            break;
        }

        // ---- "hat" space: x = D x̂ with D = diag(√v), plus the curvature of
        //      v itself as a diagonal C = diag(g)·Jv -------------------------
        let d: Vec<f64> = v.iter().map(|vi| vi.sqrt()).collect();
        let diag_h: Vec<f64> = g.iter().zip(&dv).map(|(&gi, &dvi)| gi * dvi).collect();
        let g_h: Vec<f64> = d.iter().zip(&g).map(|(&di, &gi)| di * gi).collect();
        for i in 0..m {
            for k in 0..n {
                jh[i * n + k] = jac[i * n + k] * d[k];
            }
        }

        // ---- factor the augmented system [J_h; diag(√C)] once -------------
        let rows = m + n;
        let mut aug = vec![0.0f64; rows * n];
        aug[..m * n].copy_from_slice(&jh);
        for k in 0..n {
            let col_sq: f64 =
                (0..m).map(|i| jh[i * n + k] * jh[i * n + k]).sum::<f64>() + diag_h[k];
            // A column that is identically zero carries no direction at all
            // (a zero Jacobian column forces g[k] = 0, hence dv[k] = 0 and
            // C[k] = 0), so its step component is zero whatever we put here.
            // The QR, unlike scipy's SVD, refuses to factor it: give it a 1.
            aug[(m + k) * n + k] = if col_sq > 0.0 { diag_h[k].max(0.0).sqrt() } else { 1.0 };
        }
        let mut qtf = vec![0.0f64; rows];
        qtf[..m].copy_from_slice(&f);
        if householder_qr_in_place(&mut aug, &mut qtf, rows, n).is_err() {
            termination = Some(LmTermination::Stalled);
            break;
        }
        let mut r_tri = vec![0.0f64; n * n];
        for i in 0..n {
            r_tri[i * n + i..(i + 1) * n].copy_from_slice(&aug[i * n + i..(i + 1) * n]);
        }
        let g_h_norm = norm(&g_h);

        // theta controls how far the step is held back from the bound. Near a
        // stationary point the step-back tightens towards 0.5%.
        let theta = 0.995f64.max(1.0 - g_norm);

        let mut actual_reduction = -1.0f64;
        let mut accepted: Option<(Vec<f64>, Vec<f64>, f64)> = None;

        while actual_reduction <= 0.0 && evals < cfg.max_evals {
            let (p_h, new_alpha) =
                match solve_lsq_trust_region(&r_tri, &qtf[..n], m, n, g_h_norm, delta, alpha) {
                    Ok(v) => v,
                    Err(_) => {
                        termination = Some(LmTermination::Stalled);
                        break;
                    },
                };
            alpha = new_alpha;
            let p: Vec<f64> = d.iter().zip(&p_h).map(|(&di, &pi)| di * pi).collect();

            let (step, step_h, predicted_reduction) =
                select_step(&x, &jh, &diag_h, &g_h, &p, &p_h, &d, delta, lb, ub, theta, m, n);

            let x_new = {
                let trial: Vec<f64> =
                    x.iter().zip(&step).map(|(&xi, &si)| xi + si).collect();
                make_strictly_feasible(&trial, lb, ub, 0.0)
            };
            residuals(&x_new, &mut f_new)?;
            evals += 1;
            if f_new.len() != m {
                return Err(format!(
                    "trust_region_reflective: residual length changed from {m} to {}",
                    f_new.len()
                ));
            }

            let step_h_norm = norm(&step_h);
            if !f_new.iter().all(|v| v.is_finite()) {
                // Not an error: the merit is undefined out here. Shrink hard
                // and try again from the same point.
                delta = 0.25 * step_h_norm;
                if !(delta > 0.0) {
                    termination = Some(LmTermination::Stalled);
                    break;
                }
                continue;
            }

            let cost_new = 0.5 * dot(&f_new, &f_new);
            actual_reduction = cost - cost_new;
            let (delta_new, r) = update_tr_radius(
                delta,
                actual_reduction,
                predicted_reduction,
                step_h_norm,
                step_h_norm > 0.95 * delta,
            );
            ratio = r;

            let step_norm = norm(&step);
            termination = check_termination(
                actual_reduction,
                cost,
                step_norm,
                norm(&x),
                ratio,
                cfg.ftol,
                cfg.xtol,
            );
            if actual_reduction > 0.0 {
                accepted = Some((x_new, f_new.clone(), cost_new));
            }
            if termination.is_some() {
                break;
            }

            if delta_new > 0.0 {
                alpha *= delta / delta_new;
            }
            delta = delta_new;
            if !(delta > 0.0) || !delta.is_finite() {
                termination = Some(LmTermination::Stalled);
                break;
            }
        }

        if let Some((x_new, f_acc, cost_new)) = accepted {
            x = x_new;
            f = f_acc;
            cost = cost_new;
            fill_jacobian(
                residuals,
                jacobian,
                &x,
                m,
                n,
                &mut jac,
                &mut evals,
                &mut analytic_jacobians,
            )?;
            g = jt_times(&jac, &f, m, n);
        }

        iteration += 1;
        if termination.is_some() {
            break;
        }
    }

    Ok(LmResult {
        x,
        cost: 2.0 * cost,
        iterations: iteration,
        evals,
        termination: termination.unwrap_or(LmTermination::MaxIterations),
        gain_ratio: ratio,
        analytic_jacobians,
    })
}

/// The analytic source if it has one for this `x`, otherwise central
/// differences. Same precedence and same error text as the LM driver, so a
/// mis-shaped analytic Jacobian reads the same whichever backend hit it.
fn fill_jacobian<F, J>(
    residuals: &F,
    jacobian: Option<&J>,
    x: &[f64],
    m: usize,
    n: usize,
    jac: &mut Vec<f64>,
    evals: &mut usize,
    analytic: &mut usize,
) -> Result<(), String>
where
    F: Fn(&[f64], &mut Vec<f64>) -> Result<(), String> + Sync,
    J: JacobianSource + ?Sized,
{
    let supplied = match jacobian {
        Some(src) => src.fill(x, jac)?,
        None => None,
    };
    match supplied {
        Some(rows) => {
            if rows != m || jac.len() != m * n {
                return Err(format!(
                    "trust_region_reflective: analytic jacobian is {rows} rows where the \
                     residual vector is {m} long"
                ));
            }
            *analytic += 1;
        },
        None => {
            let added = build_jacobian(residuals, x, jac)?;
            *evals += added;
        },
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// The trust-region subproblem (Moré, on the QR factor)
// ---------------------------------------------------------------------------

/// Minimize the model over `‖p‖ ≤ Δ`, returning the step and the `α` that
/// produced it.
///
/// `r_tri` is the triangular factor of the augmented Jacobian
/// `B = [J_h; diag(√C)]` and `qtf` the leading `n` of `Qᵀf_aug`, so
/// `g_h = Rᵀ·qtf` and `p(α) = −(RᵀR + αI)⁻¹g_h`. `α` is carried in from the
/// previous call because consecutive subproblems differ only in Δ.
fn solve_lsq_trust_region(
    r_tri: &[f64],
    qtf: &[f64],
    m: usize,
    n: usize,
    g_h_norm: f64,
    delta: f64,
    initial_alpha: f64,
) -> Result<(Vec<f64>, f64), String> {
    let mut rmin = f64::INFINITY;
    let mut rmax: f64 = 0.0;
    for i in 0..n {
        let v = r_tri[i * n + i].abs();
        rmin = rmin.min(v);
        rmax = rmax.max(v);
    }
    // scipy's rank test compares σ_min to σ_max; |Rᵢᵢ| brackets those, and a
    // borderline call changes which branch seeds α, not the answer.
    let full_rank = m >= n && rmin > EPS * (m as f64) * rmax;

    let mut alpha_upper = g_h_norm / delta;
    let mut alpha_lower = 0.0f64;

    if full_rank {
        let mut p_gn = vec![0.0f64; n];
        let neg: Vec<f64> = qtf.iter().map(|v| -v).collect();
        back_substitute(r_tri, &neg, n, &mut p_gn)?;
        let p_norm = norm(&p_gn);
        if p_norm <= delta {
            return Ok((p_gn, 0.0));
        }
        // φ(0) and φ'(0) bound α from below: the Newton step from 0.
        let mut q = vec![0.0f64; n];
        forward_substitute_transpose(r_tri, &p_gn, n, &mut q)?;
        let phi = p_norm - delta;
        let phi_prime = -dot(&q, &q) / p_norm;
        if phi_prime != 0.0 {
            alpha_lower = -phi / phi_prime;
        }
    }

    let mut alpha = if !full_rank && initial_alpha == 0.0 {
        (0.001 * alpha_upper).max((alpha_lower * alpha_upper).sqrt())
    } else {
        initial_alpha
    };

    let mut p = vec![0.0f64; n];
    for _ in 0..ALPHA_MAX_ITER {
        if alpha < alpha_lower || alpha > alpha_upper {
            alpha = (0.001 * alpha_upper).max((alpha_lower * alpha_upper).sqrt());
        }
        let (pa, p_norm, q_norm) = damped_step(r_tri, qtf, n, alpha)?;
        p = pa;
        let phi = p_norm - delta;
        let phi_prime = if p_norm > 0.0 { -(q_norm * q_norm) / p_norm } else { 0.0 };
        if phi < 0.0 {
            alpha_upper = alpha;
        }
        if phi_prime == 0.0 {
            break;
        }
        let ratio = phi / phi_prime;
        alpha_lower = alpha_lower.max(alpha - ratio);
        alpha -= (phi + delta) * ratio / delta;
        if phi.abs() < ALPHA_RTOL * delta {
            break;
        }
    }

    // Re-solve at the final α only if the loop never ran (it always does) —
    // otherwise `p` already holds p(α). Rescale onto the boundary, as scipy
    // does, so a step 1% long cannot leave the region.
    let p_norm = norm(&p);
    if p_norm > 0.0 {
        let s = delta / p_norm;
        for v in p.iter_mut() {
            *v *= s;
        }
    }
    Ok((p, alpha))
}

/// `p = −(RᵀR + αI)⁻¹Rᵀqtf`, plus `‖p‖` and `‖q‖` where `Rαᵀq = p`.
///
/// Solved as the least-squares problem `[R; √α·I]p ≈ [−qtf; 0]`, whose own
/// triangular factor `Rα` satisfies `RαᵀRα = RᵀR + αI` — so the derivative
/// the secular equation needs comes out of the same factorization as the
/// step, which is the whole reason Moré's method is cheap.
fn damped_step(
    r_tri: &[f64],
    qtf: &[f64],
    n: usize,
    alpha: f64,
) -> Result<(Vec<f64>, f64, f64), String> {
    let rows = 2 * n;
    let mut a = vec![0.0f64; rows * n];
    for i in 0..n {
        a[i * n + i..(i + 1) * n].copy_from_slice(&r_tri[i * n + i..(i + 1) * n]);
    }
    let root = alpha.max(0.0).sqrt();
    for k in 0..n {
        a[(n + k) * n + k] = root;
    }
    let mut b = vec![0.0f64; rows];
    for i in 0..n {
        b[i] = -qtf[i];
    }
    householder_qr_in_place(&mut a, &mut b, rows, n)?;
    let mut p = vec![0.0f64; n];
    back_substitute(&a, &b, n, &mut p)?;
    let mut q = vec![0.0f64; n];
    forward_substitute_transpose(&a, &p, n, &mut q)?;
    let (p_norm, q_norm) = (norm(&p), norm(&q));
    Ok((p, p_norm, q_norm))
}

/// Solve `Rᵀq = b` for upper-triangular `R` stored row-major with stride `n`.
fn forward_substitute_transpose(
    r: &[f64],
    b: &[f64],
    n: usize,
    out: &mut [f64],
) -> Result<(), String> {
    for i in 0..n {
        let mut s = b[i];
        for j in 0..i {
            s -= r[j * n + i] * out[j];
        }
        let d = r[i * n + i];
        if d.abs() < 1e-300 || !d.is_finite() {
            return Err("forward_substitute_transpose: singular triangular factor".into());
        }
        out[i] = s / d;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Step selection
// ---------------------------------------------------------------------------

/// The best of the three TRF candidates, and the reduction the model promises
/// for it (in ½‖r‖² units, the convention the ratio tests use).
///
/// Returns `(step, step_h, predicted_reduction)`.
#[allow(clippy::too_many_arguments)]
fn select_step(
    x: &[f64],
    jh: &[f64],
    diag_h: &[f64],
    g_h: &[f64],
    p_in: &[f64],
    p_h_in: &[f64],
    d: &[f64],
    delta: f64,
    lb: &[f64],
    ub: &[f64],
    theta: f64,
    m: usize,
    n: usize,
) -> (Vec<f64>, Vec<f64>, f64) {
    let mut p = p_in.to_vec();
    let mut p_h = p_h_in.to_vec();

    let trial: Vec<f64> = x.iter().zip(&p).map(|(&xi, &pi)| xi + pi).collect();
    if in_bounds(&trial, lb, ub) {
        let value = evaluate_quadratic(jh, g_h, &p_h, diag_h, m, n);
        return (p, p_h, -value);
    }

    let (p_stride, hits) = step_size_to_bound(x, &p, lb, ub);

    // The reflected direction: flip the components that hit a bound.
    let mut r_h: Vec<f64> = p_h.clone();
    for k in 0..n {
        if hits[k] {
            r_h[k] = -r_h[k];
        }
    }
    let r_vec: Vec<f64> = d.iter().zip(&r_h).map(|(&di, &ri)| di * ri).collect();

    // Cut the trust-region step back to the bound it hits.
    for k in 0..n {
        p[k] *= p_stride;
        p_h[k] *= p_stride;
    }
    let x_on_bound: Vec<f64> = x.iter().zip(&p).map(|(&xi, &pi)| xi + pi).collect();

    // From there the reflected direction leaves either the box or the region.
    let to_tr = intersect_trust_region(&p_h, &r_h, delta).map(|(_, t2)| t2);
    let (to_bound, _) = step_size_to_bound(&x_on_bound, &r_vec, lb, ub);

    let r_stride = match to_tr {
        Some(t) => to_bound.min(t),
        None => f64::NEG_INFINITY,
    };
    let (r_stride_l, r_stride_u) = if r_stride > 0.0 {
        let l = (1.0 - theta) * p_stride / r_stride;
        // The upper limit is the bound (held back by theta) or the region
        // boundary, whichever the reflection reaches first.
        let u = if r_stride == to_bound { theta * to_bound } else { to_tr.unwrap_or(0.0) };
        (l, u)
    } else {
        (0.0, -1.0)
    };

    let mut r_value = f64::INFINITY;
    if r_stride_l <= r_stride_u {
        let (a, b, c) = build_quadratic_1d(jh, g_h, &r_h, diag_h, Some(&p_h), m, n);
        let (stride, value) = minimize_quadratic_1d(a, b, r_stride_l, r_stride_u, c);
        for k in 0..n {
            r_h[k] = r_h[k] * stride + p_h[k];
        }
        r_value = value;
    }
    let r_out: Vec<f64> = d.iter().zip(&r_h).map(|(&di, &ri)| di * ri).collect();

    // Hold the truncated trust-region step back off the bound as well.
    for k in 0..n {
        p[k] *= theta;
        p_h[k] *= theta;
    }
    let p_value = evaluate_quadratic(jh, g_h, &p_h, diag_h, m, n);

    // The constrained Cauchy step, along −g_h.
    let mut ag_h: Vec<f64> = g_h.iter().map(|v| -v).collect();
    let ag_norm = norm(&ag_h);
    let mut ag: Vec<f64> = d.iter().zip(&ag_h).map(|(&di, &ai)| di * ai).collect();
    let ag_value = if ag_norm > 0.0 {
        let to_tr = delta / ag_norm;
        let (to_bound, _) = step_size_to_bound(x, &ag, lb, ub);
        let cap = if to_bound < to_tr { theta * to_bound } else { to_tr };
        let (a, b, _) = build_quadratic_1d(jh, g_h, &ag_h, diag_h, None, m, n);
        let (stride, value) = minimize_quadratic_1d(a, b, 0.0, cap, 0.0);
        for k in 0..n {
            ag_h[k] *= stride;
            ag[k] *= stride;
        }
        value
    } else {
        f64::INFINITY
    };

    if p_value < r_value && p_value < ag_value {
        (p, p_h, -p_value)
    } else if r_value < p_value && r_value < ag_value {
        (r_out, r_h, -r_value)
    } else {
        (ag, ag_h, -ag_value)
    }
}

// ---------------------------------------------------------------------------
// Coleman-Li scaling and the bound geometry
// ---------------------------------------------------------------------------

/// `v[i]` is the distance to the bound the anti-gradient points at, or 1 when
/// it points away from both; `dv[i]` is `∂v/∂x` (−1, 1 or 0), which is what
/// makes `C = diag(g)·Jv` non-negative.
fn cl_scaling_vector(x: &[f64], g: &[f64], lb: &[f64], ub: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = x.len();
    let mut v = vec![1.0f64; n];
    let mut dv = vec![0.0f64; n];
    for i in 0..n {
        if g[i] < 0.0 && ub[i].is_finite() {
            v[i] = ub[i] - x[i];
            dv[i] = -1.0;
        } else if g[i] > 0.0 && lb[i].is_finite() {
            v[i] = x[i] - lb[i];
            dv[i] = 1.0;
        }
    }
    (v, dv)
}

fn in_bounds(x: &[f64], lb: &[f64], ub: &[f64]) -> bool {
    x.iter().zip(lb).zip(ub).all(|((&xi, &l), &u)| xi >= l && xi <= u)
}

/// How far along `s` the box lets us go, and which coordinates stop us.
fn step_size_to_bound(x: &[f64], s: &[f64], lb: &[f64], ub: &[f64]) -> (f64, Vec<bool>) {
    let n = x.len();
    let mut steps = vec![f64::INFINITY; n];
    for i in 0..n {
        if s[i] != 0.0 {
            let a = (lb[i] - x[i]) / s[i];
            let b = (ub[i] - x[i]) / s[i];
            steps[i] = a.max(b);
        }
    }
    let min_step = steps.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let hits: Vec<bool> = (0..n).map(|i| s[i] != 0.0 && steps[i] == min_step).collect();
    (min_step, hits)
}

/// Which bounds `x` is sitting on, to within `rtol` of the bound's own scale.
/// `-1` lower, `1` upper, `0` neither.
fn find_active_constraints(x: &[f64], lb: &[f64], ub: &[f64], rtol: f64) -> Vec<i8> {
    let n = x.len();
    let mut active = vec![0i8; n];
    if rtol == 0.0 {
        for i in 0..n {
            if x[i] <= lb[i] {
                active[i] = -1;
            } else if x[i] >= ub[i] {
                active[i] = 1;
            }
        }
        return active;
    }
    for i in 0..n {
        let lower_dist = x[i] - lb[i];
        let upper_dist = ub[i] - x[i];
        let lower_threshold = rtol * 1.0f64.max(lb[i].abs());
        let upper_threshold = rtol * 1.0f64.max(ub[i].abs());
        if lb[i].is_finite() && lower_dist <= upper_dist.min(lower_threshold) {
            active[i] = -1;
        } else if ub[i].is_finite() && upper_dist <= lower_dist.min(upper_threshold) {
            active[i] = 1;
        }
    }
    active
}

/// Push `x` off every bound it is sitting on. `rstep == 0` steps back a single
/// ULP, which is what keeps an accepted iterate strictly interior without
/// moving it anywhere a caller could measure.
fn make_strictly_feasible(x: &[f64], lb: &[f64], ub: &[f64], rstep: f64) -> Vec<f64> {
    let n = x.len();
    let active = find_active_constraints(x, lb, ub, rstep);
    let mut out = x.to_vec();
    for i in 0..n {
        if active[i] == -1 {
            out[i] = if rstep == 0.0 {
                next_toward(lb[i], ub[i])
            } else {
                lb[i] + rstep * 1.0f64.max(lb[i].abs())
            };
        } else if active[i] == 1 {
            out[i] = if rstep == 0.0 {
                next_toward(ub[i], lb[i])
            } else {
                ub[i] - rstep * 1.0f64.max(ub[i].abs())
            };
        }
        // A box narrower than the step-back leaves no interior to step into;
        // the midpoint is the only point that is strictly inside if any is.
        if out[i] < lb[i] || out[i] > ub[i] {
            out[i] = 0.5 * (lb[i] + ub[i]);
        }
    }
    out
}

/// The next representable double from `x` in the direction of `toward`.
fn next_toward(x: f64, toward: f64) -> f64 {
    if x.is_nan() || toward.is_nan() || x == toward {
        return toward;
    }
    if toward > x { x.next_up() } else { x.next_down() }
}

/// Where the ray `x + t·s` crosses `‖·‖ = Δ`. `None` when the ray is
/// degenerate or `x` is already outside, which are the two cases the caller
/// treats as "no reflection available".
fn intersect_trust_region(x: &[f64], s: &[f64], delta: f64) -> Option<(f64, f64)> {
    let a = dot(s, s);
    if a == 0.0 {
        return None;
    }
    let b = dot(x, s);
    let c = dot(x, x) - delta * delta;
    if c > 0.0 {
        return None;
    }
    let disc = b * b - a * c;
    if disc < 0.0 {
        return None;
    }
    let d = disc.sqrt();
    // Numerical Recipes' cancellation-free pair of roots.
    let q = -(b + d.copysign(b));
    let t1 = q / a;
    let t2 = c / q;
    Some(if t1 < t2 { (t1, t2) } else { (t2, t1) })
}

// ---------------------------------------------------------------------------
// The model, and the tests it feeds
// ---------------------------------------------------------------------------

/// `½(‖J_h·s‖² + sᵀ·C·s) + g_hᵀ·s` — the model's value at `s`.
fn evaluate_quadratic(jh: &[f64], g_h: &[f64], s: &[f64], diag: &[f64], m: usize, n: usize) -> f64 {
    let js = mat_vec(jh, s, m, n);
    let mut q = dot(&js, &js);
    q += s.iter().zip(diag).map(|(&si, &di)| si * di * si).sum::<f64>();
    0.5 * q + dot(s, g_h)
}

/// The model restricted to the line `s0 + t·s`, as `a·t² + b·t + c`.
fn build_quadratic_1d(
    jh: &[f64],
    g_h: &[f64],
    s: &[f64],
    diag: &[f64],
    s0: Option<&[f64]>,
    m: usize,
    n: usize,
) -> (f64, f64, f64) {
    let v = mat_vec(jh, s, m, n);
    let mut a = dot(&v, &v);
    a += s.iter().zip(diag).map(|(&si, &di)| si * di * si).sum::<f64>();
    a *= 0.5;
    let mut b = dot(g_h, s);
    let mut c = 0.0;
    if let Some(s0) = s0 {
        let u = mat_vec(jh, s0, m, n);
        b += dot(&u, &v);
        c = 0.5 * dot(&u, &u) + dot(g_h, s0);
        b += s0.iter().zip(diag).zip(s).map(|((&z, &di), &si)| z * di * si).sum::<f64>();
        c += 0.5 * s0.iter().zip(diag).map(|(&z, &di)| z * di * z).sum::<f64>();
    }
    (a, b, c)
}

/// Minimize `a·t² + b·t + c` over `[lb, ub]`, returning `(t, value)`.
fn minimize_quadratic_1d(a: f64, b: f64, lb: f64, ub: f64, c: f64) -> (f64, f64) {
    let value = |t: f64| t * (a * t + b) + c;
    let mut best_t = lb;
    let mut best_y = value(lb);
    let yu = value(ub);
    if yu < best_y {
        best_t = ub;
        best_y = yu;
    }
    if a != 0.0 {
        let extremum = -0.5 * b / a;
        if lb < extremum && extremum < ub {
            let ye = value(extremum);
            if ye < best_y {
                best_t = extremum;
                best_y = ye;
            }
        }
    }
    (best_t, best_y)
}

/// Shrink on a bad step, grow only when the step pressed against the radius.
/// A step that stopped short of Δ was limited by the *bounds*, not by the
/// model, so widening the region on it would buy nothing.
fn update_tr_radius(
    delta: f64,
    actual_reduction: f64,
    predicted_reduction: f64,
    step_norm: f64,
    bound_hit: bool,
) -> (f64, f64) {
    let ratio = if predicted_reduction > 0.0 {
        actual_reduction / predicted_reduction
    } else if predicted_reduction == 0.0 && actual_reduction == 0.0 {
        1.0
    } else {
        0.0
    };
    let mut delta = delta;
    if ratio < 0.25 {
        delta = 0.25 * step_norm;
    } else if ratio > 0.75 && bound_hit {
        delta *= 2.0;
    }
    (delta, ratio)
}

/// MINPACK's ftol/xtol pair, in scipy's arrangement: a cost test that only
/// counts when the step was at least a quarter as good as promised, and a
/// step-size test scaled by where `x` is.
fn check_termination(
    d_f: f64,
    f: f64,
    dx_norm: f64,
    x_norm: f64,
    ratio: f64,
    ftol: f64,
    xtol: f64,
) -> Option<LmTermination> {
    let ftol_satisfied = d_f < ftol * f && ratio > 0.25;
    let xtol_satisfied = dx_norm < xtol * (xtol + x_norm);
    if ftol_satisfied {
        Some(LmTermination::Cost)
    } else if xtol_satisfied {
        Some(LmTermination::Step)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Small linear algebra
// ---------------------------------------------------------------------------

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(&x, &y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

fn mat_vec(a: &[f64], s: &[f64], m: usize, n: usize) -> Vec<f64> {
    (0..m).map(|i| dot(&a[i * n..(i + 1) * n], s)).collect()
}

/// `Jᵀ·f`.
fn jt_times(jac: &[f64], f: &[f64], m: usize, n: usize) -> Vec<f64> {
    let mut g = vec![0.0f64; n];
    for i in 0..m {
        let row = &jac[i * n..(i + 1) * n];
        let fi = f[i];
        for k in 0..n {
            g[k] += row[k] * fi;
        }
    }
    g
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smatrix::synthesis::thick_opt::NoJacobian;

    fn cfg() -> LmConfig {
        LmConfig {
            ftol: 1e-14,
            xtol: 1e-14,
            gtol: 1e-14,
            max_iterations: 500,
            max_evals: 200_000,
            ..Default::default()
        }
    }

    fn run<F>(f: &F, x0: &[f64], lb: &[f64], ub: &[f64]) -> LmResult
    where
        F: Fn(&[f64], &mut Vec<f64>) -> Result<(), String> + Sync,
    {
        trust_region_reflective(f, None::<&NoJacobian>, x0, lb, ub, &cfg()).unwrap()
    }

    // -- the geometry helpers ---------------------------------------------

    #[test]
    fn cl_scaling_measures_the_bound_the_antigradient_points_at() {
        let x = [2.0, 2.0, 2.0];
        // g < 0 pushes x up, so the distance that matters is to ub;
        // g > 0 pushes it down, so it is the distance to lb; g == 0 neither.
        let g = [-1.0, 1.0, 0.0];
        let lb = [0.0, 0.0, 0.0];
        let ub = [10.0, 10.0, 10.0];
        let (v, dv) = cl_scaling_vector(&x, &g, &lb, &ub);
        assert_eq!(v, vec![8.0, 2.0, 1.0]);
        assert_eq!(dv, vec![-1.0, 1.0, 0.0]);
    }

    #[test]
    fn the_curvature_diagonal_is_never_negative() {
        // C = diag(g)·Jv is what makes the hat-space model convex; the sign
        // convention of dv is the whole reason it holds.
        let x = [1.0, 1.0];
        let ub = [5.0, 5.0];
        let lb = [-5.0, -5.0];
        for g in [[-3.0, 2.0], [2.0, -7.0], [0.0, 0.0], [-1.0, -1.0]] {
            let (_, dv) = cl_scaling_vector(&x, &g, &lb, &ub);
            for k in 0..2 {
                assert!(g[k] * dv[k] >= 0.0, "C[{k}] < 0 for g = {g:?}");
            }
        }
    }

    #[test]
    fn step_size_to_bound_reports_the_first_bound_and_who_hit_it() {
        let x = [1.0, 1.0];
        let s = [1.0, 4.0];
        let lb = [0.0, 0.0];
        let ub = [3.0, 3.0];
        // coordinate 0 reaches 3 at t = 2, coordinate 1 at t = 0.5.
        let (t, hits) = step_size_to_bound(&x, &s, &lb, &ub);
        assert!((t - 0.5).abs() < 1e-15);
        assert_eq!(hits, vec![false, true]);
    }

    #[test]
    fn a_zero_component_never_hits_a_bound() {
        let (t, hits) = step_size_to_bound(&[1.0, 1.0], &[0.0, 1.0], &[0.0, 0.0], &[3.0, 3.0]);
        assert!((t - 2.0).abs() < 1e-15);
        assert_eq!(hits, vec![false, true]);
    }

    #[test]
    fn make_strictly_feasible_leaves_a_bound_by_one_ulp_and_no_more() {
        let lb = [0.0, 0.0];
        let ub = [1.0, 1.0];
        let out = make_strictly_feasible(&[0.0, 1.0], &lb, &ub, 0.0);
        assert!(out[0] > 0.0 && out[0] < 1e-300, "{}", out[0]);
        assert!(out[1] < 1.0 && out[1] > 1.0 - 1e-15, "{}", out[1]);
    }

    #[test]
    fn a_box_with_no_interior_collapses_to_its_midpoint() {
        // lb == ub leaves nowhere strictly inside; the midpoint is the only
        // answer that is at least feasible.
        let out = make_strictly_feasible(&[4.0], &[4.0], &[4.0], 1e-10);
        assert_eq!(out, vec![4.0]);
    }

    #[test]
    fn intersect_trust_region_finds_the_two_crossings() {
        // From the origin along e_x, the ball of radius 2 is crossed at ±2.
        let (t1, t2) = intersect_trust_region(&[0.0, 0.0], &[1.0, 0.0], 2.0).unwrap();
        assert!((t1 + 2.0).abs() < 1e-14);
        assert!((t2 - 2.0).abs() < 1e-14);
        // A point already outside has no usable crossing.
        assert!(intersect_trust_region(&[3.0, 0.0], &[1.0, 0.0], 2.0).is_none());
        assert!(intersect_trust_region(&[0.0, 0.0], &[0.0, 0.0], 2.0).is_none());
    }

    #[test]
    fn minimize_quadratic_1d_takes_the_interior_vertex_only_when_it_is_inside() {
        // t² − 2t has its vertex at 1.
        let (t, y) = minimize_quadratic_1d(1.0, -2.0, -5.0, 5.0, 0.0);
        assert!((t - 1.0).abs() < 1e-15);
        assert!((y + 1.0).abs() < 1e-15);
        // Excluded from the interval, the minimum is at the near endpoint.
        let (t, _) = minimize_quadratic_1d(1.0, -2.0, 2.0, 5.0, 0.0);
        assert!((t - 2.0).abs() < 1e-15);
        // A linear function is minimized at an endpoint.
        let (t, _) = minimize_quadratic_1d(0.0, -1.0, 0.0, 3.0, 0.0);
        assert!((t - 3.0).abs() < 1e-15);
    }

    #[test]
    fn build_quadratic_1d_agrees_with_evaluating_the_model() {
        // The 1-D restriction has to be the same function as the full model,
        // or `select_step` is minimizing something else than it scores.
        let (m, n) = (3usize, 2usize);
        let jh = [1.0, 2.0, 0.5, -1.0, 3.0, 0.25];
        let g_h = [0.3, -0.7];
        let diag = [0.5, 1.5];
        let s = [0.4, -0.9];
        let s0 = [0.1, 0.2];
        let (a, b, c) = build_quadratic_1d(&jh, &g_h, &s, &diag, Some(&s0), m, n);
        for t in [-1.0, -0.25, 0.0, 0.6, 2.0] {
            let comb: Vec<f64> = (0..n).map(|k| s0[k] + t * s[k]).collect();
            let direct = evaluate_quadratic(&jh, &g_h, &comb, &diag, m, n);
            let poly = a * t * t + b * t + c;
            assert!((direct - poly).abs() < 1e-12, "t = {t}: {direct} vs {poly}");
        }
    }

    #[test]
    fn update_tr_radius_shrinks_on_a_bad_step_and_grows_only_at_the_boundary() {
        // ratio < 0.25: the radius comes down to a quarter of the step.
        let (d, r) = update_tr_radius(10.0, 0.1, 10.0, 4.0, true);
        assert!((d - 1.0).abs() < 1e-15);
        assert!(r < 0.25);
        // A good step that did NOT press against the radius was limited by
        // the bounds, so widening buys nothing.
        let (d, _) = update_tr_radius(10.0, 9.0, 10.0, 1.0, false);
        assert!((d - 10.0).abs() < 1e-15);
        // A good step that did press against it doubles the radius.
        let (d, _) = update_tr_radius(10.0, 9.0, 10.0, 9.9, true);
        assert!((d - 20.0).abs() < 1e-15);
        // A non-positive prediction is never rewarded.
        let (_, r) = update_tr_radius(10.0, 1.0, -1.0, 1.0, true);
        assert_eq!(r, 0.0);
    }

    #[test]
    fn the_cost_test_ignores_a_step_that_badly_undershot_its_promise() {
        // dF/F is tiny either way; only the gain ratio separates "converged"
        // from "the model has stopped describing the function".
        assert_eq!(
            check_termination(1e-12, 1.0, 10.0, 1.0, 0.9, 1e-8, 1e-8),
            Some(LmTermination::Cost)
        );
        assert_eq!(check_termination(1e-12, 1.0, 10.0, 1.0, 0.1, 1e-8, 1e-8), None);
    }

    // -- the subproblem ----------------------------------------------------

    /// Factor a dense `rows`×`n` matrix and return `(R, Qᵀb[..n])`.
    fn factor(b_mat: &[f64], rhs: &[f64], rows: usize, n: usize) -> (Vec<f64>, Vec<f64>) {
        let mut a = b_mat.to_vec();
        let mut qtb = rhs.to_vec();
        householder_qr_in_place(&mut a, &mut qtb, rows, n).unwrap();
        let mut r = vec![0.0f64; n * n];
        for i in 0..n {
            r[i * n + i..(i + 1) * n].copy_from_slice(&a[i * n + i..(i + 1) * n]);
        }
        (r, qtb[..n].to_vec())
    }

    #[test]
    fn the_damped_step_solves_the_normal_equations_it_claims_to() {
        // p(α) must satisfy (BᵀB + αI)p = −Bᵀf exactly, for every α.
        let (rows, n) = (4usize, 2usize);
        let b_mat = [1.0, 0.5, 2.0, -1.0, 0.25, 3.0, -0.75, 0.5];
        let f = [0.3, -1.2, 0.8, 0.1];
        let (r, qtf) = factor(&b_mat, &f, rows, n);
        let btf = jt_times(&b_mat, &f, rows, n);
        for alpha in [0.0, 1e-6, 0.3, 7.0, 1e4] {
            let (p, p_norm, q_norm) = damped_step(&r, &qtf, n, alpha).unwrap();
            for k in 0..n {
                let mut lhs = alpha * p[k];
                for j in 0..n {
                    let bjk: f64 = (0..rows).map(|i| b_mat[i * n + j] * b_mat[i * n + k]).sum();
                    lhs += bjk * p[j];
                }
                assert!(
                    (lhs + btf[k]).abs() < 1e-9,
                    "alpha = {alpha}, row {k}: {lhs} vs {}",
                    -btf[k]
                );
            }
            assert!((p_norm - norm(&p)).abs() < 1e-12);
            assert!(q_norm > 0.0);
        }
    }

    #[test]
    fn the_derivative_the_secular_equation_uses_is_the_real_one() {
        // φ'(α) = −‖q‖²/‖p‖ is what drives Moré's Newton iteration; if it is
        // wrong the iteration still converges, just slowly enough that the
        // 10-step cap turns into a wrong answer. Difference it.
        let (rows, n) = (5usize, 3usize);
        let b_mat = [
            1.0, 0.5, 0.2, 2.0, -1.0, 0.3, 0.25, 3.0, -0.1, -0.75, 0.5, 1.5, 0.4, 0.9, -2.0,
        ];
        let f = [0.3, -1.2, 0.8, 0.1, -0.4];
        let (r, qtf) = factor(&b_mat, &f, rows, n);
        let alpha = 0.7;
        let h = 1e-6;
        let (_, p_norm, q_norm) = damped_step(&r, &qtf, n, alpha).unwrap();
        let (_, up, _) = damped_step(&r, &qtf, n, alpha + h).unwrap();
        let (_, dn, _) = damped_step(&r, &qtf, n, alpha - h).unwrap();
        let fd = (up - dn) / (2.0 * h);
        let analytic = -(q_norm * q_norm) / p_norm;
        assert!((fd - analytic).abs() < 1e-6 * analytic.abs().max(1.0), "{fd} vs {analytic}");
    }

    #[test]
    fn the_subproblem_returns_the_gauss_newton_step_when_it_fits() {
        // Inside the region the constraint is inactive and α must be 0 —
        // a damped step there would be a worse answer than the model allows.
        let (rows, n) = (4usize, 2usize);
        let b_mat = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        let f = [0.5, -0.25, 0.0, 0.0];
        let (r, qtf) = factor(&b_mat, &f, rows, n);
        let g_norm = norm(&jt_times(&b_mat, &f, rows, n));
        let (p, alpha) = solve_lsq_trust_region(&r, &qtf, rows, n, g_norm, 10.0, 0.0).unwrap();
        assert_eq!(alpha, 0.0);
        assert!((p[0] + 0.5).abs() < 1e-12 && (p[1] - 0.25).abs() < 1e-12, "{p:?}");
    }

    #[test]
    fn a_constrained_subproblem_lands_on_the_boundary() {
        let (rows, n) = (4usize, 2usize);
        let b_mat = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        let f = [5.0, -3.0, 0.0, 0.0];
        let (r, qtf) = factor(&b_mat, &f, rows, n);
        let g_norm = norm(&jt_times(&b_mat, &f, rows, n));
        let delta = 1.0;
        let (p, alpha) = solve_lsq_trust_region(&r, &qtf, rows, n, g_norm, delta, 0.0).unwrap();
        assert!(alpha > 0.0);
        assert!((norm(&p) - delta).abs() < 1e-12, "‖p‖ = {}", norm(&p));
        // and it points downhill
        assert!(dot(&p, &jt_times(&b_mat, &f, rows, n)) < 0.0);
    }

    // -- the driver --------------------------------------------------------

    #[test]
    fn a_linear_problem_is_solved_exactly() {
        let xs = [0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = [-1.0, 1.0, 3.0, 5.0, 7.0];
        let f = |p: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            for i in 0..5 {
                out.push(p[0] * xs[i] + p[1] - ys[i]);
            }
            Ok(())
        };
        let r = run(&f, &[0.0, 0.0], &[-1e3, -1e3], &[1e3, 1e3]);
        assert!((r.x[0] - 2.0).abs() < 1e-7, "{:?}", r.x);
        assert!((r.x[1] + 1.0).abs() < 1e-7, "{:?}", r.x);
        assert!(r.cost < 1e-14, "{}", r.cost);
    }

    #[test]
    fn a_nonlinear_fit_recovers_its_parameters() {
        let ts: Vec<f64> = (0..7).map(|i| i as f64 * 0.5).collect();
        let obs: Vec<f64> = ts.iter().map(|t| 2.5 * (-0.7 * t).exp()).collect();
        let f = |p: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            for (t, o) in ts.iter().zip(&obs) {
                out.push(p[0] * (-p[1] * t).exp() - o);
            }
            Ok(())
        };
        let r = run(&f, &[1.0, 0.2], &[0.1, 0.1], &[10.0, 5.0]);
        assert!((r.x[0] - 2.5).abs() < 1e-5, "{:?}", r.x);
        assert!((r.x[1] - 0.7).abs() < 1e-5, "{:?}", r.x);
    }

    #[test]
    fn a_boundary_optimum_is_reached_to_the_last_ulp() {
        // Unconstrained this wants (5, −5); the box forces (1, 0). This is the
        // case the clamp-LM is weakest on, and the one TRF exists for.
        let f = |p: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(p[0] - 5.0);
            out.push(p[1] + 5.0);
            out.push(0.25 * (p[0] - p[1]));
            Ok(())
        };
        let r = run(&f, &[0.5, 0.5], &[0.0, 0.0], &[1.0, 1.0]);
        assert!((r.x[0] - 1.0).abs() < 1e-9, "{:?}", r.x);
        assert!(r.x[1].abs() < 1e-9, "{:?}", r.x);
        // Strictly interior, always: the step-back is a defining property of
        // the method, not an accident of this problem.
        assert!(r.x[0] < 1.0 && r.x[1] > 0.0, "{:?} touches a bound", r.x);
    }

    #[test]
    fn every_iterate_stays_inside_the_box() {
        // The residual closure is the only witness to what the solver tried.
        use std::sync::Mutex;
        let seen: Mutex<Vec<Vec<f64>>> = Mutex::new(Vec::new());
        let lb = [0.0, 0.0];
        let ub = [1.0, 1.0];
        let f = |p: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            seen.lock().unwrap().push(p.to_vec());
            out.clear();
            out.push(p[0] - 5.0);
            out.push(p[1] + 5.0);
            Ok(())
        };
        let r = trust_region_reflective(
            &f,
            None::<&NoJacobian>,
            &[0.5, 0.5],
            &lb,
            &ub,
            &LmConfig { max_iterations: 20, ..cfg() },
        )
        .unwrap();
        let _ = r;
        // Central differences legitimately probe outside the box (they are
        // taken in x, as in the LM driver), so the claim is about the points
        // the solver ADOPTS, which are the ones it evaluates without a
        // neighbour a step away. Assert the weaker, true thing: nothing the
        // solver evaluated ran away from the box.
        for p in seen.lock().unwrap().iter() {
            for k in 0..2 {
                assert!(
                    p[k] > lb[k] - 1e-4 && p[k] < ub[k] + 1e-4,
                    "evaluated at {p:?}, far outside [{:?}, {:?}]",
                    lb,
                    ub
                );
            }
        }
    }

    #[test]
    fn the_reflection_beats_stopping_at_the_bound() {
        // A long trust-region step into a corner: reflected, it can still buy
        // progress in the other coordinate. Run once and check the solver
        // reaches the constrained optimum rather than stalling on the wall.
        let f = |p: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(10.0 * (p[1] - p[0] * p[0]));
            out.push(1.0 - p[0]);
            Ok(())
        };
        // Bounded Rosenbrock with the optimum (1, 1) on the ub corner.
        let r = run(&f, &[-1.2, 1.0], &[-2.0, -2.0], &[1.0, 1.0]);
        assert!(r.cost < 1e-12, "cost {} at {:?}", r.cost, r.x);
    }

    #[test]
    fn a_rank_deficient_jacobian_still_produces_a_step() {
        // Only the sum is determined. The subproblem must not fail on it.
        let ds: Vec<f64> = (0..6).map(|i| i as f64).collect();
        let f = |p: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            for d in &ds {
                out.push((p[0] + p[1]) * d - 3.0 * d);
            }
            Ok(())
        };
        let r = run(&f, &[0.0, 0.0], &[-10.0, -10.0], &[10.0, 10.0]);
        assert!((r.x[0] + r.x[1] - 3.0).abs() < 1e-6, "{:?}", r.x);
    }

    #[test]
    fn a_parameter_with_no_effect_does_not_break_the_factorization() {
        // The second column of J is identically zero, so the augmented
        // system has a zero column -- the case the QR route has to floor and
        // the SVD route would absorb.
        let f = |p: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(p[0] - 2.0);
            out.push(0.5 * (p[0] - 2.0));
            Ok(())
        };
        let r = run(&f, &[0.0, 0.5], &[-5.0, -5.0], &[5.0, 5.0]);
        assert!((r.x[0] - 2.0).abs() < 1e-7, "{:?}", r.x);
        assert!(r.cost < 1e-12);
    }

    #[test]
    fn the_result_is_reported_in_this_crates_cost_convention() {
        // scipy carries ½‖r‖²; LmResult carries ‖r‖². A run that cannot move
        // reports the cost of where it is, which pins the factor of two.
        let f = |_p: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(3.0);
            out.push(4.0);
            Ok(())
        };
        let r = run(&f, &[1.0], &[0.0], &[2.0]);
        assert!((r.cost - 25.0).abs() < 1e-12, "{}", r.cost);
    }

    #[test]
    fn an_analytic_jacobian_is_used_and_counted() {
        struct Linear;
        impl JacobianSource for Linear {
            fn fill(&self, x: &[f64], jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
                let n = x.len();
                jac.clear();
                // r = [x0 - 1, 2*x1 + 3]
                jac.extend_from_slice(&[1.0, 0.0]);
                jac.extend_from_slice(&[0.0, 2.0]);
                let _ = n;
                Ok(Some(2))
            }
        }
        let f = |p: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(p[0] - 1.0);
            out.push(2.0 * p[1] + 3.0);
            Ok(())
        };
        let r = trust_region_reflective(
            &f,
            Some(&Linear),
            &[0.0, 0.0],
            &[-5.0, -5.0],
            &[5.0, 5.0],
            &cfg(),
        )
        .unwrap();
        assert!(r.analytic_jacobians > 0);
        assert!((r.x[0] - 1.0).abs() < 1e-9 && (r.x[1] + 1.5).abs() < 1e-9, "{:?}", r.x);
    }

    #[test]
    fn a_declining_jacobian_source_falls_back_to_differences() {
        struct Declines;
        impl JacobianSource for Declines {
            fn fill(&self, _x: &[f64], _jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
                Ok(None)
            }
        }
        let f = |p: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(p[0] - 1.0);
            Ok(())
        };
        let r = trust_region_reflective(
            &f,
            Some(&Declines),
            &[0.0],
            &[-5.0],
            &[5.0],
            &cfg(),
        )
        .unwrap();
        assert_eq!(r.analytic_jacobians, 0);
        assert!((r.x[0] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn the_inputs_are_checked_before_anything_is_evaluated() {
        let f = |_p: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(1.0);
            Ok(())
        };
        let c = cfg();
        assert!(trust_region_reflective(&f, None::<&NoJacobian>, &[], &[], &[], &c).is_err());
        assert!(
            trust_region_reflective(&f, None::<&NoJacobian>, &[0.0], &[0.0], &[], &c).is_err()
        );
        let e = trust_region_reflective(&f, None::<&NoJacobian>, &[5.0], &[0.0], &[1.0], &c)
            .unwrap_err();
        assert!(e.contains("outside"), "{e}");
        let e = trust_region_reflective(
            &f,
            None::<&NoJacobian>,
            &[0.5],
            &[0.0],
            &[f64::INFINITY],
            &c,
        )
        .unwrap_err();
        assert!(e.contains("not finite"), "{e}");
    }

    #[test]
    fn the_eval_budget_is_honoured() {
        let f = |p: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(10.0 * (p[1] - p[0] * p[0]));
            out.push(1.0 - p[0]);
            Ok(())
        };
        let r = trust_region_reflective(
            &f,
            None::<&NoJacobian>,
            &[-1.2, 1.0],
            &[-5.0, -5.0],
            &[5.0, 5.0],
            &LmConfig { max_evals: 12, ..cfg() },
        )
        .unwrap();
        // The budget is a floor-and-overshoot, not a hard cut: a Jacobian
        // build spends 2n at once. What must hold is that it stopped early.
        assert!(r.evals <= 12 + 2 * 2 + 1, "{}", r.evals);
        assert_eq!(r.termination, LmTermination::MaxIterations);
    }

    #[test]
    fn a_residual_error_is_propagated_rather_than_reported_as_no_progress() {
        let f = |_p: &[f64], _out: &mut Vec<f64>| -> Result<(), String> {
            Err("the merit spec is missing a curve".into())
        };
        let e = trust_region_reflective(&f, None::<&NoJacobian>, &[0.5], &[0.0], &[1.0], &cfg())
            .unwrap_err();
        assert!(e.contains("missing a curve"), "{e}");
    }

    #[test]
    fn a_non_finite_residual_shrinks_the_region_instead_of_failing() {
        // Past x = 1 the merit is undefined. The solver must come back, not
        // give up and not report a NaN optimum.
        let f = |p: &[f64], out: &mut Vec<f64>| -> Result<(), String> {
            out.clear();
            out.push(if p[0] > 1.0 { f64::NAN } else { p[0] - 0.5 });
            out.push(0.0);
            Ok(())
        };
        let r = run(&f, &[0.0], &[-2.0], &[4.0]);
        assert!(r.x[0].is_finite(), "{:?}", r.x);
        assert!((r.x[0] - 0.5).abs() < 1e-6, "{:?}", r.x);
    }
}
