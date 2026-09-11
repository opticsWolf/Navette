//! Landscape scanning, local minimization (Nelder–Mead) and field-profile
//! extraction on top of the coherent-block solver.
//!
//! The landscape maps trial effective index `n_eff` to |1/r|²; its minima
//! locate guided eigenmodes. [`char_func`] is the complex objective,
//! [`char_func_xy`] the real-valued adapter for the minimizer, and
//! [`reflection_coefficient_helper`] the shared forward solve.
use num_complex::Complex64;

// Reuse the shared low-level primitives from the crate
use crate::smatrix::coherent_block::solve_coherent_block_fields_inner;

// -----------------------------------------------------------------------------
// Helper: compute the complex reflection coefficient for a given n_eff
// -----------------------------------------------------------------------------
/// Complex reflection coefficient of the stack for trial index `n_eff` at
/// wavelength `lam` [nm] and polarization `pol` (0 = s, 1 = p).
///
/// Thin forward solve shared by the landscape scanner and the minimizer;
/// `inv_n` is the precomputed reciprocal-index stack (shared read-only).
#[inline]
pub fn reflection_coefficient_helper(
    n_stack: &[Complex64],
    inv_n: &[Complex64],
    thicknesses: &[f64],
    rough_types: &[i32],
    rough_vals: &[f64],
    lam: f64,
    n_eff: Complex64,
    pol: i32,
) -> Complex64 {
    // inv_n is computed once per wavelength by the caller and shared read-only
    // across all grid points / simplex evaluations (and across rayon threads).
    let (r_front, _, _, _, _, _, _, _) = solve_coherent_block_fields_inner(
        0,
        n_stack.len() - 1,
        n_stack,
        inv_n,
        thicknesses,
        rough_vals,
        rough_types,
        lam,
        n_eff,
        pol,
    );
    r_front
}

// -----------------------------------------------------------------------------
// Characteristic function: |1 / r(n_eff)|²
// -----------------------------------------------------------------------------
/// Characteristic function |1/r(n_eff)|² whose minima are the guided modes.
///
/// Returns 1e30 at exact zeros of |r| instead of dividing by zero.
#[inline]
pub fn char_func(
    n_stack: &[Complex64],
    inv_n: &[Complex64],
    thicknesses: &[f64],
    rough_types: &[i32],
    rough_vals: &[f64],
    lam: f64,
    n_eff: Complex64,
    pol: i32,
) -> f64 {
    let r = reflection_coefficient_helper(
        n_stack, inv_n, thicknesses, rough_types, rough_vals, lam, n_eff, pol,
    );
    let abs_r = r.norm();
    if !r.is_finite() || abs_r < 1e-15 {
        // A non-finite r means the forward solve overflowed (a huge trial
        // index, or a propagation exponent past the double range). Reporting
        // it as maximally bad keeps both the scanner and the minimizer honest;
        // returning NaN would make `find_minima`'s comparisons undefined.
        1e30
    } else {
        (1.0 / abs_r).powi(2)
    }
}

/// Physical box for a trial effective index: `3 x max |n_i|` over the stack.
///
/// Every guided mode of a passive stack has `min|n| <= |n_eff| <= max|n|`, so
/// 3x is generous even for leaky and substrate-side modes. It exists because
/// `|1/r|^2` has a spurious attractor at large `|n_eff|`: an unbounded
/// Nelder-Mead run on a seed with no mode nearby sprinted to
/// `n_eff = -2.0e8 + 8.4e6j` and reported `val = 1.3e-217` -- a confident
/// answer, twelve orders "better" than a real surface-plasmon pole, and
/// physically meaningless.
///
/// Note the bound when interpreting a rejection: on a stack whose highest
/// index is 3.5, trial indices past ~10.4 are refused outright.
#[inline]
pub fn n_eff_bound(n_stack: &[Complex64]) -> f64 {
    3.0 * n_stack.iter().map(|n| n.norm()).fold(0.0_f64, f64::max)
}

/// Real‑valued wrapper for minimisation (expects [Re, Im] slice)
#[allow(clippy::too_many_arguments)]
pub fn char_func_xy(
    xy: &[f64],
    n_stack: &[Complex64],
    inv_n: &[Complex64],
    thicknesses: &[f64],
    rough_types: &[i32],
    rough_vals: &[f64],
    lam: f64,
    pol: i32,
) -> f64 {
    let n_eff = Complex64::new(xy[0], xy[1]);
    // R3.3: the box lives here, in the minimizer's adapter, and NOT in
    // `char_func` itself. The landscape scanner also calls `char_func`, over a
    // range the user chose explicitly; painting 1e30 across the part of that
    // range outside the box would corrupt a diagnostic the caller asked for.
    // The minimizer is the one that needs walls -- it picks its own points.
    let bound = n_eff_bound(n_stack);
    if !xy[0].is_finite()
        || !xy[1].is_finite()
        || xy[0].abs() > bound
        || xy[1].abs() > bound
    {
        return 1e30;
    }
    char_func(n_stack, inv_n, thicknesses, rough_types, rough_vals, lam, n_eff, pol)
}

// -----------------------------------------------------------------------------
// Coarse landscape scan (parallel, GIL‑free)
