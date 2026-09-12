// SPDX-License-Identifier: LGPL-3.0-or-later
//! Pure-Rust shared primitives for the whole crate: constants, fast complex
//! kernels, the roughness form factor, the three Redheffer star products, and
//! the non-uniform spectral differentiation operator.
//!
//! This module is the SINGLE SOURCE OF TRUTH for everything that used to be
//! duplicated across func_0/func_1/func_2/func_4 and the needle operator.
//! It deliberately carries no pyo3 / numpy dependencies so it can be compiled
//! standalone (e.g. into the needle verification crate via `#[path]`).
//!
//! Numerical behaviour is byte-identical to the code that was moved here;
//! nothing was re-derived or "cleaned up".

use num_complex::Complex64;
use num_complex::ComplexFloat;

/// √3 (roughness/field-profile closed forms).
pub const SQRT3: f64 = 1.73205080757;
/// Denominator regularization floor shared by the field-amplitude star product.
pub const LOG_MIN: f64 = 1e-100;
/// p-pol cosθ magnitude guard (matches coherent_block).
pub const EPS_COS: f64 = 1e-12;
/// Machine epsilon used by the intensity/cross star regularization (inv := 0).
pub const DBL_EPS: f64 = 2.22e-16;
/// Speed of light in nm/fs — ω = 2π·c/λ with λ in nm gives rad/fs.
pub const C_NM_PER_FS: f64 = 299.792458;

/// Propagation phase of the differential-phase reference: the phase a plane
/// wave accumulates crossing an equivalent layer of incidence medium —
/// thickness `total_d`, real index `n_inc_re` — at incidence angle
/// `theta_inc_deg` (degrees, in the incidence medium, so no Snell step):
/// `passes · 2π · n · d · cosθ / λ`.
///
/// `wavelength` and `total_d` must share length units (nm in this crate).
/// `passes` is 1 for transmitted (single traversal, `PDts`/`PDtp`) and 2
/// for a reflection round trip. The incidence medium is assumed lossless;
/// pass `Re(n)` — any extinction only attenuates, it never phases.
/// Pure function (no solves); cost is one `cos` per call — hoist per angle
/// in hot loops via [`reference_wavenumber`].
///
/// Sign convention: `+kD` matches this crate's forward-propagation phase
/// (pinned by the `solver_propagation_sign_matches_reference` test — an
/// all-matched slab of D simulated by the solver has `arg(tf) = +kD`).
/// That is the conjugate of Macleod/`e^{+iωt}` textbooks; the crate is
/// self-consistent (all phase demands, needle `P_PHI` and GD/GDD share it),
/// so differential and absolute demands agree with each other exactly —
/// only textbook-imported target numbers need conjugating.
#[inline]
pub fn reference_phase(
    wavelength: f64,
    n_inc_re: f64,
    theta_inc_deg: f64,
    total_d: f64,
    passes: f64,
) -> f64 {
    passes * reference_wavenumber(wavelength, n_inc_re, theta_inc_deg) * total_d
}

/// Axial wavenumber part of [`reference_phase`]: `2π · n · cosθ / λ`
/// (radians per unit thickness; the caller scales by `passes · D`).
/// The needle gain correction for differential
/// demands is built from this: `dM/dD = Σ −2·kz·w·Δ/tol²` over folded points
/// (exact pre-fold data), a position-independent shift of `P(z)` that never
/// moves the needle site `argmax` — only the predicted-gain bookkeeping.
#[inline]
pub fn reference_wavenumber(wavelength: f64, n_inc_re: f64, theta_inc_deg: f64) -> f64 {
    2.0 * std::f64::consts::PI * n_inc_re * theta_inc_deg.to_radians().cos() / wavelength
}

#[inline(always)]
/// Construct a complex value without importing the trait soup at call sites.
pub fn cplx(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

/// Fast algebraic principal complex square root (`re >= 0` branch).
///
/// `num_complex`'s general-case `sqrt` routes through `to_polar`/`from_polar`,
/// i.e. an `atan2` plus a `sin`/`cos` — three transcendentals to take a root.
/// This computes the same principal value from `|z|` and two real `sqrt`s,
/// using the larger component first to avoid cancellation. Differs from the
/// polar result by at most ~1 ULP, and `cos θ` is taken every interface so the
/// saving compounds. Callers normalize the sign afterwards via
/// [`forward_branch`], so the branch convention matches the reference exactly.
///
/// `|z|` is taken via `hypot`, which is correctly rounded and avoids the
/// `a*a + b*b` intermediate (a couple of ULP of error plus over/underflow
/// risk). The earlier naive form was chosen to mirror the reference's
/// `fastmath` magnitude; with `fastmath` disabled this accurate form both
/// matches the stricter reference better and lowers the component-level
/// rounding that the `arg()` in `Delta` amplifies by `1/|cross|`. Cost is one
/// libm `hypot` instead of one `sqrt`-of-a-sum — no `atan2`, so the hot-loop
/// transcendental budget is unchanged in practice.
#[inline(always)]
pub fn csqrt_fast(z: Complex64) -> Complex64 {
    let a = z.re;
    let b = z.im;
    if a == 0.0 && b == 0.0 {
        return cplx(0.0, 0.0);
    }
    let m = a.hypot(b); // |z|, correctly-rounded (accurate, no a^2+b^2 intermediate)
    if a >= 0.0 {
        let re = ((m + a) * 0.5).sqrt();
        Complex64::new(re, b / (2.0 * re))
    } else {
        let im = ((m - a) * 0.5).sqrt().copysign(b);
        Complex64::new(b / (2.0 * im), im)
    }
}

/// Pick the forward branch of `cos theta`: the wave that does not grow in +z.
///
/// SINGLE SOURCE OF TRUTH for the branch cut. [`csqrt_fast`] returns the
/// principal root; this decides whether to keep it. Five call sites carried
/// this four-line test copied out — both `coherent_block` solvers, front
/// interface and interior, plus `needle_operator::cos_from_nsin` — and it is
/// subtle enough to deserve one home.
///
/// The rule is `Im(cos) >= 0`, bit-for-bit what those five copies did. For a
/// real incident index that is exactly the physical condition `Im(kz) >= 0`
/// with `kz = n*cos`: for real `n > 0` the two tests share a sign. Absorbing
/// *layers* are handled correctly and never even reach the flip — `nsin/n` has
/// negative imaginary part, so `1 - (nsin/n)^2` has positive imaginary part,
/// so `Im cos > 0` already.
///
/// `n` is taken and ignored. It is here because the obvious "improvement" is
/// to test `Im(n*cos)` instead, and the parameter is where the reason not to
/// lives.
///
/// **The rule is not extended to cover a complex incident index, and that is
/// not an oversight** (R3.4). With `n0` complex the transverse wavevector
/// `kx = k0*n0*sin(theta)` is complex, the wave is inhomogeneous, and "does
/// not grow in +z" stops being the right criterion at all: a transparent layer
/// under an absorbing ambient legitimately grows along z while decaying along
/// x. Which root is forward then depends on the inhomogeneity of the incident
/// wave, and a real angle of incidence does not specify it. The input is
/// under-determined, not merely awkward to normalize — so the ambiguity is
/// removed at the source instead of a root being picked here: [`sanitize_incident_index`]
/// zeroes `Im(n0)` and reports it, and the stack is then solved with a
/// transparent ambient of index `Re(n0)`, for which this simple rule is
/// exactly right. (0.6.21 refused the stack outright; 0.6.26 corrected and
/// warned at the `ScatterMatrix` door; 0.6.27 applies the same rule at
/// `DesignStack::from_design`, which is where the synthesis surface builds
/// its stacks — that door had been open the whole time.)
///
/// Nothing reaches this function with a complex ambient unless the caller went
/// around those doors, which the native `Solver` still permits by design (R3.1).
///
/// Measured on the permissive native path, which is what such a caller gets:
/// for a complex ambient the test is not just wrong, it is *undecidable*. `nsin/n0`
/// should be exactly real; what comes back carries a rounding residue of order
/// 1e-31 whose sign depends on how `1/n0` rounded, and flipping on it inverts
/// `Re cos` from +0.985 to -0.985. A 2-layer stack at 10 degrees gave `Rs`
/// alternating between 0.0024 and 417 as `k_ambient` moved 1e-16 -> 1e-2, with
/// no monotonicity, because the flip tracked rounding rather than physics.
#[inline(always)]
pub fn forward_branch(cos_theta: Complex64, _n: Complex64) -> Complex64 {
    if cos_theta.im < 0.0 {
        -cos_theta
    } else {
        cos_theta
    }
}

/// The shared explanation carried by every "ambient absorption dropped"
/// warning. The Python `ScatterMatrix` door emits the same sentences from its
/// own copy of the rule; `test_both_doors_explain_it_the_same_way` pins the
/// two together so they cannot drift apart.
pub const AMBIENT_DROP_EXPLANATION: &str =
    "Its absorption has been dropped -- the stack is solved with a transparent \
     ambient of index Re(n[0]). Reflectance is not defined against an absorbing \
     ambient: R = |r|^2 stops being an energy ratio (R + T climbs past 1), and at \
     oblique incidence the transverse wavevector goes complex, so a real angle of \
     incidence no longer identifies one incident wave. The results are therefore \
     the transparent-ambient values; attenuation along the path through the \
     ambient is not included. Put the absorbing medium on the substrate side \
     (fully supported) if you need it carried.";

/// Strip absorption from an incident-medium index vector.
///
/// Returns `None` when `nk` is already transparent -- one scan, no allocation,
/// which is the entire cost on the common path. Otherwise returns the
/// corrected vector (`Im(n) == 0` throughout) plus a warning naming what was
/// dropped, for the caller to surface however it surfaces warnings.
///
/// This is the disarming half of [`forward_branch`]. That branch rule is
/// exactly right when the ambient is real, and a tripwire when it is not:
/// swept over ambients 1.0 / 1.52 / 2.35 / 4.0, every angle from 0 to 89.9
/// degrees in 0.1 steps, against transparent, weakly absorbing, strongly
/// absorbing and metal-like layers -- 32400 cases, total internal reflection
/// included -- a real ambient never flips it once. The same sweep with `k` on
/// the ambient flips 10680 of 14256, deciding on an `Im(cos)` of order 1e-20.
/// Calling this at assembly time is what makes the rule downstream correct,
/// rather than correct-by-assumption.
///
/// Dropping `Im(n0)` is a correction, not an approximation: with it gone the
/// transverse wavevector is real again, every layer takes the standard branch,
/// `R = |r|^2` is a true energy ratio and `R + T + A == 1` holds. What is lost
/// is attenuation along the path *through* the ambient, a geometry-dependent
/// factor that a semi-infinite half-space does not define in the first place.
/// See R3.4, and `navette.smatrix._sanitize_incident_medium` for the Python
/// door, which applies the same rule to `ScatterMatrix`.
///
/// COLD PATH BY CONSTRUCTION. The one engine call site is
/// `DesignStack::from_design`, which runs once per stack assembly. It must
/// never be called from `solver_arrays` or anything below it -- those run per
/// merit evaluation, thousands of times per design run.
pub fn sanitize_incident_index(nk: &[Complex64]) -> Option<(Vec<Complex64>, String)> {
    let first = nk.iter().position(|z| z.im != 0.0)?;
    let n_bad = nk[first..].iter().filter(|z| z.im != 0.0).count();
    let k_max = nk.iter().map(|z| z.im.abs()).fold(0.0_f64, f64::max);
    let fixed: Vec<Complex64> = nk.iter().map(|z| Complex64::new(z.re, 0.0)).collect();
    let msg = format!(
        "the incident medium (layer 0) is absorbing at {} of {} wavelengths \
         (first at index {}: {}; largest |Im(n)| = {:e}). {}",
        n_bad,
        nk.len(),
        first,
        nk[first],
        k_max,
        AMBIENT_DROP_EXPLANATION,
    );
    Some((fixed, msg))
}

/// Fast complex exponential. Same formula as `num_complex`
/// (`e^re · (cos im, sin im)`) but via a single `sin_cos`, which shares the
/// argument reduction between sine and cosine.
#[inline(always)]
pub fn cexp_fast(z: Complex64) -> Complex64 {
    let e = z.re.exp();
    let (s, c) = z.im.sin_cos();
    Complex64::new(e * c, e * s)
}

/// Névot-Croce (roughness type 5) amplitude factors `(reflection, transmission)`.
///
/// SINGLE SOURCE OF TRUTH for the type-5 branch. Every interface builder —
/// `coherent_block::solve_pol_specialized`,
/// `coherent_block::solve_coherent_block_fields_dual`,
/// `solver::field_prof` and `needle_operator::interface_matrix` — must call
/// this rather than open-coding the exponentials, so a fix cannot land in
/// some sites and miss others (it did: see R1.1).
///
/// * reflection  `f  = exp(−2·kz1·kz2·σ²)`   — correlated Debye-Waller damping
/// * transmission `ga = exp(+((kz1−kz2)·σ)²/2)` — ≥ 1, → 1 at low contrast
///
/// The two are conjugate: to first order in σ² the reflected loss and the
/// transmitted gain cancel exactly, so `R + T = 1` for a lossless interface.
///
/// SPECULAR ONLY. This conserves energy among the coherent beams (graded
/// interface picture) and models no diffuse scatter. It is only
/// *perturbatively* unitary — valid for `|kz·σ| ≪ 1`. Outside that regime the
/// transmission factor overshoots without bound (e.g. n=1→4.28 at σ=20 nm,
/// λ=550 nm gives T ≈ 1.08 and a negative residual absorptance). See
/// `RoughnessType` docs and `docs/plans/scatter_loss_plan.md`.
#[inline(always)]
pub fn nevot_croce_factors(kz1: Complex64, kz2: Complex64, sigma: f64) -> (Complex64, Complex64) {
    let f = (-2.0 * kz1 * kz2 * sigma * sigma).exp();
    let d = (kz1 - kz2) * sigma;
    let ga = (d * d * 0.5).exp();
    (f, ga)
}

/// Roughness form factor W(q). No Python/PyO3 overhead so it can be called
/// millions of times from hot loops and inlined by the compiler.
#[inline(always)]
pub fn w_function_inner(q: Complex64, rough_type: i32) -> Complex64 {
    match rough_type {
        0 => Complex64::new(1.0, 0.0),
        1 => {
            let val = q * SQRT3;
            if val.norm() < 1e-9 {
                Complex64::new(1.0, 0.0)
            } else {
                val.sin() / val
            }
        }
        2 => q.cos(),
        3 => {
            let denom = Complex64::new(1.0, 0.0) + (q * q) * 0.5;
            Complex64::new(1.0, 0.0) / denom
        }
        4 => (-(q * q) * 0.5).exp(),
        _ => Complex64::new(1.0, 0.0),
    }
}

/// Complex (field-amplitude) Redheffer star product. Pure Rust, inlined.
#[inline(always)]
pub fn redheffer_product_complex_field_inner(
    r_a_front: Complex64,
    t_a_back: Complex64,
    t_a_fwd: Complex64,
    r_a_back: Complex64,
    r_b_front: Complex64,
    t_b_back: Complex64,
    t_b_fwd: Complex64,
    r_b_back: Complex64,
) -> (Complex64, Complex64, Complex64, Complex64) {
    let mut denom = Complex64::new(1.0, 0.0) - r_a_back * r_b_front;
    if denom.abs() < LOG_MIN {
        // Phase-preserving regularization (matches the reference implementation).
        let phase = denom / (denom.abs() + 1e-300);
        denom = Complex64::new(LOG_MIN, 0.0) * phase + 1e-300;
    }
    let inv_denom = denom.recip();

    let s_r_front = r_a_front + t_a_back * r_b_front * t_a_fwd * inv_denom;
    let s_t_back = t_a_back * t_b_back * inv_denom;
    let s_t_fwd = t_b_fwd * t_a_fwd * inv_denom;
    let s_r_back = r_b_back + t_b_fwd * r_a_back * t_b_back * inv_denom;

    (s_r_front, s_t_back, s_t_fwd, s_r_back)
}

/// Real-valued (intensity) Redheffer star product. Pure Rust, inlined.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub fn redheffer_product_real_inner(
    ra_rf: f64,
    ra_tb: f64,
    ra_tf: f64,
    ra_rb: f64,
    rb_rf: f64,
    rb_tb: f64,
    rb_tf: f64,
    rb_rb: f64,
) -> (f64, f64, f64, f64) {
    let denom = 1.0 - ra_rb * rb_rf;
    let inv_denom = if denom.abs() < DBL_EPS { 0.0 } else { 1.0 / denom };

    let rf = ra_rf + ra_tb * rb_rf * ra_tf * inv_denom;
    let tb = ra_tb * rb_tb * inv_denom;
    let tf = rb_tf * ra_tf * inv_denom;
    let rb = rb_rb + rb_tf * ra_rb * rb_tb * inv_denom;

    (rf, tb, tf, rb)
}

/// Complex (coherency) Redheffer star product over the p-s cross-amplitudes
/// C = (p-field)·conj(s-field). Structurally identical to
/// `redheffer_product_real_inner`; the denominator `1 - C_Ab·C_Bf` is the
/// n=m term of the incoherent multiple-reflection geometric series in the
/// cross channel (different bounce orders are mutually incoherent). On the
/// diagonal (C = |field|²) this collapses exactly onto the real product, which
/// is why R/T are unaffected by the coherency mode.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub fn redheffer_product_cross_inner(
    a_cf: Complex64,
    a_db: Complex64,
    a_df: Complex64,
    a_cb: Complex64,
    b_cf: Complex64,
    b_db: Complex64,
    b_df: Complex64,
    b_cb: Complex64,
) -> (Complex64, Complex64, Complex64, Complex64) {
    let denom = Complex64::new(1.0, 0.0) - a_cb * b_cf;
    let inv = if denom.norm() < DBL_EPS {
        Complex64::new(0.0, 0.0)
    } else {
        denom.recip()
    };

    let cf = a_cf + a_db * b_cf * a_df * inv;
    let db = a_db * b_db * inv;
    let df = b_df * a_df * inv;
    let cb = b_cb + b_df * a_cb * b_db * inv;

    (cf, db, df, cb)
}

/// Second-order non-uniform central-difference derivative dy/dx
/// (numpy.gradient semantics; one-sided at endpoints). Single source of truth
/// for the dispersion post-passes of both the core engine and the needle
/// operator's spectral sensitivity chain.
pub fn grad_nonuniform(y: &[f64], x: &[f64]) -> Vec<f64> {
    let n = y.len();
    let mut d = vec![0.0; n];
    if n < 2 {
        return d;
    }
    d[0] = (y[1] - y[0]) / (x[1] - x[0]);
    d[n - 1] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
    for i in 1..n - 1 {
        let hd = x[i] - x[i - 1];
        let hs = x[i + 1] - x[i];
        d[i] = (-hs / (hd * (hd + hs))) * y[i - 1]
            + ((hs - hd) / (hd * hs)) * y[i]
            + (hd / (hs * (hd + hs))) * y[i + 1];
    }
    d
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_phase_hand_calc() {
        // λ = 500 nm, n = 1, θ = 0°, D = 100 nm, 1 pass:
        // 2π·1·100·1/500 = 2π/5 = 1.2566370614...
        let p = reference_phase(500.0, 1.0, 0.0, 100.0, 1.0);
        assert!((p - 2.0 * std::f64::consts::PI / 5.0).abs() < 1e-12, "p={p}");
        // Oblique: θ = 60° halves the axial projection.
        let po = reference_phase(500.0, 1.0, 60.0, 100.0, 1.0);
        assert!((po - p / 2.0).abs() < 1e-12, "po={po}");
        // Two passes double; zero thickness kills the reference.
        assert!((reference_phase(500.0, 1.0, 0.0, 100.0, 2.0) - 2.0 * p).abs() < 1e-12);
        assert_eq!(reference_phase(500.0, 1.5, 10.0, 0.0, 1.0), 0.0);
        // Wavenumber × D × passes reconstructs the phase.
        let kz = reference_wavenumber(500.0, 1.5, 10.0);
        assert!((1.0 * kz * 100.0 - reference_phase(500.0, 1.5, 10.0, 100.0, 1.0)).abs() < 1e-15);
    }

    /// The exact expression the five call sites carried before R3.4
    /// consolidated them. Kept here as the oracle: `forward_branch` must be
    /// this, bit for bit, forever.
    fn old_inline_rule(c: Complex64) -> Complex64 {
        if c.im < 0.0 {
            -c
        } else {
            c
        }
    }

    /// `cos` from `nsin` the way the sites compute it, without the branch.
    fn raw_cos(nsin: Complex64, n: Complex64) -> Complex64 {
        let r0 = nsin / n;
        csqrt_fast(Complex64::new(1.0, 0.0) - r0 * r0)
    }

    #[test]
    fn forward_branch_is_the_old_inline_rule_bit_for_bit() {
        // A deterministic sweep over the quadrants, including the signed
        // zeros and the exactly-real and exactly-imaginary axes, since the
        // rule branches on `im < 0.0` and -0.0 is NOT < 0.0.
        let parts = [-3.25_f64, -1.0, -1e-300, -0.0, 0.0, 1e-300, 1.0, 3.25];
        for &re in &parts {
            for &im in &parts {
                let c = Complex64::new(re, im);
                for &n in &[
                    Complex64::new(1.0, 0.0),
                    Complex64::new(1.52, 0.0),
                    Complex64::new(2.35, 0.8),
                    Complex64::new(-1.0, -1.0),
                ] {
                    let got = forward_branch(c, n);
                    let want = old_inline_rule(c);
                    assert_eq!(got.re.to_bits(), want.re.to_bits(), "re c={c} n={n}");
                    assert_eq!(got.im.to_bits(), want.im.to_bits(), "im c={c} n={n}");
                }
            }
        }
    }

    #[test]
    fn forward_branch_ignores_the_index_argument() {
        // `n` is taken and ignored on purpose (see the doc comment). Pin it,
        // so that "improving" the rule to test `Im(n*cos)` cannot slip in
        // unnoticed: for n = 2.35 + 0.8i the two rules genuinely disagree.
        let c = Complex64::new(0.7, 0.1);
        let n_real = Complex64::new(1.52, 0.0);
        let n_cplx = Complex64::new(2.35, 0.8);
        assert_eq!(forward_branch(c, n_real), forward_branch(c, n_cplx));
        // And there really is a case where the two rules disagree, so the
        // assertion above is not vacuous: for cos = -1 + 0.1i the current rule
        // keeps it (Im cos > 0) while an Im(n*cos) rule would flip it
        // (2.35*0.1 + 0.8*(-1) < 0).
        let c2 = Complex64::new(-1.0, 0.1);
        assert!(c2.im > 0.0);
        assert!((n_cplx * c2).im < 0.0);
        assert_eq!(forward_branch(c2, n_cplx), c2);
    }

    #[test]
    fn forward_branch_is_identity_for_a_real_ambient() {
        // Propagating: 30 deg from n0 = 1.0 into n = 1.5. `cos` is real and
        // positive, nothing to flip.
        let nsin = Complex64::new(30.0_f64.to_radians().sin(), 0.0);
        let c = raw_cos(nsin, Complex64::new(1.5, 0.0));
        assert_eq!(c.im, 0.0);
        assert!(c.re > 0.0);
        assert_eq!(forward_branch(c, Complex64::new(1.5, 0.0)), c);

        // Evanescent: total internal reflection, 60 deg from n0 = 1.5 into
        // n = 1.0. `cos` is purely imaginary and the principal root already
        // has the decaying sign.
        let nsin_tir = Complex64::new(1.5 * 60.0_f64.to_radians().sin(), 0.0);
        let c_tir = raw_cos(nsin_tir, Complex64::new(1.0, 0.0));
        assert!(c_tir.re.abs() < 1e-15, "c_tir={c_tir}");
        assert!(c_tir.im > 0.0, "c_tir={c_tir}");
        assert_eq!(forward_branch(c_tir, Complex64::new(1.0, 0.0)), c_tir);

        // Exactly critical: nsin == n, so `cos` is exactly zero. +0.0 is not
        // < 0.0, so the branch leaves it alone rather than producing -0.0.
        let c_crit = raw_cos(Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0));
        assert_eq!(c_crit, Complex64::new(0.0, 0.0));
        let out = forward_branch(c_crit, Complex64::new(1.0, 0.0));
        assert_eq!(out.re.to_bits(), 0.0_f64.to_bits());
        assert_eq!(out.im.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn sanitize_incident_index_leaves_a_real_ambient_alone() {
        // The common path: no allocation, no warning, nothing to say.
        assert!(sanitize_incident_index(&[Complex64::new(1.0, 0.0); 4]).is_none());
        assert!(sanitize_incident_index(&[]).is_none());
        // -0.0 is not absorption. It compares equal to 0.0 and it does not
        // trip `forward_branch` either (`-0.0 < 0.0` is false), so warning
        // about it would be noise. Pinned because the obvious "tidy up the
        // sign" edit would start warning on every such grid.
        assert!(sanitize_incident_index(&[Complex64::new(1.0, -0.0)]).is_none());
    }

    #[test]
    fn sanitize_incident_index_flattens_and_says_what_it_dropped() {
        let nk = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.05),
            Complex64::new(1.5, -0.2),
        ];
        let (fixed, msg) = sanitize_incident_index(&nk).expect("should have fired");
        assert_eq!(
            fixed,
            vec![
                Complex64::new(1.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(1.5, 0.0),
            ]
        );
        // Counts, not "the row is bad": a single stray wavelength in a
        // 200-point grid reads very differently from all 200.
        assert!(msg.contains("2 of 3 wavelengths"), "{msg}");
        assert!(msg.contains("first at index 1"), "{msg}");
        assert!(msg.contains("2e-1"), "{msg}"); // largest |Im(n)| = 0.2
        assert!(msg.contains(AMBIENT_DROP_EXPLANATION), "{msg}");
        assert!(msg.is_ascii(), "the message must survive a cp1252 console");
        // Idempotent: the corrected ambient is the one the branch rule below
        // is provably safe for, so a second pass has nothing to do.
        assert!(sanitize_incident_index(&fixed).is_none());
    }

    #[test]
    fn a_sanitized_ambient_never_reaches_the_flip() {
        // The pairing this whole item rests on, measured rather than argued:
        // with the ambient real, `forward_branch` is the identity across the
        // full angle range -- total internal reflection included, which the
        // test below does not cover (it only ever goes from n=1 into n=2.35,
        // so r0 < 1 always and the evanescent regime is never entered).
        let mut checked = 0usize;
        for &n0 in &[1.52_f64, 1.0, 2.35, 4.0] {
            for i in 0..=899 {
                let deg = i as f64 * 0.1;
                let nsin = Complex64::new(n0 * deg.to_radians().sin(), 0.0);
                for &(nr, k) in &[
                    (1.0_f64, 0.0_f64),
                    (1.0, 1e-12),
                    (1.0, 0.05),
                    (1.38, 0.0),
                    (1.38, 0.5),
                    (2.35, 0.0),
                    (2.35, 3.0),
                    (0.15, 3.5),
                ] {
                    let n = Complex64::new(nr, k);
                    let c = raw_cos(nsin, n);
                    checked += 1;
                    assert_eq!(forward_branch(c, n), c, "n0={n0} deg={deg} n={n} c={c}");
                }
            }
        }
        assert_eq!(checked, 28800);
    }

    #[test]
    fn absorbing_layers_under_a_real_ambient_never_reach_the_flip() {
        // The claim in the doc comment: with a real ambient, `1 - (nsin/n)^2`
        // has positive imaginary part for any absorbing layer, so the
        // principal square root is already forward and the flip never fires.
        // Checked across the angle range and three decades of k.
        //
        // Note what this does NOT establish: the ambient here is n = 1 and
        // the layer n = 2.35, so `r0 < 1` and the evanescent regime is never
        // entered. `a_sanitized_ambient_never_reaches_the_flip` above covers
        // that, and it is the one that speaks for the supported path.
        for &deg in &[0.0_f64, 10.0, 30.0, 45.0, 60.0, 75.0, 89.0] {
            let nsin = Complex64::new(1.0 * deg.to_radians().sin(), 0.0);
            for &k in &[1e-6_f64, 1e-3, 0.05, 0.5, 3.0] {
                let n = Complex64::new(2.35, k);
                let c = raw_cos(nsin, n);
                assert!(c.im >= 0.0, "deg={deg} k={k} c={c}");
                assert_eq!(forward_branch(c, n), c, "deg={deg} k={k}");
            }
        }
    }
}
