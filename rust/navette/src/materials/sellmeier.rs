// SPDX-License-Identifier: LGPL-3.0-or-later
//! Sellmeier dispersion model.
//!
//! n²(λ) = 1 + Σᵢ Bᵢ·λ² / (λ² − Cᵢ)   with λ in µm.  The third term is only
//! included when B3 ≠ 0 (matching the Python `np.where(B3 != 0.0, …, 0.0)`).
//! k = 0 (or Urbach tail for the `_urbach` variant).
//!
//! **Domain.** A Sellmeier fit is only a function between its poles. At
//! λ² = Cᵢ the i-th term is infinite; on the short-wavelength side of the
//! nearest pole n² turns negative and `sqrt` answers NaN. Neither is a
//! numerical accident to be smoothed over — it is the model being evaluated
//! outside the range it was fitted on, and a published coefficient set carries
//! that range in the paper it came from, not in the numbers. The kernels here
//! return what the arithmetic gives (inf / NaN, faithful to the Python
//! reference); [`sellmeier_domain_check`] turns that into an error naming the
//! wavelength, and both callers that face a user run it. See R6.2.

use ndarray::{Array1, ArrayView1};
use num_complex::Complex64;

use crate::materials::common::{map_nk, urbach_k};
use crate::materials::units::wl_um2;

/// Real refractive index from the (up to) three-term Sellmeier equation.
#[allow(clippy::too_many_arguments)]
#[inline]
pub fn sellmeier_n(
    l2: f64, // λ² in µm²
    b1: f64,
    c1: f64,
    b2: f64,
    c2: f64,
    b3: f64,
    c3: f64,
) -> f64 {
    let term1 = b1 * l2 / (l2 - c1);
    let term2 = b2 * l2 / (l2 - c2);
    // Faithful to the Python branch: contributes only when B3 is nonzero.
    let term3 = if b3 != 0.0 { b3 * l2 / (l2 - c3) } else { 0.0 };
    let n_sq = 1.0 + term1 + term2 + term3;
    n_sq.sqrt()
}

/// Reject a wavelength grid that leaves the Sellmeier model's domain.
///
/// Runs on the *result*, so the common case costs one pass looking for a
/// non-finite real part and nothing else; the coefficient arithmetic that
/// works out which pole was hit only happens once something is already wrong.
///
/// Checks `n` only. The Urbach tail is a separate model on `k` with its own
/// parameters, and an overflow there is not a Sellmeier pole.
///
/// The message is ASCII on purpose: it reaches Python as a `ValueError`, and a
/// cp1252 console cannot print what it cannot encode.
pub fn sellmeier_domain_check(
    wavelength_nm: ArrayView1<f64>,
    nk: &Array1<Complex64>,
    b1: f64,
    c1: f64,
    b2: f64,
    c2: f64,
    b3: f64,
    c3: f64,
) -> Result<(), String> {
    let Some(i) = nk.iter().position(|v| !v.re.is_finite()) else {
        return Ok(());
    };
    let w = wavelength_nm[i];
    let l2 = wl_um2(w);
    // Which term(s) the grid walked into. B3 = 0 drops the third term entirely,
    // so its C3 is not a pole at all (see `sellmeier_n`).
    let terms: [(usize, f64, f64); 3] = [(1, b1, c1), (2, b2, c2), (3, b3, c3)];
    let mut poles: Vec<String> = Vec::new();
    for (idx, b, c) in terms {
        if idx == 3 && b == 0.0 {
            continue;
        }
        if c > 0.0 {
            poles.push(format!(
                "C{idx}={c} um^2 (resonance at {:.4} nm)",
                c.sqrt() * 1000.0
            ));
        }
    }
    let at_pole = terms
        .iter()
        .any(|&(idx, b, c)| !(idx == 3 && b == 0.0) && l2 == c);
    let cause = if at_pole {
        "sits exactly on a pole (lambda^2 == C), so n^2 is infinite"
    } else if nk[i].re.is_nan() {
        "is on the short side of a resonance, where n^2 < 0 and sqrt gives NaN"
    } else {
        "makes n^2 non-finite"
    };
    let where_ = if poles.is_empty() {
        "none of its C coefficients are positive, so it has no real resonance".to_string()
    } else {
        format!("its poles are {}", poles.join(", "))
    };
    Err(format!(
        "Sellmeier: wavelength {w:.6} nm (lambda^2 = {l2:.6} um^2) {cause}, and \
     {where_}. Restrict the grid to one side of the nearest resonance, or use \
     coefficients fitted for this range -- a non-finite index propagates as NaN \
     through every layer of the stack with nothing left to say where it started."
    ))
}

/// Sellmeier complex refractive index (k = 0).
#[allow(clippy::too_many_arguments)]
pub fn sellmeier_nk(
    wavelength_nm: ArrayView1<f64>,
    b1: f64,
    c1: f64,
    b2: f64,
    c2: f64,
    b3: f64,
    c3: f64,
) -> Array1<Complex64> {
    map_nk(wavelength_nm, |w| {
        Complex64::new(sellmeier_n(wl_um2(w), b1, c1, b2, c2, b3, c3), 0.0)
    })
}

/// Sellmeier dispersion with an Urbach absorption tail on k.
#[allow(clippy::too_many_arguments)]
pub fn sellmeier_urbach_nk(
    wavelength_nm: ArrayView1<f64>,
    b1: f64,
    c1: f64,
    b2: f64,
    c2: f64,
    b3: f64,
    c3: f64,
    alpha0: f64,
    eu: f64,
    lambda_g: f64,
) -> Array1<Complex64> {
    map_nk(wavelength_nm, |w| {
        let n = sellmeier_n(wl_um2(w), b1, c1, b2, c2, b3, c3);
        let k = urbach_k(w, alpha0, eu, lambda_g);
        Complex64::new(n, k)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    // BK7 (Schott), the coefficient set the parity tests use.
    const BK7: [f64; 6] = [
        1.03961212,
        0.00600069867,
        0.231792344,
        0.0200179144,
        1.01046945,
        103.560653,
    ];

    fn check(wl: &Array1<f64>) -> Result<(), String> {
        let out = sellmeier_nk(wl.view(), BK7[0], BK7[1], BK7[2], BK7[3], BK7[4], BK7[5]);
        sellmeier_domain_check(
            wl.view(),
            &out,
            BK7[0],
            BK7[1],
            BK7[2],
            BK7[3],
            BK7[4],
            BK7[5],
        )
    }

    #[test]
    fn a_grid_inside_the_fit_is_accepted() {
        // BK7's shortest pole is at sqrt(0.00600069867) um = 77.46 nm; the visible
        // is far above it and n is an ordinary 1.5.
        let wl = array![380.0, 550.0, 780.0, 2000.0];
        assert!(check(&wl).is_ok());
        let out = sellmeier_nk(wl.view(), BK7[0], BK7[1], BK7[2], BK7[3], BK7[4], BK7[5]);
        assert!((out[1].re - 1.5185).abs() < 1e-3, "n(550) = {}", out[1].re);
    }

    #[test]
    fn the_domain_edge_is_below_the_first_pole_not_at_it() {
        // Worth pinning because it is counter-intuitive: n^2 stays positive for a
        // stretch below the 77.46 nm resonance (the other two terms hold it up),
        // so 50 nm evaluates to a perfectly finite n = 0.47 -- nonsense physically,
        // but not a domain error, and this guard does not pretend otherwise. The
        // arithmetic only breaks once term 1 dominates, at 70.7 nm and shorter.
        let out = sellmeier_nk(
            array![50.0].view(),
            BK7[0],
            BK7[1],
            BK7[2],
            BK7[3],
            BK7[4],
            BK7[5],
        );
        assert!(
            out[0].re.is_finite() && out[0].re < 1.0,
            "n(50 nm) = {}",
            out[0].re
        );
        assert!(check(&array![50.0]).is_ok());
    }

    #[test]
    fn a_wavelength_inside_the_first_resonance_is_refused() {
        // 70 nm: lambda^2 = 0.0049 um^2, just inside C1 = 0.0060, where term 1
        // goes large and negative, n^2 < 0 and sqrt gives NaN.
        let e = check(&array![550.0, 70.0]).expect_err("NaN must be refused");
        assert!(e.contains("70"), "{e}");
        assert!(e.contains("NaN"), "{e}");
        // The message has to name where the edge is, or the user cannot act on it.
        assert!(
            e.contains("77.4"),
            "the nearest resonance must be named: {e}"
        );
    }

    #[test]
    fn landing_exactly_on_a_pole_says_so() {
        // Reaching a pole exactly needs lambda^2 to hit C bit-for-bit, which a
        // round trip through `sqrt(C)*1000` does not give -- it lands a few ulp
        // short and shows up as the NaN case instead. Constructed the other way
        // round: C2 = 0.25 um^2 is exactly `wl_um2(500 nm)`.
        let wl = array![500.0];
        let out = sellmeier_nk(wl.view(), 1.0, 0.01, 0.3, 0.25, 0.0, 0.0);
        assert!(out[0].re.is_infinite(), "n = {}", out[0].re);
        let e = sellmeier_domain_check(wl.view(), &out, 1.0, 0.01, 0.3, 0.25, 0.0, 0.0)
            .expect_err("a pole must be refused");
        assert!(e.contains("pole"), "{e}");
        assert!(e.contains("infinite"), "{e}");
    }

    #[test]
    fn the_message_is_ascii() {
        // It travels to Python as a ValueError, and printing it on a cp1252
        // console must not raise a UnicodeEncodeError of its own.
        let e = check(&array![70.0]).err().unwrap();
        assert!(e.is_ascii(), "non-ASCII in: {e}");
    }

    #[test]
    fn a_dropped_third_term_is_not_a_pole() {
        // B3 = 0 removes the term (see `sellmeier_n`), so C3 is not a resonance
        // and a grid sitting on it is perfectly fine.
        let c3: f64 = 0.25;
        let wl = array![c3.sqrt() * 1000.0];
        let out = sellmeier_nk(wl.view(), 1.0, 0.01, 0.3, 0.05, 0.0, c3);
        assert!(out[0].re.is_finite(), "n = {}", out[0].re);
        assert!(sellmeier_domain_check(wl.view(), &out, 1.0, 0.01, 0.3, 0.05, 0.0, c3).is_ok());
    }

    #[test]
    fn the_urbach_variant_is_checked_on_n_not_k() {
        let wl = array![550.0];
        // A huge Urbach tail is a large k, not a domain error.
        let out = sellmeier_urbach_nk(wl.view(), 1.0, 0.01, 0.3, 0.05, 0.0, 0.0, 1e9, 0.06, 380.0);
        assert!(out[0].im > 0.0);
        assert!(sellmeier_domain_check(wl.view(), &out, 1.0, 0.01, 0.3, 0.05, 0.0, 0.0).is_ok());
    }
}
