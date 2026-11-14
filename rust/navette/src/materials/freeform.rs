// SPDX-License-Identifier: LGPL-3.0-or-later
//! Free-form optical constants for materials no closed-form model fits
//! (`docs/plans/general_fit_plan.md` §4.6).
//!
//! * [`point_wise_nk`]: n and k at given wavelengths, linear in between and
//!   held at the ends. With one point per data wavelength every value is a
//!   fit parameter; it needs several independent measurements per
//!   wavelength and a smoothness prior to be determined.
//! * [`spline_nk`]: n and k on a coarse knot grid through
//!   [`UniInterpolator`] (pchip by default, makima or linear), held at the
//!   ends. Far fewer parameters, smooth by construction. pchip keeps a
//!   non-negative k non-negative; makima may undershoot between knots.
//! * [`bspline_eps_nk`]: ε₂ as a cubic B-spline in photon energy with
//!   non-negative coefficients, ε₁ from the shared FFT Kramers-Kronig
//!   transform ([`kk`]), so the model is KK-consistent by construction.
//!   Coefficient `i` sits near `energies[i]`; ε₂ is zero outside the first
//!   and last energy, so the end coefficients are usually kept at zero.
//!
//! k and ε₂ are read as `max(value, 0)`: a knot fitted onto its zero bound
//! is probed just below it by finite differences, and makima may dip
//! under zero between knots; the solver refuses a negative k (gain).

use ndarray::{Array1, ArrayView1};
use num_complex::Complex64;

use super::kk;
use super::units::energy_ev;
use crate::interpolate::UniInterpolator;

/// Strictly increasing, finite, at least `min` points.
fn check_grid(what: &str, x: &[f64], min: usize) -> Result<(), String> {
    if x.len() < min {
        return Err(format!(
            "{what} needs at least {min} points, got {}",
            x.len()
        ));
    }
    if x.iter().any(|v| !v.is_finite()) {
        return Err(format!("{what} must be finite"));
    }
    if x.windows(2).any(|w| w[1] <= w[0]) {
        return Err(format!("{what} must be strictly increasing"));
    }
    Ok(())
}

fn check_values(what: &str, v: &[f64], n: usize) -> Result<(), String> {
    if v.len() != n {
        return Err(format!("{what} has {} values for {n} points", v.len()));
    }
    if v.iter().any(|x| !x.is_finite()) {
        return Err(format!("{what} must be finite"));
    }
    Ok(())
}

/// n and k at `points` (nm), linear in between and held at the ends;
/// `k` absent is zero.
pub fn point_wise_nk(
    wavelength_nm: ArrayView1<f64>,
    points: &[f64],
    n: &[f64],
    k: Option<&[f64]>,
) -> Result<Array1<Complex64>, String> {
    check_grid("PointWise wavelengths", points, 2)?;
    check_values("PointWise n", n, points.len())?;
    if let Some(k) = k {
        check_values("PointWise k", k, points.len())?;
    }
    Ok(wavelength_nm.mapv(|w| {
        let kv = k.map_or(0.0, |k| kk::interp(w, points, k).max(0.0));
        Complex64::new(kk::interp(w, points, n), kv)
    }))
}

/// n and k on `knots` (nm) by `method` (`pchip`, `makima`, `linear`),
/// held at the ends; `k` absent is zero.
pub fn spline_nk(
    wavelength_nm: ArrayView1<f64>,
    knots: &[f64],
    n: &[f64],
    k: Option<&[f64]>,
    method: &str,
) -> Result<Array1<Complex64>, String> {
    check_grid("Spline wavelengths", knots, 2)?;
    if !matches!(method, "pchip" | "makima" | "linear") {
        return Err(format!(
            "Spline method {method:?}: use \"pchip\", \"makima\" or \"linear\""
        ));
    }
    check_values("Spline n", n, knots.len())?;
    let rows = if let Some(k) = k {
        check_values("Spline k", k, knots.len())?;
        2
    } else {
        1
    };
    // The interpolator speaks numpy's ndarray (0.17).
    let mut flat = n.to_vec();
    if let Some(k) = k {
        flat.extend_from_slice(k);
    }
    let y =
        ndarray017::Array2::from_shape_vec((rows, knots.len()), flat).map_err(|e| e.to_string())?;
    let it = UniInterpolator::new(
        ndarray017::Array1::from_vec(knots.to_vec()),
        y,
        rows > 1,
        method,
        false,
        3,
        "clamp",
    )?;
    let wl: Vec<f64> = wavelength_nm.to_vec();
    let v = it.evaluate(&wl, 0, None)?;
    Ok((0..wl.len())
        .map(|i| Complex64::new(v[(0, i)], if rows > 1 { v[(1, i)].max(0.0) } else { 0.0 }))
        .collect())
}

/// The clamped cubic knot vector for one coefficient per energy: the end
/// energies four times, the inner ones but the first and last once (the
/// not-a-knot placement).
fn knot_vector(energies: &[f64]) -> Vec<f64> {
    let k = energies.len();
    let mut t = vec![energies[0]; 4];
    t.extend_from_slice(&energies[2..k - 2]);
    t.extend(std::iter::repeat_n(energies[k - 1], 4));
    t
}

/// The cubic B-spline with knots `t` and coefficients `c` at `x` (de Boor);
/// zero outside `[t[3], t[len c]]`.
fn de_boor(t: &[f64], c: &[f64], x: f64) -> f64 {
    let k = c.len();
    if x < t[3] || x > t[k] {
        return 0.0;
    }
    // The span: t[i] <= x < t[i + 1], the last one closed.
    let i = (t.partition_point(|&v| v <= x) - 1).clamp(3, k - 1);
    let mut d = [c[i - 3], c[i - 2], c[i - 1], c[i]];
    for r in 1..=3 {
        for j in (r..=3).rev() {
            let lo = t[j + i - 3];
            let hi = t[j + 1 + i - r];
            let a = if hi > lo { (x - lo) / (hi - lo) } else { 0.0 };
            d[j] = (1.0 - a) * d[j - 1] + a * d[j];
        }
    }
    d[3]
}

/// ε₂ as a cubic B-spline over `energies` (eV) with coefficients `eps2`,
/// ε₁ = `eps_inf` − Hilbert{ε₂} on the KK grid; returns n + ik.
pub fn bspline_eps_nk(
    wavelength_nm: ArrayView1<f64>,
    energies: &[f64],
    eps2: &[f64],
    eps_inf: f64,
) -> Result<Array1<Complex64>, String> {
    check_grid("BSplineEps energies", energies, 4)?;
    check_values("BSplineEps eps2", eps2, energies.len())?;
    let (lo, hi) = (kk::GRID_MIN, kk::GRID_MAX);
    if energies[0] < lo || energies[energies.len() - 1] > hi {
        return Err(format!(
            "BSplineEps energies must lie in the KK grid [{lo}, {hi}] eV"
        ));
    }
    let t = knot_vector(energies);
    let e_full = kk::energy_grid();
    let eps2_at = |e: f64| de_boor(&t, eps2, e).max(0.0);
    let eps2_full: Vec<f64> = e_full.iter().map(|&e| eps2_at(e)).collect();
    let eps1_full = kk::kk_fft(&eps2_full, eps_inf);
    let mut out = Vec::with_capacity(wavelength_nm.len());
    for &w in wavelength_nm {
        let e = energy_ev(w);
        if !(lo..=hi).contains(&e) {
            return Err(format!(
                "BSplineEps: {w} nm ({e:.4} eV) is outside the KK grid [{lo}, {hi}] eV"
            ));
        }
        let e1 = kk::interp(e, &e_full, &eps1_full);
        out.push(Complex64::new(e1, eps2_at(e)).sqrt());
    }
    Ok(Array1::from(out))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wl() -> Array1<f64> {
        Array1::linspace(300.0, 900.0, 61)
    }

    #[test]
    fn point_wise_is_exact_at_its_points_and_held_outside() {
        let p = [400.0, 600.0, 800.0];
        let v = point_wise_nk(
            ArrayView1::from(&[400.0, 500.0, 600.0, 200.0, 1000.0]),
            &p,
            &[2.0, 1.8, 1.7],
            Some(&[0.1, 0.0, 0.0]),
        )
        .unwrap();
        assert_eq!(v[0], Complex64::new(2.0, 0.1));
        assert_eq!(v[1], Complex64::new(1.9, 0.05));
        assert_eq!(v[2], Complex64::new(1.8, 0.0));
        assert_eq!(v[3], v[0]);
        assert_eq!(v[4], Complex64::new(1.7, 0.0));
        // A knot just under zero reads as zero.
        let v = point_wise_nk(
            ArrayView1::from(&[600.0]),
            &p,
            &[2.0, 1.8, 1.7],
            Some(&[0.1, -1e-9, 0.0]),
        )
        .unwrap();
        assert_eq!(v[0].im, 0.0);
        assert!(point_wise_nk(wl().view(), &p, &[2.0], None).is_err());
        assert!(point_wise_nk(wl().view(), &[400.0, 400.0], &[2.0, 2.0], None).is_err());
    }

    #[test]
    fn a_spline_passes_its_knots_and_pchip_keeps_k_non_negative() {
        let knots = [300.0, 400.0, 500.0, 700.0, 900.0];
        let n = [2.4, 2.2, 2.1, 2.05, 2.02];
        let k = [0.3, 0.05, 0.0, 0.0, 0.0];
        let at = spline_nk(ArrayView1::from(&knots), &knots, &n, Some(&k), "pchip").unwrap();
        for (i, z) in at.iter().enumerate() {
            assert!((z.re - n[i]).abs() < 1e-12 && (z.im - k[i]).abs() < 1e-12);
        }
        let v = spline_nk(wl().view(), &knots, &n, Some(&k), "pchip").unwrap();
        assert!(v.iter().all(|z| z.im >= 0.0));
        assert!(spline_nk(wl().view(), &knots, &n, None, "cubic").is_err());
        assert!(spline_nk(wl().view(), &knots, &n, None, "makima").is_ok());
    }

    #[test]
    fn the_b_spline_is_a_partition_of_unity_and_kk_consistent() {
        // Equal coefficients: ε₂ is that constant inside, zero outside.
        let e = [1.0, 2.0, 3.0, 4.5, 6.0, 8.0];
        let t = knot_vector(&e);
        for x in [1.0, 1.7, 3.3, 7.9, 8.0] {
            assert!((de_boor(&t, &[2.0; 6], x) - 2.0).abs() < 1e-12, "{x}");
        }
        assert_eq!(de_boor(&t, &[2.0; 6], 0.5), 0.0);
        assert_eq!(de_boor(&t, &[2.0; 6], 8.5), 0.0);
        // Non-negative coefficients give non-negative ε₂ (the convex hull).
        let c = [0.0, 0.5, 3.0, 1.0, 0.2, 0.0];
        assert!((0..400).all(|i| de_boor(&t, &c, 0.5 + 0.02 * i as f64) >= 0.0));
        // An absorption band in the UV raises n below it (normal dispersion).
        let below = Array1::linspace(350.0, 900.0, 56);
        let band = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let v = bspline_eps_nk(below.view(), &band, &[0.0, 2.0, 4.0, 4.0, 2.0, 0.0], 1.0).unwrap();
        assert!(v.iter().all(|z| z.re > 1.0));
        let n: Vec<f64> = v.iter().map(|z| z.re).collect();
        assert!(n.windows(2).all(|w| w[0] >= w[1]), "n falls toward the red");
        assert!(bspline_eps_nk(wl().view(), &[4.0, 5.0, 6.0], &[0.0; 3], 1.0).is_err());
    }
}
