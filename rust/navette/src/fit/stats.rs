// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::stats — covariance and correlations at a solution.
//!
//! `C = (JᵀJ)⁻¹ · χ²_red`, with `JᵀJ = RᵀR` from the LM's own Householder
//! QR of J (never the normal equations: the thickness ↔ index degeneracy
//! this is meant to expose is exactly where `JᵀJ` loses half its digits).
//! J is taken in the solver's normalised coordinates and mapped back to
//! the parameters by their `dx/du`.

use crate::lsq::lm::{back_substitute, householder_qr_in_place};

/// `(JᵀJ)⁻¹` of a row-major `m × n` Jacobian (`m >= n`), row-major `n × n`.
/// A column that does not move the residuals refuses, naming its index.
pub(crate) fn inverse_normal_matrix(jac: &[f64], m: usize, n: usize) -> Result<Vec<f64>, String> {
    if m < n {
        return Err(format!(
            "{m} residuals cannot determine {n} free parameters"
        ));
    }
    let mut a = jac.to_vec();
    let mut b = vec![0.0; m];
    householder_qr_in_place(&mut a, &mut b, m, n).map_err(|_| {
        "a free parameter does not move the residuals (rank-deficient Jacobian)".to_string()
    })?;
    // R⁻¹ column by column, then (RᵀR)⁻¹ = R⁻¹ R⁻ᵀ.
    let mut rinv = vec![0.0; n * n]; // row-major
    let mut e = vec![0.0; n];
    let mut col = vec![0.0; n];
    for j in 0..n {
        e.fill(0.0);
        e[j] = 1.0;
        back_substitute(&a, &e, n, &mut col).map_err(|_| {
            "singular Jacobian: the free parameters are not independent".to_string()
        })?;
        for i in 0..n {
            rinv[i * n + j] = col[i];
        }
    }
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            c[i * n + k] = (0..n).map(|j| rinv[i * n + j] * rinv[k * n + j]).sum();
        }
    }
    Ok(c)
}

/// Correlation matrix of a covariance matrix.
pub(crate) fn correlation(cov: &[f64], n: usize) -> Vec<f64> {
    let sd: Vec<f64> = (0..n).map(|i| cov[i * n + i].sqrt()).collect();
    let mut out = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            out[i * n + k] = cov[i * n + k] / (sd[i] * sd[k]);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse_of_a_known_normal_matrix() {
        // J = [[1, 0], [0, 2], [1, 1]] → JᵀJ = [[2, 1], [1, 5]], inverse / 9.
        let jac = [1.0, 0.0, 0.0, 2.0, 1.0, 1.0];
        let c = inverse_normal_matrix(&jac, 3, 2).unwrap();
        let want = [5.0 / 9.0, -1.0 / 9.0, -1.0 / 9.0, 2.0 / 9.0];
        for (got, want) in c.iter().zip(want) {
            assert!((got - want).abs() < 1e-14, "{c:?}");
        }
        let r = correlation(&c, 2);
        assert!((r[1] + 1.0 / 10f64.sqrt()).abs() < 1e-14);
        assert!((r[0] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn a_dead_column_refuses() {
        let jac = [1.0, 0.0, 2.0, 0.0, 3.0, 0.0];
        assert!(
            inverse_normal_matrix(&jac, 3, 2)
                .unwrap_err()
                .contains("does not move")
        );
        assert!(
            inverse_normal_matrix(&jac, 1, 2)
                .unwrap_err()
                .contains("cannot determine")
        );
    }
}
