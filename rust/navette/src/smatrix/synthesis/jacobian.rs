//! Navette -- Rust Rewrite of Numba-optimized thin-film optical solver
//!
//! synthesis::jacobian — the analytic Jacobian of the thickness optimizer.
//!
//! ```text
//!     J[i,k] = Σ_terms  ∂r_i/∂(curve value) · ∂(curve value)/∂d_k
//! ```
//!
//! The left factor is [`MeritSensitivity`] (R4.5 increment i, `merit.rs`):
//! which simulated curve values residual row `i` reads, and with what
//! derivative. The right factor is [`CurveDeposits`]: how each of those
//! values moves when a film's thickness moves, which
//! `SmatrixContext::simulate_with_deposits` reads off the same solver sweep
//! that produced the curves.
//!
//! Splitting it in two is what makes it testable. Each half has an exact,
//! independent finite-difference oracle — `residuals()` differenced in curve
//! space for the left, `simulate()` differenced in thickness space for the
//! right — so a disagreement in the product has only one place left to be.
//!
//! **Determinism.** `J` is materialized m×n and every element accumulates
//! over its own row's terms in a fixed order: no cross-point reduction, no
//! parallel accumulation, nothing whose order depends on scheduling (§13).

use crate::smatrix::synthesis::merit::{CurveId, MeritSensitivity};

/// The four front intensity channels the coherent sweep produces, in the
/// order [`CurveDeposits`] stores them. Everything else is either derived
/// from these (absorption: A = 1 − R − T) or not simulated at all, in which
/// case the residual pass has already failed on the missing curve.
pub const DEPOSIT_CHANNELS: [CurveId; 4] =
    [CurveId::Rs, CurveId::Rp, CurveId::Ts, CurveId::Tp];

/// Storage slot for a curve, or `None` when this pass carries no deposits
/// for it.
pub fn deposit_channel(id: CurveId) -> Option<usize> {
    match id {
        CurveId::Rs => Some(0),
        CurveId::Rp => Some(1),
        CurveId::Ts => Some(2),
        CurveId::Tp => Some(3),
        _ => None,
    }
}

/// ∂(simulated curve value)/∂(film thickness), per channel, grid point and
/// optimized film.
///
/// A "grid point" is the flat angle-major index the curves themselves use,
/// `point = angle_row · n_wav + wavelength`; a "parameter" is a position in
/// the optimizer's parameter vector, i.e. the k-th optimize-flagged film.
#[derive(Clone, Debug)]
pub struct CurveDeposits {
    n_points: usize,
    n_par: usize,
    /// `data[channel][point · n_par + par]`, channels as [`DEPOSIT_CHANNELS`].
    data: [Vec<f64>; 4],
}

impl CurveDeposits {
    /// Build from one row per grid point, each row laid out
    /// `[channel · n_par + par]` (4·n_par entries) — the shape the parallel
    /// sweep collects.
    pub fn from_point_rows(rows: Vec<Vec<f64>>, n_par: usize) -> Result<Self, String> {
        let n_points = rows.len();
        let stride = 4 * n_par;
        if rows.iter().any(|r| r.len() != stride) {
            return Err(format!(
                "curve deposits: every point row must hold 4·n_par = {stride} entries"
            ));
        }
        let mut data = [
            vec![0.0; n_points * n_par],
            vec![0.0; n_points * n_par],
            vec![0.0; n_points * n_par],
            vec![0.0; n_points * n_par],
        ];
        for (p, row) in rows.iter().enumerate() {
            for (c, block) in data.iter_mut().enumerate() {
                block[p * n_par..(p + 1) * n_par]
                    .copy_from_slice(&row[c * n_par..(c + 1) * n_par]);
            }
        }
        Ok(CurveDeposits { n_points, n_par, data })
    }

    pub fn n_points(&self) -> usize {
        self.n_points
    }

    pub fn n_par(&self) -> usize {
        self.n_par
    }

    /// ∂(curve value at `point`)/∂(thickness of parameter `par`).
    pub fn get(&self, channel: usize, point: usize, par: usize) -> f64 {
        self.data[channel][point * self.n_par + par]
    }

    /// One parameter's whole column for a channel — the slice the assembly
    /// walks per residual term.
    fn column(&self, channel: usize, point: usize) -> &[f64] {
        &self.data[channel][point * self.n_par..(point + 1) * self.n_par]
    }
}

/// Assemble the m×n Jacobian, row-major, from the two halves.
///
/// `n_wav` is the simulated grid's wavelength count — the stride the
/// sensitivity's `(angle_row, wavelength)` addresses flatten by, the same one
/// `SimCurves` uses.
///
/// Refuses rather than guesses: a spec with an uncovered row (a phase target
/// or a color demand) has no analytic Jacobian here and must go to the
/// finite-difference path. Ask with [`MeritSensitivity::is_complete`] first —
/// the error is for a caller that did not.
pub fn assemble_jacobian(
    sens: &MeritSensitivity,
    dep: &CurveDeposits,
    n_wav: usize,
    jac: &mut Vec<f64>,
) -> Result<(), String> {
    if let Some(&row) = sens.uncovered.first() {
        return Err(format!(
            "analytic jacobian: residual row {row} has no curve sensitivity \
             (phase target or color demand) — use the finite-difference path"
        ));
    }
    let m = sens.rows.len();
    let n = dep.n_par();
    jac.clear();
    jac.resize(m * n, 0.0);
    for (i, row) in sens.rows.iter().enumerate() {
        let out = &mut jac[i * n..(i + 1) * n];
        for term in row {
            let channel = deposit_channel(term.curve).ok_or_else(|| {
                format!(
                    "analytic jacobian: no deposits for curve {:?} (row {i})",
                    term.curve
                )
            })?;
            let point = term.angle_row * n_wav + term.wavelength;
            if point >= dep.n_points() {
                return Err(format!(
                    "analytic jacobian: row {i} reads grid point {point}, but the \
                     deposits cover {} points",
                    dep.n_points()
                ));
            }
            let d = term.d_residual;
            for (o, &g) in out.iter_mut().zip(dep.column(channel, point)) {
                *o += d * g;
            }
        }
    }
    Ok(())
}
