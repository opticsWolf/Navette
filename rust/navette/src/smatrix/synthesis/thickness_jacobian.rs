// SPDX-License-Identifier: LGPL-3.0-or-later
//! Navette -- Rust Rewrite of Numba-optimized thin-film optical solver
//!
//! synthesis::thickness_jacobian — the analytic Jacobian of the thickness optimizer.
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

use crate::smatrix::synthesis::merit::{CurveId, CurveTerm, MeritSensitivity, TermPart};

/// The four front intensity channels the coherent sweep produces, in the
/// order [`CurveDeposits`] stores them. Absorption (A = 1 − R − T) reads
/// them through its companions and the unpolarized curves through their
/// s/p halves ([`CurveId::unpolarized_parts`], split at assembly). Back
/// curves have their own block, [`BACK_DEPOSIT_CHANNELS`].
pub const DEPOSIT_CHANNELS: [CurveId; 4] = [CurveId::Rs, CurveId::Rp, CurveId::Ts, CurveId::Tp];

/// The complex amplitudes a phase demand reads, whose real and imaginary
/// parts follow the intensity channels when the spec asks for phase:
/// slots 4, 5 (Re, Im of `Ts`), 6, 7 (`Tp`), 8, 9 (`Rs`) and 10, 11 (`Rp`).
pub const PHASE_DEPOSIT_CHANNELS: [CurveId; 4] =
    [CurveId::Ts, CurveId::Tp, CurveId::Rs, CurveId::Rp];

/// The back-incidence curves, present when the spec reads any: their
/// intensities in slots 12-15, then the Re, Im of their complex amplitudes
/// in slots 16-23 (16, 17 for `RBs`, …, 22, 23 for `TBp`). The front phase
/// slots before them are always filled when these are.
pub const BACK_DEPOSIT_CHANNELS: [CurveId; 4] =
    [CurveId::RBs, CurveId::RBp, CurveId::TBs, CurveId::TBp];

/// First slot of the back intensities; the back amplitudes follow at
/// `BACK_SLOT + 4`.
pub const BACK_SLOT: usize = 4 + 2 * PHASE_DEPOSIT_CHANNELS.len();

/// Slot count of a sweep that carries the back channels.
pub const ALL_SLOTS: usize = BACK_SLOT + 3 * BACK_DEPOSIT_CHANNELS.len();

/// Storage slot for a curve, or `None` when this pass carries no deposits
/// for it.
pub fn deposit_channel(id: CurveId) -> Option<usize> {
    match id {
        CurveId::Rs => Some(0),
        CurveId::Rp => Some(1),
        CurveId::Ts => Some(2),
        CurveId::Tp => Some(3),
        _ => BACK_DEPOSIT_CHANNELS
            .iter()
            .position(|&c| c == id)
            .map(|k| BACK_SLOT + k),
    }
}

/// First slot of a curve's complex amplitude (Re; Im follows), or `None`.
fn amplitude_channel(id: CurveId) -> Option<usize> {
    match PHASE_DEPOSIT_CHANNELS.iter().position(|&c| c == id) {
        Some(k) => Some(4 + 2 * k),
        None => BACK_DEPOSIT_CHANNELS
            .iter()
            .position(|&c| c == id)
            .map(|k| BACK_SLOT + 4 + 2 * k),
    }
}

/// Storage slot for one sensitivity term, or `None` when no deposits exist
/// for it. A [`TermPart::TotalThickness`] term has no slot: ∂D/∂d is 1 for
/// every film and the assembly adds it directly.
pub(crate) fn term_channel(term: &CurveTerm) -> Option<usize> {
    match term.part {
        TermPart::Value => deposit_channel(term.curve),
        TermPart::Re => amplitude_channel(term.curve),
        TermPart::Im => amplitude_channel(term.curve).map(|c| c + 1),
        TermPart::TotalThickness => None,
    }
}

/// ∂(simulated curve value)/∂(film thickness), per channel, grid point and
/// optimized film.
///
/// A "grid point" is the flat angle-major index the curves themselves use,
/// `point = angle_row · n_wav + wavelength`; a "parameter" is a position in
/// the optimizer's parameter vector, i.e. the k-th optimize-flagged film.
///
/// Four channels ([`DEPOSIT_CHANNELS`]) for an intensity-only spec, twelve
/// when the amplitudes of [`PHASE_DEPOSIT_CHANNELS`] follow, and
/// [`ALL_SLOTS`] when the spec reads a back curve
/// ([`BACK_DEPOSIT_CHANNELS`]).
#[derive(Clone, Debug)]
pub struct CurveDeposits {
    n_points: usize,
    n_par: usize,
    /// `data[channel][point · n_par + par]`, channels as [`term_channel`].
    data: Vec<Vec<f64>>,
}

impl CurveDeposits {
    /// Build from one row per grid point, each row laid out
    /// `[channel · n_par + par]` (4·n_par entries, 4 + 2k channels with
    /// the first k of [`PHASE_DEPOSIT_CHANNELS`], or [`ALL_SLOTS`]) — the
    /// shape the parallel sweep collects.
    pub fn from_point_rows(rows: Vec<Vec<f64>>, n_par: usize) -> Result<Self, String> {
        let n_points = rows.len();
        let max_ch = 4 + 2 * PHASE_DEPOSIT_CHANNELS.len();
        let n_ch = match rows.first() {
            Some(r)
                if n_par > 0
                    && r.len() % n_par == 0
                    && (4..=max_ch)
                        .step_by(2)
                        .chain([ALL_SLOTS])
                        .any(|c| c * n_par == r.len()) =>
            {
                r.len() / n_par
            }
            _ => 4,
        };
        let stride = n_ch * n_par;
        if rows.iter().any(|r| r.len() != stride) {
            return Err(format!(
                "curve deposits: every point row must hold {n_ch}·n_par = {stride} entries"
            ));
        }
        let mut data = vec![vec![0.0; n_points * n_par]; n_ch];
        for (p, row) in rows.iter().enumerate() {
            for (c, block) in data.iter_mut().enumerate() {
                block[p * n_par..(p + 1) * n_par].copy_from_slice(&row[c * n_par..(c + 1) * n_par]);
            }
        }
        Ok(CurveDeposits {
            n_points,
            n_par,
            data,
        })
    }

    pub fn n_points(&self) -> usize {
        self.n_points
    }

    pub fn n_par(&self) -> usize {
        self.n_par
    }

    /// 4 (intensities), or more with the phase amplitudes and the back
    /// channels after them.
    pub fn n_channels(&self) -> usize {
        self.data.len()
    }

    /// True when every term of `sens` has a deposit channel here (after
    /// the unpolarized split). A false answer is a decline, not an error:
    /// the caller differences that iteration.
    pub fn covers(&self, sens: &MeritSensitivity) -> bool {
        let has = |t: &CurveTerm| term_channel(t).is_some_and(|c| c < self.data.len());
        sens.rows.iter().flatten().all(|t| {
            t.part == TermPart::TotalThickness
                || match (t.part, t.curve.unpolarized_parts()) {
                    (TermPart::Value, Some((a, b))) => {
                        has(&CurveTerm { curve: a, ..*t }) && has(&CurveTerm { curve: b, ..*t })
                    }
                    _ => has(t),
                }
        })
    }

    /// `out += d_residual · ∂(curve)/∂θ` for one non-thickness term. An
    /// unpolarized intensity is `½(X_s + X_p)` by the evaluator's contract,
    /// so its column is the half-sum of the two polarized columns.
    fn add_term(
        &self,
        out: &mut [f64],
        term: &CurveTerm,
        row: usize,
        n_wav: usize,
    ) -> Result<(), String> {
        let d = term.d_residual;
        match (term.part, term.curve.unpolarized_parts()) {
            (TermPart::Value, Some((a, b))) => {
                let half = 0.5 * d;
                for part in [a, b] {
                    let col = self.term_column(
                        &CurveTerm {
                            curve: part,
                            ..*term
                        },
                        row,
                        n_wav,
                    )?;
                    for (o, &g) in out.iter_mut().zip(col) {
                        *o += half * g;
                    }
                }
            }
            _ => {
                for (o, &g) in out.iter_mut().zip(self.term_column(term, row, n_wav)?) {
                    *o += d * g;
                }
            }
        }
        Ok(())
    }

    /// The deposit column a term reads, or the reason there is none.
    fn term_column(&self, term: &CurveTerm, row: usize, n_wav: usize) -> Result<&[f64], String> {
        let channel = term_channel(term)
            .filter(|&c| c < self.data.len())
            .ok_or_else(|| {
                format!(
                    "analytic jacobian: no deposits for {:?} of curve {:?} (row {row})",
                    term.part, term.curve
                )
            })?;
        let point = term.angle_row * n_wav + term.wavelength;
        if point >= self.n_points {
            return Err(format!(
                "analytic jacobian: row {row} reads grid point {point}, but the \
                 deposits cover {} points",
                self.n_points
            ));
        }
        Ok(self.column(channel, point))
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
/// Refuses rather than guesses: a spec with an uncovered row (a color
/// demand) has no analytic Jacobian here and must go to the
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
             (color demand) — use the finite-difference path"
        ));
    }
    let m = sens.rows.len();
    let n = dep.n_par();
    jac.clear();
    jac.resize(m * n, 0.0);
    for (i, row) in sens.rows.iter().enumerate() {
        let out = &mut jac[i * n..(i + 1) * n];
        for term in row {
            let d = term.d_residual;
            if term.part == TermPart::TotalThickness {
                // ∂D/∂d = 1 for every film.
                for o in out.iter_mut() {
                    *o += d;
                }
                continue;
            }
            dep.add_term(out, term, i, n_wav)?;
        }
    }
    Ok(())
}

/// F1.6: the parameter-map assembly. `params` is the same list
/// `refine_thicknesses_report` built (each entry either one row or one
/// scalable span's bulk rows with frozen fractions); the deposits are per
/// FILM ROW, laid out in parameter order, so the assembly contracts each
/// residual's per-row columns into per-parameter columns.
///
/// Bit-exactness, asserted not assumed (the no-span twin):
///
/// - A `Row` parameter is a one-element group with weight `1.0`; its
///   output column is bitwise the row column `assemble_jacobian` would
///   have written (the same additions in the same order — `1.0 * g` is
///   `g` to the last bit, and here the multiply is skipped entirely).
/// - A `Span` parameter's column is `sum_r phi_r * col_r` accumulated in
///   the span's row order — exactly the closed form the Jacobian twin
///   checks against the row columns.
pub(crate) fn assemble_jacobian_mapped(
    sens: &MeritSensitivity,
    dep: &CurveDeposits,
    n_wav: usize,
    params: &[crate::smatrix::synthesis::thickness_refine::ThicknessParam],
    jac: &mut Vec<f64>,
) -> Result<(), String> {
    if let Some(&row) = sens.uncovered.first() {
        return Err(format!(
            "analytic jacobian: residual row {row} has no curve sensitivity \
             (color demand) — use the finite-difference path"
        ));
    }
    let m = sens.rows.len();
    let n_out = params.len();
    let n_flat = dep.n_par();
    // Each parameter's base position in the FLAT deposit list: the
    // deposits are laid out in parameter order (each parameter's rows
    // consecutively), so the map is positional, not by film row.
    let mut bases = Vec::with_capacity(params.len() + 1);
    let mut run = 0usize;
    for p in params {
        bases.push(run);
        match p {
            crate::smatrix::synthesis::thickness_refine::ThicknessParam::Row(_) => run += 1,
            crate::smatrix::synthesis::thickness_refine::ThicknessParam::Span { rows, .. } => {
                run += rows.len()
            }
        }
    }
    if run != n_flat {
        return Err(format!(
            "analytic jacobian: the parameter map covers {} row(s), but the \
             deposits cover {n_flat} - the flat list and the map disagree",
            run
        ));
    }
    jac.clear();
    jac.resize(m * n_out, 0.0);
    // One residual row at a time: contract the deposits into per-ROW
    // columns (exactly `assemble_jacobian`'s inner loop, same additions in
    // the same order), then apply the parameter map.
    let mut scratch = vec![0.0f64; n_flat];
    for (i, row) in sens.rows.iter().enumerate() {
        for v in scratch.iter_mut() {
            *v = 0.0;
        }
        for term in row {
            let d = term.d_residual;
            if term.part == TermPart::TotalThickness {
                // ∂D/∂(row thickness) = 1 for every film row; a span's
                // fractions then sum it to 1 per span parameter.
                for s in scratch.iter_mut() {
                    *s += d;
                }
                continue;
            }
            dep.add_term(&mut scratch, term, i, n_wav)?;
        }
        let out = &mut jac[i * n_out..(i + 1) * n_out];
        for (p, par) in params.iter().enumerate() {
            match par {
                crate::smatrix::synthesis::thickness_refine::ThicknessParam::Row(_) => {
                    out[p] = scratch[bases[p]]
                }
                crate::smatrix::synthesis::thickness_refine::ThicknessParam::Span {
                    rows,
                    fractions,
                } => {
                    let mut acc = 0.0f64;
                    for (k, &w) in fractions.iter().enumerate() {
                        acc += w * scratch[bases[p] + k];
                    }
                    let _ = rows; // rows.len() == fractions.len(), checked by construction
                    out[p] = acc;
                }
            }
        }
    }
    Ok(())
}
