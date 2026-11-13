// SPDX-License-Identifier: LGPL-3.0-or-later
//! synthesis::thickness_refine — the thickness parameters of the LM
//! refinement inside needle synthesis, and their analytic Jacobian.
//!
//! The specific path: one parameter per optimize-flagged film, or one
//! total per scalable span. The solvers it hands them to live in
//! [`crate::lsq`] and know nothing about thicknesses.

use std::sync::atomic::{AtomicBool, Ordering};

use crate::lsq::lm::JacobianSource;
use crate::smatrix::synthesis::evaluator::SmatrixContext;
use crate::smatrix::synthesis::structure::DesignStack;
use crate::smatrix::synthesis::thickness_jacobian::assemble_jacobian_mapped;

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
pub(crate) enum ThicknessParam {
    Row(usize),
    Span {
        rows: Vec<usize>,
        fractions: Vec<f64>,
    },
}

impl ThicknessParam {
    /// The parameter's starting value: the row's thickness, or the span
    /// total D0 the fractions were frozen at.
    pub(crate) fn value0(&self, stack: &DesignStack) -> f64 {
        match self {
            ThicknessParam::Row(i) => stack.films()[*i].d_nm,
            ThicknessParam::Span { rows, .. } => rows.iter().map(|&r| stack.films()[r].d_nm).sum(),
        }
    }
}

/// The parameter list for one `refine_thicknesses` call.
///
/// Span parameters first-class: each scalable span contributes ONE
/// parameter (its bulk rows in row order, fractions frozen from the
/// current row thicknesses). Every other optimize-flagged row is a `Row`
/// parameter. A stack with no scalable span produces exactly today's
/// `opt_indices` list, as `Row` params.
pub(crate) fn thickness_params(stack: &DesignStack) -> Vec<ThicknessParam> {
    let films = stack.films();
    let mut params = Vec::new();
    for sp in stack.spans() {
        if stack.span_is_scalable(sp) {
            let rows: Vec<usize> = (sp.bulk_start..sp.end).collect();
            let d0: f64 = rows.iter().map(|&r| films[r].d_nm).sum();
            let fractions: Vec<f64> = rows.iter().map(|&r| films[r].d_nm / d0).collect();
            params.push(ThicknessParam::Span { rows, fractions });
        } else {
            for r in sp.bulk_start..sp.end {
                if films[r].optimize {
                    params.push(ThicknessParam::Row(r));
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
pub(crate) fn apply_thickness_params(
    stack: &mut DesignStack,
    params: &[ThicknessParam],
    x: &[f64],
) -> Result<(), String> {
    for (j, p) in params.iter().enumerate() {
        match p {
            ThicknessParam::Row(i) => stack.set_thickness(*i, x[j])?,
            ThicknessParam::Span { rows, fractions } => {
                for (&r, &f) in rows.iter().zip(fractions.iter()) {
                    stack.set_thickness(r, f * x[j])?;
                }
            }
        }
    }
    Ok(())
}

/// The analytic Jacobian of one `refine_thicknesses` problem: rebuild the
/// stack at `x`, take one sweep with deposits, and multiply the two halves.
///
/// Declines (`Ok(None)`) when the merit spec has rows whose dependence the
/// sensitivity pass does not carry — color demands. That is
/// a property of the spec, not of `x`, so a declining source declines every
/// iteration and the run is simply the finite-difference run. The first
/// decline is remembered, so the later calls skip the deposits sweep whose
/// answer they already know.
pub(crate) struct ThicknessJacobian {
    pub(crate) ctx: SmatrixContext,
    base_stack: DesignStack,
    params: Vec<ThicknessParam>,
    n_wav: usize,
    pub(crate) declined: AtomicBool,
}

impl ThicknessJacobian {
    pub(crate) fn new(
        ctx: SmatrixContext,
        base_stack: DesignStack,
        params: Vec<ThicknessParam>,
        n_wav: usize,
    ) -> Self {
        Self {
            ctx,
            base_stack,
            params,
            n_wav,
            declined: AtomicBool::new(false),
        }
    }
}

impl JacobianSource for ThicknessJacobian {
    fn fill(&self, x: &[f64], jac: &mut Vec<f64>) -> Result<Option<usize>, String> {
        // F2.2: a multi-environment spec is an uncovered spec for now.
        // `curve_sensitivity` walks the targets of ONE simulation and
        // knows nothing about `env_idx`; routing the deposits by design
        // -film name and summing them across environments is F2.3. Until
        // then the joint run is a finite-difference run, which is the
        // documented decline, not a silent single-environment Jacobian.
        // A stack with an incoherent film is covered: its deposits go
        // through the cascade's adjoint weights (W6).
        if self.ctx.envs.is_some() || self.declined.load(Ordering::Relaxed) {
            return Ok(None);
        }
        let mut st = self.base_stack.clone();
        apply_thickness_params(&mut st, &self.params, x)?;
        // The deposits are per FILM ROW (the flat list: each parameter's
        // rows in parameter order), then the mapped assembly contracts
        // row columns into parameter columns with the frozen weights.
        let mut flat_rows = Vec::new();
        for p in &self.params {
            match p {
                ThicknessParam::Row(i) => flat_rows.push(*i),
                ThicknessParam::Span { rows, .. } => flat_rows.extend_from_slice(rows),
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
