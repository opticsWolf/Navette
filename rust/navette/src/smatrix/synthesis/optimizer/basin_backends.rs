// SPDX-License-Identifier: LGPL-3.0-or-later
//! Basin adapters. See docs/basin-backends.md for tolerance and budget semantics.

use super::{IntervalMap, MappedJacobian, OptimizerBackend, OptimizerResult};
use crate::smatrix::synthesis::thick_opt::{
    JacobianSource, LmConfig, LmTermination, build_jacobian,
};
use basin::{
    BoundedFiniteDiff, BoxConstraints, CostFunction, DenseMatrix, Executor, Jacobian,
    LevenbergMarquardtQr, Method, NativeConvergenceTest, Residual, TerminationReason,
    TrustRegionReflective,
};
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

#[derive(Default)]
struct Counts {
    residuals: AtomicUsize,
    analytic: AtomicUsize,
    rows: AtomicUsize,
}

#[derive(Clone)]
struct ResidualProblem<'a> {
    residuals: &'a (dyn Fn(&[f64], &mut Vec<f64>) -> Result<(), String> + Sync),
    counts: Arc<Counts>,
}

impl Residual for ResidualProblem<'_> {
    type Param = Vec<f64>;
    type Output = Vec<f64>;
    type Error = String;

    fn residual(&self, x: &Vec<f64>) -> Result<Vec<f64>, String> {
        self.counts.residuals.fetch_add(1, Ordering::Relaxed);
        let mut r = Vec::new();
        (self.residuals)(x, &mut r)?;
        if r.is_empty() {
            return Err("empty residual vector".into());
        }
        let rows =
            self.counts
                .rows
                .compare_exchange(0, r.len(), Ordering::Relaxed, Ordering::Relaxed);
        if rows.is_err_and(|rows| rows != r.len()) {
            return Err("inconsistent residual length".into());
        }
        // Basin handles non-finite trial values by rejecting the step. Do not
        // turn these into callback errors, which mean the run must abort.
        Ok(r)
    }
}

impl CostFunction for ResidualProblem<'_> {
    type Param = Vec<f64>;
    type Output = f64;
    type Error = String;

    fn cost(&self, x: &Vec<f64>) -> Result<f64, String> {
        Ok(0.5 * self.residual(x)?.iter().map(|r| r * r).sum::<f64>())
    }
}

struct Problem<'a, J: ?Sized> {
    raw: ResidualProblem<'a>,
    jacobian: Option<&'a J>,
    bounds: Option<(Vec<f64>, Vec<f64>)>,
}

impl<J: ?Sized> Residual for Problem<'_, J> {
    type Param = Vec<f64>;
    type Output = Vec<f64>;
    type Error = String;

    fn residual(&self, x: &Vec<f64>) -> Result<Vec<f64>, String> {
        self.raw.residual(x)
    }
}

impl<J: ?Sized> CostFunction for Problem<'_, J> {
    type Param = Vec<f64>;
    type Output = f64;
    type Error = String;

    fn cost(&self, x: &Vec<f64>) -> Result<f64, String> {
        self.raw.cost(x)
    }
}

impl<J: JacobianSource + ?Sized> Jacobian for Problem<'_, J> {
    type Jacobian = DenseMatrix;

    fn jacobian(&self, x: &Vec<f64>) -> Result<DenseMatrix, String> {
        let mut values = Vec::new();
        if let Some(source) = self.jacobian
            && let Some(rows) = source.fill(x, &mut values)?
        {
            if rows != self.raw.counts.rows.load(Ordering::Relaxed)
                || rows.checked_mul(x.len()) != Some(values.len())
            {
                return Err("inconsistent Jacobian shape".into());
            }
            self.raw.counts.analytic.fetch_add(1, Ordering::Relaxed);
            return Ok(DenseMatrix::from_row_slice(rows, x.len(), &values));
        }

        if let Some((lower, upper)) = &self.bounds {
            // Ordinary central differences can leave the physical domain at
            // an active bound. This adapter selects feasible one-sided probes.
            BoundedFiniteDiff::new(self.raw.clone(), lower.clone(), upper.clone())
                .jacobian_method(Method::Central)
                .jacobian(x)
        } else {
            let residuals = |x: &[f64], out: &mut Vec<f64>| {
                *out = self.raw.residual(&x.to_vec())?;
                Ok(())
            };
            build_jacobian(&residuals, x, &mut values)?;
            Ok(DenseMatrix::from_row_slice(
                values.len() / x.len(),
                x.len(),
                &values,
            ))
        }
    }
}

impl<J: ?Sized> BoxConstraints for Problem<'_, J> {
    fn lower(&self) -> &Vec<f64> {
        &self.bounds.as_ref().expect("TRF supplies bounds").0
    }

    fn upper(&self) -> &Vec<f64> {
        &self.bounds.as_ref().expect("TRF supplies bounds").1
    }
}

fn map_termination(
    reason: TerminationReason,
    native_tests: &[NativeConvergenceTest],
) -> Result<LmTermination, String> {
    match reason {
        TerminationReason::SolverConverged => {
            // Basin can report several tests in an unspecified order. Prefer
            // gradient, then cost, then step when Navette needs one reason.
            if native_tests.iter().any(|test| matches!(test,
                NativeConvergenceTest::AbsoluteGradient
                | NativeConvergenceTest::GradientOrthogonality
                | NativeConvergenceTest::RobustGradientOrthogonality
                | NativeConvergenceTest::AbsoluteScaledGradient
            )) {
                Ok(LmTermination::Gradient)
            } else if native_tests.contains(&NativeConvergenceTest::RelativeModelReduction) {
                Ok(LmTermination::Cost)
            } else if native_tests.iter().any(|test| matches!(test,
                NativeConvergenceTest::RelativeTrialStep
                | NativeConvergenceTest::RelativeTrustRadius
            )) {
                Ok(LmTermination::Step)
            } else {
                Ok(LmTermination::Converged)
            }
        }
        TerminationReason::GradientTolerance
        | TerminationReason::RelativeGradientTolerance
        | TerminationReason::ProjectedGradientTolerance => Ok(LmTermination::Gradient),
        TerminationReason::CostTolerance | TerminationReason::RelativeCostTolerance => Ok(LmTermination::Cost),
        TerminationReason::ParamTolerance | TerminationReason::RelativeParamTolerance => Ok(LmTermination::Step),
        TerminationReason::MaxIter | TerminationReason::MaxEvaluations => Ok(LmTermination::MaxIterations),
        TerminationReason::NumericalNoProgress => Ok(LmTermination::Stalled),
        TerminationReason::SolverFailed => Err("solver failed (invalid initial evaluations, non-finite derivatives, or exhausted numerical progress)".into()),
        other => Err(format!("unexpected termination reason: {other:?}")),
    }
}

fn validate(x0: &[f64], lb: &[f64], ub: &[f64], cfg: &LmConfig) -> Result<(), String> {
    if x0.is_empty() {
        return Err("empty parameter vector".into());
    }
    if lb.len() != x0.len() || ub.len() != x0.len() {
        return Err("bound length mismatch".into());
    }
    for (i, ((&x, &lower), &upper)) in x0.iter().zip(lb).zip(ub).enumerate() {
        if !x.is_finite() || !lower.is_finite() || !upper.is_finite() || lower > upper {
            return Err(format!(
                "invalid parameter or finite bounds at coordinate {i}"
            ));
        }
        if cfg.backend == OptimizerBackend::BasinTrf && (x < lower || x > upper) {
            return Err(format!("initial parameter {i} is outside its bounds"));
        }
    }
    for (name, value) in [("ftol", cfg.ftol), ("xtol", cfg.xtol), ("gtol", cfg.gtol)] {
        if !value.is_finite() || value < 0.0 {
            return Err(format!("{name} must be finite and nonnegative"));
        }
    }
    if cfg.backend == OptimizerBackend::BasinLmQr
        && (!cfg.lambda_init.is_finite() || cfg.lambda_init <= 0.0)
    {
        return Err("lambda_init must be finite and positive".into());
    }
    Ok(())
}

pub(super) fn run<F, J>(
    residuals: &F,
    jacobian: Option<&J>,
    x0: &[f64],
    lb: &[f64],
    ub: &[f64],
    cfg: &LmConfig,
) -> Result<OptimizerResult, String>
where
    F: Fn(&[f64], &mut Vec<f64>) -> Result<(), String> + Sync,
    J: JacobianSource + ?Sized,
{
    let execute = || -> Result<OptimizerResult, String> {
        validate(x0, lb, ub, cfg)?;
        let counts = Arc::new(Counts::default());
        let budget_counts = Arc::clone(&counts);
        let max_evals = cfg.max_evals;
        // Basin's residual counter cannot see finite-difference probes inside
        // our Jacobian adapter. Count them at the callback boundary instead.
        let budget = move |_: &basin::NllsState<Vec<f64>>| {
            (budget_counts.residuals.load(Ordering::Relaxed) >= max_evals)
                .then_some(TerminationReason::MaxEvaluations)
        };

        let (result, x, termination) = if cfg.backend == OptimizerBackend::BasinLmQr {
            let map = IntervalMap::new(lb, ub)?;
            let mapped_r = |u: &[f64], out: &mut Vec<f64>| residuals(&map.to_bounded(u), out);
            let mapped_j = jacobian.map(|inner| MappedJacobian { inner, map: &map });
            let problem = Problem {
                raw: ResidualProblem {
                    residuals: &mapped_r,
                    counts: Arc::clone(&counts),
                },
                jacobian: mapped_j.as_ref(),
                bounds: None,
            };
            let solver = LevenbergMarquardtQr::<Vec<f64>, DenseMatrix>::new()
                .with_absolute_gradient_tolerance(cfg.gtol)
                .with_gradient_orthogonality_tolerance(cfg.gtol_scale_invariant.then_some(cfg.gtol))
                .with_relative_model_reduction_tolerance(cfg.ftol)
                .with_relative_step_tolerance(cfg.xtol)
                .with_tau(cfg.lambda_init);
            let result = Executor::from_start(problem, solver, map.to_unbounded(x0))
                .max_iter(cfg.max_iterations as u64)
                .stop_when(budget)
                .run_with_solver()?;
            let x = map.to_bounded(result.param());
            let termination = map_termination(result.reason, result.native_convergence_tests())?;
            (result.into_result(), x, termination)
        } else {
            let problem = Problem {
                raw: ResidualProblem {
                    residuals,
                    counts: Arc::clone(&counts),
                },
                jacobian,
                bounds: Some((lb.to_vec(), ub.to_vec())),
            };
            let solver = TrustRegionReflective::new()
                .with_absolute_scaled_gradient_tolerance(cfg.gtol)
                .with_relative_cost_change_tolerance(cfg.ftol)
                .with_relative_step_tolerance(cfg.xtol);
            let result = Executor::from_start(problem, solver, x0.to_vec())
                .max_iter(cfg.max_iterations as u64)
                .stop_when(budget)
                .run_with_solver()?;
            let x = result.param().clone();
            let termination = map_termination(result.reason, result.native_convergence_tests())?;
            (result.into_result(), x, termination)
        };
        let cost = 2.0 * result.cost();
        if !cost.is_finite() || x.iter().any(|v| !v.is_finite()) {
            return Err("non-finite solver result".into());
        }
        Ok(OptimizerResult {
            x,
            cost,
            iterations: result.iter() as usize,
            evals: counts.residuals.load(Ordering::Relaxed),
            termination,
            gain_ratio: f64::NAN,
            analytic_jacobians: counts.analytic.load(Ordering::Relaxed),
            backend: cfg.backend,
        })
    };
    execute().map_err(|e| format!("{}: {e}", cfg.backend.as_str()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn termination_reasons_preserve_convergence_and_failure() {
        for (reason, expected) in [
            (TerminationReason::SolverConverged, LmTermination::Converged),
            (
                TerminationReason::RelativeCostTolerance,
                LmTermination::Cost,
            ),
            (
                TerminationReason::RelativeParamTolerance,
                LmTermination::Step,
            ),
            (
                TerminationReason::MaxEvaluations,
                LmTermination::MaxIterations,
            ),
            (
                TerminationReason::NumericalNoProgress,
                LmTermination::Stalled,
            ),
        ] {
            assert_eq!(map_termination(reason, &[]).unwrap(), expected);
        }
        assert!(map_termination(TerminationReason::SolverFailed, &[]).is_err());
        assert!(map_termination(TerminationReason::UserRequested, &[]).is_err());
    }

    #[test]
    fn native_tests_identify_convergence_with_stable_precedence() {
        use NativeConvergenceTest::*;
        for (tests, expected) in [
            (vec![AbsoluteGradient], LmTermination::Gradient),
            (vec![GradientOrthogonality], LmTermination::Gradient),
            (vec![RobustGradientOrthogonality], LmTermination::Gradient),
            (vec![AbsoluteScaledGradient], LmTermination::Gradient),
            (vec![RelativeModelReduction], LmTermination::Cost),
            (vec![RelativeTrialStep], LmTermination::Step),
            (vec![RelativeTrustRadius], LmTermination::Step),
            (vec![NoFreeParameters], LmTermination::Converged),
            (
                vec![RelativeTrialStep, RelativeModelReduction],
                LmTermination::Cost,
            ),
            (
                vec![RelativeModelReduction, RelativeTrialStep],
                LmTermination::Cost,
            ),
            (
                vec![RelativeModelReduction, AbsoluteGradient],
                LmTermination::Gradient,
            ),
        ] {
            assert_eq!(
                map_termination(TerminationReason::SolverConverged, &tests).unwrap(),
                expected,
                "{tests:?}"
            );
        }
    }

    #[test]
    fn native_diagnostics_do_not_override_other_termination_reasons() {
        let tests = &[NativeConvergenceTest::AbsoluteGradient];
        for (reason, expected) in [
            (
                TerminationReason::NumericalNoProgress,
                LmTermination::Stalled,
            ),
            (
                TerminationReason::MaxEvaluations,
                LmTermination::MaxIterations,
            ),
            (
                TerminationReason::RelativeCostTolerance,
                LmTermination::Cost,
            ),
            (
                TerminationReason::RelativeParamTolerance,
                LmTermination::Step,
            ),
        ] {
            assert_eq!(map_termination(reason, tests).unwrap(), expected);
        }
        assert!(map_termination(TerminationReason::SolverFailed, tests).is_err());
    }
}
