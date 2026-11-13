// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::global — a population search over a stage's box (plan §6).
//!
//! A `global` stage hands χ² of its free parameters to one of basin's
//! derivative-free global solvers, then polishes the best point it found
//! with the stage's least squares. The search runs on the box scaled to
//! `[0, 1]ⁿ`, so one step size fits parameters of every unit (CMA-ES has a
//! single σ); a point the model cannot evaluate costs `+∞`-like, which every
//! method rejects. Behind the `opt-basin` feature.
// Without the feature only the refusal is reachable.
#![cfg_attr(not(feature = "opt-basin"), allow(dead_code, unused_imports))]

use serde::{Deserialize, Serialize};

use super::{Chi2Fn, Progress, Watcher};

/// Which basin solver a `global` stage runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GlobalMethod {
    /// Differential evolution, `DE/rand/1/bin`.
    #[default]
    De,
    /// Differential evolution with a Nelder-Mead polish of the best
    /// member every generation (memetic).
    DeInject,
    /// CMA-ES with adaptive-penalty bound handling, started at the
    /// stage's current point.
    CmaEs,
}

/// The search's best point (free order) and the χ² evaluations it spent.
pub(crate) struct Found {
    pub(crate) x: Vec<f64>,
    pub(crate) evals: usize,
}

/// What every method needs: χ² on the box, a start and a budget.
pub(crate) struct Search<'a, 'w> {
    pub(crate) chi2: &'a Chi2Fn,
    pub(crate) lo: Vec<f64>,
    pub(crate) hi: Vec<f64>,
    /// The stage's current point, inside the box.
    pub(crate) x0: Vec<f64>,
    pub(crate) budget: usize,
    pub(crate) seed: u64,
    pub(crate) population: Option<usize>,
    pub(crate) watch: Option<(&'a Watcher<'w>, usize)>,
}

/// Returned for a point the model cannot evaluate: large and finite, so
/// population sorting and CMA-ES's fitness ranking stay well defined.
const REFUSED_CHI2: f64 = 1e100;

impl Search<'_, '_> {
    fn x_of(&self, z: &[f64]) -> Vec<f64> {
        z.iter()
            .zip(self.lo.iter().zip(&self.hi))
            .map(|(&z, (&lo, &hi))| (lo + z * (hi - lo)).clamp(lo, hi))
            .collect()
    }

    pub(crate) fn run(&self, method: GlobalMethod) -> Result<Found, String> {
        if self.budget == 0 {
            return Err("global: budget must be at least 1".into());
        }
        if self
            .lo
            .iter()
            .zip(&self.hi)
            .any(|(lo, hi)| !(lo.is_finite() && hi.is_finite() && lo < hi))
        {
            return Err("global: every free parameter needs a finite box".into());
        }
        imp::run(self, method)
    }
}

#[cfg(not(feature = "opt-basin"))]
mod imp {
    use super::*;

    pub(super) fn run(_: &Search<'_, '_>, _: GlobalMethod) -> Result<Found, String> {
        Err("a global stage is not compiled in: rebuild with `--features opt-basin`".into())
    }
}

#[cfg(feature = "opt-basin")]
mod imp {
    use std::sync::Mutex;
    use std::sync::atomic::{AtomicUsize, Ordering};

    use basin::{
        BasicPopulationState, BoundedCmaEs, BoxConstraints, CmaEsState, CostFunction, De, DeInject,
        DenseMatrix, Executor, NelderMead,
    };

    use super::*;

    /// χ² on the unit box, counted, best tracked and reported.
    struct Unit<'s, 'a, 'w> {
        s: &'s Search<'a, 'w>,
        lower: Vec<f64>,
        upper: Vec<f64>,
        evals: &'s AtomicUsize,
        best: &'s Mutex<(f64, Vec<f64>)>,
    }

    impl CostFunction for Unit<'_, '_, '_> {
        type Param = Vec<f64>;
        type Output = f64;
        type Error = String;

        fn cost(&self, z: &Vec<f64>) -> Result<f64, String> {
            if let Some((w, _)) = self.s.watch {
                w.check()?;
            }
            let x = self.s.x_of(z);
            let chi2 = self.s.chi2.chi2(&x).unwrap_or(REFUSED_CHI2);
            let evals = self.evals.fetch_add(1, Ordering::Relaxed) + 1;
            let mut best = self.best.lock().unwrap_or_else(|e| e.into_inner());
            if chi2 < best.0 {
                *best = (chi2, x);
                if let Some((w, stage)) = self.s.watch {
                    w.tell(&Progress::Improved { stage, evals, chi2 });
                }
            }
            Ok(chi2)
        }
    }

    impl BoxConstraints for Unit<'_, '_, '_> {
        fn lower(&self) -> &Vec<f64> {
            &self.lower
        }

        fn upper(&self) -> &Vec<f64> {
            &self.upper
        }
    }

    pub(super) fn run(s: &Search<'_, '_>, method: GlobalMethod) -> Result<Found, String> {
        let n = s.lo.len();
        let evals = AtomicUsize::new(0);
        let z0: Vec<f64> =
            s.x0.iter()
                .zip(s.lo.iter().zip(&s.hi))
                .map(|(&x, (&lo, &hi))| ((x - lo) / (hi - lo)).clamp(0.0, 1.0))
                .collect();
        // The start counts: a search that finds nothing better keeps it.
        let start = s.chi2.chi2(&s.x0).unwrap_or(REFUSED_CHI2);
        let best = Mutex::new((start, s.x0.clone()));
        let problem = Unit {
            s,
            lower: vec![0.0; n],
            upper: vec![1.0; n],
            evals: &evals,
            best: &best,
        };
        let budget = s.budget as u64;
        let de = || {
            let de = De::new(s.seed);
            match s.population {
                Some(p) => de.with_pop_size(p),
                None => de,
            }
        };
        match method {
            GlobalMethod::De => Executor::new(
                problem,
                de(),
                BasicPopulationState::<Vec<f64>>::with_size(1),
            )
            .max_iter(u64::MAX)
            .max_cost_evals(budget)
            .run()
            .map(|_| ()),
            GlobalMethod::DeInject => Executor::new(
                problem,
                DeInject::with_inner_solver(de(), NelderMead::adaptive())
                    .with_k(1)
                    .with_inner_max_iter(20),
                BasicPopulationState::<Vec<f64>>::with_size(1),
            )
            .max_iter(u64::MAX)
            .max_cost_evals(budget)
            .run()
            .map(|_| ()),
            GlobalMethod::CmaEs => {
                let solver = BoundedCmaEs::<Vec<f64>, DenseMatrix>::new(s.seed);
                let solver = match s.population {
                    Some(p) => solver.with_lambda(p),
                    None => solver,
                };
                Executor::new(
                    problem,
                    solver,
                    CmaEsState::<Vec<f64>, DenseMatrix>::new(z0, 0.3),
                )
                .max_iter(u64::MAX)
                .max_cost_evals(budget)
                .run()
                .map(|_| ())
            }
        }
        .map_err(|e| format!("global {method:?}: {e}"))?;
        let (chi2, x) = best.into_inner().unwrap_or_else(|e| e.into_inner());
        if !(chi2 < REFUSED_CHI2) {
            return Err("global: no point in the box could be evaluated".into());
        }
        Ok(Found {
            x,
            evals: evals.into_inner(),
        })
    }
}
