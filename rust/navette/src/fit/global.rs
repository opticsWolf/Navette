// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::global — a population search over a stage's box (plan §3.8.3).
//!
//! A `global` stage hands χ² of its free parameters to a derivative-free
//! global solver, then polishes the best point it found with the stage's
//! least squares. The search runs on the box scaled to `[0, 1]ⁿ`, so one
//! step size fits parameters of every unit (CMA-ES has a single σ); a
//! point the model cannot evaluate costs a huge finite χ², which every
//! method rejects.
//!
//! `de` and `de_best` are Navette's own differential evolution
//! ([`crate::lsq::de`]), in every build. `basin_de`, `de_inject` and
//! `cma_es` run basin's solvers and need the `opt-basin` feature; they
//! stay as references for the native DE.
// Without the feature the basin arms only refuse.
#![cfg_attr(not(feature = "opt-basin"), allow(dead_code, unused_imports))]

use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use serde::{Deserialize, Serialize};

use super::{Chi2Fn, Progress, Watcher};
use crate::lsq::de::{DeConfig, DeMutation, minimize};

/// Which solver a `global` stage runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GlobalMethod {
    /// Navette's differential evolution, `DE/rand/1/bin`.
    #[default]
    De,
    /// Navette's differential evolution, `DE/best/1/bin`: greedier.
    DeBest,
    /// basin's `DE/rand/1/bin` (`opt-basin`), the reference for `de`.
    BasinDe,
    /// basin's memetic DE: a Nelder-Mead polish of the best member every
    /// generation (`opt-basin`).
    DeInject,
    /// basin's CMA-ES with adaptive-penalty bound handling, started at the
    /// stage's current point (`opt-basin`).
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

/// χ² on the unit box, counted, best tracked and reported.
struct Tracker<'s, 'a, 'w> {
    s: &'s Search<'a, 'w>,
    evals: AtomicUsize,
    best: Mutex<(f64, Vec<f64>)>,
}

impl<'s, 'a, 'w> Tracker<'s, 'a, 'w> {
    fn new(s: &'s Search<'a, 'w>) -> Self {
        // The start counts: a search that finds nothing better keeps it.
        let start = s.chi2.chi2(&s.x0).unwrap_or(REFUSED_CHI2);
        Self {
            s,
            evals: AtomicUsize::new(0),
            best: Mutex::new((start, s.x0.clone())),
        }
    }

    fn cost(&self, z: &[f64]) -> Result<f64, String> {
        if let Some((w, _)) = self.s.watch {
            w.check()?;
        }
        let x = self.s.x_of(z);
        let chi2 = self.s.chi2.chi2(&x).unwrap_or(REFUSED_CHI2);
        let evals = self.evals.fetch_add(1, Ordering::Relaxed) + 1;
        let mut best = self.best.lock().unwrap_or_else(|e| e.into_inner());
        // A tie keeps the lexicographically smaller point, so the order in
        // which parallel evaluations arrive cannot change the result.
        if chi2 < best.0 {
            *best = (chi2, x);
            if let Some((w, stage)) = self.s.watch {
                w.tell(&Progress::Searched { stage, evals, chi2 });
            }
        } else if chi2 == best.0 && x < best.1 {
            best.1 = x;
        }
        Ok(chi2)
    }

    fn found(self) -> Result<Found, String> {
        let (chi2, x) = self.best.into_inner().unwrap_or_else(|e| e.into_inner());
        if !(chi2 < REFUSED_CHI2) {
            return Err("global: no point in the box could be evaluated".into());
        }
        Ok(Found {
            x,
            evals: self.evals.into_inner(),
        })
    }
}

impl Search<'_, '_> {
    fn x_of(&self, z: &[f64]) -> Vec<f64> {
        z.iter()
            .zip(self.lo.iter().zip(&self.hi))
            .map(|(&z, (&lo, &hi))| (lo + z * (hi - lo)).clamp(lo, hi))
            .collect()
    }

    fn z0(&self) -> Vec<f64> {
        self.x0
            .iter()
            .zip(self.lo.iter().zip(&self.hi))
            .map(|(&x, (&lo, &hi))| ((x - lo) / (hi - lo)).clamp(0.0, 1.0))
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
        let t = Tracker::new(self);
        match method {
            GlobalMethod::De | GlobalMethod::DeBest => {
                let n = self.lo.len();
                let cfg = DeConfig {
                    mutation: if method == GlobalMethod::De {
                        DeMutation::Rand1
                    } else {
                        DeMutation::Best1
                    },
                    pop_size: self.population,
                    budget: self.budget,
                    seed: self.seed,
                    ..DeConfig::default()
                };
                minimize(
                    |z| t.cost(z),
                    &vec![0.0; n],
                    &vec![1.0; n],
                    Some(&self.z0()),
                    &cfg,
                )
                .map_err(|e| format!("global {method:?}: {e}"))?;
            }
            _ => basin_run(self, &t, method)?,
        }
        t.found()
    }
}

#[cfg(not(feature = "opt-basin"))]
fn basin_run(_: &Search<'_, '_>, _: &Tracker<'_, '_, '_>, m: GlobalMethod) -> Result<(), String> {
    Err(format!(
        "global {m:?}: basin's solvers are not compiled in: rebuild with `--features opt-basin`, or use `de`"
    ))
}

#[cfg(feature = "opt-basin")]
fn basin_run(
    s: &Search<'_, '_>,
    t: &Tracker<'_, '_, '_>,
    method: GlobalMethod,
) -> Result<(), String> {
    use basin::{
        BasicPopulationState, BoundedCmaEs, BoxConstraints, CmaEsState, CostFunction, De, DeInject,
        DenseMatrix, Executor, NelderMead,
    };

    struct Unit<'t, 's, 'a, 'w> {
        t: &'t Tracker<'s, 'a, 'w>,
        lower: Vec<f64>,
        upper: Vec<f64>,
    }

    impl CostFunction for Unit<'_, '_, '_, '_> {
        type Param = Vec<f64>;
        type Output = f64;
        type Error = String;

        fn cost(&self, z: &Vec<f64>) -> Result<f64, String> {
            self.t.cost(z)
        }
    }

    impl BoxConstraints for Unit<'_, '_, '_, '_> {
        fn lower(&self) -> &Vec<f64> {
            &self.lower
        }

        fn upper(&self) -> &Vec<f64> {
            &self.upper
        }
    }

    let n = s.lo.len();
    let problem = Unit {
        t,
        lower: vec![0.0; n],
        upper: vec![1.0; n],
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
        GlobalMethod::De | GlobalMethod::DeBest => unreachable!("native"),
        GlobalMethod::BasinDe => Executor::new(
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
                CmaEsState::<Vec<f64>, DenseMatrix>::new(s.z0(), 0.3),
            )
            .max_iter(u64::MAX)
            .max_cost_evals(budget)
            .run()
            .map(|_| ())
        }
    }
    .map_err(|e| format!("global {method:?}: {e}"))
}
