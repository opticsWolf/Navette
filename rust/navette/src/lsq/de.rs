// SPDX-License-Identifier: LGPL-3.0-or-later
//! Differential evolution over a finite box, in the default build.
//!
//! Synchronous `DE/rand/1/bin` or `DE/best/1/bin` (Storn & Price): every
//! trial of a generation is built from the frozen current population, the
//! generation is evaluated in parallel and collected in order, and a trial
//! replaces its target when it is no worse (`<=`). The mutation scale `F`
//! is drawn once per generation from `dither` (SciPy's dithering). A donor
//! coordinate outside the box, or not finite, is replaced by a uniform draw
//! inside its bounds (reinit-per-coordinate repair, basin's and DEoptim's
//! choice: clipping would pile donors up on the faces).
//!
//! Deterministic: the RNG is consumed on one thread in a fixed order and
//! the costs are collected in population order, so the same seed gives the
//! same bits whatever the thread count.
//!
//! The population starts with `x0` (when given) and uniform draws. A run
//! stops after the generation that would overrun `budget`, or once the
//! population has collapsed onto one point.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

/// How a donor is built.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeMutation {
    /// `x[r0] + F (x[r1] - x[r2])`: explores, the robust default.
    Rand1,
    /// `best + F (x[r0] - x[r1])`: converges faster, more greedy.
    Best1,
}

/// Settings of one run.
#[derive(Debug, Clone)]
pub struct DeConfig {
    pub mutation: DeMutation,
    /// Population size; `None` is `max(10 n, 8)`.
    pub pop_size: Option<usize>,
    /// `F` is drawn uniformly from `[dither.0, dither.1)` each generation;
    /// equal ends fix it.
    pub dither: (f64, f64),
    /// Crossover probability.
    pub cr: f64,
    /// Most cost evaluations; the initial population counts.
    pub budget: usize,
    pub seed: u64,
    /// The population has collapsed when every coordinate's spread is at
    /// most `xtol` times its box width.
    pub xtol: f64,
}

impl Default for DeConfig {
    fn default() -> Self {
        Self {
            mutation: DeMutation::Rand1,
            pop_size: None,
            dither: (0.5, 1.0),
            cr: 0.9,
            budget: 10_000,
            seed: 0,
            xtol: 1e-10,
        }
    }
}

/// The best member found.
#[derive(Debug, Clone)]
pub struct DeResult {
    pub x: Vec<f64>,
    pub cost: f64,
    pub evals: usize,
    pub generations: usize,
}

/// Minimise `cost` over the box `[lo, hi]`. A cost error stops the run
/// and is returned (a watcher asking to stop, typically).
pub fn minimize<F>(
    cost: F,
    lo: &[f64],
    hi: &[f64],
    x0: Option<&[f64]>,
    cfg: &DeConfig,
) -> Result<DeResult, String>
where
    F: Fn(&[f64]) -> Result<f64, String> + Sync,
{
    let n = lo.len();
    if n == 0 || hi.len() != n {
        return Err("de: the box needs at least one coordinate and matching ends".into());
    }
    if lo
        .iter()
        .zip(hi)
        .any(|(l, h)| !(l.is_finite() && h.is_finite() && l < h))
    {
        return Err("de: every coordinate needs a finite box with lo < hi".into());
    }
    let np = cfg.pop_size.unwrap_or((10 * n).max(8));
    if np < 4 {
        return Err(format!(
            "de: a population of {np} is too small (at least 4)"
        ));
    }
    let (fmin, fmax) = cfg.dither;
    if !(fmin > 0.0 && fmin <= fmax && fmax.is_finite()) {
        return Err(format!("de: dither ({fmin}, {fmax}) is not 0 < min <= max"));
    }
    if !(0.0..=1.0).contains(&cfg.cr) {
        return Err(format!(
            "de: crossover probability {} is outside [0, 1]",
            cfg.cr
        ));
    }
    let mut rng = StdRng::seed_from_u64(cfg.seed);
    let mut pop: Vec<Vec<f64>> = (0..np)
        .map(|_| (0..n).map(|j| rng.random_range(lo[j]..=hi[j])).collect())
        .collect();
    if let Some(x0) = x0 {
        if x0.len() != n {
            return Err("de: the start has the wrong length".into());
        }
        pop[0] = x0
            .iter()
            .zip(lo.iter().zip(hi))
            .map(|(&x, (&l, &h))| x.clamp(l, h))
            .collect();
    }
    let evaluate = |xs: &[Vec<f64>]| -> Result<Vec<f64>, String> {
        xs.par_iter()
            .map(|x| cost(x).map(|c| if c.is_nan() { f64::INFINITY } else { c }))
            .collect()
    };
    let mut costs = evaluate(&pop)?;
    let mut evals = np;
    let mut generations = 0;
    let best_of = |costs: &[f64]| {
        (0..costs.len())
            .min_by(|&a, &b| costs[a].total_cmp(&costs[b]))
            .unwrap()
    };
    while evals + np <= cfg.budget && !collapsed(&pop, lo, hi, cfg.xtol) {
        let f = if fmin < fmax {
            rng.random_range(fmin..fmax)
        } else {
            fmin
        };
        let best = best_of(&costs);
        let trials: Vec<Vec<f64>> = (0..np)
            .map(|i| {
                let r = distinct(&mut rng, np, i);
                let base = match cfg.mutation {
                    DeMutation::Rand1 => &pop[r[0]],
                    DeMutation::Best1 => &pop[best],
                };
                let (a, b) = match cfg.mutation {
                    DeMutation::Rand1 => (&pop[r[1]], &pop[r[2]]),
                    DeMutation::Best1 => (&pop[r[0]], &pop[r[1]]),
                };
                let j_rand = rng.random_range(0..n);
                (0..n)
                    .map(|j| {
                        let mut d = base[j] + f * (a[j] - b[j]);
                        if !(d.is_finite() && lo[j] <= d && d <= hi[j]) {
                            d = rng.random_range(lo[j]..=hi[j]);
                        }
                        if j == j_rand || rng.random::<f64>() < cfg.cr {
                            d
                        } else {
                            pop[i][j]
                        }
                    })
                    .collect()
            })
            .collect();
        let trial_costs = evaluate(&trials)?;
        evals += np;
        generations += 1;
        for (i, (t, c)) in trials.into_iter().zip(trial_costs).enumerate() {
            if c <= costs[i] {
                pop[i] = t;
                costs[i] = c;
            }
        }
    }
    let best = best_of(&costs);
    Ok(DeResult {
        x: pop.swap_remove(best),
        cost: costs[best],
        evals,
        generations,
    })
}

/// Three indices, pairwise distinct and all different from `i`.
fn distinct(rng: &mut StdRng, np: usize, i: usize) -> [usize; 3] {
    let mut r = [usize::MAX; 3];
    for k in 0..3 {
        r[k] = loop {
            let c = rng.random_range(0..np);
            if c != i && !r[..k].contains(&c) {
                break c;
            }
        };
    }
    r
}

fn collapsed(pop: &[Vec<f64>], lo: &[f64], hi: &[f64], xtol: f64) -> bool {
    (0..lo.len()).all(|j| {
        let (mn, mx) = pop
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(mn, mx), x| {
                (mn.min(x[j]), mx.max(x[j]))
            });
        mx - mn <= xtol * (hi[j] - lo[j])
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rastrigin(x: &[f64]) -> Result<f64, String> {
        Ok(10.0 * x.len() as f64
            + x.iter()
                .map(|v| v * v - 10.0 * (2.0 * std::f64::consts::PI * v).cos())
                .sum::<f64>())
    }

    fn rosenbrock(x: &[f64]) -> Result<f64, String> {
        Ok(x.windows(2)
            .map(|w| 100.0 * (w[1] - w[0] * w[0]).powi(2) + (1.0 - w[0]).powi(2))
            .sum())
    }

    #[test]
    fn finds_the_global_minimum_of_rastrigin() {
        let (lo, hi) = (vec![-5.12; 4], vec![5.12; 4]);
        for seed in 0..5 {
            let cfg = DeConfig {
                budget: 40_000,
                seed,
                ..DeConfig::default()
            };
            let r = minimize(rastrigin, &lo, &hi, None, &cfg).unwrap();
            assert!(r.cost < 1e-6, "seed {seed}: {}", r.cost);
            assert!(r.evals <= cfg.budget);
        }
    }

    #[test]
    fn best1_solves_rosenbrock() {
        let (lo, hi) = (vec![-2.0; 3], vec![2.0; 3]);
        let cfg = DeConfig {
            mutation: DeMutation::Best1,
            budget: 30_000,
            ..DeConfig::default()
        };
        let r = minimize(rosenbrock, &lo, &hi, None, &cfg).unwrap();
        assert!(r.cost < 1e-8, "{}", r.cost);
    }

    #[test]
    fn same_seed_same_bits_whatever_the_threads() {
        let (lo, hi) = (vec![-5.12; 3], vec![5.12; 3]);
        let cfg = DeConfig {
            budget: 3000,
            seed: 9,
            ..DeConfig::default()
        };
        let a = minimize(rastrigin, &lo, &hi, None, &cfg).unwrap();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let b = pool.install(|| minimize(rastrigin, &lo, &hi, None, &cfg).unwrap());
        assert_eq!(a.cost.to_bits(), b.cost.to_bits());
        assert!(
            a.x.iter()
                .zip(&b.x)
                .all(|(p, q)| p.to_bits() == q.to_bits())
        );
    }

    #[test]
    fn the_start_is_kept_when_nothing_beats_it() {
        let (lo, hi) = (vec![-1.0; 2], vec![1.0; 2]);
        let cfg = DeConfig {
            budget: 20,
            pop_size: Some(20),
            ..DeConfig::default()
        };
        let r = minimize(
            |x| Ok(x[0] * x[0] + x[1] * x[1]),
            &lo,
            &hi,
            Some(&[0.0, 0.0]),
            &cfg,
        )
        .unwrap();
        assert_eq!((r.cost, r.evals, r.generations), (0.0, 20, 0));
    }

    #[test]
    fn bad_settings_name_themselves() {
        let ok = DeConfig::default();
        let f = |_: &[f64]| Ok(0.0);
        assert!(
            minimize(f, &[1.0], &[0.0], None, &ok)
                .unwrap_err()
                .contains("lo < hi")
        );
        let small = DeConfig {
            pop_size: Some(3),
            ..ok.clone()
        };
        assert!(
            minimize(f, &[0.0], &[1.0], None, &small)
                .unwrap_err()
                .contains("too small")
        );
        let cr = DeConfig { cr: 1.5, ..ok };
        assert!(
            minimize(f, &[0.0], &[1.0], None, &cr)
                .unwrap_err()
                .contains("crossover")
        );
    }

    #[test]
    fn a_cost_error_stops_the_run() {
        let e = minimize(
            |_| Err("stopped".to_string()),
            &[0.0],
            &[1.0],
            None,
            &DeConfig::default(),
        )
        .unwrap_err();
        assert_eq!(e, "stopped");
    }
}
