// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::explore — which strategy works best (plan §3.8.1, §3.8.2).
//!
//! - [`generate`] turns parameter sets into strategies ([`StrategyGen`]):
//!   every ordering of the sets, regroupings of them into blocks, or
//!   cumulative release. Each generated strategy is named after its
//!   arrangement (`thicknesses > groups | polish`), so a ranking reads as
//!   "which order works".
//! - [`tournament`] runs strategies on one problem (measured data, no
//!   truth) and ranks them by final χ².
//! - [`draw_truth`] makes a known truth from the model's own error laws,
//!   and [`benchmark`] runs every strategy on every drawn truth and ranks
//!   them by how well they recover it.
//!
//! **Truths live in parameter space.** `expand` with errors on perturbs
//! every traversal entry, including the index of each layer separately,
//! and most of that has no address: a fit could never recover it, so the
//! recovery score would measure the model's blind spots, not the
//! strategy. [`draw_truth`] draws one value per declared free parameter
//! instead, from the error law its group carries (`Group::error_mask`
//! switches each channel on):
//!
//! | Parameter | Draw |
//! |---|---|
//! | `layer[..].thickness`, `.roughness`, `.interface`, `.inh_delta` | that channel's law on the value, as `expand` draws it |
//! | `group[G].n_factor`, `.k_factor`, `.thick_factor` | the n, k or thickness law's *relative* part, as a factor |
//! | `group[G].*_summand` | that channel's law's *absolute* part, added |
//! | anything else (material parameters, gradients) | held at its start |
//!
//! The draw is clamped to the parameter's bounds, ties follow, and the
//! problem's datasets are simulated on the truth with seeded Gaussian
//! noise at each point's σ. Everything is deterministic per seed; the
//! wall-clock times reported are the only thing that is not.

use std::time::{Duration, Instant};

use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand_distr::{Distribution, Normal};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::structure::enums::ErrorMask;
use crate::structure::{ErrorParams, Group};

use super::model::FitModel;
use super::params::{GroupField, LayerField, ParamAddr, default_bounds};
use super::sets::ParamSet;
use super::strategy::{Stage, Strategy, StrategyReport, run_strategy};
use super::{FitProblem, evaluate};

/// Orderings run in full up to this many; above it they are sampled.
pub const DEFAULT_ORDERING_CAP: usize = 24;

/// How to make strategies (plan §3.8.1, decision 6). Set names are
/// parameter sets, user or built-in (`thicknesses`, `groups`, ...).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum StrategyGen {
    /// Hand-written strategies, as they are.
    Explicit(Vec<Strategy>),
    /// Every order of `sets`, one stage per set, then a stage with all of
    /// them free when `polish`. `sample: (count, seed)` takes a seeded
    /// sample of the orders instead.
    Orderings {
        sets: Vec<String>,
        #[serde(default)]
        sample: Option<(usize, u64)>,
        #[serde(default = "yes")]
        polish: bool,
    },
    /// One strategy per variant; a variant is a list of blocks, a block
    /// the sets freed together: `[[["dispersion_H", "dispersion_L"],
    /// ["thicknesses"]], [["dispersion_H"], ["dispersion_L"],
    /// ["thicknesses"]]]`.
    Groupings {
        variants: Vec<Vec<Vec<String>>>,
        #[serde(default = "yes")]
        polish: bool,
    },
    /// In each order of `sets`: free the first, then the first two, and
    /// so on; the last stage frees all.
    Cumulative {
        sets: Vec<String>,
        #[serde(default)]
        sample: Option<(usize, u64)>,
    },
}

fn yes() -> bool {
    true
}

/// Strategies from [`generate`], with what was sampled.
#[derive(Debug, Clone, PartialEq)]
pub struct Generated {
    pub strategies: Vec<Strategy>,
    /// One line per generator that sampled its orders.
    pub notes: Vec<String>,
}

/// Expand `gens` into strategies. Orderings and cumulative generators
/// with more than `cap` orders and no `sample` take a seeded (seed 0)
/// sample of `cap`, and say so in `notes`.
pub fn generate(gens: &[StrategyGen], cap: usize) -> Result<Generated, String> {
    let mut strategies = Vec::new();
    let mut notes = Vec::new();
    for g in gens {
        match g {
            StrategyGen::Explicit(list) => strategies.extend(list.iter().cloned()),
            StrategyGen::Orderings {
                sets,
                sample,
                polish,
            } => {
                for order in orders(sets, *sample, cap, "orderings", &mut notes)? {
                    let mut stages: Vec<Stage> = order
                        .iter()
                        .map(|s| stage(s, std::slice::from_ref(s)))
                        .collect();
                    let mut name = order.join(" > ");
                    if *polish {
                        stages.push(stage("polish", sets));
                        name.push_str(" | polish");
                    }
                    strategies.push(Strategy { name, stages });
                }
            }
            StrategyGen::Groupings { variants, polish } => {
                for variant in variants {
                    if variant.is_empty() || variant.iter().any(|b| b.is_empty()) {
                        return Err("groupings: an empty variant or block".into());
                    }
                    let names: Vec<String> = variant.iter().map(|b| b.join("+")).collect();
                    let mut stages: Vec<Stage> = variant
                        .iter()
                        .zip(&names)
                        .map(|(b, n)| stage(n, b))
                        .collect();
                    let mut name = names.join(" > ");
                    if *polish {
                        let all: Vec<String> = variant.iter().flatten().cloned().collect();
                        stages.push(stage("polish", &all));
                        name.push_str(" | polish");
                    }
                    strategies.push(Strategy { name, stages });
                }
            }
            StrategyGen::Cumulative { sets, sample } => {
                for order in orders(sets, *sample, cap, "cumulative", &mut notes)? {
                    let stages = (1..=order.len())
                        .map(|k| stage(&format!("+{}", order[k - 1]), &order[..k]))
                        .collect();
                    let name = order
                        .iter()
                        .enumerate()
                        .map(|(i, s)| if i == 0 { s.clone() } else { format!("+{s}") })
                        .collect::<Vec<_>>()
                        .join(" > ");
                    strategies.push(Strategy { name, stages });
                }
            }
        }
    }
    let mut seen = std::collections::HashSet::new();
    if let Some(s) = strategies.iter().find(|s| !seen.insert(s.name.clone())) {
        return Err(format!("two strategies are named {:?}", s.name));
    }
    Ok(Generated { strategies, notes })
}

/// An `lsq` stage freeing the named sets.
fn stage(name: &str, sets: &[String]) -> Stage {
    let free: Vec<String> = sets.iter().map(|s| format!("set:{s}")).collect();
    Stage {
        free,
        ..Stage::new(name, &[])
    }
}

/// The orders of `sets` a generator runs.
fn orders(
    sets: &[String],
    sample: Option<(usize, u64)>,
    cap: usize,
    what: &str,
    notes: &mut Vec<String>,
) -> Result<Vec<Vec<String>>, String> {
    if sets.is_empty() {
        return Err(format!("{what}: no sets"));
    }
    let mut unique = sets.to_vec();
    unique.sort();
    unique.dedup();
    if unique.len() != sets.len() {
        return Err(format!("{what}: a set is named twice in {sets:?}"));
    }
    let total = (1..=sets.len()).fold(1usize, |a, k| a.saturating_mul(k));
    let (count, seed) = match sample {
        Some((k, seed)) if k < total => (k, seed),
        Some(_) => (total, 0),
        None if total > cap => {
            notes.push(format!(
                "{what} over {sets:?}: {total} orders, a seeded sample of {cap} (seed 0)"
            ));
            (cap, 0)
        }
        None => (total, 0),
    };
    let pick = |idx: &[usize]| idx.iter().map(|&i| sets[i].clone()).collect();
    if count == total {
        // Every order, lexicographic in the given positions.
        let mut idx: Vec<usize> = (0..sets.len()).collect();
        let mut all = vec![pick(&idx)];
        while next_permutation(&mut idx) {
            all.push(pick(&idx));
        }
        return Ok(all);
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    while out.len() < count {
        let mut idx: Vec<usize> = (0..sets.len()).collect();
        idx.shuffle(&mut rng);
        if seen.insert(idx.clone()) {
            out.push(pick(&idx));
        }
    }
    Ok(out)
}

fn next_permutation(v: &mut [usize]) -> bool {
    let Some(i) = (1..v.len()).rev().find(|&i| v[i - 1] < v[i]) else {
        return false;
    };
    let j = (i..v.len()).rev().find(|&j| v[j] > v[i - 1]).unwrap();
    v.swap(i - 1, j);
    v[i..].reverse();
    true
}

// --- tournament on one problem --------------------------------------------

/// One strategy's result in a [`tournament`].
#[derive(Debug, Clone)]
pub struct TournamentEntry {
    pub name: String,
    /// 1 is best.
    pub rank: usize,
    /// Final χ² on every dataset; infinite when the strategy failed.
    pub chi2: f64,
    pub chi2_reduced: f64,
    pub bic: f64,
    /// RMS of `std_error / |value|` over the parameters with an error bar.
    pub rel_uncertainty: f64,
    pub evals: usize,
    /// Wall clock; not deterministic, not ranked on.
    pub time: Duration,
    pub error: Option<String>,
    pub report: Option<StrategyReport>,
}

/// Run every strategy on `problem` (in parallel, collected in order) and
/// rank by final χ², then evaluations, then name; failures last.
pub fn tournament(
    problem: &FitProblem,
    strategies: &[Strategy],
    sets: &[ParamSet],
) -> Vec<TournamentEntry> {
    let mut out: Vec<TournamentEntry> = strategies
        .par_iter()
        .map(|s| {
            let t0 = Instant::now();
            let r = run_strategy(problem, s, sets);
            let time = t0.elapsed();
            match r {
                Ok(r) => TournamentEntry {
                    name: s.name.clone(),
                    rank: 0,
                    chi2: r.report.chi2,
                    chi2_reduced: r.report.chi2_reduced,
                    bic: r.report.bic,
                    rel_uncertainty: rel_uncertainty(&r),
                    evals: r.evals,
                    time,
                    error: None,
                    report: Some(r),
                },
                Err(e) => TournamentEntry {
                    name: s.name.clone(),
                    rank: 0,
                    chi2: f64::INFINITY,
                    chi2_reduced: f64::NAN,
                    bic: f64::NAN,
                    rel_uncertainty: f64::NAN,
                    evals: 0,
                    time,
                    error: Some(e),
                    report: None,
                },
            }
        })
        .collect();
    out.sort_by(|a, b| {
        a.chi2
            .total_cmp(&b.chi2)
            .then(a.evals.cmp(&b.evals))
            .then(a.name.cmp(&b.name))
    });
    for (i, e) in out.iter_mut().enumerate() {
        e.rank = i + 1;
    }
    out
}

fn rel_uncertainty(r: &StrategyReport) -> f64 {
    let rel: Vec<f64> = r
        .report
        .parameters
        .iter()
        .filter_map(|p| {
            p.std_error
                .map(|s| s / p.value.abs().max(f64::MIN_POSITIVE))
        })
        .collect();
    if rel.is_empty() {
        f64::NAN
    } else {
        (rel.iter().map(|x| x * x).sum::<f64>() / rel.len() as f64).sqrt()
    }
}

// --- truths and the benchmark ---------------------------------------------

/// A drawn truth and the measurements it gives.
#[derive(Debug, Clone)]
pub struct Truth {
    pub seed: u64,
    pub model: FitModel,
    /// Every free parameter's true value, in problem order.
    pub values: Vec<(String, f64)>,
    /// `problem` with its datasets' values simulated on the truth plus
    /// noise; its model and starts are the problem's own.
    pub problem: FitProblem,
}

/// Draw a truth for `problem` (module docs): free parameters drawn from
/// their groups' error laws, then every dataset simulated with noise.
pub fn draw_truth(problem: &FitProblem, seed: u64) -> Result<Truth, String> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut start = problem.clone();
    start.parameters.iter_mut().for_each(|p| p.fixed = true);
    let start = evaluate(&start)?.model;

    let mut at_truth = problem.clone();
    for p in at_truth.parameters.iter_mut() {
        if p.fixed {
            continue;
        }
        let target = start.resolve(&p.addr)?;
        let x0 = start.read(&target)?;
        let x = draw(&p.addr, x0, start.group_of(&target), &mut rng);
        let (dlb, dub) = default_bounds(&p.addr, x0);
        let (lb, ub) = (p.lb.unwrap_or(dlb), p.ub.unwrap_or(dub));
        p.value = Some(x.clamp(lb, ub));
        p.fixed = true;
    }
    let model = evaluate(&at_truth)?.model;
    let values = problem
        .parameters
        .iter()
        .filter(|p| !p.fixed)
        .map(|p| Ok((p.addr.to_string(), model.get(&p.addr)?)))
        .collect::<Result<Vec<_>, String>>()?;

    // Simulate on every point, windows aside, then add noise.
    let mut full = problem.clone();
    full.model = model.clone();
    full.data.iter_mut().for_each(|d| d.window = None);
    full.nk_data.iter_mut().for_each(|d| d.window = None);
    let sim = full.simulate()?;
    let mut measured = problem.clone();
    let (spectra, nk) = sim.split_at(problem.data.len());
    for (d, s) in measured.data.iter_mut().zip(spectra) {
        d.values = s
            .iter()
            .enumerate()
            .map(|(i, v)| v + noise(d.sigma_at(i), &mut rng))
            .collect();
    }
    for (d, s) in measured.nk_data.iter_mut().zip(nk) {
        let n = d.wavelengths.len();
        for (c, col) in d.columns.iter_mut().zip(s.chunks(n)) {
            c.values = col
                .iter()
                .enumerate()
                .map(|(i, v)| v + noise(c.sigma_at(i), &mut rng))
                .collect();
        }
    }
    Ok(Truth {
        seed,
        model,
        values,
        problem: measured,
    })
}

fn noise(sigma: f64, rng: &mut StdRng) -> f64 {
    Normal::new(0.0, sigma).map_or(0.0, |n| n.sample(rng))
}

/// One parameter's true value (module docs' table).
fn draw(addr: &ParamAddr, x0: f64, group: Option<&Group>, rng: &mut StdRng) -> f64 {
    let Some(g) = group else {
        return x0;
    };
    let on = |m: ErrorMask| g.error_mask[m as usize] != 0;
    match addr {
        ParamAddr::Layer { field, .. } => match field {
            LayerField::Thickness if on(ErrorMask::Thickness) => g.thickness_error(x0, rng),
            LayerField::InhDelta if on(ErrorMask::InhDelta) => g.inh_delta_error(x0, rng),
            LayerField::Roughness if on(ErrorMask::Roughness) => g.sr_roughness_error(x0, rng),
            LayerField::InterfaceThickness if on(ErrorMask::Interface) => {
                g.interface_error(x0, rng)
            }
            _ => x0,
        },
        ParamAddr::Group { field, .. } => {
            let (mask, law, params, factor) = match field {
                GroupField::NFactor => (ErrorMask::NReal, g.n_error_type, &g.n_error_params, true),
                GroupField::KFactor => (ErrorMask::NImag, g.k_error_type, &g.k_error_params, true),
                GroupField::ThickFactor => (
                    ErrorMask::Thickness,
                    g.thickness_error_type,
                    &g.thickness_error_params,
                    true,
                ),
                GroupField::ThickSummand => (
                    ErrorMask::Thickness,
                    g.thickness_error_type,
                    &g.thickness_error_params,
                    false,
                ),
                GroupField::InhDeltaSummand => (
                    ErrorMask::InhDelta,
                    g.inh_delta_error_type,
                    &g.inh_delta_error_params,
                    false,
                ),
                GroupField::RoughnessSummand => (
                    ErrorMask::Roughness,
                    g.roughness_error_type,
                    &g.roughness_error_params,
                    false,
                ),
                GroupField::InterfaceSummand => (
                    ErrorMask::Interface,
                    g.interface_error_type,
                    &g.interface_error_params,
                    false,
                ),
            };
            if !on(mask) {
                return x0;
            }
            if factor {
                x0 * Group::apply_error(1.0, law, &relative_part(params), rng)
            } else {
                x0 + Group::apply_error(0.0, law, params, rng)
            }
        }
        ParamAddr::Material { .. } | ParamAddr::Data { .. } => x0,
    }
}

/// `params` with the absolute laws zeroed.
fn relative_part(params: &ErrorParams) -> ErrorParams {
    ErrorParams {
        abs_mean_delta_g: 0.0,
        abs_std_dev: 0.0,
        abs_mean_delta_h: 0.0,
        abs_variance: 0.0,
        ..params.clone()
    }
}

/// One strategy on one truth.
#[derive(Debug, Clone)]
pub struct RunScore {
    pub seed: u64,
    pub chi2_reduced: f64,
    /// RMS over the free parameters of `(fit − truth) / |truth|`
    /// (absolute where the truth is 0); infinite when the run failed.
    pub recovery: f64,
    /// Largest `|fit − truth| / std_error`.
    pub max_z: f64,
    /// Parameters within 2σ of their truth.
    pub within_2sigma: usize,
    pub evals: usize,
    /// Wall clock; not deterministic, not ranked on.
    pub time: Duration,
    pub error: Option<String>,
}

/// One strategy over every truth.
#[derive(Debug, Clone)]
pub struct BenchmarkEntry {
    pub name: String,
    /// 1 is best.
    pub rank: usize,
    /// Mean `recovery` over the runs that finished.
    pub mean_recovery: f64,
    pub worst_recovery: f64,
    pub mean_chi2_reduced: f64,
    pub mean_evals: f64,
    /// Runs that finished with every parameter within 2σ.
    pub successes: usize,
    pub failures: usize,
    /// One per seed, in seed order.
    pub runs: Vec<RunScore>,
}

#[derive(Debug, Clone)]
pub struct BenchmarkReport {
    pub seeds: Vec<u64>,
    pub n_free: usize,
    /// Ranked: fewest failures, then lowest mean recovery error, then
    /// fewest evaluations, then name.
    pub entries: Vec<BenchmarkEntry>,
}

/// Every strategy on a truth drawn for each seed ([`draw_truth`]), from
/// the problem's starts; ranked by recovery of the truth. Runs are
/// independent and run in parallel, collected in order: the report is
/// the same for the same request, times aside.
pub fn benchmark(
    problem: &FitProblem,
    strategies: &[Strategy],
    sets: &[ParamSet],
    seeds: &[u64],
) -> Result<BenchmarkReport, String> {
    if strategies.is_empty() || seeds.is_empty() {
        return Err("benchmark: give at least one strategy and one seed".into());
    }
    let truths: Vec<Truth> = seeds
        .par_iter()
        .map(|&s| draw_truth(problem, s).map_err(|e| format!("truth for seed {s}: {e}")))
        .collect::<Result<_, _>>()?;
    let pairs: Vec<(usize, usize)> = (0..strategies.len())
        .flat_map(|i| (0..truths.len()).map(move |k| (i, k)))
        .collect();
    let scores: Vec<RunScore> = pairs
        .par_iter()
        .map(|&(i, k)| score(&truths[k], &strategies[i], sets))
        .collect();

    let n_free = truths[0].values.len();
    let mut entries: Vec<BenchmarkEntry> = strategies
        .iter()
        .zip(scores.chunks(truths.len()))
        .map(|(s, runs)| {
            let done: Vec<&RunScore> = runs.iter().filter(|r| r.error.is_none()).collect();
            let mean = |f: &dyn Fn(&RunScore) -> f64| {
                if done.is_empty() {
                    f64::INFINITY
                } else {
                    done.iter().map(|r| f(r)).sum::<f64>() / done.len() as f64
                }
            };
            BenchmarkEntry {
                name: s.name.clone(),
                rank: 0,
                mean_recovery: mean(&|r| r.recovery),
                worst_recovery: done
                    .iter()
                    .map(|r| r.recovery)
                    .fold(f64::NEG_INFINITY, f64::max),
                mean_chi2_reduced: mean(&|r| r.chi2_reduced),
                mean_evals: mean(&|r| r.evals as f64),
                successes: done.iter().filter(|r| r.within_2sigma == n_free).count(),
                failures: runs.len() - done.len(),
                runs: runs.to_vec(),
            }
        })
        .collect();
    entries.sort_by(|a, b| {
        a.failures
            .cmp(&b.failures)
            .then(a.mean_recovery.total_cmp(&b.mean_recovery))
            .then(a.mean_evals.total_cmp(&b.mean_evals))
            .then(a.name.cmp(&b.name))
    });
    for (i, e) in entries.iter_mut().enumerate() {
        e.rank = i + 1;
    }
    Ok(BenchmarkReport {
        seeds: seeds.to_vec(),
        n_free,
        entries,
    })
}

fn score(truth: &Truth, strategy: &Strategy, sets: &[ParamSet]) -> RunScore {
    let t0 = Instant::now();
    let r = run_strategy(&truth.problem, strategy, sets);
    let time = t0.elapsed();
    let r = match r {
        Ok(r) => r,
        Err(e) => {
            return RunScore {
                seed: truth.seed,
                chi2_reduced: f64::NAN,
                recovery: f64::INFINITY,
                max_z: f64::NAN,
                within_2sigma: 0,
                evals: 0,
                time,
                error: Some(e),
            };
        }
    };
    let mut sq = 0.0;
    let mut max_z: f64 = 0.0;
    let mut within = 0;
    for (addr, t) in &truth.values {
        let p = r.report.parameter(addr).expect("a declared parameter");
        let d = p.value - t;
        let rel = if *t == 0.0 { d } else { d / t.abs() };
        sq += rel * rel;
        if let Some(se) = p.std_error {
            let z = d.abs() / se;
            max_z = max_z.max(z);
            if z <= 2.0 {
                within += 1;
            }
        } else {
            max_z = f64::INFINITY;
        }
    }
    RunScore {
        seed: truth.seed,
        chi2_reduced: r.report.chi2_reduced,
        recovery: (sq / truth.values.len().max(1) as f64).sqrt(),
        max_z,
        within_2sigma: within,
        evals: r.evals,
        time,
        error: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(gens: &[StrategyGen], cap: usize) -> Vec<String> {
        generate(gens, cap)
            .unwrap()
            .strategies
            .into_iter()
            .map(|s| s.name)
            .collect()
    }

    fn sets(s: &[&str]) -> Vec<String> {
        s.iter().map(|x| x.to_string()).collect()
    }

    #[test]
    fn orderings_groupings_and_cumulative_name_their_arrangement() {
        let o = StrategyGen::Orderings {
            sets: sets(&["a", "b", "c"]),
            sample: None,
            polish: true,
        };
        assert_eq!(
            names(&[o], 24),
            [
                "a > b > c | polish",
                "a > c > b | polish",
                "b > a > c | polish",
                "b > c > a | polish",
                "c > a > b | polish",
                "c > b > a | polish"
            ]
        );
        let g = StrategyGen::Groupings {
            variants: vec![
                vec![sets(&["a", "b"]), sets(&["c"])],
                vec![sets(&["a"]), sets(&["b"]), sets(&["c"])],
            ],
            polish: false,
        };
        let made = generate(&[g], 24).unwrap();
        assert_eq!(made.strategies[0].name, "a+b > c");
        assert_eq!(made.strategies[0].stages[0].free, ["set:a", "set:b"]);
        assert_eq!(made.strategies[1].stages.len(), 3);
        let c = StrategyGen::Cumulative {
            sets: sets(&["a", "b"]),
            sample: None,
        };
        let made = generate(&[c], 24).unwrap();
        assert_eq!(made.strategies[0].name, "a > +b");
        assert_eq!(made.strategies[0].stages[1].free, ["set:a", "set:b"]);
        assert_eq!(made.strategies[1].name, "b > +a");
    }

    #[test]
    fn many_orders_are_sampled_and_say_so() {
        let five = sets(&["a", "b", "c", "d", "e"]);
        let o = StrategyGen::Orderings {
            sets: five.clone(),
            sample: None,
            polish: false,
        };
        let made = generate(std::slice::from_ref(&o), 24).unwrap();
        assert_eq!(made.strategies.len(), 24);
        assert!(made.notes[0].contains("120 orders"), "{:?}", made.notes);
        assert_eq!(made, generate(&[o], 24).unwrap());
        let s = StrategyGen::Orderings {
            sets: five,
            sample: Some((7, 3)),
            polish: false,
        };
        let made = generate(&[s], 24).unwrap();
        assert_eq!(made.strategies.len(), 7);
        assert!(made.notes.is_empty());
        let dup = StrategyGen::Orderings {
            sets: sets(&["a", "a"]),
            sample: None,
            polish: false,
        };
        assert!(generate(&[dup], 24).unwrap_err().contains("twice"));
    }

    #[test]
    fn generators_round_trip_through_json() {
        let text = r#"[
            {"orderings": {"sets": ["thicknesses", "groups"]}},
            {"groupings": {"variants": [[["a", "b"], ["c"]]], "polish": false}},
            {"cumulative": {"sets": ["a", "b"], "sample": [1, 9]}}
        ]"#;
        let gens: Vec<StrategyGen> = serde_json::from_str(text).unwrap();
        assert_eq!(
            gens[0],
            StrategyGen::Orderings {
                sets: sets(&["thicknesses", "groups"]),
                sample: None,
                polish: true
            }
        );
        let back: Vec<StrategyGen> =
            serde_json::from_str(&serde_json::to_string(&gens).unwrap()).unwrap();
        assert_eq!(back, gens);
    }
}
