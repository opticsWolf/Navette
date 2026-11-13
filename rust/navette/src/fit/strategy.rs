// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::strategy — staged fitting (plan §3.8).
//!
//! A [`Strategy`] is an ordered list of [`Stage`]s; each says what is
//! free, which data is scored, and how to solve. Every stage starts from
//! the previous one's end state, and a stage is one ordinary [`fit`] of
//! the problem with everything outside its free set held, so stages add no
//! second solver. Stages are plain data and round-trip through JSON:
//!
//! ```json
//! {"name": "thickness first", "stages": [
//!   {"name": "thicknesses", "free": ["kind=thickness"], "window": [500, 900]},
//!   {"name": "indices", "free": ["group[*].n_factor"]},
//!   {"name": "polish", "free": ["*"], "repeat": {"until_no_gain": {"rel": 1e-3, "max": 3}}}
//! ]}
//! ```
//!
//! Stage kinds: `lsq` (default), `scan` (a 1-D grid over one parameter,
//! best point kept, then `lsq`), and `multi_start` (seeded uniform or
//! Latin-hypercube starts in the bounds, `lsq` from each, best kept). The
//! global stage is phase 6.

use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::lsq::backends::LsqBackend;

use super::model::FitModel;
use super::sets::{ParamSet, select};
use super::{FitProblem, FitReport, evaluate, fit};

/// An ordered list of stages.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Strategy {
    pub name: String,
    pub stages: Vec<Stage>,
}

impl Strategy {
    pub fn from_json(text: &str) -> Result<Self, String> {
        serde_json::from_str(text).map_err(|e| format!("strategy: {e}"))
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(self).expect("a strategy is plain data")
    }
}

/// One stage: what is free, what is scored, how to solve.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Stage {
    pub name: String,
    /// Selectors (module docs of `fit::sets`); the union is free, the rest
    /// of the declared parameters held at their current values. A
    /// parameter declared `fixed` is never freed.
    pub free: Vec<String>,
    /// Dataset ids scored (spectra and nk data); `None` scores all.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub datasets: Option<Vec<String>>,
    /// Per dataset id, a factor on its weight.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub weights: Vec<(String, f64)>,
    /// Wavelength window (nm) for every scored dataset, over its own.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub window: Option<(f64, f64)>,
    #[serde(default)]
    pub kind: StageKind,
    #[serde(default, skip_serializing_if = "StageLsq::is_empty")]
    pub lsq: StageLsq,
    #[serde(default)]
    pub repeat: Repeat,
    #[serde(default)]
    pub accept: Accept,
}

impl Stage {
    /// An `lsq` stage freeing `free`, everything else default.
    pub fn new(name: impl Into<String>, free: &[&str]) -> Self {
        Self {
            name: name.into(),
            free: free.iter().map(|s| s.to_string()).collect(),
            datasets: None,
            weights: Vec::new(),
            window: None,
            kind: StageKind::Lsq,
            lsq: StageLsq::default(),
            repeat: Repeat::Once,
            accept: Accept::Always,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum StageKind {
    #[default]
    Lsq,
    /// `points` values of `param` (a selector naming one free parameter)
    /// evenly over its bounds, ends included; the best is kept, then `lsq`.
    Scan { param: String, points: usize },
    /// The current point plus `starts` sampled ones, `lsq` from each; the
    /// lowest χ² wins (the earliest on a tie). An open bound side samples
    /// to `x ± max(|x|/2, 1e-3)`.
    MultiStart {
        starts: usize,
        #[serde(default)]
        sampler: Sampler,
        #[serde(default)]
        seed: u64,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Sampler {
    #[default]
    Uniform,
    /// Latin hypercube: one start in each of `starts` strata per parameter.
    Lhs,
}

/// Overrides of the problem's solver settings for one stage.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StageLsq {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub backend: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_iterations: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_evals: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ftol: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub xtol: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gtol: Option<f64>,
}

impl StageLsq {
    fn is_empty(&self) -> bool {
        *self == Self::default()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum Repeat {
    #[default]
    Once,
    /// Run again while a run lowers χ² by more than `rel` (relative), at
    /// most `max` runs.
    UntilNoGain { rel: f64, max: usize },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Accept {
    #[default]
    Always,
    /// Keep the stage's result only if it lowered the stage's χ²; else
    /// roll back to the stage's start.
    IfBetter,
}

/// One stage after it ran.
#[derive(Debug, Clone)]
pub struct StageReport {
    pub name: String,
    /// The addresses the stage freed.
    pub free: Vec<String>,
    /// The stage's χ² (its datasets, window and weights) at its start.
    pub chi2_start: f64,
    /// Least-squares runs: repeats, multi-start starts.
    pub runs: usize,
    /// `false` when `IfBetter` rolled the stage back.
    pub accepted: bool,
    /// The kept run's report (on the stage's data).
    pub report: FitReport,
}

/// Result of [`run_strategy`].
#[derive(Debug, Clone)]
pub struct StrategyReport {
    pub name: String,
    pub stages: Vec<StageReport>,
    /// Every dataset at the end state, with the covariance of every
    /// parameter some stage freed.
    pub report: FitReport,
}

/// Run `strategy` on `problem`. `sets` are the user parameter sets
/// (built-ins need no declaration).
pub fn run_strategy(
    problem: &FitProblem,
    strategy: &Strategy,
    sets: &[ParamSet],
) -> Result<StrategyReport, String> {
    if strategy.stages.is_empty() {
        return Err(format!("strategy {:?}: no stages", strategy.name));
    }
    // Starts and ties written once; every stage reads the model.
    let mut model = evaluate(&held(problem, &vec![false; problem.parameters.len()]))?.model;
    let mut ever = vec![false; problem.parameters.len()];
    let mut stages = Vec::new();
    for stage in &strategy.stages {
        let bad = |e: String| format!("stage {:?}: {e}", stage.name);
        let free = select(&stage.free, &problem.parameters, sets).map_err(bad)?;
        let free: Vec<bool> = free
            .iter()
            .zip(&problem.parameters)
            .map(|(f, p)| *f && !p.fixed)
            .collect();
        if !free.iter().any(|f| *f) {
            return Err(bad("frees nothing (every match is declared fixed)".into()));
        }
        let report = run_stage(problem, stage, &free, sets, &model).map_err(bad)?;
        for (e, f) in ever.iter_mut().zip(&free) {
            *e |= f;
        }
        if report.accepted {
            model = report.report.model.clone();
        }
        stages.push(report);
    }
    let mut end = held(problem, &ever);
    end.model = model;
    Ok(StrategyReport {
        name: strategy.name.clone(),
        stages,
        report: evaluate(&end)?,
    })
}

/// `problem` with every parameter not flagged `free` held, starts read
/// from the model.
fn held(problem: &FitProblem, free: &[bool]) -> FitProblem {
    let mut p = problem.clone();
    for (q, f) in p.parameters.iter_mut().zip(free) {
        q.fixed = !f;
    }
    p
}

fn run_stage(
    problem: &FitProblem,
    stage: &Stage,
    free: &[bool],
    sets: &[ParamSet],
    model: &FitModel,
) -> Result<StageReport, String> {
    let mut sp = held(problem, free);
    sp.model = model.clone();
    for q in &mut sp.parameters {
        q.value = None;
    }
    scope_data(&mut sp, stage)?;
    apply_lsq(&mut sp, &stage.lsq)?;
    let start = evaluate(&sp)?;
    let addrs: Vec<String> = sp
        .parameters
        .iter()
        .zip(free)
        .filter(|(_, f)| **f)
        .map(|(p, _)| p.addr.to_string())
        .collect();

    let (mut best, mut runs) = match &stage.kind {
        StageKind::Lsq => (fit(&sp)?, 1),
        StageKind::Scan { param, points } => {
            let pick = select(std::slice::from_ref(param), &sp.parameters, sets)?;
            let hits: Vec<usize> = (0..pick.len()).filter(|&i| pick[i] && free[i]).collect();
            let [j] = hits[..] else {
                return Err(format!(
                    "scan {param:?} must name exactly one free parameter, names {}",
                    hits.len()
                ));
            };
            if *points < 2 {
                return Err(format!("scan needs at least 2 points, got {points}"));
            }
            let (lb, ub) = (start.parameters[j].lb, start.parameters[j].ub);
            if !(lb.is_finite() && ub.is_finite()) {
                return Err(format!(
                    "scan {param:?}: bounds [{lb}, {ub}] are open; give finite bounds"
                ));
            }
            let mut probe = sp.clone();
            for q in &mut probe.parameters {
                q.fixed = true;
            }
            let chi2: Vec<f64> = (0..*points)
                .into_par_iter()
                .map(|i| {
                    let x = lb + (ub - lb) * i as f64 / (*points - 1) as f64;
                    let mut p = probe.clone();
                    p.parameters[j].value = Some(x);
                    evaluate(&p).map(|r| r.chi2).unwrap_or(f64::INFINITY)
                })
                .collect();
            let i = argmin(&chi2).ok_or("every scan point failed to evaluate")?;
            sp.parameters[j].value = Some(lb + (ub - lb) * i as f64 / (*points - 1) as f64);
            (fit(&sp)?, 1)
        }
        StageKind::MultiStart {
            starts,
            sampler,
            seed,
        } => {
            let boxes: Vec<(usize, f64, f64)> = (0..free.len())
                .filter(|&i| free[i])
                .map(|i| {
                    let p = &start.parameters[i];
                    let half = (0.5 * p.value.abs()).max(1e-3);
                    let lo = if p.lb.is_finite() {
                        p.lb
                    } else {
                        p.value - half
                    };
                    let hi = if p.ub.is_finite() {
                        p.ub
                    } else {
                        p.value + half
                    };
                    (i, lo.max(p.lb), hi.min(p.ub))
                })
                .collect();
            let points = sample(&boxes, *starts, *sampler, *seed);
            let mut candidates = vec![sp.clone()];
            for pt in points {
                let mut p = sp.clone();
                for (&(i, _, _), x) in boxes.iter().zip(pt) {
                    p.parameters[i].value = Some(x);
                }
                candidates.push(p);
            }
            let reports: Vec<Option<FitReport>> =
                candidates.par_iter().map(|p| fit(p).ok()).collect();
            let chi2: Vec<f64> = reports
                .iter()
                .map(|r| r.as_ref().map_or(f64::INFINITY, |r| r.chi2))
                .collect();
            let i = argmin(&chi2).ok_or("every start failed")?;
            let n = reports.len();
            (reports.into_iter().nth(i).flatten().unwrap(), n)
        }
    };

    if let Repeat::UntilNoGain { rel, max } = stage.repeat {
        let mut runs_here = 1;
        while runs_here < max.max(1) {
            let mut again = sp.clone();
            again.model = best.model.clone();
            // From where the last run ended, not from a scan's grid point.
            for q in &mut again.parameters {
                q.value = None;
            }
            let next = fit(&again)?;
            runs += 1;
            runs_here += 1;
            let gain = (best.chi2 - next.chi2) / best.chi2.max(f64::MIN_POSITIVE);
            let keep = next.chi2 < best.chi2;
            if keep {
                best = next;
            }
            if !(gain > rel) {
                break;
            }
        }
    }

    let accepted = match stage.accept {
        Accept::Always => true,
        Accept::IfBetter => best.chi2 < start.chi2,
    };
    if !accepted {
        best = start.clone();
    }
    Ok(StageReport {
        name: stage.name.clone(),
        free: addrs,
        chi2_start: start.chi2,
        runs,
        accepted,
        report: best,
    })
}

/// Index of the smallest finite value, the first on a tie.
fn argmin(v: &[f64]) -> Option<usize> {
    let mut best: Option<usize> = None;
    for (i, x) in v.iter().enumerate() {
        if x.is_finite() && best.is_none_or(|b| *x < v[b]) {
            best = Some(i);
        }
    }
    best
}

/// `n` points in the boxes, seeded.
fn sample(boxes: &[(usize, f64, f64)], n: usize, sampler: Sampler, seed: u64) -> Vec<Vec<f64>> {
    let mut rng = StdRng::seed_from_u64(seed);
    match sampler {
        Sampler::Uniform => (0..n)
            .map(|_| {
                boxes
                    .iter()
                    .map(|&(_, lo, hi)| lo + (hi - lo) * rng.random::<f64>())
                    .collect()
            })
            .collect(),
        Sampler::Lhs => {
            // Per parameter, one draw in each stratum, strata shuffled.
            let columns: Vec<Vec<f64>> = boxes
                .iter()
                .map(|&(_, lo, hi)| {
                    let mut strata: Vec<usize> = (0..n).collect();
                    strata.shuffle(&mut rng);
                    strata
                        .into_iter()
                        .map(|s| lo + (hi - lo) * (s as f64 + rng.random::<f64>()) / n as f64)
                        .collect()
                })
                .collect();
            (0..n)
                .map(|k| columns.iter().map(|c| c[k]).collect())
                .collect()
        }
    }
}

/// The stage's datasets, window and weights.
fn scope_data(p: &mut FitProblem, stage: &Stage) -> Result<(), String> {
    let known: Vec<String> = p
        .data
        .iter()
        .map(|d| d.id.clone())
        .chain(p.nk_data.iter().map(|d| d.id.clone()))
        .collect();
    let name_ok = |id: &String| {
        if known.contains(id) {
            Ok(())
        } else {
            Err(format!("no dataset {id:?}; the problem has {known:?}"))
        }
    };
    if let Some(ids) = &stage.datasets {
        ids.iter().try_for_each(name_ok)?;
        p.data.retain(|d| ids.contains(&d.id));
        p.nk_data.retain(|d| ids.contains(&d.id));
    }
    for (id, factor) in &stage.weights {
        name_ok(id)?;
        if !(*factor > 0.0 && factor.is_finite()) {
            return Err(format!("weight factor {factor} for {id:?} is not positive"));
        }
        p.data
            .iter_mut()
            .filter(|d| d.id == *id)
            .for_each(|d| d.weight *= factor);
        p.nk_data
            .iter_mut()
            .filter(|d| d.id == *id)
            .for_each(|d| d.weight *= factor);
    }
    if let Some(w) = stage.window {
        p.data.iter_mut().for_each(|d| d.window = Some(w));
        p.nk_data.iter_mut().for_each(|d| d.window = Some(w));
    }
    Ok(())
}

fn apply_lsq(p: &mut FitProblem, s: &StageLsq) -> Result<(), String> {
    if let Some(b) = &s.backend {
        p.lsq.backend = LsqBackend::parse(b)?;
    }
    if let Some(v) = s.max_iterations {
        p.lsq.max_iterations = v;
    }
    if let Some(v) = s.max_evals {
        p.lsq.max_evals = v;
    }
    if let Some(v) = s.ftol {
        p.lsq.ftol = v;
    }
    if let Some(v) = s.xtol {
        p.lsq.xtol = v;
    }
    if let Some(v) = s.gtol {
        p.lsq.gtol = v;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stages_round_trip_through_json() {
        let text = r#"{"name": "s", "stages": [
            {"name": "a", "free": ["kind=thickness"], "window": [500, 900]},
            {"name": "b", "free": ["*"], "kind": {"scan": {"param": "layer[m:1].thickness", "points": 9}},
             "accept": "if_better"},
            {"name": "c", "free": ["*"], "kind": {"multi_start": {"starts": 8, "sampler": "lhs", "seed": 3}},
             "repeat": {"until_no_gain": {"rel": 0.001, "max": 3}}, "lsq": {"backend": "trf"}}
        ]}"#;
        let s = Strategy::from_json(text).unwrap();
        assert_eq!(s.stages[0].kind, StageKind::Lsq);
        assert_eq!(s.stages[1].accept, Accept::IfBetter);
        assert_eq!(
            s.stages[2].repeat,
            Repeat::UntilNoGain { rel: 0.001, max: 3 }
        );
        assert_eq!(Strategy::from_json(&s.to_json()).unwrap(), s);
        assert!(
            Strategy::from_json(
                r#"{"name": "s", "stages": [{"name": "a", "free": [], "colour": 1}]}"#
            )
            .unwrap_err()
            .contains("colour")
        );
    }

    #[test]
    fn lhs_puts_one_start_in_each_stratum() {
        let pts = sample(&[(0, 0.0, 1.0), (1, 10.0, 20.0)], 5, Sampler::Lhs, 7);
        for (d, (lo, hi)) in [(0.0, 1.0), (10.0, 20.0)].into_iter().enumerate() {
            let mut strata: Vec<usize> = pts
                .iter()
                .map(|p| ((p[d] - lo) / (hi - lo) * 5.0).floor() as usize)
                .collect();
            strata.sort();
            assert_eq!(strata, [0, 1, 2, 3, 4]);
        }
        assert_eq!(
            pts,
            sample(&[(0, 0.0, 1.0), (1, 10.0, 20.0)], 5, Sampler::Lhs, 7)
        );
    }
}
