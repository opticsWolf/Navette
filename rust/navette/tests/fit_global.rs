// SPDX-License-Identifier: LGPL-3.0-or-later
//! Fit engine, phases 6 and 7 (`docs/plans/general_fit_plan.md`): the
//! `global` stage. The phase-4 coating (H L H L, four thicknesses and two
//! index factors free) from starts drawn at random over a 0.5-1.8×
//! thickness and ±10 % index box: every global method followed by its
//! polish reaches the noise floor from every start. Least squares alone
//! does too there (two curves over 400-900 nm pin that stack); the same
//! design four times thicker, read through one curve, is rugged enough
//! to catch it, and differential evolution still is not caught.
//!
//! Navette's own DE runs in every build; basin's methods need
//! `opt-basin`, and with it the native DE is checked against basin's on
//! success rate and evaluations to a given chi2 over many seeds.
// The comparison helpers serve the basin tests only.
#![cfg_attr(not(feature = "opt-basin"), allow(dead_code, unused_imports))]

use std::collections::{BTreeMap, HashMap};

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;

use navette::fit::{
    FitModel, FitProblem, GlobalMethod, Parameter, Progress, SpectrumData, Stage, StageKind,
    Strategy, Watcher, draw_truth, evaluate, run_strategy, run_strategy_watched,
};
use navette::smatrix::angle::Radians;
use navette::structure::enums::{ErrorMask, LayerType};
use navette::structure::{Entry, ErrorParams, Group, Layer, MaterialSpec, Structure};

fn spec(model: &str, params: serde_json::Value) -> Entry {
    let params: BTreeMap<String, serde_json::Value> = serde_json::from_value(params).unwrap();
    Entry::Spec(MaterialSpec::new(model, params))
}

fn half_space(material: &str, role: LayerType) -> Layer {
    Layer {
        layer_type: role,
        ..Layer::film(0.0, material)
    }
}

fn erring(name: &str) -> Group {
    let mut g = Group::new(name);
    g.error_mask[ErrorMask::Thickness as usize] = 1;
    g.error_mask[ErrorMask::NReal as usize] = 1;
    g.thickness_error_params = ErrorParams {
        abs_std_dev: 0.0,
        abs_variance: 0.0,
        rel_std_dev: 0.02,
        ..ErrorParams::thickness()
    };
    g.n_error_params = ErrorParams {
        abs_std_dev: 0.0,
        abs_variance: 0.0,
        rel_std_dev: 0.01,
        ..ErrorParams::index()
    };
    g
}

/// The phase-4 benchmark coating, thicknesses times `scale`:
/// air | H L H L | glass.
fn coating(scale: f64) -> FitModel {
    let layers = vec![
        half_space("air", LayerType::Ambient),
        Layer::film(scale * NOMINAL[0], "TiO2"),
        Layer::film(scale * NOMINAL[1], "SiO2"),
        Layer::film(scale * NOMINAL[2], "TiO2"),
        Layer::film(scale * NOMINAL[3], "SiO2"),
        half_space("glass", LayerType::Substrate),
    ];
    let groups = HashMap::from([
        ("TiO2".to_string(), erring("TiO2")),
        ("SiO2".to_string(), erring("SiO2")),
    ]);
    let materials = HashMap::from([
        ("air".to_string(), spec("Konstant", json!({"n": 1.0}))),
        ("glass".to_string(), spec("Konstant", json!({"n": 1.52}))),
        (
            "TiO2".to_string(),
            spec("Cauchy", json!({"A": 2.2, "B": 0.02, "C": 0.001})),
        ),
        (
            "SiO2".to_string(),
            spec("Cauchy", json!({"A": 1.45, "B": 0.0036, "C": 0.0})),
        ),
    ]);
    FitModel::from_structure(&Structure::new(layers, groups), "main", materials).unwrap()
}

const NOMINAL: [f64; 4] = [60.0, 95.0, 60.0, 95.0];

/// Two curves (R at 8°, Rp at 45°) or, `one_curve`, R alone.
fn problem(scale: f64, one_curve: bool) -> FitProblem {
    let wl: Vec<f64> = (0..101).map(|i| 400.0 + 5.0 * i as f64).collect();
    let curve = |id: &str, ch: &str, deg: f64| {
        SpectrumData::new(
            id,
            ch,
            Radians::from_degrees(deg).unwrap(),
            wl.clone(),
            vec![0.0; wl.len()],
            0.002,
        )
    };
    let mut parameters: Vec<Parameter> = NOMINAL
        .iter()
        .enumerate()
        .map(|(i, d)| {
            Parameter::parse(&format!("layer[main:{}].thickness", i + 1))
                .unwrap()
                .bounds(0.5 * scale * d, 1.8 * scale * d)
        })
        .collect();
    for g in ["TiO2", "SiO2"] {
        parameters.push(
            Parameter::parse(&format!("group[{g}].n_factor"))
                .unwrap()
                .bounds(0.9, 1.1),
        );
    }
    let mut data = vec![curve("R", "R_avg", 8.0)];
    if !one_curve {
        data.push(curve("Rp", "Rp", 45.0));
    }
    FitProblem::new(coating(scale), data, parameters)
}

/// The truth's problem with starts drawn uniformly over the box.
fn random_start(truth_problem: &FitProblem, scale: f64, seed: u64) -> FitProblem {
    let mut p = truth_problem.clone();
    let mut rng = StdRng::seed_from_u64(1000 + seed);
    for (i, q) in p.parameters.iter_mut().enumerate() {
        let (lo, hi) = if i < 4 {
            (0.5 * scale * NOMINAL[i], 1.8 * scale * NOMINAL[i])
        } else {
            (0.9, 1.1)
        };
        q.value = Some(rng.random_range(lo..hi));
    }
    p
}

fn strategy(kind: StageKind) -> Strategy {
    let mut stage = Stage::new("all", &["*"]);
    stage.kind = kind;
    Strategy {
        name: format!("{:?}", stage.kind),
        stages: vec![stage],
    }
}

fn global_seeded(method: GlobalMethod, budget: usize, seed: u64) -> Strategy {
    strategy(StageKind::Global {
        method,
        budget,
        seed,
        population: None,
    })
}

fn global(method: GlobalMethod, budget: usize) -> Strategy {
    global_seeded(method, budget, 3)
}

const SEEDS: [u64; 6] = [1, 2, 3, 4, 5, 6];

/// One seed's truth draw: the problem and chi2 at the truth.
fn truth_and_floor(scale: f64, one_curve: bool, seed: u64) -> (FitProblem, f64) {
    let truth = draw_truth(&problem(scale, one_curve), seed).unwrap();
    let mut at_truth = truth.problem.clone();
    at_truth.model = truth.model.clone();
    for q in &mut at_truth.parameters {
        q.fixed = true;
        q.value = None;
    }
    let floor = evaluate(&at_truth).unwrap().chi2;
    (random_start(&truth.problem, scale, seed), floor)
}

/// Per seed, χ² of `strategy` from the random start over χ² at the
/// truth: at most 1 means the global minimum (the truth's basin, noise
/// moving its floor); far above 1, a wrong minimum.
fn vs_truth(scale: f64, one_curve: bool, strategy: &Strategy) -> Vec<f64> {
    SEEDS
        .iter()
        .map(|&seed| {
            let (p, floor) = truth_and_floor(scale, one_curve, seed);
            let r = run_strategy(&p, strategy, &[]).unwrap();
            r.report.chi2 / floor
        })
        .collect()
}

/// Found the global minimum (to the polish's tolerance).
fn solved(ratio: &f64) -> bool {
    *ratio <= 1.0 + 1e-6
}

#[test]
fn differential_evolution_solves_the_phase_4_coating_from_every_start() {
    let de = vs_truth(1.0, false, &global(GlobalMethod::De, 3000));
    eprintln!("de {de:?}");
    assert!(de.iter().all(solved), "{de:?}");
}

#[test]
fn greedy_de_too() {
    let c = vs_truth(1.0, false, &global(GlobalMethod::DeBest, 3000));
    eprintln!("de_best {c:?}");
    assert!(c.iter().all(solved), "{c:?}");
}

#[cfg(feature = "opt-basin")]
#[test]
fn basin_de_memetic_de_and_cma_es_too() {
    for method in [
        GlobalMethod::BasinDe,
        GlobalMethod::DeInject,
        GlobalMethod::CmaEs,
    ] {
        let c = vs_truth(1.0, false, &global(method, 3000));
        eprintln!("{method:?} {c:?}");
        assert!(c.iter().all(solved), "{method:?}: {c:?}");
    }
}

#[test]
fn the_same_seed_gives_the_same_bits() {
    let truth = draw_truth(&problem(1.0, false), 2).unwrap();
    let p = random_start(&truth.problem, 1.0, 2);
    let s = global(GlobalMethod::De, 600);
    let a = run_strategy(&p, &s, &[]).unwrap();
    let b = run_strategy(&p, &s, &[]).unwrap();
    assert_eq!(a.report.chi2.to_bits(), b.report.chi2.to_bits());
    assert_eq!(a.evals, b.evals);
}

#[test]
fn a_rugged_stack_catches_least_squares_but_not_de() {
    let lsq = vs_truth(4.0, true, &strategy(StageKind::Lsq));
    let de = vs_truth(4.0, true, &global(GlobalMethod::De, 6000));
    eprintln!("rugged: lsq {lsq:?}\nrugged: de {de:?}");
    assert!(lsq.iter().filter(|c| **c > 3.0).count() >= 2, "{lsq:?}");
    assert!(de.iter().all(solved), "{de:?}");
}

/// One run's outcome on the problem at `scale` / `one_curve`: solved after
/// the polish, the evaluations the global search spent before its chi2
/// first fell to `factor` times the truth's (`None`: never), and the best
/// chi2 the search itself reached, over the truth's.
fn outcome(
    method: GlobalMethod,
    (scale, one_curve): (f64, bool),
    seed: u64,
    budget: usize,
    factor: f64,
) -> (bool, Option<usize>, f64) {
    let (p, floor) = truth_and_floor(scale, one_curve, seed);
    let seen = std::sync::Mutex::new((None::<usize>, f64::INFINITY));
    let w = Watcher::new(|e| {
        if let Progress::Searched { evals, chi2, .. } = e {
            let mut s = seen.lock().unwrap();
            s.1 = s.1.min(*chi2);
            if *chi2 <= factor * floor {
                s.0 = Some(s.0.map_or(*evals, |v| v.min(*evals)));
            }
        }
        Ok(())
    });
    let s = global_seeded(method, budget, 100 + seed);
    let r = run_strategy_watched(&p, &s, &[], Some(&w)).unwrap();
    drop(w);
    let (reached, best) = seen.into_inner().unwrap();
    (solved(&(r.report.chi2 / floor)), reached, best / floor)
}

/// Over `seeds`: runs solved, median evaluations to `factor` (a run that
/// never got there counts as the whole budget) and the median of the
/// search's best chi2 over the truth's.
fn stats(
    method: GlobalMethod,
    problem: (f64, bool),
    budget: usize,
    factor: f64,
) -> (usize, usize, f64) {
    let runs: Vec<_> = (1..=16)
        .map(|s| outcome(method, problem, s, budget, factor))
        .collect();
    let solved = runs.iter().filter(|r| r.0).count();
    let mut evals: Vec<usize> = runs.iter().map(|r| r.1.unwrap_or(budget)).collect();
    evals.sort_unstable();
    let mut best: Vec<f64> = runs.iter().map(|r| r.2).collect();
    best.sort_by(f64::total_cmp);
    eprintln!("{method:?} {problem:?}: {solved}/16 solved, evals {evals:?}, best {best:?}");
    (solved, evals[8], best[8])
}

/// The plan's acceptance check for the native DE, against basin's DE on
/// 16 seeds each: on the phase-4 coating, no more evaluations (median, 25 %
/// slack for the different random streams) to reach ten times the truth's
/// chi2; on the rugged stack, at least basin's success rate and a search
/// result no worse (median, same slack).
#[cfg(feature = "opt-basin")]
#[test]
fn native_de_is_not_worse_than_basin_de() {
    let (ns, ne, _) = stats(GlobalMethod::De, (1.0, false), 3000, 10.0);
    let (bs, be, _) = stats(GlobalMethod::BasinDe, (1.0, false), 3000, 10.0);
    assert!(ns >= bs, "phase 4: native {ns} solved vs basin {bs}");
    assert!(
        ne as f64 <= 1.25 * be as f64,
        "phase 4: native {ne} evals vs basin {be}"
    );
    let (ns, _, nb) = stats(GlobalMethod::De, (4.0, true), 4000, 10.0);
    let (bs, _, bb) = stats(GlobalMethod::BasinDe, (4.0, true), 4000, 10.0);
    assert!(ns >= bs, "rugged: native {ns} solved vs basin {bs}");
    assert!(nb <= 1.25 * bb, "rugged: native best {nb} vs basin {bb}");
}
