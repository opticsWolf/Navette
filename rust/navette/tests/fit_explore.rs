// SPDX-License-Identifier: LGPL-3.0-or-later
//! Fit engine, phase 4 (`docs/plans/general_fit_plan.md` §3.8.1, §3.8.2):
//! truths drawn from the groups' error laws, the strategy tournament and
//! the recovery benchmark.
//!
//! The contract: a truth is a deterministic function of its seed and
//! moves only the parameters whose group switches that error channel on;
//! a benchmark is the same report for the same request (times aside) and
//! ranks a strategy that cannot express the truth below one that can.

use std::collections::{BTreeMap, HashMap};

use serde_json::json;

use navette::fit::{
    FitModel, FitProblem, Parameter, SpectrumData, Strategy, StrategyGen, benchmark, draw_truth,
    evaluate, generate, tournament,
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

/// A group with 2 % relative thickness and 1 % relative index errors.
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

/// air | H L H L | glass; `low_errs` switches the SiO2 group's errors.
fn coating(low_errs: bool) -> FitModel {
    let layers = vec![
        half_space("air", LayerType::Ambient),
        Layer::film(60.0, "TiO2"),
        Layer::film(95.0, "SiO2"),
        Layer::film(60.0, "TiO2"),
        Layer::film(95.0, "SiO2"),
        half_space("glass", LayerType::Substrate),
    ];
    let low = if low_errs {
        erring("SiO2")
    } else {
        Group::new("SiO2")
    };
    let groups = HashMap::from([
        ("TiO2".to_string(), erring("TiO2")),
        ("SiO2".to_string(), low),
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

const SIGMA: f64 = 0.002;

fn problem(low_errs: bool) -> FitProblem {
    let wl: Vec<f64> = (0..101).map(|i| 400.0 + 5.0 * i as f64).collect();
    let curve = |id: &str, ch: &str, deg: f64| {
        SpectrumData::new(
            id,
            ch,
            Radians::from_degrees(deg).unwrap(),
            wl.clone(),
            vec![0.0; wl.len()],
            SIGMA,
        )
    };
    let mut parameters: Vec<Parameter> = (1..=4)
        .map(|i| Parameter::parse(&format!("layer[main:{i}].thickness")).unwrap())
        .collect();
    parameters.push(Parameter::parse("group[TiO2].n_factor").unwrap());
    parameters.push(Parameter::parse("group[SiO2].n_factor").unwrap());
    parameters.push(Parameter::parse("material[TiO2].A").unwrap().fixed());
    FitProblem::new(
        coating(low_errs),
        vec![curve("R", "R_avg", 8.0), curve("Rp", "Rp", 45.0)],
        parameters,
    )
}

#[test]
fn a_truth_is_its_seed_and_moves_only_switched_channels() {
    let p = problem(false);
    let a = draw_truth(&p, 7).unwrap();
    let b = draw_truth(&p, 7).unwrap();
    let c = draw_truth(&p, 8).unwrap();
    assert_eq!(a.values.len(), 6, "fixed parameters are not drawn");
    for ((n, x), (_, y)) in a.values.iter().zip(&b.values) {
        assert_eq!(x.to_bits(), y.to_bits(), "{n}");
    }
    assert_eq!(a.problem.data, b.problem.data);
    assert_ne!(a.values, c.values);

    let value = |addr: &str| a.values.iter().find(|(n, _)| n == addr).unwrap().1;
    // TiO2's group draws; SiO2's has every channel off.
    assert_ne!(value("layer[main:1].thickness"), 60.0);
    assert_ne!(value("layer[main:3].thickness"), 60.0);
    assert!((value("layer[main:1].thickness") / 60.0 - 1.0).abs() < 0.1);
    assert_ne!(value("group[TiO2].n_factor"), 1.0);
    assert!((value("group[TiO2].n_factor") - 1.0).abs() < 0.05);
    assert_eq!(value("layer[main:2].thickness"), 95.0);
    assert_eq!(value("group[SiO2].n_factor"), 1.0);

    // The measurements carry the truth plus noise at the stated sigma.
    let mut at_truth = a.problem.clone();
    at_truth.model = a.model.clone();
    at_truth.parameters.iter_mut().for_each(|q| {
        q.fixed = true;
        q.value = None;
    });
    let r = evaluate(&at_truth).unwrap();
    assert!(
        (r.chi2 / r.n_residuals as f64 - 1.0).abs() < 0.3,
        "chi2/m {}",
        r.chi2 / r.n_residuals as f64
    );
}

fn strategies() -> Vec<Strategy> {
    let mut list = generate(
        &[StrategyGen::Orderings {
            sets: vec!["thicknesses".into(), "groups".into()],
            sample: None,
            polish: true,
        }],
        24,
    )
    .unwrap()
    .strategies;
    list.push(
        Strategy::from_json(
            r#"{"name": "all at once", "stages": [{"name": "all", "free": ["*"]}]}"#,
        )
        .unwrap(),
    );
    // The systematic-error model: no per-layer thicknesses, so it cannot
    // express this truth.
    list.push(
        Strategy::from_json(
            r#"{"name": "groups only", "stages": [{"name": "g", "free": ["set:groups"]}]}"#,
        )
        .unwrap(),
    );
    list
}

#[test]
fn a_benchmark_ranks_recovery_and_repeats_itself() {
    let p = problem(true);
    let seeds = [1, 2, 3];
    let a = benchmark(&p, &strategies(), &[], &seeds).unwrap();
    let b = benchmark(&p, &strategies(), &[], &seeds).unwrap();
    assert_eq!(a.n_free, 6);
    assert_eq!(a.entries.len(), 4);
    for (i, (x, y)) in a.entries.iter().zip(&b.entries).enumerate() {
        assert_eq!(x.rank, i + 1);
        assert_eq!(x.name, y.name);
        for (r, s) in x.runs.iter().zip(&y.runs) {
            assert_eq!(r.recovery.to_bits(), s.recovery.to_bits(), "{}", x.name);
            assert_eq!(r.evals, s.evals);
        }
        assert_eq!(x.runs.iter().map(|r| r.seed).collect::<Vec<_>>(), seeds);
    }
    let last = a.entries.last().unwrap();
    assert_eq!(last.name, "groups only");
    let best = &a.entries[0];
    // Recovery is noise-limited for the strategies that can express the
    // truth (thickness and index trade off), and one that cannot fails χ².
    assert!(best.mean_recovery < 0.7 * last.mean_recovery);
    assert_eq!(best.successes, 3, "{best:#?}");
    assert_eq!(last.successes, 0);
    assert!((best.mean_chi2_reduced - 1.0).abs() < 0.3);
    assert!(last.mean_chi2_reduced > 3.0);
    // Ties in recovery go to the cheaper strategy.
    assert!(a.entries[0].mean_evals <= a.entries[1].mean_evals);
}

#[test]
fn a_tournament_ranks_by_chi2_and_puts_failures_last() {
    let truth = draw_truth(&problem(true), 5).unwrap();
    let mut list = strategies();
    list.push(
        Strategy::from_json(
            r#"{"name": "broken", "stages": [{"name": "x", "free": ["set:nope"]}]}"#,
        )
        .unwrap(),
    );
    let ranked = tournament(&truth.problem, &list, &[]);
    assert_eq!(ranked.len(), 5);
    let last = ranked.last().unwrap();
    assert_eq!(last.name, "broken");
    assert!(last.error.as_deref().unwrap().contains("no set \"nope\""));
    for w in ranked.windows(2) {
        assert!(w[0].chi2 <= w[1].chi2);
    }
    assert_eq!(ranked[3].name, "groups only");
    let best = &ranked[0];
    assert!(best.evals > 0 && best.rel_uncertainty > 0.0);
    assert_eq!(best.report.as_ref().unwrap().name, best.name);
}
