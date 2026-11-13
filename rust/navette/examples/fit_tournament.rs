// SPDX-License-Identifier: LGPL-3.0-or-later
//! Which fit strategy recovers a 20-layer coating best?
//! (`docs/plans/general_fit_plan.md` §3.8.2, phase 4.)
//!
//! A 20-layer TiO2/SiO2 edge filter whose groups carry deposition error
//! laws (2 % thickness, 1 % index, relative). For each seed a truth is
//! drawn from those laws, R and T are simulated with noise, and every
//! strategy (hand-written and generated from parameter sets) is fitted
//! from the nominal design. The table ranks them by how well they recover
//! the truth. The benchmark runs twice and must give the same bits.
//!
//! ```text
//! cargo run --release -p navette --example fit_tournament [K]
//! ```

use std::collections::{BTreeMap, HashMap};
use std::process::ExitCode;
use std::time::Instant;

use serde_json::json;

use navette::fit::{
    BenchmarkReport, DEFAULT_ORDERING_CAP, FitModel, FitProblem, Parameter, SpectrumData, Strategy,
    StrategyGen, benchmark, generate,
};
use navette::smatrix::angle::Radians;
use navette::structure::enums::{ErrorMask, LayerType};
use navette::structure::{Entry, ErrorParams, Group, Layer, MaterialSpec, Structure};

const LAYERS: usize = 20;

fn spec(model: &str, params: serde_json::Value) -> Entry {
    let params: BTreeMap<String, serde_json::Value> = serde_json::from_value(params).unwrap();
    Entry::Spec(MaterialSpec::new(model, params))
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

/// air | (H L)^10 | glass, quarter waves at 700 nm stretched ±20 % along
/// the stack (a chirped edge), so no two layers are alike.
fn coating() -> FitModel {
    let mut layers = vec![Layer {
        layer_type: LayerType::Ambient,
        ..Layer::film(0.0, "air")
    }];
    for i in 0..LAYERS {
        let (material, n) = if i % 2 == 0 {
            ("TiO2", 2.3)
        } else {
            ("SiO2", 1.46)
        };
        let stretch = 0.8 + 0.4 * i as f64 / (LAYERS - 1) as f64;
        let d = (700.0 / (4.0 * n) * stretch * 10.0_f64).round() / 10.0;
        layers.push(Layer::film(d, material));
    }
    layers.push(Layer {
        layer_type: LayerType::Substrate,
        ..Layer::film(0.0, "glass")
    });
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

fn problem() -> FitProblem {
    let wl: Vec<f64> = (0..151).map(|i| 400.0 + 4.0 * i as f64).collect();
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
    let mut parameters: Vec<Parameter> = (1..=LAYERS)
        .map(|i| Parameter::parse(&format!("layer[main:{i}].thickness")).unwrap())
        .collect();
    parameters.push(Parameter::parse("group[TiO2].n_factor").unwrap());
    parameters.push(Parameter::parse("group[SiO2].n_factor").unwrap());
    FitProblem::new(
        coating(),
        vec![
            curve("R@8", "R_avg", 8.0),
            curve("T@0", "T_avg", 0.0),
            curve("Rp@45", "Rp", 45.0),
        ],
        parameters,
    )
}

fn strategies() -> (Vec<Strategy>, Vec<String>) {
    let explicit: Vec<Strategy> = serde_json::from_str(
        r#"[
        {"name": "all at once", "stages": [{"name": "all", "free": ["*"]}]},
        {"name": "groups only", "stages": [{"name": "g", "free": ["set:groups"]}]},
        {"name": "thickness window, then all", "stages": [
            {"name": "t", "free": ["set:thicknesses"], "window": [700, 1000]},
            {"name": "all", "free": ["*"]}]}
    ]"#,
    )
    .unwrap();
    let sets = || vec!["thicknesses".to_string(), "groups".to_string()];
    let made = generate(
        &[
            StrategyGen::Explicit(explicit),
            StrategyGen::Orderings {
                sets: sets(),
                sample: None,
                polish: true,
            },
            StrategyGen::Cumulative {
                sets: sets(),
                sample: None,
            },
        ],
        DEFAULT_ORDERING_CAP,
    )
    .unwrap();
    (made.strategies, made.notes)
}

fn print(report: &BenchmarkReport) {
    println!(
        "{:>4}  {:<30} {:>10} {:>10} {:>9} {:>6} {:>8} {:>9}",
        "rank", "strategy", "recovery", "worst", "chi2_red", "in 2sd", "failed", "evals"
    );
    for e in &report.entries {
        println!(
            "{:>4}  {:<30} {:>10.3e} {:>10.3e} {:>9.3} {:>4}/{} {:>8} {:>9.0}",
            e.rank,
            e.name,
            e.mean_recovery,
            e.worst_recovery,
            e.mean_chi2_reduced,
            e.successes,
            report.seeds.len(),
            e.failures,
            e.mean_evals
        );
    }
}

fn main() -> ExitCode {
    let k: u64 = std::env::args()
        .nth(1)
        .map(|a| a.parse().expect("K is a count of seeds"))
        .unwrap_or(4);
    let seeds: Vec<u64> = (1..=k).collect();
    let problem = problem();
    let (strategies, notes) = strategies();
    notes.iter().for_each(|n| println!("note: {n}"));
    println!(
        "{LAYERS}-layer coating, {} free parameters, {} strategies, seeds {seeds:?}\n",
        problem.parameters.len(),
        strategies.len()
    );

    let t0 = Instant::now();
    let a = benchmark(&problem, &strategies, &[], &seeds).unwrap();
    let took = t0.elapsed();
    print(&a);
    println!("\n{:.2} s", took.as_secs_f64());

    let b = benchmark(&problem, &strategies, &[], &seeds).unwrap();
    let same = a.entries.iter().zip(&b.entries).all(|(x, y)| {
        x.name == y.name
            && x.runs
                .iter()
                .zip(&y.runs)
                .all(|(r, s)| r.recovery.to_bits() == s.recovery.to_bits() && r.evals == s.evals)
    });
    if same {
        println!("deterministic: a second run gave the same bits");
        ExitCode::SUCCESS
    } else {
        println!("NOT deterministic: a second run differed");
        ExitCode::FAILURE
    }
}
