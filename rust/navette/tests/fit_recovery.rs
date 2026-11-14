// SPDX-License-Identifier: LGPL-3.0-or-later
//! Fit engine, phase 1 (`docs/plans/general_fit_plan.md`): recovery of a
//! known truth through the public API.
//!
//! The contract: truth → forward → seeded noise → fit from a wrong start →
//! every free parameter within 2σ of its truth, σ the report's own
//! standard error. Plus the invariants that make the numbers trustworthy:
//! a fully fixed problem changes nothing, a noiseless fit started at the
//! truth stays there, and the same request gives the same bits.

use std::collections::{BTreeMap, HashMap};

use rand::SeedableRng;
use rand::rngs::StdRng;
use rand_distr::{Distribution, Normal};
use serde_json::json;

use navette::fit::{FitModel, FitProblem, Parameter, SpectrumData, fit};
use navette::lsq::lm::LmTermination;
use navette::smatrix::angle::Radians;
use navette::structure::enums::{BlockKind, LayerType};
use navette::structure::{Architect, Entry, Group, Layer, MaterialSpec, Structure};

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

fn wavelengths() -> Vec<f64> {
    (0..121).map(|i| 400.0 + 5.0 * i as f64).collect()
}

/// air | TiO2 80 | SiO2 120 | TiO2 60 | glass, groups on both films.
fn hlh_model() -> FitModel {
    let layers = vec![
        half_space("air", LayerType::Ambient),
        Layer::film(80.0, "TiO2"),
        Layer::film(120.0, "SiO2"),
        Layer::film(60.0, "TiO2"),
        half_space("glass", LayerType::Substrate),
    ];
    let groups = HashMap::from([
        ("TiO2".to_string(), Group::new("TiO2")),
        ("SiO2".to_string(), Group::new("SiO2")),
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

/// Truth curves for `channels` at `angles_deg`, plus seeded Gaussian noise.
fn measure(model: &FitModel, channels: &[(&str, f64)], sigma: f64, seed: u64) -> Vec<SpectrumData> {
    let template: Vec<SpectrumData> = channels
        .iter()
        .map(|(ch, deg)| {
            SpectrumData::new(
                format!("{ch}@{deg}"),
                *ch,
                Radians::from_degrees(*deg).unwrap(),
                wavelengths(),
                vec![0.0; wavelengths().len()],
                sigma,
            )
        })
        .collect();
    let truth = FitProblem::new(model.clone(), template.clone(), vec![])
        .simulate()
        .unwrap();
    let mut rng = StdRng::seed_from_u64(seed);
    let noise = Normal::new(0.0, sigma).unwrap();
    template
        .into_iter()
        .zip(truth)
        .map(|(mut d, t)| {
            d.values = t.iter().map(|v| v + noise.sample(&mut rng)).collect();
            d
        })
        .collect()
}

fn within_two_sigma(report: &navette::fit::FitReport, truth: &[(&str, f64)]) {
    for (addr, want) in truth {
        let p = report.parameter(addr).unwrap();
        let se = p.std_error.unwrap();
        assert!(
            (p.value - want).abs() <= 2.0 * se,
            "{addr}: fitted {} ± {se}, truth {want} (start {})",
            p.value,
            p.start
        );
        assert!(!p.at_bound, "{addr} ended on a bound");
    }
}

#[test]
fn three_layer_stack_recovers_thickness_n_factor_and_cauchy_a() {
    let mut truth = hlh_model();
    let truths = [
        ("layer[main:2].thickness", 124.0),
        ("group[SiO2].n_factor", 1.02),
        ("material[TiO2].A", 2.25),
    ];
    for (addr, v) in truths {
        truth
            .set(&navette::fit::ParamAddr::parse(addr).unwrap(), v)
            .unwrap();
    }
    let data = measure(
        &truth,
        &[("R_avg", 8.0), ("Ts", 45.0), ("Rp", 60.0)],
        0.002,
        7,
    );
    let parameters = vec![
        Parameter::parse("layer[main:2].thickness")
            .unwrap()
            .start(115.0),
        Parameter::parse("group[SiO2].n_factor")
            .unwrap()
            .start(0.99),
        Parameter::parse("material[TiO2].A").unwrap().start(2.15),
    ];
    let report = fit(&FitProblem::new(hlh_model(), data, parameters)).unwrap();

    assert!(
        matches!(
            report.termination,
            LmTermination::Gradient | LmTermination::Step | LmTermination::Cost
        ),
        "{:?}",
        report.termination
    );
    // Noise at the stated sigma: chi2_red ~ 1.
    assert!(
        (report.chi2_reduced - 1.0).abs() < 0.2,
        "chi2_red {}",
        report.chi2_reduced
    );
    within_two_sigma(&report, &truths);
    assert_eq!(report.n_free, 3);
    assert_eq!(report.n_residuals, 3 * 121);
    let corr = report.correlation.as_ref().unwrap();
    for i in 0..3 {
        assert!((corr[i * 3 + i] - 1.0).abs() < 1e-12);
    }
    // The fitted model carries the fitted values.
    let a = navette::fit::ParamAddr::parse("material[TiO2].A").unwrap();
    assert_eq!(
        report.model.get(&a).unwrap(),
        report.parameter("material[TiO2].A").unwrap().value
    );
}

#[test]
fn psi_delta_of_a_sio2_film_on_si() {
    let layers = vec![
        half_space("air", LayerType::Ambient),
        Layer::film(100.0, "SiO2"),
        half_space("Si", LayerType::Substrate),
    ];
    let materials = HashMap::from([
        ("air".to_string(), spec("Konstant", json!({"n": 1.0}))),
        (
            "Si".to_string(),
            spec("Konstant", json!({"n": 3.88, "k": 0.02})),
        ),
        (
            "SiO2".to_string(),
            spec("Cauchy", json!({"A": 1.45, "B": 0.0036, "C": 0.0})),
        ),
    ]);
    let model =
        FitModel::from_structure(&Structure::new(layers, HashMap::new()), "film", materials)
            .unwrap();
    let mut truth = model.clone();
    let truths = [
        ("layer[film:#SiO2_1].thickness", 103.5),
        ("material[SiO2].A", 1.462),
    ];
    for (addr, v) in truths {
        truth
            .set(&navette::fit::ParamAddr::parse(addr).unwrap(), v)
            .unwrap();
    }
    let channels: Vec<(&str, f64)> = [65.0, 70.0, 75.0]
        .iter()
        .flat_map(|&a| [("Psi_R", a), ("Delta_R", a)])
        .collect();
    let data = measure(&truth, &channels, 0.05f64.to_radians(), 11);
    let parameters = vec![
        Parameter::parse("layer[film:#SiO2_1].thickness").unwrap(),
        Parameter::parse("material[SiO2].A").unwrap(),
    ];
    let report = fit(&FitProblem::new(model, data, parameters)).unwrap();
    within_two_sigma(&report, &truths);
    assert!(
        (report.chi2_reduced - 1.0).abs() < 0.2,
        "{}",
        report.chi2_reduced
    );
    for d in &report.datasets {
        assert!(
            d.rms < 3.0 * 0.05f64.to_radians(),
            "{}: rms {}",
            d.id,
            d.rms
        );
    }
}

#[test]
fn every_parameter_fixed_changes_nothing() {
    let model = hlh_model();
    let data = measure(&model, &[("R_avg", 8.0)], 0.002, 1);
    let parameters = vec![Parameter::parse("layer[main:1].thickness").unwrap().fixed()];
    let report = fit(&FitProblem::new(model.clone(), data, parameters)).unwrap();
    assert_eq!(report.termination, LmTermination::NoFreeParameters);
    assert_eq!(report.model, model);
    assert_eq!(report.parameters[0].value, 80.0);
    assert!(report.covariance.is_none());
}

#[test]
fn a_noiseless_fit_started_at_the_truth_stays_there() {
    let model = hlh_model();
    let mut data = measure(&model, &[("R_avg", 8.0), ("Tp", 50.0)], 0.002, 1);
    let clean = FitProblem::new(model.clone(), data.clone(), vec![])
        .simulate()
        .unwrap();
    for (d, c) in data.iter_mut().zip(clean) {
        d.values = c;
    }
    let parameters = vec![
        Parameter::parse("layer[main:2].thickness").unwrap(),
        Parameter::parse("material[TiO2].A").unwrap(),
    ];
    let report = fit(&FitProblem::new(model, data, parameters)).unwrap();
    assert_eq!(report.chi2, 0.0);
    assert_eq!(
        report.parameter("layer[main:2].thickness").unwrap().value,
        120.0
    );
    assert_eq!(report.parameter("material[TiO2].A").unwrap().value, 2.2);
}

#[test]
fn the_same_problem_gives_the_same_bits() {
    let model = hlh_model();
    let data = measure(&model, &[("R_avg", 8.0), ("Ts", 45.0)], 0.002, 3);
    let parameters = vec![
        Parameter::parse("layer[main:1].thickness")
            .unwrap()
            .start(85.0),
        Parameter::parse("group[TiO2].n_factor").unwrap().log(),
    ];
    let problem = FitProblem::new(model, data, parameters);
    let a = fit(&problem).unwrap();
    let b = fit(&problem).unwrap();
    for (p, q) in a.parameters.iter().zip(&b.parameters) {
        assert_eq!(p.value.to_bits(), q.value.to_bits(), "{}", p.addr);
    }
    assert_eq!(a.chi2.to_bits(), b.chi2.to_bits());
}

#[test]
fn a_repeated_shared_structure_is_one_address() {
    let pair = Structure::new(
        vec![Layer::film(100.0, "TiO2"), Layer::film(150.0, "SiO2")],
        HashMap::new(),
    );
    // Unmarked STACK ends: a marked one must hold both half-spaces.
    let mut arch = Architect::new();
    arch.add_structure(
        Structure::new(vec![Layer::film(20.0, "MgF2")], HashMap::new()),
        false,
        1,
        "top",
        BlockKind::Stack,
    )
    .unwrap();
    arch.add_structure(pair, true, 3, "pair", BlockKind::Films)
        .unwrap();
    arch.add_structure(
        Structure::new(vec![Layer::film(30.0, "Al2O3")], HashMap::new()),
        false,
        1,
        "bottom",
        BlockKind::Stack,
    )
    .unwrap();
    let model = FitModel::from_architect(&arch, HashMap::new()).unwrap();
    assert_eq!(
        model.entries("default").unwrap(),
        arch.iter_entries().as_slice()
    );
    assert_eq!(model.structure_names(), ["top", "pair", "bottom"]);

    let mut m = model.clone();
    let addr = navette::fit::ParamAddr::parse("layer[pair:#SiO2_1].thickness").unwrap();
    m.set(&addr, 140.0).unwrap();
    let sio2: Vec<f64> = m
        .entries("default")
        .unwrap()
        .iter()
        .filter(|(l, _)| l.material == "SiO2")
        .map(|(l, _)| l.thickness)
        .collect();
    assert_eq!(sio2, [140.0; 3]);
}

#[test]
fn guard_rails_name_the_address() {
    let model = hlh_model();
    let data = measure(&model, &[("R_avg", 8.0)], 0.002, 1);
    let refuse =
        |p: Vec<Parameter>| fit(&FitProblem::new(model.clone(), data.clone(), p)).unwrap_err();
    let e = refuse(vec![Parameter::parse("layer[main:9].thickness").unwrap()]);
    assert!(e.contains("layer[main:9]") && e.contains("outside"), "{e}");
    let e = refuse(vec![Parameter::parse("layer[stack:1].thickness").unwrap()]);
    assert!(e.contains("no structure \"stack\""), "{e}");
    let e = refuse(vec![Parameter::parse("material[TiO2].D").unwrap()]);
    assert!(
        e.contains("material[TiO2].D") && e.contains("no parameter"),
        "{e}"
    );
    let e = refuse(vec![Parameter::parse("group[Ta2O5].n_factor").unwrap()]);
    assert!(e.contains("no group"), "{e}");
    let e = refuse(vec![
        Parameter::parse("layer[main:1].gradient.f_end").unwrap(),
    ]);
    assert!(e.contains("no gradient"), "{e}");
    let e = refuse(vec![
        Parameter::parse("layer[main:1].thickness").unwrap(),
        Parameter::parse("layer[main:#TiO2_1].thickness").unwrap(),
    ]);
    assert!(e.contains("named twice"), "{e}");
    let e = refuse(vec![
        Parameter::parse("layer[main:1].thickness")
            .unwrap()
            .start(500.0)
            .bounds(0.0, 200.0),
    ]);
    assert!(e.contains("outside [0, 200]"), "{e}");
    let e = refuse(vec![
        Parameter::parse("group[TiO2].thick_summand").unwrap().log(),
    ]);
    assert!(e.contains("positive start"), "{e}");
}
