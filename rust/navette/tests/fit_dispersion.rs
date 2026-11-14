// SPDX-License-Identifier: LGPL-3.0-or-later
//! Fit engine, phase 2 (`docs/plans/general_fit_plan.md` §4): dispersion
//! models against measured optical constants.
//!
//! The contract is phase 1's: truth → `MaterialSpec::evaluate` → seeded
//! noise → fit from a wrong start → every free parameter within 2σ of its
//! truth. Plus a Bruggeman fraction whose finite-difference derivative is
//! stable across step sizes (Newton's tolerance must not leak into the
//! Jacobian), a joint fit of stack and nk data in one residual vector, and
//! a model ranking of a table.

use std::collections::{BTreeMap, HashMap};

use num_complex::Complex64;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand_distr::{Distribution, Normal};
use serde_json::json;

use navette::fit::{
    FitModel, FitProblem, FitReport, NkColumn, NkData, NkQuantity, Parameter, SpectrumData, fit,
    fit_dispersion, table_to_model,
};
use navette::lsq::lm::LmConfig;
use navette::smatrix::angle::Radians;
use navette::structure::enums::LayerType;
use navette::structure::{Entry, Layer, MaterialSpec, Structure};

fn spec(model: &str, params: serde_json::Value) -> MaterialSpec {
    let params: BTreeMap<String, serde_json::Value> = serde_json::from_value(params).unwrap();
    MaterialSpec::new(model, params)
}

fn grid(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| lo + (hi - lo) * i as f64 / (n - 1) as f64)
        .collect()
}

/// `truth` on `wl`, seeded noise on each requested quantity.
fn measure(
    truth: &MaterialSpec,
    wl: &[f64],
    quantities: &[(NkQuantity, f64)],
    seed: u64,
) -> NkData {
    let nk = truth.evaluate(wl).unwrap();
    let mut rng = StdRng::seed_from_u64(seed);
    let columns = quantities
        .iter()
        .map(|&(q, sigma)| {
            let noise = Normal::new(0.0, sigma).unwrap();
            let values = nk
                .iter()
                .map(|&z| q.of(z) + noise.sample(&mut rng))
                .collect();
            NkColumn::new(q, values, sigma)
        })
        .collect();
    NkData::new("measured", "film", wl.to_vec(), columns)
}

fn within_two_sigma(report: &FitReport, truth: &[(&str, f64)]) {
    for (addr, want) in truth {
        let p = report
            .parameter(addr)
            .unwrap_or_else(|| panic!("{addr} not in the report"));
        let se = p
            .std_error
            .unwrap_or_else(|| panic!("{addr}: no standard error ({:?})", report.covariance_note));
        assert!(
            (p.value - want).abs() <= 2.0 * se,
            "{addr}: fitted {} ± {se}, truth {want} (start {}; {:?} after {} iterations, chi2_red {})",
            p.value,
            p.start,
            report.termination,
            report.iterations,
            report.chi2_reduced
        );
        assert!(!p.at_bound, "{addr} ended on a bound");
    }
    assert!(
        (0.6..1.5).contains(&report.chi2_reduced),
        "chi2_red {}",
        report.chi2_reduced
    );
}

#[test]
fn cauchy_recovers_its_own_data() {
    let truth = spec("Cauchy", json!({"A": 1.46, "B": 0.00354, "C": 0.0001}));
    let start = spec("Cauchy", json!({"A": 1.5, "B": 0.003, "C": 0.0}));
    let data = measure(
        &truth,
        &grid(350.0, 1100.0, 151),
        &[(NkQuantity::N, 1e-4)],
        11,
    );
    let (fitted, report) = fit_dispersion(&start, vec![data], None, &LmConfig::default()).unwrap();
    assert_eq!(report.n_free, 3);
    within_two_sigma(
        &report,
        &[
            ("material[film].A", 1.46),
            ("material[film].B", 0.00354),
            ("material[film].C", 0.0001),
        ],
    );
    assert_eq!(
        fitted.params["A"].as_f64(),
        Some(report.parameters[0].value)
    );
}

#[test]
fn sellmeier_recovers_its_own_data() {
    // Fused silica (Malitson), C in µm².
    let truth = spec(
        "Sellmeier",
        json!({"B1": 0.6961663, "C1": 0.0046791, "B2": 0.4079426, "C2": 0.0135121,
               "B3": 0.8974794, "C3": 97.934}),
    );
    let data = measure(
        &truth,
        &grid(250.0, 2500.0, 226),
        &[(NkQuantity::N, 1e-5)],
        2,
    );
    // The two UV poles trade off against each other almost exactly; the
    // second is held at its truth and the other five are fitted.
    let start = spec(
        "Sellmeier",
        json!({"B1": 0.68, "C1": 0.0048, "B2": 0.42, "C2": 0.0135121, "B3": 0.88, "C3": 100.0}),
    );
    let params = ["B1", "C1", "B2", "B3", "C3"]
        .iter()
        .map(|k| Parameter::parse(&format!("material[film].{k}")).unwrap())
        .collect();
    let (_, report) =
        fit_dispersion(&start, vec![data], Some(params), &LmConfig::default()).unwrap();
    within_two_sigma(
        &report,
        &[
            ("material[film].B1", 0.6961663),
            ("material[film].C1", 0.0046791),
            ("material[film].B2", 0.4079426),
            ("material[film].B3", 0.8974794),
            ("material[film].C3", 97.934),
        ],
    );
}

#[test]
fn tauc_lorentz_recovers_its_own_data() {
    let truth = spec(
        "TaucLorentz",
        json!({"Eg": 3.2, "osc": [[120.0, 4.4, 1.2]], "epsilon_inf": 2.1}),
    );
    let start = spec(
        "TaucLorentz",
        json!({"Eg": 3.1, "osc": [[100.0, 4.6, 1.0]], "epsilon_inf": 2.0}),
    );
    let data = measure(
        &truth,
        &grid(250.0, 1000.0, 151),
        &[(NkQuantity::N, 1e-3), (NkQuantity::K, 1e-3)],
        3,
    );
    let (_, report) = fit_dispersion(&start, vec![data], None, &LmConfig::default()).unwrap();
    assert_eq!(report.n_free, 5);
    within_two_sigma(
        &report,
        &[
            ("material[film].Eg", 3.2),
            ("material[film].osc[0][0]", 120.0),
            ("material[film].osc[0][1]", 4.4),
            ("material[film].osc[0][2]", 1.2),
            ("material[film].epsilon_inf", 2.1),
        ],
    );
}

fn porous(fraction: f64) -> MaterialSpec {
    spec(
        "Bruggeman",
        json!({
            "host": {"model": "Cauchy", "params": {"A": 2.3, "B": 0.02, "C": 0.0}},
            "inclusion": {"model": "Konstant", "params": {"n": 1.0}},
            "fraction": fraction,
        }),
    )
}

#[test]
fn bruggeman_fraction_recovers_and_constituents_stay_held() {
    let data = measure(
        &porous(0.3),
        &grid(400.0, 1000.0, 121),
        &[(NkQuantity::N, 1e-4)],
        4,
    );
    let (fitted, report) =
        fit_dispersion(&porous(0.1), vec![data], None, &LmConfig::default()).unwrap();
    // A mixer's default free parameters are its own: the fraction only.
    assert_eq!(report.n_free, 1);
    within_two_sigma(&report, &[("material[film].fraction", 0.3)]);
    assert_eq!(fitted.params["host"], porous(0.3).params["host"]);
    let p = report.parameter("material[film].fraction").unwrap();
    assert_eq!((p.lb, p.ub), (0.0, 1.0));
}

/// `dn/df` of the closed-form two-phase Bruggeman (lossless constituents):
/// ε = (b + sqrt(b² + 8 εi εh)) / 4, b = (3f − 1) εi + (2 − 3f) εh.
fn bruggeman_dn_df(ei: f64, eh: f64, f: f64) -> f64 {
    let b = (3.0 * f - 1.0) * ei + (2.0 - 3.0 * f) * eh;
    let s = (b * b + 8.0 * ei * eh).sqrt();
    let eps = (b + s) / 4.0;
    let db = 3.0 * (ei - eh);
    let deps = (db + b * db / s) / 4.0;
    deps / (2.0 * eps.sqrt())
}

#[test]
fn bruggeman_fd_derivative_is_stable_across_steps() {
    // Newton stops at |Δε| <= 1e-9: the derivative must not see that.
    let wl = [450.0, 633.0, 900.0];
    let f = 0.3;
    let host = spec("Cauchy", json!({"A": 2.3, "B": 0.02, "C": 0.0}))
        .evaluate(&wl)
        .unwrap();
    for (j, &w) in wl.iter().enumerate() {
        let want = bruggeman_dn_df(1.0, host[j].re * host[j].re, f);
        for h in [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8] {
            let up = porous(f + h).evaluate(&[w]).unwrap()[0].re;
            let dn = porous(f - h).evaluate(&[w]).unwrap()[0].re;
            let got = (up - dn) / (2.0 * h);
            let tol = if h >= 1e-6 { 1e-6 } else { 1e-5 };
            assert!(
                ((got - want) / want).abs() < tol,
                "{w} nm, h = {h}: FD {got}, closed form {want}"
            );
        }
    }
}

#[test]
fn eps_columns_recover_like_nk() {
    let truth = spec(
        "TaucLorentz",
        json!({"Eg": 3.2, "osc": [[120.0, 4.4, 1.2]], "epsilon_inf": 2.1}),
    );
    let start = spec(
        "TaucLorentz",
        json!({"Eg": 3.1, "osc": [[100.0, 4.6, 1.0]], "epsilon_inf": 2.0}),
    );
    let data = measure(
        &truth,
        &grid(250.0, 1000.0, 151),
        &[(NkQuantity::Eps1, 5e-3), (NkQuantity::Eps2, 5e-3)],
        15,
    );
    let (_, report) = fit_dispersion(&start, vec![data], None, &LmConfig::default()).unwrap();
    within_two_sigma(
        &report,
        &[
            ("material[film].Eg", 3.2),
            ("material[film].osc[0][1]", 4.4),
        ],
    );
}

#[test]
fn joint_fit_of_a_spectrum_and_an_nk_table() {
    // Ambient | film 150 nm | substrate; Cauchy film. The film's A comes
    // from the nk table, its thickness from Rs: one residual vector.
    let film = spec("Cauchy", json!({"A": 1.9, "B": 0.015, "C": 0.0}));
    let materials = HashMap::from([
        (
            "air".to_string(),
            Entry::Spec(spec("Konstant", json!({"n": 1.0}))),
        ),
        (
            "glass".to_string(),
            Entry::Spec(spec("Konstant", json!({"n": 1.52}))),
        ),
        ("film".to_string(), Entry::Spec(film.clone())),
    ]);
    let layers = vec![
        Layer {
            layer_type: LayerType::Ambient,
            ..Layer::film(0.0, "air")
        },
        Layer::film(150.0, "film"),
        Layer {
            layer_type: LayerType::Substrate,
            ..Layer::film(0.0, "glass")
        },
    ];
    let model =
        FitModel::from_structure(&Structure::new(layers, HashMap::new()), "coat", materials)
            .unwrap();

    let wl = grid(400.0, 1000.0, 121);
    let truth_sim = FitProblem::new(
        model.clone(),
        vec![SpectrumData::new(
            "R",
            "Rs",
            Radians::ZERO,
            wl.clone(),
            vec![0.0; 121],
            1.0,
        )],
        vec![],
    )
    .simulate()
    .unwrap();
    let mut rng = StdRng::seed_from_u64(6);
    let noise = Normal::new(0.0, 1e-4).unwrap();
    let r: Vec<f64> = truth_sim[0]
        .iter()
        .map(|v| v + noise.sample(&mut rng))
        .collect();
    let spectrum = SpectrumData::new("R", "Rs", Radians::ZERO, wl.clone(), r, 1e-4);
    let table = measure(&film, &wl, &[(NkQuantity::N, 1e-3)], 7);

    let mut problem = FitProblem::new(
        model,
        vec![spectrum],
        vec![
            Parameter::parse("layer[coat:1].thickness")
                .unwrap()
                .start(140.0),
            Parameter::parse("material[film].A").unwrap().start(1.85),
        ],
    );
    problem.nk_data = vec![table];
    let report = fit(&problem).unwrap();
    assert_eq!(report.datasets.len(), 2);
    assert_eq!(report.datasets[1].id, "measured");
    assert_eq!(report.n_residuals, 242);
    within_two_sigma(
        &report,
        &[
            ("layer[coat:1].thickness", 150.0),
            ("material[film].A", 1.9),
        ],
    );
}

#[test]
fn table_to_model_ranks_the_generating_model_first() {
    // A table sampled from a Sellmeier: Sellmeier beats Cauchy and Konstant.
    let source = spec(
        "Sellmeier",
        json!({"B1": 0.6961663, "C1": 0.0046791, "B2": 0.4079426, "C2": 0.0135121,
               "B3": 0.8974794, "C3": 97.934}),
    );
    let wl = grid(220.0, 2400.0, 110);
    let nk: Vec<Complex64> = source.evaluate(&wl).unwrap();
    let table = spec(
        "Table",
        json!({"n_data": [wl, nk.iter().map(|z| z.re).collect::<Vec<_>>()]}),
    );
    let candidates = vec![
        spec("Konstant", json!({"n": 1.5})),
        spec("Cauchy", json!({"A": 1.45, "B": 0.004, "C": 0.0})),
        spec(
            "Sellmeier",
            json!({"B1": 0.7, "C1": 0.005, "B2": 0.4, "C2": 0.0135121, "B3": 0.9, "C3": 100.0}),
        ),
    ];
    let ranked = table_to_model(&table, candidates, &LmConfig::default()).unwrap();
    let order: Vec<&str> = ranked.iter().map(|c| c.model.as_str()).collect();
    assert_eq!(order, ["Sellmeier", "Cauchy", "Konstant"], "{ranked:#?}");
    assert!(ranked.iter().all(|c| c.error.is_none()));
    assert!(ranked[0].bic < ranked[1].bic && ranked[1].bic < ranked[2].bic);
}

#[test]
fn guard_rails() {
    let start = spec("Cauchy", json!({"A": 1.5, "B": 0.003, "C": 0.0}));
    let lsq = LmConfig::default();
    assert!(
        fit_dispersion(&start, vec![], None, &lsq)
            .unwrap_err()
            .contains("no nk data")
    );

    let wl = grid(400.0, 800.0, 5);
    let a = NkData::new(
        "a",
        "film",
        wl.clone(),
        vec![NkColumn::new(NkQuantity::N, vec![1.5; 5], 1e-3)],
    );
    let mut b = a.clone();
    b.material = "other".into();
    assert!(
        fit_dispersion(&start, vec![a.clone(), b], None, &lsq)
            .unwrap_err()
            .contains("one spec, one material")
    );

    let mut twice = a.clone();
    twice.columns.push(twice.columns[0].clone());
    assert!(twice.check().unwrap_err().contains("given twice"));

    let mut short = a.clone();
    short.columns[0].values.pop();
    assert!(short.check().unwrap_err().contains("4 N values for 5"));

    // A start outside the per-model bounds names the parameter.
    let bad = spec(
        "Bruggeman",
        json!({
            "host": {"model": "Konstant", "params": {"n": 2.0}},
            "inclusion": {"model": "Konstant", "params": {"n": 1.0}},
            "fraction": 1.2,
        }),
    );
    let err = fit_dispersion(&bad, vec![a.clone()], None, &lsq).unwrap_err();
    assert!(
        err.contains("material[film].fraction") && err.contains("outside"),
        "{err}"
    );

    assert!(
        table_to_model(&start, vec![start.clone()], &lsq)
            .unwrap_err()
            .contains("Table spec")
    );
}
