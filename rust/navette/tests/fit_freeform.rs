// SPDX-License-Identifier: LGPL-3.0-or-later
//! Fit engine, phase 7 (`docs/plans/general_fit_plan.md` §4.6, §3.5.4):
//! free-form optical constants.
//!
//! A film of unknown material is measured as R and T at a few angles, and
//! its n, k are recovered as a spline on coarse knots, point by point with
//! a smoothness prior, and as a Kramers-Kronig-consistent B-spline in ε₂.

use std::collections::{BTreeMap, HashMap};

use rand::SeedableRng;
use rand::rngs::StdRng;
use rand_distr::{Distribution, Normal};
use serde_json::{Value, json};

use navette::fit::{
    FitModel, FitProblem, FitReport, FitSection, ParamAddr, Parameter, SMOOTH_ID, Smoothness,
    SpectrumData, fit,
};
use navette::smatrix::angle::Radians;
use navette::structure::enums::LayerType;
use navette::structure::{Entry, Group, Layer, MaterialSpec, Structure};

fn spec(model: &str, params: Value) -> MaterialSpec {
    let params: BTreeMap<String, Value> = serde_json::from_value(params).unwrap();
    MaterialSpec::new(model, params)
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

/// air | `d` nm of `film` | glass.
fn model(film: MaterialSpec, d: f64) -> FitModel {
    let materials = HashMap::from([
        (
            "air".to_string(),
            Entry::Spec(spec("Konstant", json!({"n": 1.0}))),
        ),
        (
            "glass".to_string(),
            Entry::Spec(spec("Cauchy", json!({"A": 1.50, "B": 0.004, "C": 0.0}))),
        ),
        ("film".to_string(), Entry::Spec(film)),
    ]);
    let layers = vec![
        half_space("air", LayerType::Ambient),
        Layer::film(d, "film"),
        half_space("glass", LayerType::Substrate),
    ];
    let groups = HashMap::from([("film".to_string(), Group::new("film"))]);
    FitModel::from_structure(&Structure::new(layers, groups), "main", materials).unwrap()
}

/// R at 8°, T at 0°, Rs and Rp at 60°, each with seeded noise.
fn measure(truth: &FitModel, sigma: f64, seed: u64) -> Vec<SpectrumData> {
    let curves = [
        ("R", "R_avg", 8.0),
        ("T", "T_avg", 0.0),
        ("Rs", "Rs", 60.0),
        ("Rp", "Rp", 60.0),
    ];
    let template: Vec<SpectrumData> = curves
        .iter()
        .map(|&(id, ch, deg)| {
            SpectrumData::new(
                id,
                ch,
                Radians::from_degrees(deg).unwrap(),
                wavelengths(),
                vec![0.0; wavelengths().len()],
                sigma,
            )
        })
        .collect();
    let clean = FitProblem::new(truth.clone(), template.clone(), vec![])
        .simulate()
        .unwrap();
    let mut rng = StdRng::seed_from_u64(seed);
    let noise = Normal::new(0.0, sigma).unwrap();
    template
        .into_iter()
        .zip(clean)
        .map(|(mut d, t)| {
            d.values = t.iter().map(|v| v + noise.sample(&mut rng)).collect();
            d
        })
        .collect()
}

fn params(addrs: &[String]) -> Vec<Parameter> {
    addrs.iter().map(|a| Parameter::parse(a).unwrap()).collect()
}

fn knots(param: &str, n: usize) -> Vec<String> {
    (0..n)
        .map(|i| format!("material[film].{param}[{i}]"))
        .collect()
}

fn within(report: &FitReport, addr: &str, want: f64, sigmas: f64) {
    let p = report.parameter(addr).unwrap();
    let se = p.std_error.unwrap();
    assert!(
        (p.value - want).abs() <= sigmas * se,
        "{addr}: fitted {} +- {se}, truth {want} (start {})",
        p.value,
        p.start
    );
}

fn chi2_red_near_one(report: &FitReport) {
    assert!(
        (report.chi2_reduced - 1.0).abs() < 0.3,
        "chi2_red {}",
        report.chi2_reduced
    );
}

/// A weakly absorbing high-index film: n and k at `wl`.
fn truth_nk(wl: f64) -> (f64, f64) {
    let um = wl / 1000.0;
    (
        2.05 + 0.025 / (um * um),
        0.002 + 0.03 * (-(wl - 400.0) / 80.0).exp(),
    )
}

#[test]
fn a_spline_film_is_recovered_knot_by_knot() {
    let x: Vec<f64> = (0..7).map(|i| 400.0 + 100.0 * i as f64).collect();
    let (n, k): (Vec<f64>, Vec<f64>) = x.iter().map(|&w| truth_nk(w)).unzip();
    let truth = model(
        spec("Spline", json!({"wavelengths": x, "n": n, "k": k})),
        420.0,
    );
    let data = measure(&truth, 0.001, 1);
    let flat = spec(
        "Spline",
        json!({"wavelengths": x, "n": vec![2.1; 7], "k": vec![0.005; 7]}),
    );
    let mut addrs = knots("n", 7);
    addrs.extend(knots("k", 7));
    addrs.push("layer[main:1].thickness".into());
    let problem = FitProblem::new(model(flat, 410.0), data, params(&addrs));
    let r = fit(&problem).unwrap();
    chi2_red_near_one(&r);
    for i in 0..7 {
        within(&r, &format!("material[film].n[{i}]"), n[i], 3.0);
        within(&r, &format!("material[film].k[{i}]"), k[i], 3.0);
    }
    within(&r, "layer[main:1].thickness", 420.0, 3.0);
}

#[test]
fn a_point_wise_film_is_recovered_with_a_smoothness_prior() {
    let x: Vec<f64> = (0..25).map(|i| 400.0 + 25.0 * i as f64).collect();
    let (n, k): (Vec<f64>, Vec<f64>) = x.iter().map(|&w| truth_nk(w)).unzip();
    let truth = model(
        spec("PointWise", json!({"wavelengths": x, "n": n, "k": k})),
        420.0,
    );
    let data = measure(&truth, 0.001, 2);
    let start = spec(
        "PointWise",
        json!({"wavelengths": x, "n": vec![2.1; 25], "k": vec![0.005; 25]}),
    );
    let mut addrs = knots("n", 25);
    addrs.extend(knots("k", 25));
    addrs.push("layer[main:1].thickness".into());
    let mut problem = FitProblem::new(model(start, 410.0), data, params(&addrs));
    problem.smoothness = vec![
        Smoothness::new("film", "n", 0.005),
        Smoothness::new("film", "k", 0.002),
    ];
    let r = fit(&problem).unwrap();
    chi2_red_near_one(&r);
    let smooth = r.datasets.iter().find(|d| d.id == SMOOTH_ID).unwrap();
    assert_eq!(smooth.n_points, 2 * 23);
    let worst = |param: &str, truth: &[f64]| {
        (0..25)
            .map(|i| {
                let p = r
                    .parameter(&format!("material[film].{param}[{i}]"))
                    .unwrap();
                (p.value - truth[i]).abs()
            })
            .fold(0.0, f64::max)
    };
    assert!(worst("n", &n) < 0.02, "n off by {}", worst("n", &n));
    assert!(worst("k", &k) < 0.01, "k off by {}", worst("k", &k));
    within(&r, "layer[main:1].thickness", 420.0, 3.0);
}

#[test]
fn a_kk_b_spline_film_is_recovered() {
    // An absorption band near 3.5 eV whose tail reaches into the data.
    let e = [1.0, 2.0, 2.75, 3.25, 3.75, 4.5, 6.0, 8.0];
    let c = [0.0, 0.0, 0.15, 0.9, 2.5, 3.0, 1.0, 0.0];
    let film = |c: &[f64], eps_inf: f64| {
        spec(
            "BSplineEps",
            json!({"energies": e, "eps2": c, "epsilon_inf": eps_inf}),
        )
    };
    let truth = model(film(&c, 3.2), 420.0);
    let data = measure(&truth, 0.001, 3);
    let start = model(film(&[0.0, 0.0, 0.3, 0.5, 2.0, 3.0, 1.0, 0.0], 3.0), 410.0);
    let addrs: Vec<String> = [
        "material[film].eps2[2]",
        "material[film].eps2[3]",
        "material[film].eps2[4]",
        "material[film].epsilon_inf",
        "layer[main:1].thickness",
    ]
    .map(String::from)
    .to_vec();
    let r = fit(&FitProblem::new(start, data, params(&addrs))).unwrap();
    chi2_red_near_one(&r);
    within(&r, "material[film].eps2[2]", c[2], 3.0);
    within(&r, "material[film].eps2[3]", c[3], 3.0);
    within(&r, "material[film].eps2[4]", c[4], 3.0);
    within(&r, "material[film].epsilon_inf", 3.2, 3.0);
    within(&r, "layer[main:1].thickness", 420.0, 3.0);
}

#[test]
fn knots_are_addressable_bounded_and_smoothness_is_checked() {
    let x = [400.0, 600.0, 800.0];
    let m = model(
        spec(
            "Spline",
            json!({"wavelengths": x, "n": vec![2.0; 3], "k": vec![0.0; 3]}),
        ),
        100.0,
    );
    let all: Vec<String> = m.addresses().iter().map(|a| a.to_string()).collect();
    for i in 0..3 {
        assert!(all.contains(&format!("material[film].n[{i}]")), "{all:?}");
        assert!(all.contains(&format!("material[film].k[{i}]")));
    }
    // The knot grid is not a parameter.
    assert!(!all.iter().any(|a| a.contains("wavelengths")));
    let addr = ParamAddr::parse("material[film].k[1]").unwrap();
    assert_eq!(m.get(&addr).unwrap(), 0.0);

    let data = measure(&m, 0.01, 4);
    let mut p = FitProblem::new(
        m.clone(),
        data.clone(),
        params(&["material[film].k[1]".to_string()]),
    );
    let r = fit(&p).unwrap();
    assert_eq!(r.parameters[0].lb, 0.0, "k is non-negative per knot");

    p.smoothness = vec![Smoothness::new("film", "wavelengths", 1.0)];
    assert!(fit(&p).unwrap_err().contains("not a free-form list"));
    p.smoothness = vec![Smoothness::new("film", "n", 0.0)];
    assert!(fit(&p).unwrap_err().contains("sigma"));
    p.smoothness = vec![Smoothness::new("glass", "n", 1.0)];
    assert!(fit(&p).is_err());

    let section: FitSection = serde_json::from_value(json!({
        "parameters": [{"addr": "material[film].n[1]"}],
        "smoothness": [{"material": "film", "param": "n", "sigma": 0.01}]
    }))
    .unwrap();
    assert_eq!(section.smoothness, [Smoothness::new("film", "n", 0.01)]);
}
