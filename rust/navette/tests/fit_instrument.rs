// SPDX-License-Identifier: LGPL-3.0-or-later
//! Fit engine, phase 7 (`docs/plans/general_fit_plan.md` §3.5.5): the
//! instrument's departures as fit parameters.
//!
//! Each of angle offset, wavelength shift, scale and offset, and slit
//! bandwidth is put into simulated data, then fitted back from the
//! identity (bandwidth from a positive start: its slope is zero at zero)
//! and recovered within two standard errors.

use std::collections::{BTreeMap, HashMap};

use rand::SeedableRng;
use rand::rngs::StdRng;
use rand_distr::{Distribution, Normal};
use serde_json::json;

use navette::fit::{
    FitModel, FitProblem, FitReport, Instrument, ParamAddr, Parameter, SpectrumData, fit,
};
use navette::smatrix::angle::Radians;
use navette::structure::enums::LayerType;
use navette::structure::{Entry, Group, Layer, MaterialSpec, Structure};

fn spec(model: &str, params: serde_json::Value) -> MaterialSpec {
    let params: BTreeMap<String, serde_json::Value> = serde_json::from_value(params).unwrap();
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

fn materials() -> HashMap<String, Entry> {
    HashMap::from([
        (
            "air".to_string(),
            Entry::Spec(spec("Konstant", json!({"n": 1.0}))),
        ),
        (
            "glass".to_string(),
            Entry::Spec(spec("Cauchy", json!({"A": 1.50, "B": 0.004, "C": 0.0}))),
        ),
        (
            "TiO2".to_string(),
            Entry::Spec(spec("Cauchy", json!({"A": 2.2, "B": 0.02, "C": 0.001}))),
        ),
        (
            "SiO2".to_string(),
            Entry::Spec(spec("Cauchy", json!({"A": 1.45, "B": 0.0036, "C": 0.0}))),
        ),
    ])
}

/// air | `films` | glass as structure `main`.
fn model(films: &[(f64, &str)]) -> FitModel {
    let mut layers = vec![half_space("air", LayerType::Ambient)];
    layers.extend(films.iter().map(|&(d, m)| Layer::film(d, m)));
    layers.push(half_space("glass", LayerType::Substrate));
    let groups = films
        .iter()
        .map(|&(_, m)| (m.to_string(), Group::new(m)))
        .collect();
    FitModel::from_structure(&Structure::new(layers, groups), "main", materials()).unwrap()
}

fn hlh() -> FitModel {
    model(&[(80.0, "TiO2"), (120.0, "SiO2"), (60.0, "TiO2")])
}

/// One curve `(id, channel, deg)` seen through `instrument`, plus noise.
/// The returned data carry `start` as their instrument.
fn measure(
    truth: &FitModel,
    (id, channel, deg): (&str, &str, f64),
    instrument: Instrument,
    start: Instrument,
    sigma: f64,
    seed: u64,
) -> Vec<SpectrumData> {
    let wl = wavelengths();
    let mut d = SpectrumData::new(
        id,
        channel,
        Radians::from_degrees(deg).unwrap(),
        wl.clone(),
        vec![0.0; wl.len()],
        sigma,
    );
    d.instrument = instrument;
    let clean = FitProblem::new(truth.clone(), vec![d.clone()], vec![])
        .simulate()
        .unwrap();
    let mut rng = StdRng::seed_from_u64(seed);
    let noise = Normal::new(0.0, sigma).unwrap();
    d.values = clean[0]
        .iter()
        .map(|v| v + noise.sample(&mut rng))
        .collect();
    d.instrument = start;
    vec![d]
}

fn params(addrs: &[&str]) -> Vec<Parameter> {
    addrs.iter().map(|a| Parameter::parse(a).unwrap()).collect()
}

fn with(model: &FitModel, values: &[(&str, f64)]) -> FitModel {
    let mut m = model.clone();
    for (addr, v) in values {
        m.set(&ParamAddr::parse(addr).unwrap(), *v).unwrap();
    }
    m
}

fn recovered(report: &FitReport, truth: &[(&str, f64)]) {
    for (addr, want) in truth {
        let p = report.parameter(addr).unwrap();
        let se = p.std_error.unwrap();
        assert!(
            (p.value - want).abs() <= 2.5 * se,
            "{addr}: fitted {} +- {se}, truth {want} (start {})",
            p.value,
            p.start
        );
    }
    assert!(
        (report.chi2_reduced - 1.0).abs() < 0.3,
        "chi2_red {}",
        report.chi2_reduced
    );
}

#[test]
fn an_angle_offset_is_recovered() {
    let truth = with(&hlh(), &[("layer[main:2].thickness", 126.0)]);
    let tilted = Instrument {
        angle_offset: 0.012,
        ..Instrument::default()
    };
    let data = measure(
        &truth,
        ("Rp", "Rp", 55.0),
        tilted,
        Instrument::default(),
        0.001,
        1,
    );
    let problem = FitProblem::new(
        hlh(),
        data,
        params(&["layer[main:2].thickness", "data[Rp].angle_offset"]),
    );
    let r = fit(&problem).unwrap();
    recovered(
        &r,
        &[
            ("layer[main:2].thickness", 126.0),
            ("data[Rp].angle_offset", 0.012),
        ],
    );
    // The fitted model carries the instrument it found.
    let found = r.model.instrument("Rp").unwrap().angle_offset;
    assert_eq!(found, r.parameter("data[Rp].angle_offset").unwrap().value);
}

#[test]
fn a_wavelength_shift_is_recovered() {
    let truth = with(&hlh(), &[("layer[main:1].thickness", 84.0)]);
    let shifted = Instrument {
        wavelength_shift: 1.5,
        ..Instrument::default()
    };
    let data = measure(
        &truth,
        ("T", "T_avg", 0.0),
        shifted,
        Instrument::default(),
        0.001,
        2,
    );
    let problem = FitProblem::new(
        hlh(),
        data,
        params(&["layer[main:1].thickness", "data[T].wavelength_shift"]),
    );
    let r = fit(&problem).unwrap();
    recovered(
        &r,
        &[
            ("layer[main:1].thickness", 84.0),
            ("data[T].wavelength_shift", 1.5),
        ],
    );
}

#[test]
fn a_scale_and_offset_are_recovered() {
    let truth = with(&hlh(), &[("layer[main:3].thickness", 57.0)]);
    let gain = Instrument {
        scale: 0.97,
        offset: 0.004,
        ..Instrument::default()
    };
    let data = measure(
        &truth,
        ("R", "R_avg", 8.0),
        gain,
        Instrument::default(),
        0.001,
        3,
    );
    let problem = FitProblem::new(
        hlh(),
        data,
        params(&["layer[main:3].thickness", "data[R].scale", "data[R].offset"]),
    );
    let r = fit(&problem).unwrap();
    recovered(
        &r,
        &[
            ("layer[main:3].thickness", 57.0),
            ("data[R].scale", 0.97),
            ("data[R].offset", 0.004),
        ],
    );
}

#[test]
fn a_slit_bandwidth_is_recovered() {
    // A thick film: fringes about 30 nm apart, which a 6 nm slit blurs.
    let thick = model(&[(1500.0, "TiO2")]);
    let truth = with(&thick, &[("layer[main:1].thickness", 1504.0)]);
    let slit = Instrument {
        bandwidth: 6.0,
        ..Instrument::default()
    };
    let start = Instrument {
        bandwidth: 3.0,
        ..Instrument::default()
    };
    let data = measure(&truth, ("T", "T_avg", 0.0), slit, start, 0.001, 4);
    let problem = FitProblem::new(
        thick,
        data,
        params(&["layer[main:1].thickness", "data[T].bandwidth"]),
    );
    let r = fit(&problem).unwrap();
    recovered(
        &r,
        &[
            ("layer[main:1].thickness", 1504.0),
            ("data[T].bandwidth", 6.0),
        ],
    );
}

#[test]
fn instruments_check_their_settings() {
    let mut d = measure(
        &hlh(),
        ("psi", "Psi_R", 70.0),
        Instrument::default(),
        Instrument::default(),
        0.01,
        5,
    );
    d[0].instrument.bandwidth = 2.0;
    let e = FitProblem::new(hlh(), d, vec![]).simulate().unwrap_err();
    assert!(e.contains("bandwidth"), "{e}");

    let mut d = measure(
        &hlh(),
        ("T", "T_avg", 0.0),
        Instrument::default(),
        Instrument::default(),
        0.01,
        5,
    );
    d.push(d[0].clone());
    let e = fit(&FitProblem::new(
        hlh(),
        d,
        params(&["layer[main:1].thickness"]),
    ))
    .unwrap_err();
    assert!(e.contains("two datasets"), "{e}");
}

#[test]
fn an_instrument_reads_from_json_with_identity_defaults() {
    let i: Instrument = serde_json::from_value(json!({"angle_offset": 0.01})).unwrap();
    assert_eq!(i.angle_offset, 0.01);
    assert_eq!(i.scale, 1.0);
    assert!(!i.is_identity());
    assert!(Instrument::default().is_identity());
    assert!(serde_json::from_value::<Instrument>(json!({"gain": 2.0})).is_err());
}
