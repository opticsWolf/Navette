// SPDX-License-Identifier: LGPL-3.0-or-later
//! Fit engine, phase 3 (`docs/plans/general_fit_plan.md` §3.8, §3.5.5):
//! staged strategies, several samples, ties, priors, instrument data.
//!
//! The four example strategies of §3.8 run from JSON and recover a seeded
//! truth; a coating and its bare substrate share the substrate's
//! parameters; every converter preset round-trips, and σ propagated
//! through tan Ψ / cos Δ and N, C, S matches a finite difference.

use std::collections::{BTreeMap, HashMap};
use std::f64::consts::PI;

use rand::SeedableRng;
use rand::rngs::StdRng;
use rand_distr::{Distribution, Normal};
use serde_json::json;

use navette::fit::{
    FitModel, FitProblem, FitReport, InstrumentConvention, NkColumn, NkData, NkQuantity, PRIOR_ID,
    ParamAddr, ParamSet, Parameter, RawDataset, RawQuantity, Segment, Side, SigmaSource,
    SpectrumData, Strategy, Tie, fit, from_solver, run_strategy, to_solver,
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

/// air | `films` | glass as structure `main`, a group per film material.
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

fn with(model: &FitModel, values: &[(&str, f64)]) -> FitModel {
    let mut m = model.clone();
    for (addr, v) in values {
        m.set(&ParamAddr::parse(addr).unwrap(), *v).unwrap();
    }
    m
}

/// Truth curves for `(id, env, channel, deg)`, plus seeded noise.
fn measure(
    truth: &FitModel,
    curves: &[(&str, &str, &str, f64)],
    sigma: f64,
    seed: u64,
) -> Vec<SpectrumData> {
    let template: Vec<SpectrumData> = curves
        .iter()
        .map(|&(id, env, ch, deg)| {
            SpectrumData::new(
                id,
                ch,
                Radians::from_degrees(deg).unwrap(),
                wavelengths(),
                vec![0.0; wavelengths().len()],
                sigma,
            )
            .on(env)
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

fn params(addrs: &[&str]) -> Vec<Parameter> {
    addrs.iter().map(|a| Parameter::parse(a).unwrap()).collect()
}

fn within_two_sigma(report: &FitReport, truth: &[(&str, f64)]) {
    for (addr, want) in truth {
        let p = report.parameter(addr).unwrap();
        let se = p.std_error.unwrap();
        assert!(
            (p.value - want).abs() <= 2.0 * se,
            "{addr}: fitted {} ± {se}, truth {want} (start {})",
            p.value,
            p.start
        );
    }
}

fn chi2_red_near_one(report: &FitReport) {
    assert!(
        (report.chi2_reduced - 1.0).abs() < 0.25,
        "chi2_red {}",
        report.chi2_reduced
    );
}

// --- the four §3.8 strategies ---------------------------------------------

#[test]
fn materials_first_from_a_witness_then_the_layer() {
    let truths = [
        ("material[TiO2].A", 2.24),
        ("material[TiO2].B", 0.024),
        ("layer[main:1].thickness", 212.0),
    ];
    let truth = with(&model(&[(200.0, "TiO2")]), &truths);
    let data = measure(&truth, &[("T", "default", "T_avg", 0.0)], 0.002, 5);
    // The witness: nk of the same material, measured on its own.
    let Entry::Spec(tio2) = &truth.materials()["TiO2"] else {
        unreachable!()
    };
    let wl = wavelengths();
    let mut rng = StdRng::seed_from_u64(6);
    let noise = Normal::new(0.0, 1e-3).unwrap();
    let n: Vec<f64> = tio2
        .evaluate(&wl)
        .unwrap()
        .iter()
        .map(|z| NkQuantity::N.of(*z) + noise.sample(&mut rng))
        .collect();
    let witness = NkData::new(
        "witness",
        "TiO2",
        wl,
        vec![NkColumn::new(NkQuantity::N, n, 1e-3)],
    );
    let mut problem = FitProblem::new(
        model(&[(200.0, "TiO2")]),
        data,
        params(&[
            "material[TiO2].A",
            "material[TiO2].B",
            "layer[main:1].thickness",
        ]),
    );
    problem.nk_data = vec![witness];

    let strategy = Strategy::from_json(
        r#"{"name": "materials first", "stages": [
            {"name": "witness", "free": ["set:material[TiO2]"], "datasets": ["witness"]},
            {"name": "thickness", "free": ["kind=thickness"], "datasets": ["T"]},
            {"name": "polish", "free": ["*"]}
        ]}"#,
    )
    .unwrap();
    let r = run_strategy(&problem, &strategy, &[]).unwrap();
    assert_eq!(r.stages.len(), 3);
    assert_eq!(r.stages[0].free, ["material[TiO2].A", "material[TiO2].B"]);
    // Stage 1 never saw the spectrum: its report has the witness only.
    assert_eq!(r.stages[0].report.datasets.len(), 1);
    assert_eq!(r.report.datasets.len(), 2);
    assert_eq!(r.report.n_free, 3);
    within_two_sigma(&r.report, &truths);
    chi2_red_near_one(&r.report);
}

#[test]
fn thickness_first_in_a_window_then_indices_then_dispersion() {
    let truths = [
        ("layer[main:1].thickness", 82.0),
        ("layer[main:2].thickness", 117.0),
        ("layer[main:3].thickness", 62.0),
        ("group[TiO2].n_factor", 1.01),
        ("material[SiO2].A", 1.46),
    ];
    let truth = with(&hlh(), &truths);
    let data = measure(
        &truth,
        &[("R", "default", "R_avg", 8.0), ("T", "default", "Ts", 45.0)],
        0.002,
        21,
    );
    let problem = FitProblem::new(
        hlh(),
        data,
        params(&[
            "layer[main:1].thickness",
            "layer[main:2].thickness",
            "layer[main:3].thickness",
            "group[TiO2].n_factor",
            "material[SiO2].A",
        ]),
    );
    let strategy = Strategy::from_json(
        r#"{"name": "thickness first", "stages": [
            {"name": "thicknesses", "free": ["kind=thickness"], "window": [600, 1000]},
            {"name": "indices", "free": ["group[*].n_factor"]},
            {"name": "dispersion", "free": ["set:materials"]},
            {"name": "polish", "free": ["*"], "accept": "if_better",
             "repeat": {"until_no_gain": {"rel": 1e-6, "max": 3}}}
        ]}"#,
    )
    .unwrap();
    let r = run_strategy(&problem, &strategy, &[]).unwrap();
    // The window scored fewer points than the full range.
    assert!(r.stages[0].report.n_residuals < r.report.n_residuals);
    for s in &r.stages[..3] {
        assert!(s.report.chi2 <= s.chi2_start, "{}", s.name);
    }
    assert!(r.stages[3].runs >= 1);
    within_two_sigma(&r.report, &truths);
    chi2_red_near_one(&r.report);
}

#[test]
fn a_group_factor_model_is_compared_on_chi2() {
    // The truth: every TiO2 layer 4 % thick, SiO2 1 % low in index.
    let truth = with(
        &hlh(),
        &[
            ("group[TiO2].thick_factor", 1.04),
            ("group[SiO2].n_factor", 0.99),
        ],
    );
    let data = measure(
        &truth,
        &[("R", "default", "R_avg", 8.0), ("T", "default", "Tp", 50.0)],
        0.002,
        3,
    );
    let run = |addrs: &[&str], json: &str| {
        let problem = FitProblem::new(hlh(), data.clone(), params(addrs));
        run_strategy(&problem, &Strategy::from_json(json).unwrap(), &[]).unwrap()
    };
    let groups = run(
        &["group[TiO2].thick_factor", "group[SiO2].n_factor"],
        r#"{"name": "groups", "stages": [
            {"name": "factors", "free": ["set:groups"]}
        ]}"#,
    );
    let indices = run(
        &["group[TiO2].n_factor", "group[SiO2].n_factor"],
        r#"{"name": "indices", "stages": [
            {"name": "factors", "free": ["set:groups"]}
        ]}"#,
    );
    chi2_red_near_one(&groups.report);
    within_two_sigma(
        &groups.report,
        &[
            ("group[TiO2].thick_factor", 1.04),
            ("group[SiO2].n_factor", 0.99),
        ],
    );
    assert!(
        indices.report.chi2 > 2.0 * groups.report.chi2,
        "{} vs {}",
        indices.report.chi2,
        groups.report.chi2
    );
    assert!(groups.report.bic < indices.report.bic);
}

#[test]
fn a_scan_finds_the_basin_a_plain_fit_misses() {
    let addr = "layer[main:1].thickness";
    let truth = with(&model(&[(500.0, "TiO2")]), &[(addr, 640.0)]);
    let data = measure(&truth, &[("R", "default", "R_avg", 8.0)], 0.002, 9);
    let problem = FitProblem::new(
        model(&[(500.0, "TiO2")]),
        data,
        vec![Parameter::parse(addr).unwrap().bounds(300.0, 900.0)],
    );
    // A thick film's fringes: from 500 nm, a plain fit stops in a wrong
    // fringe order.
    let plain = fit(&problem).unwrap();
    assert!(plain.chi2_reduced > 5.0, "plain fit {}", plain.chi2_reduced);

    let strategy = Strategy::from_json(
        r#"{"name": "scan", "stages": [
            {"name": "scan", "free": ["kind=thickness"],
             "kind": {"scan": {"param": "layer[main:1].thickness", "points": 241}}}
        ]}"#,
    )
    .unwrap();
    let r = run_strategy(&problem, &strategy, &[]).unwrap();
    within_two_sigma(&r.report, &[(addr, 640.0)]);
    chi2_red_near_one(&r.report);

    // Multi-start reaches the same basin; the same seed gives the same bits.
    let multi = Strategy::from_json(
        r#"{"name": "multi", "stages": [
            {"name": "starts", "free": ["*"],
             "kind": {"multi_start": {"starts": 24, "sampler": "lhs", "seed": 4}}}
        ]}"#,
    )
    .unwrap();
    let a = run_strategy(&problem, &multi, &[]).unwrap();
    let b = run_strategy(&problem, &multi, &[]).unwrap();
    assert_eq!(a.stages[0].runs, 25);
    within_two_sigma(&a.report, &[(addr, 640.0)]);
    assert_eq!(
        a.report.parameters[0].value.to_bits(),
        b.report.parameters[0].value.to_bits()
    );
}

#[test]
fn user_sets_and_rollback() {
    let truth = with(&hlh(), &[("layer[main:2].thickness", 124.0)]);
    let data = measure(&truth, &[("R", "default", "R_avg", 8.0)], 0.002, 2);
    let problem = FitProblem::new(
        hlh(),
        data,
        vec![
            Parameter::parse("layer[main:2].thickness")
                .unwrap()
                .bounds(20.0, 400.0),
            Parameter::parse("group[TiO2].n_factor").unwrap(),
        ],
    );
    let sets = [ParamSet::new("L", &["layer[main:2].*"])];
    // Stage 2 scans two far points and lands in a worse fringe order:
    // `if_better` rolls it back to stage 1's end.
    let strategy = Strategy::from_json(
        r#"{"name": "s", "stages": [
            {"name": "L", "free": ["set:L"]},
            {"name": "far", "free": ["set:L"], "accept": "if_better",
             "kind": {"scan": {"param": "set:L", "points": 2}}}
        ]}"#,
    )
    .unwrap();
    let r = run_strategy(&problem, &strategy, &sets).unwrap();
    assert_eq!(r.stages[0].free, ["layer[main:2].thickness"]);
    assert!(!r.stages[1].accepted);
    let end = r.stages[0]
        .report
        .parameter("layer[main:2].thickness")
        .unwrap()
        .value;
    assert_eq!(
        r.report.parameter("layer[main:2].thickness").unwrap().value,
        end
    );
    // The never-freed group factor kept its start and has no error bar.
    let g = r.report.parameter("group[TiO2].n_factor").unwrap();
    assert_eq!(g.value, 1.0);
    assert!(g.std_error.is_none());

    let bad = Strategy::from_json(r#"{"name": "s", "stages": [{"name": "x", "free": ["set:H"]}]}"#)
        .unwrap();
    let e = run_strategy(&problem, &bad, &sets).unwrap_err();
    assert!(
        e.contains("stage \"x\"") && e.contains("no set \"H\""),
        "{e}"
    );
}

// --- several samples ------------------------------------------------------

#[test]
fn a_coating_and_its_bare_substrate_share_the_substrate() {
    let mut start = model(&[(150.0, "TiO2")]);
    start
        .add_environment(
            "bare",
            vec![
                Segment::Fixed(half_space("air", LayerType::Ambient)),
                Segment::Fixed(half_space("glass", LayerType::Substrate)),
            ],
        )
        .unwrap();
    assert_eq!(start.environment_names(), ["default", "bare"]);
    assert_eq!(start.entries("bare").unwrap().len(), 2);
    assert!(
        start
            .entries("cover")
            .unwrap_err()
            .contains("no environment")
    );

    let truths = [
        ("layer[main:1].thickness", 156.0),
        ("material[TiO2].A", 2.23),
        ("material[glass].A", 1.515),
    ];
    let truth = with(&start, &truths);
    let data = measure(
        &truth,
        &[
            ("coat", "default", "R_avg", 8.0),
            ("bare_s", "bare", "Rs", 45.0),
            ("bare_p", "bare", "Rp", 60.0),
        ],
        0.001,
        17,
    );
    let addrs = [
        "layer[main:1].thickness",
        "material[TiO2].A",
        "material[glass].A",
    ];
    let joint = fit(&FitProblem::new(
        start.clone(),
        data.clone(),
        params(&addrs),
    ))
    .unwrap();
    within_two_sigma(&joint, &truths);
    chi2_red_near_one(&joint);
    // One glass, both samples: the bare sample pins it, so the joint
    // error is smaller than the coating alone gives.
    let coat_only = fit(&FitProblem::new(
        start.clone(),
        data[..1].to_vec(),
        params(&addrs),
    ))
    .unwrap();
    let se = |r: &FitReport| r.parameter("material[glass].A").unwrap().std_error.unwrap();
    assert!(
        se(&joint) < 0.5 * se(&coat_only),
        "{} vs {}",
        se(&joint),
        se(&coat_only)
    );

    // A dataset on an unknown sample is refused by name.
    let mut stray = data.clone();
    stray[1].env = "cover".into();
    let e = fit(&FitProblem::new(start, stray, params(&addrs))).unwrap_err();
    assert!(e.contains("cover"), "{e}");
}

// --- ties and priors ------------------------------------------------------

#[test]
fn a_tie_moves_its_followers_with_the_driver() {
    let truth = with(
        &hlh(),
        &[
            ("layer[main:1].thickness", 84.0),
            ("layer[main:3].thickness", 84.0),
        ],
    );
    let data = measure(&truth, &[("R", "default", "R_avg", 8.0)], 0.002, 8);
    let mut problem = FitProblem::new(hlh(), data, params(&["layer[main:1].thickness"]));
    problem.ties = vec![
        Tie::parse("layer[main:1].thickness")
            .unwrap()
            .equal("layer[main:3].thickness")
            .unwrap(),
    ];
    let r = fit(&problem).unwrap();
    within_two_sigma(&r, &[("layer[main:1].thickness", 84.0)]);
    let get = |a: &str| r.model.get(&ParamAddr::parse(a).unwrap()).unwrap();
    assert_eq!(
        get("layer[main:3].thickness"),
        get("layer[main:1].thickness")
    );

    // A follower that is also a parameter is refused.
    problem.parameters = params(&["layer[main:1].thickness", "layer[main:3].thickness"]);
    let e = fit(&problem).unwrap_err();
    assert!(e.contains("layer[main:3].thickness"), "{e}");

    // A scaled follower: SiO2 at 1.5 × the first TiO2.
    problem.parameters = params(&["layer[main:1].thickness"]);
    problem.ties = vec![
        Tie::parse("layer[main:1].thickness")
            .unwrap()
            .follow("layer[main:2].thickness", 1.5, 0.0)
            .unwrap(),
    ];
    let r = fit(&problem).unwrap();
    let get = |a: &str| r.model.get(&ParamAddr::parse(a).unwrap()).unwrap();
    assert!((get("layer[main:2].thickness") - 1.5 * get("layer[main:1].thickness")).abs() < 1e-9);
}

#[test]
fn a_prior_pulls_as_hard_as_its_sigma() {
    let data = measure(&hlh(), &[("R", "default", "R_avg", 8.0)], 0.002, 4);
    let addr = "layer[main:2].thickness";
    let fit_with = |sigma: f64| {
        let p = Parameter::parse(addr).unwrap().prior(130.0, sigma);
        fit(&FitProblem::new(hlh(), data.clone(), vec![p])).unwrap()
    };
    let loose = fit_with(1e3);
    let tight = fit_with(1e-3);
    let at = |r: &FitReport| r.parameter(addr).unwrap().value;
    assert!((at(&loose) - 120.0).abs() < 0.5, "{}", at(&loose));
    assert!((at(&tight) - 130.0).abs() < 0.01, "{}", at(&tight));
    let prior = tight.datasets.iter().find(|d| d.id == PRIOR_ID).unwrap();
    assert_eq!(prior.n_points, 1);
    assert_eq!(tight.n_residuals, wavelengths().len() + 1);
}

// --- the instrument converter ---------------------------------------------

fn solver_psi_delta(delta_lo: f64, delta_hi: f64) -> [SpectrumData; 2] {
    let wl = wavelengths();
    let n = wl.len();
    let t = |i: usize| i as f64 / (n - 1) as f64;
    let angle = Radians::from_degrees(70.0).unwrap();
    let psi = (0..n).map(|i| 0.1 + 1.3 * t(i)).collect();
    let delta = (0..n)
        .map(|i| delta_lo + (delta_hi - delta_lo) * t(i))
        .collect();
    [
        SpectrumData::new("e.psi", "Psi_R", angle, wl.clone(), psi, 1.0),
        SpectrumData::new("e.delta", "Delta_R", angle, wl, delta, 1.0),
    ]
}

fn close(a: &[f64], b: &[f64], tol: f64, what: &str) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        assert!((x - y).abs() <= tol * (1.0 + y.abs()), "{what}: {x} vs {y}");
    }
}

#[test]
fn every_preset_round_trips() {
    let names = InstrumentConvention::preset_names();
    assert!(names.len() >= 10);
    let rs = SpectrumData::new(
        "r",
        "Rs",
        Radians::from_degrees(45.0).unwrap(),
        wavelengths(),
        (0..121).map(|i| 0.05 + 0.004 * i as f64).collect(),
        1.0,
    );
    for name in &names {
        let conv = InstrumentConvention::preset(name).unwrap();
        // tan Ψ / cos Δ keeps Δ in [0, π] only.
        let tan_cos = name.contains("tan_cos");
        let (lo, hi) = if tan_cos { (0.05, 3.1) } else { (-3.1, 3.1) };
        let ell = solver_psi_delta(lo, hi);
        let raw = from_solver("e", &[&ell[0], &ell[1]], &conv).unwrap();
        let back = to_solver(&raw, &conv).unwrap();
        assert_eq!(back.len(), 2, "{name}");
        for (b, want) in back.iter().zip(&ell) {
            assert_eq!(b.id, want.id, "{name}");
            assert_eq!(b.channel, want.channel, "{name}");
            close(&b.wavelengths, &want.wavelengths, 1e-12, name);
            close(&b.values, &want.values, 1e-9, name);
            assert!((b.angle.get() - want.angle.get()).abs() < 1e-12, "{name}");
            assert_eq!(b.source.as_deref(), Some(conv.describe().as_str()));
        }
        let raw = from_solver("r", &[&rs], &conv).unwrap();
        let back = to_solver(&raw, &conv).unwrap();
        close(&back[0].values, &rs.values, 1e-12, name);
        close(&back[0].wavelengths, &rs.wavelengths, 1e-12, name);
    }
    assert!(
        InstrumentConvention::preset("vendor")
            .unwrap_err()
            .contains("presets")
    );
}

#[test]
fn conventions_read_the_same_measurement_alike() {
    // Δ = -40° as each convention writes it.
    let deg = |d: f64| d.to_radians();
    let read = |name: &str, columns: Vec<f64>, x: f64| {
        let conv = InstrumentConvention::preset(name).unwrap();
        let q = RawQuantity::Ellipsometry(Side::R);
        let raw = RawDataset::new(
            "m",
            q,
            70.0,
            vec![x],
            columns.into_iter().map(|c| vec![c]).collect(),
        );
        let d = to_solver(&raw, &conv).unwrap();
        (d[0].values[0], d[1].values[0], d[0].wavelengths[0])
    };
    let want = (deg(30.0), deg(-40.0), 620.0);
    let ev = 1239.841984 / 620.0;
    let s2 = deg(60.0).sin();
    let cases = [
        read("nm_deg", vec![30.0, -40.0], 620.0),
        read("nm_deg_delta360", vec![30.0, 320.0], 620.0),
        read("nm_deg_opposite_delta", vec![30.0, 40.0], 620.0),
        read("um_deg", vec![30.0, -40.0], 0.62),
        read("ev_deg", vec![30.0, -40.0], ev),
        read(
            "ev_deg_ncs",
            vec![
                deg(60.0).cos(),
                s2 * deg(-40.0).cos(),
                s2 * deg(-40.0).sin(),
            ],
            ev,
        ),
        read(
            "ev_deg_ncs_opposite_delta",
            vec![deg(60.0).cos(), s2 * deg(40.0).cos(), s2 * deg(40.0).sin()],
            ev,
        ),
    ];
    for (i, got) in cases.iter().enumerate() {
        assert!((got.0 - want.0).abs() < 1e-12, "case {i}: psi {}", got.0);
        assert!((got.1 - want.1).abs() < 1e-12, "case {i}: delta {}", got.1);
        assert!((got.2 - want.2).abs() < 1e-9, "case {i}: wl {}", got.2);
    }
    // tan Ψ / cos Δ cannot tell -40° from +40°.
    let t = read(
        "nm_deg_tan_cos",
        vec![deg(30.0).tan(), deg(-40.0).cos()],
        620.0,
    );
    assert!((t.1 - deg(40.0)).abs() < 1e-12);

    // Percent and wavenumbers on a spectrum.
    let conv = InstrumentConvention::preset("wavenumber_deg_percent").unwrap();
    let raw = RawDataset::new(
        "t",
        RawQuantity::Spectrum("T_avg".into()),
        0.0,
        vec![1e7 / 500.0],
        vec![vec![91.5]],
    );
    let d = to_solver(&raw, &conv).unwrap();
    assert!((d[0].values[0] - 0.915).abs() < 1e-15);
    assert!((d[0].wavelengths[0] - 500.0).abs() < 1e-9);
}

/// The converter's σ against a central difference of its own map.
#[test]
fn sigma_propagates_through_the_map() {
    let deg = |d: f64| d.to_radians();
    for (name, columns) in [
        ("nm_deg_tan_cos", vec![deg(25.0).tan(), deg(70.0).cos()]),
        (
            "ev_deg_ncs",
            vec![
                deg(50.0).cos(),
                deg(50.0).sin() * deg(110.0).cos(),
                deg(50.0).sin() * deg(110.0).sin(),
            ],
        ),
    ] {
        let sig = [1e-3, 2e-3, 3e-3];
        let conv = InstrumentConvention::preset(name)
            .unwrap()
            .with_sigma(SigmaSource::Column);
        let raw = |cols: &[f64]| {
            let mut r = RawDataset::new(
                "m",
                RawQuantity::Ellipsometry(Side::R),
                70.0,
                vec![2.0],
                cols.iter().map(|c| vec![*c]).collect(),
            );
            r.sigma = sig[..cols.len()].iter().map(|s| vec![*s]).collect();
            r
        };
        let at = |cols: &[f64]| {
            let d = to_solver(&raw(cols), &conv).unwrap();
            [d[0].values[0], d[1].values[0]]
        };
        let got = to_solver(&raw(&columns), &conv).unwrap();
        for (out, g) in got.iter().enumerate() {
            let mut var = 0.0;
            for c in 0..columns.len() {
                let h = 1e-6;
                let (mut up, mut dn) = (columns.clone(), columns.clone());
                up[c] += h;
                dn[c] -= h;
                let d = (at(&up)[out] - at(&dn)[out]) / (2.0 * h);
                var += (d * sig[c]).powi(2);
            }
            let want = var.sqrt();
            let have = g.sigma[0];
            assert!(
                (have - want).abs() < 1e-6 * want,
                "{name} output {out}: sigma {have} vs {want}"
            );
        }
    }
}

#[test]
fn converted_data_fits_and_says_how_it_was_read() {
    // SiO2 on glass seen by an instrument writing Δ in [0°, 360°).
    let truth = with(
        &model(&[(100.0, "SiO2")]),
        &[("layer[main:1].thickness", 104.0)],
    );
    let solver = measure(
        &truth,
        &[
            ("e.psi", "default", "Psi_R", 60.0),
            ("e.delta", "default", "Delta_R", 60.0),
        ],
        1e-9,
        1,
    );
    let conv = InstrumentConvention::preset("nm_deg_delta360")
        .unwrap()
        .with_sigma(SigmaSource::Constant(0.05));
    let raw = from_solver("e", &[&solver[0], &solver[1]], &conv).unwrap();
    assert!(raw.columns[1].iter().all(|d| (0.0..360.0).contains(d)));
    let data = to_solver(&raw, &conv).unwrap();
    assert!((data[0].sigma[0] - deg_rad(0.05)).abs() < 1e-15);
    let r = fit(&FitProblem::new(
        model(&[(100.0, "SiO2")]),
        data,
        params(&["layer[main:1].thickness"]),
    ))
    .unwrap();
    assert!((r.parameter("layer[main:1].thickness").unwrap().value - 104.0).abs() < 1e-6);
    for d in &r.datasets {
        let s = d.source.as_deref().unwrap();
        assert!(s.contains("\"delta_range\":\"positive\""), "{s}");
    }
}

fn deg_rad(d: f64) -> f64 {
    d * PI / 180.0
}

// --- global stages ----------------------------------------------------------

#[cfg(not(feature = "opt-basin"))]
#[test]
fn a_basin_stage_without_basin_says_so_and_de_runs() {
    let data = measure(&hlh(), &[("R", "default", "R_avg", 8.0)], 0.002, 3);
    let p = FitProblem::new(hlh(), data, params(&["layer[main:1].thickness"]));
    let stage = |method: &str| {
        Strategy::from_json(&format!(
            r#"{{"name": "g", "stages": [{{"name": "g", "free": ["*"], "kind": {{"global": {{"method": "{method}", "budget": 50}}}}}}]}}"#
        ))
        .unwrap()
    };
    let e = run_strategy(&p, &stage("basin_de"), &[]).unwrap_err();
    assert!(e.contains("opt-basin"), "{e}");
    run_strategy(&p, &stage("de"), &[]).unwrap();
}

#[test]
fn analytic_and_differenced_jacobians_fit_alike() {
    let truth = with(
        &hlh(),
        &[
            ("layer[main:2].thickness", 124.0),
            ("material[TiO2].A", 2.23),
        ],
    );
    let data = measure(
        &truth,
        &[("R", "default", "R_avg", 8.0), ("T", "default", "Tp", 45.0)],
        0.002,
        9,
    );
    let mut p = FitProblem::new(
        hlh(),
        data,
        params(&[
            "layer[main:1].thickness",
            "layer[main:2].thickness",
            "group[TiO2].thick_factor",
            "material[TiO2].A",
        ]),
    );
    let a = fit(&p).unwrap();
    p.lsq.jacobian = navette::lsq::lm::JacobianMode::Fd;
    let f = fit(&p).unwrap();
    assert!(
        (a.chi2 / f.chi2 - 1.0).abs() < 1e-6,
        "{} vs {}",
        a.chi2,
        f.chi2
    );
    for (x, y) in a.parameters.iter().zip(&f.parameters) {
        let (sx, sy) = (x.std_error.unwrap(), y.std_error.unwrap());
        assert!(
            (x.value - y.value).abs() < 1e-3 * sx,
            "{}: {} vs {}",
            x.addr,
            x.value,
            y.value
        );
        assert!((sx / sy - 1.0).abs() < 1e-4, "{}: se {sx} vs {sy}", x.addr);
    }
}
