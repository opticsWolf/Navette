// SPDX-License-Identifier: LGPL-3.0-or-later
//! Where a fit spends its time (`docs/plans/general_fit_plan.md` phase 6).
//!
//! Each case is a fit request (a program document with a `fit` section)
//! whose data was simulated from a perturbed truth with noise; the fit
//! starts from the nominal model. Each case runs twice, `fd` (central
//! differences over the whole residual vector) and `an` (the fit's own
//! Jacobian, the default). Per run:
//!
//! - `eval`: one forward evaluation (`FitProblem::simulate`), best of 5;
//! - `jac`: one Jacobian, measured as `evaluate` (two evaluations, the
//!   Jacobian of the covariance, the report) minus two `eval`s;
//! - `jac/eval`: how many evaluations a Jacobian costs (`fd`: `2·n`);
//! - the fit itself: wall time (best of 3), iterations, evaluations
//!   (`fd` counts its Jacobians' evaluations, `an` does not), χ²_red;
//! - `speedup`: the `fd` fit's time over the `an` one's.
//!
//! ```text
//! cargo run --release -p navette --example fit_perf [--json PATH]
//! ```

use std::process::ExitCode;
use std::time::{Duration, Instant};

use serde_json::{Value, json};

use navette::fit::{FitProblem, evaluate, fit, request_from_json};
use navette::lsq::lm::JacobianMode;

const LAYERS: usize = 20;

/// A deterministic standard normal (xorshift + Box-Muller): no RNG crate
/// in an example.
struct Noise(u64);

impl Noise {
    fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn normal(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

fn wl(n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| 400.0 + 600.0 * i as f64 / (n - 1) as f64)
        .collect()
}

/// The program's materials; `k` scales the truth's Cauchy A.
fn materials(k: f64) -> Value {
    json!([
        {"name": "air", "model": "Konstant", "params": {"n": 1.0}},
        {"name": "glass", "model": "Konstant", "params": {"n": 1.52 * (1.0 + 0.3 * (k - 1.0))}},
        {"name": "TiO2", "model": "Cauchy", "params": {"A": 2.2 * k, "B": 0.02, "C": 0.001}},
        {"name": "SiO2", "model": "Cauchy", "params": {"A": 1.45 * k, "B": 0.0036, "C": 0.0}},
    ])
}

/// A chirped (H L)^10 edge; the truth's thicknesses are off by ±`err`.
fn coat(err: f64) -> Value {
    let layers: Vec<Value> = (0..LAYERS)
        .map(|i| {
            let (m, n) = if i % 2 == 0 {
                ("TiO2", 2.3)
            } else {
                ("SiO2", 1.46)
            };
            let stretch = 0.8 + 0.4 * i as f64 / (LAYERS - 1) as f64;
            let d = (700.0 / (4.0 * n) * stretch * 10.0_f64).round() / 10.0;
            let sign = if i % 3 == 0 { 1.0 } else { -1.0 };
            json!({"material_code": m, "thickness_nm": d * (1.0 + sign * err)})
        })
        .collect();
    json!({"label": "coat", "layers": layers, "groups": []})
}

fn film(label: &str, material: &str, d: f64) -> Value {
    json!({"label": label, "layers": [{"material_code": material, "thickness_nm": d}],
           "groups": []})
}

fn env(name: &str, stack: &[&str]) -> Value {
    let stack: Vec<Value> = stack.iter().map(|s| json!({"structure": s})).collect();
    json!({"name": name, "ambient": "air", "substrate": "glass", "stack": stack})
}

fn spectrum(id: &str, env: &str, channel: &str, deg: f64, n: usize) -> Value {
    json!({"id": id, "env": env, "channel": channel, "angle_deg": deg,
           "wavelengths": wl(n), "values": vec![0.0; n], "sigma": 0.002})
}

fn nk(id: &str, material: &str, n: usize) -> Value {
    json!({"id": id, "material": material, "wavelengths": wl(n),
           "columns": [{"quantity": "n", "values": vec![0.0; n], "sigma": 0.002}]})
}

/// Several structures need an architect; its blocks name them (each once,
/// labelled as the structure), and the environments stack them.
fn program(materials: Value, structures: Vec<Value>, fit: Value) -> String {
    let blocks: Vec<Value> = structures
        .iter()
        .map(|s| {
            json!({"structure": s["label"], "label": s["label"], "repeat_count": 1,
                       "kind": 0})
        })
        .collect();
    json!({"schema_version": 4, "kind": "program",
           "sections": {"materials": materials, "structures": structures,
                        "architect": {"blocks": blocks}, "fit": fit}})
    .to_string()
}

struct Case {
    name: &'static str,
    /// The truth's program and the nominal one, alike but for the values.
    truth: String,
    nominal: String,
}

fn coating_case(name: &'static str, witnesses: bool, nk_data: bool) -> Case {
    let make = |err: f64, k: f64| {
        let mut structures = vec![coat(err)];
        let mut envs = vec![env("coat", &["coat"])];
        let mut data = vec![
            spectrum("R@8", "coat", "R_avg", 8.0, 151),
            spectrum("T@0", "coat", "T_avg", 0.0, 151),
            spectrum("Rp@45", "coat", "Rp", 45.0, 151),
        ];
        let mut parameters = vec![
            json!({"select": "kind=thickness"}),
            json!({"addr": "material[TiO2].A", "lb": 2.0, "ub": 2.5}),
            json!({"addr": "material[SiO2].A", "lb": 1.3, "ub": 1.6}),
        ];
        if witnesses {
            structures.push(film("witH", "TiO2", 300.0 * (1.0 + err)));
            structures.push(film("witL", "SiO2", 400.0 * (1.0 - err)));
            envs.push(env("witH", &["witH"]));
            envs.push(env("witL", &["witL"]));
            envs.push(env("bare", &[]));
            for (e, n) in [("witH", 101), ("witL", 101)] {
                data.push(spectrum(&format!("R {e}"), e, "R_avg", 8.0, n));
                data.push(spectrum(&format!("T {e}"), e, "T_avg", 0.0, n));
            }
            data.push(spectrum("T bare", "bare", "T_avg", 0.0, 51));
            parameters.push(json!({"addr": "material[glass].n", "lb": 1.4, "ub": 1.7}));
        }
        let mut fit = json!({"environments": envs, "data": data, "parameters": parameters});
        if nk_data {
            fit["nk_data"] = json!([nk("nk TiO2", "TiO2", 121), nk("nk SiO2", "SiO2", 121)]);
            fit["parameters"].as_array_mut().unwrap().extend([
                json!({"addr": "material[TiO2].B", "lb": 0.0, "ub": 0.05}),
                json!({"addr": "material[SiO2].B", "lb": 0.0, "ub": 0.02}),
            ]);
        }
        program(materials(k), structures, fit)
    };
    Case {
        name,
        truth: make(0.01, 1.004),
        nominal: make(0.0, 1.0),
    }
}

fn dispersion_case(name: &'static str, model: &str, nominal: Value, truth: Value) -> Case {
    let make = |params: &Value| {
        let mats = json!([{"name": "m", "model": model, "params": params}]);
        let fit = json!({"nk_data": [nk("nk", "m", 300)],
                         "parameters": [{"select": "material[m].*"}]});
        json!({"schema_version": 4, "kind": "program",
               "sections": {"materials": mats, "fit": fit}})
        .to_string()
    };
    Case {
        name,
        truth: make(&truth),
        nominal: make(&nominal),
    }
}

fn cases() -> Vec<Case> {
    vec![
        coating_case("coat20", false, false),
        coating_case("coat20+witnesses", true, false),
        coating_case("coat20+nk", false, true),
        dispersion_case(
            "cauchy nk",
            "Cauchy",
            json!({"A": 1.9, "B": 0.01, "C": 0.0}),
            json!({"A": 1.95, "B": 0.012, "C": 0.0002}),
        ),
        dispersion_case(
            "sellmeier nk",
            "Sellmeier",
            json!({"B1": 0.7, "C1": 0.005, "B2": 0.4, "C2": 0.0135, "B3": 0.9, "C3": 100.0}),
            json!({"B1": 0.696, "C1": 0.00468, "B2": 0.408, "C2": 0.0135, "B3": 0.897,
                   "C3": 97.9}),
        ),
    ]
}

/// The nominal problem with the truth's simulated data plus noise.
fn problem(case: &Case) -> FitProblem {
    let truth = request_from_json(&case.truth).unwrap().problem;
    let sims = truth.simulate().unwrap();
    let mut p = request_from_json(&case.nominal).unwrap().problem;
    let mut noise = Noise(0x9E37_79B9_7F4A_7C15);
    let (spectra, nk) = sims.split_at(p.data.len());
    for (d, sim) in p.data.iter_mut().zip(spectra) {
        for (v, s) in d.values.iter_mut().zip(sim) {
            *v = s + 0.002 * noise.normal();
        }
    }
    for (d, sim) in p.nk_data.iter_mut().zip(nk) {
        let mut it = sim.iter();
        for c in &mut d.columns {
            for v in &mut c.values {
                *v = it.next().unwrap() + 0.002 * noise.normal();
            }
        }
    }
    p
}

fn best_of(k: usize, mut f: impl FnMut()) -> Duration {
    (0..k)
        .map(|_| {
            let t = Instant::now();
            f();
            t.elapsed()
        })
        .min()
        .unwrap()
}

fn main() -> ExitCode {
    let json_out = std::env::args().skip_while(|a| a != "--json").nth(1);
    println!(
        "{:<18} {:>3} {:>4} {:>6} {:>8} {:>8} {:>8} {:>8} {:>4} {:>6} {:>8} {:>7}",
        "case",
        "jac",
        "n",
        "m",
        "eval ms",
        "jac ms",
        "jac/eval",
        "fit ms",
        "it",
        "evals",
        "chi2_red",
        "speedup"
    );
    let mut rows = Vec::new();
    let ms = |d: Duration| d.as_secs_f64() * 1e3;
    for case in cases() {
        let mut p = problem(&case);
        let t_eval = best_of(5, || {
            p.simulate().unwrap();
        });
        let mut fd_fit = None;
        for (label, mode) in [("fd", JacobianMode::Fd), ("an", JacobianMode::Analytic)] {
            p.lsq.jacobian = mode;
            let t_evaluate = best_of(3, || {
                evaluate(&p).unwrap();
            });
            let t_jac = t_evaluate.saturating_sub(2 * t_eval);
            let t_fit = best_of(3, || {
                fit(&p).unwrap();
            });
            let r = fit(&p).unwrap();
            let speedup = match fd_fit {
                None => {
                    fd_fit = Some(t_fit);
                    String::new()
                }
                Some(fd) => format!("{:.2}x", fd.as_secs_f64() / t_fit.as_secs_f64()),
            };
            println!(
                "{:<18} {:>3} {:>4} {:>6} {:>8.3} {:>8.3} {:>8.1} {:>8.1} {:>4} {:>6} {:>8.3} {:>7}",
                case.name,
                label,
                r.n_free,
                r.n_residuals,
                ms(t_eval),
                ms(t_jac),
                t_jac.as_secs_f64() / t_eval.as_secs_f64(),
                ms(t_fit),
                r.iterations,
                r.evals,
                r.chi2_reduced,
                speedup
            );
            rows.push(json!({
                "case": case.name, "jacobian": label, "n_free": r.n_free, "m": r.n_residuals,
                "eval_ms": ms(t_eval), "jac_ms": ms(t_jac), "fit_ms": ms(t_fit),
                "iterations": r.iterations, "evals": r.evals, "chi2_reduced": r.chi2_reduced,
                "chi2_bits": format!("{:016x}", r.chi2.to_bits()),
            }));
        }
    }
    if let Some(path) = json_out {
        std::fs::write(&path, serde_json::to_string_pretty(&rows).unwrap()).unwrap();
        println!("written {path}");
    }
    ExitCode::SUCCESS
}
