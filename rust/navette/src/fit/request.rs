// SPDX-License-Identifier: LGPL-3.0-or-later
//! The fit as data (`docs/plans/general_fit_plan.md` §3.10.1).
//!
//! A fit request is a program document (`config`) whose `fit` section says
//! what to fit: the model is the program's own architect (or its one
//! structure) and material library, so a request and the program file a
//! fit is saved in are one format.
//!
//! ```json
//! {"schema_version": 4, "kind": "program", "sections": {
//!   "materials": [...], "groups": [...], "structures": [...],
//!   "fit": {
//!     "environments": [{"name": "bare", "ambient": "air", "substrate": "glass",
//!                       "stack": []}],
//!     "data": [{"id": "R", "channel": "R_avg", "angle_deg": 8,
//!               "wavelengths": [...], "values": [...], "sigma": 0.002}],
//!     "nk_data": [{"id": "ellips", "material": "TiO2", "wavelengths": [...],
//!                  "columns": [{"quantity": "n", "values": [...], "sigma": 0.01}]}],
//!     "parameters": [{"addr": "layer[main:0].thickness", "lb": 50, "ub": 80},
//!                    {"select": "kind=thickness"}],
//!     "ties": [{"driver": "layer[a:1].thickness",
//!               "targets": [{"target": "layer[b:1].thickness"}]}],
//!     "sets": [{"name": "H", "members": ["tag=H"]}],
//!     "strategy": {"name": "...", "stages": [...]},
//!     "coherence_mode": 0,
//!     "lsq": {"backend": "trf", "max_iterations": 200}
//!   }}}
//! ```
//!
//! - `environments`: extra samples, each a chain of `stack` items
//!   (`{"structure": s, "inverted": false, "repeat": 1}` or
//!   `{"layer": material, "thickness": d}`) between an optional `ambient`
//!   and `substrate` half-space. The program's own chain is `"default"`;
//!   program structures are films only, so a sample is usually declared
//!   here.
//! - `data`: solver-convention spectra; exactly one of `angle_deg` and
//!   `angle_rad`; `sigma` one number or one per point; `env` defaults to
//!   the first declared environment (`"default"` when none is). Instrument
//!   data goes through [`to_solver`](super::convert::to_solver) first.
//! - `parameters`: `addr` names one scalar; `select` expands a selector
//!   over [`FitModel::addresses`] (addresses declared elsewhere in the list
//!   are skipped). The other keys (`value`, `lb`, `ub`, `scale`, `log`,
//!   `fixed`, `tags`, `prior: [p0, sigma]`) apply to every expanded one.
//! - `strategy`: stages to run; without one, everything free is fitted in
//!   one stage.

use std::collections::HashSet;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::config::{LoadedProgram, load_program_json};
use crate::lsq::backends::LsqBackend;
use crate::lsq::lm::LmConfig;
use crate::smatrix::angle::Radians;
use crate::structure::Layer;
use crate::structure::enums::LayerType;

use super::data::{NkColumn, NkData, NkQuantity, SpectrumData};
use super::model::{FitModel, Segment};
use super::params::{Parameter, Tie, Transform};
use super::sets::{ParamSet, select};
use super::strategy::{StageLsq, Strategy, StrategyReport, run_strategy_watched};
use super::{FitProblem, FitReport, Progress, Watcher, evaluate, fit_watched};

/// A program's `fit` section.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FitSection {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub environments: Vec<EnvironmentSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub data: Vec<DataSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub nk_data: Vec<NkSpec>,
    pub parameters: Vec<ParamSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub ties: Vec<TieSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub sets: Vec<ParamSet>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strategy: Option<Strategy>,
    #[serde(default)]
    pub coherence_mode: i32,
    #[serde(default)]
    pub lsq: StageLsq,
}

/// One extra sample.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EnvironmentSpec {
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ambient: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub substrate: Option<String>,
    #[serde(default)]
    pub stack: Vec<StackItem>,
}

/// A structure of the model, or one held film.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StackItem {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub structure: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub layer: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub thickness: Option<f64>,
    #[serde(default)]
    pub inverted: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub repeat: Option<usize>,
}

/// One number for every point, or one per point.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Sigma {
    One(f64),
    Each(Vec<f64>),
}

impl Sigma {
    fn to_vec(&self) -> Vec<f64> {
        match self {
            Sigma::One(s) => vec![*s],
            Sigma::Each(v) => v.clone(),
        }
    }
}

/// One measured spectrum in solver conventions.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DataSpec {
    pub id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub env: Option<String>,
    pub channel: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub angle_deg: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub angle_rad: Option<f64>,
    pub wavelengths: Vec<f64>,
    pub values: Vec<f64>,
    pub sigma: Sigma,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub weight: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub window: Option<(f64, f64)>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
}

/// One column of an nk dataset: `n`, `k`, `eps1` or `eps2`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NkColumnSpec {
    pub quantity: String,
    pub values: Vec<f64>,
    pub sigma: Sigma,
}

/// Measured optical constants of one material.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NkSpec {
    pub id: String,
    pub material: String,
    pub wavelengths: Vec<f64>,
    pub columns: Vec<NkColumnSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub weight: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub window: Option<(f64, f64)>,
}

/// One parameter (`addr`) or a selector over the model (`select`).
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ParamSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub addr: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub select: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub value: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lb: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ub: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scale: Option<f64>,
    #[serde(default)]
    pub log: bool,
    #[serde(default)]
    pub fixed: bool,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub tags: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub prior: Option<(f64, f64)>,
}

/// `target = factor * driver + offset` for each target.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TieSpec {
    pub driver: String,
    pub targets: Vec<TieTarget>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TieTarget {
    pub target: String,
    #[serde(default = "one")]
    pub factor: f64,
    #[serde(default)]
    pub offset: f64,
}

fn one() -> f64 {
    1.0
}

/// A request, built: the problem, the strategy to run on it (if any) and
/// the user parameter sets.
#[derive(Debug, Clone)]
pub struct FitRequest {
    pub problem: FitProblem,
    pub strategy: Option<Strategy>,
    pub sets: Vec<ParamSet>,
}

/// Build the request a program document describes; it must have a `fit`
/// section.
pub fn request_from_json(text: &str) -> Result<FitRequest, String> {
    let raw: Value =
        serde_json::from_str(text).map_err(|e| format!("fit request: invalid JSON: {e}"))?;
    // Materials are specs, evaluated on each dataset's own wavelengths;
    // the loader's grid only checks array entries, so any grid will do.
    let grid = data_grid(&raw);
    let program = load_program_json(text, &grid)?;
    request_from_program(&program)
}

/// Every data wavelength, sorted and deduplicated; `[550]` when none.
fn data_grid(raw: &Value) -> Vec<f64> {
    let fit = &raw["sections"]["fit"];
    let mut w: Vec<f64> = ["data", "nk_data"]
        .iter()
        .filter_map(|k| fit[k].as_array())
        .flatten()
        .filter_map(|d| d["wavelengths"].as_array())
        .flatten()
        .filter_map(Value::as_f64)
        .filter(|x| x.is_finite())
        .collect();
    w.sort_by(f64::total_cmp);
    w.dedup();
    if w.is_empty() { vec![550.0] } else { w }
}

/// Build the request from a loaded program.
pub fn request_from_program(program: &LoadedProgram) -> Result<FitRequest, String> {
    let section = program
        .fit
        .as_ref()
        .ok_or("fit request: the program has no 'fit' section")?;
    let has_stack = program.architect.is_some() || !program.structures.is_empty();
    let mut model = if has_stack {
        FitModel::from_program(program)?
    } else {
        let materials = program
            .materials
            .as_ref()
            .ok_or("fit request: the program has neither a stack nor materials")?
            .entries()
            .clone();
        FitModel::from_materials(materials)
    };
    for env in &section.environments {
        model.add_environment(&env.name, segments(env)?)?;
    }

    // A program's structures are films only, so its own chain is rarely a
    // sample: data that names no environment is measured on the first one
    // the section declares, or on the program's chain when it declares none.
    let fallback = section
        .environments
        .first()
        .map_or(super::DEFAULT_ENV, |e| e.name.as_str());
    let data = section
        .data
        .iter()
        .map(|d| spectrum(d, fallback))
        .collect::<Result<Vec<_>, _>>()?;
    let nk_data = section
        .nk_data
        .iter()
        .map(nk)
        .collect::<Result<Vec<_>, _>>()?;
    let parameters = parameters(&section.parameters, &model, &section.sets)?;
    let ties = section
        .ties
        .iter()
        .map(|t| {
            t.targets.iter().try_fold(Tie::parse(&t.driver)?, |tie, x| {
                tie.follow(&x.target, x.factor, x.offset)
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut problem = FitProblem::new(model, data, parameters);
    problem.nk_data = nk_data;
    problem.ties = ties;
    problem.coherence_mode = section.coherence_mode;
    problem.lsq = lm_config(&section.lsq)?;
    Ok(FitRequest {
        problem,
        strategy: section.strategy.clone(),
        sets: section.sets.clone(),
    })
}

fn segments(env: &EnvironmentSpec) -> Result<Vec<Segment>, String> {
    let bad = |why: String| format!("environment {:?}: {why}", env.name);
    let half = |material: &str, role: LayerType| {
        Segment::Fixed(Layer {
            layer_type: role,
            ..Layer::film(0.0, material)
        })
    };
    let mut out = Vec::new();
    if let Some(a) = &env.ambient {
        out.push(half(a, LayerType::Ambient));
    }
    for item in &env.stack {
        match (&item.structure, &item.layer) {
            (Some(s), None) => {
                if item.thickness.is_some() {
                    return Err(bad(format!("structure {s:?} takes no thickness")));
                }
                out.push(Segment::Structure {
                    name: s.clone(),
                    inverted: item.inverted,
                    repeat: item.repeat.unwrap_or(1),
                });
            }
            (None, Some(m)) => {
                if item.inverted || item.repeat.is_some() {
                    return Err(bad(format!(
                        "layer {m:?} takes only a thickness, not inverted or repeat"
                    )));
                }
                let d = item
                    .thickness
                    .ok_or_else(|| bad(format!("layer {m:?} needs a thickness")))?;
                out.push(Segment::Fixed(Layer::film(d, m)));
            }
            _ => {
                return Err(bad(
                    "each stack item names exactly one of 'structure' and 'layer'".into(),
                ));
            }
        }
    }
    if let Some(s) = &env.substrate {
        out.push(half(s, LayerType::Substrate));
    }
    Ok(out)
}

fn angle(id: &str, deg: Option<f64>, rad: Option<f64>) -> Result<Radians, String> {
    let bad = |why: &str| format!("dataset {id:?}: {why}");
    match (deg, rad) {
        (Some(d), None) => Radians::from_degrees(d).map_err(|e| bad(&e.to_string())),
        (None, Some(r)) => Radians::new(r).map_err(|e| bad(&e.to_string())),
        _ => Err(bad("give exactly one of angle_deg and angle_rad")),
    }
}

pub(crate) fn spectrum(d: &DataSpec, fallback: &str) -> Result<SpectrumData, String> {
    let mut s = SpectrumData::new(
        d.id.clone(),
        d.channel.clone(),
        angle(&d.id, d.angle_deg, d.angle_rad)?,
        d.wavelengths.clone(),
        d.values.clone(),
        1.0,
    );
    s.sigma = d.sigma.to_vec();
    s.env = d.env.clone().unwrap_or_else(|| fallback.to_string());
    if let Some(w) = d.weight {
        s.weight = w;
    }
    s.window = d.window;
    s.source = d.source.clone();
    s.check()?;
    Ok(s)
}

pub(crate) fn nk(d: &NkSpec) -> Result<NkData, String> {
    let columns = d
        .columns
        .iter()
        .map(|c| {
            let quantity = match c.quantity.as_str() {
                "n" => NkQuantity::N,
                "k" => NkQuantity::K,
                "eps1" => NkQuantity::Eps1,
                "eps2" => NkQuantity::Eps2,
                other => {
                    return Err(format!(
                        "nk dataset {:?}: quantity {other:?}; one of n, k, eps1, eps2",
                        d.id
                    ));
                }
            };
            Ok(NkColumn {
                quantity,
                values: c.values.clone(),
                sigma: c.sigma.to_vec(),
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let mut out = NkData::new(
        d.id.clone(),
        d.material.clone(),
        d.wavelengths.clone(),
        columns,
    );
    if let Some(w) = d.weight {
        out.weight = w;
    }
    out.window = d.window;
    out.check()?;
    Ok(out)
}

fn dressed(mut p: Parameter, s: &ParamSpec) -> Parameter {
    p.value = s.value;
    p.lb = s.lb;
    p.ub = s.ub;
    p.scale = s.scale;
    if s.log {
        p.transform = Transform::Log;
    }
    p.fixed = s.fixed;
    p.tags = s.tags.clone();
    p.prior = s.prior;
    p
}

/// Declarations to parameters: addresses as given, selectors expanded over
/// the model's addresses not declared by address anywhere in the list.
pub(crate) fn parameters(
    specs: &[ParamSpec],
    model: &FitModel,
    sets: &[ParamSet],
) -> Result<Vec<Parameter>, String> {
    let mut named = HashSet::new();
    for s in specs {
        if let Some(a) = &s.addr {
            let canonical = Parameter::parse(a)?.addr.to_string();
            if !named.insert(canonical) {
                return Err(format!("parameter {a:?} is declared twice"));
            }
        }
    }
    let mut taken = named.clone();
    let mut out = Vec::new();
    for s in specs {
        match (&s.addr, &s.select) {
            (Some(a), None) => out.push(dressed(Parameter::parse(a)?, s)),
            (None, Some(sel)) => {
                let pool: Vec<Parameter> = model
                    .addresses()
                    .into_iter()
                    .filter(|a| !taken.contains(&a.to_string()))
                    .map(Parameter::new)
                    .collect();
                let picked = select(std::slice::from_ref(sel), &pool, sets)
                    .map_err(|e| format!("parameters: {e}"))?;
                for (p, hit) in pool.into_iter().zip(picked) {
                    if hit {
                        taken.insert(p.addr.to_string());
                        out.push(dressed(p, s));
                    }
                }
            }
            _ => {
                return Err(
                    "each parameter declaration has exactly one of 'addr' and 'select'".into(),
                );
            }
        }
    }
    Ok(out)
}

/// The default solver settings with `s`'s overrides.
pub(crate) fn lm_config(s: &StageLsq) -> Result<LmConfig, String> {
    let mut c = LmConfig::default();
    if let Some(b) = &s.backend {
        c.backend = LsqBackend::parse(b)?;
    }
    if let Some(v) = s.max_iterations {
        c.max_iterations = v;
    }
    if let Some(v) = s.max_evals {
        c.max_evals = v;
    }
    if let Some(v) = s.ftol {
        c.ftol = v;
    }
    if let Some(v) = s.xtol {
        c.xtol = v;
    }
    if let Some(v) = s.gtol {
        c.gtol = v;
    }
    Ok(c)
}

// ---------------------------------------------------------------------------
// Reports as JSON
// ---------------------------------------------------------------------------

/// `NaN` and infinities become `null`.
pub(crate) fn num(x: f64) -> Value {
    if x.is_finite() { json!(x) } else { Value::Null }
}

/// A [`FitReport`] as JSON (the model left out).
pub fn report_json(r: &FitReport) -> Value {
    json!({
        "parameters": r.parameters.iter().map(|p| json!({
            "addr": p.addr,
            "start": num(p.start),
            "value": num(p.value),
            "std_error": p.std_error.map(num),
            "lb": num(p.lb),
            "ub": num(p.ub),
            "fixed": p.fixed,
            "at_bound": p.at_bound,
        })).collect::<Vec<_>>(),
        "datasets": r.datasets.iter().map(|d| json!({
            "id": d.id,
            "source": d.source,
            "n_points": d.n_points,
            "chi2": num(d.chi2),
            "rms": num(d.rms),
        })).collect::<Vec<_>>(),
        "chi2": num(r.chi2),
        "chi2_reduced": num(r.chi2_reduced),
        "aic": num(r.aic),
        "bic": num(r.bic),
        "n_residuals": r.n_residuals,
        "n_free": r.n_free,
        "covariance": r.covariance.as_ref().map(|c| c.iter().map(|x| num(*x)).collect::<Vec<_>>()),
        "correlation": r.correlation.as_ref().map(|c| c.iter().map(|x| num(*x)).collect::<Vec<_>>()),
        "covariance_note": r.covariance_note,
        "termination": format!("{:?}", r.termination),
        "iterations": r.iterations,
        "evals": r.evals,
        "backend": r.backend.as_str(),
    })
}

/// A [`StrategyReport`] as JSON: its stages and its end report.
pub fn strategy_report_json(r: &StrategyReport) -> Value {
    json!({
        "name": r.name,
        "evals": r.evals,
        "stages": r.stages.iter().map(|s| json!({
            "name": s.name,
            "free": s.free,
            "chi2_start": num(s.chi2_start),
            "chi2": num(s.report.chi2),
            "runs": s.runs,
            "evals": s.evals,
            "iterations": s.report.iterations,
            "termination": format!("{:?}", s.report.termination),
            "accepted": s.accepted,
        })).collect::<Vec<_>>(),
        "report": report_json(&r.report),
    })
}

/// Measured and simulated curves of `problem` on `model`, at each
/// dataset's kept points: what a plot of the fit draws.
pub fn curves_json(problem: &FitProblem, model: &FitModel) -> Result<Value, String> {
    let mut at = problem.clone();
    at.model = model.clone();
    for q in &mut at.parameters {
        q.fixed = true;
        q.value = None;
    }
    let sims = at.simulate()?;
    let mut out = Vec::new();
    for (d, sim) in at.data.iter().zip(&sims) {
        let kept = d.kept();
        out.push(json!({
            "id": d.id,
            "kind": "spectrum",
            "env": d.env,
            "channel": d.channel,
            "angle_deg": d.angle.get().to_degrees(),
            "wavelengths": kept.iter().map(|&i| d.wavelengths[i]).collect::<Vec<_>>(),
            "measured": kept.iter().map(|&i| d.values[i]).collect::<Vec<_>>(),
            "sigma": kept.iter().map(|&i| d.sigma_at(i)).collect::<Vec<_>>(),
            "simulated": sim,
        }));
    }
    for (d, sim) in at.nk_data.iter().zip(&sims[at.data.len()..]) {
        let kept = d.kept();
        let k = kept.len();
        for (j, c) in d.columns.iter().enumerate() {
            let quantity = match c.quantity {
                NkQuantity::N => "n",
                NkQuantity::K => "k",
                NkQuantity::Eps1 => "eps1",
                NkQuantity::Eps2 => "eps2",
            };
            out.push(json!({
                "id": d.id,
                "kind": "nk",
                "material": d.material,
                "quantity": quantity,
                "wavelengths": kept.iter().map(|&i| d.wavelengths[i]).collect::<Vec<_>>(),
                "measured": kept.iter().map(|&i| c.values[i]).collect::<Vec<_>>(),
                "sigma": kept.iter().map(|&i| c.sigma_at(i)).collect::<Vec<_>>(),
                "simulated": sim.get(j * k..(j + 1) * k).unwrap_or_default(),
            }));
        }
    }
    Ok(Value::Array(out))
}

/// Run the fit a program document describes and report it as JSON:
/// `{"strategy": <strategy report or null>, "report": ..., "curves": ...}`.
/// Without a strategy, every free parameter is fitted in one stage named
/// `fit`, and the watcher sees that one stage.
pub fn fit_program(text: &str, watch: Option<&Watcher<'_>>) -> Result<Value, String> {
    let req = request_from_json(text)?;
    let (strategy, report) = match &req.strategy {
        Some(s) => {
            let r = run_strategy_watched(&req.problem, s, &req.sets, watch)?;
            (Some(strategy_report_json(&r)), r.report)
        }
        None => {
            let free: Vec<String> = req
                .problem
                .parameters
                .iter()
                .filter(|p| !p.fixed)
                .map(|p| p.addr.to_string())
                .collect();
            if let Some(w) = watch {
                w.tell(&Progress::StageStart {
                    stage: 0,
                    name: "fit".into(),
                    free,
                });
            }
            let r = fit_watched(&req.problem, watch.map(|w| (w, 0)));
            if let Some(w) = watch {
                w.check()?;
            }
            let r = r?;
            if let Some(w) = watch {
                w.tell(&Progress::StageEnd {
                    stage: 0,
                    name: "fit".into(),
                    chi2: r.chi2,
                    accepted: true,
                });
                w.check()?;
            }
            (None, r)
        }
    };
    Ok(json!({
        "strategy": strategy,
        "report": report_json(&report),
        "curves": curves_json(&req.problem, &report.model)?,
    }))
}

/// The report and curves at the start values, nothing moved.
pub fn evaluate_program(text: &str) -> Result<Value, String> {
    let req = request_from_json(text)?;
    let report = evaluate(&req.problem)?;
    Ok(json!({
        "strategy": Value::Null,
        "report": report_json(&report),
        "curves": curves_json(&req.problem, &report.model)?,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bare(fit: Value) -> String {
        json!({
            "schema_version": 4,
            "kind": "program",
            "sections": {
                "materials": [
                    {"name": "air", "model": "Konstant", "params": {"n": 1.0}},
                    {"name": "glass", "model": "Konstant", "params": {"n": 1.52}},
                    {"name": "TiO2", "model": "Cauchy",
                     "params": {"A": 2.2, "B": 0.02, "C": 0.001}},
                ],
                "structures": [{"label": "main", "layers": [
                    {"material_code": "TiO2", "thickness_nm": 80.0},
                ], "groups": []}],
                "fit": fit,
            }
        })
        .to_string()
    }

    /// `fit` plus the one sample every test measures.
    fn program(mut fit: Value) -> String {
        fit["environments"] = json!([{"name": "sample", "ambient": "air",
                                      "substrate": "glass",
                                      "stack": [{"structure": "main"}]}]);
        bare(fit)
    }

    fn data(values: Vec<f64>) -> Value {
        let wl: Vec<f64> = (0..41).map(|i| 400.0 + 10.0 * i as f64).collect();
        json!([{"id": "R", "channel": "R_avg", "angle_deg": 8.0,
                "wavelengths": wl, "values": values, "sigma": 0.002}])
    }

    #[test]
    fn select_expands_and_skips_named_addresses() {
        let text = program(json!({
            "data": data(vec![0.1; 41]),
            "parameters": [
                {"select": "material[TiO2].*", "lb": 0.0, "ub": 5.0},
                {"addr": "material[TiO2].A", "lb": 2.0, "ub": 2.5},
                {"select": "kind=thickness"},
            ],
        }));
        let req = request_from_json(&text).unwrap();
        let names: Vec<String> = req
            .problem
            .parameters
            .iter()
            .map(|p| p.addr.to_string())
            .collect();
        assert_eq!(
            names,
            [
                "material[TiO2].B",
                "material[TiO2].C",
                "material[TiO2].A",
                "layer[main:0].thickness"
            ]
        );
        assert_eq!(req.problem.parameters[2].lb, Some(2.0));
        assert_eq!(req.problem.parameters[0].ub, Some(5.0));
    }

    #[test]
    fn a_program_fits_its_own_simulation_back() {
        // Simulate at 92 nm, fit from 80 nm.
        let mut truth: FitRequest = request_from_json(&program(json!({
            "data": data(vec![0.0; 41]),
            "parameters": [{"addr": "layer[main:0].thickness", "value": 92.0, "fixed": true}],
        })))
        .unwrap();
        truth.problem.parameters[0].fixed = true;
        let sim = truth.problem.simulate().unwrap().remove(0);
        let text = program(json!({
            "data": data(sim),
            "parameters": [{"addr": "layer[main:0].thickness", "lb": 50.0, "ub": 120.0}],
        }));
        let seen = std::sync::Mutex::new(Vec::new());
        let w = Watcher::new(|p: &Progress| {
            seen.lock().unwrap().push(p.clone());
            Ok(())
        });
        let out = fit_program(&text, Some(&w)).unwrap();
        let d = out["report"]["parameters"][0]["value"].as_f64().unwrap();
        assert!((d - 92.0).abs() < 1e-3, "{d}");
        assert_eq!(out["curves"][0]["simulated"].as_array().unwrap().len(), 41);
        drop(w);
        let seen = seen.into_inner().unwrap();
        assert!(matches!(seen.first(), Some(Progress::StageStart { .. })));
        assert!(matches!(seen.last(), Some(Progress::StageEnd { .. })));
        assert!(seen.iter().any(|p| matches!(p, Progress::Improved { .. })));

        // A callback that refuses stops the fit with its message.
        let w = Watcher::new(|p: &Progress| match p {
            Progress::Improved { .. } => Err("enough".into()),
            _ => Ok(()),
        });
        let e = fit_program(&text, Some(&w)).unwrap_err();
        assert!(e.contains("enough"), "{e}");
    }

    #[test]
    fn bad_declarations_name_themselves() {
        let e = request_from_json(&program(json!({
            "data": data(vec![0.1; 41]),
            "parameters": [{"addr": "layer[main:0].thickness", "select": "*"}],
        })))
        .unwrap_err();
        assert!(e.contains("exactly one of 'addr' and 'select'"), "{e}");
        let e = request_from_json(&program(json!({
            "data": [{"id": "R", "channel": "R_avg", "wavelengths": [500.0],
                      "values": [0.1], "sigma": 0.01}],
            "parameters": [],
        })))
        .unwrap_err();
        assert!(e.contains("angle_deg and angle_rad"), "{e}");
        let e = request_from_json(&program(json!({"parameters": [], "bogus": 1}))).unwrap_err();
        assert!(e.contains("bogus"), "{e}");
    }
}
