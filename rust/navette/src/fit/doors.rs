// SPDX-License-Identifier: LGPL-3.0-or-later
//! JSON doors over the fit API: what `navette.fit` calls (plan §3.10.2).
//!
//! Each door takes JSON text and returns a JSON value, so the Python side
//! is a thin layer of dataclasses and the shapes are tested here, without
//! Python. A program document with a `fit` section ([`super::request`]) is
//! the problem everywhere; strategies are a list whose items are either a
//! strategy (`{"name", "stages"}`) or a generator
//! (`{"orderings": {...}}`, [`StrategyGen`]).

use std::collections::{BTreeMap, HashMap};

use serde_json::{Value, json};

use crate::structure::{Entry, MaterialSpec};

use super::convert::{InstrumentConvention, RawDataset, from_solver, to_solver};
use super::data::SpectrumData;
use super::explore::{StrategyGen, benchmark, draw_truth, generate, tournament};
use super::model::FitModel;
use super::request::{
    DataSpec, NkSpec, ParamSpec, Sigma, curves_json, lm_config, nk, num, parameters, report_json,
    request_from_json, strategy_report_json,
};
use super::strategy::{StageLsq, Strategy};
use super::{FitProblem, dispersion};

fn parse<T: serde::de::DeserializeOwned>(text: &str, what: &str) -> Result<T, String> {
    serde_json::from_str(text).map_err(|e| format!("{what}: {e}"))
}

/// Strategies from a list of strategies and generators, with the notes of
/// what was sampled.
pub fn strategies_from_json(
    text: &str,
    cap: usize,
) -> Result<(Vec<Strategy>, Vec<String>), String> {
    let items: Vec<Value> = parse(text, "strategies")?;
    let mut gens = Vec::new();
    for (i, item) in items.into_iter().enumerate() {
        if item.get("stages").is_some() {
            let s: Strategy =
                serde_json::from_value(item).map_err(|e| format!("strategies[{i}]: {e}"))?;
            gens.push(StrategyGen::Explicit(vec![s]));
        } else {
            gens.push(serde_json::from_value(item).map_err(|e| format!("strategies[{i}]: {e}"))?);
        }
    }
    let made = generate(&gens, cap)?;
    Ok((made.strategies, made.notes))
}

/// `{"strategies": [...], "notes": [...]}`.
pub fn generate_json(text: &str, cap: usize) -> Result<Value, String> {
    let (strategies, notes) = strategies_from_json(text, cap)?;
    Ok(json!({"strategies": strategies, "notes": notes}))
}

/// Every strategy on the program's problem, ranked by χ²:
/// `{"entries": [...], "notes": [...]}`.
pub fn tournament_json(program: &str, strategies: &str, cap: usize) -> Result<Value, String> {
    let req = request_from_json(program)?;
    let (list, notes) = strategies_from_json(strategies, cap)?;
    let entries: Vec<Value> = tournament(&req.problem, &list, &req.sets)
        .iter()
        .map(|e| {
            json!({
                "name": e.name,
                "rank": e.rank,
                "chi2": num(e.chi2),
                "chi2_reduced": num(e.chi2_reduced),
                "bic": num(e.bic),
                "rel_uncertainty": num(e.rel_uncertainty),
                "evals": e.evals,
                "time_s": e.time.as_secs_f64(),
                "error": e.error,
                "report": e.report.as_ref().map(strategy_report_json),
            })
        })
        .collect();
    Ok(json!({"entries": entries, "notes": notes}))
}

/// Every strategy on a truth drawn per seed, ranked by recovery.
pub fn benchmark_json(
    program: &str,
    strategies: &str,
    seeds: &[u64],
    cap: usize,
) -> Result<Value, String> {
    let req = request_from_json(program)?;
    let (list, notes) = strategies_from_json(strategies, cap)?;
    let b = benchmark(&req.problem, &list, &req.sets, seeds)?;
    let entries: Vec<Value> = b
        .entries
        .iter()
        .map(|e| {
            json!({
                "name": e.name,
                "rank": e.rank,
                "mean_recovery": num(e.mean_recovery),
                "worst_recovery": num(e.worst_recovery),
                "mean_chi2_reduced": num(e.mean_chi2_reduced),
                "mean_evals": num(e.mean_evals),
                "successes": e.successes,
                "failures": e.failures,
                "runs": e.runs.iter().map(|r| json!({
                    "seed": r.seed,
                    "chi2_reduced": num(r.chi2_reduced),
                    "recovery": num(r.recovery),
                    "max_z": num(r.max_z),
                    "within_2sigma": r.within_2sigma,
                    "evals": r.evals,
                    "time_s": r.time.as_secs_f64(),
                    "error": r.error,
                })).collect::<Vec<_>>(),
            })
        })
        .collect();
    Ok(json!({"seeds": b.seeds, "n_free": b.n_free, "entries": entries, "notes": notes}))
}

/// A spectrum as the request's `data` item. The default environment is
/// left unnamed, so the item lands where unnamed data does (a request's
/// first declared sample).
fn data_spec(d: &SpectrumData) -> DataSpec {
    DataSpec {
        id: d.id.clone(),
        env: (d.env != super::DEFAULT_ENV).then(|| d.env.clone()),
        channel: d.channel.clone(),
        angle_deg: None,
        angle_rad: Some(d.angle.get()),
        wavelengths: d.wavelengths.clone(),
        values: d.values.clone(),
        sigma: match d.sigma.as_slice() {
            [s] => Sigma::One(*s),
            v => Sigma::Each(v.to_vec()),
        },
        weight: (d.weight != 1.0).then_some(d.weight),
        window: d.window,
        source: d.source.clone(),
    }
}

/// A truth drawn from the program's error laws:
/// `{"seed", "values": [[addr, value], ...], "data": [...]}`, the data
/// in the request's own shape, ready to replace the section's.
pub fn truth_json(program: &str, seed: u64) -> Result<Value, String> {
    let req = request_from_json(program)?;
    let t = draw_truth(&req.problem, seed)?;
    Ok(json!({
        "seed": t.seed,
        "values": t.values,
        "data": t.problem.data.iter().map(data_spec).collect::<Vec<_>>(),
    }))
}

/// A preset name, or a convention as a JSON object.
fn convention(text: &str) -> Result<InstrumentConvention, String> {
    if text.trim_start().starts_with('{') {
        parse(text, "instrument convention")
    } else {
        InstrumentConvention::preset(text)
    }
}

/// An instrument record to solver-convention `data` items.
pub fn to_solver_json(raw: &str, conv: &str) -> Result<Value, String> {
    let raw: RawDataset = parse(raw, "raw dataset")?;
    let data = to_solver(&raw, &convention(conv)?)?;
    Ok(json!(data.iter().map(data_spec).collect::<Vec<_>>()))
}

/// Solver-convention `data` items (one spectrum, or a Psi and a Delta) as
/// the instrument would write them.
pub fn from_solver_json(id: &str, data: &str, conv: &str) -> Result<Value, String> {
    let specs: Vec<DataSpec> = parse(data, "data")?;
    let data = specs
        .iter()
        .map(|d| super::request::spectrum(d, super::DEFAULT_ENV))
        .collect::<Result<Vec<_>, _>>()?;
    let refs: Vec<&SpectrumData> = data.iter().collect();
    let raw = from_solver(id, &refs, &convention(conv)?)?;
    serde_json::to_value(raw).map_err(|e| e.to_string())
}

/// One material spec against nk data:
/// `{"spec", "report", "curves"}`. `parameters` (declarations, as in a
/// request) `None` frees every number of the spec.
pub fn fit_dispersion_json(
    spec: &str,
    nk_data: &str,
    parameters_json: Option<&str>,
    lsq: &str,
) -> Result<Value, String> {
    let spec: MaterialSpec = parse(spec, "spec")?;
    let data: Vec<NkSpec> = parse(nk_data, "nk_data")?;
    let data = data.iter().map(nk).collect::<Result<Vec<_>, _>>()?;
    let lsq: StageLsq = parse(lsq, "lsq")?;
    let name = data
        .first()
        .ok_or("fit_dispersion: no nk data")?
        .material
        .clone();
    let model_of = |s: &MaterialSpec| {
        FitModel::from_materials(HashMap::from([(name.clone(), Entry::Spec(s.clone()))]))
    };
    let params = match parameters_json {
        None => None,
        Some(text) => {
            let decl: Vec<ParamSpec> = parse(text, "parameters")?;
            Some(parameters(&decl, &model_of(&spec), &[])?)
        }
    };
    let (fitted, report) =
        dispersion::fit_dispersion(&spec, data.clone(), params, &lm_config(&lsq)?)?;
    let at = FitProblem::dispersion(model_of(&fitted), data, Vec::new());
    Ok(json!({
        "spec": fitted,
        "report": report_json(&report),
        "curves": curves_json(&at, &at.model)?,
    }))
}

/// Rank model candidates for a `Table` spec by BIC: a list of
/// `{"model", "spec", "aic", "bic", "error", "report"}`, best first.
pub fn table_to_model_json(table: &str, candidates: &str, lsq: &str) -> Result<Value, String> {
    let table: MaterialSpec = parse(table, "table")?;
    let candidates: Vec<MaterialSpec> = parse(candidates, "candidates")?;
    let lsq: StageLsq = parse(lsq, "lsq")?;
    let ranked = dispersion::table_to_model(&table, candidates, &lm_config(&lsq)?)?;
    Ok(json!(
        ranked
            .iter()
            .map(|c| json!({
                "model": c.model,
                "spec": c.spec,
                "aic": num(c.aic),
                "bic": num(c.bic),
                "error": c.error,
                "report": c.report.as_ref().map(report_json),
            }))
            .collect::<Vec<_>>()
    ))
}

/// The preset conventions by name.
pub fn conventions_json() -> Value {
    let all: BTreeMap<String, Value> = InstrumentConvention::preset_names()
        .into_iter()
        .filter_map(|n| {
            let c = InstrumentConvention::preset(&n).ok()?;
            Some((n, serde_json::to_value(c).ok()?))
        })
        .collect();
    json!(all)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strategies_mix_explicit_and_generated() {
        let (list, notes) = strategies_from_json(
            r#"[{"name": "all", "stages": [{"name": "a", "free": ["*"]}]},
                {"orderings": {"sets": ["thicknesses", "groups"]}}]"#,
            24,
        )
        .unwrap();
        let names: Vec<&str> = list.iter().map(|s| s.name.as_str()).collect();
        assert_eq!(
            names,
            [
                "all",
                "thicknesses > groups | polish",
                "groups > thicknesses | polish"
            ]
        );
        assert!(notes.is_empty());
        let e = strategies_from_json(r#"[{"bogus": {}}]"#, 24).unwrap_err();
        assert!(e.contains("strategies[0]"), "{e}");
    }

    #[test]
    fn a_convention_is_a_preset_or_an_object() {
        let all = conventions_json();
        let name = all.as_object().unwrap().keys().next().unwrap().clone();
        assert!(convention(&name).is_ok());
        assert!(convention(&all[&name].to_string()).is_ok());
        assert!(convention("no such preset").is_err());
    }
}
