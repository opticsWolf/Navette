// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::dispersion — fitting a material spec to measured optical constants.
//!
//! Not a separate optimizer: both doors build an ordinary [`FitProblem`]
//! over a materials-only [`FitModel`] and call [`fit`].
//!
//! - [`fit_dispersion`]: one spec against nk data; returns the fitted spec
//!   and the report. A fitted spec goes back into a library with
//!   `SpecProvider::upsert`.
//! - [`table_to_model`]: fit each candidate model to a `Table` material on
//!   the table's own grid and rank the candidates by BIC.

use std::collections::HashMap;

use serde_json::Value;

use crate::lsq::lm::LmConfig;
use crate::structure::providers::Entry;
use crate::structure::specs::MaterialSpec;

use super::data::{NkColumn, NkData, NkQuantity};
use super::model::FitModel;
use super::params::{ParamAddr, Parameter, PathSeg};
use super::{FitProblem, FitReport, fit};

/// The mixers: their default free parameters are their own scalars, not
/// their constituents'.
const MIXERS: [&str; 7] = [
    "Bruggeman",
    "MaxwellGarnett",
    "Looyenga",
    "Lichtenecker",
    "MoriTanaka",
    "PowerLaw",
    "Roughness",
];

/// Fit `spec` to `data`, every dataset naming the same material (that name
/// is the spec's in the addresses: `material[<name>].C1`).
///
/// `parameters: None` frees every number of the spec: each scalar and each
/// oscillator entry, except a `Table`'s data and, for a mixer (`fraction`,
/// `L`, `alpha`), its constituents, which stay held unless named. Bounds
/// default to the per-model table (`fraction` in `[0, 1]`, widths
/// positive, ...).
pub fn fit_dispersion(
    spec: &MaterialSpec,
    data: Vec<NkData>,
    parameters: Option<Vec<Parameter>>,
    lsq: &LmConfig,
) -> Result<(MaterialSpec, FitReport), String> {
    let name = data
        .first()
        .ok_or("fit_dispersion: no nk data")?
        .material
        .clone();
    if let Some(d) = data.iter().find(|d| d.material != name) {
        return Err(format!(
            "fit_dispersion: nk dataset {:?} names material {:?}, the first names {name:?}; \
             one spec, one material",
            d.id, d.material
        ));
    }
    let parameters = match parameters {
        Some(p) => p,
        None => spec_parameters(&name, spec),
    };
    let model =
        FitModel::from_materials(HashMap::from([(name.clone(), Entry::Spec(spec.clone()))]));
    let mut problem = FitProblem::dispersion(model, data, parameters);
    problem.lsq = lsq.clone();
    let report = fit(&problem)?;
    let Some(Entry::Spec(fitted)) = report.model.materials().get(&name) else {
        unreachable!("the model holds the spec it was built with")
    };
    Ok((fitted.clone(), report))
}

/// Every free-by-default number of `spec` (see [`fit_dispersion`]), as
/// parameters of `material[name]`.
pub(crate) fn spec_parameters(name: &str, spec: &MaterialSpec) -> Vec<Parameter> {
    let mixer = MIXERS.contains(&spec.model.as_str());
    let mut paths = Vec::new();
    for (key, value) in &spec.params {
        let skip = match spec.model.as_str() {
            "Table" => key == "n_data" || key == "k_data",
            _ => mixer && value.is_object(),
        };
        if !skip {
            numbers(value, vec![PathSeg::Key(key.clone())], &mut paths);
        }
    }
    paths
        .into_iter()
        .map(|path| {
            Parameter::new(ParamAddr::Material {
                material: name.to_string(),
                path,
            })
        })
        .collect()
}

/// Paths of every number in `value` (lists by index, oscillator dicts by
/// key).
fn numbers(value: &Value, path: Vec<PathSeg>, out: &mut Vec<Vec<PathSeg>>) {
    match value {
        Value::Number(_) => out.push(path),
        Value::Array(items) => {
            for (i, item) in items.iter().enumerate() {
                let mut p = path.clone();
                p.push(PathSeg::Index(i));
                numbers(item, p, out);
            }
        }
        Value::Object(map) => {
            for (k, item) in map {
                let mut p = path.clone();
                p.push(PathSeg::Key(k.clone()));
                numbers(item, p, out);
            }
        }
        _ => {}
    }
}

/// One candidate of [`table_to_model`].
#[derive(Debug, Clone)]
pub struct CandidateFit {
    pub model: String,
    /// The candidate after the fit; `None` when it failed.
    pub spec: Option<MaterialSpec>,
    pub report: Option<FitReport>,
    /// NaN when the fit failed.
    pub aic: f64,
    pub bic: f64,
    /// Why the fit failed.
    pub error: Option<String>,
}

/// Fit each of `candidates` (start specs; every number free, as in
/// [`fit_dispersion`]) to the n (and k, when the table has it) of the
/// `Table` spec `table` on the table's own wavelengths, and rank them by
/// BIC, best first; failed candidates last, in their given order.
///
/// A table carries no uncertainties, so every point weighs the same
/// (σ = 1 for n and k): the AIC/BIC forms estimate the noise scale
/// themselves, so the ranking does not depend on that choice, only on the
/// equal weighting of n and k.
pub fn table_to_model(
    table: &MaterialSpec,
    candidates: Vec<MaterialSpec>,
    lsq: &LmConfig,
) -> Result<Vec<CandidateFit>, String> {
    if table.model != "Table" {
        return Err(format!(
            "table_to_model: a Table spec is needed, got {}",
            table.model
        ));
    }
    if candidates.is_empty() {
        return Err("table_to_model: no candidates".into());
    }
    let grid: Vec<f64> = table
        .params
        .get("n_data")
        .and_then(|v| v.get(0))
        .and_then(Value::as_array)
        .and_then(|w| w.iter().map(Value::as_f64).collect())
        .ok_or("table_to_model: the table has no n_data wavelengths")?;
    let nk = table.evaluate(&grid)?;
    let mut columns = vec![NkColumn::new(
        NkQuantity::N,
        nk.iter().map(|z| z.re).collect(),
        1.0,
    )];
    if table.params.contains_key("k_data") {
        columns.push(NkColumn::new(
            NkQuantity::K,
            nk.iter().map(|z| z.im).collect(),
            1.0,
        ));
    }
    let data = NkData::new("table", "candidate", grid, columns);

    let mut ranked: Vec<CandidateFit> = candidates
        .into_iter()
        .map(|spec| {
            let model = spec.model.clone();
            match fit_dispersion(&spec, vec![data.clone()], None, lsq) {
                Ok((fitted, report)) => CandidateFit {
                    model,
                    spec: Some(fitted),
                    aic: report.aic,
                    bic: report.bic,
                    report: Some(report),
                    error: None,
                },
                Err(e) => CandidateFit {
                    model,
                    spec: None,
                    report: None,
                    aic: f64::NAN,
                    bic: f64::NAN,
                    error: Some(e),
                },
            }
        })
        .collect();
    // Stable: failures (NaN) last, in their given order.
    ranked.sort_by(|a, b| match (a.bic.is_nan(), b.bic.is_nan()) {
        (false, false) => a.bic.total_cmp(&b.bic),
        (x, y) => x.cmp(&y),
    });
    Ok(ranked)
}
