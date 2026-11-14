// SPDX-License-Identifier: LGPL-3.0-or-later
//! `_fit`: the fit engine's JSON doors (`navette::fit::doors`,
//! `navette::fit::request`) for `navette.fit`.
//!
//! Every door takes JSON text and returns JSON text; the Python package
//! owns the dataclasses on either side. The GIL is released while a fit
//! runs; a progress callback re-acquires it per event, and an exception
//! raised in it stops the fit with that message.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use navette::fit::doors::{
    benchmark_json, conventions_json, fit_dispersion_json, from_solver_json, generate_json,
    table_to_model_json, to_solver_json, tournament_json, truth_json,
};
use navette::fit::request::{evaluate_program, fit_program};
use navette::fit::{Progress, Watcher, DEFAULT_ORDERING_CAP};

fn text(v: Result<serde_json::Value, String>) -> PyResult<String> {
    v.map(|v| v.to_string()).map_err(PyValueError::new_err)
}

fn event<'py>(py: Python<'py>, p: &Progress) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    match p {
        Progress::StageStart { stage, name, free } => {
            d.set_item("event", "stage_start")?;
            d.set_item("stage", stage)?;
            d.set_item("name", name)?;
            d.set_item("free", free.clone())?;
        }
        Progress::Improved { stage, evals, chi2 } => {
            d.set_item("event", "improved")?;
            d.set_item("stage", stage)?;
            d.set_item("evals", evals)?;
            d.set_item("chi2", chi2)?;
        }
        Progress::Searched { stage, evals, chi2 } => {
            d.set_item("event", "searched")?;
            d.set_item("stage", stage)?;
            d.set_item("evals", evals)?;
            d.set_item("chi2", chi2)?;
        }
        Progress::StageEnd {
            stage,
            name,
            chi2,
            accepted,
        } => {
            d.set_item("event", "stage_end")?;
            d.set_item("stage", stage)?;
            d.set_item("name", name)?;
            d.set_item("chi2", chi2)?;
            d.set_item("accepted", accepted)?;
        }
    }
    Ok(d)
}

/// Run the fit a program document's `fit` section describes; returns
/// `{"strategy", "report", "curves"}` as JSON text. `callback(event)` gets
/// a dict per stage start, lower chi2 and stage end; raising in it stops
/// the fit.
#[pyfunction]
#[pyo3(signature = (program_json, callback=None))]
fn fit_program_json(
    py: Python<'_>,
    program_json: &str,
    callback: Option<Py<PyAny>>,
) -> PyResult<String> {
    let program = program_json.to_string();
    let out = py.detach(move || {
        let watcher = callback.map(|cb| {
            Watcher::new(move |p: &Progress| {
                Python::attach(|py| {
                    let d = event(py, p).map_err(|e| e.to_string())?;
                    cb.call1(py, (d,)).map(|_| ()).map_err(|e| e.to_string())
                })
            })
        });
        fit_program(&program, watcher.as_ref())
    });
    text(out)
}

/// The report and curves at the start values, nothing moved.
#[pyfunction]
fn evaluate_program_json(py: Python<'_>, program_json: &str) -> PyResult<String> {
    let program = program_json.to_string();
    text(py.detach(move || evaluate_program(&program)))
}

/// Strategies from a list of strategies and generators.
#[pyfunction]
#[pyo3(signature = (strategies_json, cap=DEFAULT_ORDERING_CAP))]
fn generate_json_text(strategies_json: &str, cap: usize) -> PyResult<String> {
    text(generate_json(strategies_json, cap))
}

/// Every strategy on the program's problem, ranked by chi2.
#[pyfunction]
#[pyo3(signature = (program_json, strategies_json, cap=DEFAULT_ORDERING_CAP))]
fn tournament_json_text(
    py: Python<'_>,
    program_json: &str,
    strategies_json: &str,
    cap: usize,
) -> PyResult<String> {
    let (p, s) = (program_json.to_string(), strategies_json.to_string());
    text(py.detach(move || tournament_json(&p, &s, cap)))
}

/// Every strategy on a truth drawn per seed, ranked by recovery.
#[pyfunction]
#[pyo3(signature = (program_json, strategies_json, seeds, cap=DEFAULT_ORDERING_CAP))]
fn benchmark_json_text(
    py: Python<'_>,
    program_json: &str,
    strategies_json: &str,
    seeds: Vec<u64>,
    cap: usize,
) -> PyResult<String> {
    let (p, s) = (program_json.to_string(), strategies_json.to_string());
    text(py.detach(move || benchmark_json(&p, &s, &seeds, cap)))
}

/// A truth drawn from the program's error laws, with its simulated data.
#[pyfunction]
fn truth_json_text(py: Python<'_>, program_json: &str, seed: u64) -> PyResult<String> {
    let p = program_json.to_string();
    text(py.detach(move || truth_json(&p, seed)))
}

/// An instrument record to solver-convention data items.
#[pyfunction]
fn to_solver_json_text(raw_json: &str, convention: &str) -> PyResult<String> {
    text(to_solver_json(raw_json, convention))
}

/// Solver-convention data items as the instrument writes them.
#[pyfunction]
fn from_solver_json_text(id: &str, data_json: &str, convention: &str) -> PyResult<String> {
    text(from_solver_json(id, data_json, convention))
}

/// The preset instrument conventions, by name.
#[pyfunction]
fn conventions_json_text() -> String {
    conventions_json().to_string()
}

/// One material spec against nk data.
#[pyfunction]
#[pyo3(signature = (spec_json, nk_data_json, parameters_json=None, lsq_json="{}"))]
fn fit_dispersion_json_text(
    py: Python<'_>,
    spec_json: &str,
    nk_data_json: &str,
    parameters_json: Option<&str>,
    lsq_json: &str,
) -> PyResult<String> {
    let (s, d, p, l) = (
        spec_json.to_string(),
        nk_data_json.to_string(),
        parameters_json.map(str::to_string),
        lsq_json.to_string(),
    );
    text(py.detach(move || fit_dispersion_json(&s, &d, p.as_deref(), &l)))
}

/// Model candidates for a Table spec, ranked by BIC.
#[pyfunction]
#[pyo3(signature = (table_json, candidates_json, lsq_json="{}"))]
fn table_to_model_json_text(
    py: Python<'_>,
    table_json: &str,
    candidates_json: &str,
    lsq_json: &str,
) -> PyResult<String> {
    let (t, c, l) = (
        table_json.to_string(),
        candidates_json.to_string(),
        lsq_json.to_string(),
    );
    text(py.detach(move || table_to_model_json(&t, &c, &l)))
}

#[pymodule]
pub fn _fit(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(fit_program_json, m)?)?;
    m.add_function(wrap_pyfunction!(evaluate_program_json, m)?)?;
    m.add_function(wrap_pyfunction!(generate_json_text, m)?)?;
    m.add_function(wrap_pyfunction!(tournament_json_text, m)?)?;
    m.add_function(wrap_pyfunction!(benchmark_json_text, m)?)?;
    m.add_function(wrap_pyfunction!(truth_json_text, m)?)?;
    m.add_function(wrap_pyfunction!(to_solver_json_text, m)?)?;
    m.add_function(wrap_pyfunction!(from_solver_json_text, m)?)?;
    m.add_function(wrap_pyfunction!(conventions_json_text, m)?)?;
    m.add_function(wrap_pyfunction!(fit_dispersion_json_text, m)?)?;
    m.add_function(wrap_pyfunction!(table_to_model_json_text, m)?)?;
    m.add("DEFAULT_ORDERING_CAP", DEFAULT_ORDERING_CAP)?;
    Ok(())
}
