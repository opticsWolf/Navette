// SPDX-License-Identifier: LGPL-3.0-or-later
//! The one place a Python angle becomes a `Radians`.
//!
//! Every door that takes an incidence angle has two keyword-only arguments,
//! one per unit (`angles_deg=` / `angles_rad=`, or `angle_deg=` /
//! `angle_rad=`). Exactly one must be given; both or neither is a
//! `ValueError` naming the pair. The conversion happens here and only
//! `Radians` goes on, so there is one conversion per unit in the whole
//! crate. A value out of [0, 90] deg / [0, pi/2] rad is a `ValueError`
//! naming the value and the unit it was read in.

use navette::smatrix::angle::{self, Radians};
use numpy::PyReadonlyArray1;
use pyo3::exceptions::PyValueError;
use pyo3::{PyErr, PyResult};

fn exactly_one(what: &str, both: bool) -> PyErr {
    PyValueError::new_err(if both {
        format!("give exactly one of {what}_deg= or {what}_rad=, not both")
    } else {
        format!("give exactly one of {what}_deg= or {what}_rad=")
    })
}

/// An incidence grid from the `{what}_deg=` / `{what}_rad=` pair.
pub(crate) fn pick_grid(
    what: &str,
    deg: Option<PyReadonlyArray1<'_, f64>>,
    rad: Option<PyReadonlyArray1<'_, f64>>,
) -> PyResult<Vec<Radians>> {
    match (deg, rad) {
        (Some(d), None) => {
            angle::degrees_slice(d.as_slice()?).map_err(|e| format!("`{what}_deg` {e}"))
        }
        (None, Some(r)) => {
            angle::radians_slice(r.as_slice()?).map_err(|e| format!("`{what}_rad` {e}"))
        }
        (d, _) => return Err(exactly_one(what, d.is_some())),
    }
    .map_err(PyValueError::new_err)
}

/// One incidence angle from the `{what}_deg=` / `{what}_rad=` pair.
pub(crate) fn pick_one(what: &str, deg: Option<f64>, rad: Option<f64>) -> PyResult<Radians> {
    match (deg, rad) {
        (Some(d), None) => Radians::from_degrees(d),
        (None, Some(r)) => Radians::new(r),
        (d, _) => return Err(exactly_one(what, d.is_some())),
    }
    .map_err(PyValueError::new_err)
}
