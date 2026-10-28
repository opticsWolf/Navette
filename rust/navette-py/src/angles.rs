// SPDX-License-Identifier: LGPL-3.0-or-later
//! The one place a Python angle becomes a `Radians`.
//!
//! Every door that takes an incidence angle converts here and passes only
//! `Radians` on, so there is one conversion per unit in the whole crate. A
//! value out of [0, 90] deg / [0, pi/2] rad is a `ValueError` naming the
//! value and the unit it was read in.

use navette::smatrix::angle::{self, Radians};
use pyo3::exceptions::PyValueError;
use pyo3::PyResult;

/// An incidence grid in degrees, or in radians with `in_radians`.
pub(crate) fn grid(angles: &[f64], in_radians: bool) -> PyResult<Vec<Radians>> {
    if in_radians {
        angle::radians_slice(angles)
    } else {
        angle::degrees_slice(angles)
    }
    .map_err(PyValueError::new_err)
}

/// One incidence angle in degrees, or in radians with `in_radians`.
pub(crate) fn one(a: f64, in_radians: bool) -> PyResult<Radians> {
    if in_radians {
        Radians::new(a)
    } else {
        Radians::from_degrees(a)
    }
    .map_err(PyValueError::new_err)
}
