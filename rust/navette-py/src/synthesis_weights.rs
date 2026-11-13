// SPDX-License-Identifier: LGPL-3.0-or-later
//! Python door for W4 weight profiles (`navette::smatrix::synthesis::weights`).
//!
//! `WeightProfile` wraps the Rust enum: constructors build it, `evaluate`
//! runs the Rust evaluation, `to_json` hands the target dataclasses the
//! document they forward. No profile maths here.

use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use navette::smatrix::synthesis::weights::WeightProfile;

/// A weight as a function of a target's axis: wavelength in nm for a
/// spectral target, angle in degrees for an angular one.
#[pyclass(
    name = "WeightProfile",
    module = "navette._smatrix",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyWeightProfile {
    pub(crate) inner: WeightProfile,
}

fn checked(p: WeightProfile) -> PyResult<PyWeightProfile> {
    p.validate().map_err(PyValueError::new_err)?;
    Ok(PyWeightProfile { inner: p })
}

#[pymethods]
impl PyWeightProfile {
    /// The same weight everywhere.
    #[staticmethod]
    fn constant(value: f64) -> PyResult<Self> {
        checked(WeightProfile::Constant { value })
    }

    /// Linear from `w_lo` at `lo` to `w_hi` at `hi`, held flat outside.
    #[staticmethod]
    fn ramp(lo: f64, hi: f64, w_lo: f64, w_hi: f64) -> PyResult<Self> {
        checked(WeightProfile::Ramp { lo, hi, w_lo, w_hi })
    }

    /// `base + (peak − base)·exp(−4 ln2 (x − center)²/fwhm²)`.
    #[staticmethod]
    #[pyo3(signature = (center, fwhm, peak, base=0.0))]
    fn gaussian(center: f64, fwhm: f64, peak: f64, base: f64) -> PyResult<Self> {
        checked(WeightProfile::Gaussian {
            center,
            fwhm,
            peak,
            base,
        })
    }

    /// Weight `w` on each closed `(lo, hi, w)` band, `outside` elsewhere.
    /// Bands may touch but not overlap.
    #[staticmethod]
    #[pyo3(signature = (bands, outside=0.0))]
    fn bands(bands: Vec<(f64, f64, f64)>, outside: f64) -> PyResult<Self> {
        checked(WeightProfile::Bands { bands, outside })
    }

    /// Tabulated weights, interpolated (`"pchip"` or `"linear"`) and held
    /// flat past the ends.
    #[staticmethod]
    #[pyo3(signature = (x, w, interp="pchip"))]
    fn table(x: Vec<f64>, w: Vec<f64>, interp: &str) -> PyResult<Self> {
        checked(WeightProfile::Table {
            x,
            w,
            interp: interp.to_string(),
        })
    }

    /// Pointwise product of profiles, e.g. an emphasis times a window.
    #[staticmethod]
    #[pyo3(signature = (*factors))]
    fn product(factors: Vec<PyWeightProfile>) -> PyResult<Self> {
        checked(WeightProfile::Product {
            factors: factors.into_iter().map(|f| f.inner).collect(),
        })
    }

    /// Build from a JSON document (the form `to_json` writes).
    #[staticmethod]
    fn from_json(text: &str) -> PyResult<Self> {
        let p: WeightProfile = serde_json::from_str(text)
            .map_err(|e| PyValueError::new_err(format!("WeightProfile: {e}")))?;
        checked(p)
    }

    /// The weight at each `x`.
    fn evaluate<'py>(
        &self,
        py: Python<'py>,
        x: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let w = self
            .inner
            .evaluate(x.as_slice()?)
            .map_err(PyValueError::new_err)?;
        Ok(PyArray1::from_vec(py, w))
    }

    /// The JSON document a target set carries.
    fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(&self.inner).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("WeightProfile.from_json({:?})", self.to_json()?))
    }
}
