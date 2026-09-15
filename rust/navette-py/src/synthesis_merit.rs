// SPDX-License-Identifier: LGPL-3.0-or-later
//! Python bindings for the synthesis merit bridge (`navette-smatrix`).
//!
//! Exposes `MeritSpec` (flat optimization targets + residual kernel),
//! `SimCurves` (simulated/derived rows on the solver grid) and the
//! `build_needle_targets` fold. The intended producer is the Python
//! converter (`TargetCollection.build_merit_spec`), which copies finished
//! ingestion products out of `TargetWeaver.export_entries`.

use std::sync::Arc;

use num_complex::Complex64;
use numpy::{PyArray, PyArray1, PyReadonlyArray1};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};

use navette::smatrix::synthesis::merit::{
    ConstraintKind, CurveId, MeritKey, MeritSpec, MeritTarget, SimCurves, SimTransform,
};
use navette::smatrix::synthesis::needle_pass::build_needle_targets as core_fold;

fn parse_curve(s: &str) -> PyResult<CurveId> {
    CurveId::from_str(s).ok_or_else(|| {
        PyValueError::new_err(format!(
            "Invalid curve id '{s}' (use Rs/Rp/Ru/Ts/Tp/Tu, As/Ap/Au, \
             RBs/RBp/RBu/TBs/TBp/TBu, ABs/ABp/ABu)"
        ))
    })
}

fn parse_kind(s: &str) -> PyResult<ConstraintKind> {
    ConstraintKind::from_str(s).ok_or_else(|| {
        PyValueError::new_err(format!(
            "Invalid kind '{s}' (use 'e', 'a', 'b', 'r' or 'c')"
        ))
    })
}

fn parse_transform(s: &str) -> PyResult<SimTransform> {
    SimTransform::from_str(s).ok_or_else(|| {
        PyValueError::new_err(format!(
            "Invalid transform '{s}' (use 'linear', 'log', 'phase' or 'complex')"
        ))
    })
}

/// Extract the scalar-or-array reference index (PD2): a float becomes a
/// length-1 broadcast row, an array passes through for the native length
/// validation (length 1 or `wavelengths.len()`; anything else refuses via
/// `SimCurves::reference_length_issue`). `value` is the raw Python
/// argument (`None` when defaulted).
fn reference_row(value: Option<&Bound<'_, PyAny>>, default: f64) -> PyResult<Arc<[f64]>> {
    match value {
        None => Ok(Arc::from([default])),
        Some(v) => {
            if let Ok(f) = v.extract::<f64>() {
                return Ok(Arc::from([f]));
            }
            let arr: PyReadonlyArray1<f64> = v.extract().map_err(|_| {
                PyTypeError::new_err("reference indices accept a float or a 1-D float array")
            })?;
            Ok(Arc::from(Vec::from(arr.as_slice()?)))
        }
    }
}

/// The scalar-or-array differential index plus its SHAPE: `Some(len)` is 1
/// for a scalar/length-1 row (the warning case), `None` means the caller
/// supplied a full per-wavelength row (trusted - the engine's own shape).
fn reference_row_len(value: Option<&Bound<'_, PyAny>>) -> PyResult<Option<usize>> {
    match value {
        None => Ok(Some(1)),
        Some(v) => {
            if v.extract::<f64>().is_ok() {
                return Ok(Some(1));
            }
            let arr: PyReadonlyArray1<f64> = v.extract().map_err(|_| {
                PyTypeError::new_err("reference indices accept a float or a 1-D float array")
            })?;
            Ok(Some(arr.as_slice()?.len()))
        }
    }
}

#[pyclass(name = "SimCurves")]
/// Simulated rows on the solver grid (see `SimCurves` in the core).
pub struct PySimCurves {
    inner: SimCurves,
}

impl PySimCurves {
    pub(crate) fn wrap(inner: SimCurves) -> Self {
        PySimCurves { inner }
    }
}

#[pymethods]
impl PySimCurves {
    #[new]
    #[pyo3(signature = (angles, wavelengths, total_d=0.0, n_front=None, n_back=None))]
    /// Empty rows on the given axes; fill with `set_curve`/`set_complex`.
    /// `total_d`/`n_front`/`n_back` are the stack metadata for
    /// differential-phase (`PDts`/`PDtp`) demands (defaults zero the
    /// reference, i.e. differential ≡ absolute). `n_front`/`n_back` take a
    /// float (length-1 broadcast row) or a per-wavelength float array
    /// (PD2) - length 1 or `len(wavelengths)`, anything else refuses.
    fn new(
        angles: PyReadonlyArray1<'_, f64>,
        wavelengths: PyReadonlyArray1<'_, f64>,
        total_d: f64,
        n_front: Option<&Bound<'_, PyAny>>,
        n_back: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let a = angles.as_slice()?;
        let w = wavelengths.as_slice()?;
        let nf = reference_row(n_front, 1.0)?;
        let nb = reference_row(n_back, 1.0)?;
        let sim = SimCurves {
            angles: Arc::from(a),
            wavelengths: Arc::from(w),
            total_d,
            n_front_re: nf,
            n_back_re: nb,
            curves: Default::default(),
            back: Default::default(),
            cplx: Default::default(),
            cplx_back: Default::default(),
        };
        if let Some(msg) = sim.reference_length_issue() {
            return Err(PyValueError::new_err(msg));
        }
        Ok(PySimCurves { inner: sim })
    }

    /// Store one intensity row (`Rs/…/ABu`; absorption ids rejected —
    /// absorptance derives from companions).
    fn set_curve(&mut self, curve_id: String, values: PyReadonlyArray1<'_, f64>) -> PyResult<()> {
        // Thin over the core setter (key rules + lengths validated there).
        let id = parse_curve(&curve_id)?;
        let arc: Arc<[f64]> = Arc::from(values.as_slice()?);
        self.inner.set_curve(id, arc).map_err(PyValueError::new_err)
    }

    /// Store one complex-amplitude row for phase demands (`Rs/Rp/Ts/Tp`,
    /// `RBs/RBp/TBs/TBp`; absorption/unpolarized rejected).
    fn set_complex(
        &mut self,
        curve_id: String,
        values: PyReadonlyArray1<'_, Complex64>,
    ) -> PyResult<()> {
        // Thin over the core setter (key rules + lengths validated there).
        let id = parse_curve(&curve_id)?;
        let arc: Arc<[Complex64]> = Arc::from(values.as_slice()?);
        self.inner
            .set_complex(id, arc)
            .map_err(PyValueError::new_err)
    }
}

#[pyclass(name = "MeritSpec")]
/// Flat optimization targets (see `MeritSpec` in the core).
pub struct PyMeritSpec {
    inner: MeritSpec,
}

impl PyMeritSpec {
    pub(crate) fn inner(&self) -> &MeritSpec {
        &self.inner
    }

    pub(crate) fn from_inner(inner: MeritSpec) -> Self {
        PyMeritSpec { inner }
    }
}

/// Compile a target set (JSON) into a native `MeritSpec` (thin over
/// `targets::compile_merit_spec`): the `build_merit_spec` path.
#[pyfunction]
pub(crate) fn compile_merit_spec(request_json: &str) -> PyResult<PyMeritSpec> {
    let set: navette::smatrix::synthesis::targets::TargetSet = serde_json::from_str(request_json)
        .map_err(|e| {
        PyValueError::new_err(format!("compile_merit_spec: invalid request: {e}"))
    })?;
    navette::smatrix::synthesis::targets::compile_merit_spec(&set)
        .map(PyMeritSpec::from_inner)
        .map_err(PyValueError::new_err)
}

#[pymethods]
impl PyMeritSpec {
    #[new]
    fn new() -> Self {
        PyMeritSpec {
            inner: MeritSpec::new(),
        }
    }

    /// Register a `(angle, curve)` demand group; returns its key index.
    fn add_key(&mut self, angle: f64, curve_id: String) -> PyResult<usize> {
        let id = parse_curve(&curve_id)?;
        Ok(self.inner.add_key(MeritKey { angle, curve: id }))
    }

    #[pyo3(signature = (key_idx, wavelengths, normalized, tolerances, kind, transform,
                        norm_factor, band=None, phase=false, differential_passes=None,
                        weight=1.0, count_norm=None, integral=false))]
    /// Append one target frame. Arrays are per-point values on `wavelengths`.
    /// `differential_passes` (None = absolute phase; 1.0 = `PDts`/`PDtp`
    /// transmitted) subtracts the equivalent-medium reference (see
    /// `SimCurves` metadata); requires `phase=true`.
    /// `weight` scales the frame's merit sum; `count_norm` (target-level
    /// point count, resolved by the converter) divides it. `integral`
    /// constrains the mean of the scaled diffs (single residual).
    #[allow(clippy::too_many_arguments)]
    fn add_target(
        &mut self,
        key_idx: u32,
        wavelengths: PyReadonlyArray1<'_, f64>,
        normalized: PyReadonlyArray1<'_, f64>,
        tolerances: PyReadonlyArray1<'_, f64>,
        kind: String,
        transform: String,
        norm_factor: f64,
        band: Option<PyReadonlyArray1<'_, f64>>,
        phase: bool,
        differential_passes: Option<f64>,
        weight: f64,
        count_norm: Option<f64>,
        integral: bool,
    ) -> PyResult<()> {
        let wl = wavelengths.as_slice()?;
        let nt = normalized.as_slice()?;
        let tol = tolerances.as_slice()?;
        let b: Vec<f64> = match band {
            Some(a) => a.as_slice()?.to_vec(),
            None => Vec::new(),
        };
        self.inner
            .add_target(MeritTarget {
                key_idx,
                // F2.1: the Python door is single-environment until
                // F2.4 opens `environments=`. Stated, not defaulted -
                // the field has no `Default`, so a later widening of
                // this door cannot forget to route.
                env_idx: 0,
                wavelengths: Arc::from(wl),
                kind: parse_kind(&kind)?,
                transform: parse_transform(&transform)?,
                norm_factor,
                normalized_targets: Arc::from(nt),
                tolerances: Arc::from(tol),
                band: Arc::from(b.as_slice()),
                phase,
                differential_passes,
                weight,
                count_norm,
                integral,
            })
            .map_err(PyValueError::new_err)
    }

    /// Scalar merit: Σ residual² + `missing_penalty` per missing key group.
    fn merit(&self, py: Python<'_>, sim: &PySimCurves, missing_penalty: f64) -> PyResult<f64> {
        warn_scalar_reference(py, &self.inner, &sim.inner)?;
        let inner = &self.inner;
        let sim_inner = &sim.inner;
        Ok(py.detach(move || inner.merit(sim_inner, missing_penalty)))
    }

    /// Fixed-length residual vector (zeros where inactive).
    fn residuals(
        &self,
        py: Python<'_>,
        sim: &PySimCurves,
    ) -> PyResult<Py<PyArray<f64, numpy::Ix1>>> {
        warn_scalar_reference(py, &self.inner, &sim.inner)?;
        let mut out = Vec::new();
        self.inner
            .residuals(&sim.inner, &mut out)
            .map_err(|id| PyValueError::new_err(format!("missing curve for demand {id:?}")))?;
        Ok(PyArray::from_vec(py, out).unbind())
    }

    /// Total residual components (active or not).
    fn n_residuals(&self) -> usize {
        self.inner.n_residuals()
    }
}

/// Reference-rotation factors for differential-phase demands (thin over
/// the core kernel): `exp(-i·ref)` per wavelength. `n_inc` takes a float
/// (length-1 broadcast) or a per-wavelength float array (PD2); a scalar
/// index on a differential rotation warns (ground rule 6 - announced,
/// never silent; the engine cannot know the medium it was not given),
/// once per call site via Python's warning dedup.
#[pyfunction]
#[pyo3(signature = (wavelengths, angle_deg, n_inc=None, total_d=0.0, passes=1.0))]
pub fn reference_rotation(
    py: Python<'_>,
    wavelengths: PyReadonlyArray1<'_, f64>,
    angle_deg: f64,
    n_inc: Option<&Bound<'_, PyAny>>,
    total_d: f64,
    passes: f64,
) -> PyResult<Py<PyArray1<Complex64>>> {
    use navette::smatrix::synthesis::merit::reference_rotation as core_rot;
    let row = reference_row(n_inc, 1.0)?;
    let n_wl = wavelengths.as_slice()?.len();
    if passes != 0.0 {
        if let Some(row_len) = reference_row_len(n_inc)? {
            if row_len == 1 && row_len != n_wl {
                let n_val = match n_inc {
                    Some(v) if v.extract::<f64>().is_ok() => {
                        format!("{}", v.extract::<f64>().unwrap())
                    }
                    _ => format!("{}", row.as_ref()[0]),
                };
                let msg = format!(
                    "reference_rotation: scalar reference index n_inc={} for a \
differential rotation (passes={}) - if the incidence medium is \
dispersive, pass the per-wavelength array (length {})",
                    n_val, passes, n_wl
                );
                let user_warning = py.get_type::<pyo3::exceptions::PyUserWarning>();
                PyErr::warn(
                    py,
                    &user_warning,
                    std::ffi::CString::new(msg.clone())?.as_c_str(),
                    1,
                )?;
            }
        }
    }
    let rot = core_rot(wavelengths.as_slice()?, angle_deg, &row, total_d, passes)
        .map_err(PyValueError::new_err)?;
    Ok(PyArray::from_vec(py, rot).into())
}

/// Apply per-wavelength rotation factors to flat rows (thin over the
/// core kernel): returns the rotated rows.
#[pyfunction]
pub(crate) fn rotate_rows(rows: Vec<Complex64>, rot: Vec<Complex64>) -> PyResult<Vec<Complex64>> {
    use navette::smatrix::synthesis::merit::rotate_rows as core_rotate;
    let mut out = rows;
    core_rotate(&mut out, &rot).map_err(PyValueError::new_err)?;
    Ok(out)
}

/// The permanent guard (PD2, the part that outlives PD1): a scalar
/// (length-1) reference row on a spec with a differential demand warns.
/// After PD1 the native path is correct by construction (it reads the
/// stack); the remaining way to hand a frozen index to a dispersive
/// medium is a hand-built `SimCurves` with a scalar row. The door cannot
/// know the stack it was not given - whether a constant medium is really
/// air - so this is a warning, never silent and never a refusal (ground
/// rule 6; air stays a legitimate scalar). Python's warning registry
/// dedupes per call site, so an LM loop sees it once.
fn warn_scalar_reference(py: Python<'_>, spec: &MeritSpec, sim: &SimCurves) -> PyResult<()> {
    if !spec.uses_differential() {
        return Ok(());
    }
    let (front_demanded, back_demanded) = spec.demanded_reference_sides();
    let nf = sim.n_front_re.as_ref();
    let nb = sim.n_back_re.as_ref();
    // Report only the sides a demand can actually read (review PD G1):
    // no label maps to the back yet, so a default scalar `n_back` on a
    // front-only spec is correct, not a hazard - warning about it would
    // cry wolf at callers who did exactly what the guard asked.
    if !(front_demanded && nf.len() == 1) && !(back_demanded && nb.len() == 1) {
        return Ok(());
    }
    let sides = [
        ("n_front", front_demanded, nf),
        ("n_back", back_demanded, nb),
    ]
    .into_iter()
    .filter(|(_, demanded, row)| *demanded && row.len() == 1)
    .map(|(name, _, row)| format!("{}={}", name, row[0]))
    .collect::<Vec<_>>()
    .join(", ");
    let msg = format!(
        "differential-phase demand with a scalar reference index ({}) - \
if the incidence or exit medium is dispersive, supply the \
per-wavelength array (SimCurves n_front/n_back or \
apply_reference_rotation)",
        sides
    );
    let user_warning = py.get_type::<pyo3::exceptions::PyUserWarning>();
    PyErr::warn(
        py,
        &user_warning,
        std::ffi::CString::new(msg)?.as_c_str(),
        1,
    )
}

/// Fold a spec into per-quantity `(targets, weights)` pairs (angle-major).
/// Returns a dict with `r/t/a/rb/tb/ab` pairs plus `phi0..phi3` (one pair
/// per S-matrix channel — emit one `P_PHI` call per used channel), plus
/// `grads_r`/`grads_t`: the Option-B color buckets (`dF/dcurve` per
/// angle-major point, weight/residual/U-curve-half already folded in).
/// They are NOT a target/weight pair — pass them straight to
/// `needle_gradient`'s `grads_r`/`grads_t`, which deposits them into the
/// `P` / `P_T` channels. All-zero when the spec has no color demands.
/// Omitting them is what made the documented Python fold →
/// `needle_gradient` flow return a silent-zero color gradient (R4.2).
#[pyfunction]
#[pyo3(signature = (spec, angles, wavelengths, sim=None))]
pub fn build_needle_targets(
    py: Python<'_>,
    spec: &PyMeritSpec,
    angles: PyReadonlyArray1<'_, f64>,
    wavelengths: PyReadonlyArray1<'_, f64>,
    sim: Option<&PySimCurves>,
) -> PyResult<Py<PyDict>> {
    let a = angles.as_slice()?.to_vec();
    let w = wavelengths.as_slice()?.to_vec();
    if let Some(s) = sim {
        warn_scalar_reference(py, &spec.inner, &s.inner)?;
    }
    let nt = py
        .detach({
            let spec_inner = &spec.inner;
            let sim_inner = sim.map(|s| &s.inner);
            move || core_fold(spec_inner, &a, &w, sim_inner)
        })
        .map_err(PyValueError::new_err)?;
    let d = PyDict::new(py);
    let pair = |py: Python<'_>, name: &str, p: (Vec<f64>, Vec<f64>)| -> PyResult<()> {
        let inner = PyDict::new(py);
        inner.set_item("targets", PyArray::from_vec(py, p.0))?;
        inner.set_item("weights", PyArray::from_vec(py, p.1))?;
        d.set_item(name, inner)?;
        Ok(())
    };
    pair(py, "r", nt.r)?;
    pair(py, "t", nt.t)?;
    pair(py, "a", nt.a)?;
    pair(py, "rb", nt.rb)?;
    pair(py, "tb", nt.tb)?;
    pair(py, "ab", nt.ab)?;
    d.set_item("grads_r", PyArray::from_vec(py, nt.grad_r))?;
    d.set_item("grads_t", PyArray::from_vec(py, nt.grad_t))?;
    for (i, p) in nt.phi.into_iter().enumerate() {
        // Exact dM/dD correction for differential demands (0.0 otherwise):
        // subtract from the assembled P_PHI gradient (see `needle_gradient`
        // `gain_shift_phi`). Uniform in z — the needle site never moves.
        let inner = PyDict::new(py);
        inner.set_item("targets", PyArray::from_vec(py, p.0))?;
        inner.set_item("weights", PyArray::from_vec(py, p.1))?;
        inner.set_item("gain_shift", nt.phi_gain_shift[i])?;
        d.set_item(format!("phi{i}"), inner)?;
    }
    Ok(d.unbind())
}
