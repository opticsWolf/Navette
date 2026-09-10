//! Aggregated PyO3 bindings for Navette — builds the single
//! `navette._navette` extension containing all five native submodules:
//!
//! - `_color` — colorimetry (over `navette::color`)
//! - `_interpolate` — univariate interpolation (over `navette::interpolate`)
//! - `_smatrix` — S-matrix thin-film engine (over `navette::smatrix`)
//! - `_spectralweave` — spectral weaving + targets (over `navette::spectralweave`)
//! - `_materials` — dispersion models (over `navette::materials`)
//!
//! No physics here: every kernel lives in the pure-Rust `navette`
//! umbrella. Wrappers own NumPy inputs, release the GIL while
//! rayon-parallel kernels run, and return NumPy.

mod color;
mod config;
mod interpolate;
mod materials;
mod smatrix;
mod structure;
mod synthesis_merit;
mod synthesis_pipeline;
mod spectralweave;
mod spectralweave_optical;
mod spectralweave_target;

use pyo3::prelude::*;
use pyo3::{wrap_pyfunction, wrap_pymodule};

use crate::color::_color;
use crate::interpolate::_interpolate;
use crate::materials::_materials;
use crate::smatrix::_smatrix;
use crate::spectralweave::_spectralweave;

// Use mimalloc instead of the system allocator. The hot paths here allocate
// many short-lived buffers (and a few large ones); the default Windows heap is
// slow under that pattern.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Cargo profile this extension was compiled with: `"release"` or `"debug"`.
///
/// Exists so a benchmark can refuse to run against an unoptimized build.
/// This is not a hypothetical: the dev venv once shipped a `maturin develop`
/// (debug) build for long enough that a whole round of committed timings —
/// and the conclusions drawn from them — had to be thrown away. A debug
/// build is several times slower and makes every comparison meaningless,
/// while looking perfectly healthy from Python.
///
/// `debug_assertions` is the profile switch Cargo actually flips (`dev` on,
/// `release` off), so this reports how the code was *optimized*, which is
/// what a timing depends on — not the literal `--profile` name.
///
/// See `validation/benches/_bench_common.py::require_release`.
#[pyfunction]
fn build_profile() -> &'static str {
    if cfg!(debug_assertions) {
        "debug"
    } else {
        "release"
    }
}

/// Aggregate all five engine submodules into `navette._navette`.
#[pymodule]
fn _navette(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(_color))?;
    m.add_wrapped(wrap_pymodule!(_interpolate))?;
    m.add_wrapped(wrap_pymodule!(_smatrix))?;
    m.add_wrapped(wrap_pymodule!(structure::_structure))?;
    m.add_wrapped(wrap_pymodule!(_spectralweave))?;
    m.add_wrapped(wrap_pymodule!(_materials))?;
    m.add_function(wrap_pyfunction!(build_profile, m)?)?;
    Ok(())
}
