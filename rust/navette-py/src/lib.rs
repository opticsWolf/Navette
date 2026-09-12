// SPDX-License-Identifier: LGPL-3.0-or-later
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

// ---------------------------------------------------------------------------
// Crate-level clippy allowances.
//
// Everything else in this crate is clippy-clean and the CI gate is blocking
// (`-D warnings`), so each of these is a deliberate decision, not a backlog.
// ---------------------------------------------------------------------------

// Every site is a NaN-rejecting validator: `!(x > 0.0)` rejects NaN, whereas
// the `x <= 0.0` that clippy's `partial_cmp` advice leads to waves it through.
// The negation IS the check. Do not "simplify" it.
#![allow(clippy::neg_cmp_op_on_partial_ord)]
// The wide kernel signatures (solver points, needle passes, merit folds) are
// what R6.1 / R6.3 of docs/remediation_plan.md exist to restructure, into
// request/parameter structs. Silenced here so the gate can go blocking now
// rather than after that refactor.
#![allow(clippy::too_many_arguments)]
// The engine addresses caches by flat angle-major arithmetic
// (`k = a * num_wavs + w`, `base = w * n_layers * 2`), so the loop variable is
// an index into several arrays at different strides, not a cursor over one.
// The iterator rewrites clippy proposes are longer and hide the layout.
#![allow(clippy::needless_range_loop)]
// Plan, cache and fold types are tuples-of-vectors by construction and are
// named at exactly one call site each; a `type` alias would add a hop without
// adding meaning.
#![allow(clippy::type_complexity)]
// The inherent `from_str` parsers take a `&str` and return `Result<_, String>`
// with the crate's own error text. Implementing `std::str::FromStr` would make
// them reachable through `"...".parse()` with an error type the rest of the
// crate does not use.
#![allow(clippy::should_implement_trait)]

mod color;
mod config;
mod interpolate;
mod materials;
mod smatrix;
mod spectralweave;
mod spectralweave_optical;
mod spectralweave_target;
mod structure;
mod synthesis_merit;
mod synthesis_pipeline;

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
