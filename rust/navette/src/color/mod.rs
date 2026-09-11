// SPDX-License-Identifier: LGPL-3.0-or-later
//! # navette_color
//!
//! Rust rewrite of the Loom unified color engine — a parity port of
//! `navette_colorengine.py`. Color batches are `&[[f64; 3]]` (stride-of-3 for
//! aggressive auto-vectorization); spectra are `&[f64]`.
//!
//! Module map:
//! - [`common`]      core constants, `mat3_mul_vec`, transfer fns, sRGB/Lab core
//! - [`matrices`]    generated matrix constants (natural form + numpy inverses)
//! - [`composites`]  chained convenience pipelines (e.g. sRGB <-> Lab) and gamut
//! - [`xyy`]           XYZ <-> xyY
//! - [`lch`]           Lab <-> LCh
//! - [`luv`]           XYZ <-> CIELUV
//! - [`oklab_xyz`]     XYZ <-> Oklab (direct XYZ matrices)
//! - [`oklab_srgb`]    sRGB <-> Oklab (legacy sRGB matrices)
//! - [`uvw1964`]       CIE 1964 U*V*W*
//! - [`ucs1960`]       CIE 1960 UCS & chromaticity
//! - [`bradford`]      Bradford chromatic adaptation
//! - [`delta_e_76`]    Delta E 76
//! - [`delta_e_94`]    Delta E 94
//! - [`delta_e_cmc`]   Delta E CMC(l:c)
//! - [`din99`]         DIN99
//! - [`spectral_srgb`] spectral pipeline (SPD x CMF -> sRGB)
//! - [`photometry`]    photometry engine
//! - [`shapes`]        shape handling & broadcasting
//! - [`delta_e_2000`]  CIEDE2000
//!
//! ## Migration record (R6.6, 0.6.22)
//!
//! These sixteen were named `func_01`..`func_16` after the port-task numbering
//! until 0.6.22. The port records under `docs/plans/color/func_NN_*.md` keep
//! their original names -- they are dated documents, not live docs -- so this
//! table is how a `func_NN` reference in a plan, a commit message or the code
//! review is resolved:
//!
//! | was | is | | was | is |
//! |---|---|---|---|---|
//! | `func_01` | [`xyy`]        | | `func_09` | [`delta_e_76`]    |
//! | `func_02` | [`lch`]        | | `func_10` | [`delta_e_94`]    |
//! | `func_03` | [`luv`]        | | `func_11` | [`delta_e_cmc`]   |
//! | `func_04` | [`oklab_xyz`]  | | `func_12` | [`din99`]         |
//! | `func_05` | [`oklab_srgb`] | | `func_13` | [`spectral_srgb`] |
//! | `func_06` | [`uvw1964`]    | | `func_14` | [`photometry`]    |
//! | `func_07` | [`ucs1960`]    | | `func_15` | [`shapes`]        |
//! | `func_08` | [`bradford`]   | | `func_16` | [`delta_e_2000`]  |
//!
//! The single-digit `func_0`..`func_5` that still appear in `smatrix/` are a
//! different numbering entirely -- the numba reference implementation's module
//! names, which the parity tests import by name. They are not covered here.
//!
//! The Python bindings live in the `navette-py` aggregator crate and expose
//! these kernels as the `navette._color` extension submodule (one of five in
//! the single `navette._navette` module built with `maturin develop` from the
//! workspace root). This crate itself is pure Rust: no pyo3, no I/O.

pub mod common;
pub mod matrices;
pub mod metrics;
pub mod composites;
pub mod tables;

pub mod xyy;
pub mod lch;
pub mod luv;
pub mod oklab_xyz;
pub mod oklab_srgb;
pub mod uvw1964;
pub mod ucs1960;
pub mod bradford;
pub mod delta_e_76;
pub mod delta_e_94;
pub mod delta_e_cmc;
pub mod din99;
pub mod spectral_srgb;
pub mod photometry;
pub mod shapes;
pub mod delta_e_2000;

/// Golden-vector parity suite (reference-engine vectors in `golden.rs`).
/// Wired as a test-only module so `cargo test -p navette-color` runs it;
/// it was previously orphaned (never compiled), which hid the D50 drift.
#[cfg(test)]
#[path = "parity.rs"]
mod parity;

pub mod prelude {
    pub use crate::color::common::{REF_WHITE_D50, REF_WHITE_D65};
    pub use crate::color::composites::*;
    pub use crate::color::xyy::{xyy_to_xyz, xyz_to_xyy};
    pub use crate::color::lch::{lab_to_lch, lch_to_lab};
    pub use crate::color::luv::{luv_to_xyz, xyz_to_luv};
    pub use crate::color::oklab_xyz::{oklab_to_xyz, xyz_to_oklab};
    pub use crate::color::oklab_srgb::{oklab_to_srgb, srgb_to_oklab};
    pub use crate::color::uvw1964::{uvw_to_xyz, white_point_uv1960, xyz_to_uvw};
    pub use crate::color::ucs1960::{ucs_to_xyz, uv1960_to_xy, uv1976_to_xy, xyz_to_ucs, xyz_to_ucs_uv};
    pub use crate::color::bradford::{adapt, calc_transform_matrix};
    pub use crate::color::delta_e_76::delta_e_76;
    pub use crate::color::delta_e_94::{delta_e_94, De94Params};
    pub use crate::color::delta_e_cmc::delta_e_cmc;
    pub use crate::color::din99::delta_e_din99;
    pub use crate::color::spectral_srgb::spectral_to_srgb;
    pub use crate::color::photometry::{PhotometryEngine, Vision};
    pub use crate::color::delta_e_2000::delta_e_2000;
}

// ============================================================================
// Python bindings (numpy-based). Enabled with `--features python`.
// Targets pyo3 0.28 + rust-numpy 0.28.
// ============================================================================
