//! Navette — unified optical thin-film engine (pure Rust, no Python).
//!
//! Single crate, six modules. Rust consumers depend on this one crate;
//! the Python bindings in `navette-py` build on it.
//!
//! ```rust,no_run
//! // S-matrix core through the umbrella:
//! let _ = navette::smatrix::core_engine::REQ_RS;
//! ```

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

/// CIE color science: spaces, conversions, delta-E, whites.
pub mod color;
/// Univariate interpolation (linear / pchip / makima).
pub mod interpolate;
/// Optical dispersion models (Cauchy … UBF, tables, EMA, KK).
pub mod materials;
/// S-matrix thin-film engine + needle synthesis pipeline.
pub mod smatrix;
/// Spectral weaving: resampling, merit, targets.
pub mod spectralweave;
/// Stack model: layers, groups, providers, expansion, architect.
pub mod structure;
/// Program documents: versioned envelopes + section assembly.
pub mod config;
