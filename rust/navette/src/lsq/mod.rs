// SPDX-License-Identifier: LGPL-3.0-or-later
//! Least-squares solvers, shared by every path that fits or refines.
//!
//! - [`lm`]: Navette's bounded Levenberg-Marquardt, its settings
//!   ([`lm::LmConfig`]) and Jacobian sources.
//! - [`trf`]: trust-region reflective.
//! - [`backends`]: which solver runs ([`backends::LsqBackend`]) and
//!   the interval map for backends without native bounds.
//! - [`pool`]: the thread pool a run executes on.
//! - [`de`]: differential evolution over a box (the fit's `global` stage).
//!
//! Nothing here knows what a parameter means: needle synthesis refines
//! thicknesses through it (`smatrix::synthesis::thickness_refine`).

pub mod backends;
#[cfg(feature = "opt-basin")]
mod basin;
pub mod de;
pub mod lm;
pub(crate) mod pool;
pub mod trf;
