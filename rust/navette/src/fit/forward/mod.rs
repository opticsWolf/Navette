// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::forward — forward models a fit drives.
//!
//! - [`stack`]: a fit model through `expand` and the general solver.
//! - `dispersion` (phase 2): material specs on nk data grids.

pub mod stack;

pub use stack::StackForward;
