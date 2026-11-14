// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::forward — forward models a fit drives.
//!
//! - [`stack`]: a fit model through `expand` and the general solver.
//! - [`dispersion`]: material specs on measured nk grids.

pub mod dispersion;
pub mod stack;

pub use dispersion::DispersionForward;
pub use stack::StackForward;
