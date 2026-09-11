// SPDX-License-Identifier: LGPL-3.0-or-later
//! Pure-Rust thin-film stack model: layers, groups, expansion, validation.
//!
//! First-class home of the structure model (`Layer`, `Group`, providers,
//! two-phase expansion, `Structure`, `Architect`). The solver
//! (`navette-smatrix`) and the Python bindings (`navette-py`) consume this
//! crate; it depends on neither. Units are nanometres everywhere; the
//! state schema is versioned (see [`version`]).

pub mod architect;
pub mod enums;
pub mod group;
pub mod expansion;
pub mod layer;
pub mod providers;
pub mod specs;
// `structure::structure` (clippy::module_inception): the crate's public path
// is `navette::structure::structure::Structure`, and every Python binding,
// doc example and saved program document names it. Renaming the inner module
// is a path break for no behaviour gain -- R6.6's rename batch is where a
// change like that belongs, if it is worth making at all.
#[allow(clippy::module_inception)]
pub mod structure;
pub mod validation;
pub mod version;
pub mod weaver;

pub use architect::{Architect, Block, SharedStructure};
pub use enums::{BlockKind, ErrorMask, ErrorType, LayerMask, LayerType, OptMask, RoughnessType};

// Shared group handles: structures hold references (like Python dicts holding
// the same Group objects), so film/material baking through one handle is
// visible through all of them. Expansion snapshot-clones per run.
pub type SharedGroup = std::rc::Rc<std::cell::RefCell<Group>>;
pub use group::{ErrorParams, Group};
pub use group::shared_group;
pub use layer::Layer;
pub use providers::{DictProvider, Entry, MaterialProvider, SpecProvider, assert_provider_grid, grids_equal};
pub use expansion::{ExpandOptions, SolverArrays, Span, expand};
pub use structure::Structure;
pub use specs::{MODELS, MaterialSpec};
pub use weaver::{InterpSettings, WeaverProvider, WovenBackend};
pub use validation::{Severity, ValidationIssue};
pub use version::{SCHEMA_VERSION, check_schema_version};
