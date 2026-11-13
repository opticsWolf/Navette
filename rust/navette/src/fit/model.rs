// SPDX-License-Identifier: LGPL-3.0-or-later
//! fit::model — the model one fit evaluation reads and writes.
//!
//! `Architect` and `Structure` share their structures and groups through
//! `Rc<RefCell<..>>`, which a column-parallel Jacobian cannot touch. A
//! [`FitModel`] is the `Send + Sync` mirror the fit owns instead: the
//! architect's traversal (`iter_entries`, one `(Layer, inverted)` per row
//! in chain order, with the structure row each came from), the merged
//! groups, and the material library by value. That is exactly what
//! [`expand`](crate::structure::expand) takes, so an evaluation is "clone,
//! write the parameters, expand" with no handles in between.
//!
//! A layer address names a row of a *structure*, not of the chain: a
//! structure shared by several blocks, or repeated, is one address, and a
//! write reaches every traversal entry it produced, as an edit through the
//! shared handle would.

use std::collections::{HashMap, HashSet};
use std::rc::Rc;

use serde_json::Value;

use crate::config::LoadedProgram;
use crate::structure::architect::{Architect, SharedStructure};
use crate::structure::gradient::{GradientMode, InhMode};
use crate::structure::providers::Entry;
use crate::structure::structure::Structure;
use crate::structure::{Group, Layer};

use super::params::{GroupField, LayerField, LayerRef, ParamAddr, PathSeg};

/// Name and row materials of one structure in the model.
#[derive(Debug, Clone, PartialEq)]
struct StructureRows {
    name: String,
    materials: Vec<String>,
}

/// Everything one evaluation reads; owned, `Send + Sync`.
#[derive(Debug, Clone, PartialEq)]
pub struct FitModel {
    entries: Vec<(Layer, bool)>,
    /// `(structure, row)` each entry came from.
    origin: Vec<(usize, usize)>,
    structures: Vec<StructureRows>,
    groups: HashMap<String, Group>,
    materials: HashMap<String, Entry>,
}

/// An address resolved against one model: where a write lands.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Target {
    Layer {
        entries: Vec<usize>,
        field: LayerField,
    },
    Group {
        group: String,
        field: GroupField,
    },
    /// `params[key]` then the JSON `pointer` below it.
    Material {
        material: String,
        key: String,
        pointer: String,
    },
}

impl FitModel {
    /// Mirror an architect. Each unique structure is named by the label of
    /// the first block that holds it; labels must tell them apart.
    pub fn from_architect(
        arch: &Architect,
        materials: HashMap<String, Entry>,
    ) -> Result<Self, String> {
        if arch.is_empty() {
            return Err("fit model: the architect is empty".into());
        }
        let issues = arch.validate_chain();
        if let Some(e) = issues.iter().find(|i| i.is_error()) {
            return Err(format!("fit model: {}", e.message));
        }
        let groups = arch.merged_groups()?;

        let mut ids: Vec<*const std::cell::RefCell<Structure>> = Vec::new();
        let mut structures: Vec<StructureRows> = Vec::new();
        let mut entries = Vec::new();
        let mut origin = Vec::new();
        for (b, block) in arch.blocks.iter().enumerate() {
            let key = Rc::as_ptr(&block.structure);
            let sid = match ids.iter().position(|k| *k == key) {
                Some(sid) => sid,
                None => {
                    let name = if block.label.is_empty() {
                        format!("block{b}")
                    } else {
                        block.label.clone()
                    };
                    if structures.iter().any(|s| s.name == name) {
                        return Err(format!(
                            "fit model: two different structures are both named {name:?}; \
                             give their blocks distinct labels"
                        ));
                    }
                    ids.push(key);
                    structures.push(rows_of(&name, &block.structure));
                    ids.len() - 1
                }
            };
            // Same order as `Architect::iter_entries`.
            let s = block.structure.borrow();
            let n = s.layers.len();
            for _ in 0..block.repeat_count {
                let rows: Box<dyn Iterator<Item = usize>> = if block.inverted {
                    Box::new((0..n).rev())
                } else {
                    Box::new(0..n)
                };
                for i in rows {
                    entries.push((s.layers[i].clone(), block.inverted));
                    origin.push((sid, i));
                }
            }
        }
        Ok(Self {
            entries,
            origin,
            structures,
            groups,
            materials,
        })
    }

    /// Mirror one structure, named `name` in addresses.
    pub fn from_structure(
        structure: &Structure,
        name: &str,
        materials: HashMap<String, Entry>,
    ) -> Result<Self, String> {
        let groups = structure
            .groups
            .iter()
            .map(|(k, g)| (k.clone(), g.borrow().clone()))
            .collect();
        Ok(Self {
            entries: structure
                .layers
                .iter()
                .map(|l| (l.clone(), false))
                .collect(),
            origin: (0..structure.layers.len()).map(|i| (0, i)).collect(),
            structures: vec![StructureRows {
                name: name.to_string(),
                materials: structure
                    .layers
                    .iter()
                    .map(|l| l.material.clone())
                    .collect(),
            }],
            groups,
            materials,
        })
    }

    /// The architect of a loaded program (or its one structure) with the
    /// program's material library.
    pub fn from_program(program: &LoadedProgram) -> Result<Self, String> {
        let materials = program
            .materials
            .as_ref()
            .ok_or("fit model: the program has no materials section")?
            .entries()
            .clone();
        match &program.architect {
            Some(arch) => Self::from_architect(arch, materials),
            None if program.structures.len() == 1 => {
                let (name, s) = program.structures.iter().next().unwrap();
                Self::from_structure(&s.borrow(), name, materials)
            }
            None => {
                Err("fit model: the program has no architect and not exactly one structure".into())
            }
        }
    }

    /// Materials only, no stack: a dispersion fit against nk data.
    pub fn from_materials(materials: HashMap<String, Entry>) -> Self {
        Self {
            entries: Vec::new(),
            origin: Vec::new(),
            structures: Vec::new(),
            groups: HashMap::new(),
            materials,
        }
    }

    /// The chain as `expand` takes it.
    pub fn entries(&self) -> &[(Layer, bool)] {
        &self.entries
    }

    pub fn groups(&self) -> &HashMap<String, Group> {
        &self.groups
    }

    pub fn materials(&self) -> &HashMap<String, Entry> {
        &self.materials
    }

    pub fn structure_names(&self) -> Vec<&str> {
        self.structures.iter().map(|s| s.name.as_str()).collect()
    }

    /// Read the scalar at `addr`.
    pub fn get(&self, addr: &ParamAddr) -> Result<f64, String> {
        self.read(&self.resolve(addr)?)
            .map_err(|e| format!("{addr}: {e}"))
    }

    /// Write the scalar at `addr`.
    pub fn set(&mut self, addr: &ParamAddr, value: f64) -> Result<(), String> {
        let target = self.resolve(addr)?;
        self.write(&target, value)
            .map_err(|e| format!("{addr}: {e}"))
    }

    /// Validation of every structure row and group, as the architect's own
    /// gate runs it; the first error refuses.
    pub fn check(&self) -> Result<(), String> {
        let mut seen = HashSet::new();
        for (e, (layer, _)) in self.entries.iter().enumerate() {
            let (sid, row) = self.origin[e];
            if !seen.insert((sid, row)) {
                continue;
            }
            let label = format!("layer[{}:{row}]", self.structures[sid].name);
            if let Some(issue) = layer.property_issues(&label).iter().find(|i| i.is_error()) {
                return Err(issue.message.clone());
            }
        }
        for (name, group) in &self.groups {
            if let Some(issue) = group.validate().iter().find(|i| i.is_error()) {
                return Err(format!("group[{name}]: {}", issue.message));
            }
        }
        Ok(())
    }

    pub(crate) fn resolve(&self, addr: &ParamAddr) -> Result<Target, String> {
        let target = match addr {
            ParamAddr::Layer {
                structure,
                layer,
                field,
            } => {
                let sid = self
                    .structures
                    .iter()
                    .position(|s| s.name == *structure)
                    .ok_or_else(|| {
                        format!(
                            "{addr}: no structure {structure:?}; the model has {:?}",
                            self.structure_names()
                        )
                    })?;
                let rows = &self.structures[sid].materials;
                let row = match layer {
                    LayerRef::Index(i) if *i < rows.len() => *i,
                    LayerRef::Index(i) => {
                        return Err(format!(
                            "{addr}: row {i} is outside {structure:?} ({} rows)",
                            rows.len()
                        ));
                    }
                    LayerRef::Material { material, nth } => rows
                        .iter()
                        .enumerate()
                        .filter(|(_, m)| *m == material)
                        .nth(nth - 1)
                        .map(|(i, _)| i)
                        .ok_or_else(|| {
                            format!("{addr}: {structure:?} has fewer than {nth} {material} rows")
                        })?,
                };
                let entries: Vec<usize> = (0..self.origin.len())
                    .filter(|e| self.origin[*e] == (sid, row))
                    .collect();
                if entries.is_empty() {
                    return Err(format!("{addr}: the row is not in the chain"));
                }
                Target::Layer {
                    entries,
                    field: *field,
                }
            }
            ParamAddr::Group { group, field } => {
                if !self.groups.contains_key(group) {
                    let mut known: Vec<&String> = self.groups.keys().collect();
                    known.sort();
                    return Err(format!(
                        "{addr}: no group {group:?}; the model has {known:?}"
                    ));
                }
                Target::Group {
                    group: group.clone(),
                    field: *field,
                }
            }
            ParamAddr::Material { material, path } => {
                let spec = match self.materials.get(material) {
                    Some(Entry::Spec(spec)) => spec,
                    Some(Entry::Array(_)) => {
                        return Err(format!(
                            "{addr}: {material} is an nk array; only a model spec has parameters"
                        ));
                    }
                    None => return Err(format!("{addr}: no material {material:?}")),
                };
                let PathSeg::Key(key) = &path[0] else {
                    unreachable!("ParamAddr::parse puts a key first")
                };
                let mut cur = spec.params.get(key).ok_or_else(|| {
                    format!(
                        "{addr}: {} has no parameter {key:?}; it has {:?}",
                        spec.model,
                        spec.params.keys().collect::<Vec<_>>()
                    )
                })?;
                let mut pointer = String::new();
                for seg in &path[1..] {
                    match seg {
                        PathSeg::Index(i) => {
                            cur = cur.get(*i).ok_or_else(|| {
                                format!("{addr}: index {i} is outside the list at {key}{pointer}")
                            })?;
                            pointer.push_str(&format!("/{i}"));
                        }
                        PathSeg::Key(k) => {
                            // A nested spec is entered by its own keys.
                            if let Some(next) = cur.get(k) {
                                cur = next;
                                pointer.push_str(&format!("/{k}"));
                            } else if let Some(next) = cur.get("params").and_then(|p| p.get(k)) {
                                cur = next;
                                pointer.push_str(&format!("/params/{k}"));
                            } else {
                                return Err(format!(
                                    "{addr}: nothing named {k:?} at {key}{pointer}"
                                ));
                            }
                        }
                    }
                }
                if !cur.is_number() {
                    return Err(format!("{addr}: lands on {cur}, not on a number"));
                }
                Target::Material {
                    material: material.clone(),
                    key: key.clone(),
                    pointer,
                }
            }
        };
        self.read(&target).map_err(|e| format!("{addr}: {e}"))?;
        Ok(target)
    }

    pub(crate) fn read(&self, target: &Target) -> Result<f64, String> {
        match target {
            Target::Layer { entries, field } => {
                let layer = &self.entries[entries[0]].0;
                read_layer(layer, *field)
            }
            Target::Group { group, field } => Ok(*group_slot_ref(&self.groups[group], *field)),
            Target::Material {
                material,
                key,
                pointer,
            } => {
                let Some(Entry::Spec(spec)) = self.materials.get(material) else {
                    unreachable!("resolved against this model")
                };
                spec.params[key]
                    .pointer(pointer)
                    .and_then(Value::as_f64)
                    .ok_or_else(|| "not a number".to_string())
            }
        }
    }

    pub(crate) fn write(&mut self, target: &Target, value: f64) -> Result<(), String> {
        if !value.is_finite() {
            return Err(format!("non-finite value {value}"));
        }
        match target {
            Target::Layer { entries, field } => {
                for &e in entries {
                    write_layer(&mut self.entries[e].0, *field, value)?;
                }
            }
            Target::Group { group, field } => {
                let g = self
                    .groups
                    .get_mut(group)
                    .expect("resolved against this model");
                *group_slot(g, *field) = value;
            }
            Target::Material {
                material,
                key,
                pointer,
            } => {
                let Some(Entry::Spec(spec)) = self.materials.get_mut(material) else {
                    unreachable!("resolved against this model")
                };
                let slot = spec
                    .params
                    .get_mut(key)
                    .and_then(|v| v.pointer_mut(pointer))
                    .ok_or("the address no longer resolves")?;
                *slot = Value::from(value);
            }
        }
        Ok(())
    }

    /// For a material target: the model of the innermost spec holding the
    /// scalar, and the last key on the way to it (`fraction`, `C1`, `osc`).
    pub(crate) fn material_slot(&self, target: &Target) -> Option<(String, String)> {
        let Target::Material {
            material,
            key,
            pointer,
        } = target
        else {
            return None;
        };
        let Some(Entry::Spec(spec)) = self.materials.get(material) else {
            return None;
        };
        let mut model = spec.model.clone();
        let mut leaf = key.clone();
        let mut cur = &spec.params[key];
        for seg in pointer.split('/').skip(1) {
            match seg.parse::<usize>() {
                Ok(i) => cur = &cur[i],
                Err(_) => {
                    if seg == "params" {
                        // Entering a nested spec: its own model rules below.
                        if let Some(m) = cur.get("model").and_then(Value::as_str) {
                            model = m.to_string();
                        }
                    } else {
                        leaf = seg.to_string();
                    }
                    cur = &cur[seg];
                }
            }
        }
        Some((model, leaf))
    }
}

fn rows_of(name: &str, s: &SharedStructure) -> StructureRows {
    StructureRows {
        name: name.to_string(),
        materials: s
            .borrow()
            .layers
            .iter()
            .map(|l| l.material.clone())
            .collect(),
    }
}

fn read_layer(layer: &Layer, field: LayerField) -> Result<f64, String> {
    Ok(match field {
        LayerField::Thickness => layer.thickness,
        LayerField::InhDelta => layer.inh_delta,
        LayerField::Roughness => layer.roughness,
        LayerField::InterfaceThickness => layer.interface_thickness,
        LayerField::InhRate | LayerField::InhCap => match layer.inh_mode {
            InhMode::RateCapped { rate, cap, .. } => {
                if field == LayerField::InhRate {
                    rate
                } else {
                    cap
                }
            }
            InhMode::Fixed => return Err("the layer's grading is not rate-capped".into()),
        },
        _ => {
            let mode = &layer
                .gradient
                .as_ref()
                .ok_or("the layer has no gradient")?
                .mode;
            match (mode, field) {
                (GradientMode::FixedSpan { f_start, .. }, LayerField::GradientFStart)
                | (GradientMode::RateCapped { f_start, .. }, LayerField::GradientFStart) => {
                    *f_start
                }
                (GradientMode::FixedSpan { f_end, .. }, LayerField::GradientFEnd) => *f_end,
                (GradientMode::RateCapped { rate, .. }, LayerField::GradientRate) => *rate,
                (GradientMode::RateCapped { f_min, .. }, LayerField::GradientFMin) => *f_min,
                (GradientMode::RateCapped { f_max, .. }, LayerField::GradientFMax) => *f_max,
                _ => return Err("the layer's gradient mode has no such field".into()),
            }
        }
    })
}

fn write_layer(layer: &mut Layer, field: LayerField, value: f64) -> Result<(), String> {
    let slot: &mut f64 = match field {
        LayerField::Thickness => &mut layer.thickness,
        LayerField::InhDelta => &mut layer.inh_delta,
        LayerField::Roughness => &mut layer.roughness,
        LayerField::InterfaceThickness => &mut layer.interface_thickness,
        LayerField::InhRate | LayerField::InhCap => match &mut layer.inh_mode {
            InhMode::RateCapped { rate, cap, .. } => {
                if field == LayerField::InhRate {
                    rate
                } else {
                    cap
                }
            }
            InhMode::Fixed => return Err("the layer's grading is not rate-capped".into()),
        },
        _ => {
            let mode = &mut layer
                .gradient
                .as_mut()
                .ok_or("the layer has no gradient")?
                .mode;
            match (mode, field) {
                (GradientMode::FixedSpan { f_start, .. }, LayerField::GradientFStart)
                | (GradientMode::RateCapped { f_start, .. }, LayerField::GradientFStart) => f_start,
                (GradientMode::FixedSpan { f_end, .. }, LayerField::GradientFEnd) => f_end,
                (GradientMode::RateCapped { rate, .. }, LayerField::GradientRate) => rate,
                (GradientMode::RateCapped { f_min, .. }, LayerField::GradientFMin) => f_min,
                (GradientMode::RateCapped { f_max, .. }, LayerField::GradientFMax) => f_max,
                _ => return Err("the layer's gradient mode has no such field".into()),
            }
        }
    };
    *slot = value;
    Ok(())
}

fn group_slot(g: &mut Group, field: GroupField) -> &mut f64 {
    match field {
        GroupField::NFactor => &mut g.n_factor,
        GroupField::KFactor => &mut g.k_factor,
        GroupField::ThickFactor => &mut g.thick_factor,
        GroupField::ThickSummand => &mut g.thick_summand,
        GroupField::InhDeltaSummand => &mut g.inh_delta_summand,
        GroupField::RoughnessSummand => &mut g.roughness_summand,
        GroupField::InterfaceSummand => &mut g.interface_summand,
    }
}

fn group_slot_ref(g: &Group, field: GroupField) -> &f64 {
    match field {
        GroupField::NFactor => &g.n_factor,
        GroupField::KFactor => &g.k_factor,
        GroupField::ThickFactor => &g.thick_factor,
        GroupField::ThickSummand => &g.thick_summand,
        GroupField::InhDeltaSummand => &g.inh_delta_summand,
        GroupField::RoughnessSummand => &g.roughness_summand,
        GroupField::InterfaceSummand => &g.interface_summand,
    }
}
