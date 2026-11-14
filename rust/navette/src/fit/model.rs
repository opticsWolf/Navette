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
//!
//! **Environments** (several samples, one model, plan §3.6): the model
//! holds one or more named chains. The architect's own chain is
//! [`DEFAULT_ENV`]; [`FitModel::add_environment`] adds others built from
//! the same structures and fixed surroundings (the bare substrate, a
//! witness film on silicon). Structures, groups and materials are shared
//! by construction, so a parameter is one value in every sample. Fixed
//! rows are surroundings, not addresses; their materials and groups are
//! still the model's.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::ops::Range;
use std::rc::Rc;

use serde_json::Value;

use crate::config::LoadedProgram;
use crate::structure::architect::{Architect, SharedStructure};
use crate::structure::enums::LayerType;
use crate::structure::gradient::{GradientMode, InhMode};
use crate::structure::providers::Entry;
use crate::structure::structure::Structure;
use crate::structure::{Group, Layer};

use super::data::{Instrument, SpectrumData};
use super::params::{GroupField, InstrumentField, LayerField, LayerRef, ParamAddr, PathSeg};

/// Name and row materials of one structure in the model.
#[derive(Debug, Clone, PartialEq)]
struct StructureRows {
    name: String,
    materials: Vec<String>,
}

/// The environment an architect's or structure's own chain is named.
pub const DEFAULT_ENV: &str = "default";

/// `(model, param)` lists that are a free-form material's values, one fit
/// parameter per knot (`material[X].n[3]`); smoothness rows run over them.
pub const FREEFORM_VALUES: [(&str, &str); 5] = [
    ("PointWise", "n"),
    ("PointWise", "k"),
    ("Spline", "n"),
    ("Spline", "k"),
    ("BSplineEps", "eps2"),
];

/// One piece of an environment's chain.
#[derive(Debug, Clone, PartialEq)]
pub enum Segment {
    /// Every row of a model structure, `repeat` times, reversed when
    /// `inverted`; its rows stay addressable.
    Structure {
        name: String,
        inverted: bool,
        repeat: usize,
    },
    /// One fixed row (ambient, substrate, a cover): held, not addressable.
    Fixed(Layer),
}

impl Segment {
    /// A structure once, not inverted.
    pub fn structure(name: impl Into<String>) -> Self {
        Segment::Structure {
            name: name.into(),
            inverted: false,
            repeat: 1,
        }
    }
}

/// Everything one evaluation reads; owned, `Send + Sync`.
#[derive(Debug, Clone, PartialEq)]
pub struct FitModel {
    /// Every environment's chain, back to back.
    entries: Vec<(Layer, bool)>,
    /// `(structure, row)` each entry came from; `None` for a fixed row.
    origin: Vec<Option<(usize, usize)>>,
    /// Name and range in `entries` of each environment.
    envs: Vec<(String, Range<usize>)>,
    structures: Vec<StructureRows>,
    groups: HashMap<String, Group>,
    materials: HashMap<String, Entry>,
    /// Each dataset's instrument by dataset id, seeded from the data when
    /// a problem compiles; a fitted model keeps its fitted values.
    instruments: BTreeMap<String, Instrument>,
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
    Data {
        dataset: String,
        field: InstrumentField,
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
                    origin.push(Some((sid, i)));
                }
            }
        }
        Ok(Self {
            envs: vec![(DEFAULT_ENV.to_string(), 0..entries.len())],
            entries,
            origin,
            structures,
            groups,
            materials,
            instruments: BTreeMap::new(),
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
            origin: (0..structure.layers.len()).map(|i| Some((0, i))).collect(),
            envs: vec![(DEFAULT_ENV.to_string(), 0..structure.layers.len())],
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
            instruments: BTreeMap::new(),
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
            envs: Vec::new(),
            structures: Vec::new(),
            groups: HashMap::new(),
            materials,
            instruments: BTreeMap::new(),
        }
    }

    /// Add a sample configuration: its chain from the model's structures
    /// and fixed rows, ambient first and substrate last.
    pub fn add_environment(&mut self, name: &str, segments: Vec<Segment>) -> Result<(), String> {
        let bad = |why: String| format!("environment {name:?}: {why}");
        if self.envs.iter().any(|(n, _)| n == name) {
            return Err(bad("the model already has it".into()));
        }
        let mut entries = Vec::new();
        let mut origin = Vec::new();
        for seg in segments {
            match seg {
                Segment::Fixed(layer) => {
                    let label = format!("{name}:fixed {}", layer.material);
                    if let Some(i) = layer.property_issues(&label).iter().find(|i| i.is_error()) {
                        return Err(bad(i.message.clone()));
                    }
                    entries.push((layer, false));
                    origin.push(None);
                }
                Segment::Structure {
                    name: s,
                    inverted,
                    repeat,
                } => {
                    let sid = self
                        .structures
                        .iter()
                        .position(|r| r.name == s)
                        .ok_or_else(|| {
                            bad(format!(
                                "no structure {s:?}; the model has {:?}",
                                self.structure_names()
                            ))
                        })?;
                    let n = self.structures[sid].materials.len();
                    // A row's current layer: the first entry it produced.
                    let layer_of = |row: usize| {
                        let e = self
                            .origin
                            .iter()
                            .position(|o| *o == Some((sid, row)))
                            .expect("every structure row is in the default chain");
                        self.entries[e].0.clone()
                    };
                    for _ in 0..repeat {
                        let rows: Box<dyn Iterator<Item = usize>> = if inverted {
                            Box::new((0..n).rev())
                        } else {
                            Box::new(0..n)
                        };
                        for row in rows {
                            entries.push((layer_of(row), inverted));
                            origin.push(Some((sid, row)));
                        }
                    }
                }
            }
        }
        let roles: Vec<LayerType> = entries.iter().map(|(l, _)| l.layer_type).collect();
        if roles.len() < 2
            || roles[0] != LayerType::Ambient
            || roles[roles.len() - 1] != LayerType::Substrate
            || roles[1..roles.len() - 1]
                .iter()
                .any(|r| *r != LayerType::Film)
        {
            return Err(bad(
                "the chain must run ambient, films, substrate (a bare sample is ambient, substrate)"
                    .into(),
            ));
        }
        let start = self.entries.len();
        self.entries.extend(entries);
        self.origin.extend(origin);
        self.envs
            .push((name.to_string(), start..self.entries.len()));
        Ok(())
    }

    /// Environment names, the default first.
    pub fn environment_names(&self) -> Vec<&str> {
        self.envs.iter().map(|(n, _)| n.as_str()).collect()
    }

    /// One environment's chain as `expand` takes it.
    pub fn entries(&self, env: &str) -> Result<&[(Layer, bool)], String> {
        self.envs
            .iter()
            .find(|(n, _)| n == env)
            .map(|(_, r)| &self.entries[r.clone()])
            .ok_or_else(|| {
                format!(
                    "no environment {env:?}; the model has {:?}",
                    self.environment_names()
                )
            })
    }

    /// Where environment `env`'s chain sits among every environment's
    /// entries (the indices a layer [`Target`] names).
    pub(crate) fn env_range(&self, env: &str) -> Option<Range<usize>> {
        self.envs
            .iter()
            .find(|(n, _)| n == env)
            .map(|(_, r)| r.clone())
    }

    pub fn groups(&self) -> &HashMap<String, Group> {
        &self.groups
    }

    pub fn materials(&self) -> &HashMap<String, Entry> {
        &self.materials
    }

    /// The instrument of dataset `id`, once a problem has seeded it.
    pub fn instrument(&self, id: &str) -> Option<&Instrument> {
        self.instruments.get(id)
    }

    /// Seed each dataset's instrument from `data`, keeping any the model
    /// already holds (a fitted model passed on to the next stage).
    pub(crate) fn seed_instruments(&mut self, data: &[SpectrumData]) {
        for d in data {
            self.instruments.entry(d.id.clone()).or_insert(d.instrument);
        }
    }

    /// Drop identity instruments that `data` would seed as the identity
    /// too: an absent entry reads as the identity. One whose dataset
    /// starts elsewhere stays, or a later seeding would restore the start.
    pub(crate) fn prune_instruments(&mut self, data: &[SpectrumData]) {
        self.instruments.retain(|id, i| {
            !i.is_identity()
                || data
                    .iter()
                    .any(|d| d.id == *id && !d.instrument.is_identity())
        });
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

    /// The scalars a selector declaration (`{"select": "kind=thickness"}`)
    /// expands over, in a stable order: every film row's thickness,
    /// structure by structure; the factors and the thickness summand of
    /// every group; every top-level number of every material spec, and
    /// every knot value of a free-form one ([`FREEFORM_VALUES`]). Other
    /// fields (grading, roughness, gradients, oscillator entries) are
    /// declared by address.
    pub fn addresses(&self) -> Vec<ParamAddr> {
        let mut out = Vec::new();
        for (sid, s) in self.structures.iter().enumerate() {
            for row in 0..s.materials.len() {
                let film = self
                    .origin
                    .iter()
                    .position(|o| *o == Some((sid, row)))
                    .is_some_and(|e| self.entries[e].0.layer_type == LayerType::Film);
                if film {
                    out.push(ParamAddr::Layer {
                        structure: s.name.clone(),
                        layer: LayerRef::Index(row),
                        field: LayerField::Thickness,
                    });
                }
            }
        }
        let mut groups: Vec<&String> = self.groups.keys().collect();
        groups.sort();
        for g in groups {
            for field in [
                GroupField::NFactor,
                GroupField::KFactor,
                GroupField::ThickFactor,
                GroupField::ThickSummand,
            ] {
                out.push(ParamAddr::Group {
                    group: g.clone(),
                    field,
                });
            }
        }
        let mut materials: Vec<(&String, &Entry)> = self.materials.iter().collect();
        materials.sort_by(|a, b| a.0.cmp(b.0));
        for (name, entry) in materials {
            if let Entry::Spec(spec) = entry {
                for (key, v) in &spec.params {
                    if v.is_number() {
                        out.push(ParamAddr::Material {
                            material: name.clone(),
                            path: vec![PathSeg::Key(key.clone())],
                        });
                    } else if let Some(list) = v.as_array()
                        && FREEFORM_VALUES.contains(&(spec.model.as_str(), key.as_str()))
                    {
                        // A free-form model's values, knot by knot.
                        out.extend((0..list.len()).map(|i| ParamAddr::Material {
                            material: name.clone(),
                            path: vec![PathSeg::Key(key.clone()), PathSeg::Index(i)],
                        }));
                    }
                }
            }
        }
        out
    }

    /// Validation of every structure row and group, as the architect's own
    /// gate runs it; the first error refuses.
    pub fn check(&self) -> Result<(), String> {
        let mut seen = HashSet::new();
        for (e, (layer, _)) in self.entries.iter().enumerate() {
            let Some((sid, row)) = self.origin[e] else {
                continue;
            };
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
                    .filter(|e| self.origin[*e] == Some((sid, row)))
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
            ParamAddr::Data { dataset, field } => {
                if !self.instruments.contains_key(dataset) {
                    return Err(format!(
                        "{addr}: no dataset {dataset:?}; the problem has {:?}",
                        self.instruments.keys().collect::<Vec<_>>()
                    ));
                }
                Target::Data {
                    dataset: dataset.clone(),
                    field: *field,
                }
            }
        };
        self.read(&target).map_err(|e| format!("{addr}: {e}"))?;
        Ok(target)
    }

    /// The group governing a resolved target, if the model has one: a
    /// layer's material's group, or the named group.
    pub(crate) fn group_of(&self, target: &Target) -> Option<&Group> {
        match target {
            Target::Layer { entries, .. } => self.groups.get(&self.entries[entries[0]].0.material),
            Target::Group { group, .. } => self.groups.get(group),
            Target::Material { .. } | Target::Data { .. } => None,
        }
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
            Target::Data { dataset, field } => {
                Ok(*instrument_slot_ref(&self.instruments[dataset], *field))
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
            Target::Data { dataset, field } => {
                let i = self
                    .instruments
                    .get_mut(dataset)
                    .expect("resolved against this model");
                if *field == InstrumentField::Bandwidth && value < 0.0 {
                    return Err(format!("bandwidth {value} is negative"));
                }
                *instrument_slot(i, *field) = value;
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

fn instrument_slot(i: &mut Instrument, field: InstrumentField) -> &mut f64 {
    match field {
        InstrumentField::AngleOffset => &mut i.angle_offset,
        InstrumentField::WavelengthShift => &mut i.wavelength_shift,
        InstrumentField::Scale => &mut i.scale,
        InstrumentField::Offset => &mut i.offset,
        InstrumentField::Bandwidth => &mut i.bandwidth,
    }
}

fn instrument_slot_ref(i: &Instrument, field: InstrumentField) -> &f64 {
    match field {
        InstrumentField::AngleOffset => &i.angle_offset,
        InstrumentField::WavelengthShift => &i.wavelength_shift,
        InstrumentField::Scale => &i.scale,
        InstrumentField::Offset => &i.offset,
        InstrumentField::Bandwidth => &i.bandwidth,
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
