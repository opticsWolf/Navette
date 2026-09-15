// SPDX-License-Identifier: LGPL-3.0-or-later
//! navette::config — versioned program documents and section assembly.
//!
//! Rust-first port of the Python `config.program` + section builders:
//! envelope gating, legacy-flat detection, prefix namespacing, and
//! dependency-order assembly (materials → groups → structures →
//! architect) into native objects. File reading is JSON-only (YAML
//! stays Python-side authoring); the Python layer thins to dict
//! handover + native calls. Context-merging (live objects) stays
//! Python-side dict assembly.

use std::collections::{BTreeMap, HashMap};

use serde::Deserialize;
use serde_json::Value;

use crate::smatrix::synthesis::design_config::{
    GroupRow, LayerRow, MaterialDef, StructureCfg, apply_row, build_group, spec_from_def,
};
use crate::structure::{
    Architect, BlockKind, Entry, Group, Layer, LayerType, MaterialProvider, SharedGroup,
    SharedStructure, SpecProvider, Structure,
};

/// Program envelope version (single canonical gate) — what this build
/// WRITES.
///
/// F2.4 bumped it to 2 for `design:` / `environments:`. The bump is not
/// because a v1 binary would refuse those sections — it would not.
/// `load_program_parts` reads the `sections` mapping with bare lookups,
/// so a v1 build handed a v2 program DROPS them and optimizes a
/// single-environment stack that looks entirely plausible. Silent drop,
/// not refusal: §1.2's failure class, one gate over. The section
/// whitelist below closes the hole for every section added after this
/// one; the bump is what protects the documents written before it.
pub const PROGRAM_SCHEMA_VERSION: u32 = 2;

/// The oldest program envelope this build READS (F2.4, mirroring F1.4's
/// state gate). A v1 program has no `design`/`environments` sections, so
/// nothing about it is ambiguous to a v2 reader — it loads unchanged.
pub const MIN_READABLE_PROGRAM_SCHEMA_VERSION: u32 = 1;

/// Section names a `kind == "program"` document may carry.
///
/// The list exists because the lookups below are `sections.get(..)`: an
/// unrecognised name used to be dropped in silence, which is exactly how
/// a v2 document loses its environments on a v1 build. Refusing an
/// unknown name means the NEXT section added is a loud failure on an old
/// build instead of a quiet one.
const PROGRAM_SECTIONS: &[&str] = &[
    "materials",
    "groups",
    "structures",
    "architect",
    "design",
    "environments",
];

/// One architect block: structure label reference + placement.
#[derive(Clone, Debug, Deserialize, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct BlockCfg {
    pub structure: String,
    #[serde(default)]
    pub inverted: bool,
    #[serde(default = "d_repeat")]
    pub repeat_count: i64,
    #[serde(default)]
    pub label: String,
    #[serde(default = "d_kind")]
    pub kind: i32,
}

fn d_repeat() -> i64 {
    1
}

fn d_kind() -> i32 {
    1
}

/// Authoring-time validation (the pydantic bound inventory, natively).
/// Called by the PyO3 config constructors; assembly assumes validity.
impl MaterialDef {
    pub fn validate(&self) -> Result<(), String> {
        const MODELS: [&str; 6] = [
            "Konstant",
            "TableMaterial",
            "Cauchy",
            "CauchyUrbach",
            "Sellmeier",
            "SellmeierUrbach",
        ];
        if !MODELS.contains(&self.model.as_str()) {
            return Err(format!("Unknown model: {}", self.model));
        }
        let num = |key: &str| -> Result<f64, String> {
            self.params
                .get(key)
                .and_then(|v| v.as_f64())
                .ok_or_else(|| format!("{}: '{}' must be a number.", self.model, key))
        };
        let pos = |key: &str| -> Result<(), String> {
            if num(key)? <= 0.0 {
                return Err(format!("{}: '{}' must be > 0.", self.model, key));
            }
            Ok(())
        };
        let nonneg = |key: &str, default: f64| -> Result<(), String> {
            let v = self
                .params
                .get(key)
                .and_then(|v| v.as_f64())
                .unwrap_or(default);
            if v < 0.0 {
                return Err(format!("{}: '{}' must be >= 0.", self.model, key));
            }
            Ok(())
        };
        match self.model.as_str() {
            "Konstant" => {
                pos("n")?;
                nonneg("k", 0.0)?;
            }
            "TableMaterial" => {
                if self.n_data.is_none() {
                    return Err("TableMaterial requires n_data".to_string());
                }
                nonneg("n_factor", 1.0)?;
                nonneg("k_factor", 1.0)?;
            }
            "Cauchy" => {
                num("A")?;
                num("B")?;
                num("C")?;
            }
            "CauchyUrbach" => {
                num("A")?;
                num("B")?;
                num("C")?;
                pos("alpha0")?;
                pos("Eu")?;
                pos("lambda_g")?;
            }
            "Sellmeier" => {
                pos("B1")?;
                pos("C1")?;
                pos("B2")?;
                pos("C2")?;
                nonneg("B3", 0.0)?;
                nonneg("C3", 0.0)?;
            }
            _ => {
                pos("B1")?;
                pos("C1")?;
                pos("B2")?;
                pos("C2")?;
                nonneg("B3", 0.0)?;
                nonneg("C3", 0.0)?;
                pos("alpha0")?;
                pos("Eu")?;
                pos("lambda_g")?;
            }
        }
        Ok(())
    }
}

impl LayerRow {
    pub fn validate(&self) -> Result<(), String> {
        if !(self.thickness_nm > 0.0) {
            return Err("thickness_nm must be > 0.".to_string());
        }
        if self.roughness_nm < 0.0 {
            return Err("roughness_nm must be >= 0.".to_string());
        }
        crate::structure::RoughnessType::try_from_i32(self.rough_type)
            .map_err(|e| format!("rough_type: {e}"))?;
        if !(0.0..=1.0).contains(&self.inh_delta) {
            return Err("inh_delta must be in [0, 1].".to_string());
        }
        if self.interface_thickness_nm < 0.0 {
            return Err("interface_thickness_nm must be >= 0.".to_string());
        }
        match self.layer_type {
            0..=2 => {}
            t => return Err(format!("layer_type must be 0, 1 or 2 (got {t}).")),
        }
        // F1.5: the gradient spec's own rule surface (one source, N5) -
        // the same checks the Layer constructor's gate runs, at the
        // config door where the layer does not exist yet.
        if let Some(g) = &self.gradient {
            for issue in g.issues() {
                if issue.is_error() {
                    return Err(issue.message);
                }
            }
        }
        Ok(())
    }
}

impl GroupRow {
    pub fn validate(&self) -> Result<(), String> {
        if self.error_mask.len() != 6 {
            return Err(format!(
                "group {:?}: error_mask needs 6 entries.",
                self.name
            ));
        }
        if self.optimization_mask.len() != 7
            || self.optimization_mask.iter().any(|x| *x != 0 && *x != 1)
        {
            return Err(format!(
                "group {:?}: optimization_mask must be 7 binary entries.",
                self.name
            ));
        }
        for (key, v) in [
            ("thickness", self.thickness_error_type),
            ("n", self.n_error_type),
            ("k", self.k_error_type),
            ("inh_delta", self.inh_delta_error_type),
            ("roughness", self.roughness_error_type),
            ("interface", self.interface_error_type),
        ] {
            crate::structure::ErrorType::try_from_i32(v).map_err(|_| {
                format!("group {:?}: {key} error type must be 0, 1 or 2.", self.name)
            })?;
        }
        Ok(())
    }
}

impl StructureCfg {
    /// Nested validation: every layer + group row.
    pub fn validate(&self) -> Result<(), String> {
        for (i, row) in self.layers.iter().enumerate() {
            row.validate().map_err(|e| format!("layers[{i}]: {e}"))?;
        }
        for row in &self.groups {
            row.validate()?;
        }
        Ok(())
    }
}

impl BlockCfg {
    pub fn validate(&self) -> Result<(), String> {
        if self.repeat_count < 1 {
            return Err("repeat_count must be >= 1.".to_string());
        }
        crate::structure::BlockKind::try_from_i32(self.kind).map_err(|e| format!("kind: {e}"))?;
        Ok(())
    }
}

/// Everything a program file restores (absent sections stay empty).
#[derive(Debug, Default)]
pub struct LoadedProgram {
    pub name: Option<String>,
    pub materials: Option<SpecProvider>,
    pub groups: HashMap<String, Group>,
    /// Named structures as shared handles: architect blocks built by
    /// [`load_program_parts`] alias these (edits propagate both ways).
    pub structures: HashMap<String, SharedStructure>,
    pub architect: Option<Architect>,
    /// F2.4: named design segments, defined once and shared by every
    /// environment that references them.
    ///
    /// Carried as the typed request data rather than as built objects:
    /// a design segment has no meaning on its own — it is a piece of a
    /// `DesignRequest`, and the thing that gives it one is the
    /// environment roster below plus the material library above.
    pub design: BTreeMap<String, crate::smatrix::synthesis::environments::DesignSegmentCfg>,
    /// F2.4: the environment roster, in evaluation order. Empty is the
    /// flat case and is what every pre-F2.4 program document is.
    pub environments: Vec<crate::smatrix::synthesis::environments::EnvironmentCfg>,
}

fn px(value: &str, prefix: Option<&str>) -> String {
    match prefix {
        Some(p) => format!("{p}{value}"),
        None => value.to_string(),
    }
}

fn code_key(def: &MaterialDef) -> &str {
    def.code.as_deref().unwrap_or(&def.name)
}

// ---------------------------------------------------------------------------
// Envelope gate (mirrors `_gate` + `load_document`, JSON form)
// ---------------------------------------------------------------------------

/// Gate + classify a parsed document.
/// Returns `(kind, name, payload)`; legacy-flat maps to
/// `("materials"|"structure", None, section-content)`.
pub fn gate_document(raw: &Value) -> Result<(String, Option<String>, Value), String> {
    let top = raw
        .as_object()
        .ok_or_else(|| "document top level must be a mapping.".to_string())?;
    if !top.contains_key("kind") {
        // Legacy flat form (no envelope to gate).
        if let Some(items) = top.get("materials") {
            return Ok(("materials".to_string(), None, items.clone()));
        }
        if let Some(layers) = top.get("layers") {
            let mut map = serde_json::Map::new();
            map.insert("label".to_string(), Value::from("stack"));
            map.insert("layers".to_string(), layers.clone());
            map.insert(
                "groups".to_string(),
                top.get("groups").cloned().unwrap_or(Value::Array(vec![])),
            );
            return Ok(("structure".to_string(), None, Value::Object(map)));
        }
        return Err("legacy document needs 'materials' or 'layers' at top level.".to_string());
    }
    let kinds = ["materials", "groups", "structure", "architect", "program"];
    let kind = top
        .get("kind")
        .and_then(|v| v.as_str())
        .ok_or_else(|| "program kind must be a string.".to_string())?;
    if !kinds.contains(&kind) {
        return Err(format!(
            "program kind {kind:?} unknown (expected one of {}).",
            kinds.join(", ")
        ));
    }
    // F2.4: a RANGE, not a point, and each direction refuses with its own
    // message - below is a document this build has outgrown, above is one
    // written by a build that knows more, and the remedy differs (F1.4's
    // contract, repeated here for the second gate).
    match top.get("schema_version").and_then(|v| v.as_u64()) {
        Some(v)
            if (MIN_READABLE_PROGRAM_SCHEMA_VERSION as u64..=PROGRAM_SCHEMA_VERSION as u64)
                .contains(&v) => {}
        Some(v) if v > PROGRAM_SCHEMA_VERSION as u64 => {
            return Err(format!(
                "program schema_version {v} unsupported (this code reads \
                 {MIN_READABLE_PROGRAM_SCHEMA_VERSION}..={PROGRAM_SCHEMA_VERSION}); \
                 refusing a program written by a newer build - upgrade navette \
                 to read it."
            ));
        }
        other => {
            return Err(format!(
                "program schema_version {other:?} unsupported (this code reads \
                 {MIN_READABLE_PROGRAM_SCHEMA_VERSION}..={PROGRAM_SCHEMA_VERSION})."
            ));
        }
    }
    let payload_keys: &[&str] = match kind {
        "materials" => &["materials"],
        "groups" => &["groups"],
        "structure" => &["label", "layers", "groups"],
        "architect" => &["structures", "blocks"],
        _ => &["sections"],
    };
    let unknown: Vec<&String> = top
        .keys()
        .filter(|k| {
            *k != "schema_version"
                && *k != "kind"
                && *k != "name"
                && !payload_keys.contains(&k.as_str())
        })
        .collect();
    if !unknown.is_empty() {
        return Err(format!("unknown top-level keys: {unknown:?}."));
    }
    if kind == "program" {
        match top.get("sections") {
            Some(Value::Object(_)) => {}
            _ => return Err("program document needs a 'sections' mapping.".to_string()),
        }
    } else if top.contains_key("sections") {
        return Err(format!(
            "standalone {kind:?} document must not carry 'sections'."
        ));
    }
    let name = top
        .get("name")
        .and_then(|v| v.as_str())
        .map(|s| s.to_string());
    if kind == "program" {
        return Ok((kind.to_string(), name, top["sections"].clone()));
    }
    if kind == "materials" || kind == "groups" {
        let key = kind;
        let items = top.get(key).cloned().unwrap_or(Value::Array(vec![]));
        return Ok((kind.to_string(), name, items));
    }
    let mut payload = serde_json::Map::new();
    for (k, v) in top {
        if k != "schema_version" && k != "kind" && k != "name" {
            payload.insert(k.clone(), v.clone());
        }
    }
    Ok((kind.to_string(), name, Value::Object(payload)))
}

// ---------------------------------------------------------------------------
// Section assembly (each usable standalone or nested)
// ---------------------------------------------------------------------------

fn from_json<T: for<'de> serde::Deserialize<'de>>(v: &Value, what: &str) -> Result<T, String> {
    serde_json::from_value(v.clone()).map_err(|e| format!("{what}: invalid section: {e}"))
}

/// `SpecProvider` from a `materials` section (prefix-aware).
pub fn load_materials(
    items: &Value,
    grid: &[f64],
    prefix: Option<&str>,
) -> Result<SpecProvider, String> {
    let items: Vec<Value> = from_json(items, "materials")
        .map_err(|_| "materials section must be a list.".to_string())?;
    let mut entries = HashMap::new();
    for item in &items {
        let mut def: MaterialDef = from_json(item, "material")?;
        def.name = px(&def.name, prefix);
        let code = px(code_key(&def), prefix);
        def.code = Some(code.clone());
        entries.insert(code, Entry::Spec(spec_from_def(&def)?));
    }
    SpecProvider::new(entries, grid.to_vec())
}

/// `{name: Group}` from a `groups` section (prefix-aware).
pub fn load_groups(items: &Value, prefix: Option<&str>) -> Result<HashMap<String, Group>, String> {
    let items: Vec<Value> =
        from_json(items, "groups").map_err(|_| "groups section must be a list.".to_string())?;
    let mut out = HashMap::new();
    for item in &items {
        let mut row: GroupRow = from_json(item, "group")?;
        row.name = px(&row.name, prefix);
        out.insert(row.name.clone(), build_group(&row)?);
    }
    Ok(out)
}

fn layer_from_row(
    row: &LayerRow,
    provider: Option<&dyn MaterialProvider>,
    prefix: Option<&str>,
) -> Result<Layer, String> {
    let code = px(&row.material_code, prefix);
    if let Some(p) = provider
        && !p.contains(&code)
    {
        // NOTE: prefixed codes only resolve against prefixed providers;
        // unprefixed fallthrough mirrors the file-wins rule.
        let unprefixed = row.material_code.as_str();
        if !p.contains(unprefixed) {
            return Err(format!("material code {code:?} not found in provider"));
        }
        let mut layer = Layer::film(row.thickness_nm, unprefixed);
        apply_layer_fields(&mut layer, row)?;
        return Ok(layer);
    }
    let mut layer = Layer::film(row.thickness_nm, &code);
    apply_layer_fields(&mut layer, row)?;
    Ok(layer)
}

fn apply_layer_fields(layer: &mut Layer, row: &LayerRow) -> Result<(), String> {
    layer.layer_type = LayerType::Film;
    apply_row(layer, row)
}

/// Native `Structure` from a `structure` section.
/// Per-section `groups` merge over `library_groups` (own wins).
/// `materials` is required (mirrors the `KeyError` when absent).
pub fn load_structure(
    payload: &Value,
    materials: Option<&dyn MaterialProvider>,
    library_groups: &HashMap<String, Group>,
    prefix: Option<&str>,
) -> Result<Structure, String> {
    let cfg: StructureCfg = from_json(payload, "structure")?;
    assemble_named(&cfg, materials, library_groups, prefix)
}

/// `{label: Structure}` (duplicate labels raise).
pub fn load_named_structures(
    items: &Value,
    materials: Option<&dyn MaterialProvider>,
    library_groups: &HashMap<String, Group>,
    prefix: Option<&str>,
) -> Result<HashMap<String, Structure>, String> {
    let items: Vec<Value> = from_json(items, "structures")
        .map_err(|_| "structures section must be a list.".to_string())?;
    let mut out = HashMap::new();
    for item in &items {
        let cfg: StructureCfg = from_json(item, "structure")?;
        let label = px(&cfg.label, prefix);
        if out.contains_key(&label) {
            return Err(format!("duplicate structure label '{label}'."));
        }
        let st = assemble_named(&cfg, materials, library_groups, prefix)?;
        out.insert(label, st);
    }
    Ok(out)
}

fn assemble_named(
    cfg: &StructureCfg,
    materials: Option<&dyn MaterialProvider>,
    library_groups: &HashMap<String, Group>,
    prefix: Option<&str>,
) -> Result<Structure, String> {
    let provider = materials
        .ok_or_else(|| "structure section needs materials: no provider given.".to_string())?;
    let mut layers = Vec::with_capacity(cfg.layers.len());
    for row in &cfg.layers {
        layers.push(layer_from_row(row, Some(provider), prefix)?);
    }
    let mut merged: HashMap<String, SharedGroup> = HashMap::new();
    for (k, g) in library_groups {
        merged.insert(k.clone(), crate::structure::group::shared_group(g.clone()));
    }
    for row in &cfg.groups {
        let mut named = row.clone();
        named.name = px(&named.name, prefix);
        merged.insert(
            named.name.clone(),
            crate::structure::group::shared_group(build_group(&named)?),
        );
    }
    Ok(Structure {
        layers,
        groups: merged,
    })
}

/// Native `Architect` over shared handles: blocks alias the map entries
/// (shared-block invariant — matches the live-composition path).
///
/// Native `Architect`: blocks reference `structures` by label.
pub fn load_architect(
    payload: &Value,
    structures: &HashMap<String, Structure>,
    prefix: Option<&str>,
) -> Result<Architect, String> {
    let blocks: Vec<Value> = payload
        .get("blocks")
        .and_then(|v| v.as_array())
        .ok_or_else(|| "architect section needs 'blocks'.".to_string())?
        .clone();
    let mut arch = Architect::new();
    for (i, raw) in blocks.iter().enumerate() {
        let b: BlockCfg = from_json(raw, "block")?;
        let target = px(&b.structure, prefix);
        let st = structures
            .get(&target)
            .ok_or_else(|| format!("architect block {i}: unknown structure label '{target}'."))?;
        if b.repeat_count < 1 {
            return Err(format!("architect block {i}: repeat_count must be >= 1."));
        }
        let kind =
            BlockKind::try_from_i32(b.kind).map_err(|e| format!("architect block {i}: {e}"))?;
        let label = if b.label.is_empty() {
            b.label.clone()
        } else {
            px(&b.label, prefix)
        };
        arch.add_structure(st.clone(), b.inverted, b.repeat_count as usize, label, kind)
            .map_err(|e| format!("architect block {i}: {e}"))?;
    }
    Ok(arch)
}

/// Shared-handle variant of [`load_architect`] for whole-program restore:
/// blocks alias the structures-map entries (edits propagate both ways).
pub fn load_architect_shared(
    payload: &Value,
    structures: &HashMap<String, SharedStructure>,
    prefix: Option<&str>,
) -> Result<Architect, String> {
    let blocks: Vec<Value> = payload
        .get("blocks")
        .and_then(|v| v.as_array())
        .ok_or_else(|| "architect section needs 'blocks'.".to_string())?
        .clone();
    let mut arch = Architect::new();
    for (i, raw) in blocks.iter().enumerate() {
        let b: BlockCfg = from_json(raw, "block")?;
        let target = px(&b.structure, prefix);
        let st = structures
            .get(&target)
            .ok_or_else(|| format!("architect block {i}: unknown structure label '{target}'."))?;
        if b.repeat_count < 1 {
            return Err(format!("architect block {i}: repeat_count must be >= 1."));
        }
        let kind =
            BlockKind::try_from_i32(b.kind).map_err(|e| format!("architect block {i}: {e}"))?;
        let label = if b.label.is_empty() {
            b.label.clone()
        } else {
            px(&b.label, prefix)
        };
        arch.add_shared(st.clone(), b.inverted, b.repeat_count as usize, label, kind)
            .map_err(|e| format!("architect block {i}: {e}"))?;
    }
    Ok(arch)
}

// ---------------------------------------------------------------------------
// Full program (dependency order; file sections win)
// ---------------------------------------------------------------------------

/// Restore a program from parsed JSON text + evaluation grid.
/// Sections load in dependency order (materials → groups → structures →
/// architect); a standalone section document loads just that part.
pub fn load_program_json(text: &str, grid: &[f64]) -> Result<LoadedProgram, String> {
    load_program_json_prefixed(text, grid, None)
}

/// Prefix-namespaced variant (multi-load collisions).
pub fn load_program_json_prefixed(
    text: &str,
    grid: &[f64],
    prefix: Option<&str>,
) -> Result<LoadedProgram, String> {
    let raw: Value =
        serde_json::from_str(text).map_err(|e| format!("program: invalid JSON: {e}"))?;
    let (kind, name, payload) = gate_document(&raw)?;
    load_program_parts(&kind, name, &payload, grid, prefix)
}

fn load_program_parts(
    kind: &str,
    name: Option<String>,
    payload: &Value,
    grid: &[f64],
    prefix: Option<&str>,
) -> Result<LoadedProgram, String> {
    let mut prog = LoadedProgram {
        name,
        ..Default::default()
    };
    let sections: BTreeMap<String, Value> = if kind == "program" {
        payload
            .as_object()
            .ok_or_else(|| "program document needs a 'sections' mapping.".to_string())?
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect()
    } else {
        [(kind.to_string(), payload.clone())].into_iter().collect()
    };

    // F2.4: every lookup below is a bare `get`, so an unrecognised name
    // used to be dropped without a word - the same silence that makes a
    // v2 program lose its environments on a v1 build. Refuse it here, the
    // way unknown TOP-LEVEL keys are already refused, so the next section
    // added fails loudly on an old build.
    if kind == "program" {
        let unknown: Vec<&str> = sections
            .keys()
            .map(|k| k.as_str())
            .filter(|k| !PROGRAM_SECTIONS.contains(k))
            .collect();
        if !unknown.is_empty() {
            return Err(format!(
                "program section(s) {} unknown (expected one of {}) - refusing \
                 rather than dropping a payload this build does not read.",
                unknown.join(", "),
                PROGRAM_SECTIONS.join(", ")
            ));
        }
    }

    if let Some(items) = sections.get("materials") {
        prog.materials = Some(load_materials(items, grid, prefix)?);
    }
    // F2.4: shape-checked HERE, where the document is, so a typo in a
    // segment reference or a stray key is reported against the file
    // rather than surfacing later as a compile error with no filename in
    // it. The SEMANTIC checks (a reference to a segment that does not
    // exist, a free variable in the surroundings, a roster that does not
    // match the demands) stay in `build_environments` - one compiler, one
    // set of refusals.
    if let Some(items) = sections.get("design") {
        prog.design = serde_json::from_value(items.clone())
            .map_err(|e| format!("program section 'design': {e}"))?;
    }
    if let Some(items) = sections.get("environments") {
        prog.environments = serde_json::from_value(items.clone())
            .map_err(|e| format!("program section 'environments': {e}"))?;
    }
    if !prog.design.is_empty() && prog.environments.is_empty() {
        return Err(
            "program section 'design' is defined but 'environments' is not - a \
             design segment is only reachable through an environment that \
             references it."
                .to_string(),
        );
    }
    if let Some(items) = sections.get("groups") {
        prog.groups = load_groups(items, prefix)?;
    }
    let provider = prog.materials.as_ref().map(|p| p as &dyn MaterialProvider);
    // Wrap named structures as shared handles FIRST so architect blocks
    // below alias the same cores (shared-block invariant).
    let share = |st: Structure| SharedStructure::new(std::cell::RefCell::new(st));

    if kind == "structure" {
        let label = payload
            .get("label")
            .and_then(|v| v.as_str())
            .unwrap_or("stack");
        prog.structures.insert(
            px(label, prefix),
            share(load_structure(payload, provider, &prog.groups, prefix)?),
        );
    } else if let Some(items) = sections.get("structures") {
        for (label, st) in load_named_structures(items, provider, &prog.groups, prefix)? {
            prog.structures.insert(label, share(st));
        }
    }

    if kind == "architect" {
        if payload.get("structures").is_none() {
            return Err("standalone architect document needs 'structures' + 'blocks'.".to_string());
        }
        for (label, st) in
            load_named_structures(&payload["structures"], provider, &prog.groups, prefix)?
        {
            prog.structures.insert(label, share(st));
        }
    }
    if let Some(arch_payload) = sections.get("architect") {
        prog.architect = Some(load_architect_shared(
            arch_payload,
            &prog.structures,
            prefix,
        )?);
    }
    Ok(prog)
}

/// Restore a program from a JSON file.
pub fn load_program_file(
    path: &str,
    grid: &[f64],
    prefix: Option<&str>,
) -> Result<LoadedProgram, String> {
    let text =
        std::fs::read_to_string(path).map_err(|e| format!("program: cannot read {path}: {e}"))?;
    load_program_json_prefixed(&text, grid, prefix)
}

// ---------------------------------------------------------------------------
// Tests (standalone: no Python)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn grid() -> Vec<f64> {
        vec![500.0, 600.0]
    }

    fn program_json() -> String {
        serde_json::json!({
          "schema_version": 1,
          "kind": "program",
          "name": "demo",
          "sections": {
            "materials": [
              {"name": "L", "code": "L", "model": "Konstant",
               "params": {"n": 1.45}},
              {"name": "H", "code": "H", "model": "Konstant",
               "params": {"n": 2.1}},
            ],
            "groups": [
              {"name": "H", "thick_factor": 1.0},
            ],
            "structures": [
              {"label": "main", "layers": [
                {"material_code": "L", "thickness_nm": 100.0},
                {"material_code": "H", "thickness_nm": 60.0},
              ], "groups": []},
            ],
            "architect": {"blocks": [
              {"structure": "main", "label": "run", "repeat_count": 2},
            ]},
          }
        })
        .to_string()
    }

    #[test]
    fn full_program_restores() {
        let prog = load_program_json(&program_json(), &grid()).unwrap();
        assert_eq!(prog.name.as_deref(), Some("demo"));
        assert!(prog.materials.as_ref().unwrap().contains("L"));
        assert!(prog.groups.contains_key("H"));
        assert_eq!(prog.structures["main"].borrow().layers.len(), 2);
        let arch = prog.architect.unwrap();
        assert_eq!(arch.blocks.len(), 1);
        assert_eq!(arch.blocks[0].repeat_count, 2);
    }

    #[test]
    fn blocks_alias_named_structures() {
        let prog = load_program_json(&program_json(), &grid()).unwrap();
        let arch = prog.architect.as_ref().unwrap();
        let shared = &arch.blocks[0].structure;
        let named = prog.structures.get("main").unwrap();
        assert!(std::rc::Rc::ptr_eq(shared, named));
        // Edits propagate both ways through the shared handle.
        named.borrow_mut().layers.pop();
        assert_eq!(arch.blocks[0].structure.borrow().layers.len(), 1);
    }

    #[test]
    fn prefix_namespaces() {
        let prog = load_program_json_prefixed(&program_json(), &grid(), Some("p_")).unwrap();
        assert!(prog.materials.as_ref().unwrap().contains("p_L"));
        assert!(prog.groups.contains_key("p_H"));
        assert!(prog.structures.contains_key("p_main"));
    }

    #[test]
    fn gate_refuses() {
        let grid = grid();
        assert!(load_program_json("{\"kind\": \"program\"}", &grid).is_err()); // no version
        // F2.4: this assertion used to read "2 refuses"; the bump makes 2
        // the version this build writes, so BOTH endpoints of the range are
        // pinned here instead - the whole point of a range gate is that
        // neither end moves by accident.
        assert!(
            load_program_json(
                "{\"schema_version\": 1, \"kind\": \"program\", \"sections\": {}}",
                &grid
            )
            .is_ok(),
            "v1 is the oldest readable program and must still load"
        );
        assert!(
            load_program_json(
                "{\"schema_version\": 2, \"kind\": \"program\", \"sections\": {}}",
                &grid
            )
            .is_ok(),
            "v2 is what this build writes"
        );
        let newer = load_program_json(
            "{\"schema_version\": 3, \"kind\": \"program\", \"sections\": {}}",
            &grid,
        )
        .unwrap_err();
        assert!(newer.contains("newer build"), "{newer}");
        let stale = load_program_json(
            "{\"schema_version\": 0, \"kind\": \"program\", \"sections\": {}}",
            &grid,
        )
        .unwrap_err();
        assert!(stale.contains("reads 1..=2"), "{stale}");
        // An unknown SECTION refuses by name rather than being dropped.
        let bogus = load_program_json(
            "{\"schema_version\": 2, \"kind\": \"program\", \"sections\": {\"bogus\": []}}",
            &grid,
        )
        .unwrap_err();
        assert!(
            bogus.contains("bogus") && bogus.contains("unknown"),
            "{bogus}"
        );
        assert!(load_program_json("{\"schema_version\": 1, \"kind\": \"nope\"}", &grid).is_err());
        assert!(
            load_program_json(
                "{\"schema_version\": 1, \"kind\": \"groups\", \"groups\": [], \"bogus\": 1}",
                &grid
            )
            .is_err()
        );
        assert!(load_program_json("{\"layers\": []}", &grid).is_err()); // legacy, no materials key
    }

    #[test]
    fn legacy_flat_materials() {
        let prog = load_program_json(
      "{\"materials\": [{\"name\": \"L\", \"model\": \"Konstant\", \"params\": {\"n\": 1.45}}]}",
      &grid(),
    )
    .unwrap();
        assert!(prog.materials.as_ref().unwrap().contains("L"));
    }

    #[test]
    fn refs_must_resolve() {
        let mut v: Value = serde_json::from_str(&program_json()).unwrap();
        v["sections"]["structures"][0]["layers"][0]["material_code"] = Value::from("X");
        assert!(load_program_json(&v.to_string(), &grid()).is_err());
        let mut v: Value = serde_json::from_str(&program_json()).unwrap();
        v["sections"]["architect"]["blocks"][0]["structure"] = Value::from("ghost");
        assert!(load_program_json(&v.to_string(), &grid()).is_err());
    }
}
