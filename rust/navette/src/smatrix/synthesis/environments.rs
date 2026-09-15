// SPDX-License-Identifier: LGPL-3.0-or-later
//! synthesis::environments — named design segments + per-environment
//! segment lists (F2.1).
//!
//! The multi-environment model in one paragraph: a design is one or more
//! **named segments**, each defined exactly once. An **environment** is an
//! ordered list of segments — either inline *fixed* layers (its own cover
//! glass, substrate, housing) or a reference to a design segment. The
//! design object is shared: every environment reads the same films, so a
//! thickness step or a needle insertion propagates everywhere by
//! construction. That is the whole point, and it is why the design is a
//! reference and not a copy (rejected alternative (d) in
//! `docs/plans/multi_environment_plan.md` §3.2).
//!
//! What this module does (F2.1): schema, validation, K assemblies, and
//! the routing table. What it deliberately does NOT do: evaluate. The
//! driver loop, `residuals_multi` and the K solves are F2.2; needle and
//! LM routing are F2.3; the Python surface is F2.4. Nothing here is
//! reachable from Python yet.
//!
//! **The absent-is-today rule.** A request with no `environments` is not
//! a one-environment request that happens to look flat — it takes the
//! flat door unchanged ([`super::design_config::build_design`]) and is
//! wrapped in a single environment afterwards. The bitwise gate is
//! therefore structural, not a promise: there is no second assembler for
//! a flat request to drift away from.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use serde::Deserialize;

use super::cycle::ContrastMap;
use super::design_config::{
    DesignRequest, LayerRow, RowAssembly, build_contrast, d_inh_delta, d_layer_type, d_true,
    eval_library, split_and_build_films,
};
use super::structure::DesignStack;
use crate::structure::Group;

// ---------------------------------------------------------------------------
// Request schema (additive-optional: absent == today's flat request)
// ---------------------------------------------------------------------------

/// One row of an environment's *fixed* surroundings.
///
/// A near-mirror of [`LayerRow`], and deliberately its own type for one
/// reason: `optimize` and `needle` default to **false** here and to
/// `true` there. That difference is what lets the compile tell "the
/// author said nothing" (force false, §3.1) apart from "the author asked
/// for a free variable in the surroundings" (refuse, §4.2). With a
/// shared type and `serde`'s `default = "d_true"`, both cases arrive as
/// `true` and the refusal would have to fire on every well-formed
/// request — so the schema carries the distinction instead of the code
/// guessing at it.
#[derive(Clone, Debug, Deserialize, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct FixedLayerRow {
    pub material_code: String,
    pub thickness_nm: f64,
    #[serde(default = "d_true")]
    pub coherent: bool,
    #[serde(default)]
    pub roughness_nm: f64,
    #[serde(default)]
    pub rough_type: i32,
    #[serde(default)]
    pub inhomogen: bool,
    #[serde(default = "d_inh_delta")]
    pub inh_delta: f64,
    #[serde(default)]
    pub interface: bool,
    #[serde(default)]
    pub interface_thickness_nm: f64,
    /// Refused when `true` (§4.2). Surroundings are not design variables.
    #[serde(default)]
    pub optimize: bool,
    /// Refused when `true` (§4.2).
    #[serde(default)]
    pub needle: bool,
    #[serde(default = "d_layer_type")]
    pub layer_type: i32,
    #[serde(default)]
    pub gradient: Option<crate::structure::gradient::GradientSpec>,
}

impl FixedLayerRow {
    /// Refuse an explicit free-variable flag, then hand back the row the
    /// assembler understands with both flags forced false.
    fn to_row(&self, label: &str) -> Result<LayerRow, String> {
        if self.optimize || self.needle {
            return Err(format!(
                "{label}: optimize/needle must be false - surroundings are \
                 not design variables"
            ));
        }
        Ok(LayerRow {
            material_code: self.material_code.clone(),
            thickness_nm: self.thickness_nm,
            coherent: self.coherent,
            roughness_nm: self.roughness_nm,
            rough_type: self.rough_type,
            inhomogen: self.inhomogen,
            inh_delta: self.inh_delta,
            interface: self.interface,
            interface_thickness_nm: self.interface_thickness_nm,
            optimize: false,
            needle: false,
            layer_type: self.layer_type,
            gradient: self.gradient.clone(),
        })
    }
}

/// One named design segment: an ordered film list, defined once.
#[derive(Clone, Debug, Deserialize, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct DesignSegmentCfg {
    pub layers: Vec<LayerRow>,
}

/// One entry in an environment's ordered stack: inline fixed layers, or
/// a reference to a design segment.
///
/// The plan writes this as a union (`{layers: [...]} | {design: "id"}`).
/// It is two `Option` fields plus an exactly-one check rather than
/// `#[serde(untagged)]` because an untagged enum reports only "data did
/// not match any variant" — it throws away *which* shape was wrong and
/// why, and every refusal in this module is required to name its
/// environment and segment.
#[derive(Clone, Debug, Default, Deserialize, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct EnvSegmentCfg {
    #[serde(default)]
    pub layers: Option<Vec<FixedLayerRow>>,
    #[serde(default)]
    pub design: Option<String>,
}

/// One environment: a name and its ordered segment list.
#[derive(Clone, Debug, Deserialize, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct EnvironmentCfg {
    pub name: String,
    pub stack: Vec<EnvSegmentCfg>,
}

/// The name a single-environment (flat) request compiles under.
///
/// §8 question 1, answered: named, not `None`-only. A demand tagged
/// `environment="default"` on a flat request must mean something, and an
/// error message that can say `known environments: "default"` beats one
/// that has to explain that there are environments but none of them have
/// names.
pub const DEFAULT_ENV: &str = "default";

// ---------------------------------------------------------------------------
// Compiled product
// ---------------------------------------------------------------------------

/// One shared design parameter's identity: which segment defined it,
/// where inside that segment, and under what (global) film name.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DesignSlot {
    /// Global film name — the cross-environment parameter identity
    /// (§2.4b). Positions differ per environment; this does not.
    pub name: Arc<str>,
    /// Defining segment id.
    pub segment: Arc<str>,
    /// Index within that segment's layer list.
    pub intra: usize,
}

/// K assemblies plus the table that maps a shared design parameter to
/// its absolute position in each of them.
///
/// The routing table is the answer to the second of the two gaps in
/// §2.4: positional identity cannot survive assemblies of different
/// length, so identity is the film NAME and the table carries the
/// positions. It is built once at compile and never consulted mid-eval
/// (§4.6).
#[derive(Clone, Debug)]
pub struct CompiledEnvironments {
    names: Vec<String>,
    stacks: Vec<DesignStack>,
    slots: Vec<DesignSlot>,
    /// `routing[slot][env]` = index into `stacks[env].spans()`.
    ///
    /// Spans, not rows: post-F0.1 a design film can expand to several
    /// solver rows, and the parameter belongs to the span. Building this
    /// row-keyed would have had to be rebuilt span-keyed after Phase A —
    /// which is the whole reason Phase B sits after gradients.
    routing: Vec<Vec<usize>>,
}

impl CompiledEnvironments {
    /// Environment names, in evaluation order. `names()[0]` is what an
    /// untagged demand resolves to (§4.2).
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// The K assembled stacks, in the same order as [`Self::names`].
    pub fn stacks(&self) -> &[DesignStack] {
        &self.stacks
    }

    /// Mutable view for the driver: F2.2 writes thicknesses back.
    pub fn stacks_mut(&mut self) -> &mut [DesignStack] {
        &mut self.stacks
    }

    /// Shared design parameters, in definition order (segment id, then
    /// position within the segment). Deterministic: the segment map is
    /// ordered.
    pub fn slots(&self) -> &[DesignSlot] {
        &self.slots
    }

    pub fn n_envs(&self) -> usize {
        self.names.len()
    }

    /// True when this compiles to exactly one environment — the `K == 1`
    /// side of §4.6's single branch.
    pub fn is_single(&self) -> bool {
        self.names.len() == 1
    }

    /// Resolve an environment name to its index.
    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }

    /// Span index of design parameter `slot` inside environment `env`.
    pub fn span_of(&self, slot: usize, env: usize) -> Option<usize> {
        self.routing.get(slot).and_then(|per| per.get(env)).copied()
    }

    /// Reverse lookup: which shared parameter (if any) owns span `span`
    /// of environment `env`. Surroundings answer `None` — they never
    /// deposit into the shared buckets (§4.4).
    pub fn slot_of_span(&self, env: usize, span: usize) -> Option<usize> {
        self.routing
            .iter()
            .position(|per| per.get(env) == Some(&span))
    }
}

// ---------------------------------------------------------------------------
// Compile
// ---------------------------------------------------------------------------

/// Auto-name for a surrounding film: `{env}.fixed[{seg}][{i}]` (§4.1).
///
/// Surroundings are auto-named only (§8 question 5, answered: auto-only
/// in v1). A user name there would address nothing — no buckets, no
/// contrast keys — and inviting one invites a `contrast` entry that
/// silently never splits.
fn fixed_name(env: &str, seg: usize, i: usize) -> String {
    format!("{env}.fixed[{seg}][{i}]")
}

/// Validate the design-segment registry and flatten it to slots.
///
/// Design film names are global (§3.1): they are the cross-environment
/// parameter identity, so two segments defining the same name is not a
/// shadowing question, it is an ambiguity, and it refuses naming both
/// segments.
fn compile_slots(design: &BTreeMap<String, DesignSegmentCfg>) -> Result<Vec<DesignSlot>, String> {
    let mut slots: Vec<DesignSlot> = Vec::new();
    let mut owner: HashMap<&str, &str> = HashMap::new();
    for (id, seg) in design {
        if seg.layers.is_empty() {
            return Err(format!(
                "design segment {id:?} is empty - a segment with no layers \
                 contributes no shared parameters"
            ));
        }
        for (i, row) in seg.layers.iter().enumerate() {
            if row.layer_type != 1 {
                return Err(format!(
                    "design segment {id:?} layer {i}: layer_type must be 1 \
                     (film) - half-spaces belong to an environment's fixed \
                     segments, not to the shared design"
                ));
            }
            let name = row.material_code.as_str();
            if let Some(prev) = owner.insert(name, id) {
                return Err(format!(
                    "design film {name:?} is defined by both segment {prev:?} \
                     and segment {id:?} - design film names are global (they \
                     ARE the cross-environment parameter identity)"
                ));
            }
            slots.push(DesignSlot {
                name: Arc::from(name),
                segment: Arc::from(id.as_str()),
                intra: i,
            });
        }
    }
    Ok(slots)
}

/// Validate one environment's segment list and flatten it to named rows.
///
/// Returns the rows in assembly order plus, for each design slot that
/// this environment carries, the authoring index the row landed at.
fn compile_env_rows(
    env: &EnvironmentCfg,
    design: &BTreeMap<String, DesignSegmentCfg>,
    slots: &[DesignSlot],
) -> Result<(Vec<(String, LayerRow)>, Vec<usize>), String> {
    let mut rows: Vec<(String, LayerRow)> = Vec::new();
    let mut refs: HashMap<&str, usize> = HashMap::new();
    // Authoring index per slot: films only, in the order the assembler
    // will see them. Half-space rows do not advance it.
    let mut film_idx: usize = 0;
    let mut slot_at: Vec<Option<usize>> = vec![None; slots.len()];

    for (s, seg) in env.stack.iter().enumerate() {
        match (&seg.layers, &seg.design) {
            (Some(_), Some(_)) | (None, None) => {
                return Err(format!(
                    "environment {:?} segment {s}: exactly one of \"layers\" \
                     or \"design\" must be set",
                    env.name
                ));
            }
            (Some(layers), None) => {
                for (i, fixed) in layers.iter().enumerate() {
                    let name = fixed_name(&env.name, s, i);
                    let row = fixed.to_row(&format!("{name:?}"))?;
                    if row.layer_type == 1 {
                        film_idx += 1;
                    }
                    rows.push((name, row));
                }
            }
            (None, Some(id)) => {
                let seg_def = design.get(id.as_str()).ok_or_else(|| {
                    let known = design
                        .keys()
                        .map(|k| format!("{k:?}"))
                        .collect::<Vec<_>>()
                        .join(", ");
                    let known = if known.is_empty() {
                        "none defined".to_string()
                    } else {
                        known
                    };
                    format!(
                        "environment {:?} segment {s}: unknown design segment \
                         {id:?} (known: {known})",
                        env.name
                    )
                })?;
                *refs.entry(id.as_str()).or_insert(0) += 1;
                for (i, row) in seg_def.layers.iter().enumerate() {
                    let name = row.material_code.clone();
                    let slot = slots
                        .iter()
                        .position(|sl| &*sl.segment == id.as_str() && sl.intra == i)
                        .expect("every design layer has a slot (compile_slots)");
                    slot_at[slot] = Some(film_idx);
                    film_idx += 1;
                    rows.push((name, row.clone()));
                }
            }
        }
    }

    // Exact-once (§4.2, §8 question 4 answered: refuse in v1). Omission
    // today is almost certainly a typo'd environment, and silent
    // non-contribution would lie about jointness.
    for id in design.keys() {
        match refs.get(id.as_str()).copied().unwrap_or(0) {
            1 => {}
            0 => {
                return Err(format!(
                    "environment {:?} does not reference design segment \
                     {id:?} - every environment must reference each design \
                     segment exactly once",
                    env.name
                ));
            }
            n => {
                return Err(format!(
                    "environment {:?} references design segment {id:?} {n} \
                     times - every environment must reference each design \
                     segment exactly once",
                    env.name
                ));
            }
        }
    }

    let at = slot_at
        .into_iter()
        .map(|o| o.expect("exact-once check fills every slot"))
        .collect();
    Ok((rows, at))
}

/// Compile a request into K assemblies plus the routing table.
///
/// A request with no `environments` is NOT re-assembled here: it goes
/// through [`build_design`](super::design_config::build_design) exactly
/// as it did before F2.1 and is wrapped in one environment named
/// [`DEFAULT_ENV`]. There is no second code path for a flat request to
/// drift from — which is what makes "old flat calls assemble
/// bitwise-identical stacks" structural rather than a claim under test.
pub fn build_environments(
    req: &DesignRequest,
    wavelengths: &[f64],
) -> Result<(CompiledEnvironments, ContrastMap, Vec<String>), String> {
    if req.environments.is_empty() {
        if !req.design.is_empty() {
            return Err(
                "\"design\" segments are defined but \"environments\" is empty \
                 - a design segment is only reachable through an environment \
                 that references it"
                    .to_string(),
            );
        }
        let (stack, cmap, warnings) = super::design_config::build_design(req, wavelengths)?;
        let slots: Vec<DesignSlot> = stack
            .spans()
            .iter()
            .map(|sp| DesignSlot {
                name: Arc::from(&*stack.films()[sp.bulk_start].material),
                segment: Arc::from(DEFAULT_ENV),
                intra: sp.logical,
            })
            .collect();
        let routing = (0..slots.len()).map(|i| vec![i]).collect();
        return Ok((
            CompiledEnvironments {
                names: vec![DEFAULT_ENV.to_string()],
                stacks: vec![stack],
                slots,
                routing,
            },
            cmap,
            warnings,
        ));
    }

    if !req.structure.layers.is_empty() {
        return Err(format!(
            "structure.layers must be empty when \"environments\" is set \
             (the shared films live in \"design\" segments, the surroundings \
             in each environment's fixed segments) - got {} rows",
            req.structure.layers.len()
        ));
    }
    if req.design.is_empty() {
        return Err(
            "\"environments\" is set but no \"design\" segments are defined \
             - there would be no shared parameters to optimize"
                .to_string(),
        );
    }

    // Environment names: non-empty and unique. They prefix every
    // auto-generated surrounding name, so a duplicate would silently
    // alias two environments' films onto one nk entry.
    let mut seen: HashSet<&str> = HashSet::new();
    for env in &req.environments {
        if env.name.is_empty() {
            return Err("environment names must be non-empty".to_string());
        }
        if !seen.insert(env.name.as_str()) {
            return Err(format!("duplicate environment name {:?}", env.name));
        }
    }

    let slots = compile_slots(&req.design)?;
    let nk_base = eval_library(&req.library, wavelengths)?;
    let groups: HashMap<String, Group> = req
        .structure
        .groups
        .iter()
        .map(|row| super::design_config::build_group(row).map(|g| (row.name.clone(), g)))
        .collect::<Result<_, String>>()?;

    // One nk table for all environments. Surrounding films are named
    // per environment, so the names cannot collide across environments
    // and one aliased table serves every assembly - no per-environment
    // copy of the spectra.
    let mut per_env_rows: Vec<Vec<(String, LayerRow)>> = Vec::with_capacity(req.environments.len());
    let mut per_env_slot_at: Vec<Vec<usize>> = Vec::with_capacity(req.environments.len());
    let mut nk = nk_base.clone();
    for env in &req.environments {
        let (rows, at) = compile_env_rows(env, &req.design, &slots)?;
        for (name, row) in &rows {
            if name != &row.material_code {
                // Auto-named surrounding: alias its name onto the
                // library spectrum its material code resolves to.
                let spectrum = nk_base.get(row.material_code.as_str()).ok_or_else(|| {
                    format!("material code {:?} not found in library", row.material_code)
                })?;
                nk.insert(Arc::from(name.as_str()), spectrum.clone());
            }
        }
        per_env_rows.push(rows);
        per_env_slot_at.push(at);
    }

    // A design film whose material code happens to spell an auto-name
    // would make two different films one nk entry. Cheap to check, and
    // the alternative is a silent alias.
    for slot in &slots {
        if slot.name.contains(".fixed[") {
            return Err(format!(
                "design film {:?} collides with the auto-naming scheme for \
                 surroundings ({{env}}.fixed[{{seg}}][{{i}}]) - rename the \
                 material",
                slot.name
            ));
        }
    }

    let cmap = build_contrast(&req.contrast, &nk)?;

    let mut stacks: Vec<DesignStack> = Vec::with_capacity(req.environments.len());
    let mut warnings: Vec<String> = Vec::new();
    for (e, env) in req.environments.iter().enumerate() {
        let borrowed: Vec<(String, &LayerRow)> = per_env_rows[e]
            .iter()
            .map(|(n, r)| (n.clone(), r))
            .collect();
        let (amb, sub, films) = split_and_build_films(
            &borrowed,
            &RowAssembly {
                nk: &nk,
                film_flags: &req.film_flags,
                per_film_flags: &req.per_film_flags,
                ambient_name: &req.ambient_name,
                substrate_name: &req.substrate_name,
                wavelengths,
                gate_label: "build_environments",
            },
        )
        .map_err(|e| format!("environment {:?}: {e}", env.name))?;

        let background: HashSet<String> = films
            .iter()
            .filter(|l| (l.inhomogen || l.gradient.is_some()) && !l.optimize && !l.needle)
            .map(|l| l.material.clone())
            .collect();
        let (stack, warns) =
            DesignStack::from_design(amb, sub, &films, &nk, &groups, wavelengths, &background)
                .map_err(|e| format!("environment {:?}: {e}", env.name))?;
        warnings.extend(warns.into_iter().map(|w| format!("{}: {w}", env.name)));
        stacks.push(stack);
    }

    // design_slot -> (env, span). `Span::logical` is the authoring film
    // index, which is what `compile_env_rows` recorded per slot.
    let mut routing: Vec<Vec<usize>> = vec![Vec::with_capacity(stacks.len()); slots.len()];
    for (e, stack) in stacks.iter().enumerate() {
        for (slot, per) in routing.iter_mut().enumerate() {
            let want = per_env_slot_at[e][slot];
            let idx = stack
                .spans()
                .iter()
                .position(|sp| sp.logical == want)
                .ok_or_else(|| {
                    format!(
                        "environment {:?}: design film {:?} lost its span at \
                         assembly (authoring index {want})",
                        req.environments[e].name, slots[slot].name
                    )
                })?;
            per.push(idx);
        }
    }

    let names = req.environments.iter().map(|e| e.name.clone()).collect();
    Ok((
        CompiledEnvironments {
            names,
            stacks,
            slots,
            routing,
        },
        cmap,
        warnings,
    ))
}

/// Resolve a demand's `environment` tag to an index.
///
/// `None` is environment 0 (§4.2), which is the one thing that keeps
/// every pre-F2.1 target set meaningful without an edit. An unknown name
/// refuses naming the demand and the environments that do exist — a
/// typo'd tag would otherwise evaluate against the wrong surroundings
/// and look like a physics result.
pub fn resolve_env(tag: Option<&str>, names: &[String], demand: &str) -> Result<u32, String> {
    let Some(tag) = tag else {
        return Ok(0);
    };
    if let Some(i) = names.iter().position(|n| n == tag) {
        return Ok(i as u32);
    }
    let known = if names.is_empty() {
        format!("{DEFAULT_ENV:?}")
    } else {
        names
            .iter()
            .map(|n| format!("{n:?}"))
            .collect::<Vec<_>>()
            .join(", ")
    };
    Err(format!(
        "{demand}: unknown environment {tag:?} (known environments: {known})"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;

    use super::super::design_config::{MaterialDef, StructureCfg};

    fn lib() -> Vec<MaterialDef> {
        [("L", 1.45), ("H", 2.1), ("G", 1.52)]
            .into_iter()
            .map(|(name, n)| MaterialDef {
                name: name.to_string(),
                code: None,
                model: "Konstant".to_string(),
                params: [("n".to_string(), Value::from(n))].into_iter().collect(),
                n_data: None,
                k_data: None,
            })
            .collect()
    }

    fn film(code: &str, d: f64) -> LayerRow {
        LayerRow {
            material_code: code.to_string(),
            thickness_nm: d,
            coherent: true,
            roughness_nm: 0.0,
            rough_type: 0,
            inhomogen: false,
            inh_delta: 0.1,
            interface: false,
            interface_thickness_nm: 0.0,
            optimize: true,
            needle: true,
            layer_type: 1,
            gradient: None,
        }
    }

    fn fixed(code: &str, d: f64, layer_type: i32) -> FixedLayerRow {
        FixedLayerRow {
            material_code: code.to_string(),
            thickness_nm: d,
            coherent: true,
            roughness_nm: 0.0,
            rough_type: 0,
            inhomogen: false,
            inh_delta: 0.1,
            interface: false,
            interface_thickness_nm: 0.0,
            optimize: false,
            needle: false,
            layer_type,
            gradient: None,
        }
    }

    fn req(design: Vec<(&str, Vec<LayerRow>)>, envs: Vec<EnvironmentCfg>) -> DesignRequest {
        DesignRequest {
            structure: StructureCfg {
                label: "t".to_string(),
                layers: vec![],
                groups: vec![],
            },
            library: lib(),
            contrast: BTreeMap::new(),
            film_flags: BTreeMap::new(),
            per_film_flags: BTreeMap::new(),
            ambient_name: "air".to_string(),
            substrate_name: "sub".to_string(),
            design: design
                .into_iter()
                .map(|(id, layers)| (id.to_string(), DesignSegmentCfg { layers }))
                .collect(),
            environments: envs,
        }
    }

    fn flat(layers: Vec<LayerRow>) -> DesignRequest {
        DesignRequest {
            structure: StructureCfg {
                label: "t".to_string(),
                layers,
                groups: vec![],
            },
            library: lib(),
            contrast: BTreeMap::new(),
            film_flags: BTreeMap::new(),
            per_film_flags: BTreeMap::new(),
            ambient_name: "air".to_string(),
            substrate_name: "sub".to_string(),
            design: BTreeMap::new(),
            environments: vec![],
        }
    }

    fn design_seg(id: &str) -> EnvSegmentCfg {
        EnvSegmentCfg {
            layers: None,
            design: Some(id.to_string()),
        }
    }

    fn fixed_seg(rows: Vec<FixedLayerRow>) -> EnvSegmentCfg {
        EnvSegmentCfg {
            layers: Some(rows),
            design: None,
        }
    }

    const WL: [f64; 3] = [500.0, 550.0, 600.0];

    #[test]
    fn k1_segmented_assembles_bitwise_identically_to_flat() {
        // The F2.1 gate. Not "close": the same `LayerSpec` values and the
        // same span partition, compared field for field.
        let films = vec![film("L", 100.0), film("H", 80.0)];
        let (flat_env, _, _) = build_environments(&flat(films.clone()), &WL).unwrap();
        let seg = req(
            vec![("coat", films)],
            vec![EnvironmentCfg {
                name: "only".to_string(),
                stack: vec![fixed_seg(vec![]), design_seg("coat"), fixed_seg(vec![])],
            }],
        );
        let (multi, _, _) = build_environments(&seg, &WL).unwrap();

        let a = &flat_env.stacks()[0];
        let b = &multi.stacks()[0];
        assert_eq!(a.films(), b.films(), "films differ");
        assert_eq!(a.spans(), b.spans(), "spans differ");
        assert_eq!(a.ambient(), b.ambient(), "ambient differs");
        assert_eq!(a.substrate(), b.substrate(), "substrate differs");
        assert_eq!(multi.n_envs(), 1);
        assert!(multi.is_single());
    }

    #[test]
    fn two_environments_share_the_design_at_different_positions() {
        let films = vec![film("L", 100.0), film("H", 80.0)];
        let seg = req(
            vec![("coat", films)],
            vec![
                EnvironmentCfg {
                    name: "bare".to_string(),
                    stack: vec![design_seg("coat"), fixed_seg(vec![fixed("G", 0.0, 2)])],
                },
                EnvironmentCfg {
                    name: "laminated".to_string(),
                    stack: vec![
                        fixed_seg(vec![fixed("G", 1000.0, 1), fixed("L", 25.0, 1)]),
                        design_seg("coat"),
                        fixed_seg(vec![fixed("G", 0.0, 2)]),
                    ],
                },
            ],
        );
        let (env, _, _) = build_environments(&seg, &WL).unwrap();
        assert_eq!(env.n_envs(), 2);
        assert!(!env.is_single());
        assert_eq!(env.slots().len(), 2, "two shared parameters");
        assert_eq!(&*env.slots()[0].name, "L");
        assert_eq!(&*env.slots()[1].name, "H");

        // Same parameter, different absolute positions - the gap §2.4(b)
        // names, closed by the table rather than by position.
        assert_eq!(env.span_of(0, 0), Some(0));
        assert_eq!(env.span_of(0, 1), Some(2));
        assert_eq!(env.span_of(1, 0), Some(1));
        assert_eq!(env.span_of(1, 1), Some(3));
        assert_eq!(env.slot_of_span(1, 2), Some(0));
        // Surroundings own no shared parameter.
        assert_eq!(env.slot_of_span(1, 0), None);
        assert_eq!(env.slot_of_span(1, 1), None);

        // The shared films are the same objects in both assemblies.
        for (slot, name) in [(0usize, "L"), (1, "H")] {
            for e in 0..2 {
                let sp = env.span_of(slot, e).unwrap();
                let row = env.stacks()[e].spans()[sp].bulk_start;
                assert_eq!(&*env.stacks()[e].films()[row].material, name);
                assert!(env.stacks()[e].films()[row].optimize);
            }
        }
        // Surroundings are frozen, whatever the request said.
        let lam = &env.stacks()[1];
        assert!(!lam.films()[0].optimize && !lam.films()[0].needle);
        assert!(!lam.films()[1].optimize && !lam.films()[1].needle);
    }

    #[test]
    fn refusals_name_the_environment_and_the_segment() {
        let films = vec![film("L", 100.0)];
        let one = |stack: Vec<EnvSegmentCfg>| {
            req(
                vec![("coat", films.clone())],
                vec![EnvironmentCfg {
                    name: "bare".to_string(),
                    stack,
                }],
            )
        };

        // Unknown design reference.
        let e = build_environments(&one(vec![design_seg("coatt")]), &WL).unwrap_err();
        assert!(e.contains("\"bare\""), "{e}");
        assert!(e.contains("segment 0"), "{e}");
        assert!(e.contains("\"coatt\""), "{e}");
        assert!(e.contains("known: \"coat\""), "{e}");

        // Omitted design segment.
        let e = build_environments(&one(vec![fixed_seg(vec![])]), &WL).unwrap_err();
        assert!(e.contains("\"bare\"") && e.contains("\"coat\""), "{e}");
        assert!(e.contains("exactly once"), "{e}");

        // Referenced twice.
        let e = build_environments(&one(vec![design_seg("coat"), design_seg("coat")]), &WL)
            .unwrap_err();
        assert!(e.contains("2 times") && e.contains("exactly once"), "{e}");

        // Both keys / neither key on a segment.
        let both = EnvSegmentCfg {
            layers: Some(vec![]),
            design: Some("coat".to_string()),
        };
        let e = build_environments(&one(vec![both]), &WL).unwrap_err();
        assert!(e.contains("exactly one of"), "{e}");
        let e = build_environments(&one(vec![EnvSegmentCfg::default()]), &WL).unwrap_err();
        assert!(e.contains("exactly one of"), "{e}");

        // A free variable in the surroundings.
        let mut greedy = fixed("G", 1000.0, 1);
        greedy.optimize = true;
        let e = build_environments(&one(vec![fixed_seg(vec![greedy]), design_seg("coat")]), &WL)
            .unwrap_err();
        assert!(e.contains("bare.fixed[0][0]"), "{e}");
        assert!(e.contains("optimize/needle must be false"), "{e}");
        assert!(e.contains("not design variables"), "{e}");
    }

    #[test]
    fn design_names_are_global_across_segments() {
        let seg = req(
            vec![
                ("front", vec![film("L", 100.0)]),
                ("back", vec![film("L", 60.0)]),
            ],
            vec![EnvironmentCfg {
                name: "e".to_string(),
                stack: vec![design_seg("front"), design_seg("back")],
            }],
        );
        let e = build_environments(&seg, &WL).unwrap_err();
        assert!(e.contains("\"L\""), "{e}");
        assert!(e.contains("\"back\"") && e.contains("\"front\""), "{e}");
    }

    #[test]
    fn flat_request_wraps_as_one_named_environment() {
        let (env, _, _) = build_environments(&flat(vec![film("L", 100.0)]), &WL).unwrap();
        assert_eq!(env.names(), [DEFAULT_ENV.to_string()]);
        assert_eq!(env.slots().len(), 1);
        assert_eq!(&*env.slots()[0].name, "L");
        assert_eq!(env.span_of(0, 0), Some(0));
    }

    #[test]
    fn a_design_without_environments_refuses() {
        let mut r = flat(vec![]);
        r.design.insert(
            "coat".to_string(),
            DesignSegmentCfg {
                layers: vec![film("L", 10.0)],
            },
        );
        let e = build_environments(&r, &WL).unwrap_err();
        assert!(e.contains("environments"), "{e}");
    }

    #[test]
    fn flat_films_and_environments_together_refuse() {
        let mut r = req(
            vec![("coat", vec![film("L", 100.0)])],
            vec![EnvironmentCfg {
                name: "e".to_string(),
                stack: vec![design_seg("coat")],
            }],
        );
        r.structure.layers = vec![film("H", 50.0)];
        let e = build_environments(&r, &WL).unwrap_err();
        assert!(e.contains("structure.layers must be empty"), "{e}");
        assert!(e.contains("got 1 rows"), "{e}");
    }

    #[test]
    fn untagged_demands_resolve_to_environment_zero() {
        let names = vec!["bare".to_string(), "laminated".to_string()];
        assert_eq!(resolve_env(None, &names, "spectral[0]").unwrap(), 0);
        assert_eq!(resolve_env(Some("bare"), &names, "spectral[0]").unwrap(), 0);
        assert_eq!(
            resolve_env(Some("laminated"), &names, "spectral[0]").unwrap(),
            1
        );
        let e = resolve_env(Some("bear"), &names, "spectral[0]").unwrap_err();
        assert!(e.contains("spectral[0]"), "{e}");
        assert!(e.contains("\"bear\""), "{e}");
        assert!(e.contains("\"bare\", \"laminated\""), "{e}");
    }
}
