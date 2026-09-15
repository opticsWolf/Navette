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

    /// F2.3: which shared parameter owns solver row `row` of environment
    /// `env`, if any.
    ///
    /// The needle scan produces sites as (row, depth-into-that-row); this
    /// is the first half of §4.4's "candidate locus translated back to
    /// (segment, intra-segment position)". Rows in fixed segments answer
    /// `None` — surroundings never deposit, and they are inadmissible
    /// hosts anyway, so this is a second lock on the same door.
    pub fn slot_of_row(&self, env: usize, row: usize) -> Option<usize> {
        let span = self.stacks[env]
            .spans()
            .iter()
            .position(|sp| row >= sp.start && row < sp.end)?;
        self.slot_of_span(env, span)
    }

    /// F2.3: the host row of shared parameter `slot` inside environment
    /// `env` — the second half of the locus translation.
    ///
    /// The bulk row, not the span start: an interface-carrying design film
    /// is a legal needle host (N1 — singleton-*bulk*), and its slice row is
    /// derived, never a host.
    pub fn host_row(&self, env: usize, slot: usize) -> Option<usize> {
        let span = *self.routing.get(slot)?.get(env)?;
        Some(self.stacks[env].spans()[span].bulk_start)
    }

    /// F2.3: split shared parameter `slot` in every environment at once.
    ///
    /// This is what makes "one insertion edits the single shared design
    /// object" true of the assemblies as well as of the model. Every
    /// environment's template is split at ITS OWN row for that slot, with
    /// the same depth and the same seed, so the K assemblies stay the same
    /// design in different surroundings rather than drifting apart after
    /// the first needle.
    ///
    /// The host becomes three parameters — top, seed, bottom — exactly as
    /// `insert_needle_seed` splits one span into three, so the slot list
    /// grows by two at `slot + 1` and every later slot's span index shifts
    /// by two. The caller splits the pipeline's design stack itself with
    /// the same `(row, depth, seed)`; `expand`'s alignment check is what
    /// catches it if the two ever disagree.
    pub fn insert_seed(
        &mut self,
        slot: usize,
        depth_into_layer_nm: f64,
        seed: &crate::smatrix::synthesis::structure::LayerSpec,
    ) -> Result<(), String> {
        if slot >= self.slots.len() {
            return Err(format!(
                "insert_seed: slot {slot} of {} does not exist",
                self.slots.len()
            ));
        }
        let hosts: Vec<usize> = (0..self.names.len())
            .map(|e| {
                self.host_row(e, slot).ok_or_else(|| {
                    format!(
                        "insert_seed: slot {slot} has no host in environment {:?}",
                        self.names[e]
                    )
                })
            })
            .collect::<Result<_, _>>()?;
        let host_spans: Vec<usize> = (0..self.names.len())
            .map(|e| self.routing[slot][e])
            .collect();

        for (e, &row) in hosts.iter().enumerate() {
            self.stacks[e]
                .insert_needle_seed(row, depth_into_layer_nm, seed.clone())
                .map_err(|msg| format!("environment {:?}: {msg}", self.names[e]))?;
        }

        // Spans after the host shifted by two in every environment.
        for per in self.routing.iter_mut() {
            for (e, s) in per.iter_mut().enumerate() {
                if *s > host_spans[e] {
                    *s += 2;
                }
            }
        }
        // The host's own entry stays put: the top portion keeps the span.
        let segment = self.slots[slot].segment.clone();
        let host_name = self.slots[slot].name.clone();
        self.slots.insert(
            slot + 1,
            DesignSlot {
                name: seed.material.clone(),
                segment: segment.clone(),
                intra: 0,
            },
        );
        self.slots.insert(
            slot + 2,
            DesignSlot {
                name: host_name,
                segment,
                intra: 0,
            },
        );
        self.routing
            .insert(slot + 1, host_spans.iter().map(|s| s + 1).collect());
        self.routing
            .insert(slot + 2, host_spans.iter().map(|s| s + 2).collect());
        // `intra` stops being the AUTHORING index the moment a needle
        // lands, so it is re-derived as the position within the segment -
        // which is what every reader after the compile actually wants.
        let mut next: HashMap<Arc<str>, usize> = HashMap::new();
        for sl in self.slots.iter_mut() {
            let n = next.entry(sl.segment.clone()).or_insert(0);
            sl.intra = *n;
            *n += 1;
        }
        Ok(())
    }

    /// F2.2: the K stacks to solve, given the shared design's CURRENT
    /// thicknesses.
    ///
    /// `design` is environment 0's stack — the object the pipeline carries
    /// and the thickness optimizer writes into. Every other environment is
    /// its compiled template with the design spans' row thicknesses copied
    /// across by [`Self::span_of`]. Copying per ROW rather than per span
    /// total is what keeps a graded design film's profile intact: the
    /// fractions never go through a divide.
    ///
    /// This is deliberately not a re-assembly. Under K > 1 the structure is
    /// frozen for F2.2 (needle insertion, cleanup removal and inflate are
    /// refused by the driver and are F2.3's), so the only thing that moves
    /// between evals is a thickness, and re-running `from_design` K times
    /// per eval would pay the 126 µs assembly (R6) for nothing. The
    /// alignment check below is what makes the assumption fail loudly
    /// instead of quietly mis-routing if it is ever violated.
    ///
    /// Surroundings are never touched: a fixed segment's rows keep the
    /// template's thicknesses, which is §3.1's "fully fixed" made
    /// operational rather than promised.
    pub fn expand(&self, design: &DesignStack) -> Result<Vec<DesignStack>, String> {
        self.check_alignment(design)?;
        let mut out = Vec::with_capacity(self.names.len());
        out.push(design.clone());
        for e in 1..self.names.len() {
            let mut st = self.stacks[e].clone();
            for slot in 0..self.slots.len() {
                let src = design.spans()[self.routing[slot][0]];
                let dst = self.stacks[e].spans()[self.routing[slot][e]];
                for (a, b) in (src.start..src.end).zip(dst.start..dst.end) {
                    st.set_thickness(b, design.films()[a].d_nm)?;
                }
            }
            out.push(st);
        }
        Ok(out)
    }

    /// The cross-environment alignment assert: the live design stack still
    /// has the span layout the compile recorded, and every slot still names
    /// the same film in every environment.
    ///
    /// Span layout, not row indices — the same object produces the same
    /// spans wherever it is embedded, and that is a stronger and cheaper
    /// statement than comparing positions across assemblies of different
    /// length. A failure here means a structural move slipped past the
    /// driver's refusal; it must not be allowed to become a silent
    /// mis-route.
    fn check_alignment(&self, design: &DesignStack) -> Result<(), String> {
        if design.spans().len() != self.stacks[0].spans().len() {
            return Err(format!(
                "environment '{}': the design stack now has {} spans, compiled \
                 with {} - a structural move under K > 1 is F2.3's, not F2.2's",
                self.names[0],
                design.spans().len(),
                self.stacks[0].spans().len()
            ));
        }
        for (slot, ds) in self.slots.iter().enumerate() {
            let src = design.spans()[self.routing[slot][0]];
            let n_src = src.end - src.start;
            for e in 0..self.names.len() {
                let dst = self.stacks[e].spans()[self.routing[slot][e]];
                if dst.end - dst.start != n_src {
                    return Err(format!(
                        "design film '{}' (segment '{}') is {} row(s) in \
                         environment '{}' and {} in '{}' - the shared design \
                         must expand identically everywhere",
                        ds.name,
                        ds.segment,
                        n_src,
                        self.names[0],
                        dst.end - dst.start,
                        self.names[e]
                    ));
                }
                if self.stacks[e].films()[dst.bulk_start].material
                    != design.films()[src.bulk_start].material
                {
                    return Err(format!(
                        "design film '{}' (segment '{}') routes to material \
                         '{}' in environment '{}' but '{}' in '{}' - the \
                         routing table and the stacks have drifted apart",
                        ds.name,
                        ds.segment,
                        design.films()[src.bulk_start].material,
                        self.names[0],
                        self.stacks[e].films()[dst.bulk_start].material,
                        self.names[e]
                    ));
                }
            }
        }
        Ok(())
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

    use super::super::design_config::DesignRequest;
    use super::super::evaluator::SmatrixContext;
    use super::super::merit::MeritSpec;

    use super::super::design_config::{MaterialDef, StructureCfg};

    fn lib() -> Vec<MaterialDef> {
        [("L", 1.45), ("H", 2.1), ("G", 1.52), ("M", 1.8)]
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

    // --- F2.2: K solves, joint merit -----------------------------------

    /// One `Rs = 0` demand at normal incidence, replicated across `k`
    /// environments. `k = 1` is the pre-F2.2 spec exactly.
    fn env_spec(k: usize) -> MeritSpec {
        use crate::smatrix::synthesis::merit::{
            ConstraintKind, CurveId, MeritKey, MeritSpec, MeritTarget, SimTransform,
        };
        let mut spec = MeritSpec::new();
        spec.set_n_envs(k).unwrap();
        let key = spec.add_key(MeritKey {
            angle: 0.0,
            curve: CurveId::Rs,
        }) as u32;
        for e in 0..k {
            spec.add_target(MeritTarget {
                key_idx: key,
                env_idx: e as u32,
                wavelengths: Arc::from(WL.as_slice()),
                kind: ConstraintKind::Exact,
                transform: SimTransform::Linear,
                norm_factor: 1.0,
                normalized_targets: Arc::from([0.0, 0.0, 0.0].as_slice()),
                tolerances: Arc::from([0.01, 0.01, 0.01].as_slice()),
                band: Arc::from([].as_slice()),
                phase: false,
                differential_passes: None,
                integral: false,
                weight: 1.0,
                count_norm: None,
            })
            .unwrap();
        }
        spec
    }

    fn ctx(spec: MeritSpec, envs: Option<Arc<CompiledEnvironments>>) -> SmatrixContext {
        SmatrixContext {
            wavls: WL.to_vec(),
            sin_theta: vec![0.0],
            spec,
            clamp_min_nm: 2.0,
            clamp_max_nm: 1000.0,
            thin_layer_policy: crate::smatrix::synthesis::config::ThinLayerPolicy::Remove,
            lm: crate::smatrix::synthesis::thick_opt::LmConfig::default(),
            clamp_accumulator: Default::default(),
            envs,
        }
    }

    /// The §4.6 (a) gate at the merit level: with K = 1 the branch that
    /// exists for K > 1 must not move a single bit. Compared on the raw
    /// bit pattern, because a tolerance here would pass exactly the drift
    /// the gate is for.
    #[test]
    fn k1_joint_merit_is_the_flat_merit_bitwise() {
        use crate::smatrix::synthesis::context::DesignContext;
        let films = vec![film("L", 100.0), film("H", 80.0)];
        let (envs, _, _) = build_environments(&flat(films), &WL).unwrap();
        let stack = envs.stacks()[0].clone();

        let flat_ctx = ctx(env_spec(1), None);
        let seg_ctx = ctx(env_spec(1), Some(Arc::new(envs)));

        let a = flat_ctx.evaluate_merit(&stack).unwrap();
        let b = seg_ctx.evaluate_merit(&stack).unwrap();
        assert_eq!(a.to_bits(), b.to_bits(), "K=1 branch moved: {a} vs {b}");
        assert!(a > 0.0, "the twin needs a merit with content, got {a}");
    }

    /// §4.6 (c): two environments with the same surroundings are the same
    /// physics twice, so the joint merit is exactly twice the single one.
    /// Exactly, not nearly — `x + x` is exact in binary floating point,
    /// and any inexactness here would mean the two solves disagree.
    #[test]
    fn two_identical_environments_double_the_merit() {
        use crate::smatrix::synthesis::context::DesignContext;
        let films = vec![film("L", 100.0), film("H", 80.0)];
        let (one, _, _) = build_environments(
            &req(
                vec![("coat", films.clone())],
                vec![EnvironmentCfg {
                    name: "a".to_string(),
                    stack: vec![design_seg("coat")],
                }],
            ),
            &WL,
        )
        .unwrap();
        let (two, _, _) = build_environments(
            &req(
                vec![("coat", films)],
                vec![
                    EnvironmentCfg {
                        name: "a".to_string(),
                        stack: vec![design_seg("coat")],
                    },
                    EnvironmentCfg {
                        name: "b".to_string(),
                        stack: vec![design_seg("coat")],
                    },
                ],
            ),
            &WL,
        )
        .unwrap();
        assert_eq!(two.n_envs(), 2);
        let stack = one.stacks()[0].clone();

        let m1 = ctx(env_spec(1), Some(Arc::new(one)))
            .evaluate_merit(&stack)
            .unwrap();
        let m2 = ctx(env_spec(2), Some(Arc::new(two)))
            .evaluate_merit(&stack)
            .unwrap();
        assert_eq!((2.0 * m1).to_bits(), m2.to_bits(), "{m1} vs {m2}");
    }

    /// The hand oracle, made exact: one shared design film behind two
    /// DIFFERENT cover sequences must solve exactly as two separate flat
    /// coatings would. Compared curve by curve on the bit pattern rather
    /// than against a 1e-12 numpy transfer matrix — the reference here is
    /// the engine's own flat door, which is the thing the multi path is
    /// claiming to reproduce.
    #[test]
    fn a_shared_film_behind_two_covers_solves_as_two_flat_coatings() {
        use crate::smatrix::synthesis::merit::CurveId;
        let design = vec![film("L", 120.0), film("H", 70.0)];
        let (envs, _, _) = build_environments(
            &req(
                vec![("coat", design.clone())],
                vec![
                    EnvironmentCfg {
                        name: "bare".to_string(),
                        stack: vec![fixed_seg(vec![]), design_seg("coat")],
                    },
                    EnvironmentCfg {
                        name: "laminated".to_string(),
                        stack: vec![
                            fixed_seg(vec![fixed("G", 300.0, 1), fixed("M", 40.0, 1)]),
                            design_seg("coat"),
                        ],
                    },
                ],
            ),
            &WL,
        )
        .unwrap();

        let base = envs.stacks()[0].clone();
        let solver = ctx(env_spec(2), Some(Arc::new(envs)));
        let sims = solver.simulate_all(&base).unwrap();

        // The two flat equivalents, built through the ordinary door.
        let mut cover_g = film("G", 300.0);
        cover_g.optimize = false;
        cover_g.needle = false;
        let mut cover_m = film("M", 40.0);
        cover_m.optimize = false;
        cover_m.needle = false;
        let flats = [
            design.clone(),
            vec![cover_g, cover_m, design[0].clone(), design[1].clone()],
        ];
        for (e, rows) in flats.into_iter().enumerate() {
            let (fe, _, _) = build_environments(&flat(rows), &WL).unwrap();
            let want = ctx(env_spec(1), None).simulate(&fe.stacks()[0]).unwrap();
            for id in [CurveId::Rs, CurveId::Rp, CurveId::Ts, CurveId::Tp] {
                let got = sims[e].curve(id).expect("curve present");
                let exp = want.curve(id).expect("curve present");
                assert_eq!(got.len(), exp.len(), "env {e} curve {id:?} length");
                for (i, (g, x)) in got.iter().zip(exp.iter()).enumerate() {
                    assert_eq!(
                        g.to_bits(),
                        x.to_bits(),
                        "env {e} curve {id:?} point {i}: {g} vs {x}"
                    );
                }
            }
        }
    }

    /// `expand` moves the shared thicknesses and nothing else: the
    /// surroundings are §3.1's "fully fixed", so a design step must not
    /// disturb one row of them.
    #[test]
    fn expand_moves_the_design_and_leaves_the_surroundings_alone() {
        let (envs, _, _) = build_environments(
            &req(
                vec![("coat", vec![film("L", 100.0)])],
                vec![
                    EnvironmentCfg {
                        name: "a".to_string(),
                        stack: vec![design_seg("coat")],
                    },
                    EnvironmentCfg {
                        name: "b".to_string(),
                        stack: vec![fixed_seg(vec![fixed("G", 300.0, 1)]), design_seg("coat")],
                    },
                ],
            ),
            &WL,
        )
        .unwrap();

        let mut design = envs.stacks()[0].clone();
        design.set_thickness(0, 137.0).unwrap();
        let out = envs.expand(&design).unwrap();
        assert_eq!(out.len(), 2);
        // Environment 0 IS the design object.
        assert_eq!(out[0].films()[0].d_nm, 137.0);
        // Environment 1: cover untouched, shared film moved.
        assert_eq!(out[1].films().len(), 2);
        assert_eq!(out[1].films()[0].d_nm, 300.0, "the cover moved");
        assert_eq!(out[1].films()[1].d_nm, 137.0, "the shared film did not");
        // And the template is not mutated by an expansion.
        assert_eq!(envs.stacks()[1].films()[1].d_nm, 100.0);
    }

    /// A structural move under K > 1 is F2.3's. If one ever reaches
    /// `expand`, it must name the drift rather than route thicknesses
    /// into the wrong spans.
    #[test]
    fn expand_refuses_a_stack_that_no_longer_matches_the_compile() {
        let (envs, _, _) = build_environments(
            &req(
                vec![("coat", vec![film("L", 100.0), film("H", 80.0)])],
                vec![
                    EnvironmentCfg {
                        name: "a".to_string(),
                        stack: vec![design_seg("coat")],
                    },
                    EnvironmentCfg {
                        name: "b".to_string(),
                        stack: vec![design_seg("coat")],
                    },
                ],
            ),
            &WL,
        )
        .unwrap();
        let (other, _, _) = build_environments(&flat(vec![film("L", 100.0)]), &WL).unwrap();
        let e = envs.expand(&other.stacks()[0]).unwrap_err();
        assert!(e.contains("spans"), "{e}");
        assert!(e.contains("F2.3"), "{e}");
    }

    /// Phase-A interaction (plan §1.4 + amendment A4): a graded film in a
    /// FIXED environment segment expands with its full profile, pinned.
    ///
    /// The background rule is what does the physics, and it is implied
    /// rather than declared — a profiled film with neither `optimize` nor
    /// `needle` is background. Surroundings are always both-false, so
    /// EVERY profiled fixed row is background, and A4's two plumbing lines
    /// in `build_environments` are what carry that. This twin is their
    /// regression test: the same film expressed as a background-pinned
    /// flat design must give bit-equal rows.
    ///
    /// Compared on `nk` and thickness rather than on the whole `LayerSpec`
    /// because the two doors name the film differently on purpose — the
    /// flat door uses the material code, the environment door auto-names
    /// `<env>.fixed[<seg>][<i>]` (§8 question 5). The profile is what is
    /// under test, and the profile is the numbers.
    #[test]
    fn a_graded_film_in_a_fixed_segment_keeps_its_profile() {
        let mut graded = fixed("M", 90.0, 1);
        graded.inhomogen = true;
        let (envs, _, _) = build_environments(
            &req(
                vec![("coat", vec![film("L", 100.0)])],
                vec![EnvironmentCfg {
                    name: "e0".to_string(),
                    stack: vec![fixed_seg(vec![graded]), design_seg("coat")],
                }],
            ),
            &WL,
        )
        .unwrap();

        let mut flat_graded = film("M", 90.0);
        flat_graded.inhomogen = true;
        flat_graded.optimize = false;
        flat_graded.needle = false;
        let (want, _, _) =
            build_environments(&flat(vec![flat_graded, film("L", 100.0)]), &WL).unwrap();

        let got = &envs.stacks()[0];
        let exp = &want.stacks()[0];
        assert!(
            got.films().len() > 2,
            "the profile was homogenized: {} rows",
            got.films().len()
        );
        assert_eq!(got.films().len(), exp.films().len(), "row count");
        assert_eq!(got.spans().len(), exp.spans().len(), "span count");
        for (i, (a, b)) in got.films().iter().zip(exp.films()).enumerate() {
            assert_eq!(a.d_nm.to_bits(), b.d_nm.to_bits(), "row {i} thickness");
            assert_eq!(a.optimize, b.optimize, "row {i} optimize");
            assert_eq!(a.needle, b.needle, "row {i} needle");
            for (w, (x, y)) in a.nk.iter().zip(b.nk.iter()).enumerate() {
                assert_eq!(x.re.to_bits(), y.re.to_bits(), "row {i} nk[{w}].re");
                assert_eq!(x.im.to_bits(), y.im.to_bits(), "row {i} nk[{w}].im");
            }
        }
        // And it is still fixed: only the design film is a parameter.
        assert_eq!(envs.slots().len(), 1);
        assert_eq!(&*envs.slots()[0].name, "L");
    }

    // --- F2.2: the joint run -------------------------------------------

    /// Thickness-only pipeline config: no needle, no cleanup, no inflate —
    /// the three structural moves F2.2 refuses under K > 1.
    fn thickness_only_cfg() -> crate::smatrix::synthesis::config::PipelineConfig {
        crate::smatrix::synthesis::config::PipelineConfig {
            max_macro_cycles: 2,
            needles_per_cycle: 0,
            enable_cleanup: false,
            enable_inflate: false,
            // F2.3: `remove` (the default) eliminates a sub-floor film on
            // the final sweep, and an elimination is the one structural
            // move the compile cannot follow. Joint runs take the policy
            // needle runs were already documented to take.
            thin_layer_policy: crate::smatrix::synthesis::config::ThinLayerPolicy::ClampUpFinal,
            ..Default::default()
        }
    }

    fn two_env_req(k: usize) -> DesignRequest {
        let envs = (0..k)
            .map(|i| EnvironmentCfg {
                name: format!("e{i}"),
                stack: vec![design_seg("coat")],
            })
            .collect();
        req(
            vec![("coat", vec![film("L", 100.0), film("H", 80.0)])],
            envs,
        )
    }

    /// §4.6 (c) at the run level: K = 2 with identical surroundings
    /// optimizes to the design K = 1 reaches.
    ///
    /// The MERIT twin above is bitwise; this one is not, and the reason is
    /// worth stating rather than hiding behind a tolerance. Two things
    /// differ along the path even though the surface's minimum is
    /// identical: the joint run's residual vector is twice as long (so
    /// LM's damping arithmetic rounds differently), and until F2.3 a
    /// multi-environment spec declines the analytic Jacobian, so K = 2
    /// converges by finite differences while K = 1 converges analytically.
    /// Two different descents onto the same minimum agree to about 1e-7
    /// nm here — a genuine disagreement about WHERE the minimum is would
    /// be orders larger.
    #[test]
    fn a_joint_run_over_identical_environments_lands_where_one_does() {
        use crate::smatrix::synthesis::driver::run_environments;
        let run = |k: usize| {
            let (envs, cmap, _) = build_environments(&two_env_req(k), &WL).unwrap();
            run_environments(
                envs,
                cmap,
                &WL,
                &[0.0],
                &env_spec(k),
                thickness_only_cfg(),
                Default::default(),
                Default::default(),
                |_, _| Ok(()),
            )
            .unwrap()
        };
        let (r1, s1, _) = run(1);
        let (r2, s2, _) = run(2);

        assert_eq!(s1.films().len(), s2.films().len());
        for (i, (a, b)) in s1.films().iter().zip(s2.films()).enumerate() {
            assert!(
                (a.d_nm - b.d_nm).abs() < 1e-5,
                "film {i}: {} vs {}",
                a.d_nm,
                b.d_nm
            );
        }
        // And the joint merit really is the doubled one, at the optimum
        // the joint run itself found.
        assert!(
            (r2.final_mf - 2.0 * r1.final_mf).abs() <= 1e-7 * r2.final_mf.abs(),
            "{} vs 2 x {}",
            r2.final_mf,
            r1.final_mf
        );
    }

    /// A spec compiled against a different roster than the design would
    /// score demands against surroundings nobody asked for.
    #[test]
    fn the_run_refuses_a_spec_that_does_not_match_the_roster() {
        use crate::smatrix::synthesis::driver::run_environments;
        let (envs, cmap, _) = build_environments(&two_env_req(2), &WL).unwrap();
        let e = run_environments(
            envs,
            cmap,
            &WL,
            &[0.0],
            &env_spec(1),
            thickness_only_cfg(),
            Default::default(),
            Default::default(),
            |_, _| Ok(()),
        )
        .unwrap_err();
        assert!(e.contains("1 environment(s)") && e.contains("2 "), "{e}");
        assert!(e.contains("e0, e1"), "{e}");
    }

    // --- F2.3: the joint needle ----------------------------------------

    fn with_contrast(mut r: DesignRequest) -> DesignRequest {
        r.contrast = [
            ("L".to_string(), "H".to_string()),
            ("H".to_string(), "L".to_string()),
        ]
        .into_iter()
        .collect();
        r
    }

    /// `two_env_req` plus the L <-> H contrast map the needle sweep needs.
    fn needle_req(k: usize) -> DesignRequest {
        with_contrast(two_env_req(k))
    }

    /// Two environments that genuinely DIFFER: the same coat bare, and
    /// behind a laminate. The design rows sit at different depths, so a
    /// locus is only shared through the routing table.
    fn two_cover_req() -> DesignRequest {
        with_contrast(req(
            vec![("coat", vec![film("L", 100.0), film("H", 80.0)])],
            vec![
                EnvironmentCfg {
                    name: "bare".to_string(),
                    stack: vec![fixed_seg(vec![]), design_seg("coat")],
                },
                EnvironmentCfg {
                    name: "laminated".to_string(),
                    stack: vec![
                        fixed_seg(vec![fixed("G", 300.0, 1), fixed("M", 40.0, 1)]),
                        design_seg("coat"),
                    ],
                },
            ],
        ))
    }

    /// One macro cycle, one needle: enough to test the insertion without
    /// letting two descents diverge over repeated re-optimizations.
    fn one_needle_cfg() -> crate::smatrix::synthesis::config::PipelineConfig {
        crate::smatrix::synthesis::config::PipelineConfig {
            max_macro_cycles: 1,
            needles_per_cycle: 1,
            enable_cleanup: false,
            enable_inflate: false,
            thin_layer_policy: crate::smatrix::synthesis::config::ThinLayerPolicy::ClampUpFinal,
            ..Default::default()
        }
    }

    /// Each environment's own merit — the same single-environment door the
    /// flat path uses, applied to that environment's assembly.
    fn env_merits(envs: &CompiledEnvironments, stack: &DesignStack) -> Vec<f64> {
        let flat = ctx(env_spec(1), None);
        envs.expand(stack)
            .unwrap()
            .iter()
            .map(|st| flat.spec.merit(&flat.simulate(st).unwrap(), 1e6))
            .collect()
    }

    /// §4.6 (c) for the needle: K = 2 with identical surroundings inserts
    /// the seed K = 1 inserts, in the same place, at the same thickness.
    ///
    /// Bit-equality is not claimed and would be wrong to claim — the joint
    /// residual vector is twice as long, so the thickness LM that follows
    /// the insertion rounds differently (the same reason F2.2's run twin is
    /// not bitwise). What must be identical is the STRUCTURE: same rows,
    /// same materials, same order.
    #[test]
    fn a_joint_needle_over_identical_environments_lands_where_one_does() {
        use crate::smatrix::synthesis::driver::run_environments;
        let run = |k: usize| {
            let (envs, cmap, _) = build_environments(&needle_req(k), &WL).unwrap();
            run_environments(
                envs,
                cmap,
                &WL,
                &[0.0],
                &env_spec(k),
                one_needle_cfg(),
                Default::default(),
                Default::default(),
                |_, _| Ok(()),
            )
            .unwrap()
        };
        let (_, s1, _) = run(1);
        let (_, s2, _) = run(2);

        assert!(
            s1.films().len() > 2,
            "no needle was inserted at all - the gate tests nothing"
        );
        assert_eq!(
            s1.films().len(),
            s2.films().len(),
            "K = 2 inserted a different number of seeds than K = 1"
        );
        for (i, (a, b)) in s1.films().iter().zip(s2.films()).enumerate() {
            assert_eq!(a.material, b.material, "row {i}");
            // Relative, and loose on purpose. The two runs insert the same
            // seed at the same place and then descend differently - K = 2
            // differences the Jacobian and carries twice the residuals, so
            // the LM that follows the insertion rounds its way to a
            // slightly different point on the same minimum (F2.2's
            // CORRECTION 8, one insertion further along). A real
            // disagreement about WHERE to put the needle shows up in the
            // material sequence above, not in the fourth decimal.
            assert!(
                (a.d_nm - b.d_nm).abs() <= 1e-3 * a.d_nm.max(1.0),
                "row {i}: {} vs {} nm",
                a.d_nm,
                b.d_nm
            );
        }
    }

    /// §4.4's promise, made operational: ONE insertion edits the single
    /// shared design object and reaches every environment. The seed lands
    /// in the design segment everywhere — at DIFFERENT rows, because the
    /// laminate sits above it, under the SAME name, because a design
    /// parameter is identified by name and not by position.
    #[test]
    fn the_joint_needle_lands_in_the_design_of_every_environment() {
        use crate::smatrix::synthesis::cycle::{NeedleCycleConfig, run_needle_cycles};
        use crate::smatrix::synthesis::pipeline::SpectralInputs;

        let (envs, cmap, _) = build_environments(&two_cover_req(), &WL).unwrap();
        let before: Vec<usize> = envs.stacks().iter().map(|s| s.films().len()).collect();
        let mut stack = envs.stacks()[0].clone();
        let mut solver = ctx(env_spec(2), Some(Arc::new(envs)));
        let spectral = SpectralInputs::from_spec(&env_spec(2), &[0.0], &WL).unwrap();
        let hist = run_needle_cycles(
            &mut solver,
            &mut stack,
            &spectral,
            &cmap,
            &NeedleCycleConfig {
                max_needles: 1,
                ..Default::default()
            },
        )
        .unwrap();

        let ins = hist
            .iter()
            .find_map(|h| h.insertion.as_ref())
            .expect("the joint sweep found no site to insert at");
        let envs = solver.envs.as_ref().unwrap();

        // The split reached every template, not just the one the pipeline
        // carries: one host row became three, K times.
        for (e, st) in envs.stacks().iter().enumerate() {
            assert_eq!(
                st.films().len(),
                before[e] + 2,
                "environment {:?} did not split",
                envs.names()[e]
            );
        }

        // The shared design object and environment 0's template are the
        // same rows - that is what makes "the pipeline carries environment
        // 0's stack" true after a structural move, not just before one.
        assert_eq!(stack.films().len(), envs.stacks()[0].films().len());
        for (i, (a, b)) in stack
            .films()
            .iter()
            .zip(envs.stacks()[0].films())
            .enumerate()
        {
            assert_eq!(a.material, b.material, "row {i}");
        }

        // Positions differ, names match.
        let host_slot = envs
            .slot_of_row(0, ins.film_idx)
            .expect("host is a design row");
        let seed_slot = host_slot + 1;
        let mut seed_rows = Vec::new();
        for e in 0..envs.n_envs() {
            let row = envs
                .host_row(e, seed_slot)
                .expect("the seed is a parameter");
            assert_eq!(
                envs.stacks()[e].films()[row].material,
                ins.material,
                "environment {:?} seed material",
                envs.names()[e]
            );
            assert!(
                envs.stacks()[e].films()[row].needle,
                "environment {:?}: the seed is not a host for the next needle",
                envs.names()[e]
            );
            seed_rows.push(row);
        }
        assert_ne!(
            seed_rows[0], seed_rows[1],
            "the laminate adds two rows, so the shared seed cannot sit at the \
             same row in both environments"
        );

        // And the routing still routes: a design step reaches every
        // environment through the table that the insertion just rewrote.
        envs.expand(&stack)
            .expect("alignment survived the insertion");
    }

    /// The gradient gate. The joint slope at a shared locus is the SUM of
    /// the per-environment analytic slopes (that is what the sweep adds
    /// up), and that sum is the real derivative of the joint merit — not
    /// environment 0's derivative wearing a joint label.
    ///
    /// The reference is built from the public single-environment doors:
    /// `build_needle_targets_env` and `run_needle_pass`, one environment at
    /// a time. The measurement is a forward difference of the joint merit
    /// across an actual insertion, so it also proves the two halves of
    /// F2.3 agree — the sweep's arithmetic and `insert_seed`'s bookkeeping.
    #[test]
    fn the_joint_slope_is_the_sum_of_the_environments_slopes() {
        use crate::smatrix::synthesis::needle_pass::{
            NeedlePassInput, build_needle_targets_env, run_needle_pass,
        };

        let (envs, cmap, _) = build_environments(&two_cover_req(), &WL).unwrap();
        let base = envs.stacks()[0].clone();
        let spec = env_spec(2);
        let solver = ctx(spec.clone(), Some(Arc::new(envs)));
        let envs = solver.envs.as_ref().unwrap();

        // The locus: 40 nm into the first design film, seeded with L's
        // contrast partner.
        const SLOT: usize = 0;
        const DEPTH: f64 = 40.0;
        let partner = cmap
            .get(&base.films()[envs.host_row(0, SLOT).unwrap()].material)
            .expect("L has a contrast partner")
            .clone();

        // Reference: per environment, analytically, through the flat doors.
        let stacks = envs.expand(&base).unwrap();
        let mut analytic = 0.0f64;
        for (e, st) in stacks.iter().enumerate() {
            let sim = solver.simulate(st).unwrap();
            let fold = build_needle_targets_env(&spec, &[0.0], &WL, Some(&sim), e as u32).unwrap();
            let sa = st.solver_arrays();
            let res = run_needle_pass(
                &NeedlePassInput {
                    n_stack_cache: &sa.n_stack_cache,
                    thicknesses: &sa.thicknesses,
                    rough_types: &sa.rough_types,
                    rough_vals: &sa.rough_vals,
                    n_layers: sa.n_layers as usize,
                    wavls: &WL,
                    sin_theta: &[0.0],
                    fold: &fold,
                    needle_n_per_wav: &partner.nk,
                    start_idx: 0,
                    end_idx: (sa.n_layers - 1) as usize,
                    calc_s: true,
                    calc_p: true,
                },
                st.films(),
                st.spans(),
                2.0,
            )
            .unwrap();
            let host = envs.host_row(e, SLOT).unwrap();
            let (i, _) = res
                .sites
                .iter()
                .enumerate()
                .find(|(_, s)| s.film_idx == host && (s.depth_into_layer_nm - DEPTH).abs() < 1e-9)
                .expect("the locus is on this environment's scan grid");
            analytic += res.p_profile[i];
            analytic += fold.phi_gain_shift.iter().sum::<f64>();
        }

        // Measurement: insert a thin seed there, everywhere at once, and
        // read the joint merit on both sides.
        let delta = 1e-4;
        let seed = crate::smatrix::synthesis::structure::LayerSpec {
            material: partner.material.clone(),
            nk: partner.nk.clone(),
            d_nm: delta,
            coherent: true,
            rough_type: 0,
            rough_val: 0.0,
            optimize: true,
            needle: true,
        };
        let mut grown = envs.as_ref().clone();
        let mut stack = base.clone();
        stack
            .insert_needle_seed(envs.host_row(0, SLOT).unwrap(), DEPTH, seed.clone())
            .unwrap();
        grown.insert_seed(SLOT, DEPTH, &seed).unwrap();

        let f0 = solver
            .spec
            .merit_multi(&solver.simulate_all(&base).unwrap(), 1e6);
        let after = ctx(spec, Some(Arc::new(grown)));
        let f1 = after
            .spec
            .merit_multi(&after.simulate_all(&stack).unwrap(), 1e6);
        let measured = (f1 - f0) / delta;

        assert!(
            (measured - analytic).abs() <= 1e-4 * analytic.abs().max(1e-12),
            "joint slope: measured {measured}, analytic sum {analytic}"
        );
    }

    /// The N1 twin at K > 1: "singleton-**bulk**", so a design film that
    /// carries an interface slice is still a legal needle host — in every
    /// environment, not only in the one the scan happens to start from.
    #[test]
    fn an_interface_carrying_design_film_hosts_a_needle_in_every_environment() {
        use crate::smatrix::synthesis::needle_pass::build_scan_sites;

        let mut iface = film("L", 100.0);
        iface.interface = true;
        iface.interface_thickness_nm = 2.0;
        let r = with_contrast(req(
            vec![("coat", vec![iface, film("H", 80.0)])],
            vec![
                // Both covers are non-empty: an interface slice is
                // emitted only where the film HAS a neighbour above it, so
                // a design film that is topmost in one environment and
                // buried in another is not the same number of rows - which
                // F2.1's alignment check refuses, correctly and for its own
                // reasons. The N1 question is whether the slice-carrying
                // film hosts a needle, so both environments give it one.
                EnvironmentCfg {
                    name: "thin cover".to_string(),
                    stack: vec![fixed_seg(vec![fixed("M", 20.0, 1)]), design_seg("coat")],
                },
                EnvironmentCfg {
                    name: "laminated".to_string(),
                    stack: vec![
                        fixed_seg(vec![fixed("G", 300.0, 1), fixed("M", 40.0, 1)]),
                        design_seg("coat"),
                    ],
                },
            ],
        ));
        let (envs, _, _) = build_environments(&r, &WL).unwrap();
        let stacks = envs.expand(&envs.stacks()[0].clone()).unwrap();
        for (e, st) in stacks.iter().enumerate() {
            let host = envs.host_row(e, 0).expect("slot 0 has a host row");
            assert_eq!(
                envs.slot_of_row(e, host),
                Some(0),
                "environment {:?}: the bulk row lost its parameter",
                envs.names()[e]
            );
            let sites = build_scan_sites(st.films(), st.spans(), 2.0);
            assert!(
                sites.iter().any(|s| s.film_idx == host),
                "environment {:?}: the interface-carrying film hosts no site",
                envs.names()[e]
            );
        }
    }

    /// The anti-§3.2a gate. A joint run improves BOTH environments, and it
    /// does not reach environment 1's number by accident: optimizing
    /// environment 0 alone — the sequential recipe §3.2a warns about —
    /// leaves environment 1 worse off than the joint run does.
    #[test]
    fn a_two_environment_run_improves_both_environments() {
        use crate::smatrix::synthesis::driver::run_environments;

        let (envs, cmap, _) = build_environments(&two_cover_req(), &WL).unwrap();
        let base = envs.stacks()[0].clone();
        let start = env_merits(&envs, &base);
        let (_, joint, _) = run_environments(
            envs,
            cmap,
            &WL,
            &[0.0],
            &env_spec(2),
            thickness_only_cfg(),
            Default::default(),
            Default::default(),
            |_, _| Ok(()),
        )
        .unwrap();

        // The sequential recipe: environment 0, on its own.
        let (solo_envs, solo_cmap, _) = build_environments(&needle_req(1), &WL).unwrap();
        let (_, solo, _) = run_environments(
            solo_envs,
            solo_cmap,
            &WL,
            &[0.0],
            &env_spec(1),
            thickness_only_cfg(),
            Default::default(),
            Default::default(),
            |_, _| Ok(()),
        )
        .unwrap();

        let (envs, _, _) = build_environments(&two_cover_req(), &WL).unwrap();
        let joint_m = env_merits(&envs, &joint);
        let solo_m = env_merits(&envs, &solo);
        for e in 0..2 {
            assert!(
                joint_m[e] < start[e],
                "environment {e}: {} did not improve on {}",
                joint_m[e],
                start[e]
            );
        }
        assert!(
            joint_m[1] < solo_m[1],
            "the joint run left environment 1 at {} while optimizing \
             environment 0 alone reaches {} - the joint run is not joint",
            joint_m[1],
            solo_m[1]
        );
    }

    /// Elimination is a removal whoever ordered it. `remove` deletes a
    /// sub-floor film on the pipeline's final sweep, so a joint run would
    /// end with K templates carrying a span the design no longer has — the
    /// refusal names the policy and the one that works.
    #[test]
    fn a_joint_run_refuses_the_removing_thin_layer_policy() {
        use crate::smatrix::synthesis::driver::run_environments;
        let cfg = crate::smatrix::synthesis::config::PipelineConfig {
            max_macro_cycles: 1,
            needles_per_cycle: 0,
            enable_cleanup: false,
            enable_inflate: false,
            thin_layer_policy: crate::smatrix::synthesis::config::ThinLayerPolicy::Remove,
            ..Default::default()
        };
        let (envs, cmap, _) = build_environments(&two_env_req(2), &WL).unwrap();
        let e = run_environments(
            envs,
            cmap,
            &WL,
            &[0.0],
            &env_spec(2),
            cfg,
            Default::default(),
            Default::default(),
            |_, _| Ok(()),
        )
        .unwrap_err();
        assert!(e.contains("thin_layer_policy 'remove'"), "{e}");
        assert!(e.contains("clamp_up_final"), "{e}");
    }

    /// The other half of the same rule: mid-run, the context pins the floor
    /// as an LM bound rather than letting a film fall through it. A joint
    /// run therefore never loses a span, whatever the optimizer wants.
    #[test]
    fn a_joint_optimization_never_loses_a_span() {
        use crate::smatrix::synthesis::context::DesignContext;
        let mut thin = film("H", 2.5);
        thin.needle = false;
        let r = req(
            vec![("coat", vec![film("L", 100.0), thin])],
            (0..2)
                .map(|i| EnvironmentCfg {
                    name: format!("e{i}"),
                    stack: vec![design_seg("coat")],
                })
                .collect(),
        );
        let (envs, _, _) = build_environments(&r, &WL).unwrap();
        let spans_before = envs.stacks()[0].spans().len();

        // The hazard, demonstrated rather than assumed: let that film fall
        // below the floor and the removing sweep deletes the span outright.
        // Under K > 1 the compile would still be carrying it.
        let mut fallen = envs.stacks()[0].clone();
        fallen.set_thickness(1, 0.5).unwrap();
        fallen
            .clamp_all_policy(
                2.0,
                1000.0,
                crate::smatrix::synthesis::config::ThinLayerPolicy::Remove,
                false,
            )
            .unwrap();
        assert_eq!(
            fallen.spans().len(),
            spans_before - 1,
            "the removing sweep no longer removes - this test has gone stale"
        );

        // The joint context does not let it get there: the floor is the
        // LM's lower bound, and the sweep clamps up.
        let mut stack = envs.stacks()[0].clone();
        let mut solver = ctx(env_spec(2), Some(Arc::new(envs)));
        solver.optimize_thicknesses(&mut stack).unwrap();
        assert_eq!(
            stack.spans().len(),
            spans_before,
            "the joint clamp eliminated a shared design parameter"
        );
        for l in stack.films() {
            assert!(l.d_nm >= 2.0 - 1e-12, "{} fell through the floor", l.d_nm);
        }
        // And the compile still routes into it.
        solver
            .envs
            .as_ref()
            .unwrap()
            .expand(&stack)
            .expect("the compile still matches the design");
    }

    /// F2.3 gave the compile `insert_seed`, not its inverse: needles run
    /// jointly now, and the two moves that REMOVE or RESHUFFLE parameters
    /// are still refused rather than optimizing environment 0 alone.
    #[test]
    fn the_run_refuses_structural_moves_while_k_is_greater_than_one() {
        use crate::smatrix::synthesis::driver::run_environments;
        for (mutate, needle) in [
            (
                Box::new(
                    |c: &mut crate::smatrix::synthesis::config::PipelineConfig| {
                        c.enable_cleanup = true
                    },
                ) as Box<dyn Fn(&mut _)>,
                "enable_cleanup",
            ),
            (
                Box::new(
                    |c: &mut crate::smatrix::synthesis::config::PipelineConfig| {
                        c.enable_inflate = true
                    },
                ),
                "enable_inflate",
            ),
        ] {
            let mut cfg = thickness_only_cfg();
            mutate(&mut cfg);
            let (envs, cmap, _) = build_environments(&two_env_req(2), &WL).unwrap();
            let e = run_environments(
                envs,
                cmap,
                &WL,
                &[0.0],
                &env_spec(2),
                cfg,
                Default::default(),
                Default::default(),
                |_, _| Ok(()),
            )
            .unwrap_err();
            assert!(e.contains(needle), "{needle}: {e}");
            assert!(e.contains("insert_seed"), "{e}");
        }
        // K = 1 keeps every one of them.
        let (envs, cmap, _) = build_environments(&two_env_req(1), &WL).unwrap();
        assert!(
            run_environments(
                envs,
                cmap,
                &WL,
                &[0.0],
                &env_spec(1),
                crate::smatrix::synthesis::config::PipelineConfig {
                    max_macro_cycles: 1,
                    ..Default::default()
                },
                Default::default(),
                Default::default(),
                |_, _| Ok(()),
            )
            .is_ok()
        );
    }

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
