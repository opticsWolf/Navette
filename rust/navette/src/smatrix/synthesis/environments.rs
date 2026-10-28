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
//! LM routing are F2.3; the Python surface is F2.4. All of it has
//! shipped: `run_design_environments` (navette-py) is the door, and
//! `run_needle(design=..., environments=[...])` the Python spelling of
//! the same request.
//!
//! **The absent-is-today rule.** A request with no `environments` is not
//! a one-environment request that happens to look flat — it takes the
//! flat door unchanged ([`super::design_config::build_design`]) and is
//! wrapped in a single environment afterwards. The bitwise gate is
//! therefore structural, not a promise: there is no second assembler for
//! a flat request to drift away from.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use num_complex::Complex64;
use serde::Deserialize;

use super::cycle::{MaterialPool, new_needle, pool_entry, pool_from_layers};
use super::design_config::{
    DesignRequest, LayerRow, RowAssembly, build_pool, d_inh_delta, d_layer_type, d_true,
    eval_library, split_and_build_films,
};
use super::needle_pass::{ScanSite, span_admissible};
use super::structure::{DesignStack, EdgeMaterials, LayerSpec, edge_allows, span_index_of_row};
use crate::structure::providers::MaterialProvider;
use crate::structure::{ExpandOptions, Group, Layer, RoughnessType};

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
///
/// `outer` and `inner` are the segment's own edge lists (pool labels):
/// the materials allowed for its first film, which faces the incoming
/// light, and for its last, which faces away from it. They are checked on
/// the segment's own top and bottom in every environment, whatever
/// surrounds it there. `None` leaves that edge free.
#[derive(Clone, Debug, Deserialize, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct DesignSegmentCfg {
    pub layers: Vec<LayerRow>,
    #[serde(default)]
    pub outer: Option<Vec<String>>,
    #[serde(default)]
    pub inner: Option<Vec<String>>,
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
///
/// `repeat` (design references only, default 1) places that many copies
/// of the segment back to back. Every copy is the same design: one
/// parameter per film, one needle landing in every copy.
///
/// `inverted` (design references only) places the segment last film
/// first, as the Architect expands an inverted block: each interface and
/// roughness stays with the plane it belongs to, so a copy is the
/// segment turned over, not its rows reversed. Its outer film still
/// carries the outer list and now faces away from the light.
#[derive(Clone, Debug, Default, Deserialize, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct EnvSegmentCfg {
    #[serde(default)]
    pub layers: Option<Vec<FixedLayerRow>>,
    #[serde(default)]
    pub design: Option<String>,
    #[serde(default)]
    pub repeat: Option<usize>,
    #[serde(default)]
    pub inverted: bool,
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
    /// `routing[slot][env]` = the indices into `stacks[env].spans()` of
    /// every copy of that film, in assembly order. `routing[slot][0][0]`
    /// is the primary: the copy in the design object (environment 0's
    /// stack) that the optimizer moves. Every other copy, environment 0's
    /// own included, is pinned (`optimize` false) and follows it.
    ///
    /// Spans, not rows: post-F0.1 a design film can expand to several
    /// solver rows, and the parameter belongs to the span. Building this
    /// row-keyed would have had to be rebuilt span-keyed after Phase A —
    /// which is the whole reason Phase B sits after gradients.
    routing: Vec<Vec<Vec<usize>>>,
    /// Edge lists per design segment id; a segment with neither list set
    /// has no entry. Empty for a flat request, whose stack carries its own.
    segment_edges: BTreeMap<Arc<str>, EdgeMaterials>,
    /// `inverted[env][span]`: the span is an entry of an inverted copy.
    /// Only design spans are ever inverted.
    inverted: Vec<Vec<bool>>,
    /// Per design film, what an inverted copy is re-emitted from.
    faces: Vec<Face>,
    /// The simulation grid the stacks were assembled on.
    wavelengths: Vec<f64>,
}

/// One design film as the expander sees it, group policy already applied:
/// its index, coherence, and the roughness and interface slice on its
/// forward top face. An inverted copy moves both to the plane they
/// belong to, which is why it is re-emitted from these rather than copied
/// row by row.
#[derive(Clone, Debug)]
struct Face {
    nk: Arc<[Complex64]>,
    coherent: bool,
    rough_val: f64,
    rough_type: i32,
    /// The interface slice's thickness as emitted, when the film has one.
    slice: Option<f64>,
}

/// A boundary locus of the shared design: the plane before film `k` of
/// design segment `seg` (an index into [`CompiledEnvironments::segments`]),
/// in the segment's forward order; `k` equal to the segment's length is
/// the plane after its last film. Every copy of the segment in every
/// environment has the locus once. Where two segments touch, the plane is
/// two loci, one per segment: a film there joins one of them.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Gap {
    pub seg: usize,
    pub k: usize,
}

/// One placement of a [`Gap`] in one environment's stack.
#[derive(Clone, Copy, Debug)]
pub(crate) struct GapPlace {
    pub env: usize,
    /// The span index the new film's span takes.
    pub at: usize,
    /// The new film's row once it is in; the rows above it stay put, so
    /// it is also the boundary's row before (`films.len()`: the substrate
    /// side).
    pub row: usize,
    /// The copy is inverted.
    pub inverted: bool,
    /// The slice opening span `at` moves into the new span, ahead of it.
    pub take_slice: bool,
}

/// Serves one inverted run's films to the expander, keyed by their
/// position in the run.
struct RunProvider<'a>(Vec<&'a [Complex64]>);

impl MaterialProvider for RunProvider<'_> {
    fn nk(&self, name: &str, _wavelengths: &[f64]) -> Result<Vec<Complex64>, String> {
        name.parse::<usize>()
            .ok()
            .and_then(|i| self.0.get(i))
            .map(|nk| nk.to_vec())
            .ok_or_else(|| format!("inverted run: no film {name:?}"))
    }

    fn contains(&self, name: &str) -> bool {
        name.parse::<usize>().is_ok_and(|i| i < self.0.len())
    }

    fn grid(&self) -> Option<&[f64]> {
        None
    }
}

fn rough_kind(code: i32) -> Result<RoughnessType, String> {
    Ok(match code {
        0 => RoughnessType::None,
        1 => RoughnessType::Linear,
        2 => RoughnessType::Step,
        3 => RoughnessType::Exponential,
        4 => RoughnessType::Gaussian,
        5 => RoughnessType::NevotCroce,
        _ => return Err(format!("unknown roughness type {code}")),
    })
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

    /// The edge lists of design segment `segment`, or `None` when both of
    /// its edges are free.
    pub fn segment_edges(&self, segment: &str) -> Option<&EdgeMaterials> {
        self.segment_edges.get(segment)
    }

    pub fn n_envs(&self) -> usize {
        self.names.len()
    }

    /// True when this compiles to exactly one environment — the `K == 1`
    /// side of §4.6's single branch.
    pub fn is_single(&self) -> bool {
        self.names.len() == 1
    }

    /// True when some design film is used more than once in some
    /// environment (a `repeat`, or a second reference).
    pub fn has_copies(&self) -> bool {
        self.routing
            .iter()
            .any(|per| per.iter().any(|copies| copies.len() > 1))
    }

    /// True when some design segment has an edge list.
    pub fn has_segment_edges(&self) -> bool {
        !self.segment_edges.is_empty()
    }

    /// True when some design segment is placed inverted somewhere.
    pub fn has_inverted(&self) -> bool {
        self.inverted.iter().any(|env| env.contains(&true))
    }

    /// True when the run must take the joint path: more than one
    /// environment, one environment with copies, or an inverted copy (its
    /// rows are re-emitted from the design at every evaluation). A single
    /// environment whose films are each used once, forward, is the flat
    /// run exactly.
    pub fn is_joint(&self) -> bool {
        !self.is_single() || self.has_copies() || self.has_inverted()
    }

    /// True when span `span` of environment `env` is an entry of an
    /// inverted copy.
    pub fn is_inverted(&self, env: usize, span: usize) -> bool {
        self.inverted
            .get(env)
            .and_then(|v| v.get(span))
            .copied()
            .unwrap_or(false)
    }

    /// True when the inverted entry at `span` has its own interface slice
    /// carved by the next entry of its run (the slice then sits at the
    /// start of that entry's span). The last entry of a run drops it.
    fn carves_next(&self, env: usize, span: usize, slot: usize) -> bool {
        self.is_inverted(env, span)
            && self.is_inverted(env, span + 1)
            && self.faces[slot].slice.is_some()
    }

    /// True when copy span `span` of design film `slot` in environment
    /// `env` emits the film's own interface slice. A forward copy emits it
    /// unless it is the first film of its stack (the incident plane is
    /// clean); an inverted copy emits it when the next entry of its run
    /// carves it. Every copy of a film agrees (checked at the compile).
    fn emits_slice(&self, env: usize, span: usize, slot: usize) -> bool {
        if self.is_inverted(env, span) {
            self.carves_next(env, span, slot)
        } else {
            self.stacks[env].spans()[span].slice
        }
    }

    /// The bulk thickness of design film `slot`, read from its primary in
    /// `design`. Every copy's bulk is the same slab: the copies agree on
    /// whether the film's own slice is emitted, and it is the same
    /// thickness wherever it is.
    fn forward_bulk(&self, design: &DesignStack, slot: usize) -> f64 {
        let p = self.routing[slot][0][0];
        design.films()[design.spans()[p].bulk_start].d_nm
    }

    /// The film's whole thickness, its own slice included when emitted:
    /// what an inverted copy is re-emitted at.
    fn film_total(&self, design: &DesignStack, slot: usize) -> f64 {
        let b = self.forward_bulk(design, slot);
        match self.faces[slot].slice {
            Some(s) if self.emits_slice(0, self.routing[slot][0][0], slot) => b + s,
            _ => b,
        }
    }

    /// The depth from the forward top of a design film's bulk of a point
    /// `depth` below the top of row `row` in environment `env`, when that
    /// row is the bulk of an inverted copy; `None` for every other row.
    pub(crate) fn forward_depth(
        &self,
        env: usize,
        row: usize,
        depth: f64,
        design: &DesignStack,
    ) -> Option<f64> {
        let span = span_index_of_row(self.stacks[env].spans(), row);
        if !self.is_inverted(env, span) {
            return None;
        }
        let (slot, _) = self.copy_of_span(env, span)?;
        Some(self.forward_bulk(design, slot) - depth)
    }

    /// The interior needle sites of `st`, environment `env`'s expanded
    /// stack. Forward rows get [`build_scan_sites`]' grid; the bulk row of
    /// an inverted copy gets the same grid measured from the film's
    /// forward top, so every copy of a film offers the same loci in the
    /// same order of depth into the forward film.
    ///
    /// [`build_scan_sites`]: super::needle_pass::build_scan_sites
    pub(crate) fn scan_sites(
        &self,
        env: usize,
        st: &DesignStack,
        design: &DesignStack,
        step: f64,
    ) -> Vec<ScanSite> {
        assert!(step > 0.0, "scan_step_nm must be positive");
        let spans = st.spans();
        let mut sites = Vec::new();
        let mut cumulative = 0.0f64;
        for (film_idx, layer) in st.films().iter().enumerate() {
            let d = layer.d_nm;
            if layer.needle && d > 0.0 && span_admissible(spans, film_idx) {
                let span = span_index_of_row(spans, film_idx);
                let mirrored = self
                    .is_inverted(env, span)
                    .then(|| self.copy_of_span(env, span))
                    .flatten();
                if let Some((slot, _)) = mirrored {
                    let b = self.forward_bulk(design, slot);
                    let n_steps = (b / step) as i32;
                    for k in (1..n_steps).rev() {
                        let pos_in_layer = b - k as f64 * step;
                        sites.push(ScanSite {
                            film_idx,
                            depth_into_layer_nm: pos_in_layer,
                            z_nm: cumulative + pos_in_layer,
                            boundary: false,
                        });
                    }
                } else {
                    let n_steps = (d / step) as i32;
                    for k in 1..n_steps {
                        let pos_in_layer = k as f64 * step;
                        sites.push(ScanSite {
                            film_idx,
                            depth_into_layer_nm: pos_in_layer,
                            z_nm: cumulative + pos_in_layer,
                            boundary: false,
                        });
                    }
                }
            }
            cumulative += d;
        }
        sites
    }

    /// Resolve an environment name to its index.
    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }

    /// Span index of design parameter `slot`'s first copy inside
    /// environment `env` (in environment 0, the primary).
    pub fn span_of(&self, slot: usize, env: usize) -> Option<usize> {
        self.spans_of(slot, env).first().copied()
    }

    /// Span indices of every copy of design parameter `slot` inside
    /// environment `env`, in assembly order.
    pub fn spans_of(&self, slot: usize, env: usize) -> &[usize] {
        self.routing
            .get(slot)
            .and_then(|per| per.get(env))
            .map_or(&[], Vec::as_slice)
    }

    /// How many places design parameter `slot` occupies, over every
    /// environment and every copy.
    pub fn n_placements(&self, slot: usize) -> usize {
        self.routing
            .get(slot)
            .map_or(0, |per| per.iter().map(Vec::len).sum())
    }

    /// Reverse lookup: which shared parameter (if any) owns span `span`
    /// of environment `env`. Surroundings answer `None` — they never
    /// deposit into the shared buckets (§4.4).
    pub fn slot_of_span(&self, env: usize, span: usize) -> Option<usize> {
        self.copy_of_span(env, span).map(|(slot, _)| slot)
    }

    /// Which shared parameter owns span `span` of environment `env`, and
    /// which of its copies there that span is.
    pub fn copy_of_span(&self, env: usize, span: usize) -> Option<(usize, usize)> {
        self.routing.iter().enumerate().find_map(|(slot, per)| {
            per.get(env)?
                .iter()
                .position(|&s| s == span)
                .map(|c| (slot, c))
        })
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
        self.copy_of_row(env, row).map(|(slot, _)| slot)
    }

    /// [`Self::slot_of_row`], plus which copy of the parameter the row
    /// belongs to.
    pub fn copy_of_row(&self, env: usize, row: usize) -> Option<(usize, usize)> {
        let span = self.stacks[env]
            .spans()
            .iter()
            .position(|sp| row >= sp.start && row < sp.end)?;
        self.copy_of_span(env, span)
    }

    /// F2.3: the host row of shared parameter `slot` inside environment
    /// `env` — the second half of the locus translation. The first copy's
    /// (in environment 0, the primary's).
    ///
    /// The bulk row, not the span start: an interface-carrying design film
    /// is a legal needle host (N1 — singleton-*bulk*), and its slice row is
    /// derived, never a host.
    pub fn host_row(&self, env: usize, slot: usize) -> Option<usize> {
        let span = self.span_of(slot, env)?;
        Some(self.stacks[env].spans()[span].bulk_start)
    }

    /// F2.3: split shared parameter `slot` in every environment at once.
    ///
    /// This is what makes "one insertion edits the single shared design
    /// object" true of the assemblies as well as of the model. Every
    /// environment's template is split at ITS OWN row for every copy of
    /// that slot, with the same depth and the same seed, so the K
    /// assemblies stay the same design in different surroundings rather
    /// than drifting apart after the first needle.
    ///
    /// The host becomes three parameters — top, seed, bottom — exactly as
    /// `insert_needle_seed` splits one span into three, so the slot list
    /// grows by two at `slot + 1` and every later span shifts by two per
    /// copy split before it. The caller splits the pipeline's design stack
    /// with [`Self::split_design`] first, from the same routing;
    /// `expand`'s alignment check is what catches it if the two ever
    /// disagree.
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
        for e in 0..self.names.len() {
            if self.routing[slot][e].is_empty() {
                return Err(format!(
                    "insert_seed: slot {slot} has no host in environment {:?}",
                    self.names[e]
                ));
            }
        }

        for e in 0..self.names.len() {
            // Last copy first: a split moves only the spans after it, so
            // the copies before keep their indices until their turn.
            let mut copies = self.routing[slot][e].clone();
            copies.reverse();
            for host in copies {
                // A template's thicknesses are stale (the design's are
                // copied in at every expansion), so an inverted copy is
                // split at its own row's mirror depth: only the shape
                // matters, and the re-emission sets the rows.
                let inv = self.is_inverted(e, host);
                let depth = if inv {
                    let d = self.stacks[e].films()[self.stacks[e].spans()[host].bulk_start].d_nm;
                    (d - depth_into_layer_nm).max(0.0)
                } else {
                    depth_into_layer_nm
                };
                split_copy(
                    &mut self.stacks[e],
                    host,
                    depth,
                    seed,
                    e == 0 && host != self.routing[slot][0][0],
                )
                .map_err(|msg| format!("environment {:?}: {msg}", self.names[e]))?;
                self.inverted[e].splice(host + 1..host + 1, [inv, inv]);
                for per in self.routing.iter_mut() {
                    for s in per[e].iter_mut() {
                        if *s > host {
                            *s += 2;
                        }
                    }
                }
            }
        }

        // The host's own entries stay put: the top portion keeps the span.
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
        let shifted = |by: usize| -> Vec<Vec<usize>> {
            self.routing[slot]
                .iter()
                .map(|copies| copies.iter().map(|s| s + by).collect())
                .collect()
        };
        let (seeds, bottoms) = (shifted(1), shifted(2));
        self.routing.insert(slot + 1, seeds);
        self.routing.insert(slot + 2, bottoms);
        // An inverted copy meets the pieces bottom first: the host's span
        // holds the bottom portion and the top portion comes after the
        // seed.
        for e in 0..self.names.len() {
            for c in 0..self.routing[slot][e].len() {
                let host = self.routing[slot][e][c];
                if self.is_inverted(e, host) {
                    self.routing[slot][e][c] = host + 2;
                    self.routing[slot + 2][e][c] = host;
                }
            }
        }
        // The top portion keeps the host's face; the bottom portion is a
        // clone of the host's bulk row (its roughness included) with no
        // slice; the seed is the seed.
        let bottom = Face {
            slice: None,
            ..self.faces[slot].clone()
        };
        self.faces.insert(
            slot + 1,
            Face {
                nk: seed.nk.clone(),
                coherent: seed.coherent,
                rough_val: seed.rough_val,
                rough_type: seed.rough_type,
                slice: None,
            },
        );
        self.faces.insert(slot + 2, bottom);
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

    /// Split every copy of shared parameter `slot` in the pipeline's design
    /// stack (environment 0's), as [`Self::insert_seed`] splits the
    /// templates. Call it before `insert_seed`: it reads the routing as it
    /// stands. The seed in the primary is a free parameter; the seeds in
    /// the other copies are pinned and follow it.
    ///
    /// `depth_into_layer_nm` is measured from the forward top of the
    /// film's bulk; an inverted copy is split at the mirror depth.
    pub fn split_design(
        &self,
        design: &mut DesignStack,
        slot: usize,
        depth_into_layer_nm: f64,
        seed: &crate::smatrix::synthesis::structure::LayerSpec,
    ) -> Result<(), String> {
        let copies = self.spans_of(slot, 0);
        let Some(&primary) = copies.first() else {
            return Err(format!(
                "split_design: slot {slot} has no host in environment {:?}",
                self.names[0]
            ));
        };
        let b = self.forward_bulk(design, slot);
        for &host in copies.iter().rev() {
            let depth = if self.is_inverted(0, host) {
                b - depth_into_layer_nm
            } else {
                depth_into_layer_nm
            };
            split_copy(design, host, depth, seed, host != primary)?;
        }
        Ok(())
    }

    /// The design segments as slot ranges, in slot order. A segment's
    /// films are consecutive slots, and a needle joins its host's segment,
    /// so they stay consecutive.
    pub(crate) fn segments(&self) -> Vec<std::ops::Range<usize>> {
        let mut out: Vec<std::ops::Range<usize>> = Vec::new();
        for (i, sl) in self.slots.iter().enumerate() {
            match out.last_mut() {
                Some(r) if self.slots[r.start].segment == sl.segment => r.end = i + 1,
                _ => out.push(i..i + 1),
            }
        }
        out
    }

    /// Every boundary locus of the design, segment by segment, top to
    /// bottom in forward order.
    pub(crate) fn gaps(&self) -> Vec<Gap> {
        self.segments()
            .iter()
            .enumerate()
            .flat_map(|(seg, r)| (0..=r.len()).map(move |k| Gap { seg, k }))
            .collect()
    }

    /// The design slots on either side of `gap` in forward order: the film
    /// above and the film below, `None` past the segment's edge.
    pub(crate) fn gap_neighbours(&self, gap: Gap) -> (Option<usize>, Option<usize>) {
        let r = &self.segments()[gap.seg];
        (
            (gap.k > 0).then(|| r.start + gap.k - 1),
            (gap.k < r.len()).then(|| r.start + gap.k),
        )
    }

    /// The edge lists `gap` answers to: the segment's outer list on its top
    /// edge, its inner list on its bottom edge.
    pub(crate) fn gap_edges(&self, gap: Gap) -> (Option<&EdgeMaterials>, bool, bool) {
        let r = &self.segments()[gap.seg];
        (
            self.segment_edges(&self.slots[r.start].segment),
            gap.k == 0,
            gap.k == r.len(),
        )
    }

    /// Where `gap` sits in each copy of its segment in environment `env`,
    /// whose stack has `spans` over `n_films` rows, in assembly order.
    ///
    /// A forward copy puts the new film before its slot `k` (after its
    /// last film when `k` is the segment's length). An inverted copy meets
    /// the films last first, so the same locus is after its slot `k`
    /// there, and the plane before the forward first film is below the
    /// copy. An inverted span that opens with a slice holds the slice its
    /// predecessor carved, and the new film takes it over: the slice stays
    /// against the film it belongs to.
    pub(crate) fn gap_places(
        &self,
        env: usize,
        spans: &[crate::structure::Span],
        n_films: usize,
        gap: Gap,
    ) -> Vec<GapPlace> {
        let r = self.segments()[gap.seg].clone();
        let n = r.len();
        (0..self.routing[r.start][env].len())
            .map(|c| {
                let span_of = |slot: usize| self.routing[slot][env][c];
                let inverted = self.is_inverted(env, span_of(r.start));
                let at = match (inverted, gap.k == n) {
                    (false, false) => span_of(r.start + gap.k),
                    (false, true) => span_of(r.end - 1) + 1,
                    (true, false) => span_of(r.start + gap.k) + 1,
                    (true, true) => span_of(r.end - 1),
                };
                let take_slice = inverted
                    && spans
                        .get(at)
                        .is_some_and(|sp| sp.slice && self.is_inverted(env, at));
                let row = match spans.get(at) {
                    Some(sp) if take_slice => sp.bulk_start,
                    Some(sp) => sp.start,
                    None => n_films,
                };
                GapPlace {
                    env,
                    at,
                    row,
                    inverted,
                    take_slice,
                }
            })
            .collect()
    }

    /// True when a film placed at `pl` leaves every film's own slice where
    /// it was. The one move that would not: a film turned into the end of an
    /// inverted run whose last entry has a slice it could not place before
    /// (the plane below it belonged to the next segment) and now could.
    pub(crate) fn gap_keeps_slices(&self, pl: &GapPlace) -> bool {
        let prev = pl.at.checked_sub(1);
        let stranded = prev.is_some_and(|p| {
            self.is_inverted(pl.env, p)
                && !self.is_inverted(pl.env, pl.at)
                && self
                    .slot_of_span(pl.env, p)
                    .is_some_and(|s| self.faces[s].slice.is_some())
        });
        !(pl.inverted && stranded)
    }

    /// The slot a film inserted at `gap` takes.
    fn gap_slot(&self, gap: Gap) -> usize {
        self.segments()[gap.seg].start + gap.k
    }

    /// Insert `film` at boundary locus `gap` in every copy in the
    /// pipeline's design stack (environment 0's), as
    /// [`Self::insert_boundary`] does in the templates. Call it first: it
    /// reads the routing as it stands. The first copy is the free one; the
    /// others are pinned and follow it.
    pub fn boundary_design(
        &self,
        design: &mut DesignStack,
        gap: Gap,
        film: &LayerSpec,
    ) -> Result<(), String> {
        let places = self.gap_places(0, design.spans(), design.films().len(), gap);
        insert_places(design, &places, film, true)
    }

    /// F2.3 for a boundary: insert `film` at locus `gap` in every copy in
    /// every environment's template, as a new design slot in the gap's
    /// segment. Routing shifts by one per copy placed before a span; the
    /// caller places the film in the design stack with
    /// [`Self::boundary_design`] first. Returns the new slot.
    ///
    /// A film turned into an inverted copy gets its face like any other:
    /// its own index and coherence, no roughness and no slice, so the
    /// plane above it in forward order is smooth wherever the film lands.
    pub fn insert_boundary(&mut self, gap: Gap, film: &LayerSpec) -> Result<usize, String> {
        let slot = self.gap_slot(gap);
        let segment = self.slots[self.segments()[gap.seg].start].segment.clone();
        let all: Vec<Vec<GapPlace>> = (0..self.names.len())
            .map(|e| {
                let st = &self.stacks[e];
                self.gap_places(e, st.spans(), st.films().len(), gap)
            })
            .collect();
        if let Some(pl) = all.iter().flatten().find(|pl| !self.gap_keeps_slices(pl)) {
            return Err(format!(
                "environment {:?}: a film at span {} would give the inverted copy \
                 above it an interface slice that has nowhere to go",
                self.names[pl.env], pl.at
            ));
        }
        let mut new_routing = Vec::with_capacity(self.names.len());
        for (e, places) in all.iter().enumerate() {
            insert_places(&mut self.stacks[e], places, film, e == 0)
                .map_err(|msg| format!("environment {:?}: {msg}", self.names[e]))?;
            // Every span at or after a placement moves down one for it.
            for per in self.routing.iter_mut() {
                for s in per[e].iter_mut() {
                    *s += places.iter().filter(|pl| pl.at <= *s).count();
                }
            }
            // In assembly order the c-th new span sits c further down, for
            // the c films placed above it.
            let placed: Vec<usize> = places.iter().enumerate().map(|(c, pl)| pl.at + c).collect();
            for (pl, &p) in places.iter().zip(&placed) {
                self.inverted[e].insert(p, pl.inverted);
            }
            new_routing.push(placed);
        }
        self.routing.insert(slot, new_routing);
        self.slots.insert(
            slot,
            DesignSlot {
                name: film.material.clone(),
                segment,
                intra: 0,
            },
        );
        self.faces.insert(
            slot,
            Face {
                nk: film.nk.clone(),
                coherent: film.coherent,
                rough_val: film.rough_val,
                rough_type: film.rough_type,
                slice: None,
            },
        );
        self.renumber_intra();
        Ok(slot)
    }

    /// `intra` as the position within the segment, re-derived after a
    /// structural move.
    fn renumber_intra(&mut self) {
        let mut next: HashMap<Arc<str>, usize> = HashMap::new();
        for sl in self.slots.iter_mut() {
            let n = next.entry(sl.segment.clone()).or_insert(0);
            sl.intra = *n;
            *n += 1;
        }
    }

    /// Copy every design film's primary thicknesses into its other copies
    /// in `design` (environment 0's stack), in place.
    ///
    /// The optimizer moves the primaries only; this is what makes the
    /// copies in the stack the pipeline carries (and returns) agree with
    /// them. A design without copies is left untouched.
    pub fn sync(&self, design: &mut DesignStack) -> Result<(), String> {
        if !self.has_copies() && !self.has_inverted() {
            return Ok(());
        }
        let src = design.clone();
        self.fill(0, design, &src)
    }

    /// Write the design's current state into `st`, environment `env`'s
    /// stack (in environment 0, a copy of the design itself).
    ///
    /// A forward copy of a forward primary takes its rows one for one;
    /// copying per ROW rather than per span total is what keeps a graded
    /// design film's profile intact. A forward copy of an inverted primary
    /// takes the forward bulk (inverted films are plain). Every inverted
    /// run is then re-emitted whole through the expander, from the films'
    /// faces at the primaries' thicknesses: an inverted entry carves its
    /// predecessor, so a run is one emission, not one per film.
    fn fill(&self, env: usize, st: &mut DesignStack, design: &DesignStack) -> Result<(), String> {
        for (slot, per) in self.routing.iter().enumerate() {
            let p = per[0][0];
            let src = design.spans()[p];
            let forward = !self.is_inverted(0, p);
            for (c, &dst_i) in per[env].iter().enumerate() {
                if (env == 0 && c == 0) || self.is_inverted(env, dst_i) {
                    continue;
                }
                let dst = st.spans()[dst_i];
                if forward {
                    for (a, b) in (src.start..src.end).zip(dst.start..dst.end) {
                        st.set_thickness(b, design.films()[a].d_nm)?;
                    }
                } else {
                    st.set_thickness(dst.bulk_start, self.forward_bulk(design, slot))?;
                }
            }
        }
        let flags = &self.inverted[env];
        let mut i = 0;
        while i < flags.len() {
            if !flags[i] {
                i += 1;
                continue;
            }
            let mut j = i;
            while j < flags.len() && flags[j] {
                j += 1;
            }
            self.reemit_run(env, st, design, i..j)?;
            i = j;
        }
        Ok(())
    }

    /// Re-emit the inverted run on spans `run` of `st` (environment `env`).
    /// A primary's bulk rows keep their thickness: they are the optimizer's
    /// parameters, and the emission would give them back only to rounding.
    fn reemit_run(
        &self,
        env: usize,
        st: &mut DesignStack,
        design: &DesignStack,
        run: std::ops::Range<usize>,
    ) -> Result<(), String> {
        let mut slots = Vec::with_capacity(run.len());
        let mut seq: Vec<(Layer, bool)> = Vec::with_capacity(run.len());
        for (k, span) in run.clone().enumerate() {
            let (slot, copy) = self.copy_of_span(env, span).ok_or_else(|| {
                format!(
                    "inverted span {span} of environment {:?} is not a design film",
                    self.names[env]
                )
            })?;
            let face = &self.faces[slot];
            let thickness = self.film_total(design, slot);
            seq.push((
                Layer {
                    material: k.to_string(),
                    thickness,
                    coherent: face.coherent,
                    rough_type: rough_kind(face.rough_type)?,
                    roughness: face.rough_val,
                    interface: face.slice.is_some(),
                    interface_thickness: face.slice.unwrap_or(0.0),
                    ..Layer::default()
                },
                true,
            ));
            slots.push(env == 0 && copy == 0);
        }
        let provider = RunProvider(
            run.clone()
                .map(|span| {
                    let slot = self.slot_of_span(env, span).expect("checked above");
                    &*self.faces[slot].nk
                })
                .collect(),
        );
        let (sa, spans) = crate::structure::expand(
            &seq,
            &provider,
            &self.wavelengths,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )?;
        let base = st.spans()[run.start].start;
        for (k, span) in run.enumerate() {
            let have = st.spans()[span];
            let want = spans[k];
            if have.end - have.start != want.end - want.start || have.slice != want.slice {
                return Err(format!(
                    "environment {:?}: an inverted copy no longer has the rows its \
                     re-emission gives ({} vs {} at span {span})",
                    self.names[env],
                    have.end - have.start,
                    want.end - want.start
                ));
            }
            for r in want.start..want.end {
                let row = base + r;
                let keep = slots[k] && row >= have.bulk_start;
                st.rewrite_row(
                    row,
                    (!keep).then_some(sa.thicknesses[r]),
                    sa.row(r),
                    !sa.incoherent[r],
                    sa.rough_types[r],
                    sa.rough_vals[r],
                );
            }
        }
        Ok(())
    }

    /// F2.2: the K stacks to solve, given the shared design's CURRENT
    /// thicknesses.
    ///
    /// `design` is environment 0's stack — the object the pipeline carries
    /// and the thickness optimizer writes into. Its copies are re-synced
    /// from their primaries ([`Self::sync`]); every other environment is
    /// its compiled template with each design film's primary row
    /// thicknesses copied into every copy by the routing. Copying per ROW
    /// rather than per span total is what keeps a graded design film's
    /// profile intact: the fractions never go through a divide.
    ///
    /// This is deliberately not a re-assembly. Under K > 1 the structure is
    /// frozen BETWEEN evaluations: the only structural move that exists is
    /// a needle insertion, and it does not happen here — it happens in
    /// [`Self::insert_seed`], which splits every environment's template at
    /// the same time, so the templates are correct again before the next
    /// expansion. Cleanup removal and inflate are refused outright by the
    /// driver (there is no inverse of `insert_seed`). So the only thing
    /// that moves between evals is a thickness, and re-running
    /// `from_design` K times per eval would pay the 126 µs assembly (R6)
    /// for nothing. The alignment check below is what makes the assumption
    /// fail loudly instead of quietly mis-routing if it is ever violated.
    ///
    /// Surroundings are never touched: a fixed segment's rows keep the
    /// template's thicknesses, which is §3.1's "fully fixed" made
    /// operational rather than promised.
    pub fn expand(&self, design: &DesignStack) -> Result<Vec<DesignStack>, String> {
        self.check_alignment(design)?;
        let mut out = Vec::with_capacity(self.names.len());
        let mut first = design.clone();
        self.fill(0, &mut first, design)?;
        out.push(first);
        for e in 1..self.names.len() {
            let mut st = self.stacks[e].clone();
            self.fill(e, &mut st, design)?;
            out.push(st);
        }
        Ok(out)
    }

    /// The cross-environment alignment assert: the live design stack still
    /// has the span layout the compile recorded, and every copy of every
    /// slot still names the same film in every environment.
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
                 with {} - under K > 1 the only structural move is \
                 `insert_seed`, which splits every environment's template \
                 at once, so a mismatch here is a driver or caller gap, \
                 not a pending feature",
                self.names[0],
                design.spans().len(),
                self.stacks[0].spans().len()
            ));
        }
        // Bulk rows, not span rows: an inverted copy's interface slices sit
        // in the neighbouring span, so only the bulk is the same everywhere.
        for (slot, ds) in self.slots.iter().enumerate() {
            let src = design.spans()[self.routing[slot][0][0]];
            let n_src = src.end - src.bulk_start;
            for e in 0..self.names.len() {
                for &c in &self.routing[slot][e] {
                    let dst = self.stacks[e].spans()[c];
                    if dst.end - dst.bulk_start != n_src {
                        return Err(format!(
                            "design film '{}' (segment '{}') is {} bulk row(s) in \
                             environment '{}' and {} in '{}' - the shared design \
                             must expand identically everywhere",
                            ds.name,
                            ds.segment,
                            n_src,
                            self.names[0],
                            dst.end - dst.bulk_start,
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
        }
        Ok(())
    }
}

/// Place `film` at every one of `places` in `stack` (one environment's
/// placements, in assembly order). With `pin`, every copy but the first is
/// pinned: it follows the first.
fn insert_places(
    stack: &mut DesignStack,
    places: &[GapPlace],
    film: &LayerSpec,
    pin: bool,
) -> Result<(), String> {
    // Last copy first: an insertion moves only the spans after it.
    for (c, pl) in places.iter().enumerate().rev() {
        stack.insert_span(pl.at, film.clone(), pl.take_slice)?;
        if pin && c > 0 {
            let sp = stack.spans()[pl.at];
            for r in sp.start..sp.end {
                stack.set_row_optimize(r, false);
            }
        }
    }
    Ok(())
}

/// Split one copy (span `host`) of a design film in `stack` with `seed`.
/// `pin` pins the new seed (`optimize` false): it is a copy in the design
/// object, not the primary, so it follows the primary's seed.
fn split_copy(
    stack: &mut DesignStack,
    host: usize,
    depth_into_layer_nm: f64,
    seed: &crate::smatrix::synthesis::structure::LayerSpec,
    pin: bool,
) -> Result<(), String> {
    let row = stack.spans()[host].bulk_start;
    stack.insert_needle_seed(row, depth_into_layer_nm, seed.clone())?;
    if pin {
        let sp = stack.spans()[host + 1];
        for r in sp.start..sp.end {
            stack.set_row_optimize(r, false);
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Compile
// ---------------------------------------------------------------------------

/// Auto-name for a surrounding film: `{env}.fixed[{seg}][{i}]` (§4.1).
///
/// Surroundings are auto-named only (§8 question 5, answered: auto-only
/// in v1). A user name there would address nothing — no buckets, no
/// needle sites — and inviting one invites the belief that it does.
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
/// Returns the rows in assembly order, the authoring indices each design
/// slot's copies landed at (in assembly order), and the NAMES of the rows
/// that came from a fixed segment.
///
/// That third value is M1's fix (review PB). `to_row` forces
/// `optimize`/`needle` false, but the assembler applies per-film
/// overrides afterwards, so it needs to know which rows those flags were
/// forced on. Reported from here, where the distinction is a fact of the
/// compile, rather than recovered downstream from the auto-naming
/// pattern — a name is presentation and a forced flag is a contract.
fn compile_env_rows(
    env: &EnvironmentCfg,
    design: &BTreeMap<String, DesignSegmentCfg>,
    slots: &[DesignSlot],
) -> Result<CompiledRows, String> {
    let mut rows: Vec<(String, LayerRow)> = Vec::new();
    let mut inverted: Vec<bool> = Vec::new();
    let mut fixed: HashSet<String> = HashSet::new();
    let mut refs: HashSet<&str> = HashSet::new();
    // Authoring indices per slot, one per copy: films only, in the order
    // the assembler will see them. Half-space rows do not advance it.
    let mut film_idx: usize = 0;
    let mut slot_at: Vec<Vec<usize>> = vec![Vec::new(); slots.len()];

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
                if seg.repeat.is_some() {
                    return Err(format!(
                        "environment {:?} segment {s}: \"repeat\" applies to a \
                         design reference - list a fixed segment again to use \
                         it twice",
                        env.name
                    ));
                }
                if seg.inverted {
                    return Err(format!(
                        "environment {:?} segment {s}: \"inverted\" applies to a \
                         design reference - list fixed layers in the order you \
                         want them",
                        env.name
                    ));
                }
                for (i, layer) in layers.iter().enumerate() {
                    let name = fixed_name(&env.name, s, i);
                    let row = layer.to_row(&format!("{name:?}"))?;
                    if row.layer_type == 1 {
                        // Films only: a half-space row placed in a fixed
                        // segment never reaches the per-film override.
                        film_idx += 1;
                        fixed.insert(name.clone());
                        inverted.push(false);
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
                let repeat = seg.repeat.unwrap_or(1);
                if repeat == 0 {
                    return Err(format!(
                        "environment {:?} segment {s}: \"repeat\" must be at \
                         least 1",
                        env.name
                    ));
                }
                refs.insert(id.as_str());
                if seg.inverted
                    && let Some(row) = seg_def
                        .layers
                        .iter()
                        .find(|r| r.inhomogen || r.gradient.is_some())
                {
                    return Err(format!(
                        "environment {:?} segment {s}: design film {:?} is graded - \
                         an inverted copy cannot carry a profile",
                        env.name, row.material_code
                    ));
                }
                let order: Vec<usize> = if seg.inverted {
                    (0..seg_def.layers.len()).rev().collect()
                } else {
                    (0..seg_def.layers.len()).collect()
                };
                for _ in 0..repeat {
                    for &i in &order {
                        let row = &seg_def.layers[i];
                        let name = row.material_code.clone();
                        inverted.push(seg.inverted);
                        let slot = slots
                            .iter()
                            .position(|sl| &*sl.segment == id.as_str() && sl.intra == i)
                            .expect("every design layer has a slot (compile_slots)");
                        slot_at[slot].push(film_idx);
                        film_idx += 1;
                        rows.push((name, row.clone()));
                    }
                }
            }
        }
    }

    // At least once (§4.2; exact-once until segments could repeat).
    // Omission is almost certainly a typo'd environment, and silent
    // non-contribution would lie about jointness. A segment used twice is
    // two copies of one design: one parameter per film, routed to both.
    for id in design.keys() {
        if !refs.contains(id.as_str()) {
            return Err(format!(
                "environment {:?} does not reference design segment {id:?} - \
                 every environment must reference each design segment at \
                 least once",
                env.name
            ));
        }
    }

    Ok(CompiledRows {
        rows,
        slot_at,
        fixed,
        inverted,
    })
}

/// One environment's segment list, flattened (see [`compile_env_rows`]).
struct CompiledRows {
    rows: Vec<(String, LayerRow)>,
    slot_at: Vec<Vec<usize>>,
    fixed: HashSet<String>,
    /// Per film (half-space rows excluded), in assembly order: the film is
    /// an entry of an inverted copy.
    inverted: Vec<bool>,
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
) -> Result<(CompiledEnvironments, MaterialPool, Vec<String>), String> {
    if req.environments.is_empty() {
        if !req.design.is_empty() {
            return Err(
                "\"design\" segments are defined but \"environments\" is empty \
                 - a design segment is only reachable through an environment \
                 that references it"
                    .to_string(),
            );
        }
        let (stack, pool, warnings) = super::design_config::build_design(req, wavelengths)?;
        let slots: Vec<DesignSlot> = stack
            .spans()
            .iter()
            .map(|sp| DesignSlot {
                name: Arc::from(&*stack.films()[sp.bulk_start].material),
                segment: Arc::from(DEFAULT_ENV),
                intra: sp.logical,
            })
            .collect();
        let routing = (0..slots.len()).map(|i| vec![vec![i]]).collect();
        // Forward throughout, so the faces are never re-emitted from; they
        // are read off the rows to keep the slot list and the faces aligned.
        let faces = stack
            .spans()
            .iter()
            .map(|sp| {
                let bulk: &LayerSpec = &stack.films()[sp.bulk_start];
                Face {
                    nk: bulk.nk.clone(),
                    coherent: bulk.coherent,
                    rough_val: bulk.rough_val,
                    rough_type: bulk.rough_type,
                    slice: sp.slice.then(|| stack.films()[sp.start].d_nm),
                }
            })
            .collect();
        let inverted = vec![vec![false; stack.spans().len()]];
        return Ok((
            CompiledEnvironments {
                names: vec![DEFAULT_ENV.to_string()],
                stacks: vec![stack],
                slots,
                routing,
                segment_edges: BTreeMap::new(),
                inverted,
                faces,
                wavelengths: wavelengths.to_vec(),
            },
            pool,
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
    let mut per_env_slot_at: Vec<Vec<Vec<usize>>> = Vec::with_capacity(req.environments.len());
    let mut per_env_fixed: Vec<HashSet<String>> = Vec::with_capacity(req.environments.len());
    let mut per_env_inverted: Vec<Vec<bool>> = Vec::with_capacity(req.environments.len());
    let mut nk = nk_base.clone();
    for env in &req.environments {
        let CompiledRows {
            rows,
            slot_at: at,
            fixed,
            inverted,
        } = compile_env_rows(env, &req.design, &slots)?;
        per_env_inverted.push(inverted);
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
        per_env_fixed.push(fixed);
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

    let pool = build_pool(req.materials.as_deref(), &nk)?;

    let mut stacks: Vec<DesignStack> = Vec::with_capacity(req.environments.len());
    let mut inverted: Vec<Vec<bool>> = Vec::with_capacity(req.environments.len());
    let mut env0_films: Vec<Layer> = Vec::new();
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
                fixed_films: &per_env_fixed[e],
                ambient_name: &req.ambient_name,
                substrate_name: &req.substrate_name,
                wavelengths,
                gate_label: "build_environments",
                // Names are unique per design film (`compile_slots`) and
                // per surrounding (`fixed_name`), so a repeated name here
                // is a copy of a repeated segment.
                copies: true,
            },
        )
        .map_err(|e| format!("environment {:?}: {e}", env.name))?;

        let background: HashSet<String> = films
            .iter()
            .filter(|l| (l.inhomogen || l.gradient.is_some()) && !l.optimize && !l.needle)
            .map(|l| l.material.clone())
            .collect();
        let orient = &per_env_inverted[e];
        if orient.len() != films.len() {
            return Err(format!(
                "environment {:?}: {} films assembled from {} compiled films",
                env.name,
                films.len(),
                orient.len()
            ));
        }
        let (stack, warns) = DesignStack::from_design_oriented(
            amb,
            sub,
            &films,
            orient,
            &nk,
            &groups,
            wavelengths,
            &background,
        )
        .map_err(|e| format!("environment {:?}: {e}", env.name))?;
        warnings.extend(warns.into_iter().map(|w| format!("{}: {w}", env.name)));
        inverted.push(stack.spans().iter().map(|sp| orient[sp.logical]).collect());
        stacks.push(stack);
        if e == 0 {
            env0_films = films;
        }
    }

    // design_slot -> (env, copies' spans). `Span::logical` is the
    // authoring film index, which is what `compile_env_rows` recorded per
    // copy.
    let mut routing: Vec<Vec<Vec<usize>>> = vec![Vec::with_capacity(stacks.len()); slots.len()];
    for (e, stack) in stacks.iter().enumerate() {
        for (slot, per) in routing.iter_mut().enumerate() {
            let copies = per_env_slot_at[e][slot]
                .iter()
                .map(|&want| {
                    stack
                        .spans()
                        .iter()
                        .position(|sp| sp.logical == want)
                        .ok_or_else(|| {
                            format!(
                                "environment {:?}: design film {:?} lost its span at \
                                 assembly (authoring index {want})",
                                req.environments[e].name, slots[slot].name
                            )
                        })
                })
                .collect::<Result<Vec<usize>, String>>()?;
            per.push(copies);
        }
    }
    // Each film's face from its authored layer through its group, by the
    // expander's own arithmetic: the same whichever way the primary faces.
    let default_group = Group::new("_default_");
    let faces: Vec<Face> = routing
        .iter()
        .enumerate()
        .map(|(slot, per)| {
            let layer = &env0_films[per_env_slot_at[0][slot][0]];
            let group = groups
                .get(layer.material.as_str())
                .unwrap_or(&default_group);
            let total = (layer.thickness * group.thick_factor + group.thick_summand).max(0.0);
            let bulk = stacks[0].spans()[per[0][0]].bulk_start;
            Face {
                nk: stacks[0].films()[bulk].nk.clone(),
                coherent: layer.coherent,
                rough_val: (layer.roughness + group.roughness_summand).max(0.0),
                rough_type: layer.rough_type as i32,
                slice: layer
                    .interface
                    .then(|| (layer.interface_thickness + group.interface_summand).min(total)),
            }
        })
        .collect();

    // One parameter per design film: in the design object (environment
    // 0's stack) every copy but the primary is pinned and follows it.
    for per in &routing {
        for &c in &per[0][1..] {
            let sp = stacks[0].spans()[c];
            for r in sp.start..sp.end {
                stacks[0].set_row_optimize(r, false);
            }
        }
    }

    // The default pool is the shared design's own layers: a surrounding's
    // material is never a needle unless the request names it.
    let pool = pool.unwrap_or_else(|| {
        let design_rows: HashSet<usize> = routing
            .iter()
            .map(|per| stacks[0].spans()[per[0][0]].bulk_start)
            .collect();
        pool_from_layers(&stacks[0], |row| design_rows.contains(&row))
    });
    let pool = match &req.insert_rules {
        Some(r) => pool.with_insert_rules(r)?,
        None => pool,
    };
    let segment_edges = segment_edge_lists(req, &pool)?;
    let names = req.environments.iter().map(|e| e.name.clone()).collect();
    let mut compiled = CompiledEnvironments {
        names,
        stacks,
        slots,
        routing,
        segment_edges,
        inverted,
        faces,
        wavelengths: wavelengths.to_vec(),
    };
    compiled.check_slices()?;
    if let Some(d) = req.seed_edges {
        compiled.seed_edges(d)?;
    }
    compiled.check_edges()?;
    Ok((compiled, pool, warnings))
}

impl CompiledEnvironments {
    /// The edge lists are per segment, checked on each segment's own first
    /// and last film in every copy in every environment. A joint run never
    /// removes a layer (removal is refused while K > 1, and the floor is an
    /// LM bound), an interior needle lands strictly inside its host, and a
    /// boundary needle on a segment edge is on that edge's list, so the
    /// segments' edges hold for the whole run.
    fn check_edges(&self) -> Result<(), String> {
        for (id, edges) in &self.segment_edges {
            let r = self
                .segments()
                .into_iter()
                .find(|r| self.slots[r.start].segment == *id)
                .expect("every segment with edge lists has films");
            for (e, stack) in self.stacks.iter().enumerate() {
                let spans = stack.spans();
                // Every copy has its own top and bottom: the n-th copy of
                // the first film pairs with the n-th copy of the last.
                let tops = &self.routing[r.start][e];
                let bottoms = &self.routing[r.end - 1][e];
                for (&t, &b) in tops.iter().zip(bottoms) {
                    edges
                        .check(stack.films(), Some(&spans[t]), Some(&spans[b]))
                        .map_err(|err| {
                            format!(
                                "environment {:?}, design segment {id:?}: {err}",
                                self.names[e]
                            )
                        })?;
                }
            }
        }
        Ok(())
    }

    /// `seed_edges` on the joint door: a segment whose first (last) film is
    /// off its outer (inner) list gets a `d_nm` film of the list's first
    /// material on that edge, in every copy, as a boundary needle would
    /// put it there. The film is a design film like any other (named
    /// `"<label>#<n>"`, optimized, a host), and the old edge film stays
    /// behind it. The edge is the segment's, wherever the segment sits.
    fn seed_edges(&mut self, d_nm: f64) -> Result<(), String> {
        if !(d_nm.is_finite() && d_nm > 0.0) {
            return Err(format!(
                "seed_edges: the seed film's thickness must be positive and finite (got {d_nm} nm)"
            ));
        }
        let ids: Vec<Arc<str>> = self.segment_edges.keys().cloned().collect();
        for id in ids {
            let edges = self.segment_edges[&id].clone();
            for (list, outer) in [(&edges.outer, true), (&edges.inner, false)] {
                let Some(first) = list.as_ref().and_then(|l| l.first()) else {
                    continue;
                };
                let segs = self.segments();
                let seg = segs
                    .iter()
                    .position(|r| self.slots[r.start].segment == id)
                    .expect("every segment with edge lists has films");
                let r = segs[seg].clone();
                let edge = if outer { r.start } else { r.end - 1 };
                let design = &self.stacks[0];
                let sp = design.spans()[self.routing[edge][0][0]];
                if edge_allows(list, design.films(), &sp) {
                    continue;
                }
                let film = new_needle(design.films(), &pool_entry(&first.0, first.1.clone()), d_nm);
                let gap = Gap {
                    seg,
                    k: if outer { 0 } else { r.len() },
                };
                self.insert_boundary(gap, &film)
                    .map_err(|e| format!("seed_edges, design segment {id:?}: {e}"))?;
            }
        }
        Ok(())
    }

    /// Every copy of a design film must agree on whether the film's own
    /// interface slice is emitted, or the copies are different films: a
    /// forward copy at the top of its stack has no plane above it to put
    /// the slice on, and neither has the last film of an inverted copy
    /// (that plane faces the next segment).
    fn check_slices(&self) -> Result<(), String> {
        for (slot, per) in self.routing.iter().enumerate() {
            if self.faces[slot].slice.is_none() {
                continue;
            }
            let want = self.emits_slice(0, per[0][0], slot);
            for (e, copies) in per.iter().enumerate() {
                if copies.iter().any(|&c| self.emits_slice(e, c, slot) != want) {
                    let (with, without) = if want {
                        (&self.names[0], &self.names[e])
                    } else {
                        (&self.names[e], &self.names[0])
                    };
                    return Err(format!(
                        "design film {:?} (segment {:?}) has an interface slice that \
                         environment {with:?} puts on a plane and environment \
                         {without:?} cannot: there the film is the first of its \
                         stack, or the last of an inverted copy, and that plane \
                         belongs to its neighbour. Give the film no interface, or \
                         place every copy where the plane is the film's own",
                        self.slots[slot].name, self.slots[slot].segment
                    ));
                }
            }
        }
        Ok(())
    }
}

/// Each design segment's edge lists, resolved against the pool.
///
/// The request-level `outer`/`inner` predate per-segment lists. With one
/// design segment they are that segment's lists; with several they would
/// not say which segment they mean, so they refuse, and so does a request
/// that sets both a request-level list and the segment's own.
fn segment_edge_lists(
    req: &DesignRequest,
    pool: &MaterialPool,
) -> Result<BTreeMap<Arc<str>, EdgeMaterials>, String> {
    let global = req.outer.is_some() || req.inner.is_some();
    if global && req.design.len() > 1 {
        return Err(format!(
            "outer/inner: with {} design segments the request-level edge \
             lists are ambiguous - set \"outer\"/\"inner\" on each design \
             segment instead",
            req.design.len()
        ));
    }
    let mut out = BTreeMap::new();
    for (id, seg) in &req.design {
        let own = seg.outer.is_some() || seg.inner.is_some();
        if global && own {
            return Err(format!(
                "outer/inner: design segment {id:?} sets its own edge lists \
                 and the request sets them too - give them in one place"
            ));
        }
        let (outer, inner) = if own {
            (seg.outer.as_deref(), seg.inner.as_deref())
        } else {
            (req.outer.as_deref(), req.inner.as_deref())
        };
        if outer.is_none() && inner.is_none() {
            continue;
        }
        let edges = pool
            .edge_materials(outer, inner)
            .map_err(|e| format!("design segment {id:?}: {e}"))?;
        out.insert(Arc::from(id.as_str()), edges);
    }
    Ok(out)
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
    #[allow(unused_imports)]
    use crate::smatrix::angle::{Radians, deg, degs};
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
            materials: Some(vec![]), // no needles unless a test asks
            insert_rules: None,
            outer: None,
            inner: None,
            seed_edges: None,
            film_flags: BTreeMap::new(),
            per_film_flags: BTreeMap::new(),
            ambient_name: "air".to_string(),
            substrate_name: "sub".to_string(),
            design: design
                .into_iter()
                .map(|(id, layers)| {
                    (
                        id.to_string(),
                        DesignSegmentCfg {
                            layers,
                            outer: None,
                            inner: None,
                        },
                    )
                })
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
            materials: Some(vec![]), // no needles unless a test asks
            insert_rules: None,
            outer: None,
            inner: None,
            seed_edges: None,
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
            repeat: None,
            inverted: false,
        }
    }

    fn fixed_seg(rows: Vec<FixedLayerRow>) -> EnvSegmentCfg {
        EnvSegmentCfg {
            layers: Some(rows),
            design: None,
            repeat: None,
            inverted: false,
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
            angle: Radians::ZERO,
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

    /// Under K > 1 the only structural move is `insert_seed`, which keeps
    /// the templates in step. If a drift ever reaches `expand` anyway, it
    /// must name the drift rather than route thicknesses into the wrong
    /// spans.
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
        assert!(e.contains("insert_seed"), "{e}");
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
                &[Radians::ZERO],
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

    /// `env_spec`, plus the roster a compiled spec would carry (M7).
    fn env_spec_named(names: &[&str]) -> MeritSpec {
        let mut spec = env_spec(names.len());
        let owned: Vec<String> = names.iter().map(|s| s.to_string()).collect();
        spec.set_env_roster(&owned).unwrap();
        spec
    }

    /// M7 (review PB): the same names in a different order are a
    /// different binding, and the count cannot see it.
    ///
    /// A demand's `env_idx` is a POSITION in the roster, resolved when
    /// the spec was compiled and gone afterwards. Before the spec kept
    /// its names, `["lam", "bare"]` against a design compiling to
    /// `["bare", "lam"]` matched on count and every demand scored
    /// against the wrong surroundings — measured on the 0.7.5 wheel as a
    /// completed run reporting 8523.676912295265 where the correct
    /// pairing gives 19876.13664037235.
    #[test]
    fn the_run_refuses_a_permuted_roster_naming_both_orders() {
        use crate::smatrix::synthesis::driver::run_environments;
        let (envs, cmap, _) = build_environments(&two_env_req(2), &WL).unwrap();
        let e = run_environments(
            envs,
            cmap,
            &WL,
            &[Radians::ZERO],
            &env_spec_named(&["e1", "e0"]),
            thickness_only_cfg(),
            Default::default(),
            Default::default(),
            |_, _| Ok(()),
        )
        .unwrap_err();
        // Both lists, in their own order, because "they differ" without
        // saying how is a message that sends the reader back to the code.
        assert!(e.contains("e1, e0"), "{e}");
        assert!(e.contains("e0, e1"), "{e}");
        assert!(e.contains("POSITION"), "{e}");
        assert!(e.is_ascii(), "{e}");
    }

    /// M7's other half: a matching roster changes NOTHING.
    ///
    /// The check is additive or it is a regression — recording names
    /// must not move a single bit of the answer, and a spec that records
    /// none must keep running on the count check alone. Both spellings
    /// of the same request, compared on the bit.
    #[test]
    fn a_matching_roster_is_bit_identical_to_an_unrostered_spec() {
        use crate::smatrix::synthesis::driver::run_environments;
        let run = |spec: &MeritSpec| {
            let (envs, cmap, _) = build_environments(&two_env_req(2), &WL).unwrap();
            run_environments(
                envs,
                cmap,
                &WL,
                &[Radians::ZERO],
                spec,
                thickness_only_cfg(),
                Default::default(),
                Default::default(),
                |_, _| Ok(()),
            )
            .unwrap()
            .0
        };
        let named = run(&env_spec_named(&["e0", "e1"]));
        let plain = run(&env_spec(2));
        assert_eq!(
            named.final_mf.to_bits(),
            plain.final_mf.to_bits(),
            "{} vs {}",
            named.final_mf,
            plain.final_mf
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
            &[Radians::ZERO],
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

    /// `r` with the default needle pool: the design's own L and H.
    fn with_pool(mut r: DesignRequest) -> DesignRequest {
        r.materials = None;
        r
    }

    /// `two_env_req` plus the L, H pool the needle sweep needs.
    fn needle_req(k: usize) -> DesignRequest {
        with_pool(two_env_req(k))
    }

    /// Two environments that genuinely DIFFER: the same coat bare, and
    /// behind a laminate. The design rows sit at different depths, so a
    /// locus is only shared through the routing table.
    fn two_cover_req() -> DesignRequest {
        with_pool(req(
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
    ///
    /// The LM gets room to converge. At the default 200 iterations both
    /// runs stop about 1 nm short of the minimum, partway down a slow
    /// valley, and where each stops depends on rounding: an optimized
    /// build put K = 1 and K = 2 0.09 nm apart (74.496 vs 74.404 nm). Both
    /// reach the same minimum (73.445 nm, 1e-6 apart) when they converge.
    #[test]
    fn a_joint_needle_over_identical_environments_lands_where_one_does() {
        use crate::smatrix::synthesis::driver::run_environments;
        use crate::smatrix::synthesis::thick_opt::LmConfig;
        let lm = LmConfig {
            max_iterations: 2000,
            ..LmConfig::default()
        };
        let run = |k: usize| {
            let (envs, cmap, _) = build_environments(&needle_req(k), &WL).unwrap();
            run_environments(
                envs,
                cmap,
                &WL,
                &[Radians::ZERO],
                &env_spec(k),
                one_needle_cfg(),
                Default::default(),
                lm.clone(),
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
            // Relative. The two runs insert the same seed at the same place
            // and then descend differently - K = 2 differences the
            // Jacobian and carries twice the residuals, so the LM that
            // follows the insertion rounds its way to a slightly different
            // point on the same minimum (F2.2's CORRECTION 8, one insertion
            // further along). Converged, that is about 1e-6; 1e-5 leaves
            // room for other platforms' rounding. A real disagreement about
            // WHERE to put the needle shows up in the material sequence
            // above.
            assert!(
                (a.d_nm - b.d_nm).abs() <= 1e-5 * a.d_nm.max(1.0),
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
    fn joint_rules_check_the_design_edges_and_narrow_the_pool() {
        // The coat is L over H, behind nothing in one environment and a
        // laminate in the other: the edges are the DESIGN's, not env 0's
        // outermost row.
        let mut r = two_cover_req();
        r.outer = Some(vec!["L".into()]);
        r.inner = Some(vec!["H".into()]);
        build_environments(&r, &WL).unwrap();

        r.outer = Some(vec!["H".into()]);
        let e = build_environments(&r, &WL).err().unwrap();
        assert!(e.contains("outer") && e.contains("\"L\""), "{e}");

        r.outer = None;
        r.insert_rules = Some(
            [("L".to_string(), vec!["Z".to_string()])]
                .into_iter()
                .collect(),
        );
        let e = build_environments(&r, &WL).err().unwrap();
        assert!(e.contains("\"Z\" is not a pool material"), "{e}");
    }

    /// Step 5a: a segment used twice is one design in two places. Every
    /// copy routes to the same parameter; in the design object the primary
    /// is free and the other copy pinned; `expand` writes the primary into
    /// every copy in every environment.
    #[test]
    fn a_repeated_segment_routes_one_film_to_every_copy() {
        use crate::smatrix::synthesis::context::DesignContext;
        let mut twice = design_seg("coat");
        twice.repeat = Some(2);
        let r = req(
            vec![("coat", vec![film("L", 100.0), film("H", 80.0)])],
            vec![
                EnvironmentCfg {
                    name: "rep".to_string(),
                    stack: vec![twice.clone()],
                },
                EnvironmentCfg {
                    name: "refs".to_string(),
                    stack: vec![
                        fixed_seg(vec![fixed("G", 300.0, 1)]),
                        design_seg("coat"),
                        fixed_seg(vec![fixed("M", 40.0, 1)]),
                        design_seg("coat"),
                    ],
                },
            ],
        );
        let (envs, _, _) = build_environments(&r, &WL).unwrap();
        assert_eq!(
            envs.slots().len(),
            2,
            "one parameter per film, not per copy"
        );
        assert_eq!(envs.spans_of(0, 0), &[0, 2]);
        assert_eq!(envs.spans_of(1, 0), &[1, 3]);
        assert_eq!(envs.spans_of(0, 1), &[1, 4]);
        assert_eq!(envs.n_placements(0), 4);
        assert_eq!(envs.copy_of_span(1, 4), Some((0, 1)));
        assert_eq!(envs.copy_of_span(1, 3), None, "a surrounding");
        assert!(envs.has_copies() && envs.is_joint());

        let d = envs.stacks()[0].films();
        assert!(d[0].optimize && d[1].optimize, "the primaries are free");
        assert!(!d[2].optimize && !d[3].optimize, "the copies follow");
        assert!(d.iter().all(|f| f.needle), "every copy is scanned");

        // One step on the primary reaches every copy everywhere.
        let mut design = envs.stacks()[0].clone();
        design.set_thickness(0, 120.0).unwrap();
        let out = envs.expand(&design).unwrap();
        assert_eq!(out[0].films()[2].d_nm, 120.0);
        assert_eq!(out[1].films()[1].d_nm, 120.0);
        assert_eq!(out[1].films()[4].d_nm, 120.0);
        envs.sync(&mut design).unwrap();
        assert_eq!(design.films()[2].d_nm, 120.0);

        // One environment with copies is still joint, and solves as the
        // same rows solved flat, bit for bit.
        let one = req(
            vec![("coat", vec![film("L", 100.0), film("H", 80.0)])],
            vec![EnvironmentCfg {
                name: "rep".to_string(),
                stack: vec![twice],
            }],
        );
        let (one, _, _) = build_environments(&one, &WL).unwrap();
        assert!(one.is_single() && one.is_joint());
        assert!(
            !build_environments(&two_env_req(1), &WL)
                .unwrap()
                .0
                .is_joint()
        );
        let stack = one.stacks()[0].clone();
        assert_eq!(stack.films().len(), 4);
        let a = ctx(env_spec(1), None).evaluate_merit(&stack).unwrap();
        let b = ctx(env_spec(1), Some(Arc::new(one)))
            .evaluate_merit(&stack)
            .unwrap();
        assert_eq!(a.to_bits(), b.to_bits(), "{a} vs {b}");
    }

    /// Step 5a at the run level: a needle in a repeated segment lands in
    /// every copy, and the optimizer keeps the copies the same design.
    #[test]
    fn a_needle_lands_in_every_copy_of_a_repeated_segment() {
        use crate::smatrix::synthesis::driver::run_environments;
        let mut twice = design_seg("coat");
        twice.repeat = Some(2);
        let r = with_pool(req(
            vec![("coat", vec![film("L", 100.0), film("H", 80.0)])],
            vec![EnvironmentCfg {
                name: "rep".to_string(),
                stack: vec![twice],
            }],
        ));
        let (envs, cmap, _) = build_environments(&r, &WL).unwrap();
        let (_, s, _) = run_environments(
            envs,
            cmap,
            &WL,
            &[Radians::ZERO],
            &env_spec(1),
            one_needle_cfg(),
            Default::default(),
            Default::default(),
            |_, _| Ok(()),
        )
        .unwrap();
        let f = s.films();
        assert!(
            f.len() > 4,
            "no needle was inserted - the test tests nothing"
        );
        assert_eq!(
            f.len() % 2,
            0,
            "{} films: the copies split differently",
            f.len()
        );
        let (a, b) = f.split_at(f.len() / 2);
        for (i, (x, y)) in a.iter().zip(b).enumerate() {
            assert_eq!(x.material, y.material, "film {i}");
            assert_eq!(x.d_nm.to_bits(), y.d_nm.to_bits(), "film {i}");
            assert!(
                x.optimize && !y.optimize,
                "film {i}: only the primary is free"
            );
        }
    }

    // --- Step 5b: inverted copies --------------------------------------

    fn inverted_seg(id: &str) -> EnvSegmentCfg {
        EnvSegmentCfg {
            inverted: true,
            ..design_seg(id)
        }
    }

    /// L over H over M, every plane carrying something: an interface
    /// slice on L and H, roughness on all three.
    fn faced_coat() -> Vec<LayerRow> {
        let mut l = film("L", 100.0);
        l.interface = true;
        l.interface_thickness_nm = 3.0;
        l.roughness_nm = 1.5;
        l.rough_type = 1;
        let mut h = film("H", 80.0);
        h.interface = true;
        h.interface_thickness_nm = 2.0;
        h.roughness_nm = 0.8;
        h.rough_type = 1;
        let mut m = film("M", 60.0);
        m.roughness_nm = 0.4;
        m.rough_type = 1;
        vec![l, h, m]
    }

    /// `faced_coat` with the films' own thicknesses `d`, turned over.
    fn turned_rows(d: [f64; 3]) -> SolverArrays {
        let c = faced_coat();
        let films: Vec<(&str, f64, Option<f64>, f64)> = c
            .iter()
            .zip(d)
            .map(|(r, d)| {
                (
                    r.material_code.as_str(),
                    d,
                    r.interface.then_some(r.interface_thickness_nm),
                    r.roughness_nm,
                )
            })
            .collect();
        turned(&films)
    }

    /// Forward films `(material, thickness, slice, roughness)` turned over
    /// and expanded the Architect's way (an inverted run, linear
    /// roughness).
    fn turned(films: &[(&str, f64, Option<f64>, f64)]) -> SolverArrays {
        let nk = eval_library(&lib(), &WL).unwrap();
        let provider = crate::structure::DictProvider::with_grid(
            nk.iter()
                .map(|(k, v)| (k.to_string(), crate::structure::Entry::Array(v.clone())))
                .collect(),
            WL.to_vec(),
        )
        .unwrap();
        let mut seq: Vec<(Layer, bool)> = films
            .iter()
            .map(|&(m, d, slice, rough)| {
                (
                    Layer {
                        material: m.to_string(),
                        thickness: d,
                        roughness: rough,
                        rough_type: if rough > 0.0 {
                            RoughnessType::Linear
                        } else {
                            RoughnessType::None
                        },
                        interface: slice.is_some(),
                        interface_thickness: slice.unwrap_or(0.0),
                        ..Layer::default()
                    },
                    true,
                )
            })
            .collect();
        seq.reverse();
        crate::structure::expand(
            &seq,
            &provider,
            &WL,
            &HashMap::new(),
            ExpandOptions::deterministic(),
        )
        .unwrap()
        .0
    }

    use crate::structure::SolverArrays;

    fn assert_rows_are(st: &DesignStack, want: &SolverArrays, tol: f64) {
        let f = st.films();
        assert_eq!(f.len(), want.n_rows(), "row count");
        for (r, row) in f.iter().enumerate() {
            assert!(
                (row.d_nm - want.thicknesses[r]).abs() <= tol,
                "row {r}: {} vs {} nm",
                row.d_nm,
                want.thicknesses[r]
            );
            assert_eq!(&*row.nk, want.row(r), "row {r}: index");
            assert_eq!(
                row.rough_type, want.rough_types[r],
                "row {r}: roughness type"
            );
            assert_eq!(row.rough_val, want.rough_vals[r], "row {r}: roughness");
        }
    }

    /// An inverted copy is the Architect's expansion of the segment turned
    /// over, bit for bit: each slice and roughness moves with its plane
    /// rather than with its row.
    #[test]
    fn an_inverted_copy_is_the_segment_turned_over() {
        let r = req(
            vec![("coat", faced_coat())],
            vec![
                EnvironmentCfg {
                    name: "front".to_string(),
                    stack: vec![design_seg("coat")],
                },
                EnvironmentCfg {
                    name: "back".to_string(),
                    stack: vec![inverted_seg("coat")],
                },
            ],
        );
        let (envs, _, _) = build_environments(&r, &WL).unwrap();
        let want = turned_rows([100.0, 80.0, 60.0]);
        assert_rows_are(&envs.stacks()[1], &want, 0.0);

        // What moved, in words. Turned over, M comes first and meets the
        // light with a clean plane; H's bulk carries M's roughness; L's
        // span opens with H's slice, carved from H. L's own slice is on
        // neither stack: L is the first film forward and the last turned,
        // and both planes are clean.
        let back = &envs.stacks()[1];
        let sp = back.spans();
        assert_eq!(sp.len(), 3);
        assert!(!sp[0].slice && back.films()[0].rough_val == 0.0);
        assert!(!sp[1].slice && back.films()[sp[1].bulk_start].rough_val == 0.4);
        assert!(sp[2].slice && back.films()[sp[2].start].d_nm == 2.0);
        assert_eq!(back.films()[sp[2].bulk_start].d_nm, 100.0);
        assert_eq!(back.films()[sp[2].bulk_start].rough_val, 0.8);
        assert!(envs.is_inverted(1, 0) && !envs.is_inverted(0, 0));
        assert!(envs.is_joint());

        // A thickness step on the primary re-emits the turned copy as a
        // fresh compile at the new thickness would.
        let mut design = envs.stacks()[0].clone();
        let l_bulk = design.spans()[0].bulk_start;
        design.set_thickness(l_bulk, 110.0).unwrap();
        let out = envs.expand(&design).unwrap();
        assert_rows_are(&out[1], &turned_rows([110.0, 80.0, 60.0]), 1e-12);
    }

    /// Lossless and between equal media, a stack reflects the same from
    /// either side: the coat turned over, behind the same glass both ways,
    /// has the forward coat's merit.
    #[test]
    fn a_turned_coat_reflects_as_the_coat_does() {
        let glass = || fixed_seg(vec![fixed("G", 0.0, 0)]);
        let base = || fixed_seg(vec![fixed("G", 0.0, 2)]);
        let r = req(
            vec![(
                "coat",
                vec![film("L", 100.0), film("H", 80.0), film("M", 60.0)],
            )],
            vec![
                EnvironmentCfg {
                    name: "front".to_string(),
                    stack: vec![glass(), design_seg("coat"), base()],
                },
                EnvironmentCfg {
                    name: "back".to_string(),
                    stack: vec![glass(), inverted_seg("coat"), base()],
                },
            ],
        );
        let (envs, _, _) = build_environments(&r, &WL).unwrap();
        let m = env_merits(&envs, &envs.stacks()[0]);
        assert!(m[0] > 0.0);
        assert!(
            (m[0] - m[1]).abs() <= 1e-12 * m[0],
            "front {} vs back {}",
            m[0],
            m[1]
        );
    }

    /// A needle splits the turned copy at the mirror depth, and the pieces
    /// take their planes with them: after inserting M 30 nm into H's bulk,
    /// the turned copy is the turned expansion of the split coat.
    #[test]
    fn a_needle_splits_an_inverted_copy_at_the_mirror_depth() {
        let r = req(
            vec![("coat", faced_coat())],
            vec![
                EnvironmentCfg {
                    name: "front".to_string(),
                    stack: vec![design_seg("coat")],
                },
                EnvironmentCfg {
                    name: "back".to_string(),
                    stack: vec![inverted_seg("coat")],
                },
            ],
        );
        let (mut envs, _, _) = build_environments(&r, &WL).unwrap();
        let mut design = envs.stacks()[0].clone();
        let seed = LayerSpec::constant("M", 1.8, 0.0, 10.0, WL.len());
        envs.split_design(&mut design, 1, 30.0, &seed).unwrap();
        envs.insert_seed(1, 30.0, &seed).unwrap();

        // Forward: H keeps its slice and 30 nm of bulk above the seed, and
        // the 48 nm below it is a clone of H's bulk, roughness included.
        let d: Vec<f64> = design.films().iter().map(|f| f.d_nm).collect();
        assert_eq!(d, [100.0, 2.0, 30.0, 10.0, 48.0, 60.0]);
        assert_eq!(envs.slots().len(), 5);

        let (h, seed_rough) = (0.8, 0.0);
        let want = turned(&[
            ("L", 100.0, Some(3.0), 1.5),
            ("H", 32.0, Some(2.0), h),
            ("M", 10.0, None, seed_rough),
            ("H", 48.0, None, h),
            ("M", 60.0, None, 0.4),
        ]);
        let out = envs.expand(&design).unwrap();
        assert_rows_are(&out[1], &want, 1e-12);
        // Routing: in the turned copy the bottom portion comes first.
        assert_eq!(envs.spans_of(3, 1), &[1]);
        assert_eq!(envs.spans_of(2, 1), &[2]);
        assert_eq!(envs.spans_of(1, 1), &[3]);
    }

    /// One environment holding the coat both ways round, either way round
    /// first: a needle run keeps the turned copy the forward copy's rows
    /// in reverse, and only the primary (the first copy) is free.
    #[test]
    fn a_needle_run_keeps_the_turned_copy_turned() {
        use crate::smatrix::synthesis::driver::run_environments;
        for turned_first in [false, true] {
            let mut stack = vec![
                fixed_seg(vec![fixed("G", 0.0, 0)]),
                design_seg("coat"),
                fixed_seg(vec![fixed("G", 40.0, 1)]),
                inverted_seg("coat"),
                fixed_seg(vec![fixed("G", 0.0, 2)]),
            ];
            if turned_first {
                stack.swap(1, 3);
            }
            let r = with_pool(req(
                vec![("coat", vec![film("L", 100.0), film("H", 80.0)])],
                vec![EnvironmentCfg {
                    name: "both".to_string(),
                    stack,
                }],
            ));
            let (envs, cmap, _) = build_environments(&r, &WL).unwrap();
            assert!(envs.is_joint());
            let (_, s, _) = run_environments(
                envs,
                cmap,
                &WL,
                &[Radians::ZERO],
                &env_spec(1),
                one_needle_cfg(),
                Default::default(),
                Default::default(),
                |_, _| Ok(()),
            )
            .unwrap();
            let f = s.films();
            assert!(
                f.len() > 5,
                "no needle was inserted - the test tests nothing"
            );
            let n = (f.len() - 1) / 2;
            let (first, rest) = f.split_at(n);
            let (glass, second) = rest.split_at(1);
            assert_eq!(glass[0].d_nm, 40.0, "the spacer stayed put");
            for (i, (x, y)) in first.iter().zip(second.iter().rev()).enumerate() {
                assert_eq!(x.material, y.material, "film {i} ({turned_first})");
                assert_eq!(
                    x.d_nm.to_bits(),
                    y.d_nm.to_bits(),
                    "film {i} ({turned_first})"
                );
                assert!(
                    x.optimize && !y.optimize,
                    "film {i}: only the primary is free"
                );
            }
        }
    }

    // --- Step 6: boundary needles -------------------------------------

    /// The same rows and the same span shape, to `tol` in thickness.
    fn assert_same_rows(got: &DesignStack, want: &DesignStack, tol: f64, what: &str) {
        let shape = |s: &DesignStack| -> Vec<(usize, usize, bool, usize)> {
            s.spans()
                .iter()
                .map(|sp| (sp.start, sp.end, sp.slice, sp.bulk_start))
                .collect()
        };
        assert_eq!(shape(got), shape(want), "{what}: spans");
        for (r, (a, b)) in got.films().iter().zip(want.films()).enumerate() {
            assert!(
                (a.d_nm - b.d_nm).abs() <= tol,
                "{what}, row {r}: {} vs {} nm",
                a.d_nm,
                b.d_nm
            );
            assert_eq!(a.nk, b.nk, "{what}, row {r}: index");
            assert_eq!(
                a.rough_type, b.rough_type,
                "{what}, row {r}: roughness type"
            );
            assert_eq!(a.rough_val, b.rough_val, "{what}, row {r}: roughness");
        }
    }

    /// The faced coat forward, and turned over twice behind a spacer.
    fn faced_envs(coat: Vec<LayerRow>) -> DesignRequest {
        req(
            vec![("coat", coat)],
            vec![
                EnvironmentCfg {
                    name: "front".to_string(),
                    stack: vec![design_seg("coat")],
                },
                EnvironmentCfg {
                    name: "back".to_string(),
                    stack: vec![
                        inverted_seg("coat"),
                        fixed_seg(vec![fixed("G", 40.0, 1)]),
                        inverted_seg("coat"),
                    ],
                },
            ],
        )
    }

    /// A film placed on a boundary locus is the film authored there: every
    /// copy, forward and turned, has the rows a fresh compile of the coat
    /// with the film in it gives, slices re-mixed and moved with their
    /// planes. The plane above L in forward order is the one locus that
    /// refuses: L's slice would gain a plane in the turned copies that the
    /// forward copy does not give it.
    #[test]
    fn a_boundary_film_is_the_film_authored_there() {
        let film_g = LayerSpec::constant("G", 1.52, 0.0, 10.0, WL.len());
        for k in 1..=3 {
            let (mut envs, _, _) = build_environments(&faced_envs(faced_coat()), &WL).unwrap();
            let gap = Gap { seg: 0, k };
            let mut design = envs.stacks()[0].clone();
            envs.boundary_design(&mut design, gap, &film_g).unwrap();
            let slot = envs.insert_boundary(gap, &film_g).unwrap();
            assert_eq!(slot, k);
            assert_eq!(envs.slots().len(), 4);
            assert_eq!(envs.n_placements(slot), 3);
            let got = envs.expand(&design).unwrap();

            let mut coat = faced_coat();
            coat.insert(k, film("G", 10.0));
            let (fresh, _, _) = build_environments(&faced_envs(coat), &WL).unwrap();
            let want = fresh.expand(&fresh.stacks()[0]).unwrap();
            for (e, (g, w)) in got.iter().zip(&want).enumerate() {
                assert_same_rows(g, w, 1e-12, &format!("k = {k}, environment {e}"));
            }
            assert_same_rows(&design, &want[0], 0.0, &format!("k = {k}, the design"));
            // Every copy but the primary follows it.
            let rows: Vec<bool> = envs
                .spans_of(slot, 0)
                .iter()
                .map(|&s| design.films()[design.spans()[s].bulk_start].optimize)
                .collect();
            assert_eq!(rows, [true]);
        }

        let (mut envs, _, _) = build_environments(&faced_envs(faced_coat()), &WL).unwrap();
        let e = envs
            .insert_boundary(Gap { seg: 0, k: 0 }, &film_g)
            .unwrap_err();
        assert!(e.contains("nowhere to go"), "{e}");
    }

    /// The cases table at a segment edge. The neighbour outside the design
    /// brings no rules and excludes a needle only when it has the needle's
    /// index everywhere; inside, the design film's index and rules apply.
    #[test]
    fn a_segment_edge_takes_a_needle_that_matters_somewhere() {
        let edge = |covers: [Option<&str>; 2]| -> Vec<String> {
            let envs: Vec<EnvironmentCfg> = covers
                .iter()
                .enumerate()
                .map(|(i, cover)| EnvironmentCfg {
                    name: format!("e{i}"),
                    stack: match cover {
                        Some(m) => vec![fixed_seg(vec![fixed(m, 300.0, 1)]), design_seg("coat")],
                        None => vec![design_seg("coat")],
                    },
                })
                .collect();
            let mut r = req(
                vec![("coat", vec![film("L", 100.0), film("H", 80.0)])],
                envs,
            );
            r.materials = Some(["L", "H", "M", "G"].map(String::from).to_vec());
            let (envs, pool, _) = build_environments(&r, &WL).unwrap();
            let stacks = envs.expand(&envs.stacks()[0]).unwrap();
            let gap = Gap { seg: 0, k: 0 };
            let places: Vec<Vec<GapPlace>> = stacks
                .iter()
                .enumerate()
                .map(|(e, st)| envs.gap_places(e, st.spans(), st.films().len(), gap))
                .collect();
            (0..pool.len())
                .filter(|&m| pool.gap_accepts(&envs, &stacks, gap, &places, m).unwrap())
                .map(|m| pool[m].material.to_string())
                .collect()
        };
        // L is the inside neighbour's index everywhere.
        assert_eq!(edge([None, None]), ["H", "M", "G"]);
        // M under an M cover does nothing there, but it is a film in air.
        assert_eq!(edge([None, Some("M")]), ["H", "M", "G"]);
        assert_eq!(edge([Some("M"), Some("M")]), ["H", "G"]);
        assert_eq!(edge([Some("M"), Some("G")]), ["H", "M", "G"]);
    }

    /// A joint boundary run over two identical environments lands where
    /// the single-environment run does: with no interior site, the needle
    /// goes on a boundary of the coat, the same one, at the same thickness.
    ///
    /// The host is held at 100 nm. A film the LM has converged gives its
    /// boundaries no first-order gain (measured: the merit change of an L
    /// film on top of a converged H film falls as the square of its
    /// thickness), so no boundary would pass the convergence gate.
    #[test]
    fn a_joint_boundary_needle_over_identical_environments_lands_where_one_does() {
        use crate::smatrix::synthesis::cycle::NeedleCycleConfig;
        use crate::smatrix::synthesis::driver::run_environments;
        use crate::smatrix::synthesis::thick_opt::LmConfig;
        let lm = LmConfig {
            max_iterations: 2000,
            ..LmConfig::default()
        };
        let run = |k: usize| {
            let mut r = two_env_req(k);
            let mut h = film("H", 100.0);
            h.optimize = false;
            r.design.get_mut("coat").unwrap().layers = vec![h];
            r.materials = Some(["L", "H"].map(String::from).to_vec());
            let (envs, pool, _) = build_environments(&r, &WL).unwrap();
            assert_eq!(envs.is_joint(), k > 1);
            run_environments(
                envs,
                pool,
                &WL,
                &[Radians::ZERO],
                &env_spec(k),
                one_needle_cfg(),
                // A step longer than any film: no interior site anywhere.
                NeedleCycleConfig {
                    boundary_needles: true,
                    scan_step_nm: 1000.0,
                    ..Default::default()
                },
                lm.clone(),
                |_, _| Ok(()),
            )
            .unwrap()
        };
        let (r1, s1, _) = run(1);
        let (r2, s2, _) = run(2);
        let ins = |r: &crate::smatrix::synthesis::pipeline::PipelineResult| {
            r.phases
                .iter()
                .flat_map(|p| &p.needle_results)
                .filter_map(|h| h.insertion.clone())
                .map(|i| (i.film_idx, i.material, i.boundary, i.seed_nm))
                .collect::<Vec<_>>()
        };
        assert!(
            !ins(&r1).is_empty(),
            "no needle was inserted - the test tests nothing"
        );
        assert!(ins(&r1).iter().all(|i| i.2), "{:?}", ins(&r1));
        assert_eq!(ins(&r1), ins(&r2));
        assert_eq!(s1.films().len(), s2.films().len());
        for (i, (a, b)) in s1.films().iter().zip(s2.films()).enumerate() {
            assert_eq!(a.material, b.material, "film {i}");
            assert!(
                (a.d_nm - b.d_nm).abs() <= 1e-5 * a.d_nm.max(1.0),
                "film {i}: {} vs {} nm",
                a.d_nm,
                b.d_nm
            );
        }
        assert!(
            (r2.final_mf - 2.0 * r1.final_mf).abs() <= 1e-7 * r2.final_mf.abs(),
            "{} vs 2 x {}",
            r2.final_mf,
            r1.final_mf
        );
    }

    /// `seed_edges` on the joint door seeds the segment's edge in every
    /// copy: on top forward, at the bottom of a turned copy.
    #[test]
    fn seed_edges_seed_a_segment_edge_in_every_copy() {
        let mut r = with_pool(req(
            vec![("coat", vec![film("L", 100.0), film("H", 80.0)])],
            vec![
                EnvironmentCfg {
                    name: "front".to_string(),
                    stack: vec![design_seg("coat")],
                },
                EnvironmentCfg {
                    name: "back".to_string(),
                    stack: vec![inverted_seg("coat")],
                },
            ],
        ));
        r.design.get_mut("coat").unwrap().outer = Some(vec!["H".into()]);
        let e = build_environments(&r, &WL).unwrap_err();
        assert!(e.contains("outer"), "{e}");

        r.seed_edges = Some(5.0);
        let (envs, _, _) = build_environments(&r, &WL).unwrap();
        let names = |st: &DesignStack| -> Vec<(String, f64)> {
            st.films()
                .iter()
                .map(|f| (f.material.to_string(), f.d_nm))
                .collect()
        };
        let front = names(&envs.stacks()[0]);
        assert_eq!(front[0], ("H#1".to_string(), 5.0));
        assert!(envs.stacks()[0].films()[0].optimize);
        let back = names(&envs.expand(&envs.stacks()[0]).unwrap()[1]);
        assert_eq!(back.last().unwrap(), &("H#1".to_string(), 5.0));
        assert_eq!(envs.slots().len(), 3);

        r.seed_edges = Some(0.0);
        let e = build_environments(&r, &WL).unwrap_err();
        assert!(e.contains("positive"), "{e}");
    }

    /// Inverting a fixed segment, or a segment with a graded film,
    /// refuses at the compile.
    #[test]
    fn an_inverted_reference_refuses_what_it_cannot_turn() {
        let env = |seg: EnvSegmentCfg| EnvironmentCfg {
            name: "e".to_string(),
            stack: vec![seg, design_seg("coat")],
        };
        let mut fixed_turned = fixed_seg(vec![fixed("G", 10.0, 1)]);
        fixed_turned.inverted = true;
        let r = req(
            vec![("coat", vec![film("L", 100.0)])],
            vec![env(fixed_turned)],
        );
        let e = build_environments(&r, &WL).unwrap_err();
        assert!(
            e.contains("\"inverted\" applies to a design reference"),
            "{e}"
        );

        let mut graded = film("H", 80.0);
        graded.inhomogen = true;
        let r = req(
            vec![("coat", vec![film("L", 100.0)]), ("grad", vec![graded])],
            vec![EnvironmentCfg {
                name: "e".to_string(),
                stack: vec![design_seg("coat"), inverted_seg("grad")],
            }],
        );
        let e = build_environments(&r, &WL).unwrap_err();
        assert!(
            e.contains("\"H\" is graded") && e.contains("inverted copy"),
            "{e}"
        );

        // L's slice sits on a plane behind the cover forward; turned, L is
        // the last film and the plane belongs to the substrate.
        let r = req(
            vec![("coat", faced_coat())],
            vec![
                EnvironmentCfg {
                    name: "covered".to_string(),
                    stack: vec![fixed_seg(vec![fixed("G", 50.0, 1)]), design_seg("coat")],
                },
                EnvironmentCfg {
                    name: "turned".to_string(),
                    stack: vec![inverted_seg("coat")],
                },
            ],
        );
        let e = build_environments(&r, &WL).unwrap_err();
        assert!(
            e.contains("design film \"L\"")
                && e.contains("\"covered\" puts on a plane")
                && e.contains("\"turned\" cannot"),
            "{e}"
        );
    }

    /// J1: each segment's lists are checked on its own first and last film
    /// in every environment. Here the segments run a-then-b in one
    /// environment and b-then-a behind a laminate in the other, so the
    /// assembly's outermost design film is a different segment's edge in
    /// each.
    #[test]
    fn edge_lists_belong_to_their_segment() {
        let two_segments = || {
            with_pool(req(
                vec![
                    ("a", vec![film("L", 100.0), film("H", 80.0)]),
                    ("b", vec![film("M", 60.0), film("G", 40.0)]),
                ],
                vec![
                    EnvironmentCfg {
                        name: "ab".to_string(),
                        stack: vec![design_seg("a"), design_seg("b")],
                    },
                    EnvironmentCfg {
                        name: "ba".to_string(),
                        stack: vec![
                            fixed_seg(vec![fixed("G", 300.0, 1)]),
                            design_seg("b"),
                            design_seg("a"),
                        ],
                    },
                ],
            ))
        };
        let set = |r: &mut DesignRequest, id: &str, outer: &str, inner: &str| {
            let seg = r.design.get_mut(id).unwrap();
            seg.outer = Some(vec![outer.into()]);
            seg.inner = Some(vec![inner.into()]);
        };

        let mut r = two_segments();
        set(&mut r, "a", "L", "H");
        set(&mut r, "b", "M", "G");
        let (envs, _, _) = build_environments(&r, &WL).unwrap();
        let a = envs.segment_edges("a").unwrap();
        assert_eq!(&*a.outer.as_ref().unwrap()[0].0, "L");
        assert_eq!(&*a.inner.as_ref().unwrap()[0].0, "H");
        assert!(envs.segment_edges("b").is_some());

        // b's top film is M wherever b sits, so a list without M refuses
        // and names the segment.
        set(&mut r, "b", "H", "G");
        let e = build_environments(&r, &WL).unwrap_err();
        assert!(
            e.contains("design segment \"b\"") && e.contains("outer") && e.contains("\"M\""),
            "{e}"
        );

        // One segment free, the other listed.
        let mut r = two_segments();
        r.design.get_mut("b").unwrap().inner = Some(vec!["G".into()]);
        let (envs, _, _) = build_environments(&r, &WL).unwrap();
        assert!(envs.segment_edges("a").is_none());
        assert!(envs.segment_edges("b").unwrap().outer.is_none());

        // Request-level lists do not say which of two segments they mean.
        let mut r = two_segments();
        r.outer = Some(vec!["L".into()]);
        let e = build_environments(&r, &WL).unwrap_err();
        assert!(e.contains("ambiguous"), "{e}");

        // With one segment they are its lists, but not on top of its own.
        let mut r = two_cover_req();
        r.outer = Some(vec!["L".into()]);
        r.design.get_mut("coat").unwrap().inner = Some(vec!["H".into()]);
        let e = build_environments(&r, &WL).unwrap_err();
        assert!(e.contains("one place"), "{e}");
        r.outer = None;
        r.design.get_mut("coat").unwrap().outer = Some(vec!["L".into()]);
        let (envs, _, _) = build_environments(&r, &WL).unwrap();
        assert!(envs.segment_edges("coat").is_some());
    }

    #[test]
    fn the_joint_needle_lands_in_the_design_of_every_environment() {
        use crate::smatrix::synthesis::cycle::{NeedleCycleConfig, run_needle_cycles};
        use crate::smatrix::synthesis::pipeline::SpectralInputs;

        let (envs, cmap, _) = build_environments(&two_cover_req(), &WL).unwrap();
        let before: Vec<usize> = envs.stacks().iter().map(|s| s.films().len()).collect();
        let mut stack = envs.stacks()[0].clone();
        let mut solver = ctx(env_spec(2), Some(Arc::new(envs)));
        let spectral = SpectralInputs::from_spec(&env_spec(2), &[Radians::ZERO], &WL).unwrap();
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
        // The needle's own name: fresh, of the winning pool label.
        let seed_name = &stack.films()[ins.film_idx + 1].material;
        assert!(
            seed_name.starts_with(&format!("{}#", ins.material)),
            "{seed_name}"
        );
        let mut seed_rows = Vec::new();
        for e in 0..envs.n_envs() {
            let row = envs
                .host_row(e, seed_slot)
                .expect("the seed is a parameter");
            assert_eq!(
                &envs.stacks()[e].films()[row].material,
                seed_name,
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

        // The locus: 40 nm into the first design film (L), seeded with
        // the pool's other material.
        const SLOT: usize = 0;
        const DEPTH: f64 = 40.0;
        let host_nk = &base.films()[envs.host_row(0, SLOT).unwrap()].nk;
        let partner = cmap
            .iter()
            .find(|m| &m.nk != host_nk)
            .expect("the pool offers L a needle")
            .clone();

        // Reference: per environment, analytically, through the flat doors.
        let stacks = envs.expand(&base).unwrap();
        let mut analytic = 0.0f64;
        for (e, st) in stacks.iter().enumerate() {
            let sim = solver.simulate(st).unwrap();
            let fold = build_needle_targets_env(&spec, &[Radians::ZERO], &WL, Some(&sim), e as u32)
                .unwrap();
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
        let r = with_pool(req(
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
            &[Radians::ZERO],
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
            &[Radians::ZERO],
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
            &[Radians::ZERO],
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
                &[Radians::ZERO],
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
                &[Radians::ZERO],
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
        assert!(e.contains("at least once"), "{e}");

        // A repeat of none, and a repeat on a fixed segment.
        let mut zero = design_seg("coat");
        zero.repeat = Some(0);
        let e = build_environments(&one(vec![zero]), &WL).unwrap_err();
        assert!(e.contains("segment 0") && e.contains("at least 1"), "{e}");
        let mut fixed_twice = fixed_seg(vec![]);
        fixed_twice.repeat = Some(2);
        let e = build_environments(&one(vec![fixed_twice, design_seg("coat")]), &WL).unwrap_err();
        assert!(e.contains("design reference"), "{e}");

        // Both keys / neither key on a segment.
        let both = EnvSegmentCfg {
            layers: Some(vec![]),
            design: Some("coat".to_string()),
            repeat: None,
            inverted: false,
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

    /// M1 (review PB): the per-film override cannot re-open the
    /// forced-false invariant on a surrounding.
    ///
    /// `to_row` refuses an explicit `optimize`/`needle` on a fixed row,
    /// but flag application runs global map -> row -> per-film override
    /// and the override is keyed by MATERIAL CODE — which a surrounding
    /// row carries just as a design row does. Before this refusal the
    /// override landed after the forced false and undid it, and under
    /// K > 1 the LM then moved environment 0's copy of the surrounding
    /// against environment-0-only residuals while every other
    /// environment kept the compiled thickness. Measured on the 0.7.5
    /// wheel: 500.0 -> 534.0298 nm, no refusal, and the shared film
    /// driven to its clamp floor.
    #[test]
    fn a_per_film_override_cannot_free_a_fixed_surrounding() {
        let two = || {
            req(
                vec![("coat", vec![film("L", 100.0), film("H", 80.0)])],
                vec![
                    EnvironmentCfg {
                        name: "bare".to_string(),
                        stack: vec![fixed_seg(vec![fixed("G", 500.0, 1)]), design_seg("coat")],
                    },
                    EnvironmentCfg {
                        name: "lam".to_string(),
                        stack: vec![fixed_seg(vec![fixed("G", 1000.0, 1)]), design_seg("coat")],
                    },
                ],
            )
        };

        for flag in ["optimize", "needle"] {
            let mut r = two();
            r.per_film_flags.insert(
                "G".to_string(),
                [(flag.to_string(), Value::Bool(true))]
                    .into_iter()
                    .collect(),
            );
            let e = build_environments(&r, &WL).unwrap_err();
            // The environment, the row, the flag, the code, and what to
            // do instead — the same five things `to_row`'s refusal names.
            assert!(e.contains("bare"), "{flag}: {e}");
            assert!(e.contains("bare.fixed[0][0]"), "{flag}: {e}");
            assert!(e.contains(flag), "{flag}: {e}");
            assert!(e.contains("per_film_flags"), "{flag}: {e}");
            assert!(e.contains("FIXED surrounding row"), "{flag}: {e}");
            assert!(e.contains("design segment"), "{flag}: {e}");
            assert!(e.is_ascii(), "{flag}: {e}");
        }

        // The invariant is about `true`, not about addressing the
        // material: `false` is a no-op and every other flag still lands.
        let mut r = two();
        r.per_film_flags.insert(
            "G".to_string(),
            [
                ("optimize".to_string(), Value::Bool(false)),
                ("needle".to_string(), Value::Bool(false)),
                ("roughness".to_string(), Value::from(1.5)),
            ]
            .into_iter()
            .collect(),
        );
        let (envs, _, _) = build_environments(&r, &WL).unwrap();
        let cover = &envs.stacks()[0].films()[0];
        assert_eq!(cover.rough_val, 1.5, "the benign override did not land");
        assert!(!cover.optimize && !cover.needle);

        // And a DESIGN film stays addressable — its name IS the
        // parameter's identity, so a flag on it is legitimate.
        let mut r = two();
        r.per_film_flags.insert(
            "L".to_string(),
            [("needle".to_string(), Value::Bool(false))]
                .into_iter()
                .collect(),
        );
        let (envs, _, _) = build_environments(&r, &WL).unwrap();
        let shared = &envs.stacks()[0].films()[1];
        assert_eq!(&*shared.material, "L");
        assert!(shared.optimize && !shared.needle);
    }

    /// M1's other half: the flat door is untouched by the new field.
    ///
    /// There every row is a `LayerRow` with both flags defaulting true,
    /// so there is no forced-false invariant to protect and a per-film
    /// override must keep working exactly as it did. This is the twin
    /// that would fail if `fixed_films` were ever populated by inferring
    /// from a name instead of being reported by the compile.
    #[test]
    fn the_flat_door_still_takes_a_per_film_override() {
        let mut r = flat(vec![film("L", 100.0), film("H", 80.0)]);
        r.per_film_flags.insert(
            "L".to_string(),
            [
                ("optimize".to_string(), Value::Bool(true)),
                ("needle".to_string(), Value::Bool(true)),
            ]
            .into_iter()
            .collect(),
        );
        let (envs, _, _) = build_environments(&r, &WL).unwrap();
        let f = &envs.stacks()[0].films()[0];
        assert!(f.optimize && f.needle);
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
                outer: None,
                inner: None,
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
