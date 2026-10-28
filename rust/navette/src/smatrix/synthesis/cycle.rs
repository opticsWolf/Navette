// SPDX-License-Identifier: LGPL-3.0-or-later
//! Navette -- Rust Rewrite of Numba-optimized thin-film optical solver
//!
//! synthesis::cycle — one needle pass: repeated analytic insertions.
//!
//! Port of NeedleSynthesizer.run() with compute_p_function replaced by the
//! analytic sweep (needle_pass). The Python convergence test measured the
//! MF drop from inserting a 1 nm TEST needle:
//!     improvement_py = MF_now − MF(test needle, 1 nm) ≈ −P(z*)·1 nm
//! (first-order, since P = dF/dδ). We generalize to the seed thickness:
//!     predicted_improvement = −P_best · δ_seed
//! and stop when it falls below `convergence_threshold`. Thresholds are in
//! merit-units-per-nm-of-seed — recalibrate when porting configs.

use std::collections::{BTreeMap, HashMap};
use std::sync::Arc;
use std::time::Instant;

use num_complex::Complex64;

use crate::smatrix::synthesis::context::DesignContext;
use crate::smatrix::synthesis::environments::{CompiledEnvironments, Gap, GapPlace};
use crate::smatrix::synthesis::needle_pass::{
    NeedlePassInput, ScanSite, build_boundary_sites, build_needle_targets_env, build_scan_sites,
    needle_pass_scan,
};
use crate::smatrix::synthesis::pipeline::SpectralInputs;
use crate::smatrix::synthesis::structure::{
    DesignStack, EdgeMaterial, EdgeMaterials, LayerSpec, edge_allows, is_edge_material,
};
use crate::structure::Span;

/// Knobs for the inner insertion loop (Python `NeedleConfig` subset).
#[derive(Clone, Debug)]
pub struct NeedleCycleConfig {
    /// Max insertions per pass (pipeline forwards `needles_per_cycle`).
    pub max_needles: usize,
    /// Stop when predicted improvement drops below this (merit units).
    pub convergence_threshold: f64,
    /// Seed layer thickness (nm).
    pub needle_seed_thickness_nm: f64,
    /// Scan resolution (nm).
    pub scan_step_nm: f64,
    /// Re-fold needle targets against the live sim each cycle (op-point
    /// masking for one-sided/banded/integral kinds + live PD gain shift).
    /// `false` keeps the static conservative fold (A/B benchmark switch).
    pub refold_per_cycle: bool,
    /// Record per-step wall time in each [`NeedleCycleResult::timing`]
    /// (bench switch). Off, no clock is read at all; on or off, the run
    /// itself is the same — the timers observe and never steer.
    pub record_timing: bool,
    /// Also scan every span boundary, the outer edges included (see
    /// `build_boundary_sites`): a needle there is a new film between two
    /// layers, and reaches films too thin for an interior site. Off, the
    /// scan and the run are exactly what they were. A joint run scans the
    /// shared design's boundary loci (see `MaterialPool::gap_accepts`).
    pub boundary_needles: bool,
}

impl Default for NeedleCycleConfig {
    fn default() -> Self {
        NeedleCycleConfig {
            max_needles: 10,
            convergence_threshold: 1e-4,
            needle_seed_thickness_nm: 5.0,
            scan_step_nm: 2.0,
            refold_per_cycle: true,
            record_timing: false,
            boundary_needles: false,
        }
    }
}

/// Wall-clock seconds one insertion cycle spent in each step, recorded only
/// under [`NeedleCycleConfig::record_timing`].
///
/// `scan` is the analytic sweep (site build + one scan per host
/// material over that host's sites, every environment on the joint arm); `refold` the
/// simulate + fold that precedes it; `merit` the convergence test's
/// `evaluate_merit`; `lm` the re-optimization after the insertion (0 on a
/// cycle that stopped without inserting). `initial_lm` is the pass's
/// opening optimization, carried on the pass's first record only.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CycleTiming {
    pub refold_s: f64,
    pub scan_s: f64,
    pub merit_s: f64,
    pub lm_s: f64,
    pub initial_lm_s: f64,
}

/// `Some(now)` only when timing is on: an untimed run reads no clock.
fn tick(on: bool) -> Option<Instant> {
    on.then(Instant::now)
}

/// Seconds since [`tick`], or 0 for an untimed run.
fn tock(t0: Option<Instant>) -> f64 {
    t0.map_or(0.0, |t| t.elapsed().as_secs_f64())
}

/// Record of one insertion cycle — mirrors Python `NeedleCycleResult`
/// with the scan-MF replaced by its analytic analog.
#[derive(Clone, Debug)]
pub struct NeedleCycleResult {
    pub cycle: usize,
    pub merit_before: f64,
    pub merit_after: f64,
    /// Most-negative P value at the chosen site (dF/dδ there).
    pub best_p: Option<f64>,
    /// −best_p · seed_thickness — the recalibrated convergence metric.
    pub predicted_improvement: Option<f64>,
    pub layer_count: usize,
    pub insertion: Option<Insertion>,
    /// Per-step wall time; `None` unless `record_timing` was set.
    pub timing: Option<CycleTiming>,
}

/// Where a needle landed. Inside a film, `film_idx` is the host it split.
/// On a `boundary`, it is the row the new film now occupies (the old row
/// there moved down one), and the depth is 0.
#[derive(Clone, Debug)]
pub struct Insertion {
    pub film_idx: usize,
    pub depth_into_layer_nm: f64,
    pub material: Arc<str>,
    pub boundary: bool,
    /// Thickness the needle went in at: the seed thickness, or a rough
    /// boundary's start thickness (see [`rough_boundary_start_nm`]).
    pub seed_nm: f64,
}

/// A boundary needle on a rough interface, and the thickness it starts at.
///
/// Névot-Croce depends on the media on both sides of an interface, so a
/// film slipped into a rough boundary never thins back to the rough
/// `A | B` it came from: `A | N` sharp plus `N | B` rough is a different
/// interface pair at any thickness. The merit jumps on insertion, and P,
/// which is first order in the thickness, cannot see the jump. These sites
/// are scored by evaluating the merit with the needle in place instead.
///
/// The needle starts where the roughness cap leaves the interface alone,
/// `d0 = max(seed, cap * sigma)`, so it is a layer that can carry the
/// roughness it sits on. Returns `None` for a smooth boundary. `b` is the
/// boundary above row `b` (`films.len()`: the substrate side); its
/// interface is solver slot `b + 1`, read from the capped arrays.
fn rough_boundary_start_nm(
    sa: &crate::smatrix::synthesis::structure::SolverArrays,
    b: usize,
    seed_nm: f64,
    cap: f64,
) -> Option<f64> {
    let k = b + 1;
    let sigma = sa.rough_vals[k];
    (sa.rough_types[k] != 0 && sigma > 0.0).then(|| seed_nm.max(cap * sigma))
}

/// The needle material pool: every material a needle may be made of.
///
/// A layer whose `needle` flag is set accepts a needle of any entry whose
/// index differs from its own; the scan scores every (site, material) pair
/// and inserts the single best. An entry's `material` is its label, and an
/// inserted needle is named `"<label>#<n>"`, unique in the stack. nk arrays
/// must be evaluated on the simulation grid; thickness and flags of an
/// entry are ignored (the needle is built fresh). An empty pool inserts
/// nothing.
///
/// [`Self::with_insert_rules`] narrows which entries go into which: a film
/// is matched to the entry of its own index (bit for bit, like the layer
/// merge), so needles, being pool materials, follow the rules too.
#[derive(Clone, Debug, Default)]
pub struct MaterialPool {
    materials: Vec<LayerSpec>,
    /// `allowed[h][m]`: may a film of entry `h` take a needle of entry
    /// `m`? `None`: every pair of different index.
    allowed: Option<Vec<Vec<bool>>>,
    /// Under rules, the row for a film whose index no entry carries (the
    /// `"*"` rule); `None` makes such a film an error.
    unmatched: Option<Vec<bool>>,
}

/// Insertion rules: host label -> the labels of the needles it accepts.
/// A host label with no key keeps the default (every entry of another
/// index); an empty list accepts nothing. The key `"*"` is the rule for a
/// needle-flagged film whose index matches no pool entry.
pub type InsertRules = BTreeMap<String, Vec<String>>;

impl std::ops::Deref for MaterialPool {
    type Target = [LayerSpec];
    fn deref(&self) -> &[LayerSpec] {
        &self.materials
    }
}

impl From<Vec<LayerSpec>> for MaterialPool {
    fn from(materials: Vec<LayerSpec>) -> Self {
        MaterialPool {
            materials,
            ..Default::default()
        }
    }
}

impl MaterialPool {
    /// The empty pool: the needle pass inserts nothing.
    pub fn new() -> Self {
        Self::default()
    }

    /// The entry whose index is `nk`, bit for bit.
    pub fn entry_of(&self, nk: &[Complex64]) -> Option<usize> {
        self.materials.iter().position(|m| *m.nk == *nk)
    }

    fn label_index(&self, label: &str, what: &str) -> Result<usize, String> {
        let mut hits = self
            .materials
            .iter()
            .enumerate()
            .filter(|(_, m)| &*m.material == label);
        match (hits.next(), hits.next()) {
            (Some((i, _)), None) => Ok(i),
            (None, _) => Err(format!(
                "{what}: {label:?} is not a pool material; the pool holds [{}]",
                self.labels().join(", ")
            )),
            (Some(_), Some(_)) => Err(format!(
                "{what}: two pool materials are labelled {label:?}; give them distinct labels"
            )),
        }
    }

    fn labels(&self) -> Vec<String> {
        self.materials
            .iter()
            .map(|m| format!("{:?}", m.material))
            .collect()
    }

    /// This pool with `rules` applied (see [`InsertRules`]). Empty rules
    /// leave it as it is. A label the pool does not hold is an error, and
    /// so is a host listed as its own needle: that is no needle at all.
    pub fn with_insert_rules(mut self, rules: &InsertRules) -> Result<Self, String> {
        if rules.is_empty() {
            return Ok(self);
        }
        let n = self.materials.len();
        let mut allowed: Vec<Vec<bool>> =
            (0..n).map(|h| (0..n).map(|m| m != h).collect()).collect();
        let mut unmatched = None;
        for (host, needles) in rules {
            let h = if host == "*" {
                None
            } else {
                Some(self.label_index(host, "insert rules (host)")?)
            };
            let mut row = vec![false; n];
            for needle in needles {
                let m = self.label_index(needle, "insert rules (needle)")?;
                if Some(m) == h {
                    return Err(format!(
                        "insert rules: {host:?} lists itself; a needle of a film's own \
                         index is no needle"
                    ));
                }
                row[m] = true;
            }
            match h {
                Some(h) => allowed[h] = row,
                None => unmatched = Some(row),
            }
        }
        self.allowed = Some(allowed);
        self.unmatched = unmatched;
        Ok(self)
    }

    /// The edge materials named by `outer` and `inner` (pool labels).
    pub fn edge_materials(
        &self,
        outer: Option<&[String]>,
        inner: Option<&[String]>,
    ) -> Result<EdgeMaterials, String> {
        let resolve =
            |l: Option<&[String]>, what: &str| -> Result<Option<Vec<EdgeMaterial>>, String> {
                l.map(|l| {
                    l.iter()
                        .map(|label| {
                            let m = &self.materials[self.label_index(label, what)?];
                            Ok((m.material.clone(), m.nk.clone()))
                        })
                        .collect()
                })
                .transpose()
            };
        Ok(EdgeMaterials {
            outer: resolve(outer, "outer edge")?,
            inner: resolve(inner, "inner edge")?,
        })
    }

    /// Does `host` accept a needle of entry `m`? Without rules, exactly when
    /// the indices differ.
    fn accepts(&self, host: &LayerSpec, m: usize) -> Result<bool, String> {
        let Some(allowed) = &self.allowed else {
            return Ok(host.nk != self.materials[m].nk);
        };
        match (self.entry_of(&host.nk), &self.unmatched) {
            (Some(h), _) => Ok(allowed[h][m]),
            (None, Some(row)) => Ok(row[m]),
            (None, None) => Err(format!(
                "insert rules: film {:?} accepts needles but its index is no pool \
                 material's, so no rule names it. Add its material to the pool, add a \
                 \"*\" rule, or clear its needle flag.",
                host.material
            )),
        }
    }

    /// The `sites` of `stack` that accept a needle of entry `m`. Filtering
    /// whole rows keeps each surviving row's steps intact.
    fn sites_for(
        &self,
        m: usize,
        stack: &DesignStack,
        sites: &[ScanSite],
    ) -> Result<Vec<ScanSite>, String> {
        let mut out = Vec::new();
        for s in sites {
            let ok = if s.boundary {
                self.boundary_accepts(stack, s.film_idx, m)?
            } else {
                self.accepts(&stack.films()[s.film_idx], m)?
            };
            if ok {
                out.push(*s);
            }
        }
        Ok(out)
    }

    /// May a needle of entry `m` go on the boundary above row `b` (`b ==
    /// films.len()`: on the substrate)?
    ///
    /// The needle becomes a film between the two neighbours, so it is
    /// checked against both:
    ///  * its index differs from either neighbour's, the adjacent row's
    ///    and a layer's own bulk row's (the ambient and the substrate
    ///    count): a needle of a neighbour's index is that neighbour grown
    ///    thicker, which the optimizer already covers;
    ///  * every neighbour that is a needle host accepts it under the
    ///    insertion rules (a closed host closes its boundaries too);
    ///  * on an outer edge it becomes the edge layer, so it must be on
    ///    that edge's list when one is set.
    fn boundary_accepts(&self, stack: &DesignStack, b: usize, m: usize) -> Result<bool, String> {
        let films = stack.films();
        let n = films.len();
        let nk = &self.materials[m].nk;
        let adjacent = [
            b.checked_sub(1).map_or(stack.ambient(), |r| &films[r]),
            films.get(b).unwrap_or(stack.substrate()),
        ];
        if adjacent.iter().any(|f| f.nk == *nk) {
            return Ok(false);
        }
        let neighbours = [b.checked_sub(1), (b < n).then_some(b)];
        for row in neighbours.into_iter().flatten() {
            let sp = stack.spans()[stack.span_index_of_row(row)];
            if !sp.is_singleton_bulk() {
                continue; // a graded profile: excluded by its adjacent row only
            }
            let layer = &films[sp.bulk_start];
            if layer.nk == *nk || (layer.needle && !self.accepts(layer, m)?) {
                return Ok(false);
            }
        }
        let edges = stack.edge_materials();
        let needle = &self.materials[m];
        let listed = |l: &Option<Vec<EdgeMaterial>>| {
            l.as_ref()
                .is_none_or(|l| l.iter().any(|e| is_edge_material(e, needle)))
        };
        Ok((b != 0 || listed(&edges.outer)) && (b != n || listed(&edges.inner)))
    }

    /// May a needle of entry `m` go on boundary locus `gap` of the shared
    /// design? `stacks` are the K expanded assemblies and `places[e]` the
    /// gap's placements in `stacks[e]`, copy by copy.
    ///
    /// The flat rules, read per placement:
    ///  * some design film beside the gap is a needle host, and no span
    ///    beside it, inside or outside, is empty in any placement (that
    ///    plane is also the next one);
    ///  * the needle differs in index from the design films beside it, in
    ///    every placement (adjacent row and bulk row), and every one of
    ///    them that is a host accepts it under the insertion rules;
    ///  * past the segment's edge the neighbour is no design parameter: it
    ///    differs per environment and brings no rules. It excludes the
    ///    needle only when it has the needle's index in every placement,
    ///    where the needle would do nothing anywhere;
    ///  * on the segment's top edge the needle is on its outer list, on its
    ///    bottom edge on its inner list;
    ///  * no film's own interface slice moves (see
    ///    [`CompiledEnvironments::gap_keeps_slices`]).
    pub(crate) fn gap_accepts(
        &self,
        envs: &CompiledEnvironments,
        stacks: &[DesignStack],
        gap: Gap,
        places: &[Vec<GapPlace>],
        m: usize,
    ) -> Result<bool, String> {
        let nk = &self.materials[m].nk;
        let (above, below) = envs.gap_neighbours(gap);
        let design: Vec<usize> = [above, below].into_iter().flatten().collect();
        let primary = &stacks[0];
        let thick =
            |st: &DesignStack, sp: &Span| st.films()[sp.start..sp.end].iter().any(|f| f.d_nm > 0.0);
        let mut host = false;
        for &slot in &design {
            let sp = primary.spans()[envs.span_of(slot, 0).expect("a design slot has a primary")];
            if !sp.is_singleton_bulk() {
                continue; // a graded profile: excluded by its adjacent row only
            }
            let layer = &primary.films()[sp.bulk_start];
            if layer.nk == *nk || (layer.needle && !self.accepts(layer, m)?) {
                return Ok(false);
            }
            host |= layer.needle && thick(primary, &sp);
        }
        if !host {
            return Ok(false);
        }
        // Per side, whether every outside neighbour so far has the needle's
        // index; a side with no outside neighbour never excludes.
        let mut outside_same = [None::<bool>; 2];
        for (e, st) in stacks.iter().enumerate() {
            let spans = st.spans();
            let films = st.films();
            for (c, pl) in places[e].iter().enumerate() {
                if !envs.gap_keeps_slices(pl) {
                    return Ok(false);
                }
                let mine: Vec<usize> = design.iter().map(|&s| envs.spans_of(s, e)[c]).collect();
                // The spans above and below the new film, and the rows that
                // will touch it.
                let sides = [
                    (pl.at.checked_sub(1), pl.row.checked_sub(1)),
                    (
                        (pl.at < spans.len()).then_some(pl.at),
                        (pl.row < films.len()).then_some(pl.row),
                    ),
                ];
                for (side, (span, row)) in sides.into_iter().enumerate() {
                    let same = match (span, row) {
                        (Some(si), Some(r)) => {
                            let sp = &spans[si];
                            if !thick(st, sp) {
                                return Ok(false);
                            }
                            films[r].nk == *nk
                                || (sp.is_singleton_bulk() && films[sp.bulk_start].nk == *nk)
                        }
                        _ => {
                            let half = if side == 0 {
                                st.ambient()
                            } else {
                                st.substrate()
                            };
                            half.nk == *nk
                        }
                    };
                    if span.is_some_and(|si| mine.contains(&si)) {
                        if same {
                            return Ok(false);
                        }
                    } else {
                        let o = outside_same[side].get_or_insert(true);
                        *o &= same;
                    }
                }
            }
        }
        if outside_same.contains(&Some(true)) {
            return Ok(false);
        }
        let (edges, top, bottom) = envs.gap_edges(gap);
        let needle = &self.materials[m];
        let listed = |l: Option<&Vec<EdgeMaterial>>| {
            l.is_none_or(|l| l.iter().any(|e| is_edge_material(e, needle)))
        };
        Ok((!top || listed(edges.and_then(|e| e.outer.as_ref())))
            && (!bottom || listed(edges.and_then(|e| e.inner.as_ref()))))
    }
}

/// Rules for one needle run beyond the pool itself: which needles go into
/// which materials, and which materials may form the edge layers.
#[derive(Clone, Debug, Default)]
pub struct NeedleRules {
    pub insert: Option<InsertRules>,
    /// Allowed materials (pool labels) of the layer facing the ambient.
    pub outer: Option<Vec<String>>,
    /// Allowed materials (pool labels) of the layer on the substrate.
    pub inner: Option<Vec<String>>,
    /// A starting design off an edge list gets a film of the list's first
    /// material, this thick (nm), on that edge instead of an error.
    /// `None`: the design must comply as given.
    pub seed_edges: Option<f64>,
}

impl NeedleRules {
    /// Apply to a run: the pool takes the insertion rules, the stack the
    /// edge materials. Its current edges must comply, or be seeded
    /// (`seed_edges`).
    pub fn apply(
        &self,
        stack: &mut DesignStack,
        pool: MaterialPool,
    ) -> Result<MaterialPool, String> {
        let pool = match &self.insert {
            Some(r) => pool.with_insert_rules(r)?,
            None => pool,
        };
        if self.outer.is_some() || self.inner.is_some() {
            let edges = pool.edge_materials(self.outer.as_deref(), self.inner.as_deref())?;
            if let Some(d) = self.seed_edges {
                seed_edges(stack, &edges, d)?;
            }
            stack.set_edge_materials(edges)?;
        }
        Ok(pool)
    }
}

/// Put a `d_nm` film of the list's first material on each edge whose
/// current layer is off its list. The film is a needle like any other
/// (named `"<label>#<n>"`, optimized, a host), and the old edge layer is
/// kept behind it.
fn seed_edges(stack: &mut DesignStack, edges: &EdgeMaterials, d_nm: f64) -> Result<(), String> {
    if !(d_nm.is_finite() && d_nm > 0.0) {
        return Err(format!(
            "seed_edges: the seed film's thickness must be positive and finite (got {d_nm} nm)"
        ));
    }
    for (list, outer) in [(&edges.outer, true), (&edges.inner, false)] {
        let Some(first) = list.as_ref().and_then(|l| l.first()) else {
            continue;
        };
        let spans = stack.spans();
        let edge = if outer { spans.first() } else { spans.last() };
        if edge.is_some_and(|sp| edge_allows(list, stack.films(), sp)) {
            continue;
        }
        let film = new_needle(stack.films(), &pool_entry(&first.0, first.1.clone()), d_nm);
        let b = if outer { 0 } else { stack.films().len() };
        stack.insert_film_at_boundary(b, film)?;
    }
    Ok(())
}

/// A film name for a new needle of `label`, used by no film in `films`.
fn needle_name(films: &[LayerSpec], label: &str) -> Arc<str> {
    (1usize..)
        .map(|n| format!("{label}#{n}"))
        .find(|c| films.iter().all(|f| &*f.material != c.as_str()))
        .map(Arc::from)
        .expect("an unused needle name exists")
}

/// A pool entry: `nk` under `label`, the needle's own fresh flags.
pub fn pool_entry(label: &str, nk: Arc<[Complex64]>) -> LayerSpec {
    LayerSpec {
        material: Arc::from(label),
        nk,
        d_nm: 0.0,
        coherent: true,
        rough_type: 0,
        rough_val: 0.0,
        optimize: true,
        needle: true,
    }
}

/// `entries` with every repeated index dropped (the first label wins):
/// two entries of one index would sweep the same needle twice.
pub fn dedup_pool(entries: impl IntoIterator<Item = LayerSpec>) -> MaterialPool {
    let mut out: Vec<LayerSpec> = Vec::new();
    for e in entries {
        if !out.iter().any(|m| m.nk == e.nk) {
            out.push(e);
        }
    }
    out.into()
}

/// A fresh needle of pool entry `m`, `d_nm` thick, named apart from every
/// film in `films`.
pub(crate) fn new_needle(films: &[LayerSpec], m: &LayerSpec, d_nm: f64) -> LayerSpec {
    LayerSpec {
        material: needle_name(films, &m.material),
        nk: m.nk.clone(),
        d_nm,
        coherent: true,
        rough_type: 0,
        rough_val: 0.0,
        optimize: true,
        needle: true,
    }
}

/// The default pool: one entry per distinct index among the stack's
/// LAYERS, labelled by the first layer that carries it.
///
/// A layer is the bulk row of a singleton-bulk span, so graded sublayers
/// and interface slices never contribute an index of their own. `keep`
/// narrows the rows considered: the joint door keeps the shared design's
/// rows, so a surrounding's material is never a needle by default.
pub fn pool_from_layers(stack: &DesignStack, keep: impl Fn(usize) -> bool) -> MaterialPool {
    dedup_pool(
        stack
            .spans()
            .iter()
            .filter(|sp| sp.is_singleton_bulk() && keep(sp.bulk_start))
            .map(|sp| {
                let l = &stack.films()[sp.bulk_start];
                pool_entry(&l.material, l.nk.clone())
            }),
    )
}

/// F2.3: one candidate locus, agreed on by every environment.
///
/// A locus is a *shared design parameter* plus a step into it, never a row
/// — row 7 of one assembly and row 3 of another are the same physical
/// place when both answer to the same slot.
#[derive(Clone, Debug)]
struct JointCandidate {
    site: JointSite,
    /// The needle material: the index into the pool whose sweep won.
    material: usize,
    /// P summed over environments.
    p: f64,
}

/// Where a joint needle goes: inside a design film, or on a boundary
/// locus as a new design film.
#[derive(Clone, Copy, Debug)]
enum JointSite {
    Interior {
        /// Index into `CompiledEnvironments::slots`.
        slot: usize,
        depth_into_layer_nm: f64,
    },
    Boundary(Gap),
}

/// What the joint sweep found: the best P-scored candidate, the gain
/// shift summed over the K folds, the number of admissible design sites
/// (zero: "stack too thin"), and the rough boundary loci left for a trial
/// as (locus, pool index, start thickness).
struct JointSweep {
    best: Option<JointCandidate>,
    gtot: f64,
    n_sites: usize,
    rough: Vec<(Gap, usize, f64)>,
}

/// F2.3: the joint analytic sweep — K scans, one candidate.
///
/// Each environment scans ITS OWN assembly with ITS OWN fold (an operating
/// point masks one-sided and banded kinds, and environment 1's operating
/// point is not environment 0's), then every site that lands in the shared
/// design is routed **by design parameter** into a shared bucket and summed.
/// The sum is the right combination because the joint merit is a sum of the
/// per-environment merits, so dF/dδ at a shared locus is the sum of the
/// per-environment slopes — the same identity F2.2's residual concatenation
/// rests on, one derivative up.
///
/// Surroundings never deposit: their rows answer `None` to `slot_of_row`,
/// which is the second lock on the door the `needle` flag already holds.
///
/// Under `boundary_needles` the boundary loci of the design are scored the
/// same way: one site per placement, in every copy in every environment,
/// summed per (material, locus). A locus whose plane is rough in some
/// placement is set aside for a trial, as on the flat arm, at the largest
/// start thickness over its placements.
///
/// Returns the winner (most negative ΣP, `None` when nothing improves), the
/// gain shift summed over the K folds for the convergence test, the number
/// of admissible design sites (zero means "stack too thin") and the rough
/// loci. Under `record_timing` the K re-folds' wall time is added to
/// `refold_s`, so the caller can split it out of the sweep's total.
fn joint_needle_sweep<C: DesignContext + ?Sized>(
    ctx: &C,
    stack: &DesignStack,
    spectral: &SpectralInputs,
    pool: &MaterialPool,
    cfg: &NeedleCycleConfig,
    refold_s: &mut f64,
) -> Result<JointSweep, String> {
    let envs = ctx
        .environments()
        .ok_or_else(|| "joint_needle_sweep: the context has no environments".to_string())?;
    let stacks = envs.expand(stack)?;

    let mut gtot = 0.0f64;
    // (pool index, slot, step into the host) -> (depth, ΣP, environments seen)
    let mut buckets: HashMap<(usize, usize, usize), (f64, f64, usize)> = HashMap::new();
    let mut n_sites = 0usize;

    // The boundary loci each material may take, P-scored or rough.
    let mut gap_sets: Vec<Vec<(Gap, Vec<Vec<GapPlace>>)>> = vec![Vec::new(); pool.len()];
    let mut rough = Vec::new();
    if cfg.boundary_needles {
        let solver: Vec<_> = stacks.iter().map(DesignStack::solver_arrays).collect();
        for gap in envs.gaps() {
            let places: Vec<Vec<GapPlace>> = stacks
                .iter()
                .enumerate()
                .map(|(e, st)| envs.gap_places(e, st.spans(), st.films().len(), gap))
                .collect();
            let d0 = places
                .iter()
                .flatten()
                .filter_map(|pl| {
                    rough_boundary_start_nm(
                        &solver[pl.env],
                        pl.row,
                        cfg.needle_seed_thickness_nm,
                        stack.roughness_cap(),
                    )
                })
                .reduce(f64::max);
            for (m, set) in gap_sets.iter_mut().enumerate() {
                if !pool.gap_accepts(envs, &stacks, gap, &places, m)? {
                    continue;
                }
                n_sites += 1;
                match d0 {
                    Some(d0) => rough.push((gap, m, d0)),
                    None => set.push((gap, places.clone())),
                }
            }
        }
    }
    // (pool index, locus) -> (ΣP, placements seen)
    let mut gap_buckets: HashMap<(usize, Gap), (f64, usize)> = HashMap::new();

    for (e, st) in stacks.iter().enumerate() {
        let mut fold = spectral.folds[e].clone();
        let t0 = tick(cfg.record_timing);
        if cfg.refold_per_cycle
            && let Ok(sim) = ctx.simulate(st)
            && let Ok(f) = build_needle_targets_env(
                &spectral.spec,
                &spectral.angles,
                &spectral.wavls,
                Some(&sim),
                e as u32,
            )
        {
            fold = f;
        }
        *refold_s += tock(t0);
        gtot += fold.phi_gain_shift.iter().sum::<f64>();

        let sa = st.solver_arrays();
        for (mi, mat) in pool.iter().enumerate() {
            let input = NeedlePassInput {
                n_stack_cache: &sa.n_stack_cache,
                thicknesses: &sa.thicknesses,
                rough_types: &sa.rough_types,
                rough_vals: &sa.rough_vals,
                n_layers: sa.n_layers as usize,
                wavls: spectral.wavls.as_slice(),
                sin_theta: spectral.sin_theta.as_slice(),
                fold: &fold,
                needle_n_per_wav: &mat.nk,
                start_idx: 0,
                end_idx: (sa.n_layers - 1) as usize,
                calc_s: true,
                calc_p: true,
            };
            // Only the rows that accept this needle. A design row has the
            // same index in every environment, so every environment keeps
            // the same rows.
            let all = envs.scan_sites(e, st, stack, cfg.scan_step_nm);
            let mut sites = pool.sites_for(mi, st, &all)?;
            let n_interior = sites.len();
            // The boundary sites follow, one per placement of each locus
            // this material may take, at the depth the locator sums to.
            let mut gap_of_site = Vec::new();
            let films = st.films();
            for (gap, places) in &gap_sets[mi] {
                for pl in &places[e] {
                    let z_nm = films[..pl.row].iter().map(|f| f.d_nm).sum();
                    sites.push(ScanSite {
                        film_idx: pl.row,
                        depth_into_layer_nm: 0.0,
                        z_nm,
                        boundary: true,
                    });
                    gap_of_site.push(*gap);
                }
            }
            let res = needle_pass_scan(&input, &sites)?;
            for (gap, &p) in gap_of_site.iter().zip(&res.p_profile[n_interior..]) {
                let entry = gap_buckets.entry((mi, *gap)).or_insert((0.0, 0));
                entry.0 += p;
                entry.1 += 1;
            }
            // The scan walks rows in order and emits steps 1..n_steps
            // inside each admissible forward row, so the running count per
            // row IS the step index — and it is the same index in every
            // environment, because a design row has the same thickness
            // everywhere by construction (`expand` copies it). An inverted
            // copy's grid runs the other way; its step index is read back
            // from the depth into the forward film, which is also the
            // depth the bucket keeps.
            let mut step: HashMap<usize, usize> = HashMap::new();
            for (site, &p) in res.sites.iter().zip(res.p_profile.iter()).take(n_interior) {
                let n = step.entry(site.film_idx).or_insert(0);
                let mut ordinal = *n;
                *n += 1;
                let Some((slot, copy)) = envs.copy_of_row(e, site.film_idx) else {
                    continue; // surroundings: inadmissible, and never shared
                };
                let mut depth = site.depth_into_layer_nm;
                if let Some(u) = envs.forward_depth(e, site.film_idx, depth, stack) {
                    depth = u;
                    ordinal = ((u / cfg.scan_step_nm).round() as usize).saturating_sub(1);
                }
                if e == 0 && copy == 0 {
                    n_sites += 1; // the design's own sites, each once
                }
                let entry = buckets
                    .entry((mi, slot, ordinal))
                    .or_insert((depth, 0.0, 0));
                entry.1 += p;
                entry.2 += 1;
            }
        }
    }

    // §4.4's alignment assert, in its span-aware form: the same design
    // object must produce the same scan grid in every place it is used, so
    // every bucket is filled once per copy in every environment. A short
    // bucket means the assemblies drifted apart, and summing it would sum
    // over different loci.
    if let Some(((m, slot, ord), (_, _, n))) = buckets
        .iter()
        .find(|((_, slot, _), v)| v.2 != envs.n_placements(*slot))
    {
        let m = &pool[*m].material;
        let k = envs.n_placements(*slot);
        return Err(format!(
            "joint needle sweep: the locus (slot {slot}, step {ord}, needle {m:?}) \
             appears in {n} of the {k} places the design film is used - the \
             assemblies no longer share the design's scan grid"
        ));
    }

    // A boundary locus is placed once per copy of its segment.
    if let Some(((m, gap), (_, n))) = gap_buckets.iter().find(|((_, gap), v)| {
        let first = envs.segments()[gap.seg].start;
        v.1 != envs.n_placements(first)
    }) {
        let m = &pool[*m].material;
        return Err(format!(
            "joint needle sweep: the boundary locus (segment {}, before film {}, \
             needle {m:?}) appears in {n} of its places - the assemblies no longer \
             share the design's boundaries",
            gap.seg, gap.k
        ));
    }

    // Deterministic order before the min: a hash map's iteration order is
    // not stable, and two loci can tie on P exactly (identical
    // environments, mirrored stacks). Per material, interior sites come
    // before boundary loci, so an exact tie goes to the interior site as on
    // the flat arm.
    let mut items: Vec<((usize, usize, usize, usize), JointSite, f64)> = buckets
        .into_iter()
        .map(|((m, slot, ord), (depth, p, _))| {
            let site = JointSite::Interior {
                slot,
                depth_into_layer_nm: depth,
            };
            ((m, 0, slot, ord), site, p)
        })
        .chain(
            gap_buckets
                .into_iter()
                .map(|((m, gap), (p, _))| ((m, 1, gap.seg, gap.k), JointSite::Boundary(gap), p)),
        )
        .collect();
    items.sort_by_key(|a| a.0);

    let mut best: Option<JointCandidate> = None;
    for ((material, ..), site, p) in items {
        if p >= 0.0 {
            continue; // not an improvement anywhere
        }
        let better = match &best {
            None => true,
            Some(b) => p < b.p,
        };
        if better {
            best = Some(JointCandidate { site, material, p });
        }
    }
    Ok(JointSweep {
        best,
        gtot,
        n_sites,
        rough,
    })
}

/// Run one full needle pass on `stack`.
///
/// Mirrors `synth.run(max_needles = cfg.max_needles)`:
/// initial optimization → [re-fold → scan → select → insert → optimize]×N.
/// Each insertion is the single best (site, material) pair over every
/// needle-flagged layer and every `pool` entry that layer accepts.
///
/// The fold refreshes against the live sim each cycle (when
/// `cfg.refold_per_cycle` and the context simulates): one-sided/banded
/// demands mask at the operating point and PD gain shifts go live, so
/// insertion tracks the merit the optimizer actually sees. Any simulate
/// or fold failure keeps the previous fold (mock contexts without a
/// solver run the static conservative fold throughout).
pub fn run_needle_cycles<C: DesignContext + ?Sized>(
    ctx: &mut C,
    stack: &mut DesignStack,
    spectral: &SpectralInputs,
    pool: &MaterialPool,
    cfg: &NeedleCycleConfig,
) -> Result<Vec<NeedleCycleResult>, String> {
    use crate::smatrix::synthesis::needle_pass::build_needle_targets;

    let wavls = &spectral.wavls;
    let sin_theta = &spectral.sin_theta;
    let mut fold = spectral.folds[0].clone();
    let mut history = Vec::new();

    // F2.3: §4.6's one branch, read once, outside every loop. `stack` is
    // the shared design object and IS environment 0's assembly; the other
    // K-1 are re-expressed from it inside the sweep.
    let multi = ctx.environments().is_some();

    // Initial optimization.
    let t0 = tick(cfg.record_timing);
    ctx.optimize_thicknesses(stack)?;
    let initial_lm_s = tock(t0);

    for cycle in 0..cfg.max_needles {
        if multi {
            if joint_cycle(ctx, stack, spectral, pool, cfg, cycle, &mut history)? {
                break;
            }
            continue;
        }
        let mut t = CycleTiming::default();
        // Refresh the fold against the live operating point.
        let t0 = tick(cfg.record_timing);
        if cfg.refold_per_cycle
            && let Ok(sim) = ctx.simulate(stack)
            && let Ok(f) = build_needle_targets(
                &spectral.spec,
                &spectral.angles,
                &spectral.wavls,
                Some(&sim),
            )
        {
            fold = f;
        }
        t.refold_s = tock(t0);
        let t0 = tick(cfg.record_timing);
        // 1. Candidate sites: every step inside a needle-flagged layer.
        //    `build_scan_sites` reads the flag and skips multi-row spans, so
        //    interface slices and graded sublayers are never sites.
        //    Boundary sites follow the interior ones, so an exact tie goes
        //    to the interior site.
        //    Boundary sites on a rough interface are set aside: P cannot
        //    score them (see `rough_boundary_start_nm`), a trial does.
        let sa = stack.solver_arrays();
        let mut sites = build_scan_sites(stack.films(), stack.spans(), cfg.scan_step_nm);
        let mut rough_sites = Vec::new();
        if cfg.boundary_needles {
            for s in build_boundary_sites(stack.films(), stack.spans()) {
                if rough_boundary_start_nm(&sa, s.film_idx, 0.0, 0.0).is_some() {
                    rough_sites.push(s);
                } else {
                    sites.push(s);
                }
            }
        }

        // 2. Analytic sweep (both polarizations summed) once per pool
        //    material, over the sites whose layer accepts it (its index
        //    differs); the global best (site, material) pair is kept.
        //    Insertion rules narrow the pairs further.
        let per_material: Vec<Vec<_>> = (0..pool.len())
            .map(|m| pool.sites_for(m, stack, &sites))
            .collect::<Result<_, _>>()?;
        let rough_per_material: Vec<Vec<_>> = (0..pool.len())
            .map(|m| pool.sites_for(m, stack, &rough_sites))
            .collect::<Result<_, _>>()?;
        if per_material
            .iter()
            .chain(&rough_per_material)
            .all(Vec::is_empty)
        {
            break; // "Stack too thin for needle insertion" (or no pool).
        }

        let mut best: Option<(Insertion, &LayerSpec, f64)> = None;
        for (mat, mat_sites) in pool.iter().zip(&per_material) {
            if mat_sites.is_empty() {
                continue;
            }
            let input = NeedlePassInput {
                n_stack_cache: &sa.n_stack_cache,
                thicknesses: &sa.thicknesses,
                rough_types: &sa.rough_types,
                rough_vals: &sa.rough_vals,
                n_layers: sa.n_layers as usize,
                wavls: wavls.as_slice(),
                sin_theta: sin_theta.as_slice(),
                fold: &fold,
                needle_n_per_wav: &mat.nk,
                start_idx: 0,
                end_idx: (sa.n_layers - 1) as usize,
                calc_s: true,
                calc_p: true,
            };
            let res = needle_pass_scan(&input, mat_sites)?;
            if let Some((site, p)) = res.best() {
                let better = match &best {
                    None => true,
                    Some((_, _, bp)) => p < *bp,
                };
                if better {
                    best = Some((
                        Insertion {
                            film_idx: site.film_idx,
                            depth_into_layer_nm: site.depth_into_layer_nm,
                            material: mat.material.clone(),
                            boundary: site.boundary,
                            seed_nm: cfg.needle_seed_thickness_nm,
                        },
                        mat,
                        p,
                    ));
                }
            }
        }

        t.scan_s = tock(t0);

        // 3. Convergence check on the predicted improvement.
        let t0 = tick(cfg.record_timing);
        let current_mf = ctx.evaluate_merit(stack)?;
        t.merit_s = tock(t0);

        // Inserting a seed of thickness δ also grows the equivalent-medium
        // reference by δ, so the true slope is P(z*) + Σ gain_shift (0 for
        // absolute-phase demand sets; the site itself is unaffected, only
        // this bookkeeping).
        let gtot: f64 = fold.phi_gain_shift.iter().sum();
        // (insertion, material, P if scored by P, improvement it promises)
        let mut best: Option<(Insertion, &LayerSpec, Option<f64>, f64)> =
            best.map(|(ins, mat, p)| {
                (
                    ins,
                    mat,
                    Some(p),
                    -(p + gtot) * cfg.needle_seed_thickness_nm,
                )
            });

        // Rough boundary sites: the merit with the needle in place, at its
        // start thickness. The trials count to the scan's time.
        let t0 = tick(cfg.record_timing);
        for (mat, mat_sites) in pool.iter().zip(&rough_per_material) {
            for site in mat_sites {
                let d0 = rough_boundary_start_nm(
                    &sa,
                    site.film_idx,
                    cfg.needle_seed_thickness_nm,
                    stack.roughness_cap(),
                )
                .expect("rough sites were picked by the same test");
                let mut trial = stack.clone();
                trial.insert_film_at_boundary(site.film_idx, new_needle(stack.films(), mat, d0))?;
                let gain = current_mf - ctx.evaluate_merit(&trial)?;
                if best.as_ref().is_none_or(|b| gain > b.3) {
                    best = Some((
                        Insertion {
                            film_idx: site.film_idx,
                            depth_into_layer_nm: 0.0,
                            material: mat.material.clone(),
                            boundary: true,
                            seed_nm: d0,
                        },
                        mat,
                        None,
                        gain,
                    ));
                }
            }
        }
        t.scan_s += tock(t0);

        let timing = cfg.record_timing.then_some(t);
        let result = match best {
            None => {
                // No improving site anywhere.
                history.push(NeedleCycleResult {
                    cycle: cycle + 1,
                    merit_before: current_mf,
                    merit_after: current_mf,
                    best_p: None,
                    predicted_improvement: None,
                    layer_count: stack.films().len(),
                    insertion: None,
                    timing,
                });
                break;
            }
            Some((ins, mat, p, predicted)) => {
                if predicted < cfg.convergence_threshold {
                    history.push(NeedleCycleResult {
                        cycle: cycle + 1,
                        merit_before: current_mf,
                        merit_after: current_mf,
                        best_p: p,
                        predicted_improvement: Some(predicted),
                        layer_count: stack.films().len(),
                        insertion: None,
                        timing,
                    });
                    break; // "Convergence reached — stopping."
                }

                // 4. Insert the needle (split the host, Python
                //    `_insert_needle`) under a fresh name: a needle is a
                //    layer of its own and accepts needles in later cycles.
                //    A boundary needle splits nothing: it is a new film. A
                //    trial-scored one goes in at its start thickness.
                let seed = new_needle(stack.films(), mat, ins.seed_nm);
                if ins.boundary {
                    stack.insert_film_at_boundary(ins.film_idx, seed)?;
                } else {
                    stack.insert_needle_seed(ins.film_idx, ins.depth_into_layer_nm, seed)?;
                }

                // 5. Re-optimize.
                let t0 = tick(cfg.record_timing);
                let new_mf = ctx.optimize_thicknesses(stack)?;
                t.lm_s = tock(t0);

                NeedleCycleResult {
                    cycle: cycle + 1,
                    merit_before: current_mf,
                    merit_after: new_mf,
                    best_p: p,
                    predicted_improvement: Some(predicted),
                    layer_count: stack.films().len(),
                    insertion: Some(ins),
                    timing: cfg.record_timing.then_some(t),
                }
            }
        };
        history.push(result);
    }

    // The opening optimization belongs to no insertion; it rides on the
    // pass's first record (lost only when the stack was too thin to scan).
    if let Some(t) = history.first_mut().and_then(|h| h.timing.as_mut()) {
        t.initial_lm_s = initial_lm_s;
    }

    Ok(history)
}

/// F2.3: one cycle of the joint arm. Returns `true` when the pass is over.
///
/// Same five steps as the K == 1 body — sweep, convergence test, insert,
/// re-optimize, record — with two differences and no others. The sweep is
/// [`joint_needle_sweep`]'s (K scans summed by design parameter), and the
/// insertion happens TWICE: once on the shared design object and once on
/// the compile, so every environment's template changes at its own row for
/// the same parameter. Nothing else about the insertion changes: the same
/// seed, the same depth, the same `insert_needle_seed` (or, on a boundary,
/// a new film of its own), the same re-optimization through the context.
/// Rough boundary loci are scored by a trial, as on the flat arm.
fn joint_cycle<C: DesignContext + ?Sized>(
    ctx: &mut C,
    stack: &mut DesignStack,
    spectral: &SpectralInputs,
    pool: &MaterialPool,
    cfg: &NeedleCycleConfig,
    cycle: usize,
    history: &mut Vec<NeedleCycleResult>,
) -> Result<bool, String> {
    let mut t = CycleTiming::default();
    let t0 = tick(cfg.record_timing);
    let sweep = joint_needle_sweep(&*ctx, stack, spectral, pool, cfg, &mut t.refold_s)?;
    // The sweep interleaves K re-folds with K scans; the scan is the rest.
    t.scan_s = tock(t0) - t.refold_s;
    if sweep.n_sites == 0 {
        return Ok(true); // "Stack too thin for needle insertion."
    }
    let t0 = tick(cfg.record_timing);
    let current_mf = ctx.evaluate_merit(stack)?;
    t.merit_s = tock(t0);

    // The gain shift is summed over the K folds for the same reason P is:
    // growing the seed grows every environment's equivalent-medium
    // reference, and the joint merit charges all of them.
    // (site, material, P if scored by P, seed thickness, improvement it promises)
    let mut best: Option<(JointSite, usize, Option<f64>, f64, f64)> = sweep.best.map(|c| {
        let predicted = -(c.p + sweep.gtot) * cfg.needle_seed_thickness_nm;
        (
            c.site,
            c.material,
            Some(c.p),
            cfg.needle_seed_thickness_nm,
            predicted,
        )
    });

    // Rough boundary loci: the joint merit with the needle in place in
    // every copy, at its start thickness. The trial runs on a copy of the
    // compile, swapped into the context for the one evaluation. The trials
    // count to the scan's time.
    let t0 = tick(cfg.record_timing);
    for &(gap, m, d0) in &sweep.rough {
        let lost = || "joint needle: the context lost its environments".to_string();
        let film = new_needle(stack.films(), &pool[m], d0);
        let mut trial = stack.clone();
        let mut envs = ctx.environments().ok_or_else(lost)?.clone();
        envs.boundary_design(&mut trial, gap, &film)?;
        envs.insert_boundary(gap, &film)?;
        let saved = std::mem::replace(ctx.environments_mut().ok_or_else(lost)?, envs);
        let merit = ctx.evaluate_merit(&trial);
        *ctx.environments_mut().ok_or_else(lost)? = saved;
        let gain = current_mf - merit?;
        if best.as_ref().is_none_or(|b| gain > b.4) {
            best = Some((JointSite::Boundary(gap), m, None, d0, gain));
        }
    }
    t.scan_s += tock(t0);

    let timing = cfg.record_timing.then_some(t);
    let stalled = |history: &mut Vec<NeedleCycleResult>, p, predicted| {
        history.push(NeedleCycleResult {
            cycle: cycle + 1,
            merit_before: current_mf,
            merit_after: current_mf,
            best_p: p,
            predicted_improvement: predicted,
            layer_count: stack.films().len(),
            insertion: None,
            timing,
        });
    };

    let Some((site, m, p, seed_nm, predicted)) = best else {
        stalled(history, None, None);
        return Ok(true); // no improving locus in any environment
    };
    if predicted < cfg.convergence_threshold {
        stalled(history, p, Some(predicted));
        return Ok(true); // "Convergence reached - stopping."
    }

    let lost = || "joint needle: the context lost its environments".to_string();
    let mat = &pool[m];
    let seed = new_needle(stack.films(), mat, seed_nm);
    // Locus -> the shared design's own row: the host's before a split, the
    // new film's after a boundary insertion. `host_row` answers the BULK
    // row, so an interface-carrying design film is still a legal host (N1).
    let row_of = |ctx: &C, slot: usize| {
        ctx.environments()
            .and_then(|e| e.host_row(0, slot))
            .ok_or_else(|| {
                format!("joint needle: design parameter {slot} has no row in environment 0")
            })
    };
    // Every copy in the design object, then every template.
    let (row0, depth, boundary) = match site {
        JointSite::Interior {
            slot,
            depth_into_layer_nm,
        } => {
            let row0 = row_of(&*ctx, slot)?;
            ctx.environments().ok_or_else(lost)?.split_design(
                stack,
                slot,
                depth_into_layer_nm,
                &seed,
            )?;
            ctx.environments_mut().ok_or_else(lost)?.insert_seed(
                slot,
                depth_into_layer_nm,
                &seed,
            )?;
            (row0, depth_into_layer_nm, false)
        }
        JointSite::Boundary(gap) => {
            ctx.environments()
                .ok_or_else(lost)?
                .boundary_design(stack, gap, &seed)?;
            let slot = ctx
                .environments_mut()
                .ok_or_else(lost)?
                .insert_boundary(gap, &seed)?;
            (row_of(&*ctx, slot)?, 0.0, true)
        }
    };

    let t0 = tick(cfg.record_timing);
    let new_mf = ctx.optimize_thicknesses(stack)?;
    t.lm_s = tock(t0);
    history.push(NeedleCycleResult {
        cycle: cycle + 1,
        merit_before: current_mf,
        merit_after: new_mf,
        best_p: p,
        predicted_improvement: Some(predicted),
        layer_count: stack.films().len(),
        insertion: Some(Insertion {
            film_idx: row0,
            depth_into_layer_nm: depth,
            material: mat.material.clone(),
            boundary,
            seed_nm,
        }),
        timing: cfg.record_timing.then_some(t),
    });
    Ok(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[allow(unused_imports)]
    use crate::smatrix::angle::{Radians, deg, degs};
    // Only the tests build NeedleTargets by hand; the pass returns them.
    use crate::smatrix::synthesis::merit::{MeritSpec, SimCurves};
    use crate::smatrix::synthesis::needle_pass::NeedleTargets;

    /// Frozen context: merit constant, optimization a no-op (keeps the
    /// analytic AR condition intact so the scan sees exact zeros).
    struct StillCtx;
    impl DesignContext for StillCtx {
        fn evaluate_merit(&self, _s: &DesignStack) -> Result<f64, String> {
            Ok(0.0)
        }
        fn simulate(&self, _s: &DesignStack) -> Result<SimCurves, String> {
            Err("mock context has no simulator".into())
        }
        fn optimize_thicknesses(&mut self, s: &mut DesignStack) -> Result<f64, String> {
            self.evaluate_merit(s)
        }
    }

    /// air | G(200 nm, n = 1.52, host) | glass: all matched, so R is
    /// exactly the single-interface value ((n−1)/(n+1))² ≈ 0.04258 and
    /// T = 1 − R (lossless) at normal incidence.
    /// (A lossless R = 0 point would NOT discriminate: T = 1 − R
    /// identically there, so every intensity gradient vanishes with R's.)
    fn glass_stack() -> (DesignStack, f64) {
        let nw = 1;
        let r0 = ((1.52_f64 - 1.0) / (1.52_f64 + 1.0)).powi(2);
        let mut ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
        ambient.optimize = false;
        ambient.needle = false;
        let mut substrate = LayerSpec::constant("sub", 1.52, 0.0, 0.0, nw);
        substrate.optimize = false;
        substrate.needle = false;
        let stack = DesignStack::with_films(
            ambient,
            substrate,
            vec![LayerSpec::constant("G", 1.52, 0.0, 200.0, nw)],
        )
        .unwrap();
        (stack, r0)
    }

    fn fold_with(bucket: &str, target: f64) -> NeedleTargets {
        let zero = || (vec![0.0f64; 1], vec![0.0f64; 1]);
        let demand = || (vec![target; 1], vec![1.0f64; 1]);
        let (r, t) = match bucket {
            "r" => (demand(), zero()),
            "t" => (zero(), demand()),
            _ => unreachable!(),
        };
        NeedleTargets {
            r,
            t,
            a: zero(),
            rb: zero(),
            tb: zero(),
            ab: zero(),
            phi: [zero(), zero(), zero(), zero()],
            phi_gain_shift: [0.0; 4],
            grad_r: vec![0.0f64; 1],
            grad_t: vec![0.0f64; 1],
        }
    }

    fn entry(label: &str, n: f64) -> LayerSpec {
        pool_entry(label, vec![Complex64::new(n, 0.0)].into())
    }

    fn pool_h() -> MaterialPool {
        vec![entry("H", 2.35)].into()
    }

    fn cycle_cfg() -> NeedleCycleConfig {
        NeedleCycleConfig {
            max_needles: 2,
            convergence_threshold: 1e-4,
            needle_seed_thickness_nm: 5.0,
            scan_step_nm: 2.0,
            refold_per_cycle: true,
            record_timing: false,
            boundary_needles: false,
        }
    }

    /// SpectralInputs wrapping a hand-built fold (empty spec — the
    /// StillCtx simulator errors, so the static fold stands).
    fn spectral_of(fold: NeedleTargets) -> SpectralInputs {
        SpectralInputs {
            wavls: vec![1000.0],
            sin_theta: vec![0.0],
            folds: vec![fold],
            spec: MeritSpec::new(),
            angles: vec![Radians::ZERO],
        }
    }

    /// Rows with needle=false are never insertion sites, whatever the
    /// pool holds. The run breaks with the stack untouched.
    #[test]
    fn non_host_rows_never_seed() {
        let nw = 1;
        let mut ambient = LayerSpec::constant("air", 1.0, 0.0, 0.0, nw);
        ambient.optimize = false;
        ambient.needle = false;
        let mut substrate = LayerSpec::constant("sub", 1.52, 0.0, 0.0, nw);
        substrate.optimize = false;
        substrate.needle = false;
        let mut host = LayerSpec::constant("G", 1.52, 0.0, 200.0, nw);
        host.needle = false; // the pool offers H, the flag refuses
        let mut stack = DesignStack::with_films(ambient, substrate, vec![host]).unwrap();
        let mut ctx = StillCtx;
        let spectral = spectral_of(fold_with("t", 2.0));
        let hist =
            run_needle_cycles(&mut ctx, &mut stack, &spectral, &pool_h(), &cycle_cfg()).unwrap();
        assert!(hist.is_empty());
        assert_eq!(stack.films().len(), 1);
        assert!((stack.films()[0].d_nm - 200.0).abs() < 1e-12);
    }

    /// A layer never accepts a needle of its own index. The default pool
    /// of a one-material stack is that material, so nothing is admissible
    /// and the run breaks with the stack untouched.
    #[test]
    fn same_index_needles_are_never_scored() {
        let (mut stack, _) = glass_stack();
        let pool = pool_from_layers(&stack, |_| true);
        assert_eq!(pool.len(), 1);
        assert_eq!(pool[0].material.as_ref(), "G");
        let mut ctx = StillCtx;
        let spectral = spectral_of(fold_with("t", 0.0));
        let hist = run_needle_cycles(&mut ctx, &mut stack, &spectral, &pool, &cycle_cfg()).unwrap();
        assert!(hist.is_empty());
        assert_eq!(stack.films().len(), 1);
    }

    /// Every pool material is a candidate in every layer of another
    /// index, needles included: a needle is a layer of its own, under a
    /// fresh name. No row is ever air (the old contrast lookup could fall
    /// back to n = 1 when a site in a film without an entry won).
    #[test]
    fn needles_are_layers_with_fresh_names() {
        let (mut stack, _) = glass_stack();
        let pool: MaterialPool = vec![entry("H", 2.35), entry("L", 1.38)].into();
        let mut ctx = StillCtx;
        let spectral = spectral_of(fold_with("t", 0.0));
        let cfg = NeedleCycleConfig {
            max_needles: 6,
            ..cycle_cfg()
        };
        let hist = run_needle_cycles(&mut ctx, &mut stack, &spectral, &pool, &cfg).unwrap();
        let n_ins = hist.iter().filter(|h| h.insertion.is_some()).count();
        assert!(n_ins >= 2, "{n_ins} insertions: degenerate");
        assert_eq!(stack.films().len(), 1 + 2 * n_ins);

        let names: Vec<&str> = stack.films().iter().map(|l| l.material.as_ref()).collect();
        let mut unique = names.clone();
        unique.sort_unstable();
        unique.dedup();
        assert_eq!(unique.len(), names.len() - n_ins, "{names:?}"); // a split host keeps its name
        for h in hist.iter().filter_map(|h| h.insertion.as_ref()) {
            assert!(["H", "L"].contains(&h.material.as_ref()));
        }
        for l in stack.films().iter().filter(|l| l.material.contains('#')) {
            assert!(l.needle && l.optimize);
            let label = l.material.split('#').next().unwrap();
            assert_eq!(
                l.nk,
                pool.iter().find(|m| &*m.material == label).unwrap().nk
            );
        }
        assert!(stack.films().iter().all(|l| l.nk[0].re > 1.0));
    }

    #[test]
    fn needle_names_skip_taken_ones() {
        let nw = 1;
        let films = vec![
            LayerSpec::constant("H#1", 2.35, 0.0, 5.0, nw),
            LayerSpec::constant("H#3", 2.35, 0.0, 5.0, nw),
        ];
        assert_eq!(needle_name(&films, "H").as_ref(), "H#2");
        assert_eq!(needle_name(&films, "L").as_ref(), "L#1");
    }

    fn rules(pairs: &[(&str, &[&str])]) -> InsertRules {
        pairs
            .iter()
            .map(|(h, ns)| (h.to_string(), ns.iter().map(|n| n.to_string()).collect()))
            .collect()
    }

    /// air | L(30) H(3) G(30) | glass (1.52).
    fn boundary_stack() -> DesignStack {
        let (s, _) = glass_stack();
        let films = vec![
            LayerSpec::constant("L", 1.38, 0.0, 30.0, 1),
            LayerSpec::constant("H", 2.35, 0.0, 3.0, 1),
            LayerSpec::constant("G", 1.52, 0.0, 30.0, 1),
        ];
        DesignStack::with_films(s.ambient().clone(), s.substrate().clone(), films).unwrap()
    }

    /// The labels a boundary above row `b` accepts.
    fn on_boundary(pool: &MaterialPool, stack: &DesignStack, b: usize) -> Vec<String> {
        (0..pool.len())
            .filter(|&m| pool.boundary_accepts(stack, b, m).unwrap())
            .map(|m| pool[m].material.to_string())
            .collect()
    }

    fn boundary_pool() -> MaterialPool {
        vec![
            entry("L", 1.38),
            entry("H", 2.35),
            entry("G", 1.52),
            entry("A", 1.0),
            entry("M", 1.7),
        ]
        .into()
    }

    /// A boundary needle differs in index from both neighbours, the ambient
    /// and the substrate included.
    #[test]
    fn a_boundary_needle_differs_from_both_neighbours() {
        let (s, pool) = (boundary_stack(), boundary_pool());
        assert_eq!(on_boundary(&pool, &s, 0), ["H", "G", "M"]); // air | L
        assert_eq!(on_boundary(&pool, &s, 1), ["G", "A", "M"]); // L | H
        assert_eq!(on_boundary(&pool, &s, 3), ["L", "H", "A", "M"]); // G | glass
    }

    /// Every neighbouring host must accept the needle: a closed host
    /// closes its boundaries too.
    #[test]
    fn a_boundary_needle_needs_both_hosts_rules() {
        let s = boundary_stack();
        let pool = boundary_pool()
            .with_insert_rules(&rules(&[("H", &[]), ("L", &["M", "A"])]))
            .unwrap();
        assert_eq!(on_boundary(&pool, &s, 0), ["M"]); // A is the ambient's index
        assert!(on_boundary(&pool, &s, 1).is_empty()); // H closed
        assert!(on_boundary(&pool, &s, 2).is_empty());
        assert_eq!(on_boundary(&pool, &s, 3), ["L", "H", "A", "M"]); // G: default
    }

    /// On an outer edge the needle becomes the edge layer: the list decides.
    #[test]
    fn a_boundary_needle_on_an_edge_is_on_its_list() {
        let mut s = boundary_stack();
        let pool = boundary_pool();
        let labels = |l: &[&str]| l.iter().map(|x| x.to_string()).collect::<Vec<_>>();
        let edges = pool
            .edge_materials(Some(&labels(&["L", "M"])), Some(&labels(&["G"])))
            .unwrap();
        s.set_edge_materials(edges).unwrap();
        assert_eq!(on_boundary(&pool, &s, 0), ["M"]);
        assert!(on_boundary(&pool, &s, 3).is_empty()); // G is the only inner, and the edge already
        assert_eq!(on_boundary(&pool, &s, 1), ["G", "A", "M"]); // inside: no list

        // A site that accepts nothing is in no material's site list, so
        // neither the scan nor a trial ever visits it.
        let sites = build_boundary_sites(s.films(), s.spans());
        assert!(sites.iter().any(|x| x.film_idx == 3));
        for m in 0..pool.len() {
            let kept = pool.sites_for(m, &s, &sites).unwrap();
            assert!(kept.iter().all(|x| x.film_idx != 3), "{}", pool[m].material);
        }
    }

    /// A film too thin for an interior site is reached only through its
    /// boundaries: off, the pass has no site; on, every needle is a new
    /// film on a boundary, and no two touching films share an index.
    #[test]
    fn boundary_needles_reach_films_too_thin_to_split() {
        let thin = || {
            let (s, _) = glass_stack();
            let g = LayerSpec::constant("G", 1.52, 0.0, 3.0, 1);
            DesignStack::with_films(s.ambient().clone(), s.substrate().clone(), vec![g]).unwrap()
        };
        // H only: an H needle never takes an H needle, so every later
        // site is a boundary too.
        let pool = pool_h();
        let spectral = spectral_of(fold_with("t", 0.0));

        let mut stack = thin();
        let hist =
            run_needle_cycles(&mut StillCtx, &mut stack, &spectral, &pool, &cycle_cfg()).unwrap();
        assert!(hist.is_empty());

        let mut stack = thin();
        let cfg = NeedleCycleConfig {
            boundary_needles: true,
            max_needles: 3,
            ..cycle_cfg()
        };
        let hist = run_needle_cycles(&mut StillCtx, &mut stack, &spectral, &pool, &cfg).unwrap();
        let ins: Vec<&Insertion> = hist.iter().filter_map(|h| h.insertion.as_ref()).collect();
        assert!(!ins.is_empty(), "nothing landed on a boundary");
        assert!(
            ins.iter()
                .all(|i| i.boundary && i.depth_into_layer_nm == 0.0)
        );
        let f = stack.films();
        assert_eq!(f.len(), 1 + ins.len());
        assert!(f.iter().any(|l| &*l.material == "G" && l.d_nm == 3.0));
        assert!(
            f.windows(2).all(|w| w[0].nk != w[1].nk),
            "touching films share an index"
        );
    }

    /// Every film lowers the merit by 0.1: a trial with a needle in place
    /// sees a gain, whatever P says.
    struct CountCtx;
    impl DesignContext for CountCtx {
        fn evaluate_merit(&self, s: &DesignStack) -> Result<f64, String> {
            Ok(1.0 - 0.1 * s.films().len() as f64)
        }
        fn simulate(&self, _s: &DesignStack) -> Result<SimCurves, String> {
            Err("mock context has no simulator".into())
        }
        fn optimize_thicknesses(&mut self, s: &mut DesignStack) -> Result<f64, String> {
            self.evaluate_merit(s)
        }
    }

    /// A 3 nm film with rough faces: both its boundaries are rough, so the
    /// needle is trial-scored and goes in at `max(seed, cap * sigma_eff)`.
    /// sigma 4 is capped to 1.5 by the 3 nm film: d0 = 3 at cap 2; with the
    /// cap off, sigma stays 4 and d0 is the seed.
    #[test]
    fn a_rough_boundary_needle_starts_at_its_start_thickness() {
        let rough = |cap: f64| {
            let (s, _) = glass_stack();
            let mut g = LayerSpec::constant("G", 1.52, 0.0, 3.0, 1);
            (g.rough_type, g.rough_val) = (5, 4.0);
            let mut sub = s.substrate().clone();
            (sub.rough_type, sub.rough_val) = (5, 4.0);
            let mut st = DesignStack::with_films(s.ambient().clone(), sub, vec![g]).unwrap();
            st.set_roughness_cap(cap).unwrap();
            st
        };
        let spectral = spectral_of(fold_with("t", 0.0));
        let cfg = NeedleCycleConfig {
            boundary_needles: true,
            max_needles: 1,
            needle_seed_thickness_nm: 2.0,
            ..cycle_cfg()
        };
        for (cap, d0) in [(2.0, 3.0), (0.0, 2.0)] {
            let mut stack = rough(cap);
            let hist =
                run_needle_cycles(&mut CountCtx, &mut stack, &spectral, &pool_h(), &cfg).unwrap();
            let h = &hist[0];
            let ins = h.insertion.as_ref().expect("no rough boundary needle");
            assert!(ins.boundary);
            assert_eq!(ins.seed_nm, d0);
            assert_eq!(h.best_p, None);
            assert!((h.predicted_improvement.unwrap() - 0.1).abs() < 1e-12);
            assert_eq!(stack.films()[ins.film_idx].d_nm, d0);
            // The rough interface stays with the film below.
            assert_eq!(stack.films()[ins.film_idx].rough_val, 0.0);
            assert_eq!(stack.films()[ins.film_idx + 1].rough_val, 4.0);
        }
    }

    /// The glass host G takes only L, and L takes nothing: every needle is
    /// L, and H, though in the pool, never lands.
    #[test]
    fn insert_rules_narrow_the_scan() {
        let (mut stack, _) = glass_stack();
        let pool = MaterialPool::from(vec![entry("G", 1.52), entry("H", 2.35), entry("L", 1.38)])
            .with_insert_rules(&rules(&[("G", &["L"]), ("L", &[])]))
            .unwrap();
        let mut ctx = StillCtx;
        let spectral = spectral_of(fold_with("t", 0.0));
        let cfg = NeedleCycleConfig {
            max_needles: 4,
            ..cycle_cfg()
        };
        let hist = run_needle_cycles(&mut ctx, &mut stack, &spectral, &pool, &cfg).unwrap();
        let ins: Vec<&str> = hist
            .iter()
            .filter_map(|h| h.insertion.as_ref().map(|i| i.material.as_ref()))
            .collect();
        assert!(!ins.is_empty(), "the rule left L admissible in G");
        assert!(ins.iter().all(|m| *m == "L"), "{ins:?}");
        assert!(stack.films().iter().all(|l| l.nk[0].re != 2.35));
    }

    /// An empty list closes a host: with G closed there is no site at all.
    #[test]
    fn an_empty_rule_closes_a_host() {
        let (mut stack, _) = glass_stack();
        let pool = MaterialPool::from(vec![entry("G", 1.52), entry("H", 2.35)])
            .with_insert_rules(&rules(&[("G", &[])]))
            .unwrap();
        let mut ctx = StillCtx;
        let spectral = spectral_of(fold_with("t", 0.0));
        let hist = run_needle_cycles(&mut ctx, &mut stack, &spectral, &pool, &cycle_cfg()).unwrap();
        assert!(hist.is_empty());
        assert_eq!(stack.films().len(), 1);
    }

    /// Under rules, a host whose index is in no pool entry needs the "*"
    /// rule; without it the run refuses rather than guessing.
    #[test]
    fn a_film_outside_the_pool_needs_the_star_rule() {
        let pool = || MaterialPool::from(vec![entry("H", 2.35), entry("L", 1.38)]);
        let mut ctx = StillCtx;
        let spectral = spectral_of(fold_with("t", 0.0));

        let (mut stack, _) = glass_stack();
        let strict = pool().with_insert_rules(&rules(&[("H", &["L"])])).unwrap();
        let e =
            run_needle_cycles(&mut ctx, &mut stack, &spectral, &strict, &cycle_cfg()).unwrap_err();
        assert!(e.contains("\"*\"") && e.contains("\"G\""), "{e}");

        let (mut stack, _) = glass_stack();
        let star = pool()
            .with_insert_rules(&rules(&[("*", &["H"]), ("H", &[])]))
            .unwrap();
        let hist = run_needle_cycles(&mut ctx, &mut stack, &spectral, &star, &cycle_cfg()).unwrap();
        assert!(
            hist.iter()
                .filter_map(|h| h.insertion.as_ref())
                .all(|i| i.material.as_ref() == "H")
        );
    }

    #[test]
    fn insert_rules_refuse_bad_labels() {
        let pool = || MaterialPool::from(vec![entry("H", 2.35), entry("L", 1.38)]);
        let e = pool()
            .with_insert_rules(&rules(&[("X", &["L"])]))
            .unwrap_err();
        assert!(e.contains("\"X\" is not a pool material"), "{e}");
        let e = pool()
            .with_insert_rules(&rules(&[("H", &["Y"])]))
            .unwrap_err();
        assert!(e.contains("\"Y\" is not a pool material"), "{e}");
        let e = pool()
            .with_insert_rules(&rules(&[("H", &["H"])]))
            .unwrap_err();
        assert!(e.contains("lists itself"), "{e}");
        // Empty rules are no rules: the pool is unchanged.
        let p = pool().with_insert_rules(&InsertRules::new()).unwrap();
        assert!(p.allowed.is_none());
    }

    /// The edge lists resolve pool labels to their indices, and the stack
    /// refuses an edge already outside its list.
    #[test]
    fn needle_rules_set_the_edges() {
        let (mut stack, _) = glass_stack();
        let pool = MaterialPool::from(vec![entry("G", 1.52), entry("H", 2.35)]);
        let bad = NeedleRules {
            outer: Some(vec!["H".into()]),
            ..Default::default()
        };
        let e = bad.apply(&mut stack, pool.clone()).unwrap_err();
        assert!(e.contains("outer") && e.contains("\"H\""), "{e}");
        assert!(stack.edge_materials().is_free());

        let good = NeedleRules {
            outer: Some(vec!["G".into(), "H".into()]),
            inner: Some(vec!["G".into()]),
            ..Default::default()
        };
        good.apply(&mut stack, pool).unwrap();
        let edges = stack.edge_materials();
        assert_eq!(edges.outer.as_ref().unwrap().len(), 2);
        assert_eq!(&*edges.inner.as_ref().unwrap()[0].0, "G");
    }

    /// `seed_edges` puts a film of the list's first material on each edge
    /// off its list and keeps the old layer behind it; an edge already on
    /// its list is left alone, and the thickness must be positive.
    #[test]
    fn seed_edges_put_the_first_listed_material_on_an_edge_off_its_list() {
        let pool = MaterialPool::from(vec![entry("G", 1.52), entry("H", 2.35), entry("L", 1.38)]);
        let rules = |d| NeedleRules {
            outer: Some(vec!["H".into(), "L".into()]),
            inner: Some(vec!["G".into()]),
            seed_edges: Some(d),
            ..Default::default()
        };
        let (mut stack, _) = glass_stack();
        rules(2.0).apply(&mut stack, pool.clone()).unwrap();
        let f: Vec<(&str, f64)> = stack
            .films()
            .iter()
            .map(|l| (l.material.as_ref(), l.d_nm))
            .collect();
        assert_eq!(f, [("H#1", 2.0), ("G", 200.0)]);
        assert!(stack.films()[0].needle && stack.films()[0].optimize);
        assert!(stack.edge_materials().outer.is_some());

        for d in [0.0, -1.0, f64::NAN] {
            let (mut stack, _) = glass_stack();
            let e = rules(d).apply(&mut stack, pool.clone()).unwrap_err();
            assert!(e.contains("seed_edges"), "{e}");
        }
    }

    #[test]
    fn satisfied_r_demand_inserts_nothing() {
        // R demanded at its exact value → residual ~1e-16 → dust profile,
        // convergence gate stops: no insertion, film count untouched.
        let (mut stack, r0) = glass_stack();
        let mut ctx = StillCtx;
        let spectral = spectral_of(fold_with("r", r0));
        let hist =
            run_needle_cycles(&mut ctx, &mut stack, &spectral, &pool_h(), &cycle_cfg()).unwrap();
        assert_eq!(hist.len(), 1);
        assert!(hist[0].insertion.is_none());
        assert!(hist[0].best_p.map(|p| p > -1e-9).unwrap_or(true));
        assert_eq!(stack.films().len(), 1);
    }

    /// Real solver context over one Rs demand at 1000 nm, normal incidence.
    fn rs_ctx(
        kind: crate::smatrix::synthesis::merit::ConstraintKind,
        target: f64,
    ) -> (
        crate::smatrix::synthesis::evaluator::SmatrixContext,
        SpectralInputs,
    ) {
        use crate::smatrix::synthesis::evaluator::SmatrixContext;
        use crate::smatrix::synthesis::merit::{MeritKey, MeritTarget, SimTransform};
        use crate::smatrix::synthesis::thick_opt::LmConfig;

        let mut spec = MeritSpec::new();
        let k = spec.add_key(MeritKey {
            angle: Radians::ZERO,
            curve: crate::smatrix::synthesis::merit::CurveId::Rs,
        });
        spec.add_target(MeritTarget {
            key_idx: k as u32,
            env_idx: 0,
            wavelengths: vec![1000.0].into(),
            kind,
            transform: SimTransform::Linear,
            norm_factor: 1.0,
            normalized_targets: vec![target].into(),
            tolerances: vec![0.05].into(),
            band: vec![].into(),
            phase: false,
            differential_passes: None,
            integral: false,
            weight: 1.0,
            count_norm: None,
        })
        .unwrap();
        let ctx = SmatrixContext {
            wavls: vec![1000.0],
            sin_theta: vec![0.0],
            spec: spec.clone(),
            clamp_min_nm: 2.0,
            clamp_max_nm: 1000.0,
            lm: LmConfig::default(),
            clamp_accumulator: crate::smatrix::synthesis::structure::ClampReport::default(),
            envs: None,
            thin_layer_policy: crate::smatrix::synthesis::config::ThinLayerPolicy::Remove,
        };
        let spectral = SpectralInputs::from_spec(&spec, &[Radians::ZERO], &[1000.0]).unwrap();
        (ctx, spectral)
    }

    #[test]
    fn live_refold_masks_satisfied_above_demand() {
        // Glass R ≈ 0.0426 against Above 0.01: satisfied (sim ≥ target)
        // at the operating point. The static conservative fold would
        // insert (Above folds active without a sim); the per-cycle live
        // re-fold masks it, so nothing is inserted. Real solver context.
        use crate::smatrix::synthesis::merit::ConstraintKind;
        let (mut ctx, spectral) = rs_ctx(ConstraintKind::Above, 0.01);
        let (mut stack, _) = glass_stack();
        let hist =
            run_needle_cycles(&mut ctx, &mut stack, &spectral, &pool_h(), &cycle_cfg()).unwrap();
        assert!(hist.iter().all(|h| h.insertion.is_none()));
        assert_eq!(stack.films().len(), 1);
    }

    /// The timers observe and never steer: the same pass with
    /// `record_timing` on ends in the same stack, bit for bit, with the
    /// same history — and only the timed run carries timings.
    #[test]
    fn record_timing_changes_nothing_but_the_timing() {
        use crate::smatrix::synthesis::merit::ConstraintKind;
        let run = |record_timing: bool| {
            // Glass R ≈ 0.0426 against R = 0: violated, so the pass
            // inserts and re-optimizes through the real LM.
            let (mut ctx, spectral) = rs_ctx(ConstraintKind::Exact, 0.0);
            let (mut stack, _) = glass_stack();
            let cfg = NeedleCycleConfig {
                record_timing,
                ..cycle_cfg()
            };
            let hist = run_needle_cycles(&mut ctx, &mut stack, &spectral, &pool_h(), &cfg).unwrap();
            (hist, stack)
        };
        let (h_off, s_off) = run(false);
        let (h_on, s_on) = run(true);

        assert!(
            h_off.iter().any(|h| h.insertion.is_some()),
            "no insertion: degenerate"
        );
        let bits =
            |s: &DesignStack| -> Vec<u64> { s.films().iter().map(|l| l.d_nm.to_bits()).collect() };
        assert_eq!(bits(&s_off), bits(&s_on));
        assert_eq!(h_off.len(), h_on.len());
        for (a, b) in h_off.iter().zip(&h_on) {
            assert_eq!(a.merit_before.to_bits(), b.merit_before.to_bits());
            assert_eq!(a.merit_after.to_bits(), b.merit_after.to_bits());
            assert_eq!(a.best_p.map(f64::to_bits), b.best_p.map(f64::to_bits));
            assert_eq!(a.layer_count, b.layer_count);
            assert_eq!(
                a.insertion
                    .as_ref()
                    .map(|i| (i.film_idx, i.depth_into_layer_nm.to_bits())),
                b.insertion
                    .as_ref()
                    .map(|i| (i.film_idx, i.depth_into_layer_nm.to_bits())),
            );
            assert!(a.timing.is_none());
            let t = b.timing.expect("timed run must carry timings");
            for v in [t.refold_s, t.scan_s, t.merit_s, t.lm_s, t.initial_lm_s] {
                assert!(v.is_finite() && v >= 0.0, "{t:?}");
            }
            if b.insertion.is_some() {
                assert!(t.scan_s > 0.0 && t.lm_s > 0.0, "{t:?}");
            }
        }
        assert!(h_on[0].timing.unwrap().initial_lm_s > 0.0);
        assert!(
            h_on[1..]
                .iter()
                .all(|h| h.timing.unwrap().initial_lm_s == 0.0)
        );
    }

    #[test]
    fn violated_t_demand_drives_insertion() {
        // Same stack, T ≈ 0.957 against a T = 0 demand: inserting H
        // disrupts the matching and lowers T, so the scan must find a
        // negative-P site and the host film must split.
        let (mut stack, _) = glass_stack();
        let mut ctx = StillCtx;
        let spectral = spectral_of(fold_with("t", 0.0));
        let hist =
            run_needle_cycles(&mut ctx, &mut stack, &spectral, &pool_h(), &cycle_cfg()).unwrap();
        // max_needles = 2 and T stays violated → two insertions
        // (1 → 3 → 5 films), both H into the original host lineage.
        assert_eq!(hist.len(), 2);
        for h in &hist {
            let ins = h.insertion.as_ref().expect("expected an insertion");
            assert!(h.best_p.unwrap() < 0.0);
            assert_eq!(ins.material.as_ref(), "H");
        }
        assert_eq!(stack.films().len(), 5);
    }
}
