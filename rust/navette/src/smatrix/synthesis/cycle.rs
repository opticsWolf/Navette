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
use crate::smatrix::synthesis::needle_pass::{
    NeedlePassInput, ScanSite, build_needle_targets_env, build_scan_sites, needle_pass_scan,
};
use crate::smatrix::synthesis::pipeline::SpectralInputs;
use crate::smatrix::synthesis::structure::{DesignStack, EdgeMaterial, EdgeMaterials, LayerSpec};

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

#[derive(Clone, Debug)]
pub struct Insertion {
    pub film_idx: usize,
    pub depth_into_layer_nm: f64,
    pub material: Arc<str>,
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

    /// The `sites` whose film accepts a needle of entry `m`. Filtering
    /// whole rows keeps each surviving row's steps intact.
    fn sites_for(
        &self,
        m: usize,
        films: &[LayerSpec],
        sites: &[ScanSite],
    ) -> Result<Vec<ScanSite>, String> {
        let mut out = Vec::new();
        for s in sites {
            if self.accepts(&films[s.film_idx], m)? {
                out.push(*s);
            }
        }
        Ok(out)
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
}

impl NeedleRules {
    /// Apply to a run: the pool takes the insertion rules, the stack the
    /// edge materials (its current edges must comply).
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
            stack.set_edge_materials(edges)?;
        }
        Ok(pool)
    }
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
fn new_needle(films: &[LayerSpec], m: &LayerSpec, d_nm: f64) -> LayerSpec {
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
    /// Index into `CompiledEnvironments::slots`.
    slot: usize,
    depth_into_layer_nm: f64,
    /// The needle material: the index into the pool whose sweep won.
    material: usize,
    /// P summed over environments.
    p: f64,
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
/// Returns the winner (most negative ΣP, `None` when nothing improves), the
/// gain shift summed over the K folds for the convergence test, and the
/// number of admissible design sites — zero means "stack too thin".
/// Under `record_timing` the K re-folds' wall time is added to `refold_s`,
/// so the caller can split it out of the sweep's total.
fn joint_needle_sweep<C: DesignContext + ?Sized>(
    ctx: &C,
    stack: &DesignStack,
    spectral: &SpectralInputs,
    pool: &MaterialPool,
    cfg: &NeedleCycleConfig,
    refold_s: &mut f64,
) -> Result<(Option<JointCandidate>, f64, usize), String> {
    let envs = ctx
        .environments()
        .ok_or_else(|| "joint_needle_sweep: the context has no environments".to_string())?;
    let k = envs.n_envs();
    let stacks = envs.expand(stack)?;

    let mut gtot = 0.0f64;
    // (pool index, slot, step into the host) -> (depth, ΣP, environments seen)
    let mut buckets: HashMap<(usize, usize, usize), (f64, f64, usize)> = HashMap::new();
    let mut n_sites = 0usize;

    for (e, st) in stacks.iter().enumerate() {
        let mut fold = spectral.folds[e].clone();
        let t0 = tick(cfg.record_timing);
        if cfg.refold_per_cycle
            && let Ok(sim) = ctx.simulate(st)
            && let Ok(f) = build_needle_targets_env(
                &spectral.spec,
                &spectral.angles_deg,
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
            let all = build_scan_sites(st.films(), st.spans(), cfg.scan_step_nm);
            let sites = pool.sites_for(mi, st.films(), &all)?;
            let res = needle_pass_scan(&input, &sites)?;
            // `build_scan_sites` walks rows in order and emits steps
            // 1..n_steps inside each admissible row, so the running count
            // per row IS the step index — and it is the same index in
            // every environment, because a design row has the same
            // thickness everywhere by construction (`expand` copies it).
            let mut step: HashMap<usize, usize> = HashMap::new();
            for (site, &p) in res.sites.iter().zip(res.p_profile.iter()) {
                let n = step.entry(site.film_idx).or_insert(0);
                let ordinal = *n;
                *n += 1;
                let Some(slot) = envs.slot_of_row(e, site.film_idx) else {
                    continue; // surroundings: inadmissible, and never shared
                };
                if e == 0 {
                    n_sites += 1;
                }
                let entry = buckets.entry((mi, slot, ordinal)).or_insert((
                    site.depth_into_layer_nm,
                    0.0,
                    0,
                ));
                entry.1 += p;
                entry.2 += 1;
            }
        }
    }

    // §4.4's alignment assert, in its span-aware form: the same design
    // object must produce the same scan grid in every environment, so every
    // bucket is filled exactly K times. A short bucket means the assemblies
    // drifted apart, and summing it would sum over different loci.
    if let Some(((m, slot, ord), (_, _, n))) = buckets.iter().find(|(_, v)| v.2 != k) {
        let m = &pool[*m].material;
        return Err(format!(
            "joint needle sweep: the locus (slot {slot}, step {ord}, needle {m:?}) \
             appears in {n} of {k} environments - the assemblies no longer \
             share the design's scan grid"
        ));
    }

    // Deterministic order before the min: a hash map's iteration order is
    // not stable, and two loci can tie on P exactly (identical
    // environments, mirrored stacks).
    let mut items: Vec<((usize, usize, usize), (f64, f64, usize))> = buckets.into_iter().collect();
    items.sort_by_key(|a| a.0);

    let mut best: Option<JointCandidate> = None;
    for ((material, slot, _), (depth, p, _)) in items {
        if p >= 0.0 {
            continue; // not an improvement anywhere
        }
        let better = match &best {
            None => true,
            Some(b) => p < b.p,
        };
        if better {
            best = Some(JointCandidate {
                slot,
                depth_into_layer_nm: depth,
                material,
                p,
            });
        }
    }
    Ok((best, gtot, n_sites))
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
                &spectral.angles_deg,
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
        let sites = build_scan_sites(stack.films(), stack.spans(), cfg.scan_step_nm);

        // 2. Analytic sweep (both polarizations summed) once per pool
        //    material, over the sites whose layer accepts it (its index
        //    differs); the global best (site, material) pair is kept.
        //    Insertion rules narrow the pairs further.
        let per_material: Vec<Vec<_>> = (0..pool.len())
            .map(|m| pool.sites_for(m, stack.films(), &sites))
            .collect::<Result<_, _>>()?;
        if per_material.iter().all(Vec::is_empty) {
            break; // "Stack too thin for needle insertion" (or no pool).
        }
        let sa = stack.solver_arrays();

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
            Some((ins, mat, p)) => {
                // Inserting a seed of thickness δ also grows the
                // equivalent-medium reference by δ, so the true slope is
                // P(z*) + Σ gain_shift (0 for absolute-phase demand sets;
                // the site itself is unaffected, only this bookkeeping).
                let gtot: f64 = fold.phi_gain_shift.iter().sum();
                let predicted = -(p + gtot) * cfg.needle_seed_thickness_nm;
                if predicted < cfg.convergence_threshold {
                    history.push(NeedleCycleResult {
                        cycle: cycle + 1,
                        merit_before: current_mf,
                        merit_after: current_mf,
                        best_p: Some(p),
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
                let seed = new_needle(stack.films(), mat, cfg.needle_seed_thickness_nm);
                stack.insert_needle_seed(ins.film_idx, ins.depth_into_layer_nm, seed)?;

                // 5. Re-optimize.
                let t0 = tick(cfg.record_timing);
                let new_mf = ctx.optimize_thicknesses(stack)?;
                t.lm_s = tock(t0);

                NeedleCycleResult {
                    cycle: cycle + 1,
                    merit_before: current_mf,
                    merit_after: new_mf,
                    best_p: Some(p),
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
/// the compile, so every environment's template splits at its own row for
/// the same parameter. Nothing else about the insertion changes: the same
/// seed, the same depth, the same `insert_needle_seed`, the same
/// re-optimization through the context.
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
    let (best, gtot, n_sites) =
        joint_needle_sweep(&*ctx, stack, spectral, pool, cfg, &mut t.refold_s)?;
    // The sweep interleaves K re-folds with K scans; the scan is the rest.
    t.scan_s = tock(t0) - t.refold_s;
    if n_sites == 0 {
        return Ok(true); // "Stack too thin for needle insertion."
    }
    let t0 = tick(cfg.record_timing);
    let current_mf = ctx.evaluate_merit(stack)?;
    t.merit_s = tock(t0);
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

    let Some(cand) = best else {
        stalled(history, None, None);
        return Ok(true); // no improving locus in any environment
    };
    // The gain shift is summed over the K folds for the same reason P is:
    // growing the seed grows every environment's equivalent-medium
    // reference, and the joint merit charges all of them.
    let predicted = -(cand.p + gtot) * cfg.needle_seed_thickness_nm;
    if predicted < cfg.convergence_threshold {
        stalled(history, Some(cand.p), Some(predicted));
        return Ok(true); // "Convergence reached - stopping."
    }

    // Locus -> the shared design's own row. `host_row` answers the BULK
    // row, so an interface-carrying design film is still a legal host (N1).
    let row0 = ctx
        .environments()
        .and_then(|e| e.host_row(0, cand.slot))
        .ok_or_else(|| {
            format!(
                "joint needle: design parameter {} has no host row in \
                 environment 0",
                cand.slot
            )
        })?;
    let mat = &pool[cand.material];
    let seed = new_needle(stack.films(), mat, cfg.needle_seed_thickness_nm);
    stack.insert_needle_seed(row0, cand.depth_into_layer_nm, seed.clone())?;
    ctx.environments_mut()
        .ok_or_else(|| "joint needle: the context lost its environments".to_string())?
        .insert_seed(cand.slot, cand.depth_into_layer_nm, &seed)?;

    let t0 = tick(cfg.record_timing);
    let new_mf = ctx.optimize_thicknesses(stack)?;
    t.lm_s = tock(t0);
    history.push(NeedleCycleResult {
        cycle: cycle + 1,
        merit_before: current_mf,
        merit_after: new_mf,
        best_p: Some(cand.p),
        predicted_improvement: Some(predicted),
        layer_count: stack.films().len(),
        insertion: Some(Insertion {
            film_idx: row0,
            depth_into_layer_nm: cand.depth_into_layer_nm,
            material: mat.material.clone(),
        }),
        timing: cfg.record_timing.then_some(t),
    });
    Ok(false)
}

#[cfg(test)]
mod tests {
    use super::*;
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
            angles_deg: vec![0.0],
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
            angle: 0.0,
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
        let spectral = SpectralInputs::from_spec(&spec, &[0.0], &[1000.0]).unwrap();
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
