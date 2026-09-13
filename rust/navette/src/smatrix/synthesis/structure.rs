// SPDX-License-Identifier: LGPL-3.0-or-later
//! Navette -- Rust Rewrite of Numba-optimized thin-film optical solver
//!
//! synthesis::structure — DesignStack / LayerSpec.
//!
//! Port of the Loom_Structure / Layer subset used by needle synthesis
//! (loom_structure.py) plus the mutation primitives used by
//! needle_synthesis.py (`_insert_needle`, `merge_adjacent_layers`) and
//! needle_pipeline.py (`clamp_all_layers`).
//!
//! Stack model (matches Python exactly):
//!     layer_list = [ambient, film_0, …, film_{N-1}, substrate]
//! Ambient and substrate are fixed: never optimized, never removed,
//! never hosts. Only `films` is mutable.
//!
//! Solver-array layout (inherited from navette/smatrix.py):
//!   n_stack_cache — flat f64, wav-major with re/im interleaved per layer:
//!                   base = w * n_layers * 2 → [Re0, Im0, Re1, Im1, …]
//!   thicknesses   — [n_layers] nm
//!   incoherent_flags / rough_types / rough_vals — [n_layers]

use std::collections::HashMap;
use std::sync::Arc;

use num_complex::Complex64;

use crate::smatrix::optics_core::cplx;
use crate::smatrix::synthesis::config::ThinLayerPolicy;
use crate::structure::Span;

// ---------------------------------------------------------------------------
// LayerSpec
// ---------------------------------------------------------------------------

/// One layer of the design stack.
///
/// `nk` holds the complex refractive index evaluated on the *fixed*
/// simulation wavelength grid (one entry per wavelength), mirroring
/// `ArrayMaterialProvider` semantics from needle_synthesis.py where nk
/// arrays are pre-interpolated onto the grid.
#[derive(Clone, Debug)]
pub struct LayerSpec {
    pub material: Arc<str>,
    /// Complex nk per simulation wavelength, len == num_wavs.
    pub nk: Arc<[Complex64]>,
    /// Physical thickness in nm.
    pub d_nm: f64,
    /// `true` → coherent propagation (incoherent_flag = 0).
    pub coherent: bool,
    /// Roughness model id passed through to the solver untouched.
    pub rough_type: i32,
    /// Roughness amplitude (Å) passed through to the solver untouched.
    pub rough_val: f64,
    /// Included in thickness optimization.
    pub optimize: bool,
    /// Admissible host for needle insertion.
    pub needle: bool,
}

impl LayerSpec {
    /// Constant-index helper for tests and synthetic designs.
    pub fn constant(material: &str, n_re: f64, n_im: f64, d_nm: f64, num_wavs: usize) -> Self {
        LayerSpec {
            material: Arc::from(material),
            nk: vec![cplx(n_re, n_im); num_wavs].into(),
            d_nm,
            coherent: true,
            rough_type: 0,
            rough_val: 0.0,
            optimize: true,
            needle: true,
        }
    }

    /// Verbatim port of Python `Layer.clone()` semantics: an independent
    /// copy (Arc clones share the immutable nk buffer — fine, it is never
    /// mutated after construction).
    pub fn cloned(&self) -> Self {
        self.clone()
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ClampReport
// ---------------------------------------------------------------------------

/// What one clamp pass did (F0.2).
///
/// Result dicts surface it only when something happened (B4): an
/// `Option<ClampReport>` that is `None` when nothing was removed and
/// nothing was capped keeps a no-span run's serialized output
/// byte-identical. `rows_removed` / `spans_capped` are the old
/// `(n_removed, n_capped)` counts; `spans_removed` names each span the
/// floor took, because a silent correct deletion and a silent wrong one
/// look identical from the outside.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ClampReport {
    /// One description per span removed whole: "material (D nm)".
    pub spans_removed: Vec<String>,
    pub spans_capped: usize,
    pub rows_removed: usize,
}

impl ClampReport {
    pub fn is_empty(&self) -> bool {
        self.spans_removed.is_empty() && self.spans_capped == 0 && self.rows_removed == 0
    }

    /// One phase can run several clamp passes (post-cleanup, post-inflate,
    /// plus every `optimize_thicknesses` sweep in between) - the phase's
    /// report is their aggregate (the evaluator site accumulates rather
    /// than reporting per call).
    pub fn merge(&mut self, other: ClampReport) {
        self.spans_removed.extend(other.spans_removed);
        self.spans_capped += other.spans_capped;
        self.rows_removed += other.rows_removed;
    }
}

// DesignStack
// ---------------------------------------------------------------------------

/// Thin-film stack: fixed ambient + films + fixed substrate.
///
/// F0.1: `spans` rides alongside `films` — row `i` of `films` belongs to
/// exactly one span of `spans` (asserted by `assert_spans_partition` at
/// every constructor and count-changing mutator). Nothing observable of
/// any existing run reads it yet; it is the bookkeeping F0.2/F1.6/F2.3
/// consult.
///
/// F1.7: `recipes` rides alongside `spans`, aligned index-for-index —
/// span `i`'s re-emission inputs (`None` for every plain film). It rides
/// the same mutators `spans` does, and `assert_spans_partition` checks
/// both lengths at once.
#[derive(Clone, Debug)]
pub struct DesignStack {
    ambient: LayerSpec,
    substrate: LayerSpec,
    films: Vec<LayerSpec>,
    spans: Vec<Span>,
    recipes: Vec<Option<SpanRecipe>>,
    num_wavs: usize,
}

impl DesignStack {
    /// Build a stack from boundary layers plus film layers.
    ///
    /// All layers must share the same nk length (the simulation grid size).
    /// Rows not born from `expand` get one singleton span each (F0.1).
    pub fn with_films(
        ambient: LayerSpec,
        substrate: LayerSpec,
        films: Vec<LayerSpec>,
    ) -> Result<Self, String> {
        let spans: Vec<Span> = (0..films.len())
            .map(|i| Span {
                start: i,
                end: i + 1,
                logical: i,
                slice: false,
                bulk_start: i,
            })
            .collect();
        let recipes: Vec<Option<SpanRecipe>> = (0..films.len()).map(|_| None).collect();
        Self::from_parts(ambient, substrate, films, spans, recipes)
    }

    /// The one internal constructor: validated row lists plus the span
    /// partition that covers them. `from_design` passes the spans it kept
    /// from `expand` plus the recipes captured from the same emission;
    /// every other door synthesizes singletons above and one `None` per
    /// film. F1.7: `from_parts` checks the recipe vector's length against
    /// the span vector (the both-lengths check B2 asks for) at every
    /// observable point.
    fn from_parts(
        ambient: LayerSpec,
        substrate: LayerSpec,
        films: Vec<LayerSpec>,
        spans: Vec<Span>,
        recipes: Vec<Option<SpanRecipe>>,
    ) -> Result<Self, String> {
        let num_wavs = ambient.nk.len();
        if substrate.nk.len() != num_wavs {
            return Err(format!(
                "substrate nk length {} != ambient nk length {}",
                substrate.nk.len(),
                num_wavs
            ));
        }
        for (i, f) in films.iter().enumerate() {
            if f.nk.len() != num_wavs {
                return Err(format!(
                    "film {} ('{}') nk length {} != ambient nk length {}",
                    i,
                    f.material,
                    f.nk.len(),
                    num_wavs
                ));
            }
        }
        if recipes.len() != spans.len() {
            return Err(format!(
                "recipes length {} != spans length {} - the recipe vector must be \
                 aligned index-for-index with the span partition",
                recipes.len(),
                spans.len()
            ));
        }
        let s = DesignStack {
            ambient,
            substrate,
            films,
            spans,
            recipes,
            num_wavs,
        };
        s.assert_spans_partition();
        Ok(s)
    }

    // -- properties ---------------------------------------------------------

    pub fn num_wavs(&self) -> usize {
        self.num_wavs
    }

    pub fn ambient(&self) -> &LayerSpec {
        &self.ambient
    }

    pub fn substrate(&self) -> &LayerSpec {
        &self.substrate
    }

    /// The film layers (excluding ambient and substrate).
    pub fn films(&self) -> &[LayerSpec] {
        &self.films
    }

    pub fn film(&self, idx: usize) -> Option<&LayerSpec> {
        self.films.get(idx)
    }

    /// Total number of solver-visible layers including boundaries.
    pub fn total_layer_count(&self) -> usize {
        self.films.len() + 2
    }

    /// Total physical thickness of the film stack (nm).
    pub fn total_thickness_nm(&self) -> f64 {
        self.films.iter().map(|f| f.d_nm).sum()
    }

    /// The span partition (F0.1). One span per contiguous run of rows;
    /// `assert_spans_partition` holds at every observable point.
    pub fn spans(&self) -> &[Span] {
        &self.spans
    }

    /// The span containing solver row `row`.
    ///
    /// Panics on an out-of-range row: callers are the crate's own mutators
    /// and gates, all of which validate the row first — an out-of-range
    /// row here is a bookkeeping bug, not a user error.
    pub fn span_of_row(&self, row: usize) -> &Span {
        &self.spans[span_index_of_row(&self.spans, row)]
    }

    /// F1.6: the one scalable-span predicate, shared by the parameter
    /// builder (evaluator), the clamp and the needle-pipeline
    /// construction refusal: the span has at least two BULK rows and
    /// every one of them is optimize-flagged (the interface slice is
    /// derived and carries optimize=false already). This is the shape
    /// `from_design` gives a profiled (`inhomogen` or `gradient`)
    /// carrier with `optimize = true` whose mode scales exactly.
    pub(crate) fn span_is_scalable(&self, sp: &Span) -> bool {
        span_is_scalable_rows(sp, &self.films)
    }

    /// R3 / §6 row 6: the one needle-host admissibility rule, shared by
    /// `insert_needle_seed` (the public door) and the scan-site filter.
    /// A host row inside a multi-row span is refused — a needle must
    /// never split a graded profile or land on a derived row. The message
    /// names the span's material and its row range. `None` = admissible.
    /// F1.6: a SCALABLE span (an optimized profiled carrier) is also
    /// refused — scaling the profile preserves it, splitting a foreign
    /// row into the middle of it does not. The rule is unchanged in
    /// form: multi-row spans were never hosts; scalability only ever
    /// ADDS reasons a span must stay whole.
    pub(crate) fn needle_host_refusal(&self, film_idx: usize) -> Option<String> {
        if film_idx >= self.films.len() {
            return None; // range is the caller's error to raise
        }
        let sp = self.span_of_row(film_idx);
        if sp.is_singleton_bulk() {
            return None;
        }
        Some(format!(
            "insert_needle_seed: film {} ('{}') lies in solver rows {}..{}, a \
             multi-row span - not an admissible needle host",
            film_idx, self.films[film_idx].material, sp.start, sp.end
        ))
    }

    /// R2: the one named invariant check — `spans` is a contiguous,
    /// ordered, non-overlapping partition of `0..films.len()`, and every
    /// span's `bulk_start` agrees with the arithmetic predicate (B9: one
    /// source, one assertion that it matches the derivation). Called at
    /// the tail of both constructors and every count-changing mutator;
    /// F1.7 extends it to the recipe vector's length (B2's both-lengths
    /// check: `recipes` is aligned index-for-index with `spans`, checked
    /// in `from_parts` and asserted here so a mutator bug cannot hide).
    pub(crate) fn assert_spans_partition(&self) {
        assert_eq!(
            self.recipes.len(),
            self.spans.len(),
            "recipes must stay aligned index-for-index with spans"
        );
        let n = self.films.len();
        let mut next = 0usize;
        for sp in &self.spans {
            assert_eq!(
                sp.start, next,
                "spans must be contiguous and ordered (gap at row {next})"
            );
            assert!(sp.end > sp.start, "span {sp:?} is empty");
            assert!(sp.end <= n, "span {sp:?} exceeds film rows ({n})");
            assert!(
                sp.bulk_start >= sp.start && sp.bulk_start <= sp.end,
                "span {sp:?} has bulk_start outside its rows"
            );
            debug_assert_eq!(
                sp.bulk_start,
                sp.start + usize::from(sp.slice),
                "Span::bulk_start disagrees with the arithmetic predicate"
            );
            next = sp.end;
        }
        assert_eq!(
            next, n,
            "spans must cover every film row exactly once (stopped at {next} of {n})"
        );
    }

    /// Build a stack by expanding design films through `navette-structure`.
    ///
    /// Entry point for the Python driver (`stack_from_layers`): films are
    /// design layers (material NAMES resolved via `nk` on `wavelengths`),
    /// groups apply thickness/nk/roughness/interface policy, and the
    /// expansion emits solver rows. Physics (thickness/nk/coherent/
    /// roughness) comes from the rows; identity and optimizer flags come
    /// from the carrier design layer — except interface-slice rows, which
    /// are derived (optimize/needle forced false, never free parameters).
    /// Graded profiles are HOMOGENIZED with a warning per film (base
    /// index, single row): the pipeline's operators assume uniform slabs,
    /// so dropping the profile loudly beats refusing (substrate/media use)
    /// and beats flattening silently. Span-aware graded optimization is
    /// future work (D2).
    pub fn from_design(
        mut ambient: LayerSpec,
        substrate: LayerSpec,
        films: &[crate::structure::Layer],
        nk: &std::collections::HashMap<std::sync::Arc<str>, Vec<num_complex::Complex64>>,
        groups: &std::collections::HashMap<String, crate::structure::Group>,
        wavelengths: &[f64],
        background: &std::collections::HashSet<String>,
    ) -> Result<(Self, Vec<String>), String> {
        use std::collections::HashMap;
        // Two things are corrected-and-announced here rather than
        // refused, and this is the only place either can be: an absorbing
        // incident medium (R3.4), and a graded/gradient film the pipeline
        // cannot carry. Both follow the same rule -- never refused, never
        // silent.
        //
        // Profiled films go one of two ways (F1.1 extends the graded
        // rule to mixture gradients; F1.6/F1.7 reshape the second way):
        // - BACKGROUND (named in `background`): expanded WITH the profile
        //   and pinned (optimize/needle forced false on the whole carrier
        //   span). True physics, fixed: needle never hosts there, LM
        //   skips the rows, nk-keyed merge and flag-guarded cleanup preserve
        //   the span. No warning - explicit opt-in, nothing dropped.
        // - otherwise, when `optimize` is true: the carrier KEEPS its
        //   profile as ONE span, ONE LM parameter (the carrier's total
        //   thickness, F1.6). This used to be gated on the mode scaling
        //   exactly (Fixed / FixedSpan); F1.7 lands the profile-refresh
        //   machinery, so the RATE modes (whose profiles read absolute
        //   depth) keep their profiles too and are rebuilt at the
        //   construction points - the one non-additive part of F1.6's
        //   licence, completed: `optimize = true` never means
        //   "homogenize me" for any profiled film.
        // - otherwise HOMOGENIZED with a warning per film: the pipeline's
        //   operators assume uniform slabs. A graded film homogenizes to
        //   its base index (the flag flip below works because the base nk
        //   exists); a gradient film has NO base nk, so it must synthesise
        //   its own uniform row - EMA at f_mid over both endpoint spectra
        //   (R4: never reach for the flag flip here, there is nothing to
        //   flip to). The row is emitted by rewriting the film's own
        //   provider entry and dropping the gradient before `expand`, so
        //   the carrier name, flags and span bookkeeping are the uniform
        //   path's, and the emitted nk is bitwise the direct EMA call.
        let mut warnings = Vec::new();

        // THE LAYER-0 GATE (R3.4, 0.6.27). Every production `DesignStack`
        // is built here -- `assemble_stack`, `design_from_config` and the
        // PyO3 `DesignStack.from_design` all funnel through this function,
        // and every other `with_films` call in the crate is a test. So one
        // check here covers the synthesis surface completely.
        //
        // Once is enough, and this is why it is here rather than nearer the
        // solver: `ambient` is private with no mutator, so nothing
        // downstream -- needle insertion, merge, clamp, thickness steps --
        // can put the absorption back. Putting the check in
        // `solver_arrays()` instead would run it per merit evaluation,
        // thousands of times a run, to re-establish something that cannot
        // have changed.
        if let Some((fixed, msg)) =
            crate::smatrix::optics_core::sanitize_incident_index(&ambient.nk)
        {
            ambient.nk = fixed.into();
            warnings.push(msg);
        }

        // The nk table is built FIRST: the gradient homogenize below
        // resolves both endpoint spectra through it, rewrites the film's
        // own entry, and only then does the provider freeze it.
        let mut entries = HashMap::new();
        for (name, values) in nk {
            entries.insert(name.to_string(), values.clone());
        }
        let flat: Vec<crate::structure::Layer> = films
            .iter()
            .map(|l| {
                let mut h = l.clone();
                if h.inhomogen || h.gradient.is_some() {
                    if background.contains(h.material.as_str()) {
                        h.optimize = false;
                        h.needle = false;
                    } else if h.optimize {
                        // F1.6 (U2/B5) + F1.7: a profiled carrier with
                        // `optimize = true` keeps its profile AND its
                        // optimize flag: the span becomes ONE LM parameter
                        // (the carrier's total thickness), so the film is
                        // neither homogenized nor pinned. The fixed modes
                        // scale exactly; the rate modes ride
                        // `refresh_profiles` (their recipes are captured
                        // below). This is the one non-additive part of
                        // F1.6's licence, completed at F1.7: `optimize =
                        // true` used to mean "homogenize me" for profiled
                        // films. The needle flag stays as authored - a
                        // multi-row span is never a needle host
                        // (needle_host_refusal), so it cannot be
                        // meaningfully true here, but forcing it would
                        // be a second silent correction and this branch
                        // refuses to correct anything.
                    } else if h.inhomogen {
                        warnings.push(format!(
                            "DesignStack::from_design: film '{}' is graded (inhomogen) - the pipeline treats it as \
                             homogeneous (base index, profile dropped); needle/thickness steps assume uniform slabs. \
                             Pin it (optimize=False, needle=False) to keep the profile as background.",
                            h.material
                        ));
                        h.inhomogen = false;
                    } else {
                        // R4 (D4.3): the gradient homogenize. Both
                        // endpoint spectra must be registered (A5 keeps
                        // the path self-contained) - the refusals name
                        // the layer + the absent material, never a
                        // half-mixture fallback.
                        let grad = h.gradient.clone().expect("checked above");
                        if let Some(msg) = grad.expansion_error(h.inhomogen) {
                            return Err(format!("layer '{}': {msg}", h.material));
                        }
                        let wav_n = wavelengths.len();
                        let Some(nk_a) = entries.get(grad.material_a.as_str()) else {
                            return Err(crate::structure::gradient::endpoint_error(
                                &h.material,
                                &grad,
                                format!(
                                    "the provider does not carry '{}' (the film's own nk rides its name)",
                                    grad.material_a
                                ),
                            ));
                        };
                        let Some(nk_b) = entries.get(grad.material_b.as_str()) else {
                            return Err(crate::structure::gradient::endpoint_error(
                                &h.material,
                                &grad,
                                format!(
                                    "the provider does not carry '{}' (nk_b must ride the film dict)",
                                    grad.material_b
                                ),
                            ));
                        };
                        if nk_a.len() != wav_n || nk_b.len() != wav_n {
                            return Err(crate::structure::gradient::endpoint_error(
                                &h.material,
                                &grad,
                                format!(
                                    "endpoint spectrum length {} / {} != grid length {wav_n}",
                                    nk_a.len(),
                                    nk_b.len()
                                ),
                            ));
                        }
                        let f_mid = grad.f_mid(h.thickness);
                        let row = crate::structure::gradient::mix_row(
                            grad.ema,
                            nk_b,
                            nk_a,
                            f_mid,
                        );
                        entries.insert(h.material.as_str().to_string(), row);
                        h.gradient = None;
                        warnings.push(format!(
                            "DesignStack::from_design: film '{}' carries a gradient ('{}' -> '{}', {}); the pipeline treats \
                             it as homogeneous (EMA at f_mid = {:.6}); pin it (optimize=False, needle=False) to keep \
                             the profile as background.",
                            h.material,
                            grad.material_a,
                            grad.material_b,
                            grad.ema.name(),
                            f_mid
                        ));
                    }
                }
                Ok(h)
            })
            .collect::<Result<Vec<_>, String>>()?;
        let num_wavs = wavelengths.len();
        if ambient.nk.len() != num_wavs || substrate.nk.len() != num_wavs {
            return Err("DesignStack::from_design: boundary nk length != grid length.".to_string());
        }
        let provider = crate::structure::DictProvider::with_grid(
            entries
                .into_iter()
                .map(|(name, values)| (name, crate::structure::Entry::Array(values)))
                .collect(),
            wavelengths.to_vec(),
        )?;
        let seq: Vec<(crate::structure::Layer, bool)> =
            flat.iter().cloned().map(|l| (l, false)).collect();
        // F0.1 / §8.15: `expand`'s backwards-rescale branch (entry `k`
        // rescaling the previous span's already-emitted rows) is reachable
        // only through `inv = true` entries. The design path builds `inv =
        // false` throughout, so emission is purely local — the property
        // F1.7's per-span refresh rests on. Assert it beside the expansion
        // (debug-only; cannot move a fingerprint).
        debug_assert!(
            seq.iter().all(|(_, inv)| !*inv),
            "from_design builds inv = false throughout; a true entry would \
             make emission reach backwards into previous spans"
        );
        let (sa, spans, sa_bulk_nk, owner_of) = crate::structure::expand_with_recipe_inputs(
            &seq,
            &provider,
            wavelengths,
            groups,
            crate::structure::ExpandOptions::deterministic(),
        )?;
        // F1.7 (B2): the recipes are captured from the SAME emission that
        // built the spans - one recipe per profiled carrier, `None` for
        // every plain film. The design path is deterministic (B3), so the
        // prev operand is the previous entry's phase-1 nk, and the owner
        // resolved to Some(k) (self-owned) or None - the only shapes this
        // path can produce (asserted beside the expansion above).
        let groups = std::sync::Arc::new(groups.clone());
        let provider: std::sync::Arc<crate::structure::DictProvider> =
            std::sync::Arc::new(provider);
        let mut recipes: Vec<Option<SpanRecipe>> = spans.iter().map(|_| None).collect();
        for (si, span) in spans.iter().enumerate() {
            let k = span.logical;
            let carrier = &flat[k];
            if carrier.inhomogen || carrier.gradient.is_some() {
                // The per-span-refresh legality (B2): the design path's
                // owner resolves to the entry itself or None - a foreign
                // owner would make emission reach backwards into
                // previous spans and refresh ill-defined.
                debug_assert!(
                    owner_of[k].is_none_or(|o| o == k),
                    "the design path's owner must be the entry itself or None"
                );
                let emitted_total: f64 = sa.thicknesses[span.start..span.end].iter().sum();
                recipes[si] = Some(SpanRecipe {
                    layer: carrier.clone(),
                    groups: std::sync::Arc::clone(&groups),
                    bulk_nk: sa_bulk_nk[k].clone(),
                    self_owned: owner_of[k].is_some(),
                    prev_eff_nk: if k > 0 {
                        Some(sa_bulk_nk[k - 1].clone())
                    } else {
                        None
                    },
                    provider: std::sync::Arc::clone(&provider),
                    wavelengths: wavelengths.to_vec(),
                    opts: crate::structure::ExpandOptions::deterministic(),
                    emitted_total,
                });
            }
        }
        let mut rows = Vec::with_capacity(sa.n_rows());
        for span in spans.iter() {
            let carrier = &seq[span.logical].0;
            for r in span.start..span.end {
                let is_slice = span.slice && r == span.start;
                let nk_row: Vec<num_complex::Complex64> = sa.row(r).to_vec();
                rows.push(LayerSpec {
                    material: std::sync::Arc::from(carrier.material.as_str()),
                    nk: nk_row.into(),
                    d_nm: sa.thicknesses[r],
                    coherent: !sa.incoherent[r],
                    rough_type: sa.rough_types[r],
                    rough_val: sa.rough_vals[r],
                    optimize: carrier.optimize && !is_slice,
                    needle: carrier.needle && !is_slice,
                });
            }
        }
        Self::from_parts(ambient, substrate, rows, spans, recipes).and_then(|mut s| {
            // F1.7 (refresh point 1): a stack leaving any construction
            // door has fresh profiles. For the rate spans this re-runs
            // their emission at the totals just emitted - the gradient
            // branch's totals are exact, so this is a bitwise no-op
            // there; the legacy branch's uniform split sums to the
            // carrier only to float precision (N7), so the re-emission
            // reconciles the total, which is the honest extent. Fixed
            // modes are visited and left alone.
            s.refresh_profiles()?;
            Ok((s, warnings))
        })
    }

    // -- mutation primitives -------------------------------------------------

    /// Split the film at `film_idx` and insert a seed layer.
    ///
    /// Verbatim port of Python `_insert_needle`: replaces films[film_idx]
    /// with `[top_portion, seed, bottom_portion]` where both portions are
    /// clones of the original host (keeping its flags) re-thickened to
    /// `depth_into_layer_nm` and `d_total − depth_into_layer_nm`.
    ///
    /// Note: like the Python code this permits degenerate zero-thickness
    /// portions when `depth_into_layer_nm` is 0 or == host thickness;
    /// cleanup removes them later.
    pub fn insert_needle_seed(
        &mut self,
        film_idx: usize,
        depth_into_layer_nm: f64,
        seed: LayerSpec,
    ) -> Result<(), String> {
        let original = self.films.get(film_idx).ok_or_else(|| {
            format!(
                "film_idx {} out of range ({} films)",
                film_idx,
                self.films.len()
            )
        })?;
        // R3: the admissibility rule lives here too, not only in the scan
        // — this function is public and `build_scan_sites` is one caller
        // of it, not the only one.
        if let Some(msg) = self.needle_host_refusal(film_idx) {
            return Err(msg);
        }
        let d_total = original.d_nm;
        let d_bot = d_total - depth_into_layer_nm;

        let mut top = original.cloned();
        top.d_nm = depth_into_layer_nm;
        let mut bot = original.cloned();
        bot.d_nm = d_bot;

        // Replace in place without touching unrelated entries.
        self.films[film_idx] = top;
        self.films.insert(film_idx + 1, seed);
        self.films.insert(film_idx + 2, bot);

        // F0.1 span bookkeeping. The host row lies in a singleton-bulk
        // span (enforced above), so the split is local: rows up to and
        // including the re-thickened top stay in the host's span; the seed
        // becomes its own singleton; the bottom portion keeps the host's
        // logical identity in a new span; later spans shift by two. For
        // the 2-row (interface) host the top and bottom portions are
        // re-thickened slice rows, so the slice flag follows whichever
        // half kept the span's leading edge.
        let old_spans = std::mem::take(&mut self.spans);
        let old_recipes = std::mem::take(&mut self.recipes);
        let mut new_spans = Vec::with_capacity(old_spans.len() + 2);
        // F1.7: the recipe vector rides the rebuild. A split host's rows
        // changed, so all three resulting spans drop their recipes (the
        // affected spans keep their current rows and are simply never
        // refreshed); untouched spans - before or shifted past - keep
        // theirs (content-only, no row indices inside).
        let mut new_recipes: Vec<Option<SpanRecipe>> = Vec::with_capacity(old_recipes.len() + 2);
        for (sid, sp) in old_spans.iter().enumerate() {
            if sp.end <= film_idx {
                new_spans.push(*sp);
                new_recipes.push(old_recipes[sid].clone());
            } else if film_idx < sp.start {
                new_spans.push(Span {
                    start: sp.start + 2,
                    end: sp.end + 2,
                    logical: sp.logical,
                    slice: sp.slice,
                    bulk_start: sp.bulk_start + 2,
                });
                new_recipes.push(old_recipes[sid].clone());
            } else {
                // The host span [s, e) with film_idx in it.
                let (s, e) = (sp.start, sp.end);
                if e - s == 1 {
                    // One-row host: top keeps the span, seed and bot are new.
                    new_spans.push(Span {
                        start: s,
                        end: film_idx + 1,
                        logical: sp.logical,
                        slice: false,
                        bulk_start: s,
                    });
                    new_spans.push(Span {
                        start: film_idx + 1,
                        end: film_idx + 2,
                        logical: film_idx + 1,
                        slice: false,
                        bulk_start: film_idx + 1,
                    });
                    new_spans.push(Span {
                        start: film_idx + 2,
                        end: film_idx + 3,
                        logical: sp.logical,
                        slice: false,
                        bulk_start: film_idx + 2,
                    });
                } else if film_idx == s {
                    // Split at the slice row of a 2-row span: the top is
                    // the re-thickened slice (no bulk rows of its own);
                    // the bottom portion and the untouched bulk row form
                    // the slice-carrying remainder.
                    new_spans.push(Span {
                        start: s,
                        end: s + 1,
                        logical: sp.logical,
                        slice: true,
                        bulk_start: s + 1,
                    });
                    new_spans.push(Span {
                        start: s + 1,
                        end: s + 2,
                        logical: film_idx + 1,
                        slice: false,
                        bulk_start: s + 1,
                    });
                    new_spans.push(Span {
                        start: s + 2,
                        end: e + 2,
                        logical: sp.logical,
                        slice: true,
                        bulk_start: s + 3,
                    });
                } else {
                    // Split at the bulk row (the only other row of a
                    // singleton-bulk span): the slice keeps its row and
                    // the top stays beside it; the bottom portion is a
                    // plain row of the host's material.
                    new_spans.push(Span {
                        start: s,
                        end: e,
                        logical: sp.logical,
                        slice: sp.slice,
                        bulk_start: sp.bulk_start,
                    });
                    new_spans.push(Span {
                        start: e,
                        end: e + 1,
                        logical: e,
                        slice: false,
                        bulk_start: e,
                    });
                    new_spans.push(Span {
                        start: e + 1,
                        end: e + 2,
                        logical: sp.logical,
                        slice: false,
                        bulk_start: e + 1,
                    });
                }
                // The host's rows changed shape (split): every result
                // drops its recipe.
                new_recipes.push(None);
                new_recipes.push(None);
                new_recipes.push(None);
            }
        }
        self.spans = new_spans;
        self.recipes = new_recipes;
        self.assert_spans_partition();
        Ok(())
    }

    /// Merge consecutive film layers of identical material AND bit-identical
    /// nk. The nk key is load-bearing: interface slices, group-scaled rows
    /// and graded-profile sublayers share material names with their bulk
    /// but must never collapse into one slab (first-row nk would win).
    /// Uniform same-material rows carry identical vectors, so they still
    /// merge. Combined thickness is the sum; the *first* layer's other
    /// properties (optimize/needle/roughness/…) are preserved. Returns the
    /// number of merges performed.
    pub fn merge_adjacent(&mut self) -> usize {
        if self.films.is_empty() {
            return 0;
        }
        // Maximal merge runs over the CURRENT rows (same material AND
        // bit-identical nk), then rebuild films and spans in one pass so
        // the two cannot disagree (F0.1; N2 — the span vector is what
        // makes the nk-keyed merge safe to keep running over graded rows).
        let old_spans = std::mem::take(&mut self.spans);
        let films = std::mem::take(&mut self.films);
        let mut runs: Vec<(usize, usize)> = Vec::with_capacity(films.len());
        let mut i = 0usize;
        while i < films.len() {
            let mut j = i + 1;
            while j < films.len()
                && films[j].material == films[i].material
                && films[j].nk == films[i].nk
            {
                j += 1;
            }
            runs.push((i, j));
            i = j;
        }
        let merge_count: usize = runs.iter().map(|(s, e)| e - s - 1).sum();
        if merge_count == 0 {
            // Nothing merges: put the originals back untouched (identical
            // content to a full rebuild, without the churn).
            self.films = films;
            self.spans = old_spans;
            return 0;
        }

        let mut merged: Vec<LayerSpec> = Vec::with_capacity(runs.len());
        // Per output row: the span owning the run's leader row, and the
        // run's old-row extent [s, e).
        let mut owners: Vec<(usize, usize, usize)> = Vec::with_capacity(runs.len());
        for (s, e) in runs {
            let mut combined_d = films[s].d_nm;
            for r in s + 1..e {
                combined_d += films[r].d_nm;
            }
            let mut result = films[s].cloned();
            result.d_nm = combined_d;
            merged.push(result);
            owners.push((span_index_of_row(&old_spans, s), s, e));
        }

        // Rebuild the partition: consecutive output rows owned by the
        // same old span form one span. A span all of whose rows merged
        // into an earlier run disappears (the merged row joined the FIRST
        // span touched — "first layer's properties win", matching the
        // films rule). The slice flag survives only when the old span's
        // slice row IS its surviving group's first row; if the slice row
        // was absorbed into an earlier run, the flag would point at a
        // row that no longer is one.
        let old_recipes = std::mem::take(&mut self.recipes);
        let mut new_spans: Vec<Span> = Vec::new();
        let mut new_recipes: Vec<Option<SpanRecipe>> = Vec::new();
        let mut p = 0usize;
        let mut g = 0usize;
        while g < owners.len() {
            let (sid, leader, _run_end) = owners[g];
            let mut h = g + 1;
            while h < owners.len() && owners[h].0 == sid {
                h += 1;
            }
            let old = &old_spans[sid];
            let slice_kept = old.slice && leader == old.start;
            new_spans.push(Span {
                start: p,
                end: p + (h - g),
                logical: old.logical,
                slice: slice_kept,
                bulk_start: p + usize::from(slice_kept),
            });
            // F1.7 recipe rule: the span's rows survived WHOLE (every old
            // row is still here, none absorbed across the boundary, the
            // leading edge intact) - an INTRA-span merge (the RateCapped
            // saturated tail legitimately merges into itself, F1.2) only
            // re-cut the rows, and refresh rebuilds them from the
            // recipe's carrier at the new total. Anything else - a
            // foreign row absorbed into this span, this span's leading
            // edge eaten by the previous one - means re-emission would
            // NOT reproduce what stands here: the recipe drops.
            let old_lo = leader;
            let old_hi = owners[h - 1].2;
            let rows_whole = old_lo == old.start && old_hi == old.end;
            new_recipes.push(if rows_whole {
                old_recipes[sid].clone()
            } else {
                None
            });
            p += h - g;
            g = h;
        }

        self.films = merged;
        self.spans = new_spans;
        self.recipes = new_recipes;
        self.assert_spans_partition();
        merge_count
    }

    /// Remove film at `film_idx`. Returns the removed layer.
    ///
    /// F0.1: the containing span shrinks with its rows; a span reduced to
    /// nothing is dropped; later spans renumber down by one. Removing a
    /// 2-row span's bulk leaves the slice row alone in its span — a
    /// derived row whose carrier is gone; that state is recorded honestly
    /// (slice flag kept, not singleton-bulk) and F0.2's clamp rules are
    /// what give it a policy.
    pub fn remove_film(&mut self, film_idx: usize) -> Result<LayerSpec, String> {
        if film_idx >= self.films.len() {
            return Err(format!(
                "film_idx {} out of range ({} films)",
                film_idx,
                self.films.len()
            ));
        }
        let old_spans = std::mem::take(&mut self.spans);
        let old_recipes = std::mem::take(&mut self.recipes);
        let mut new_spans = Vec::with_capacity(old_spans.len());
        let mut new_recipes: Vec<Option<SpanRecipe>> = Vec::with_capacity(old_recipes.len());
        for (sid, sp) in old_spans.into_iter().enumerate() {
            if sp.end <= film_idx {
                new_spans.push(sp);
                new_recipes.push(old_recipes[sid].clone());
            } else if film_idx < sp.start {
                new_spans.push(Span {
                    start: sp.start - 1,
                    end: sp.end - 1,
                    logical: sp.logical,
                    slice: sp.slice,
                    bulk_start: sp.bulk_start - 1,
                });
                new_recipes.push(old_recipes[sid].clone());
            } else if sp.end - sp.start > 1 {
                // Containing span with rows surviving the removal: the
                // row set changed, so the recipe drops (the surviving
                // rows keep their current values and are never
                // refreshed).
                let slice_survives = sp.slice && film_idx != sp.start;
                new_spans.push(Span {
                    start: sp.start,
                    end: sp.end - 1,
                    logical: sp.logical,
                    slice: slice_survives,
                    bulk_start: sp.start + usize::from(slice_survives),
                });
                new_recipes.push(None);
            }
            // else: the removed row was the span's only row — dropped.
        }
        self.spans = new_spans;
        self.recipes = new_recipes;
        let removed = self.films.remove(film_idx);
        self.assert_spans_partition();
        Ok(removed)
    }

    /// Enforce [min_nm, max_nm] on every film layer.
    ///
    /// F0.2: the floor and the cap are span quantities. The comparison
    /// reads `D`, the span's slice-inclusive total - expansion carves the
    /// slice out of the carrier, so the slice IS part of the authored
    /// thickness (B1): excluded from scaling, never from measuring. A
    /// slice row is never a floor or a cap candidate on its own, in every
    /// branch; it leaves the stack only when its carrier span does.
    ///
    /// - `D < min_nm` -> the whole span is removed, all rows in one
    ///   operation, and the report names it (a silent correct deletion
    ///   and a silent wrong one look identical otherwise).
    /// - `D > max_nm` on a multi-row span -> refused, not rescaled: a
    ///   pinned profile above the manufacturing ceiling is an authoring
    ///   error, and row-by-row capping would warp it. The pipeline
    ///   pre-checks at `NeedlePipeline::new`, so a run never hits this
    ///   mid-flight except across a merge that grew a span.
    /// - A one-row span is row and span in one (`D == d`): removal and
    ///   capping behave exactly as before - a no-span run is
    ///   bit-identical, full stop.
    ///
    /// Returns the [`ClampReport`] - the old `(n_removed, n_capped)` grew
    /// into it. Verbatim port lineage:
    /// `ClampedNeedleSynthesizer.clamp_all_layers`.
    /// F0.3: `clamp_up` selects the floor's behaviour - `false` removes a
    /// sub-minimum span (today, and every in-run pass under
    /// `Remove`/`ClampUpFinal`), `true` clamps a surviving SINGLETON-BULK
    /// span up to `clamp_min_nm` instead (the final pass under
    /// `ClampUpFinal`, every pass under `ClampUpAlways`). Singleton-bulk is
    /// N1's "one physical layer", not a row count (C2): a plain
    /// interface-carrying film is two rows and still one layer - the slice
    /// stays (B1: it belongs to the authored thickness) and the bulk row
    /// takes the remainder, so the film is clamped instead of deleted. A
    /// multi-row span is removed whole under every policy; the exceptions
    /// are F1.6's branches above this one (a scalable span scales to the
    /// floor, fractions preserved). Clamp-ups are neither removals nor
    /// caps, so they do not appear in the report.
    pub fn clamp_all(
        &mut self,
        min_nm: f64,
        max_nm: f64,
        clamp_up: bool,
    ) -> Result<ClampReport, String> {
        // Legacy/Python surface: the bool form. `true` means the
        // strongest clamp-up posture (rows clamped up as soon as they
        // are seen, spans scaled up); `false` is plain `Remove`.
        // ClampUpFinal reaches the sweep through `clamp_all_policy`.
        let policy = if clamp_up {
            ThinLayerPolicy::ClampUpAlways
        } else {
            ThinLayerPolicy::Remove
        };
        self.clamp_all_policy(min_nm, max_nm, policy, clamp_up)
    }

    /// The policy-aware clamp (F0.2/F0.3/F1.6). `final_pass` reproduces
    /// F0.3's pass structure exactly: DURING the run only ClampUpAlways
    /// clamps anything up; the FINAL pass substitutes clamping for
    /// removal under both clamp-up policies. The span rules ride the
    /// same structure:
    ///
    /// - the ceiling REFUSAL narrows to PINNED multi-row spans (F0.2); a
    ///   SCALABLE span (every bulk row optimize-flagged — the shape
    ///   `from_design` gives a profiled carrier with `optimize = true`,
    ///   F1.6) is CAPPED to the ceiling as a bound, fractions preserved;
    /// - an under-thickness SCALABLE span under an active clamp-up is
    ///   SCALED to the floor, fractions preserved (F0.3's deferral
    ///   lands here: no longer "removed whole under every policy");
    ///   otherwise it is removed whole as before;
    /// - every other rule is unchanged, and a stack with no scalable
    ///   span behaves exactly as before.
    pub fn clamp_all_policy(
        &mut self,
        min_nm: f64,
        max_nm: f64,
        policy: ThinLayerPolicy,
        final_pass: bool,
    ) -> Result<ClampReport, String> {
        debug_assert!(min_nm >= 0.0 && max_nm > min_nm);
        let clamp_up_rows = if final_pass {
            policy != ThinLayerPolicy::Remove
        } else {
            policy == ThinLayerPolicy::ClampUpAlways
        };
        // Refusals are checked BEFORE any mutation, so the stack is never
        // left half-clamped behind an error. Scalable spans are exempt:
        // their total is an LM bound (ub), not an authoring error, and
        // the sweep below caps them.
        for sp in &self.spans {
            if sp.end - sp.start <= 1 || self.span_is_scalable(sp) {
                continue;
            }
            let d: f64 = self.films[sp.start..sp.end].iter().map(|l| l.d_nm).sum();
            if d > max_nm {
                return Err(format!(
                    "clamp_all: span '{}' is {:.1} nm thick, above the {:.1} nm ceiling - \
                     refusing rather than rescaling a profile",
                    self.films[sp.start].material, d, max_nm
                ));
            }
        }

        let old = std::mem::take(&mut self.films);
        let old_spans = std::mem::take(&mut self.spans);
        let old_recipes = std::mem::take(&mut self.recipes);
        let mut surviving: Vec<LayerSpec> = Vec::with_capacity(old.len());
        let mut new_spans: Vec<Span> = Vec::with_capacity(old_spans.len());
        let mut new_recipes: Vec<Option<SpanRecipe>> = Vec::with_capacity(old_recipes.len());
        let mut report = ClampReport::default();
        let mut p = 0usize;

        for (sid, sp) in old_spans.iter().enumerate() {
            let d: f64 = old[sp.start..sp.end].iter().map(|l| l.d_nm).sum();
            // The span's scalability is read from the OLD rows (the
            // predicate only consults optimize flags, which the rebuild
            // preserves per output row).
            let scalable = span_is_scalable_rows(sp, &old);
            if d < min_nm {
                // F0.3 + F1.6: a scalable multi-row span under an active
                // clamp-up posture is SCALED to the floor, fractions
                // preserved (the profile survives; this is the deferral
                // landing). Under `Remove`, and for every non-scalable
                // span, the F0.3 rule stands: a singleton-bulk span clamps
                // up, multi-row spans are removed whole. The clamp-up
                // predicate is N1's "one physical layer", not a row
                // count (C2): a plain interface-carrying film is two rows
                // and still one layer - the slice stays (B1: it belongs
                // to the authored thickness) and the bulk row takes the
                // rest, so the film is clamped instead of deleted.
                if scalable && sp.end - sp.start > 1 && clamp_up_rows {
                    let f = min_nm / d;
                    let mut rows = old[sp.start..sp.end].to_vec();
                    for r in rows.iter_mut() {
                        r.d_nm *= f;
                    }
                    let e = p + rows.len();
                    surviving.extend(rows);
                    new_spans.push(Span {
                        start: p,
                        end: e,
                        logical: sp.logical,
                        slice: sp.slice,
                        bulk_start: p + usize::from(sp.slice),
                    });
                    // The span survived with all its rows (scaled) - the
                    // recipe rides 1:1 and refresh rebuilds the profile
                    // at the new floor.
                    new_recipes.push(old_recipes[sid].clone());
                    p = e;
                    continue;
                }
                if clamp_up_rows && sp.is_singleton_bulk() {
                    let mut thin = old[sp.start..sp.end].to_vec();
                    let slice_rows = usize::from(sp.slice);
                    let slice_total: f64 = thin[..slice_rows].iter().map(|l| l.d_nm).sum();
                    // The clamp sets the SLICE-INCLUSIVE total to the
                    // floor (F0.2's rule reads the span's total); the
                    // bulk row takes the remainder. A one-row span is
                    // the slice_rows == 0 case of this, bitwise the old
                    // `thin[0].d_nm = min_nm`.
                    thin[slice_rows].d_nm = min_nm - slice_total;
                    let e = p + thin.len();
                    surviving.extend(thin);
                    new_spans.push(Span {
                        start: p,
                        end: e,
                        logical: sp.logical,
                        slice: sp.slice,
                        bulk_start: p + slice_rows,
                    });
                    new_recipes.push(old_recipes[sid].clone());
                    p = e;
                    continue;
                }
                report
                    .spans_removed
                    .push(format!("{} ({:.1} nm)", old[sp.start].material, d));
                report.rows_removed += sp.end - sp.start;
                continue;
            }
            let mut rows = old[sp.start..sp.end].to_vec();
            if rows.len() == 1 && rows[0].d_nm > max_nm {
                // One-row span: D == d, so this is today's row cap.
                rows[0].d_nm = max_nm;
                report.spans_capped += 1;
            } else if rows.len() > 1 && scalable && d > max_nm {
                // F1.6: a scalable span above the ceiling is capped as a
                // bound (fractions preserved), not refused - the refusal
                // above fired only for pinned spans.
                let f = max_nm / d;
                for r in rows.iter_mut() {
                    r.d_nm *= f;
                }
                report.spans_capped += 1;
            }
            let e = p + rows.len();
            surviving.extend(rows);
            new_spans.push(Span {
                start: p,
                end: e,
                logical: sp.logical,
                slice: sp.slice,
                bulk_start: p + usize::from(sp.slice),
            });
            new_recipes.push(old_recipes[sid].clone());
            p = e;
        }

        self.films = surviving;
        self.spans = new_spans;
        self.recipes = new_recipes;
        self.assert_spans_partition();
        Ok(report)
    }

    /// Test-only: force a row's optimize flag. The F0.2 span-exemption
    /// twin needs an optimize=true graded span, which `from_design` never
    /// produces before F1.6 (non-background graded films homogenize).
    #[cfg(test)]
    pub(crate) fn set_row_optimize_for_test(&mut self, row: usize, optimize: bool) {
        self.films[row].optimize = optimize;
    }

    /// Set the thickness of film `film_idx` (used by the thickness optimizer).
    pub fn set_thickness(&mut self, film_idx: usize, d_nm: f64) -> Result<(), String> {
        let f = self
            .films
            .get_mut(film_idx)
            .ok_or_else(|| format!("film_idx {} out of range", film_idx))?;
        f.d_nm = d_nm;
        Ok(())
    }

    // -- profile refresh (F1.7, U3/U4) ---------------------------------------

    /// Rebuild every rate-mode span's profile from its current total
    /// thickness.
    ///
    /// A rate-type profile is a function of absolute thickness -
    /// `f(z) = f_start + rate * z / ref_thickness` (gradient), `delta =
    /// min(rate * D / ref, cap)` (legacy) - so scaling the span makes its
    /// stored nk stale, and the span is re-emitted here. Fixed-mode spans
    /// are visited and left alone: their profile was never a function of
    /// the total (U3, verified in-tree).
    ///
    /// **Construction points only (U4, and it governs this item):**
    /// `from_design`, the top of each macro cycle (before the budget
    /// check and needle scan), and after the final clamp pass - and
    /// nowhere else. Never inside an LM round, never inside the
    /// Jacobian, never inside a needle scan or a cleanup trial: the row
    /// count of a span is frozen for the duration of one LM solve (the
    /// residual must not change length mid-solve, and a live profile
    /// would put step discontinuities at every `ceil()` boundary). The
    /// refresh-count twin pins the count mechanically; a fourth call
    /// site appearing later makes it fail loudly.
    ///
    /// Row count may change HERE and only here, so
    /// `assert_spans_partition` runs immediately after, and the F1.6
    /// `fractions` are re-captured on the next `optimize_thicknesses`
    /// call rather than carried across a refresh. The per-span legality
    /// (re-emitting span k cannot disturb span k-1) rests on the design
    /// path's `inv = false` emission - the debug assert below records
    /// it.
    pub(crate) fn refresh_profiles(&mut self) -> Result<(), String> {
        if REFRESH_COUNTING.with(std::cell::Cell::get) {
            REFRESH_CALLS.with(|c| c.set(c.get() + 1));
        }
        for si in 0..self.spans.len() {
            let Some(recipe) = self.recipes[si].as_ref() else {
                continue;
            };
            // B3: refresh never re-rolls the error ensemble. Unreachable
            // through `from_design` today (it always expands
            // deterministic); if a tolerancing path ever grows spans,
            // this fails loudly at the door instead of quietly changing
            // a per-run ensemble into a per-cycle one.
            if recipe.opts.apply_errors {
                return Err(format!(
                    "refresh_profiles: span '{}' was built with apply_errors = true - \
                     refreshing would re-roll the error ensemble; span recipes are \
                     construction-only (deterministic) today",
                    self.films[self.spans[si].start].material
                ));
            }
            if !recipe.is_rate_mode() {
                continue; // fixed modes: bitwise no-op by design
            }
            let sp = self.spans[si];
            // The span's current slice-inclusive extent (B1): the
            // interface carve draws from the carrier's thickness, so the
            // slice IS part of the extent being refreshed.
            let total: f64 = self.films[sp.start..sp.end].iter().map(|l| l.d_nm).sum();
            // A span whose rows still sum (bitwise) to the measured total
            // of the last emission has not moved: its profile already
            // describes this extent exactly, so refresh is a bitwise
            // no-op. This is what makes the from_design call site a true
            // no-op for freshly built spans instead of an ulp-chasing
            // re-emission, and it is one-step convergent after a move
            // (the marker is re-measured from the new rows).
            if total == recipe.emitted_total {
                continue;
            }
            // The sublayer-count and delta rules read the CARRIER's
            // thickness (at construction it equals the emitted extent:
            // identity groups on this path); refresh feeds the current
            // extent through both channels.
            let mut emitted_layer = recipe.layer.clone();
            emitted_layer.thickness = total;
            let block = crate::structure::emit_standalone(
                &emitted_layer,
                &recipe.groups,
                &recipe.bulk_nk,
                total,
                recipe.self_owned,
                recipe.prev_eff_nk.clone(),
                recipe.opts,
                &*recipe.provider,
                &recipe.wavelengths,
            )?;
            let carrier = &recipe.layer;
            let material: Arc<str> = Arc::from(carrier.material.as_str());
            let n_old = sp.end - sp.start;
            let n_new = block.thicknesses.len();
            let mut new_rows = Vec::with_capacity(n_new);
            for (r, d) in block.thicknesses.iter().enumerate() {
                let is_slice = block.span.slice && r == 0;
                new_rows.push(LayerSpec {
                    material: material.clone(),
                    nk: block.nk[r].clone().into(),
                    d_nm: *d,
                    coherent: block.coherent[r],
                    rough_type: block.rough_types[r],
                    rough_val: block.rough_vals[r],
                    optimize: carrier.optimize && !is_slice,
                    needle: carrier.needle && !is_slice,
                });
            }
            let delta = n_new as isize - n_old as isize;
            self.films.splice(sp.start..sp.end, new_rows);
            // This span's record, re-based onto the stack (its logical
            // identity is stable across mutators; the emission's local
            // start is 0).
            self.spans[si] = Span {
                start: sp.start,
                end: sp.start + n_new,
                logical: sp.logical,
                slice: block.span.slice,
                bulk_start: sp.start + usize::from(block.span.slice),
            };
            // Later spans shift with the row count change; their recipes
            // are content-only (no row indices) and follow untouched.
            if delta != 0 {
                for later in &mut self.spans[si + 1..] {
                    later.start = (later.start as isize + delta) as usize;
                    later.end = (later.end as isize + delta) as usize;
                    later.bulk_start = (later.bulk_start as isize + delta) as usize;
                }
            }
            // Re-measure the marker from the new rows (one-step
            // convergence: the next refresh sees an unmoved span).
            let new_total: f64 = self.films[sp.start..sp.start + n_new]
                .iter()
                .map(|l| l.d_nm)
                .sum();
            if let Some(r) = self.recipes[si].as_mut() {
                r.emitted_total = new_total;
            }
        }
        self.assert_spans_partition();
        Ok(())
    }

    // -- solver-array materialization ----------------------------------------

    /// Materialize flat solver arrays in the layouts consumed by
    /// `core_engine` / `needle_engine`.
    pub fn solver_arrays(&self) -> SolverArrays {
        let n_layers = self.total_layer_count();
        let nw = self.num_wavs;

        let mut n_stack_cache = vec![0.0f64; nw * n_layers * 2];
        let mut thicknesses = vec![0.0f64; n_layers];
        let mut incoherent_flags = vec![0i32; n_layers];
        let mut rough_types = vec![0i32; n_layers];
        let mut rough_vals = vec![0.0f64; n_layers];

        for (slot, layer) in std::iter::once(&self.ambient)
            .chain(self.films.iter())
            .chain(std::iter::once(&self.substrate))
            .enumerate()
        {
            thicknesses[slot] = layer.d_nm;
            incoherent_flags[slot] = i32::from(!layer.coherent);
            rough_types[slot] = layer.rough_type;
            rough_vals[slot] = layer.rough_val;
            for w in 0..nw {
                let base = w * n_layers * 2 + slot * 2;
                n_stack_cache[base] = layer.nk[w].re;
                n_stack_cache[base + 1] = layer.nk[w].im;
            }
        }

        SolverArrays {
            n_stack_cache,
            thicknesses,
            incoherent_flags,
            rough_types,
            rough_vals,
            n_layers: n_layers as i32,
        }
    }
}

// ---------------------------------------------------------------------------
// Span helpers
// ---------------------------------------------------------------------------

/// Index of the span containing `row` in a sorted, contiguous partition.
/// Row `row` must be covered by `spans` (callers validate); the binary
/// search exploits the ordering `assert_spans_partition` guarantees.
pub(crate) fn span_index_of_row(spans: &[Span], row: usize) -> usize {
    spans.partition_point(|s| s.end <= row)
}

/// F1.6: the scalability predicate over an explicit row list — the form
/// the clamp needs (it judges the OLD rows while rebuilding the stack).
pub(crate) fn span_is_scalable_rows(sp: &Span, films: &[LayerSpec]) -> bool {
    let bulk = sp.end - sp.bulk_start;
    bulk >= 2 && (sp.bulk_start..sp.end).all(|r| films[r].optimize)
}

// ---------------------------------------------------------------------------
// SpanRecipe + profile refresh (F1.7)
// ---------------------------------------------------------------------------

/// What one span's re-emission needs (B2) - everything one iteration of
/// `expand`'s emission loop reads, captured at the door that built the
/// span. `None` for every plain film, so a stack with no profiled layer
/// pays one `Vec` of `None` and nothing else.
///
/// Per-span refresh is legal because the design path builds `inv =
/// false` throughout: `expand`'s backwards-rescale branch (entry `k`
/// rescaling the PREVIOUS span's already-emitted rows) is unreachable
/// here, so re-emitting span k cannot disturb span k-1. That is a
/// property of the CALLER, not of the algorithm - `refresh_profiles`
/// records it as a `debug_assert`.
///
/// The recipe rides `DesignStack::recipes`, aligned index-for-index with
/// `spans` (never a map keyed by `logical` - a second index that
/// renumbers on `remove_film` is exactly the failure mode F0.1 exists
/// to prevent). Any mutator that alters a span's ROW SET drops its
/// recipe to `None` (the affected span keeps its current rows and is
/// simply never refreshed); spans that only shift keep theirs.
#[derive(Clone)]
pub(crate) struct SpanRecipe {
    /// The carrier as the door rewrote it (flags, interface, roughness,
    /// gradient spec, inh mode). `thickness` is the CONSTRUCTION value,
    /// kept for the record; refresh feeds the span's current extent
    /// through `emit_standalone` directly.
    pub(crate) layer: crate::structure::Layer,
    /// The groups map the emission resolves policy through (the
    /// carrier's own group and, for gradients, each endpoint's).
    pub(crate) groups: std::sync::Arc<HashMap<String, crate::structure::Group>>,
    /// The entry's phase-1 spectrum (provider-resolved, group-scaled).
    pub(crate) bulk_nk: Vec<num_complex::Complex64>,
    /// Whether the span carries its own interface flag (the design
    /// path's owner resolves to `Some(k)` - self-owned - or `None`).
    pub(crate) self_owned: bool,
    /// The PREVIOUS entry's terminal spectrum - the slice-mix "prev"
    /// operand. `None` for the stack's first span.
    pub(crate) prev_eff_nk: Option<Vec<num_complex::Complex64>>,
    /// The run's provider (gradient endpoints resolve through it).
    /// Concrete: the design path always builds a `DictProvider`, and a
    /// trait object here would drag `!Send`/`!Sync` implementors
    /// (`SpecProvider`'s RefCell cache, the PyO3 shelf) into a stack
    /// type the LM bounds to `Send + Sync`.
    pub(crate) provider: std::sync::Arc<crate::structure::DictProvider>,
    /// The run's wavelength grid (the gradient sublayer rule reads it).
    pub(crate) wavelengths: Vec<f64>,
    /// The options the emission ran under (B3: refresh refuses
    /// `apply_errors: true` rather than re-rolling the ensemble).
    pub(crate) opts: crate::structure::ExpandOptions,
    /// The extent the rows currently describe, MEASURED as the sum of the
    /// span's rows at the last emission. A span whose current sum still
    /// equals this (bitwise) has not moved since its profile was built -
    /// refresh skips it, which is what makes refresh at construction a
    /// bitwise no-op and the rate span's staleness detection exact
    /// rather than heuristic. Updated after every re-emission to the new
    /// rows' measured sum, so the next refresh skips it too (one-step
    /// convergence, no ulp chasing).
    pub(crate) emitted_total: f64,
}

impl std::fmt::Debug for SpanRecipe {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SpanRecipe")
            .field("layer", &self.layer)
            .field("bulk_nk.len", &self.bulk_nk.len())
            .field("self_owned", &self.self_owned)
            .field("prev_eff_nk", &self.prev_eff_nk.as_ref().map(|v| v.len()))
            .field("wavelengths", &self.wavelengths)
            .field("opts", &self.opts)
            .finish_non_exhaustive()
    }
}

impl SpanRecipe {
    /// Whether the recipe's profile is a function of the span's absolute
    /// thickness (U3) - the spans `refresh_profiles` rebuilds. The rate
    /// modes: `InhMode::RateCapped` (legacy single-material, delta grows
    /// with thickness) and `GradientMode::RateCapped` (f grows with
    /// depth). The fixed modes are scale-free (U3, verified) and are
    /// left alone.
    fn is_rate_mode(&self) -> bool {
        use crate::structure::gradient::{GradientMode, InhMode};
        if let Some(grad) = &self.layer.gradient {
            matches!(grad.mode, GradientMode::RateCapped { .. })
        } else if self.layer.inhomogen {
            matches!(self.layer.inh_mode, InhMode::RateCapped { .. })
        } else {
            // A rate mode on a non-profiled layer is inert (F1.3's
            // warning covers the authoring side); nothing to refresh.
            false
        }
    }
}

thread_local! {
    /// F1.7's refresh-count instrumentation (the U4 gate's counter).
    /// Counts `refresh_profiles` CALLS while a test holds the counting
    /// flag on. Thread-local, so parallel tests cannot pollute each
    /// other's counts, and off by default so no production or test path
    /// pays for it.
    static REFRESH_CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static REFRESH_COUNTING: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Turn the refresh counter on (the refresh-count twin's harness).
#[cfg(test)]
pub(crate) fn refresh_counting_on() {
    REFRESH_COUNTING.with(|c| c.set(true));
}

/// Turn the refresh counter off.
#[cfg(test)]
pub(crate) fn refresh_counting_off() {
    REFRESH_COUNTING.with(|c| c.set(false));
}

/// Zero the refresh counter.
#[cfg(test)]
pub(crate) fn refresh_count_reset() {
    REFRESH_CALLS.with(|c| c.set(0));
}

/// The refresh counter's value since the last reset.
#[cfg(test)]
pub(crate) fn refresh_count() -> usize {
    REFRESH_CALLS.with(|c| c.get())
}

// ---------------------------------------------------------------------------
// SolverArrays
// ---------------------------------------------------------------------------

/// Flat solver input arrays (mirrors loom_structure.SolverArrays).
#[derive(Clone, Debug)]
pub struct SolverArrays {
    /// Wav-major, re/im interleaved: base = w * n_layers * 2.
    pub n_stack_cache: Vec<f64>,
    /// Per-layer thickness (nm), including ambient + substrate.
    pub thicknesses: Vec<f64>,
    /// 1 = incoherent spacer, 0 = coherent.
    pub incoherent_flags: Vec<i32>,
    pub rough_types: Vec<i32>,
    pub rough_vals: Vec<f64>,
    pub n_layers: i32,
}

// ---------------------------------------------------------------------------
// Tests — mirror Python behaviors verbatim
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smatrix::synthesis::evaluator::{apply_params, build_params};
    use std::collections::{HashMap, HashSet};

    const NW: usize = 8; // simulation wavelengths

    fn air(d: f64) -> LayerSpec {
        let mut l = LayerSpec::constant("air", 1.0, 0.0, d, NW);
        l.optimize = false;
        l.needle = false;
        l
    }

    fn sub(d: f64) -> LayerSpec {
        let mut l = LayerSpec::constant("sub", 1.52, 0.0, d, NW);
        l.optimize = false;
        l.needle = false;
        l
    }

    fn stack(films: Vec<LayerSpec>) -> DesignStack {
        DesignStack::with_films(air(0.0), sub(0.0), films).unwrap()
    }

    /// F1.1 (R4): a non-background gradient film homogenizes to ONE row
    /// that is bitwise the direct EMA call at f_mid over both endpoint
    /// spectra; the warning names the mixture and f_mid. A background
    /// gradient film expands WITH the full profile, pinned.
    #[test]
    fn gradient_homogenize_is_bitwise_the_ema_and_background_expands() {
        use crate::materials::MixRule;
        use crate::structure::gradient::{GradientSpec, mix_row};
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.35, 0.01); NW]);
        nk.insert(Arc::from("L"), vec![Complex64::new(1.46, 0.0); NW]);
        let mut film = crate::structure::Layer::film(100.0, "H");
        film.gradient = Some(GradientSpec::fixed_span("H", "L", 0.0, 1.0));
        // F1.6: optimize=true on a FixedSpan gradient now KEEPS the
        // profile (the scalable-span branch), so the homogenize path is
        // exercised with the stable posture that always homogenized:
        // optimize=false + needle=true (not background - background
        // requires needle=false too).
        film.optimize = false;
        film.needle = true;

        // Non-background: one row, bitwise the EMA at f_mid, warning
        // names the mixture and f_mid.
        let (stack, warns) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            std::slice::from_ref(&film),
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        assert_eq!(stack.films().len(), 1);
        assert!((stack.films()[0].d_nm - 100.0).abs() < 1e-12);
        let want = mix_row(
            MixRule::Bruggeman {
                max_iter: 100,
                tol: 1e-9,
            },
            &nk[&Arc::from("L")],
            &nk[&Arc::from("H")],
            0.5,
        );
        assert_eq!(&stack.films()[0].nk[..], &want[..]);
        assert_eq!(warns.len(), 1, "{warns:?}");
        assert!(
            warns[0].contains("'H'") && warns[0].contains("'L'"),
            "{}",
            warns[0]
        );
        assert!(warns[0].contains("f_mid = 0.5"), "{}", warns[0]);
        assert!(warns[0].contains("Bruggeman"), "{}", warns[0]);

        // Background: the FULL profile expands, pinned, silent.
        let mut bg = HashSet::new();
        bg.insert("H".to_string());
        let (bstack, warns) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            std::slice::from_ref(&film),
            &nk,
            &HashMap::new(),
            &wl,
            &bg,
        )
        .unwrap();
        assert!(warns.is_empty());
        assert!(bstack.films().len() > 1, "profile expanded, not flattened");
        let d: f64 = bstack.films().iter().map(|f| f.d_nm).sum();
        assert!((d - 100.0).abs() < 1e-9);
        assert!(bstack.films().iter().all(|f| !f.optimize && !f.needle));
        assert!(!bstack.span_of_row(0).is_singleton_bulk());

        // Provider door: a missing inclusion spectrum refuses, naming
        // the absent material and the layer - never a half-mixture.
        let mut partial = HashMap::new();
        partial.insert(Arc::from("H"), vec![Complex64::new(2.35, 0.01); NW]);
        let err = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            std::slice::from_ref(&film),
            &partial,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap_err();
        assert!(err.contains("layer 'H'"), "{err}");
        assert!(err.contains("'L'"), "{err}");
    }

    /// D1: design films expand through the structure model (interface
    /// slices appear as derived rows; graded homogenizes, background pins).
    #[test]
    fn from_design_expands_interfaces_and_refuses_grading() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.35, 0.0); NW]);
        let mut films = vec![
            crate::structure::Layer::film(20.0, "H"),
            crate::structure::Layer::film(50.0, "H"),
        ];
        films[1].interface = true;
        films[1].interface_thickness = 5.0;
        let (stack, warns) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &films,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        assert!(warns.is_empty());
        // Plain film + slice + bulk (ambient/substrate flank outside films).
        assert_eq!(stack.films().len(), 3);
        assert!((stack.films()[1].d_nm - 5.0).abs() < 1e-12);
        assert!((stack.films()[2].d_nm - 45.0).abs() < 1e-12);
        // Slice rows are derived: never free, never hosts.
        assert!(!stack.films()[1].optimize && !stack.films()[1].needle);
        assert!(stack.films()[2].optimize && stack.films()[2].needle);
        assert_eq!(stack.films()[1].material.as_ref(), "H");
        // Graded input homogenizes with a warning (base index, one row).
        let mut graded = vec![crate::structure::Layer::film(50.0, "H")];
        graded[0].inhomogen = true;
        graded[0].inh_delta = 0.2;
        // F1.6: the stable homogenize posture (optimize=false +
        // needle=true - not background, which requires needle=false).
        graded[0].optimize = false;
        graded[0].needle = true;
        let (gstack, warns) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &graded,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        assert_eq!(warns.len(), 1);
        assert!(warns[0].contains("homogeneous"));
        assert_eq!(gstack.films().len(), 1);
        assert!((gstack.films()[0].d_nm - 50.0).abs() < 1e-12);
        // Background graded: profile expands (11 sublayers), all pinned.
        let bg: std::collections::HashSet<String> = ["H".to_string()].into_iter().collect();
        let (bstack, bwarns) =
            DesignStack::from_design(air(0.0), sub(0.0), &graded, &nk, &HashMap::new(), &wl, &bg)
                .unwrap();
        assert!(bwarns.is_empty());
        assert_eq!(bstack.films().len(), 11);
        assert!(bstack.films().iter().all(|f| !f.optimize && !f.needle));
        // Pinned span survives merge (nk-keyed) and was never a host.
        let mut mstack = bstack;
        assert_eq!(mstack.merge_adjacent(), 0);
        assert!((mstack.total_thickness_nm() - 50.0).abs() < 1e-12);
    }

    /// R3.4 (0.6.27): the layer-0 gate lives here, because this is the only
    /// place a production `DesignStack` is built.
    #[test]
    fn from_design_drops_ambient_absorption_once_and_says_so() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.35, 0.3); NW]);
        let films = vec![crate::structure::Layer::film(50.0, "H")];

        let mut amb = air(0.0);
        amb.nk = vec![Complex64::new(1.33, 0.02); NW].into();
        let mut absorbing_sub = sub(0.0);
        absorbing_sub.nk = vec![Complex64::new(1.52, 0.1); NW].into();

        let (stack, warns) = DesignStack::from_design(
            amb,
            absorbing_sub,
            &films,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();

        assert_eq!(warns.len(), 1, "{warns:?}");
        assert!(
            warns[0].contains("incident medium (layer 0) is absorbing"),
            "{}",
            warns[0]
        );
        assert!(
            warns[0].contains(&format!("{NW} of {NW} wavelengths")),
            "{}",
            warns[0]
        );

        // Layer 0 is flattened; the real part is kept, not replaced.
        assert!(stack.ambient().nk.iter().all(|z| z.im == 0.0));
        assert!(
            stack
                .ambient()
                .nk
                .iter()
                .all(|z| (z.re - 1.33).abs() < 1e-15)
        );

        // Nothing else is touched. An absorbing substrate and absorbing films
        // are ordinary physics -- the gate is about layer 0 only, and a gate
        // that quietly flattened the rest would be far worse than no gate.
        assert!(
            stack
                .substrate()
                .nk
                .iter()
                .all(|z| (z.im - 0.1).abs() < 1e-15)
        );
        assert!(
            stack.films()[0]
                .nk
                .iter()
                .all(|z| (z.im - 0.3).abs() < 1e-15)
        );

        // Once is enough: the ambient is private with no mutator, so the
        // stack operations the pipeline runs thousands of times cannot put
        // the absorption back. This is the whole argument for gating at
        // assembly instead of per merit evaluation.
        let mut s = stack;
        s.insert_needle_seed(0, 25.0, h(0.0)).unwrap();
        s.merge_adjacent();
        s.clamp_all(0.5, 500.0, false).unwrap();
        s.set_thickness(0, 12.0).unwrap();
        assert!(s.ambient().nk.iter().all(|z| z.im == 0.0));

        // A single stray wavelength is counted, not rounded up to the row.
        let mut spotty = air(0.0);
        let mut grid = vec![Complex64::new(1.0, 0.0); NW];
        grid[3] = Complex64::new(1.0, 1e-12);
        spotty.nk = grid.into();
        let (_, warns) = DesignStack::from_design(
            spotty,
            sub(0.0),
            &films,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        assert_eq!(warns.len(), 1);
        assert!(
            warns[0].contains(&format!("1 of {NW} wavelengths")),
            "{}",
            warns[0]
        );
        assert!(warns[0].contains("first at index 3"), "{}", warns[0]);

        // And a transparent ambient stays silent -- the common path.
        let (_, warns) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &films,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        assert!(warns.is_empty(), "{warns:?}");
    }

    fn h(d: f64) -> LayerSpec {
        LayerSpec::constant("H", 2.35, 0.0, d, NW)
    }

    fn l(d: f64) -> LayerSpec {
        LayerSpec::constant("L", 1.46, 0.0, d, NW)
    }

    #[test]
    fn basic_properties() {
        let s = stack(vec![h(100.0), l(50.0)]);
        assert_eq!(s.film_count_public(), 2);
        assert_eq!(s.total_layer_count(), 4);
        assert!((s.total_thickness_nm() - 150.0).abs() < 1e-12);
        assert_eq!(s.num_wavs(), NW);
    }

    #[test]
    fn mismatched_nk_lengths_rejected() {
        let mut bad = h(10.0);
        bad.nk = vec![cplx(2.35, 0.0); NW - 1].into();
        assert!(DesignStack::with_films(air(0.0), sub(0.0), vec![bad]).is_err());
    }

    #[test]
    fn insert_needle_seed_splits_host() {
        // Python: [top=depth, seed, bot=d_total-depth]; host flags preserved.
        let mut s = stack(vec![h(100.0)]);
        let seed = LayerSpec::constant("L", 1.46, 0.0, 5.0, NW);
        s.insert_needle_seed(0, 60.0, seed).unwrap();

        let f = s.films();
        assert_eq!(f.len(), 3);
        assert_eq!(f[0].material.as_ref(), "H");
        assert!((f[0].d_nm - 60.0).abs() < 1e-12);
        assert_eq!(f[1].material.as_ref(), "L");
        assert!((f[1].d_nm - 5.0).abs() < 1e-12);
        assert_eq!(f[2].material.as_ref(), "H");
        assert!((f[2].d_nm - 40.0).abs() < 1e-12);
        // total preserved
        assert!((s.total_thickness_nm() - 105.0).abs() < 1e-12);
        // host flags kept on fragments
        assert!(f[0].optimize && f[0].needle);
        assert!(f[2].optimize && f[2].needle);
    }

    #[test]
    fn insert_out_of_range_errors() {
        let mut s = stack(vec![h(100.0)]);
        assert!(s.insert_needle_seed(3, 5.0, l(1.0)).is_err());
    }

    #[test]
    fn merge_adjacent_sums_and_keeps_first_props() {
        // H(50) L(30) H(20) H(40) L(60) → merges: one pair of H's.
        // First H keeps its own props; give the first fragment different
        // flags to verify "keeps FIRST layer's properties".
        let mut first_h = h(20.0);
        first_h.optimize = false;
        let second_h = h(40.0);

        let mut s = stack(vec![h(50.0), l(30.0), first_h, second_h, l(60.0)]);
        let n = s.merge_adjacent();
        assert_eq!(n, 1);

        let f = s.films();
        assert_eq!(f.len(), 4);
        assert_eq!(f[2].material.as_ref(), "H");
        assert!((f[2].d_nm - 60.0).abs() < 1e-12);
        // first-of-pair properties win
        assert!(!f[2].optimize);
        assert!(f[2].needle);
        assert!((s.total_thickness_nm() - 200.0).abs() < 1e-12);
    }

    /// nk-keyed merge: same material with different nk (interface slice,
    /// group scaling, graded sublayers) must NEVER collapse.
    #[test]
    fn merge_adjacent_respects_nk() {
        // Same name, perturbed nk → no merge.
        let a = h(10.0);
        let mut b = h(10.0);
        b.nk = std::sync::Arc::from(vec![Complex64::new(2.5, 0.0); NW]);
        let mut s = stack(vec![a, b]);
        assert_eq!(s.merge_adjacent(), 0);
        assert_eq!(s.films().len(), 2);
        // Slice-like row (same name, mixed nk, derived flags) + bulk → kept.
        let mut slice = h(4.0);
        slice.optimize = false;
        slice.needle = false;
        slice.nk = std::sync::Arc::from(vec![Complex64::new(1.9, 0.0); NW]);
        let mut s = stack(vec![slice, h(46.0)]);
        assert_eq!(s.merge_adjacent(), 0);
        assert!(!s.films()[0].optimize);
    }

    #[test]
    fn merge_full_run_collapses() {
        let mut s = stack(vec![h(1.0), h(2.0), h(3.0)]);
        assert_eq!(s.merge_adjacent(), 2);
        assert_eq!(s.film_count_public(), 1);
        assert!((s.films()[0].d_nm - 6.0).abs() < 1e-12);
    }

    #[test]
    fn clamp_removes_thin_caps_thick() {
        let mut s = stack(vec![h(1.0), l(5000.0), h(50.0), l(2.0 - 1e-9)]);
        let rep = s.clamp_all(2.0, 1000.0, false).unwrap();
        assert_eq!((rep.rows_removed, rep.spans_capped), (2, 1));
        let f = s.films();
        assert_eq!(f.len(), 2);
        assert!((f[0].d_nm - 1000.0).abs() < 1e-12);
        assert!((f[1].d_nm - 50.0).abs() < 1e-12);
    }

    #[test]
    fn clamp_boundary_values_survive() {
        // Exactly-at-boundary layers survive untouched (strict < and >).
        let mut s = stack(vec![h(2.0), l(1000.0)]);
        let rep = s.clamp_all(2.0, 1000.0, false).unwrap();
        assert_eq!((rep.rows_removed, rep.spans_capped), (0, 0));
        assert_eq!(s.film_count_public(), 2);
    }

    #[test]
    fn remove_film_returns_layer() {
        let mut s = stack(vec![h(10.0), l(20.0)]);
        let removed = s.remove_film(0).unwrap();
        assert_eq!(removed.material.as_ref(), "H");
        assert_eq!(s.film_count_public(), 1);
        assert!(s.remove_film(5).is_err());
    }

    #[test]
    fn solver_arrays_layout_matches_smatrix_convention() {
        // wav-major re/im interleave: base = w * n_layers * 2 + slot * 2
        let s = stack(vec![h(10.0)]);
        let sa = s.solver_arrays();
        assert_eq!(sa.n_layers, 3);
        assert_eq!(sa.n_stack_cache.len(), NW * 3 * 2);
        assert_eq!(sa.thicknesses.len(), 3);

        // ambient slot 0, film slot 1, substrate slot 2
        assert!((sa.thicknesses[0]).abs() < 1e-15);
        assert!((sa.thicknesses[1] - 10.0).abs() < 1e-12);
        assert!((sa.thicknesses[2]).abs() < 1e-15);

        for w in 0..NW {
            let base = w * 3 * 2;
            assert!((sa.n_stack_cache[base] - 1.0).abs() < 1e-12); // air Re
            assert!(sa.n_stack_cache[base + 1].abs() < 1e-12); // air Im
            assert!((sa.n_stack_cache[base + 2] - 2.35).abs() < 1e-12); // H Re
            assert!((sa.n_stack_cache[base + 4] - 1.52).abs() < 1e-12); // sub Re
        }
        assert!(sa.incoherent_flags.iter().all(|&v| v == 0));
    }

    #[test]
    fn incoherent_flag_and_roughness_passthrough() {
        let mut sp = LayerSpec::constant("spacer", 1.45, 0.0, 1_000_000.0, NW);
        sp.coherent = false;
        sp.rough_type = 2;
        sp.rough_val = 7.5;
        let s = stack(vec![sp]);
        let sa = s.solver_arrays();
        assert_eq!(sa.incoherent_flags, vec![0, 1, 0]);
        assert_eq!(sa.rough_types[1], 2);
        assert!((sa.rough_vals[1] - 7.5).abs() < 1e-12);
    }

    // ------------------------------------------------------------------
    // F0.1 — span bookkeeping on the mutators and both constructors
    // ------------------------------------------------------------------

    /// Span summary: (start, end, logical, slice, bulk_start).
    fn span_rows(s: &DesignStack) -> Vec<(usize, usize, usize, bool, usize)> {
        s.spans()
            .iter()
            .map(|sp| (sp.start, sp.end, sp.logical, sp.slice, sp.bulk_start))
            .collect()
    }

    /// from_design keeps the spans expand emitted (they died at the door
    /// before F0.1). The 2-row span here is N1's interface case.
    #[test]
    fn from_design_keeps_spans() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.0, 0.0); NW]);
        let mut films = vec![crate::structure::Layer::film(20.0, "H")];
        films.push(crate::structure::Layer {
            interface: true,
            interface_thickness: 1.0,
            ..crate::structure::Layer::film(50.0, "H")
        });
        let (stack, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &films,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        assert_eq!(stack.films().len(), 3);
        assert_eq!(
            span_rows(&stack),
            vec![(0, 1, 0, false, 0), (1, 3, 1, true, 2)]
        );
        // The 2-row span is still singleton-bulk (N1).
        assert!(stack.span_of_row(1).is_singleton_bulk());
        assert!(stack.span_of_row(2).is_singleton_bulk());
    }

    /// N1 twin — the one that would fail if anyone rewrote the predicate
    /// as `end - start > 1`. An interface-carrying plain film with
    /// needle = true: needle sites are still generated inside its bulk,
    /// and the site set is bit-identical to the flag-only scan (the
    /// before/after claim of F0.1). The bulk row is also still a
    /// thin-removal candidate (the same filter `remove_thin_layers`
    /// applies; that function is untouched by this item).
    #[test]
    fn n1_interface_film_bulk_is_still_a_host() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.0, 0.0); NW]);
        let films = vec![
            crate::structure::Layer::film(100.0, "H"),
            crate::structure::Layer {
                interface: true,
                interface_thickness: 1.0,
                ..crate::structure::Layer::film(50.0, "H")
            },
        ];
        let (stack, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &films,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        // Row layout: 0 = plain bulk, 1 = slice (needle=false), 2 = bulk.
        let with_spans = crate::smatrix::synthesis::needle_pass::build_scan_sites(
            stack.films(),
            stack.spans(),
            10.0,
        );
        let flag_only =
            crate::smatrix::synthesis::needle_pass::build_scan_sites(stack.films(), &[], 10.0);
        assert_eq!(
            with_spans, flag_only,
            "singleton-bulk scan must be unchanged"
        );
        // Sites inside the interface film's bulk (row 2, d = 49): interior
        // multiples of 10 below 49 -> k in 1..4.
        assert!(
            with_spans
                .iter()
                .any(|s| s.film_idx == 2 && s.depth_into_layer_nm == 10.0)
        );
        assert!(
            with_spans.iter().all(|s| s.film_idx != 1),
            "no sites on the slice"
        );
        // Thin-removal candidate: same predicate remove_thin_layers uses.
        assert!(stack.films()[2].optimize && stack.films()[2].d_nm < 2.0 * 25.0);
    }

    /// R3 twin: the refusal exists at the public door itself, names the
    /// material and the row range, and is asserted directly (not only
    /// through the scan). Singleton-bulk hosts still insert.
    #[test]
    fn needle_refusal_at_the_mutator_names_span_and_range() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.0, 0.0); NW]);
        let graded = vec![crate::structure::Layer {
            inhomogen: true,
            inh_delta: 0.2,
            ..crate::structure::Layer::film(50.0, "H")
        }];
        let bg: HashSet<String> = ["H".to_string()].into_iter().collect();
        let (mut stack, _) =
            DesignStack::from_design(air(0.0), sub(0.0), &graded, &nk, &HashMap::new(), &wl, &bg)
                .unwrap();
        // 11 pinned rows: never a host, and the error says which span.
        assert_eq!(stack.films().len(), 11);
        let err = stack.insert_needle_seed(3, 5.0, h(5.0)).unwrap_err();
        assert!(err.contains("multi-row span"), "{err}");
        assert!(err.contains("rows 0..11"), "{err}");
        assert!(err.contains("'H'"), "{err}");
        // Nothing changed.
        assert_eq!(stack.films().len(), 11);
        // A singleton-bulk host still inserts (span rule, not flag rule).
        let plain = vec![crate::structure::Layer::film(50.0, "H")];
        let (mut s2, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &plain,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        s2.insert_needle_seed(0, 25.0, h(5.0)).unwrap();
        assert_eq!(s2.films().len(), 3);
    }

    /// N2 merge twins: intra-span (graded rows that became identical),
    /// cross-span (two adjacent carriers), and the slice-absorption case.
    #[test]
    fn merge_bookkeeping_intra_and_cross_span() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.0, 0.0); NW]);
        // Intra-span: delta = 0 makes every sublayer's factors exactly 1.0,
        // so all rows carry bit-identical nk and merge to one row inside
        // ONE span.
        let zero = vec![crate::structure::Layer {
            inhomogen: true,
            inh_delta: 0.0,
            ..crate::structure::Layer::film(50.0, "H")
        }];
        let bg: HashSet<String> = ["H".to_string()].into_iter().collect();
        let (gstack, _) =
            DesignStack::from_design(air(0.0), sub(0.0), &zero, &nk, &HashMap::new(), &wl, &bg)
                .unwrap();
        let n_before = gstack.films().len();
        assert!(n_before > 1, "the case needs multiple rows");
        let mut g = gstack;
        assert_eq!(g.merge_adjacent(), n_before - 1);
        assert_eq!(span_rows(&g), vec![(0, 1, 0, false, 0)]);
        assert!((g.total_thickness_nm() - 50.0).abs() < 1e-12);

        // Cross-span: two adjacent carriers merge into the FIRST span;
        // the second is dropped.
        let two = vec![
            crate::structure::Layer::film(20.0, "H"),
            crate::structure::Layer::film(30.0, "H"),
        ];
        let (mut s2, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &two,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        assert_eq!(span_rows(&s2).len(), 2);
        assert_eq!(s2.merge_adjacent(), 1);
        assert_eq!(span_rows(&s2), vec![(0, 1, 0, false, 0)]);
        assert!((s2.total_thickness_nm() - 50.0).abs() < 1e-12);
    }

    /// The contrived-but-reachable absorption case: an interface slice
    /// whose mixed nk is bit-identical to its neighbours' merges away,
    /// and the span bookkeeping records honestly what happened to the
    /// rows. n = 1.0 makes the looyenga eps round trip exact (cbrt of 1
    /// is 1, cubed is 1, sqrt of 1 is 1), so the slice row's nk is
    /// bit-identical to the plain rows and the nk-keyed merge fires.
    #[test]
    fn merge_absorbing_a_slice_drops_the_flag_it_no_longer_carries() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(1.0, 0.0); NW]);
        let films = vec![
            crate::structure::Layer::film(20.0, "H"),
            crate::structure::Layer {
                interface: true,
                interface_thickness: 1.0,
                ..crate::structure::Layer::film(50.0, "H")
            },
        ];
        let (mut stack, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &films,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        // rows: 0 = H20, 1 = slice (nk = mix(H,H) = H), 2 = bulk.
        // One maximal run [0, 3): all three merge into row 0's identity.
        assert_eq!(stack.merge_adjacent(), 2);
        assert_eq!(stack.films().len(), 1);
        assert!((stack.films()[0].d_nm - 70.0).abs() < 1e-12);
        // The first span survives with the merged row; the second is gone.
        assert_eq!(span_rows(&stack), vec![(0, 1, 0, false, 0)]);
    }

    /// insert bookkeeping: the host span splits locally; the seed is its
    /// own singleton; the bottom keeps the host's logical identity; the
    /// slice flag follows whichever half kept the leading edge.
    #[test]
    fn insert_splits_the_span_like_the_rows() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.0, 0.0); NW]);
        let films = vec![
            crate::structure::Layer::film(20.0, "H"),
            crate::structure::Layer {
                interface: true,
                interface_thickness: 1.0,
                ..crate::structure::Layer::film(50.0, "H")
            },
        ];
        let (base, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &films,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();

        // Split at the BULK row (2 = span end - 1): slice and top stay
        // together, seed alone, bottom alone.
        let mut s = base.clone();
        s.insert_needle_seed(2, 10.0, l(5.0)).unwrap();
        assert_eq!(s.films().len(), 5);
        assert_eq!(
            span_rows(&s),
            vec![
                (0, 1, 0, false, 0),
                (1, 3, 1, true, 2),
                (3, 4, 3, false, 3),
                (4, 5, 1, false, 4),
            ]
        );
        assert!((s.films()[2].d_nm - 10.0).abs() < 1e-12);
        assert!((s.films()[4].d_nm - 39.0).abs() < 1e-12);

        // Split at the SLICE row (1 = span start): the re-thickened slice
        // stays the leading edge; the bottom portion joins the bulk in a
        // slice-carrying span.
        let mut s = base;
        s.insert_needle_seed(1, 10.0, l(5.0)).unwrap();
        assert_eq!(s.films().len(), 5);
        assert_eq!(
            span_rows(&s),
            vec![
                (0, 1, 0, false, 0),
                (1, 2, 1, true, 2),
                (2, 3, 2, false, 2),
                (3, 5, 1, true, 4),
            ]
        );
    }

    /// remove bookkeeping: singletons drop whole; multi-row spans shrink;
    /// the slice and bulk halves record what survived.
    #[test]
    fn remove_updates_the_span_partition() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.0, 0.0); NW]);
        let films = vec![
            crate::structure::Layer::film(20.0, "H"),
            crate::structure::Layer {
                interface: true,
                interface_thickness: 1.0,
                ..crate::structure::Layer::film(50.0, "H")
            },
        ];
        let (base, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &films,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();

        // Remove the bulk: the slice row is alone in its span (derived
        // row whose carrier is gone — recorded, not hidden).
        let mut s = base.clone();
        s.remove_film(2).unwrap();
        assert_eq!(span_rows(&s), vec![(0, 1, 0, false, 0), (1, 2, 1, true, 2)]);
        assert!(!s.span_of_row(1).is_singleton_bulk());

        // Remove the slice: the bulk row is a plain singleton again.
        let mut s = base.clone();
        s.remove_film(1).unwrap();
        assert_eq!(
            span_rows(&s),
            vec![(0, 1, 0, false, 0), (1, 2, 1, false, 1)]
        );
        assert!(s.span_of_row(1).is_singleton_bulk());

        // Remove the only row of a singleton: its span is dropped.
        let mut s = base;
        s.remove_film(0).unwrap();
        assert_eq!(span_rows(&s), vec![(0, 2, 1, true, 1)]);
    }

    /// clamp bookkeeping: rows are still judged one by one (the
    /// span-quantity floor/cap is F0.2), and the partition records what
    /// survived. The B1 defect is deliberately still visible here —
    /// F0.2 licence item 7 is the fix.
    #[test]
    fn clamp_records_surviving_rows_per_span() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("SiO2"), vec![Complex64::new(1.46, 0.0); NW]);
        nk.insert(Arc::from("TiO2"), vec![Complex64::new(2.35, 0.0); NW]);
        let films = vec![
            crate::structure::Layer::film(100.0, "SiO2"),
            crate::structure::Layer {
                interface: true,
                interface_thickness: 1.0,
                ..crate::structure::Layer::film(50.0, "TiO2")
            },
        ];
        let (mut stack, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &films,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        // rows: [100][1 slice][49 bulk]. F0.2 licence item 7 (B1): the
        // 1 nm slice is never a floor candidate on its own - it survives,
        // and the carrier keeps the nanometre the slice was carved from.
        // Assert on BOTH rows and the total: a film-count assertion alone
        // passes for the wrong reason if the slice merged into the bulk.
        let rep = stack.clamp_all(2.0, 1000.0, false).unwrap();
        assert_eq!((rep.rows_removed, rep.spans_capped), (0, 0));
        assert!(rep.is_empty());
        assert_eq!(stack.films().len(), 3);
        assert_eq!(
            span_rows(&stack),
            vec![(0, 1, 0, false, 0), (1, 3, 1, true, 2)]
        );
        assert!((stack.films()[1].d_nm - 1.0).abs() < 1e-12);
        assert!((stack.films()[2].d_nm - 49.0).abs() < 1e-12);
        assert!((stack.total_thickness_nm() - 150.0).abs() < 1e-12);

        // Control: interface_thickness = 3.0 holds identically in both
        // directions.
        let films3 = vec![
            crate::structure::Layer::film(100.0, "SiO2"),
            crate::structure::Layer {
                interface: true,
                interface_thickness: 3.0,
                ..crate::structure::Layer::film(50.0, "TiO2")
            },
        ];
        let (mut c3, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &films3,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        let rep3 = c3.clamp_all(2.0, 1000.0, false).unwrap();
        assert!(rep3.is_empty());
        assert!((c3.total_thickness_nm() - 150.0).abs() < 1e-12);

        // A fully-deleted span disappears from the partition.
        let graded = vec![crate::structure::Layer {
            inhomogen: true,
            inh_delta: 0.2,
            ..crate::structure::Layer::film(50.0, "TiO2")
        }];
        let bg: HashSet<String> = ["TiO2".to_string()].into_iter().collect();
        let (mut g, _) =
            DesignStack::from_design(air(0.0), sub(0.0), &graded, &nk, &HashMap::new(), &wl, &bg)
                .unwrap();
        let n_rows = g.films().len();
        assert!(n_rows > 1);
        let rep = g.clamp_all(100.0, 1000.0, false).unwrap();
        assert_eq!(rep.rows_removed, n_rows);
        assert_eq!(rep.spans_removed.len(), 1);
        assert!(
            rep.spans_removed[0].starts_with("TiO2"),
            "{}",
            rep.spans_removed[0]
        );
        assert!(g.films().is_empty());
        assert!(g.spans().is_empty());
    }

    // ------------------------------------------------------------------
    // F0.2 — the span-quantity floor and cap (licence twins)
    // ------------------------------------------------------------------

    /// Licence item 1, measured (A2): a background-pinned 5 nm graded
    /// film (delta = 0.1 -> four 1.25 nm sublayers) under the default
    /// floor. Before F0.2 the film vanished row by row; now it survives
    /// whole. Assert film count AND total thickness, not row count.
    #[test]
    fn f02_floor_leaves_a_sub_floor_graded_span_whole() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("TiO2"), vec![Complex64::new(2.35, 0.0); NW]);
        let graded = vec![crate::structure::Layer {
            inhomogen: true,
            inh_delta: 0.1,
            ..crate::structure::Layer::film(5.0, "TiO2")
        }];
        let bg: HashSet<String> = ["TiO2".to_string()].into_iter().collect();
        let (mut stack, _) =
            DesignStack::from_design(air(0.0), sub(0.0), &graded, &nk, &HashMap::new(), &wl, &bg)
                .unwrap();
        assert_eq!(stack.films().len(), 4);
        let rep = stack.clamp_all(2.0, 1000.0, false).unwrap();
        assert!(rep.is_empty(), "{rep:?}");
        assert_eq!(stack.films().len(), 4, "the film must survive whole");
        assert!((stack.total_thickness_nm() - 5.0).abs() < 1e-12);
    }

    /// Licence item 1, reported: the floor above the span total removes
    /// the span AS A UNIT and the report names it. The assertion is on
    /// the report - a silent correct deletion and a silent wrong one
    /// look identical from the outside.
    #[test]
    fn f02_floor_above_a_span_removes_it_whole_and_names_it() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("TiO2"), vec![Complex64::new(2.35, 0.0); NW]);
        let graded = vec![crate::structure::Layer {
            inhomogen: true,
            inh_delta: 0.1,
            ..crate::structure::Layer::film(5.0, "TiO2")
        }];
        let bg: HashSet<String> = ["TiO2".to_string()].into_iter().collect();
        let (mut stack, _) =
            DesignStack::from_design(air(0.0), sub(0.0), &graded, &nk, &HashMap::new(), &wl, &bg)
                .unwrap();
        let rep = stack.clamp_all(6.0, 1000.0, false).unwrap();
        assert_eq!(rep.rows_removed, 4);
        assert_eq!(rep.spans_removed.len(), 1);
        assert!(
            rep.spans_removed[0].starts_with("TiO2 (5.0 nm)"),
            "{}",
            rep.spans_removed[0]
        );
        assert!(stack.films().is_empty());
    }

    /// Licence item 2, measured: a PINNED 1000 nm graded span
    /// (delta = 0.5 -> 57 rows) under a 300 nm ceiling refuses -
    /// message names the span material and both numbers. F1.6 narrows
    /// this to pinned spans; a scalable span is capped instead
    /// (f16_scalable_span_above_ceiling_is_capped_not_refused).
    #[test]
    fn f02_cap_refuses_a_multirow_span_above_the_ceiling() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("TiO2"), vec![Complex64::new(2.35, 0.0); NW]);
        let graded = vec![crate::structure::Layer {
            inhomogen: true,
            inh_delta: 0.5,
            ..crate::structure::Layer::film(1000.0, "TiO2")
        }];
        let bg: HashSet<String> = ["TiO2".to_string()].into_iter().collect();
        let (mut stack, _) =
            DesignStack::from_design(air(0.0), sub(0.0), &graded, &nk, &HashMap::new(), &wl, &bg)
                .unwrap();
        assert_eq!(stack.films().len(), 57, "the plan's measured row count");
        let err = stack.clamp_all(2.0, 300.0, false).unwrap_err();
        // The exact message (C3): deterministic format! output over
        // pinned inputs - pin it whole, so a lost line-continuation can
        // never regrow a whitespace run.
        assert_eq!(
            err,
            "clamp_all: span 'TiO2' is 1000.0 nm thick, above the 300.0 nm ceiling - \
             refusing rather than rescaling a profile"
        );
        // Nothing changed behind the error.
        assert_eq!(stack.films().len(), 57);
        assert!((stack.total_thickness_nm() - 1000.0).abs() < 1e-12);
        // Control: under a 1500 nm ceiling the same span passes.
        let (mut stack2, _) =
            DesignStack::from_design(air(0.0), sub(0.0), &graded, &nk, &HashMap::new(), &wl, &bg)
                .unwrap();
        let rep = stack2.clamp_all(2.0, 1500.0, false).unwrap();
        assert!(rep.is_empty());
    }

    // ------------------------------------------------------------------
    // F1.6 - the clamp's span branches
    // ------------------------------------------------------------------

    /// Licence item 4, landed: an under-thickness SCALABLE span (the
    /// optimized profiled carrier F1.6's from_design branch produces)
    /// under ClampUpFinal is SCALED to the floor, fractions preserved -
    /// inverting F0.3's deferral deliberately. The scale-up is a
    /// clamp-up, so it is not a ClampReport entry (F0.3's rule).
    #[test]
    fn f16_under_thickness_scalable_span_is_scaled_to_the_floor() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("TiO2"), vec![Complex64::new(2.35, 0.0); NW]);
        let graded = vec![crate::structure::Layer {
            inhomogen: true,
            inh_delta: 0.1,
            ..crate::structure::Layer::film(5.0, "TiO2")
        }];
        // NOT background: the carrier keeps optimize=true (F1.6) and the
        // span is scalable (all bulk rows free).
        let (mut stack, warns) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &graded,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        assert!(warns.is_empty());
        let rows: Vec<usize> = (0..stack.films().len()).collect();
        let d0: f64 = rows.iter().map(|&r| stack.films()[r].d_nm).sum();
        let fractions: Vec<f64> = rows.iter().map(|&r| stack.films()[r].d_nm / d0).collect();
        assert!(stack.span_is_scalable(stack.spans().first().unwrap()));

        // ClampUpFinal, final pass: scaled to the floor, not removed.
        let rep = stack
            .clamp_all_policy(6.0, 1000.0, ThinLayerPolicy::ClampUpFinal, true)
            .unwrap();
        assert!(rep.is_empty(), "scale-ups are not clamp-report entries");
        assert_eq!(stack.films().len(), rows.len(), "no row removed");
        let d1: f64 = rows.iter().map(|&r| stack.films()[r].d_nm).sum();
        assert!(
            (d1 - 6.0).abs() < 1e-9,
            "the span total IS the floor ({d1})"
        );
        for (&r, &f) in rows.iter().zip(&fractions) {
            let f_now = stack.films()[r].d_nm / d1;
            assert!((f_now - f).abs() <= 1e-12 * f.abs().max(1e-12));
        }

        // Under Remove the same span is removed whole (F0.3 unchanged).
        let (mut stack2, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &graded,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        let rep = stack2
            .clamp_all_policy(6.0, 1000.0, ThinLayerPolicy::Remove, true)
            .unwrap();
        assert_eq!(rep.spans_removed.len(), 1);
        assert!(stack2.films().is_empty());
    }

    /// C2 (review finding 2, measured): an interface-carrying PLAIN film
    /// is ONE physical layer (N1) and clamps up instead of being deleted.
    /// The slice stays bitwise (B1: it belongs to the authored
    /// thickness) and the bulk row takes `min - slice`, so the span's
    /// slice-inclusive total lands on the floor. The review's case:
    /// 0.8 nm TiO2 with a 0.2 nm interface under a 2.0 nm floor was
    /// removed whole; now it is clamped and kept.
    #[test]
    fn c2_interface_film_clamps_up_instead_of_being_deleted() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.35, 0.0); NW]);
        nk.insert(Arc::from("L"), vec![Complex64::new(1.46, 0.0); NW]);
        // The interface slice is carved by the flag owner, which resolves
        // to Some(k) only from the second entry on - a lead film is
        // required (F1.6's interface twin documents the same).
        let mut lead = crate::structure::Layer::film(20.0, "L");
        lead.optimize = false;
        lead.needle = false;
        let mut carrier = crate::structure::Layer::film(0.8, "H");
        carrier.interface = true;
        carrier.interface_thickness = 0.2;
        let (mut stack, warns) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &[lead.clone(), carrier.clone()],
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        assert!(warns.is_empty(), "a plain film is not profiled: no warning");
        let sp = stack.spans()[1];
        assert!(sp.slice, "the carrier carried its interface");
        assert!(sp.is_singleton_bulk(), "one physical layer, two rows");
        let slice_d = stack.films()[sp.start].d_nm;
        assert_eq!(slice_d, 0.2);
        let bulk_before = stack.films()[sp.bulk_start].d_nm;
        assert_eq!(
            bulk_before,
            0.8 - slice_d,
            "the authored 0.8 includes the slice"
        );

        // ClampUpFinal: clamped, kept, not reported (F0.3's rule).
        let rep = stack
            .clamp_all_policy(2.0, 300.0, ThinLayerPolicy::ClampUpFinal, true)
            .unwrap();
        assert!(rep.is_empty(), "clamp-ups are not clamp-report entries");
        assert_eq!(stack.films().len(), 3, "no row removed");
        assert_eq!(stack.films()[sp.start].d_nm, slice_d, "slice bitwise");
        let bulk = stack.films()[sp.bulk_start].d_nm;
        assert_eq!(bulk, 2.0 - slice_d, "the bulk takes the remainder");
        assert!(
            (bulk + slice_d - 2.0).abs() < 1e-9,
            "the span total IS the floor (2.0 - 0.2 need not round-trip bitwise)"
        );
        stack.assert_spans_partition();

        // ClampUpAlways: the same clamp through the every-pass posture.
        let (mut stack2, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &[lead.clone(), carrier.clone()],
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        let rep = stack2
            .clamp_all_policy(2.0, 300.0, ThinLayerPolicy::ClampUpAlways, true)
            .unwrap();
        assert!(rep.is_empty());
        assert_eq!(stack2.films().len(), 3);
        assert_eq!(stack2.films()[sp.start].d_nm, slice_d);
        assert_eq!(stack2.films()[sp.bulk_start].d_nm, 2.0 - slice_d);

        // Control: under Remove the same span is removed whole, reported
        // (F0.2's rule, unchanged) - rows 2 (slice + bulk).
        let (mut stack3, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &[lead, carrier],
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        let rep = stack3
            .clamp_all_policy(2.0, 300.0, ThinLayerPolicy::Remove, true)
            .unwrap();
        assert_eq!(rep.spans_removed.len(), 1);
        assert_eq!(rep.rows_removed, 2, "the whole carrier goes, slice too");
        assert_eq!(stack3.films().len(), 1, "only the lead survives");
    }

    /// C2 control: a plain interface-carrying film ABOVE the floor is
    /// untouched - the clamp-up branch is not reachable for a healthy
    /// span, and the slice rows never move through clamp either way.
    #[test]
    fn c2_interface_film_above_the_floor_is_untouched() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.35, 0.0); NW]);
        nk.insert(Arc::from("L"), vec![Complex64::new(1.46, 0.0); NW]);
        let mut lead = crate::structure::Layer::film(20.0, "L");
        lead.optimize = false;
        lead.needle = false;
        let mut carrier = crate::structure::Layer::film(50.0, "H");
        carrier.interface = true;
        carrier.interface_thickness = 0.2;
        let (mut stack, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &[lead, carrier],
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        let rep = stack
            .clamp_all_policy(2.0, 300.0, ThinLayerPolicy::ClampUpFinal, true)
            .unwrap();
        assert!(rep.is_empty());
        assert_eq!(stack.films().len(), 3);
        assert_eq!(stack.films()[1].d_nm, 0.2);
        assert_eq!(stack.films()[2].d_nm, 50.0 - 0.2);
    }

    /// A SCALABLE span above the ceiling is CAPPED as a bound (fractions
    /// preserved), not refused - the F0.2 refusal narrowed to pinned
    /// spans, whose twin above is unchanged.
    #[test]
    fn f16_scalable_span_above_ceiling_is_capped_not_refused() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("TiO2"), vec![Complex64::new(2.35, 0.0); NW]);
        let graded = vec![crate::structure::Layer {
            inhomogen: true,
            inh_delta: 0.5,
            ..crate::structure::Layer::film(1000.0, "TiO2")
        }];
        let (mut stack, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &graded,
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        assert_eq!(stack.films().len(), 57);
        let rep = stack
            .clamp_all_policy(2.0, 300.0, ThinLayerPolicy::Remove, false)
            .unwrap();
        assert_eq!(rep.spans_capped, 1);
        let d: f64 = (0..stack.films().len())
            .map(|r| stack.films()[r].d_nm)
            .sum();
        assert!(
            (d - 300.0).abs() < 1e-9,
            "the span total IS the ceiling ({d})"
        );
        // Uniform scaling: every row kept its share of the span.
        let f0 = stack.films()[0].d_nm;
        for f in stack.films() {
            assert!((f.d_nm / f0 - 1.0).abs() < 1e-9);
        }
    }

    /// The interface twin: a scalable span's slice row is an interface
    /// property, not part of the layer's thickness - it is NOT scaled,
    /// and the span partition still holds after the writeback.
    #[test]
    fn f16_interface_slice_is_not_scaled() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.35, 0.0); NW]);
        nk.insert(Arc::from("L"), vec![Complex64::new(1.46, 0.0); NW]);
        // A lead entry: the interface slice is carved by the FLAG OWNER,
        // which resolves to Some(k) only from the second entry on.
        // Pinned (optimize=false) so it is not a second LM parameter.
        let mut lead = crate::structure::Layer::film(20.0, "L");
        lead.optimize = false;
        lead.needle = false;
        let mut carrier = crate::structure::Layer::film(100.0, "H");
        carrier.inhomogen = true;
        carrier.inh_delta = 0.2;
        carrier.interface = true;
        carrier.interface_thickness = 5.0;
        let (mut stack, warns) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &[lead, carrier],
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        assert!(warns.is_empty());
        let sp = stack.spans()[1];
        assert!(sp.slice, "the carrier carried its interface");
        assert!(stack.span_is_scalable(&sp));
        // The slice row IS the span's first row (bulk_start = start + 1).
        let slice_row = sp.start;
        let slice_d = stack.films()[slice_row].d_nm;
        assert_eq!(slice_d, 5.0);

        // The 40% move through the writeback: only bulk rows scale.
        let params = build_params(&stack);
        assert_eq!(params.len(), 1);
        apply_params(&mut stack, &params, &[1.4 * 95.0]).unwrap();
        assert_eq!(stack.films()[slice_row].d_nm, slice_d, "slice unchanged");
        let bulk_d: f64 = (sp.bulk_start..sp.end).map(|r| stack.films()[r].d_nm).sum();
        assert!((bulk_d - 1.4 * 95.0).abs() < 1e-9, "the bulk total moved");
        stack.assert_spans_partition();
    }

    /// A scalable span is still not a needle host (the rule is unchanged
    /// in form; scalability only ever adds reasons a span stays whole).
    #[test]
    fn f16_scalable_span_is_still_not_a_needle_host() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.35, 0.0); NW]);
        let mut carrier = crate::structure::Layer::film(100.0, "H");
        carrier.inhomogen = true;
        carrier.inh_delta = 0.2;
        let (stack, _) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &[carrier],
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        let refusal = stack.needle_host_refusal(0);
        assert!(refusal.is_some(), "{refusal:?}");
        assert!(refusal.unwrap().contains("multi-row span"));
    }

    /// A one-row span above the ceiling is still CAPPED, not refused:
    /// D == d for a one-row span, and the no-span path is bit-identical,
    /// full stop.
    #[test]
    fn f02_one_row_cap_unchanged() {
        let mut s = stack(vec![h(5000.0)]);
        let rep = s.clamp_all(2.0, 1000.0, false).unwrap();
        assert_eq!(rep.spans_capped, 1);
        assert!(rep.spans_removed.is_empty());
        assert!((s.films()[0].d_nm - 1000.0).abs() < 1e-12);
    }

    // small helper so tests read like the Python property name
    impl DesignStack {
        fn film_count_public(&self) -> usize {
            self.films().len()
        }
    }

    // ------------------------------------------------------------------
    // F1.7 - profile refresh for the rate modes (U3/U4)
    // ------------------------------------------------------------------

    /// A gradient RateCapped carrier on the F1.6 posture (optimize=true:
    /// since F1.7 the rate profile is kept, not homogenized). `d` is the
    /// carrier's authored thickness; `f_start`/`rate` set the ramp.
    fn rate_gradient_stack(d: f64, rate: f64, with_interface: bool) -> DesignStack {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.35, 0.01); NW]);
        nk.insert(Arc::from("L"), vec![Complex64::new(1.46, 0.0); NW]);
        // A lead film: the interface slice is carved by the flag owner,
        // which resolves to Some(k) only from the second entry on; it
        // also gives the slice-mix a real "prev" operand.
        let mut lead = crate::structure::Layer::film(20.0, "L");
        lead.optimize = false;
        lead.needle = false;
        let mut carrier = crate::structure::Layer::film(d, "H");
        carrier.gradient = Some(crate::structure::GradientSpec::rate_capped(
            "H", "L", 0.0, rate, 100.0, 0.0, 1.0,
        ));
        if with_interface {
            carrier.interface = true;
            carrier.interface_thickness = 5.0;
        }
        let (stack, warns) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            &[lead, carrier],
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        assert!(
            warns.is_empty(),
            "the F1.6/F1.7 posture must not warn: {warns:?}"
        );
        stack
    }

    /// A legacy RateCapped carrier on the F1.6 posture.
    fn rate_legacy_stack(d: f64, rate: f64) -> DesignStack {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.35, 0.01); NW]);
        let mut carrier = crate::structure::Layer::film(d, "H");
        carrier.inhomogen = true;
        carrier.inh_mode = crate::structure::InhMode::RateCapped {
            rate,
            ref_thickness: 100.0,
            cap: 0.3,
        };
        let (stack, warns) = DesignStack::from_design(
            air(0.0),
            sub(0.0),
            std::slice::from_ref(&carrier),
            &nk,
            &HashMap::new(),
            &wl,
            &HashSet::new(),
        )
        .unwrap();
        assert!(
            warns.is_empty(),
            "the F1.6/F1.7 posture must not warn: {warns:?}"
        );
        stack
    }

    /// The stack's film rows as comparable bytes: material, nk (le bytes
    /// per component), thickness, and every flag. Two stacks are bitwise
    /// equal iff this is equal.
    fn films_bytes(s: &DesignStack) -> Vec<(String, Vec<u8>, f64, bool, i32, f64, bool, bool)> {
        s.films()
            .iter()
            .map(|f| {
                let mut nk_bytes = Vec::with_capacity(f.nk.len() * 16);
                for z in f.nk.iter() {
                    nk_bytes.extend_from_slice(&z.re.to_le_bytes());
                    nk_bytes.extend_from_slice(&z.im.to_le_bytes());
                }
                (
                    f.material.to_string(),
                    nk_bytes,
                    f.d_nm,
                    f.coherent,
                    f.rough_type,
                    f.rough_val,
                    f.optimize,
                    f.needle,
                )
            })
            .collect()
    }

    /// Refresh-correctness twin, gradient engine (full bitwise): a rate
    /// span scaled by hand from 100 nm to 200 nm, then refreshed, is
    /// BITWISE equal to the same span expanded from scratch at 200 nm.
    /// The gradient branch's totals are exact (the last sublayer absorbs
    /// the remainder, N7), so the refreshed extent is exactly the
    /// hand-set total and the comparison is meaningful to the last bit.
    /// Refresh and construction are the same code (B2) - this is a
    /// regression pin, not a load-bearing proof.
    #[test]
    fn f17_refresh_rebuilds_the_gradient_rate_span_bitwise() {
        let at_100 = rate_gradient_stack(100.0, 0.3, false);
        let si = 1; // the carrier span (lead film is span 0)
        let sp100 = at_100.spans()[si];
        assert_eq!(
            sp100.end - sp100.start,
            6,
            "count(100 nm) = ceil(100/min(20, 400/(10*2.35)))"
        );
        // Hand-scale: double every row (bitwise exact on the gradient
        // branch's binary-exact split) and pin the total.
        let mut scaled = at_100.clone();
        for r in sp100.start..sp100.end {
            let d = scaled.films()[r].d_nm * 2.0;
            scaled.set_thickness(r, d).unwrap();
        }
        let total: f64 = (sp100.start..sp100.end)
            .map(|r| scaled.films()[r].d_nm)
            .sum();
        assert_eq!(total, 200.0, "the hand-scale lands exactly on 200");
        scaled.refresh_profiles().unwrap();

        let at_200 = rate_gradient_stack(200.0, 0.3, false);
        assert_eq!(
            films_bytes(&scaled),
            films_bytes(&at_200),
            "refresh(100 -> 200) == construction(200), bitwise"
        );
        let sp200 = scaled.spans()[si];
        assert_eq!(
            sp200.end - sp200.start,
            12,
            "count(200 nm) doubles: the ceil re-derived at the refreshed extent"
        );
        assert_eq!(
            (sp100.logical, sp100.slice, sp100.bulk_start - sp100.start),
            (sp200.logical, sp200.slice, sp200.bulk_start - sp200.start),
            "span bookkeeping shape survives the refresh"
        );
        scaled.assert_spans_partition();
    }

    /// The same twin with an INTERFACE slice on the rate carrier: the
    /// slice row is re-carved at the new total and its mix nk is bitwise
    /// the construction's (the recipe's "prev" operand feeds the same
    /// looyenga_mix). The slice row's thickness is never scaled (F1.6).
    #[test]
    fn f17_refresh_recarves_the_interface_slice_bitwise() {
        let at_100 = rate_gradient_stack(100.0, 0.3, true);
        let si = 1;
        let sp = at_100.spans()[si];
        assert!(sp.slice, "the carrier carried its interface");
        let slice_d = at_100.films()[sp.start].d_nm;
        assert_eq!(slice_d, 5.0);
        let mut scaled = at_100.clone();
        for r in sp.start..sp.end {
            let d = scaled.films()[r].d_nm * 2.0;
            scaled.set_thickness(r, d).unwrap();
        }
        scaled.refresh_profiles().unwrap();
        let at_200 = rate_gradient_stack(200.0, 0.3, true);
        assert_eq!(
            films_bytes(&scaled),
            films_bytes(&at_200),
            "refresh with an interface slice == construction, bitwise"
        );
        let sp2 = scaled.spans()[si];
        assert!(sp2.slice);
        assert_eq!(
            scaled.films()[sp2.start].d_nm,
            5.0,
            "the slice is re-carved, not scaled"
        );
    }

    /// Refresh-correctness twin, legacy engine: the NK rows are bitwise
    /// the construction's (the ramp reads delta = rate*D/ref, and the
    /// refreshed D is the span's current extent); the THICKNESS rows
    /// agree to 1e-9 relative rather than bitwise, because the legacy
    /// branch's uniform split sums to the carrier only to float
    /// precision (N7) - the refreshed D is the float sum, within an ulp
    /// of the hand-set total. The plan's flat "bitwise" is corrected
    /// here in the same spirit as F1.6's d_r/D correction.
    #[test]
    fn f17_refresh_matches_construction_on_the_legacy_rate_span() {
        let at_100 = rate_legacy_stack(100.0, 0.25);
        let sp = at_100.spans()[0];
        let n100 = sp.end - sp.start;
        assert_eq!(n100, 16, "count(100 nm, rate 0.25): the F1.3 arithmetic");
        let mut scaled = at_100.clone();
        for r in sp.start..sp.end {
            let d = scaled.films()[r].d_nm * 2.0;
            scaled.set_thickness(r, d).unwrap();
        }
        scaled.refresh_profiles().unwrap();
        let at_200 = rate_legacy_stack(200.0, 0.25);
        assert_eq!(
            scaled.films().len(),
            at_200.films().len(),
            "the re-derived count matches construction"
        );
        for (a, b) in scaled.films().iter().zip(at_200.films().iter()) {
            assert_eq!(a.nk, b.nk, "the nk ramp is bitwise the construction's");
            assert_eq!(a.material, b.material);
            assert_eq!(a.coherent, b.coherent);
            assert!((a.rough_val - b.rough_val).abs() < 1e-12);
            let rel = ((a.d_nm - b.d_nm) / b.d_nm).abs();
            assert!(
                rel < 1e-9,
                "thickness within float-precision of construction ({rel})"
            );
        }
        scaled.assert_spans_partition();
    }

    /// Fixed-mode no-op twin (U3): refresh on Fixed and FixedSpan spans
    /// is bitwise identity, over a randomized differential in the style
    /// of `test_differential.py`'s seeded loops (300 seeds). The plan
    /// names that file; the synthesis-level randomized differential runs
    /// Rust-side because `refresh_profiles` is deliberately NOT public
    /// API (a fourth, user-called refresh point is exactly the risk the
    /// count twin guards against).
    #[test]
    fn f17_fixed_mode_refresh_is_bitwise_identity() {
        use rand::{Rng, SeedableRng, rngs::StdRng};
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let mut nk = HashMap::new();
        nk.insert(Arc::from("H"), vec![Complex64::new(2.35, 0.01); NW]);
        nk.insert(Arc::from("L"), vec![Complex64::new(1.46, 0.0); NW]);
        for seed in 0..300u64 {
            let mut rng = StdRng::seed_from_u64(seed);
            let pick = |rng: &mut StdRng, lo: f64, hi: f64| {
                lo + (hi - lo) * f64::from(rng.random::<u32>()) / f64::from(u32::MAX)
            };
            let n_films = 1 + (seed % 3) as usize;
            let mut films = Vec::with_capacity(n_films);
            for j in 0..n_films {
                let d = pick(&mut rng, 40.0, 280.0);
                let mut carrier = crate::structure::Layer::film(d, "H");
                match (seed + j as u64) % 3 {
                    0 => {
                        carrier.inhomogen = true;
                        carrier.inh_delta = pick(&mut rng, 0.05, 0.3);
                    }
                    1 => {
                        carrier.gradient = Some(crate::structure::GradientSpec::fixed_span(
                            "H",
                            "L",
                            pick(&mut rng, 0.0, 0.3),
                            pick(&mut rng, 0.7, 1.0),
                        ));
                    }
                    _ => {} // plain film
                }
                films.push(carrier);
            }
            let (stack, _) = DesignStack::from_design(
                air(0.0),
                sub(0.0),
                &films,
                &nk,
                &HashMap::new(),
                &wl,
                &HashSet::new(),
            )
            .unwrap();
            let before_films = films_bytes(&stack);
            let before_spans = stack.spans().to_vec();
            let mut refreshed = stack.clone();
            refreshed.refresh_profiles().unwrap();
            assert_eq!(
                films_bytes(&refreshed),
                before_films,
                "fixed-mode refresh is bitwise identity (seed {seed})"
            );
            assert_eq!(
                refreshed.spans(),
                &before_spans[..],
                "fixed-mode refresh leaves the partition alone (seed {seed})"
            );
        }
    }

    /// B3 errors-on refusal twin: a recipe carrying `apply_errors: true`
    /// is refused by name, instead of re-rolling the ensemble.
    /// Unreachable through `from_design` (it always expands
    /// deterministic); the twin constructs one directly, so the door is
    /// proved shut before anything opens it.
    #[test]
    fn f17_recipe_with_errors_is_refused() {
        let wl: Vec<f64> = (0..NW).map(|i| 400.0 + i as f64 * 50.0).collect();
        let films = vec![LayerSpec::constant("H", 2.35, 0.01, 100.0, NW)];
        let spans = vec![Span {
            start: 0,
            end: 1,
            logical: 0,
            slice: false,
            bulk_start: 0,
        }];
        let recipe = SpanRecipe {
            layer: crate::structure::Layer::film(100.0, "H"),
            groups: Arc::new(HashMap::new()),
            bulk_nk: films[0].nk.to_vec(),
            self_owned: false,
            prev_eff_nk: None,
            provider: Arc::new(crate::structure::DictProvider::new()),
            wavelengths: wl.clone(),
            opts: crate::structure::ExpandOptions {
                apply_errors: true,
                seed: None,
            },
            emitted_total: 100.0,
        };
        let mut stack =
            DesignStack::from_parts(air(0.0), sub(0.0), films, spans, vec![Some(recipe)]).unwrap();
        let err = match stack.refresh_profiles() {
            Err(e) => e,
            Ok(()) => panic!("expected the errors-on refusal"),
        };
        assert!(err.contains("apply_errors"), "{err}");
        assert!(err.contains("re-roll"), "{err}");
    }
}
