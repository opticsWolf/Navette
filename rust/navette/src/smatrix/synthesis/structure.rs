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

use std::sync::Arc;

use num_complex::Complex64;

use crate::smatrix::optics_core::cplx;
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
#[derive(Clone, Debug)]
pub struct DesignStack {
    ambient: LayerSpec,
    substrate: LayerSpec,
    films: Vec<LayerSpec>,
    spans: Vec<Span>,
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
        Self::from_parts(ambient, substrate, films, spans)
    }

    /// The one internal constructor: validated row lists plus the span
    /// partition that covers them. `from_design` passes the spans it kept
    /// from `expand`; every other door synthesizes singletons above.
    fn from_parts(
        ambient: LayerSpec,
        substrate: LayerSpec,
        films: Vec<LayerSpec>,
        spans: Vec<Span>,
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
        let s = DesignStack {
            ambient,
            substrate,
            films,
            spans,
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

    /// R3 / §6 row 6: the one needle-host admissibility rule, shared by
    /// `insert_needle_seed` (the public door) and the scan-site filter.
    /// A host row inside a multi-row span is refused — a needle must
    /// never split a graded profile or land on a derived row. The message
    /// names the span's material and its row range. `None` = admissible.
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
    /// F1.7 extends it to the recipe vector's length.
    pub(crate) fn assert_spans_partition(&self) {
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
        // rule to mixture gradients):
        // - BACKGROUND (named in `background`): expanded WITH the profile
        //   and pinned (optimize/needle forced false on the whole carrier
        //   span). True physics, fixed: needle never hosts there, LM
        //   skips the rows, nk-keyed merge and flag-guarded cleanup preserve
        //   the span. No warning - explicit opt-in, nothing dropped.
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
        let (sa, spans) = crate::structure::expand(
            &seq,
            &provider,
            wavelengths,
            groups,
            crate::structure::ExpandOptions::deterministic(),
        )?;
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
        Self::from_parts(ambient, substrate, rows, spans).map(|s| (s, warnings))
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
        let mut new_spans = Vec::with_capacity(old_spans.len() + 2);
        for sp in &old_spans {
            if sp.end <= film_idx {
                new_spans.push(*sp);
            } else if film_idx < sp.start {
                new_spans.push(Span {
                    start: sp.start + 2,
                    end: sp.end + 2,
                    logical: sp.logical,
                    slice: sp.slice,
                    bulk_start: sp.bulk_start + 2,
                });
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
            }
        }
        self.spans = new_spans;
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
        // Per output row: the span owning the run's leader row, and that
        // leader's original row index.
        let mut owners: Vec<(usize, usize)> = Vec::with_capacity(runs.len());
        for (s, e) in runs {
            let mut combined_d = films[s].d_nm;
            for r in s + 1..e {
                combined_d += films[r].d_nm;
            }
            let mut result = films[s].cloned();
            result.d_nm = combined_d;
            merged.push(result);
            owners.push((span_index_of_row(&old_spans, s), s));
        }

        // Rebuild the partition: consecutive output rows owned by the
        // same old span form one span. A span all of whose rows merged
        // into an earlier run disappears (the merged row joined the FIRST
        // span touched — "first layer's properties win", matching the
        // films rule). The slice flag survives only when the old span's
        // slice row IS its surviving group's first row; if the slice row
        // was absorbed into an earlier run, the flag would point at a
        // row that no longer is one.
        let mut new_spans: Vec<Span> = Vec::new();
        let mut p = 0usize;
        let mut g = 0usize;
        while g < owners.len() {
            let (sid, leader) = owners[g];
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
            p += h - g;
            g = h;
        }

        self.films = merged;
        self.spans = new_spans;
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
        let mut new_spans = Vec::with_capacity(old_spans.len());
        for sp in old_spans {
            if sp.end <= film_idx {
                new_spans.push(sp);
            } else if film_idx < sp.start {
                new_spans.push(Span {
                    start: sp.start - 1,
                    end: sp.end - 1,
                    logical: sp.logical,
                    slice: sp.slice,
                    bulk_start: sp.bulk_start - 1,
                });
            } else if sp.end - sp.start > 1 {
                // Containing span with rows surviving the removal.
                let slice_survives = sp.slice && film_idx != sp.start;
                new_spans.push(Span {
                    start: sp.start,
                    end: sp.end - 1,
                    logical: sp.logical,
                    slice: slice_survives,
                    bulk_start: sp.start + usize::from(slice_survives),
                });
            }
            // else: the removed row was the span's only row — dropped.
        }
        self.spans = new_spans;
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
    /// `Remove`/`ClampUpFinal`), `true` sets a surviving ONE-ROW span to
    /// `clamp_min_nm` instead (the final pass under `ClampUpFinal`, every
    /// pass under `ClampUpAlways`). A multi-row span is removed whole
    /// under every policy: clamping a span up is F1.6's scale operation
    /// and does not exist yet. Clamp-ups are neither removals nor caps,
    /// so they do not appear in the report.
    pub fn clamp_all(
        &mut self,
        min_nm: f64,
        max_nm: f64,
        clamp_up: bool,
    ) -> Result<ClampReport, String> {
        debug_assert!(min_nm >= 0.0 && max_nm > min_nm);
        // Refusals are checked BEFORE any mutation, so the stack is never
        // left half-clamped behind an error.
        for sp in &self.spans {
            if sp.end - sp.start <= 1 {
                continue;
            }
            let d: f64 = self.films[sp.start..sp.end].iter().map(|l| l.d_nm).sum();
            if d > max_nm {
                return Err(format!(
                    "clamp_all: span '{}' is {:.1} nm thick, above the {:.1} nm ceiling -                      refusing rather than rescaling a profile",
                    self.films[sp.start].material, d, max_nm
                ));
            }
        }

        let old = std::mem::take(&mut self.films);
        let old_spans = std::mem::take(&mut self.spans);
        let mut surviving: Vec<LayerSpec> = Vec::with_capacity(old.len());
        let mut new_spans: Vec<Span> = Vec::with_capacity(old_spans.len());
        let mut report = ClampReport::default();
        let mut p = 0usize;

        for sp in &old_spans {
            let d: f64 = old[sp.start..sp.end].iter().map(|l| l.d_nm).sum();
            if d < min_nm {
                // F0.3: clamp-up applies to one-row spans only; a
                // multi-row span's profile is scaled whole by F1.6 or,
                // until then, removed whole under every policy (the
                // ThinLayerPolicy doc comment states the deferral).
                if clamp_up && sp.end - sp.start == 1 {
                    let mut thin = old[sp.start..sp.end].to_vec();
                    thin[0].d_nm = min_nm;
                    surviving.extend(thin);
                    new_spans.push(Span {
                        start: p,
                        end: p + 1,
                        logical: sp.logical,
                        slice: sp.slice,
                        bulk_start: p + usize::from(sp.slice),
                    });
                    p += 1;
                    continue;
                }
                report
                    .spans_removed
                    .push(format!("{} ({:.1} nm)", old[sp.start].material, d));
                report.rows_removed += sp.end - sp.start;
                continue;
            }
            let mut rows = old[sp.start..sp.end].to_vec();
            if sp.end - sp.start == 1 && rows[0].d_nm > max_nm {
                // One-row span: D == d, so this is today's row cap.
                rows[0].d_nm = max_nm;
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
            p = e;
        }

        self.films = surviving;
        self.spans = new_spans;
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

    /// Licence item 2, measured: a pinned 1000 nm graded span
    /// (delta = 0.5 -> 57 rows) under a 300 nm ceiling refuses -
    /// message names the span material and both numbers.
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
        assert!(err.contains("'TiO2'"), "{err}");
        assert!(err.contains("1000.0"), "{err}");
        assert!(err.contains("300.0"), "{err}");
        assert!(err.contains("refusing"), "{err}");
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
}
