# Implementation plan — second audit: amendments & recommendations

STATUS: audit of `docs/implementation_plan.md` **revision 3** (commit
`baaad26`) against the same tree as the first audit. `baaad26` touches only
`docs/`, so the code base is still exactly `72a2d4d` (0.6.32) and every
citation below — like every citation in revision 3 — was checked against that
tree. Revision 2's audit (`docs/implementation_plan_amendment.md`) is
superseded as a worklist (its items are folded into the plan as `A*`/`N*`
blocks, several of them corrected); this file is the record of the second
pass. Findings here use the prefix **B\*** to stay disjoint from the plan's
own `A*`/`N*`/`U*` tags.

Method: every new claim revision 3 makes (the `N*` items, the `U*`
requirements, the corrected `A*` items, and the measured tables) was traced
to the tree. The ledger is §1. The findings are §2.

---

## 0. Verdict

Revision 3 is accurate at a level that changes what an audit is for. Of its
~40 new checkable claims, **every line-number citation and every measured
number checks out exactly** (§1) — including the F0.1 clamp-arithmetic table,
which reproduces `sub_layer_count`'s integer boundaries correctly at every
row, and the F1.6 analytic-Jacobian cost claim, which matches the deposits
doc comment verbatim. The two corrections it makes to the first audit (A1's
mechanism, A6's scope) are themselves verified: the program payload does live
inside `sections` with silent section-ignore (`config.rs:592-624`), and
`MixRule` (`mod.rs:54`) is the right selector with all six kernels.

The second pass found **one more live bug of the same class the plan already
documents (B1), two genuine design gaps in F1.7 (B2, B3), one licence
self-contradiction in F0.2 (B4), one scope gap in F1.6's wording (B5), and
four small corrections (B6–B9, B10)**. Nothing re-sequences the ladder;
nothing invalidates a gate. The items are, in order:

| # | Finding | Class | Severity |
|---|---|---|---|
| B1 | `clamp_all` deletes sub-floor **interface-slice rows** today — A2's bug one row type over, which F0.2 fixes by implication but its own text says it must not | new live bug + F0.2 self-contradiction | high |
| B2 | F1.7's `refresh_profiles` has nothing to refresh from: no span carries its profile recipe | design gap | high |
| B3 | F1.7 does not say what happens to the legacy δ error draw on refresh (re-draw per cycle changes the Monte-Carlo ensemble; the refresh-correctness twin as written only holds with draws off) | design gap | medium |
| B4 | F0.2's licence "bit-identical, full stop" is violated by its own item 6 unless the clamp report is conditionally emitted | licence contradiction | medium |
| B5 | F1.6's scalable set says "graded carrier" while its own table includes `GradientMode::FixedSpan` — the gradient-carrier case is unaddressed in the change text and licence | scope gap | medium |
| B6 | F1.3's `sub_layer_count` caller list misses the caller inside `property_issues` itself (`layer.rs:227`) — the one place the count is user-visible text | small | low |
| B7 | "evaluator.rs:373 is the last statement of `optimize_thicknesses_report`" is wrong — and the truth is good news for F0.3: the reported merit is *already* evaluated after the clamp | small | low |
| B8 | F1.7's "`optimize_thicknesses` … runs seven times per macro cycle" conflates seven call *sites* with seven invocations | small | low |
| B9 | The tree already computes what the singleton-bulk predicate derives: `bulk_spans` (`expansion.rs:158,313`) dies at the same door as `spans` | recommendation | — |
| B10 | Two one-line citation drifts (`expansion.rs:257` vs the `:258` branch; `layer.rs:127` vs the `:128` read) | trivial | — |

---

## 1. Claims ledger (revision 3's new claims, verified)

| Plan claim | Verdict | Evidence |
|---|---|---|
| **N1**: a span includes its interface-slice row; `start` is captured before the slice is pushed; `slice: emitted_slice` | ✓ exact | `expansion.rs:209` (`let start = col_thick.len();`, before the "Plane slice" block at `:214`); `:307-312` (`spans.push(Span { start, end: col_thick.len(), logical: k, slice: emitted_slice })`) |
| N1: the bulk row carries `needle = carrier.needle && !is_slice`, `is_slice` only at `r == span.start` | ✓ exact | `structure.rs:266-268` |
| **N2**: `merge_adjacent` keys on `material == && nk ==`; `cleanup_design` runs it twice | ✓ exact | `structure.rs:346-347`; `cleanup.rs:110` (pass 1), `:116` (pass 3) |
| N2: all sublayers of one carrier share the carrier's material name | ✓ exact | `structure.rs:270` |
| **N3**: `inflate_design`'s fall-through selects all films; `round_to_qwot`'s clamp-up is `(round_half_even(ratio) * step).max(step)`; neither is pipeline-reachable | ✓ exact, and stronger than stated | `inflate.rs:135` (`_ => (0..n_films).collect(),`), `:186`; neither function appears in `_smatrix.pyi`, `navette-py`, or `pipeline.rs` — `inflate_design` is pipeline-gated (`pipeline.rs:30, 201`), `round_to_qwot` has no caller at all outside tests |
| **N4**: `Layer`'s serde is hand-written; `serialize_map(Some(13))`; "unknown keys ignored" | ✓ exact | `layer.rs:335` (`impl Serialize`), `:338` (`Some(13)`), `:356-357` (`impl Deserialize`, doc "version-checked, unknown keys ignored") |
| **N5**: `validation.rs` is 82 lines with only the issue type and gate; `property_issues` at `layer.rs:155`; it takes no provider | ✓ exact | `wc -l` = 82; `layer.rs:155` (`pub fn property_issues(&self, label: &str)`) |
| **N6**: README `:139`/`:170` are historical, `:190` is the bump site | ✓ exact | `:139` "skips the tree-wide reformat commit (0.6.32)"; `:170` "`cargo fmt --all --check` (since 0.6.32)"; `:190` wheel filename |
| **N7**: legacy divides uniformly (`step_t = thickness / sub`); gradient branch absorbs the remainder | ✓ exact | `expansion.rs:273` |
| **N8**: `run_needle`'s `layers` is a required positional | ✓ exact | `pipeline.py:188-189` |
| **N9**: state keys silently ignored; `set_properties` warns; `LayerRow` refuses; program sections silently ignored; program top-level refuses | ✓ exact, all four doors | `layer.rs:357` (doc); `layer.rs:252` (returns issues); `design_config.rs:74` (`#[serde(deny_unknown_fields)]`); `config.rs:592-624` (bare `sections.get(...)` — no section whitelist); `config.rs:297-307` (top-level refuse) |
| **U1**: the only clamp-up in the tree is `round_to_qwot`, no pipeline caller; `needle_seed_thickness_nm = 5.0` at `cycle.rs:47` | ✓ exact | `inflate.rs:186`; `cycle.rs:47` |
| **U2**: whole-layer scaling precedent is `Group::thick_factor` at `expansion.rs:128` | ✓ exact | `:128` (`layer.thickness * group.thick_factor + group.thick_summand`) |
| **U3**: `InhMode::Fixed`'s factors depend on δ and `i/(sub-1)` only; thickness enters solely as `step_t` | ✓ exact | `expansion.rs:265-267`, `:273` |
| **U5**: four rules see solver rows — floor, cap, layer budget, inflate | ✓ exact | `structure.rs:388-394` (both branches); `pipeline.rs:102` (`films().len() >= max_film_layers`); `inflate.rs:135` |
| F0.1's clamp table (3/5/6/10/10/20 nm × δ) and F1.1's <6 nm figures | ✓ all nine rows reproduce `sub_layer_count` exactly | recomputed against `layer.rs:126-133`; the tree's own test pins two of them (`layer.rs:427-434`: `graded(1000.0, 0.5) == 57`, `graded(100.0, 0.1) == 11`) |
| F0.2's "seven call sites" of `optimize_thicknesses` | ✓ exact | `cleanup.rs:87, :121`; `cycle.rs:107, :249`; `inflate.rs:148, :191`; `pipeline.rs:271` (evaluator hits are `#[cfg(test)]`) |
| F0.2: `check_budgets` at `:102`, defaults `max_film_layers: 40` (`config.rs:84`), `clamp_min_nm: 2.0` (`:88`), `stagnation_oscillation_ratio: 0.75` (`:100`), window 5 (`:98`); `PipelinePhaseResult.layer_count` exists | ✓ exact | all six verified |
| F0.2: `clamp_all`'s return value discarded at all four call sites | ✓ exact | `pipeline.rs:192, 213, 273`; `evaluator.rs:373` — all statement-position |
| F0.3: `lb` at `evaluator.rs:341`; `StagnationOscillation` variant and its `"STAGNATION_OSCILLATION"` string; `clamp_all`'s doc comment states remove-not-clamp | ✓ exact | `evaluator.rs:341`; `config.rs:18, 31`; `structure.rs:377-380` |
| F0.3: `round_to_qwot`'s `.max(step)` is the only clamp-up | ✓ exact | `inflate.rs:186` |
| F1.1 (A9.1 corrected): **five** setters call `gate_layer` | ✓ exact | `rust/navette-py/src/structure.rs` — `:239` (thickness), `:271` (inhomogen), `:294` (inh_delta), `:308` (roughness), `:330` (interface_thickness); the other seven setters and the two enum-checking setters (`set_rough_type`, `set_layer_type`) do not |
| F1.1 (A6 corrected): six kernels; `MixRule` at `mod.rs:54` with parameterized variants; Python dispatches all six at `__init__.py:232` | ✓ exact | `ema.rs:23, 37, 52, 70, 87, 136` (`general_power_law` at `:52`); `MixRule { Bruggeman { max_iter, tol }, MaxwellGarnett, Looyenga, Lichtenecker, MoriTanaka { l }, PowerLaw { alpha } }` derives only `Clone, Copy, Debug`; `materials/__init__.py:232-233` |
| F1.1 (N5): provider-existence refusal needs the provider | ✓ | `property_issues(&self, label)` has no provider parameter |
| F1.1 (A5): `_FILM_DEFAULTS` at `pipeline.py:85`; unknown-flag refusal fires twice; `ArrayFilm` one nk; provider entries film-name → nk at `driver.rs:124` | ✓ exact | `pipeline.py:85-88`, `:101-103`, `:112-114`, `:117`; `driver.rs:32-45`, `:124` |
| F1.3: `sub_layer_count` reads `self.inh_delta`; callers at `expansion.rs:257`, `structure.rs:272`, `architect.rs:565` | ✓ for the three named — **but see B6: there is a fourth** | `layer.rs:128`; `structure/structure.rs:272-273`; `architect.rs:565-566` |
| F1.6: `opt_indices` at `evaluator.rs:310`; residual closure at `:327`; `lb` at `:341`; deposits cost doc at `:58`; `assemble_jacobian` at `jacobian.rs:120` | ✓ exact | all five verified |
| F1.6: `optimize = true` on a graded film today means homogenize (background predicate at `driver.rs:130`, flatten at `structure.rs:236`) | ✓ exact | both verified |
| F1.7: refresh point 2 at the cycle top before the pre-flight (`pipeline.rs:139`); `DepositJacobian::fill` at `evaluator.rs:394`; the `Ok(None)` decline at `:406` | ✓ exact | all three verified |
| F2.2 (A8): `run()` at `pipeline.rs:122`, calls `clamp_all` three times; `evaluate_merit` at `evaluator.rs:284` | ✓ exact | verified |
| F2.4 (A1 corrected): program payload lives inside `sections`; unknown sections silently ignored; unknown top-level keys refuse | ✓ exact | `config.rs:295` (`_ => &["sections"]`), `:297-307`, `:580-633`; `program.py:45`, and the module docstring (`:7`) indeed claims "refused when missing/stale/future" |
| F2.4 (N8): `layers` positional; DesignStack mutators PyO3-exposed | ✓ | `pipeline.py:188`; `_smatrix.pyi:450-455` |
| F1.4 (A3): the policy comment in two places; `test_stale_schema_versions_refused` at `:51`; `test_state_fingerprint` at `:96`; `_accepted_version` at `test_request_bits.py:232` asserting `len(ok) == 1` | ✓ exact | `types.py:113-118`; `test_roundtrip.py:70-76, 51-66, 96-98`; `test_request_bits.py:232-240` |
| Ground rule 5 (R1): no in-tree bit-exact needle-run pin | ✓ consistent with two independent searches | the pin's "scratchpad harness outside the repo" half is the author's knowledge, accepted on authority — the in-tree absence is verified |

Everything else revision 3 asserts that was checked in the first audit
(version sites, harness counts, the `A*`-corrected items it adopts) remains
verified; the tree has not moved.

---

## 2. Findings

### B1 — `clamp_all` deletes sub-floor interface-slice rows: the same live-bug class as A2, one row type over — and F0.2's own text currently preserves it

**Evidence.** The chain is short and every link is in the tree:

1. Interface-slice rows are *films*. `from_design` pushes every row of every
   span, the slice included (`structure.rs:263-275`), with
   `d_nm = sa.thicknesses[r]` — the interface thickness — and flags zeroed
   (`optimize: carrier.optimize && !is_slice`).
2. `clamp_all` is unguarded over all films (`structure.rs:388-394`): the
   removal branch is `if layer.d_nm < min_nm { n_removed += 1; continue; }`.
   There is no slice check.
3. Validation permits legal sub-floor slices. `layer.rs:189` refuses only
   `interface_thickness >= thickness`; a 1.5 nm slice on a 10 nm film is a
   legal authoring state, and the default `interface_thickness` of 0.0 is
   also below any floor when `interface = true` is set with an explicit
   small width.
4. So today, on the design path, a film declared with `interface = true` and
   a slice thinner than `clamp_min_nm` loses its slice row at the **first**
   clamp of the run — after cleanup, after inflate, after every thickness
   optimization (`pipeline.rs:192, 213, 273`; `evaluator.rs:373`). The
   interface physics disappears silently; `clamp_all`'s discarded return
   value (F0.2 rule 4) means nobody is told. It is A2's bug with the other
   derived row type, and it is Python-reachable too (`_smatrix.pyi:455`
   exposes `clamp_all`).

**And F0.2's current text would keep it.** F0.2 says: *"For a singleton-bulk
span (every stack that exists today) row and span are the same object and
nothing changes."* An interface-carrying plain film is singleton-**bulk**
with **two** rows — the plan's own N1 finding — so for it "row and span are
the same object" is false, and "nothing changes" preserves the deletion.
(The parenthetical is also false in the other direction: pinned graded
spans exist today and are multi-row — the plan's own A2 analysis.) The
pieces of the fix are already in F0.2 — *"the slice is not user thickness
and never was"* — they are just not connected to the singleton-bulk case.

**Amendment (three sentences to F0.2):**

1. **Change:** the slice row is never a floor or cap candidate, in every
   branch, singleton-bulk or not. It is not user thickness; the span
   comparison already excludes it, and the per-row branches must too.
2. **Licence:** add item 7 — "a slice row thinner than `clamp_min_nm` is
   preserved (previously: deleted at the first clamp)". This is a behaviour
   change on stacks that use no new feature, which is exactly why it must be
   enumerated rather than slipped into the bulk rule.
3. **Twin:** a design film with `interface = true`, `interface_thickness =
   1.0`, under `clamp_all(2.0, 1000.0)` — before: the slice row is gone;
   after: it survives with the bulk row. Assert on both rows.

### B2 — `refresh_profiles` has nothing to refresh from

**Evidence.** F1.7's contract is: re-derive `n_sub`, re-evaluate `f(z)`,
re-run the EMA kernel, rebuild the span's rows in place. But what F1.1/F0.1
put on the stack is: `LayerSpec` rows (material, nk, d_nm, flags) and
`Vec<Span>` (four integer/bool fields — `expansion.rs:31-38`). **Neither
carries the recipe**: the profile spec (`InhMode::RateCapped { rate,
ref_thickness, cap }` or `GradientMode::RateCapped` + `GradientSpec`), the
resolved endpoint spectra `nk_a(λ)`/`nk_b(λ)` (or the legacy base nk and δ),
the group summand, or the frozen error draw. `DesignStack` has no provider
either. "Re-run the EMA kernel" has no inputs on the object being refreshed.

**Amendment.** F1.7 (or F1.6, which already introduces span-carried state)
adds per-span profile provenance, and says where:

- Recommended shape: a side table keyed by logical index —
  `profiles: BTreeMap<usize, SpanProfile>` — populated where the recipe is
  already in hand (`from_design`/`expand`, where `bulk_start` is computed at
  `expansion.rs:256/:313`), maintained by the same five mutators F0.1
  already bookkeeps, and carrying the resolved spectra so refresh needs no
  provider (the grid is fixed per run). This keeps `Span` `Copy` (B9's
  cheaper cousin of the same observation) and costs nothing for plain
  stacks.
- The recipe must include the **frozen error contribution** (B3) and the
  group summand, or the refreshed profile is not the same function the
  constructor ran.

### B3 — the δ error draw on refresh is unspecified, and one twin as written cannot pass

**Evidence.** The legacy branch draws the δ error at expansion
(`expansion.rs:260-264`: `current_delta = group.inh_delta_error(current_delta,
rng.rng())` under `opts.apply_errors`). `delta_nominal` on refresh is
`(delta_layer(D) + summand) · 0.5` — the `delta_layer` part moves with `D`,
which is the point of F1.7. What happens to the **drawn** part is not
stated anywhere in F1.7:

- If refresh re-runs the branch *with* `apply_errors`, every rate-mode span
  re-draws its δ error at every cycle top — the Monte-Carlo ensemble
  becomes per-cycle instead of per-run, which is a statistics change the
  plan's own discipline (F1.3: "the draw itself is not clamped") exists to
  guard.
- If refresh skips the draw, the drawn contribution must be **carried
  frozen** from construction — which is B2's recipe field.

The plan's own twin collides with this: *"a rate-mode span scaled by hand
from 100 nm to 200 nm, then refreshed, is bitwise equal to the same span
expanded from scratch at 200 nm"* — a from-scratch expansion at 200 nm
**draws a fresh δ error** under `apply_errors`, so bitwise equality can only
hold with draws off (`ExpandOptions::deterministic`) on both sides.

**Amendment.** One decision, stated in F1.7's "What refresh does": **the
draw is frozen at construction and carried in the recipe; refresh recomputes
`delta_layer` and re-applies the frozen draw.** The refresh-correctness twin
gains "(runs deterministic on both sides; the drawn-error case is pinned by
the recipe contract, not by this twin)". This also closes the loop with
F1.3's statistical gate, which pins the mean across the clamp boundary at
construction — per run, not per cycle.

### B4 — F0.2's licence contradicts its own item 6 unless the report is conditional

**Evidence.** The licence opens: *"A run whose stack contains **no**
multi-row span must be bit-identical, full stop."* Item 6 then licenses:
*"The result dict gains a clamp report that was never there"* — and the
change text threads `ClampReport` into `PipelinePhaseResult` at all three
pipeline sites unconditionally. An unconditionally-threaded field changes
the report shape of **every** run, no-span runs included — the licence
violates itself on the same page, and every golden/parity consumer of the
result dict sees a new key at 0.6.34.

**Amendment.** Make the field `Option<ClampReport>` (or emit `None` when
nothing was removed and nothing was capped). No-span runs then serialize
byte-identically, the licence's "full stop" holds literally, and the
deletion-report twin is unchanged (it asserts `Some`). One line in the
change text; the six-item licence stays six.

### B5 — F1.6's scalable set: the text says "graded", the table includes gradients

**Evidence.** The scale-free table at F1.6 lists `GradientMode::FixedSpan`
→ **F1.6**. But the Change section says *"A graded carrier with `optimize =
true`"*, the licence says *"a graded film with `optimize = true`"*, and
§8.12 says *"the graded carrier"*. In this plan's own vocabulary (§D2)
"graded" and "gradient" are distinct engines. Between F1.1 and F1.6 a
gradient carrier with `optimize = true` homogenizes (F1.1's R4 path); F1.6
must say whether it flips to scalable — the table says yes, the change text
never does. An implementer reading "graded" as `inhomogen`-only ships a
FixedSpan row of the table as dead code, and the homogenize-warning removal
(licence item 2) silently does not cover gradients.

**Amendment.** Wherever F1.6 names the scalable set, say: *"any profiled
carrier — `inhomogen` or `gradient` — with `optimize = true`, whose mode is
in the scale-free set of the table."* The licence's first sentence and §8.12
get the same two words. The pinned path (`optimize = false, needle = false`)
is untouched either way.

### B6 — there is a fifth `sub_layer_count` site, and it is the one users can read

**Evidence.** F1.3 says the function "has three other callers —
`expansion.rs:257`, `structure.rs:272` and `architect.rs:565`, the latter
two for row-count prediction. All four must see the same `delta_layer`."
There is a **fourth caller** the list omits: `property_issues` itself
(`layer.rs:227`), where the `inh_delta == 0` note interpolates the count
into its message ("expands to {n} identical sub-layers"). It is the one
site where the count becomes user-visible text — under `InhMode::RateCapped`
a stale count there ships as a diagnostic that contradicts the solver, which
is worse than a wrong prediction in a cold struct.

**Amendment.** F1.3's caller list gains `layer.rs:227` (five sites, "all must
see the same `delta_layer`"), and the sentence notes that this site is a
message, so its twin is a message-equality check, not a row-count check.

### B7 — "last statement" is wrong, and the truth already satisfies F0.3's ordering

**Evidence.** F0.2 says *"evaluator.rs:373 is the last statement of
`optimize_thicknesses_report`"*. It is not: the clamp at `:373` is the last
**mutation**, and the reported merit is evaluated *after* it
(`:375`, `self.evaluate_merit(stack).map(|mf| (mf, Some(res)))`). The same
ordering already exists at the run's final clamp
(`pipeline.rs:271-274`: optimize → clamp → `evaluate_merit` → `final_mf`).

**Why this matters in the good direction.** F0.3's gate — *"final_mf equals a
fresh merit evaluation of the returned stack … the assertion that catches
reporting the pre-clamp number"* — describes an invariant the tree **already
maintains** at both reporting sites. The twin is still worth writing, but it
should be framed for what it is: a pin on an existing ordering, so that
F0.3's implementation cannot silently invert it — not the introduction of an
ordering that does not exist yet. (The one place the plan's concern is
genuinely live is F0.3's `ClampUpFinal` end-of-run clamp-up itself, which
must sit *before* the final merit evaluation for the same reason.)

**Amendment.** F0.2's sentence becomes "…is the last mutation of
`optimize_thicknesses_report`; the reported merit is evaluated after it
(`:375`), and F0.3's twin pins that ordering." F0.3's gate keeps the
assertion and gains the framing.

### B8 — "seven times per macro cycle" vs seven call sites

**Evidence.** F1.7's cost argument says `optimize_thicknesses` "itself runs
seven times per macro cycle". There are seven call **sites**
(`cleanup.rs:87, :121`; `cycle.rs:107, :249`; `inflate.rs:148, :191`;
`pipeline.rs:271`), but the per-cycle invocation count is variable and
larger: `cycle.rs:249` fires per needle insertion and `cleanup.rs:87` fires
per cleanup-trial iteration. F0.2's own wording — "seven call sites, dozens
of times per macro cycle" — is the correct one, and it makes F1.7's 1000×
argument **stronger**, not weaker.

**Amendment.** F1.7's sentence reads "…which itself runs dozens of times per
macro cycle across seven call sites". No gate changes.

### B9 — the tree already computes the bulk ranges the predicate derives (recommendation)

**Evidence.** `expand` keeps an internal `bulk_spans: Vec<(usize, usize)>`
(`expansion.rs:158`), pushes one entry per span (`:313`), consumes it for
the interface slice's owner lookup (`:231`) — and returns only
`(SolverArrays, Vec<Span>)` (`:103`). Every caller discards it at the same
door as the spans (`from_design` at `structure.rs:256`; the architect sites
at `architect.rs:299, 339, 391`). So the quantity F0.1 defines
arithmetically (`end − start − slice == 1`) and F0.2 needs concretely
(`Σ d_r` over bulk rows) is **already computed and already dropped**.

**Recommendation.** F0.1 may carry the bulk range alongside the span
(a field on `Span`, or the parallel vector returned as a third element) —
then F0.2's span bulk thickness and F1.6's `φ_r` fractions read it instead
of re-deriving it, and the arithmetic predicate remains as the
`debug_assert` cross-check. Not a correction: the derivation is correct and
free. But it is the same information at the same door twice, and the plan's
own §1.3 logic ("the information is already computed; F0.1 is plumbing, not
derivation") applies to `bulk_spans` exactly as it does to `spans`.

### B10 — citation drift (two one-liners)

1. F1.1's "Where" cites "beside the `inhomogen` branch at
   [expansion.rs:257](rust/navette/src/structure/expansion.rs:257)" — the
   branch is at `:258`; `:257` is `let sub = layer.sub_layer_count();`,
   which F1.3 cites correctly one section later.
2. F1.3's "sub_layer_count reads `self.inh_delta` directly today
   (layer.rs:127)" — the read is at `:128`; `:127` is the
   `inhomogen && thickness > 0.0` guard.

The only two off-by-ones found in revision 3; recorded because the plan's
citation discipline is otherwise exact and these are cheap to fix now while
the tree is frozen.

---

## 3. Note on the one unverifiable claim

Ground rule 5 (R1) states the needle-run fingerprint "is carried by a
scratchpad harness that lives outside the repo". The in-tree half of that
claim (no bit-exact needle pin in `validation/`) is verified by two
independent searches; the scratchpad half is the author's knowledge and
cannot be checked from the tree. It is accepted on authority — and F0.1's
committed pin makes the question moot the moment it lands.

---

## 4. Adoption

- Add B1–B5 as AMENDED blocks under F0.2 (B1, B4), F1.7 (B2, B3), and F1.6
  (B5); B6–B8 as one-line corrections under F1.3, F0.2/F0.3, and F1.7; B9
  as a recommendation under F0.1; B10 as two character fixes.
- No re-sequencing: the ladder stays `0.6.33 → 0.6.47`. B1 and B4 are
  changes to F0.2's *content* (one extra rule and licence item; one
  `Option`), not to its position; B2/B3 are F1.7 design decisions that fit
  the item as written.
- The §0.4 ledger gains a `B*` block with the same Verdict format; §6 gains
  no new rows (B1 belongs under the existing row-2 interaction; B2/B3 under
  row 10; B4 is a licence detail, not an interaction).
- No version bump; docs-only.

---

## 5. Postscript — status after revision 4

Revision 4 (the +397/−75 rewrite of `docs/implementation_plan.md`) adopted
the findings above and corrected three of them. Each correction was
re-verified against the tree before being accepted here. The record:

| Finding | Status after revision 4 | Evidence |
|---|---|---|
| **B1** | **Confirmed by reproduction on the installed 0.6.32 release build, and worse than this file recorded.** The audit said the interface physics disappears; it also **shortens the film**. Expansion carves the slice out of the carrier (`layer_thickness -= t_interface`, [expansion.rs:227-231](rust/navette/src/structure/expansion.rs:227) — mechanism verified in source, consistent with the measured 150.0 → 149.0; control at 3.0 nm holds 150.0). That falsifies the plan's "the slice is not user thickness and never was": the floor and the cap must compare the slice-inclusive total (`D`), while F1.6 scales only the bulk (`D_bulk`). Two different numbers the plan was about to conflate. | adopted as F0.2 defect 5 with the measured numbers, the `D`/`D_bulk` table, the never-a-candidate rule, and licence item 7 (six → seven everywhere) |
| **B2** | **Confirmed, and the fix shape recommended here was wrong.** A side map keyed by `logical` renumbers under `remove_film` and `insert_needle_seed` — the exact failure F0.1 exists to prevent. The gap is also larger than a profile spec: the stack carries no provider, no wavelength values (only the count), no groups — refresh needs everything one iteration of `expand`'s emission loop reads (the carrier `Layer`, the resolved `Group`, the endpoint spectra on the run's grid). | superseded: the recipe rides `spans` index-for-index (`Vec<Option<SpanRecipe>>`, `None` for plain films), and the missing body is supplied by factoring `expand`'s per-entry emission into `emit_entry`, called by both construction and refresh — open item 15, resolved below |
| **B3** | **Withdrawn.** No error draw ever happens on the synthesis path. The only `DesignStack`-producing paths are `from_design` (expands with `ExpandOptions::deterministic()`, [structure.rs:261](rust/navette/src/smatrix/synthesis/structure.rs:261) — verified) and `with_films` (no `expand` call at all — verified); the two `apply_errors: true` sites (`structure/structure.rs:240`, [architect.rs:344](rust/navette/src/structure/architect.rs:344) — verified) return arrays and never build a stack. The refresh-correctness twin is sound as written, not accidentally sound. The defensive refusal (a recipe records its `ExpandOptions`; refresh refuses `apply_errors: true`) is retained so it stays that way. | B3 dropped; the decision recorded in F1.7 |
| **B6** | **Corrected upward.** A fifth reader exists and is public API: the PyO3 getter at [navette-py/src/structure.rs:360-362](rust/navette-py/src/structure.rs:360) (`Layer.sub_layer_count`, confirmed live: returns 11 for a 100 nm δ=0.1 layer). Five readers, two of them user-visible — the advisory message and the getter. A stale value lands in a user's hands, not in a cold struct. | adopted into F1.3's five-reader table |
| **B4, B5, B7, B8, B9, B10** | Adopted as written. B9 as a decision to carry the bulk range (third return element or a `bulk_start` field; `Span` stays `Copy` either way); B7 as the reframed pin on an ordering the tree already keeps; B4 as `Option<ClampReport>`; B5 as the profiled-carrier wording plus the gradient half of licence item 1 (a gradient has no base index; R4 synthesizes the row from an EMA at f_mid). | — |

Spot-checked revision 4 claims that depended on this file's citations — the
emission loop at [expansion.rs:163](rust/navette/src/structure/expansion.rs:163),
the `oi != k` backwards rescale at [:231-238](rust/navette/src/structure/expansion.rs:231),
the `inv = false` construction at
[structure.rs:255](rust/navette/src/smatrix/synthesis/structure.rs:255), and
the `Span`/`LayerSpec`/`DesignStack` field tables — all resolve.

Open item 15 (the `emit_entry` extraction) was resolved by the audit: see
the RESOLVED block in `docs/implementation_plan.md` §8.15. Summary of the
call: factored, and the extraction lands in F0.1 — gate-match (F0.1 is the
only item where any difference is a defect; F1.7's differential gate could
mask a refactor bug that stays inside licence bounds), one restructure under
one gate (F0.1 already restructures the same loop for the bulk-range carry,
and F1.1/F1.6's later branch additions then edit `emit_entry` instead of a
monolith), and failure mode (the improvised alternative was the
`logical`-keyed map).
