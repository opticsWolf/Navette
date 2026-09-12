# Implementation plan — pre-flight audit: amendments & recommendations

STATUS: audit of `docs/implementation_plan.md` (written at commit `c01e718`)
read against the tree it declares as its base (`72a2d4d`, 0.6.32). The diff
`72a2d4d..HEAD` touches only `docs/`, so every code citation below applies to
the plan's own base commit. This document does not supersede anything; it is
a set of corrections to fold into the plan as CORRECTIONS blocks (A-items)
and non-blocking advice (R-items), in the same ritual the plan itself
prescribes for its predecessors.

Method: every file:line citation in the plan was opened and read; every
"already exists" / "already safe" claim was traced to the code; the gates'
preconditions (fixtures, benches, sync tests) were inventoried. The three
source plans (`implementation_plan.md`, `plans/multi_environment_plan.md`,
`plans/gradient_layers_plan.md`) were read in full.

---

## 0. Verdict

The plan is unusually accurate. Every line-number citation checked lands
exactly (§1 below). Its five corrections of the source plans (§1.1–1.5) all
hold against the tree. The sequencing argument (F0.1 before gradients before
multi-environment), the version ladder, and the gate philosophy survive the
audit unchanged.

What the audit found is different in kind: **omissions and feasibility gaps,
not contradictions.**

| # | Finding | Class | Severity |
|---|---|---|---|
| A1 | The program-schema gate (`config.rs`) is the same equality trap §1.2 diagnosed in `version.rs` — and F2.4 must bump it, which refuses every existing program file | blind spot (same class as §1.2) | high |
| A2 | F0.1's mutator list omits `clamp_all` — a fifth row-count-changing mutator called per-cycle in the hot pipeline | omission | high |
| A3 | F1.4 covers only the Rust half of the state-schema gate; the Python mirror constant, its policy comment, and the cross-language sync test are unaddressed | omission | medium |
| A4 | §1.4's "zero lines of new code" is off by one line: the background predicate tests `l.inhomogen`, which a gradient film never carries | small factual correction | low |
| A5 | A gradient design film's `material_b` has no nk anywhere in the synthesis design path — the decision cannot wait for F1.5 | feasibility gap | medium |
| A6 | Mori-Tanaka is in-tree, Python-bound, and named by §D9.4 — but excluded from the plan's `EmaModel` enum without a stated reason | inconsistency | low |
| A7 | The F2.2 <1% bench gate needs a baseline that no existing benchmark produces | feasibility gap | medium |
| A8 | F2.2/F2.3's "Where" cites `merit.rs` signatures; the actual call sites live in `pipeline.rs` and `evaluator.rs` | precision | low |
| A9 | Minor precision: "four setters", the F1.2 message punctuation, `Span` derives | wording | trivial |

---

## 1. Claims ledger (what was verified and holds)

| Plan claim | Verdict | Evidence |
|---|---|---|
| §1.2: the schema gate is an equality, "v1 is the baseline, there is no past" | ✓ exact | `rust/navette/src/structure/version.rs:8,16-19` (`Some(v) if v == SCHEMA_VERSION`); doc comment verbatim |
| §1.3: `Span { start, end, logical, slice }` exists; `from_design` builds rows from spans, then discards them | ✓ exact | `expansion.rs:31-38`; `structure.rs:256-269` (rows built per span, `Self::with_films(...)` drops the vector) |
| §1.4: background set = `inhomogen && !optimize && !needle` | ✓ exact for inhomogen | `driver.rs:128-134` (one-line caveat → A4) |
| §1.4: background-pinned graded films homogenize/keep-profile behavior | ✓ | `structure.rs:225-236` (homogenize flips `inhomogen=false` with warning; background contains → keeps profile) |
| §1.5: thin-removal filters on `optimize && d_nm < threshold`, comment names pinned graded sublayers | ✓ exact | `cleanup.rs:52-58` |
| §1.5: needle host admissibility is `layer.needle` | ✓ exact | `needle_pass.rs:62-70` (`build_scan_sites`, `if layer.needle && d > 0.0`) |
| F1.1: `layer.rs` is 606 lines; no gradient field, no `inh_mode`, no `gradient.rs` | ✓ | `wc -l` = 606; `Layer` fields `layer.rs:62-87`; `structure/` has no gradient module |
| F1.1: inhomogen branch and `δ_nominal = (inh_delta + inh_delta_summand) · 0.5` | ✓ exact | `expansion.rs:258,259`; summand lives in `group.rs:90` |
| F1.1: the four EMA kernels are in-tree | ✓ | `ema.rs` — `lichtenecker:23`, `looyenga:37`, `maxwell_garnett:70`, `bruggeman:136` (Mori-Tanaka too → A6) |
| F1.3: `sub_layer_count()` is the `thickness^0.4` rule | ✓ | `layer.rs:126-133` (δ-dependent factor present — consistent with the plan's "pure function of thickness, no circularity" note) |
| F1.5: `LayerRow` carries `inhomogen` at the cited line; `ArrayFilm` exists | ✓ | `design_config.rs:85`; `driver.rs:32-45` |
| F2.2: `merit`/`residuals` are single-`&SimCurves`, at the cited lines | ✓ exact | `merit.rs:736`, `merit.rs:755` |
| F2.2: single-env call sites are thin | ✓ (but elsewhere than merit.rs → A8) | `evaluator.rs:284` (`spec.merit(&sim, 1e6)`), `evaluator.rs:333` (LM residual closure) |
| F2.4: `run_needle` at the cited line; `_dump` pass-through; three demand kinds; backside curves | ✓ | `pipeline.py:188`; `spectralweave/target.py:99,157,252`; `SpectralTarget:34 / AngularTarget:138 / ColorTarget:198`; `SimCurves::back` `merit.rs:321-323` |
| Ground rule 1: seven version sites in five files | ✓ | `pyproject.toml:14`; `Cargo.toml:9,29`; `Cargo.lock:465,492`; `__about__.py:12`; `README.md:190` (139/170 are prose history, not bump sites) |
| Ground rule 3: the battery's named artifacts exist | ✓ | 10 review harnesses in `validation/review/`; `tools/check_{exposure,cie_sync,pyi_sync,toolchain}.py`; `opt-minpack-lm`/`opt-argmin` features (`rust/navette/Cargo.toml:24,29`); `build_profile()` (`__init__.py:53`) |
| Ground rule 5 predecessors: randomized differential suite exists | ✓ | `validation/regression/structure/test_differential.py:95,119` (`test_random_stacks_bit_identical`, `test_random_architects_bit_identical`) — anchor naming → R1 |
| F1.4: `test_roundtrip.py` / `test_program.py` exist | ✓ | `validation/regression/structure/`, `validation/regression/config/` (ownership mix-up → A3) |
| F0.1: merge is nk-keyed; distinct sublayer nk never collapse | ✓ | `structure.rs:332-360` (merges on `material == && nk ==`) |
| F1.5: "Python `Layer`, validated natively" | ✓ | `Layer` is the pyo3 class (`rust/navette-py/src/structure.rs:185`), re-exported from `navette.structure.models` |
| Remediation plan closed, every row struck | ✓ | `docs/remediation_plan.md` master table: R1.1–R4.5 all struck with DONE versions |

---

## 2. Amendments (fold into the plan as CORRECTIONS blocks)

### A1 — F2.4 must give the program schema the §1.2 treatment

**Plan section:** F2.4 (0.6.42), cross-plan interaction row 4.

**Evidence.** The program envelope gate is a second equality gate, built the
same way as the state gate §1.2 diagnosed:

- `config.rs:26` — `pub const PROGRAM_SCHEMA_VERSION: u32 = 1;`
- `config.rs:282-284` — `Some(v) if v == PROGRAM_SCHEMA_VERSION as u64 => {}`,
  else `"program schema_version {other:?} unsupported"` — the identical shape
  as `version.rs:16`.
- Unknown top-level keys **refuse**, they are not ignored: `config.rs:293-305`
  collects `unknown` keys and errors; `config.rs:733-740` pins it
  (`"bogus": 1` → err; `schema_version: 2` → err).

**Why it bites.** `design:`/`environments:` are new top-level sections in a
program document. Under the v1 gate they refuse as unknown keys, so F2.4
*necessarily* bumps `PROGRAM_SCHEMA_VERSION` — and the moment it does, every
v1 program file on disk refuses to load. `plans/multi_environment_plan.md`
§4.1's "additive-optional → old calls parse unchanged" is true for the payload
structs and false for the envelope gate. This is §1.2's own failure class,
reproduced one gate over, and the plan's §1.2 correction was written for the
state gate only.

**Amendment.** F2.4 gains the F1.4 treatment, in both homes:

- Split the constant: `PROGRAM_SCHEMA_VERSION = 2` (what this code writes) +
  `MIN_READABLE_PROGRAM_SCHEMA_VERSION = 1`; `gate_document` accepts the
  range, refuses below (stale) and above (newer build) with distinct
  messages. Both homes move together: `config.rs:26` and
  `src/navette/config/program.py:45`.
- Extend `validation/smoke/test_request_bits.py` — it already syncs both
  version pairs (`:12-13`; state pair `:245-257`, program pair `:268-289`);
  the program test must probe both range endpoints.
- The existing refusal test at `config.rs:727-731` (`schema_version: 2`
  refuses) flips meaning after the bump: a v2 program must now be accepted
  and a v3 refused. Update it deliberately, not as drive-by fallout.
- Cross-plan interaction row 4 gains: "the program envelope bump repeats the
  state-gate pattern — same range, same two-sided messages, same sync test."

**Gates to add to F2.4.** A v1 program fixture loads and assembles
bit-identically; a `schema_version: 3` program refuses naming the newer
build; both endpoints of the program range asserted in
`test_request_bits.py`.

### A2 — `clamp_all` is a fifth row-count-changing mutator; F0.1 must cover it

**Plan section:** F0.1 (0.6.33).

**Evidence.** F0.1's change list maintains spans across "insert_needle_seed,
merge_adjacent, remove_film, and remove_thin_layers". There is a fifth:

- `structure.rs:381` — `clamp_all(min_nm, max_nm)`: removes **every** film
  below `min_nm` (no `optimize`/`needle` guard at all — the loop is
  unconditional) and caps the rest. Row count changes.
- It is hot: called after cleanup and after each LM step, per cycle
  (`pipeline.rs:192, 213, 273`) and after thickness optimization
  (`evaluator.rs:373`).
- Default `clamp_min_nm = 2.0` (`config.rs:88`). Gradient sublayers at the
  F1.1 max-step floor stay above 2 nm today, so nothing breaks now — but the
  protection is coincidence, which is exactly the coupling F0.1 exists to
  retire.

**Amendment.**

- F0.1's mutator list gains `clamp_all`: span maintenance + the
  `debug_assert` partition check. (Its cap branch never changes counts — no
  span action there.)
- The non-singleton exemption extends to `clamp_all`'s **removal branch**,
  alongside the existing exemption in `remove_thin_layers`. Without it, the
  first optimizable gradient span (F0.1's own motivating future) is eaten at
  the LM boundary every cycle, not at cleanup.
- F0.1's twin batch gains a clamp twin: an optimizable thin-sublayer span
  under a pipeline config whose `clamp_min_nm` exceeds the sublayer
  thickness survives the clamp.

### A3 — F1.4 must name the Python half of the state-schema gate

**Plan section:** F1.4 (0.6.37).

**Evidence.** F1.4's change list names `version.rs` only. The gate has a
second, non-thin home and two pinned consumers:

- `src/navette/structure/types.py:119` — `SCHEMA_VERSION = 1` is a second
  literal (the Python `check_schema_version` at `:123-127` is thin over the
  native gate — that half is fine; the constant is not).
- The constant's doc comment states a policy the bump contradicts:
  "purely additive keys are safe without a bump (readers ignore unknown
  keys) — the fingerprint test in
  `validation/regression/structure/test_roundtrip.py` enforces this decision
  on every key-set change." F1.4's fields are additive-with-defaults and the
  plan bumps anyway (the newer-writer hazard justifies it) — so the comment
  must be rewritten the same commit, or the code contradicts its own
  documented policy the week it lands.
- `validation/smoke/test_request_bits.py:245-257`
  (`test_state_schema_version_matches_the_native_gate`) asserts the Python
  constant equals what the native gate accepts. A range gate turns "what the
  native accepts" from a point into an interval — the test must probe both
  endpoints, and the Python side needs the second constant
  (`MIN_READABLE_SCHEMA_VERSION`) to do it.

**Also.** F1.4's gate sentence "test_roundtrip.py and test_program.py
extended rather than replaced" mixes schemas: `test_program.py`
(`validation/regression/config/`) exercises the *program* schema. The
state-schema item extends `test_roundtrip.py` only; the `test_program.py`
half moves to F2.4, where it belongs (A1's item).

### A4 — §1.4's "zero new code" claim is off by one line

**Plan section:** §1.4 (What changed), F2.2's Phase-A interaction twin.

**Evidence.** The background predicate at `driver.rs:130` selects
`l.inhomogen && !l.optimize && !l.needle`. A gradient film carries
`gradient: Some(..)` with `inhomogen = false` (the two are mutually
exclusive by §D2 validation), so the predicate never selects a gradient film
until it learns the field. And `assemble_stack`'s per-field copy
(`driver.rs:115-120`) must copy `gradient` into the `Layer` at all, or
`from_design` cannot see it.

**Amendment.** §1.4's sentence becomes: "a graded film inside a fixed
environment segment keeps its full profile with zero new physics code; a
gradient film needs two one-line plumbing touches when the field exists
(F1.1): the background predicate gains the `gradient.is_some()` disjunct,
and `assemble_stack` copies the field." The physics claim (§D-rule) is
unchanged; F2.2's twin still pins it.

### A5 — two nk spectra must reach the design-path provider; decide at F1.1, not F1.5

**Plan section:** F1.1 (0.6.34), F1.5 (0.6.38).

**Evidence.** The structure path is fine: `Layer.gradient` names
`material_a`/`material_b` and a `MaterialProvider` resolves both — §D2's
"same provider" rule covers it. The synthesis design path is not:

- `pipeline.py:91-121` (`_film_dicts`) evaluates **one** material per film
  into `nk` (`_eval_nk(mat, wl)`, `:117`).
- `ArrayFilm` (`driver.rs:32-45`) carries one `nk: Vec<Complex64>`.
- `assemble_stack` registers provider entries as film-name → film-nk only
  (`driver.rs:244-251`); there is no materials library in this path.

A gradient design film's `material_b` therefore has no nk anywhere in the
design path — and F1.1's own scope (expansion branch, homogenize path in
`from_design`) consumes it from the first commit that exercises a gradient
film, not at F1.5.

**Amendment.** Resolve as an open decision at F1.1 (new slot in §8),
recommended resolution:

- `ArrayFilm` gains `gradient: Option<GradientJson>`; the JSON carries
  `material_b` **and** `nk_b` (the Python side evaluates both materials at
  the door, in `_film_dicts`); `assemble_stack` injects
  `material_b → nk_b` into the provider entries when a gradient is present.
  This keeps the path self-contained (its existing convention: every nk is
  film-supplied) and keeps §D2's same-provider refusal as the gate.
- The design film keeps its **own name** (it is the contrast/needle
  identity, multi-env plan §2.3); `material_a`/`material_b` are provider
  keys only, never the film's identity.

F1.5 mirrors the decision in `LayerRow`. Refusal when `nk_b`/`material_b` is
absent — never fall back to `material_a`'s nk (a half-mixture is worse than
a refusal; §D2's own words).

### A6 — Mori-Tanaka: in-tree, bound, and unaccounted for

**Plan section:** F1.1 (0.6.34); `§D2` of the gradient plan (authoritative).

**Evidence.** `§D2`'s `EmaModel` enum — restated in F1.1 — lists Bruggeman,
MaxwellGarnett, Looyenga, Lichtenecker. But `mori_tanaka` is in-tree
(`ema.rs:87`; spec variant with shape factor `l` at
`materials/mod.rs:66`) and Python-bound (`_materials.pyi:56`
`ema_mori_tanaka(n_i, n_h, f, l=1/3)`; `materials/__init__.py:246-248`
handles `"MoriTanaka"`). §D9.4's host-asymmetry warning names Mori-Tanaka in
the same breath as Maxwell-Garnett, so users will expect it here.

**Amendment.** Either choice is safe — a later variant is a new serde enum
name, and an older binary meeting it hits the newer-build refusal (the F1.4
range gate earning its keep). Resolve as a one-liner open decision at F1.1
(§8 already has the slot). Recommended: v1 ships the four, with a stated
exclusion sentence in the `GradientSpec` docstring — "Mori-Tanaka is
in-tree and bound for homogeneous EMA composites; it joins `EmaModel` as a
fifth variant when a real case needs it (its shape factor `l` needs a spec
field first)" — so the omission reads as a decision, not an oversight.

### A7 — the <1% bench gate needs a benchmark that does not exist yet

**Plan section:** F2.1 (0.6.39, baseline recording), F2.2 gate (b).

**Evidence.** Inventory of `validation/benches/`:

- `synthesis/` — `bench_refold.py` only: refold-vs-static comparison with a
  simulate+fold vs LM cost breakdown. The closest existing proxy, but not a
  plain eval loop.
- `smatrix/` — `bench_backside_speed.py`, `bench_core_engine_scaling.py`.
- `structure/` — `bench_grid_assert.py`.

No harness times assemble → simulate → merit of the design pipeline in
isolation. F2.1 says "Baseline recorded here: pre-segment eval timings from
`validation/benches/`" — as things stand there is nothing to record.

**Amendment.** F2.1 adds `validation/benches/synthesis/bench_eval.py` (one
fixed design problem, timed assemble+simulate+merit loop, release-gated via
`_bench_common.require_release` like `bench_refold.py`'s preamble) as part
of the item, and its results become the F2.2 baseline. This is the plan's
own logic ("a baseline taken after the change is not a baseline") applied
one step earlier: the baseline's *harness* must exist before the baseline.

### A8 — name the real change surface of F2.2/F2.3

**Plan section:** F2.2 (0.6.40), F2.3 (0.6.41).

**Evidence.** The cited `merit.rs:736/755` signatures are correctly quoted
and correctly kept. But the single-environment call sites the multi-variant
must route through live elsewhere:

- `evaluator.rs:282-285` — `DesignContext::evaluate_merit` calls
  `self.spec.merit(&sim, 1e6)`; `:333` — the thickness-LM residual closure
  clones **one** base stack, sets thicknesses, simulates, calls
  `spec.residuals(&sim, out)`.
- `pipeline.rs:63-97` — `NeedlePipeline` owns exactly one `DesignStack` and
  one context; `run()` (`:119`) is the macro-loop.

K>1 therefore means: `NeedlePipeline` holds the shared design object plus
per-environment surroundings and rebuilds K stacks at run start and after
each insertion (multi-env §4.6's assembly-time rule); the LM residual
closure becomes span-aware across K stacks; `evaluate_merit` becomes the
routing point for `residuals_multi`. None of it touches `merit.rs` — the
plan's "no merit formula change" contract holds — but sizing the item from
its `merit.rs` citation alone underestimates where the lines go.

**Amendment.** F2.2/F2.3 "Where" sections gain: `pipeline.rs` (NeedlePipeline
fields + run loop), `evaluator.rs` (DesignContext impl + LM residual
closure), `driver.rs` (K assembly + the one K-branch). No re-sequencing;
this is why F2.3 is effort-L. Re-check the multi-env plan's "~350 Rust
lines" estimate at F2.2 against this surface.

### A9 — minor precision (three one-liners)

1. **"Four setters" (F1.1).** The native `Layer` binding exposes twelve
   setters (`rust/navette-py/src/structure.rs:230-355`). The R3.5 numeric
   gate covers the numeric ones — thickness, inh_delta, roughness,
   interface_thickness — which is what "four" means. Reword to "the four
   numeric-property setters" so the count is decodable at review.
2. **F1.2 refusal message punctuation.** F1.2's quoted message is ASCII
   (`"...are the same slope - give one."`) — correct per ground rule 7.
   §D0's source version uses an em-dash. Implement F1.2's ASCII form; do
   not copy §D0's punctuation.
3. **`Span` derives (F0.1).** `Span` already derives `Clone, Copy,
   PartialEq, Eq` (`expansion.rs:29-38`) — storing `Vec<Span>` on
   `DesignStack` needs no derive work. F0.1's data-model note can say so
   and stay focused on the bookkeeping, which is where its own risk
   assessment already correctly puts the item.

---

## 3. Recommendations (non-blocking)

**R1 — name the fingerprint anchors in the plan.** Ground rule 5 makes two
fingerprints hard gates, but only one has a nameable home:
`validation/regression/structure/test_differential.py:95`
(`test_random_stacks_bit_identical`; also `:119` for architects). No single
named test pins a bit-exact needle run — the coverage is spread over the
review harnesses (`fd_step1`, `fd_rchannel`, `weaver_race`) and smoke tests.
Give each item's "fingerprint unmoved" line the exact test(s) that prove it,
and if no needle-run anchor is named by F0.1, add one small explicit pin
(a fixed design run, hashed spectra) rather than leaving the gate
unmechanized.

**R2 — one named partition helper.** Implement F0.1's `debug_assert` as a
single `assert_spans_partition(&self)` (contiguous, ordered, no gaps, no
overlap) called at the tail of every count-changing mutator **and both
constructors**. A named helper survives A2's fifth mutator and any future
one; four independently-written asserts do not.

**R3 — the needle-refusal twin should exercise the mutator directly.**
`insert_needle_seed` (`structure.rs:296`) is public and does no
admissibility check of its own — the `needle` flag gates
`build_scan_sites` (`needle_pass.rs:62-70`). F0.1's refusal should exist at
the mutator as well as the scan-site filter, and the twin should call the
mutator directly, not only through the scan.

**R4 — gradient homogenize cannot reuse the flag flip.** Today's homogenize
mutates a clone (`h.inhomogen = false`, `structure.rs:236`) and re-expands
uniform — reusable because a base nk exists. A gradient film has no base
nk: homogenize must synthesize the row itself (EMA at `f_mid` over both
provider spectra, which A5 guarantees are present). One sentence in F1.1's
"Where" so nobody reaches for the existing flip.

**R5 — commit the v1 state fixture before touching `version.rs`.** The plan
requires the fixture; sequence it as F1.4's first commit so the range
gate's acceptance test has its oracle from the first line of the change,
and the refusal-side test (v3 naming the newer build) is written against
the un-bumped gate before it changes.

**R6 — record the F2.1 baseline in-repo.** A results-file convention
already exists (`validation/benches/smatrix/results/`). Record F2.1's
baseline under `validation/benches/synthesis/results/` with the commit hash
and date, not only in the item text — the F2.2 gate then measures against
an artifact, not a prose number.

**R7 — progress log records adopted amendment IDs.** When an item lands,
its CORRECTIONS block lists which A/R numbers it adopted. Keeps the plan
file, this audit, and reality reconcilable at review time.

---

## 4. Risk-table additions (plan §7)

| Risk | Mitigation |
|---|---|
| F2.4's program bump refuses every existing v1 program file | A1 range gate in both homes; `test_request_bits` endpoint probes; v1 fixture |
| `clamp_all` deletes a future optimizable gradient span after every LM step | A2 exemption in the removal branch + the clamp twin |
| A gradient design film's `material_b` has no nk on the design path | A5 option 1 (nk_b at the door); absent-material refusal, never a fallback |
| F2.2's <1% gate measures against a baseline that does not exist | A7 `bench_eval.py` added in F2.1; R6 records it as an artifact |
| The plan file's own claims drift from reality mid-series | R7 amendment IDs in the progress log |

---

## 5. Adoption

- Add A1–A9 as CORRECTIONS blocks under the affected sections of
  `docs/implementation_plan.md` (F0.1, F1.1, F1.4, F1.5, F2.1, F2.2, F2.3,
  F2.4, §1.4, §7, §8). Strike nothing; the source plans stay untouched.
- No version bump: this is a docs-only change; the ladder starts at F0.1
  (0.6.33) as planned.
- The two open decisions this audit adds (A5's nk routing, A6's
  Mori-Tanaka) slot into the plan's §8 with their recommended resolutions,
  to be closed at F1.1.
