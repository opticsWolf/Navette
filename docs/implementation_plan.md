# Navette Feature Implementation Plan

Derived from `docs/plans/gradient_layers_plan.md` and
`docs/plans/multi_environment_plan.md`, both written against the 0.4.x tree
and never started. This document supersedes the phasing, version reservations
and progression sections of both; the physics, schema and rejected-alternative
sections of the source plans remain authoritative and are cited by section
(`§D*` = gradients, `§*` = multi-environment).

Branch: `dev_feature`. Base: `72a2d4d` (0.6.32, `main` == `dev_remedy`).

Status of the predecessor: `docs/remediation_plan.md` is **closed** — every
row of its master table is struck. This plan is the successor workload and
follows the same ritual.

**Revision 2 (pre-flight audit folded in).** `docs/implementation_plan_amendment.md`
audited revision 1 against the tree; that audit was then itself re-checked
against the tree, line by line. Both passes are folded in here as **AMENDED**
blocks under the items they touch, with an `A*` (from the amendment) or `N*`
(found while re-checking it) tag. Three of the amendment's findings were
wrong or under-scoped and are corrected rather than adopted verbatim (A1, A2,
A6, A9.1). The amendment file stays as the record of the first pass; this
document is the one to work from. §0.4 is the ledger.

**Revision 3 (two requested features, and one binding constraint).** Four
directed requirements arrived after revision 2 and are folded in here with a
`U*` tag: a **clamp-up-to-minimum** option beside today's remove-only floor
(U1 → F0.3); **scalable graded spans** (U2 → F1.6); a **recalculation loop**
for the rate-type profiles whose index is a function of absolute thickness
(U3 → F1.7); and the constraint that binds the last two — **re-slicing a
graded span runs at construction only, never inside an LM round or any other
inner loop** (U4 → F1.7, and it turns out to buy more than it costs).
Revision 2's F0.1 also splits: the bookkeeping (F0.1) is separated from the
behaviour changes it enables (F0.2), so exactly one release in the series is
allowed to move a number. The ladder is now `0.6.33 → 0.6.47`.

---

## 0. How to use this plan

### 0.1 Ground rules

1. **One item = one feature = one `0.0.1` version bump.** Seven version
   sites in five files move together
   (`pyproject.toml:14`, `Cargo.toml:9,29`, `src/navette/__about__.py:12`,
   `Cargo.lock:465,492`, `README.md:190`). Then: CHANGELOG entry, item
   marked **DONE** here with a CORRECTIONS block if reality differed from
   the design, master-table row struck, commit, push to `dev_feature`.

   **N6 — do not `sed` the README by file.** `README.md` carries three
   `0.6.32` strings and only one is a bump site. `:190` is the wheel
   filename and moves. `:139` ("the tree-wide reformat commit (0.6.32)")
   and `:170` ("`cargo fmt --all --check` (since 0.6.32)") are historical
   facts and must not move; a blanket
   `sed -i 's/0\.6\.32/0.6.33/g' README.md` silently rewrites both and
   makes the file lie about its own history. Bump the line, not the file.
2. **Release build only.** `maturin develop --release`;
   `navette.build_profile()` must print `release` before any timing claim.
3. **Verification battery** (all green before a commit):
   - `python -m pytest validation -q`
   - `cargo test --workspace`
   - `cargo clippy --workspace --all-targets -- -D warnings`
   - `cargo fmt --all --check`
   - `python tools/check_exposure.py`, `check_cie_sync.py`,
     `check_pyi_sync.py`, `check_toolchain.py`
   - the ten `validation/review/*.py` harnesses, each exit 0
   - feature-gated: `cargo test -p navette --features opt-minpack-lm`,
     `--features opt-argmin`, and both
4. **Toolchain skew is a known trap.** `tools/check_toolchain.py` exists
   because clippy reports only the lints its own version knows; CI was red
   for 17 consecutive pushes while every local run was honestly green. Run
   the gates against the newest installed toolchain
   (`cargo +<newest> clippy ...`) before claiming green, and wait for CI
   before fast-forwarding `main`.
5. **Bit-exactness is a contract.** Two fingerprints must not move for any
   item that is not explicitly a physics change:
   - random-stack solver fingerprint
   - needle-run fingerprint

   Every item below names which of the two it is allowed to move. **No item
   in this plan is allowed to move either one on a stack that uses no new
   feature.** That is the whole safety argument: both features are opt-in by
   a field that defaults to absent.

   **R1 — only one of the two has an in-tree anchor.** The random-stack
   fingerprint is
   `validation/regression/structure/test_differential.py:95`
   (`test_random_stacks_bit_identical`, with `:119` for architects). There
   is **no** in-tree test that pins a bit-exact needle run; the coverage is
   spread over `validation/review/{fd_step1,fd_rchannel,weaver_race}.py`
   and the smoke tests, and the fingerprint itself is carried by a
   scratchpad harness that lives outside the repo. A gate that only exists
   in a scratchpad is not a gate. **F0.1 commits the needle-run pin**
   (fixed design, fixed seeds, hashed spectra) as its first commit, before
   touching `DesignStack`, and every "fingerprint unmoved" line below names
   the test that proves it.
6. **Behaviour changes are announced, never silent.** Corrections at
   assembly (`from_design`) are warnings with the offending name in them;
   contradictions are refusals naming both sides and the alternative.
7. **Python-surface error messages stay ASCII** (cp1252 consoles).
   Clippy `-D warnings` rejects `///` doc comments on statements.

   *Reality check:* the tree does not currently keep this — the homogenize
   warning at
   [structure.rs:230](rust/navette/src/smatrix/synthesis/structure.rs:230)
   contains an em-dash, which is cp1252-safe but not ASCII. New messages
   are ASCII; do not "fix" the existing one as drive-by, and do not treat
   its punctuation as precedent (see A9.2 under F1.2).

### 0.2 Priority / risk / effort

| Field | Meaning |
|---|---|
| **Priority** | P0 = blocks a later item; P1 = user-visible feature; P2 = plumbing the features need; P3 = docs |
| **Risk** | probability the change itself breaks something already shipped |
| **Effort** | S < 1 h · M = hours · L = 1–2 days · XL = multi-day |

### 0.3 Master item table

| ID | Item | Version | Priority | Risk | Effort | Source |
|---|---|---|---|---|---|---|
| F0.1 | Span provenance on `DesignStack` — the bookkeeping, and no behaviour change | 0.6.33 | **P0** | M | **L** | §D4.1–2, corrected §2 |
| F0.2 | Span-level pipeline accounting — floor, cap, layer budget, inflate, reported counts | 0.6.34 | **P0** | **L** (the only item licensed to move a number) | M | **U5**, A2, N3 |
| F0.3 | `ThinLayerPolicy` — clamp up to the minimum instead of removing | 0.6.35 | P1 | M (the LM lower bound couples to it) | M | **U1** |
| F1.1 | Gradient data model + `FixedSpan` expansion + homogenize path | 0.6.36 | P1 | M (new expansion branch) | L | §D2–D3, §D4.3 |
| F1.2 | Gradient `RateCapped` mode — thickness-relative slope with caps | 0.6.37 | P1 | M (saturation meets `merge_adjacent`) | M | §D0(b), §D2 |
| F1.3 | `InhMode::RateCapped` — thickness-relative single-material drift | 0.6.38 | P1 | M (legacy path must stay bitwise) | M | §D0(b), §D2 |
| F1.6 | One thickness parameter per graded span — the scale-free profiles | 0.6.39 | P1 | M (the LM parameter list stops being a row list) | L | **U2** |
| F1.7 | Profile refresh for the rate modes — at construction points only | 0.6.40 | P1 | M (a refresh in the wrong place costs 1000×) | M | **U3, U4** |
| F1.4 | Schema v2 + a readable-version **range**, not a point | 0.6.41 | **P0** | M (every state file reads through this gate) | M | §D5 + correction §1.2 |
| F1.5 | `design_config` rows + Python `Layer.gradient` surface | 0.6.42 | P1 | S | M | §D5 |
| F2.1 | Environment segment schema + compile/validation + `bench_eval.py` | 0.6.43 | P2 | S | M | §4.1–4.2 |
| F2.2 | K assemblies, K solves, `residuals_multi` — joint merit | 0.6.44 | P1 | M (driver loop) | L | §4.3, §4.6 |
| F2.3 | Needle + LM joint: locus translation, name-routed fold sum | 0.6.45 | P1 | **L** (the hard one — fold routing) | L | §4.4 |
| F2.4 | Python `environments=` / `design=` surface + program sections + program schema range | 0.6.46 | P1 | **M** (second schema gate) | **L** | §4.5 |
| F3.1 | Docs, worked examples, exposure re-audit, release | 0.6.47 | P3 | S | M | §7-S5, §D7 |

Fifteen items, `0.6.33 → 0.6.47`. Tag `v0.7.0` at F3.1 if the minor marker is
wanted; the version number itself stays on the `0.0.1` ladder throughout.

**Read the Version column, not the ID.** F1.6 and F1.7 are new and sit
*between* F1.3 and F1.4 in the ladder. The ID order and the build order
disagree on purpose: the IDs were assigned when the items were written,
and renumbering eight existing items to restore the coincidence would
invalidate every cross-reference in this document to buy nothing.

Why the two new Phase A items land **before** F1.4: F1.4 is the single schema
bump of the phase, and it has to capture the complete Phase A field set in one
release. Scalable spans add fields (§8.12). Bumping the schema before those
fields exist would mean two bumps in one series, which §6 row 4 forbids.

Rows that moved in revision 3: F0.1 splits into F0.1 + F0.2 (risk L→M for the
bookkeeping half — all of the L moved to F0.2); F0.3, F1.6 and F1.7 are new;
F1.4 through F3.1 shift four rungs up the ladder.

Four rows moved in revision 2: F0.1 risk/effort M→L (N1, N3 and the live
`clamp_all` bug), F1.2 risk S→M (N2), F1.4 effort S–M→M (A3, N4), F2.4
risk S→M and effort M→L (A1).

### 0.4 Audit ledger

Every finding, where it landed, and whether it survived re-checking.
`A*` = from `docs/implementation_plan_amendment.md`. `N*` = found while
re-checking that amendment against the tree. `U*` = a directed requirement
from the 2026-09-12 review — not a finding to be verified but a requirement
to be met, so the Verdict column records what checking it against the tree
*revealed*, not whether it is adopted.

| ID | Finding | Verdict | Lands in |
|---|---|---|---|
| A1 | Program schema gate is a second equality gate; F2.4 must give it the §1.2 treatment | **conclusion right, mechanism wrong** — new sections are silently *ignored*, not refused (N9) | F2.4, §6 row 4 |
| A2 | `clamp_all` is a fifth row-count mutator F0.1 omits | **right, and worse than stated** — the hazard is live today, not future (N-arith) | F0.1 |
| A3 | F1.4 covers only the Rust half of the state gate | **right, and understated** — two tests invert, not one, and the range breaks a shared helper | F1.4 |
| A4 | §1.4's "zero new code" is off by one line | verified exact | §1.4, F1.1, F2.2 |
| A5 | A gradient design film's `material_b` has no nk on the design path | verified exact | F1.1, F1.5, §8.8 |
| A6 | Mori-Tanaka is in-tree but excluded from `EmaModel` | **right, under-scoped** — six kernels, and `MixRule` already exists | F1.1, §8.9 |
| A7 | The <1% bench gate has no baseline harness | verified exact | F2.1 |
| A8 | F2.2/F2.3's real surface is `pipeline.rs`/`evaluator.rs`, not `merit.rs` | verified exact | F2.2, F2.3 |
| A9.1 | "four setters" is imprecise | **right, but the corrected count is five** | F1.1 |
| A9.2 | F1.2's message punctuation is ASCII; §D0's is not | verified | F1.2, ground rule 7 |
| A9.3 | `Span` already derives what it needs | verified exact | F0.1 |
| R1 | Name the fingerprint anchors; the needle one has none | verified | ground rule 5, F0.1 |
| R2 | One named `assert_spans_partition` helper | adopted | F0.1 |
| R3 | The needle refusal belongs at the mutator too | verified — `insert_needle_seed` checks nothing | F0.1 |
| R4 | Gradient homogenize cannot reuse the `inhomogen = false` flip | verified | F1.1 |
| R5 | Commit the v1 fixture before touching `version.rs` | adopted | F1.4 |
| R6 | Record the F2.1 baseline as an in-repo artifact | adopted | F2.1 |
| R7 | Progress log records adopted amendment IDs | adopted | §9 |
| **N1** | A span **includes its interface slice row** — F0.1's "non-singleton" rule as written moves the needle fingerprint | new, high | F0.1 |
| **N2** | `merge_adjacent` collapses a saturated `RateCapped` tail; F0.1's "tests unchanged" claim fails at F1.2 | new, high | F0.1, F1.2 |
| **N3** | `inflate_design` re-thicknesses **every** film with no flag guard | new, medium | F0.1 |
| **N4** | `Layer`'s serde is hand-written — F1.4's "serde defaults" mechanism does not exist | new, medium | F1.4, §8.10 |
| **N5** | Provider-existence refusals cannot live in `property_issues`; and the "Where" names the wrong file | new, medium | F1.1 |
| **N6** | Ground rule 1's blanket README `sed` rewrites two historical facts | new, low | ground rule 1 |
| **N7** | Gradient and legacy sublayer thickness conventions differ on purpose | new, low | F1.1 |
| **N8** | `run_needle`'s `layers` is a required positional — F2.4 is a signature change | new, low | F2.4 |
| **N9** | Unknown *section* names are ignored; unknown *top-level* keys refuse; `LayerRow` refuses | new, medium | F2.4, §1.2 |
| **U1** | Removing a sub-minimum layer and setting it *to* the minimum are different operations | **only removal exists** — the one clamp-up in the tree ([inflate.rs:186](rust/navette/src/smatrix/synthesis/inflate.rs:186)) quantizes a different quantity and has no pipeline caller | F0.3 |
| **U2** | A graded span must be scalable as one physical layer | reachable; the whole-layer scaling precedent already exists (`Group::thick_factor`, [expansion.rs:128](rust/navette/src/structure/expansion.rs:128)) | F1.6 |
| **U3** | Only the rate-type profiles need a recalculation after a scale | **verified in-tree** — `InhMode::Fixed`'s factors are a function of δ and `i/(sub-1)` only ([expansion.rs:265](rust/navette/src/structure/expansion.rs:265)); thickness enters solely through `step_t` | F1.6, F1.7 |
| **U4** | The recalculation runs at **construction** only, never in an LM round | binding, and it is also a *correctness* requirement, not only a cost one (F1.7) | F1.7, §6 row 10 |
| **U5** | A graded layer is one physical layer; per-layer rules must never see its solver rows | four rules see them today: floor, cap, layer budget, inflate | F0.1, F0.2 |

---

## 1. What changed since the source plans were written

Both plans were written on 2026-09-06/07 against a 0.4.x tree. Read against
the tree at `72a2d4d`, five of their statements no longer hold: three are
wrong or unsafe (§1.1, §1.2, §1.5) and two describe work that is already
partly done (§1.3, §1.4). Corrected here, not in the source documents, which
stay as the record of the original reasoning.

### 1.1 Version reservations are dead

Gradients reserved 0.4.28–0.4.31; multi-environment reserved 0.4.32–0.4.35.
Every one of those numbers shipped as something else. Renumbered above.

### 1.2 The schema bump as written would refuse every existing state file

`§D5` says "`SCHEMA_VERSION 1 → 2` with default-`None` migration ... every v1
file parses unchanged and reproduces bit-identical results with zero user
action". That is **false against `rust/navette/src/structure/version.rs`**:

```rust
Some(v) if v == SCHEMA_VERSION => Ok(()),
Some(v) => Err(...refusing a stale state...)
```

The gate is an equality, deliberately — its own doc comment says "v1 is the
baseline, there is no past". Bumping the constant to 2 would make every v1
state file on disk refuse to load, which is the opposite of the plan's intent.

Serde defaults alone are not the answer either: a v2 state carrying a gradient,
fed to a v1-reading binary, would drop the gradient and solve a different
stack in silence. That is exactly the failure class this codebase refuses.

**AMENDED (N9) — "in silence" is true of exactly one of the three doors, and
the plan should say which.** The tree has three unknown-key policies and the
argument only applies to the first:

- **State files** (`Layer::deserialize`,
  [layer.rs:357](rust/navette/src/structure/layer.rs:357)): hand-written, doc
  comment "unknown keys ignored" — **silent drop**. This is the door §1.2 is
  about, and the argument holds exactly.
- **`Layer::set_properties`**
  ([layer.rs:252](rust/navette/src/structure/layer.rs:252)): unknown keys
  become returned **warnings**. Loud, not silent.
- **Design-config requests** (`LayerRow`,
  [design_config.rs:74](rust/navette/src/smatrix/synthesis/design_config.rs:74)):
  `#[serde(deny_unknown_fields)]` — **refuses**. Loud, and already correct.

So the range gate is bought for the state path, and the other two doors need
nothing. Keep the resolution; drop the implication that all three leak.

**Resolution (F1.4, promoted to P0 and sequenced before the Python surface):**
split the single constant into a readable **range**.

- `SCHEMA_VERSION = 2` — what this code writes.
- `MIN_READABLE_SCHEMA_VERSION = 1` — the oldest it accepts.
- `check_schema_version` accepts `MIN_READABLE ..= SCHEMA_VERSION`, refuses
  below (stale) and above (written by a newer build) with distinct messages.
  The above-range message is the important one: it is what stops a v1 binary
  from silently dropping a gradient.
- Fields added in v2 reproduce v1 behaviour exactly when absent
  (`gradient: None`, `inh_mode: Fixed`, `shape: Linear`) — but **not by
  `#[serde(default)]`**; see N4 under F1.4.
- Pinned by a test that loads a v1 state fixture committed to the tree and
  asserts bit-identical expansion, plus a test that a `schema_version: 3`
  state refuses.

### 1.3 `Span` exists and is already returned — but is dropped at the door

`§D4.1` says "`Span` already exists in `expansion.rs`" — correct.
`crate::structure::expand` returns `(SolverArrays, Vec<Span>)` and
`Span { start, end, logical, slice }` is exactly the provenance needed.

What the plan does not say is where it dies:
[structure.rs:264](rust/navette/src/smatrix/synthesis/structure.rs:264) —
`from_design` iterates `spans` to build rows, then throws the vector away.
`DesignStack` keeps only a flat `Vec<LayerSpec>`. So the information is
already computed; F0.1 is plumbing, not derivation.

**AMENDED (N1) — but not as cheap as revision 1 claimed.** "Front-loaded
because it is cheap" was the wrong reason to sequence it first. It is
front-loaded because everything else depends on it; it is not cheap, because
a span is not what the plan assumed (see F0.1). Effort M→L.

**AMENDED (revision 3) — and it is two items, not one.** F0.1 is the plumbing:
the span vector survives `from_design`, the mutators maintain it, one assert
proves the partition, one predicate answers "is this span one physical
layer". Nothing observable changes. F0.2 is what the plumbing buys: the
floor, the cap, the layer budget and the reported counts all become span
quantities, and the live `clamp_all` bug is fixed. The split is what lets
F0.1 keep an *unconditional* "both fingerprints unmoved" gate — a single
item could not, because part of it is required to change behaviour, and a
gate with an allowed-delta list is a much weaker gate.

### 1.4 The `optimize/needle = false` rule already routes gradients correctly

`§4.2` (multi-env) forces `optimize = false, needle = false` on every
fixed-segment film. `assemble_stack`
([driver.rs:130](rust/navette/src/smatrix/synthesis/driver.rs:130)) computes
the background set as `inhomogen && !optimize && !needle`, and `from_design`
expands background-pinned graded films **with** their profile.

Consequence: **a graded film inside a fixed environment segment is
automatically background-expanded and keeps its full profile.** A graded
cover glass in a surrounding stack does the physically right thing without
new physics code. Stated here so nobody adds a special case for it, and
pinned by a twin in F2.2.

**AMENDED (A4) — "zero lines of new code" is off by two.** The claim is true
of the *physics* and false of the *plumbing*, because the predicate at
[driver.rs:130](rust/navette/src/smatrix/synthesis/driver.rs:130) tests
`l.inhomogen`, and a gradient film carries `gradient: Some(..)` with
`inhomogen = false` (the two are mutually exclusive by §D2 validation). Two
one-line touches are required at F1.1, when the field exists:

1. the background predicate gains the `|| l.gradient.is_some()` disjunct;
2. `assemble_stack`'s per-field copy
   ([driver.rs:115–120](rust/navette/src/smatrix/synthesis/driver.rs:115))
   copies `gradient` into the `Layer` at all — it is a field-by-field copy,
   not a struct move, so an uncopied field is simply invisible downstream.

The §D-rule itself is unchanged and F2.2's twin still pins it.

### 1.5 The two "unprotected span" holes in §D4.1–2 are already closed

The gradient plan presents thin-removal and needle insertion as live hazards
for profile spans. Against this tree they are not: `remove_thin_layers`
filters candidates on `l.optimize`, and needle host admissibility is
`layer.needle`. Background-pinned spans carry both as `false`, so both
operations already skip them.

**AMENDED (A2, N3) — closed for those two, wide open for three others.** The
audit's re-check found that the pipeline mutates the stack through five more
entry points, and only the two §D4.1–2 names are flag-guarded:

| Operation | Guard | Verdict |
|---|---|---|
| `remove_thin_layers` ([cleanup.rs:60](rust/navette/src/smatrix/synthesis/cleanup.rs:60)) | `l.optimize` | safe |
| needle host selection ([needle_pass.rs:73](rust/navette/src/smatrix/synthesis/needle_pass.rs:73)) | `layer.needle` | safe |
| `insert_needle_seed` ([structure.rs:296](rust/navette/src/smatrix/synthesis/structure.rs:296)) | **none** — public, checks nothing (R3) | reachable |
| `merge_adjacent` ([structure.rs:332](rust/navette/src/smatrix/synthesis/structure.rs:332)) | nk-keyed only | collapses saturated tails (N2) |
| `clamp_all` ([structure.rs:381](rust/navette/src/smatrix/synthesis/structure.rs:381)) | **none** — unconditional removal | **live bug** (A2) |
| `inflate_design` ([inflate.rs:135](rust/navette/src/smatrix/synthesis/inflate.rs:135)) | **none** — re-thicknesses every film | distorts profiles (N3) |
| `round_to_qwot` ([inflate.rs:182](rust/navette/src/smatrix/synthesis/inflate.rs:182)) | **none** | not pipeline-reachable today |

So §D4.1–2 was right about the danger and wrong about *which two operations*
carry it. The plan is not wrong that the protection is a coincidence — it is
wrong that the coincidence currently holds. F0.1 is re-scoped accordingly.

---

## 2. Phase 0 — foundations

### F0.1 — span provenance on `DesignStack`, the bookkeeping half (0.6.33)

**Revision 1 called this "not a bug today — a coupling." Revision 2 downgrades
that: there is a bug today.** The coupling argument still stands and is still
the reason the item exists; it is no longer the only reason.

#### Scope, after the revision-3 split (U5)

F0.1 is the **bookkeeping half and nothing else**: `DesignStack` carries its
spans, the row-count mutators maintain them, one named assert proves the
partition, and one predicate answers "is this span one physical layer".
Every observable of every existing run is bit-identical afterwards. That is
the gate, and it is only reachable because the behaviour changes were lifted
out into **F0.2** (0.6.34) — the item that is *licensed* to change what a run
produces, with an explicit list of what it may change.

The analysis that follows stays here, because it is the reason both items
exist. Where a paragraph describes a behaviour change it is marked
**→ F0.2**, and that is where it is implemented and gated.

#### The live bug (A2, with corrected arithmetic) — fixed at → F0.2

`clamp_all` ([structure.rs:381](rust/navette/src/smatrix/synthesis/structure.rs:381))
removes **every** film below `min_nm` with no `optimize`/`needle` guard — the
loop is unconditional. It is hot: pipeline.rs
[:192](rust/navette/src/smatrix/synthesis/pipeline.rs:192),
[:213](rust/navette/src/smatrix/synthesis/pipeline.rs:213),
[:273](rust/navette/src/smatrix/synthesis/pipeline.rs:273) and
[evaluator.rs:373](rust/navette/src/smatrix/synthesis/evaluator.rs:373) — after
cleanup, after inflate, after every thickness optimization, and once more at
the end of the run. Default `clamp_min_nm = 2.0`
([config.rs:88](rust/navette/src/smatrix/synthesis/config.rs:88)).

The amendment stated that gradient sublayers "stay above 2 nm today, so
nothing breaks now." **That is wrong in both tenses.** Run the real
`sub_layer_count` rule
([layer.rs:126](rust/navette/src/structure/layer.rs:126),
`int(ceil(t^0.4) * (1 + (δ/0.1)·0.5)) + 1`):

| Film | δ | Sublayers | Each | vs 2.0 nm floor |
|---|---|---|---|---|
| 3 nm graded | 0.1 | 4 | 0.750 nm | **all deleted** |
| 5 nm graded | 0.1 | 4 | 1.250 nm | **all deleted** |
| 6 nm graded | 0.2 | 7 | 0.857 nm | **all deleted** |
| 10 nm graded | 0.2 | 7 | 1.429 nm | **all deleted** |
| 10 nm graded | 0.1 | 5 | 2.000 nm | survives |
| 20 nm graded | 0.2 | 9 | 2.222 nm | survives |

A background-pinned graded film thinner than ~10 nm (δ = 0.1) or ~20 nm
(δ = 0.2) is **erased entirely** by the first `clamp_all` of the run — not
degraded, not warned about: gone, and the stack solves without it. That is
already true at `72a2d4d` for the `inhomogen` feature that ships today.

F1.1 makes it worse, not better: the `clamp(…, 3, 64)` floor in the new
sublayer rule means any gradient film thinner than **6 nm** has sublayers
under 2 nm (3 nm → 1.00 nm, 5 nm → 1.67 nm). A 5 nm graded interlayer is an
entirely ordinary thing to author.

#### The coupling (unchanged, still the structural reason)

A span is protected only as a side effect of being *entirely frozen*, and
even that only covers two of the seven mutators (§1.5). Two things break the
argument outright:

- **The obvious next feature.** A gradient whose *total* thickness is a free
  parameter — the first thing anyone asks for with a rugate filter — needs
  `optimize = true` on its rows. The moment that exists, cleanup eats the
  sublayers one at a time, each removal looking locally like an improvement.
- **F2.3.** Multi-environment fold routing maps a shared design parameter to
  its position in each of K assemblies. With profiles, a parameter owns a
  *span*, not a row. A flat `Vec<LayerSpec>` cannot express that, and nothing
  in the current structure can reconstruct it after `from_design` returns.

#### AMENDED (N1) — a span is not what revision 1 assumed

**This is the finding that would have broken the item.** Revision 1's rule was
"non-singleton span ⇒ protected", justified by "no stack without a profile has
a non-singleton span." Read
[expansion.rs:209](rust/navette/src/structure/expansion.rs:209) and
[expansion.rs:307](rust/navette/src/structure/expansion.rs:307): `start` is
captured **before** the interface-slice row is pushed, and the span's `end` is
after the bulk. So **one span contains the interface slice row *and* the bulk
row(s)**, and `slice: emitted_slice` records which.

Consequence: an ordinary film with `interface = true` and no grading at all
is a **two-row span**. Under revision 1's rule as written:

- needle host selection would refuse to host in it — but today the bulk row
  carries `needle = carrier.needle && !is_slice`, and `is_slice` is only true
  for `r == span.start`, so the bulk row **is** an admissible host today;
- `remove_thin_layers` would exempt its bulk row, which is removable today.

Both are behaviour changes on stacks that use no new feature. **F0.1's own
gate — "both fingerprints unmoved" — would fail, and it would fail for the
right reason.**

**The predicate is bulk-row count, not row count.** A span is *singleton* iff

```
span.end - span.start - (span.slice as usize) == 1
```

which leaves interface-carrying films singleton (unchanged behaviour) and
makes graded and gradient spans non-singleton (the intended protection).
`Span` needs no new field — `slice` already carries what is required. Say
"singleton-bulk", never "singleton", everywhere in this item and in F2.3.

#### AMENDED (N2) — `merge_adjacent`'s nk key stops holding at F1.2

Revision 1: "Existing `merge_adjacent` tests unchanged (merge is nk-keyed;
distinct sublayer nk never collapse)." True for `FixedSpan` linear, **false
for `RateCapped`**, and F1.2's own gate says so: its saturation proof is that
the tail rows are *bitwise equal to pure `material_b`* — which makes them
bitwise equal to **each other**. Every sublayer of one carrier shares the
carrier's material name
([structure.rs:270](rust/navette/src/smatrix/synthesis/structure.rs:270)), and
`merge_adjacent` keys on `material == && nk ==`
([structure.rs:346-347](rust/navette/src/smatrix/synthesis/structure.rs:347)). So
the saturated tail merges, and `cleanup_design` calls `merge_adjacent`
unconditionally, twice per cleanup
([cleanup.rs:110](rust/navette/src/smatrix/synthesis/cleanup.rs:110),
[:116](rust/navette/src/smatrix/synthesis/cleanup.rs:116)).

The merge is optically a no-op (identical nk, summed thickness), so this is
not a physics hazard — it is a **bookkeeping** hazard, and bookkeeping is the
entire content of this item. F0.1's span maintenance must handle:

- **intra-span merges** — the span shrinks, stays contiguous, keeps `logical`;
- **cross-span merges** — two adjacent carriers with the same material name
  and bit-identical nk. Rare but reachable. Decide and state it: the merged
  row joins the **first** span (matching `merge_adjacent`'s existing
  "first layer's properties win" rule), and the second span is dropped.

The same reasoning kills a degenerate case: a `FixedSpan` with
`f_start == f_end` produces N bit-identical rows that merge to one, silently
turning a "gradient" into a homogeneous film. Refuse it at validation in F1.1
(naming both endpoints), rather than letting the merge quietly be right.

#### AMENDED (N3) — `inflate_design` is a thickness mutator with no guard

[inflate.rs:135](rust/navette/src/smatrix/synthesis/inflate.rs:135) selects
`(0..n_films).collect()` when `max_layers` is `None` and calls `set_thickness`
on each. Row count is unchanged, so the span partition survives — but the
*profile* does not: every sublayer grows by its own `addon_qwot · λ₀/(4n)`,
which is a different amount per sublayer because the nk differ. A pinned
profile comes out of an inflate pass with a distorted thickness distribution.
It is `enable_inflate`-gated, so it is not on the default path, and interface
slices are distorted the same way today.

`round_to_qwot` ([inflate.rs:182](rust/navette/src/smatrix/synthesis/inflate.rs:182))
has the identical shape and is not reachable from the pipeline.

**Decision for this item:** in scope, minimally. `inflate_design` gains the
same singleton-bulk exemption on its *selection* list — a non-singleton-bulk
span is not an inflate candidate. That is one filter on one `collect()`, it
cannot move either fingerprint (no existing stack has a non-singleton-bulk
span outside the already-pinned graded case, whose inflation is a bug), and
leaving it out would mean F0.1 protects a span from four operations and hands
it to a fifth. `round_to_qwot` gets the same filter for symmetry and a note
that it has no live caller.

#### Change

- `DesignStack` gains `spans: Vec<Span>` alongside `films`, populated in
  `from_design` from the vector it already receives and currently discards.
  Non-`from_design` constructors (`with_films`, used by tests and by the
  PyO3 `DesignStack.__init__`,
  [_smatrix.pyi:428](src/navette/_smatrix.pyi:428)) synthesise singleton
  spans — one row, one span, `slice: false` — so the invariant "every row
  belongs to exactly one span" holds everywhere.
- **Five** mutations that change row count maintain it: `insert_needle_seed`,
  `merge_adjacent`, `remove_film`, `remove_thin_layers`, and **`clamp_all`**
  (A2). This is the actual work of the item; the data model is trivial and
  the bookkeeping is not.
- **One predicate, added and unit-tested here, consulted everywhere else.**
  `Span::is_singleton_bulk(&self) -> bool` (the N1 arithmetic) and
  `DesignStack::span_of_row(&self, row) -> &Span`. Adding them changes
  nothing on its own, which is exactly why they belong in the bit-exact
  half; F0.2, F0.3, F1.6 and F2.3 are all callers.
- **→ F0.2 — not in this item.** The three behaviour changes the predicate
  enables (`remove_thin_layers` and `clamp_all` treating a span as one unit,
  A2; `inflate_design` / `round_to_qwot` excluding non-singleton-bulk spans
  from their selection lists, N3) are the whole content of F0.2. They are
  kept out of here so this item's fingerprint gate can be absolute.
- Needle host selection refuses non-singleton-bulk spans, naming the span's
  material and row range — likewise alongside the `needle` flag. **And the
  refusal exists at `insert_needle_seed` itself** (R3): that function is
  public, does no admissibility check of its own, and `build_scan_sites` is
  only one of its callers. One admissibility predicate, two call sites (§6
  row 6).
- **R2:** one named `assert_spans_partition(&self)` — contiguous, ordered,
  no gaps, no overlap, covers `0..films.len()` — called at the tail of every
  count-changing mutator **and both constructors**. A named helper survives
  the fifth mutator and the sixth; five independently-written asserts do not.

#### Gates

- **Fingerprints unmoved, both**, and named (R1):
  - random-stack: `validation/regression/structure/test_differential.py:95`
    and `:119`.
  - needle-run: **the new pin, committed first** (ground rule 5) — a fixed
    design, fixed seeds, hashed spectra, in
    `validation/regression/synthesis/`. Recorded before the first line of
    `DesignStack` changes, so it has an oracle rather than a hope.
- **N1 twin — the one that would have caught the mistake.** A plain film with
  `interface = true`, `needle = true`, thickness well above the scan step:
  needle sites are still generated inside its bulk, and its bulk row is still
  a thin-removal candidate. Before and after F0.1, bit-identical. This twin
  exists specifically to fail if anyone rewrites the predicate as
  `end - start > 1`.
- `assert_spans_partition` after every mutation. The assert is the proof that
  the bookkeeping is right — reasoning about five mutators is not.
- **Predicate unit tests:** `is_singleton_bulk` on a plain film, a plain film
  with `interface = true`, a graded film, and a graded film with
  `interface = true`. Four cases, no pipeline involved, and the second is the
  N1 case that would otherwise have been discovered by a moved fingerprint.
- **Needle-refusal twin (R3):** message names material + row range, asserted
  through `insert_needle_seed` **directly**, not only through the scan.
- **→ F0.2:** the thin-removal, clamp and inflate twins travel with the
  behaviour they test.
- **Merge twins (N2):** intra-span merge of bit-identical rows keeps the
  partition and merges into one span; cross-span merge lands in the first
  span. Existing `merge_adjacent` tests stay green unchanged — but do not
  claim, as revision 1 did, that they are sufficient.

**Risk.** M (was L, before the split). It still touches the needle host path
and the cleanup loop, but it no longer alters any outcome, so the fingerprint
gate here is unconditional: **any** movement in either fingerprint means the
change is wrong, full stop. There is no allowed-delta list for this item, and
that absence is the entire point of splitting it.

---

### F0.2 — span-level pipeline accounting (0.6.34)

**U5: a graded layer is one physical layer.** Four rules in the pipeline
disagree with that today, because each of them counts, compares or selects
*solver rows*. This item makes all four span-aware. It is the one release in
the series licensed to change what a run produces, and the licence is
enumerated below — anything outside that list is a defect.

#### The four rules, and what each does wrong today

Every number here was measured on the shipped 0.6.32 release build, not
derived from reading the source.

**1. The floor deletes the layer.** `clamp_all`
([structure.rs:381](rust/navette/src/smatrix/synthesis/structure.rs:381))
compares each *row* against `clamp_min_nm` and removes it. A 3-layer design
carrying a pinned 10 nm graded interlayer, run through
`clamp_all(2.0, 1000.0)`: **7 rows removed**, total thickness 220.0 → 210.0 nm,
the interlayer gone. Not thinned — gone, and the stack solves without it.

**2. The cap never fires.** The same function caps a row at `clamp_max_nm`.
A plain 1000 nm film under `clamp_max_nm = 300` comes back at 300.0 nm,
correctly. A *graded* 1000 nm film under the same cap comes back at
**1000.0 nm**, because no individual sublayer exceeds 300. The manufacturing
ceiling is silently absent on exactly the layers whose expansion made them
invisible to it.

Note the asymmetry this creates and how bad it is: the floor over-fires
(deletes a layer that is merely subdivided) and the ceiling under-fires
(ignores a layer that is genuinely too thick). Both errors come from the same
mistake and both are fixed by the same change.

**3. The layer budget counts rows.** `check_budgets`
([pipeline.rs:102](rust/navette/src/smatrix/synthesis/pipeline.rs:102))
compares `films().len()` against `max_film_layers`, default 40
([config.rs:84](rust/navette/src/smatrix/synthesis/config.rs:84)). A single
1000 nm graded film at δ = 0.5 expands to **57 rows**, so the pre-flight check
of macro cycle 1 returns `LAYER_BUDGET_REACHED` and the run terminates before
doing any work at all. The budget is a *manufacturability* limit — it means
physical layers, it always did, and one graded film is one physical layer.

**4. Nobody is told.** `clamp_all` returns `(n_removed, n_capped)`. All four
call sites discard it —
[pipeline.rs:192](rust/navette/src/smatrix/synthesis/pipeline.rs:192),
[:213](rust/navette/src/smatrix/synthesis/pipeline.rs:213),
[:273](rust/navette/src/smatrix/synthesis/pipeline.rs:273) and
[evaluator.rs:373](rust/navette/src/smatrix/synthesis/evaluator.rs:373). The
deletion in case 1 produces no warning, no counter and no entry in the phase
result. The user's design silently loses a layer and the report says nothing.

And it is hot: `evaluator.rs:373` is the last statement of
`optimize_thicknesses_report`, so the clamp runs after **every** thickness
optimization — seven call sites, dozens of times per macro cycle — not at the
three pipeline sites alone.

#### Change

- **The floor and the cap become span quantities.** `clamp_all` iterates
  spans, not rows. For a singleton-bulk span (every stack that exists today)
  row and span are the same object and nothing changes. For a multi-row span
  the comparison is against the span's **bulk** thickness, `Σ d_r` over
  `start..end` excluding the interface slice — the slice is not user
  thickness and never was.
  - `D_span < min_nm` → the **whole span** is removed, all rows in one
    operation, and it is reported (below). Not row-by-row, which is what
    produces the "removed 7 rows" nonsense above.
  - `D_span > max_nm` → **refused, not rescaled.** See the next bullet.
- **Too thick is a refusal, and it lands at `NeedlePipeline::new`**
  ([pipeline.rs:74](rust/navette/src/smatrix/synthesis/pipeline.rs:74)), which
  already holds both the stack and the config and already calls
  `cfg.validated()`. A pinned profile that exceeds the manufacturing ceiling
  is an authoring error: the user asked for a 1000 nm graded film and a 300 nm
  machine. Silently squeezing the profile into 300 nm answers a question
  nobody asked, and doing it row-by-row (the only thing `clamp_all` can do)
  warps the profile as well. Refuse at the door, name the span, name both
  numbers.

  Once F1.6 exists, a span the user has explicitly marked scalable *is*
  rescalable and the cap becomes a bound on its LM parameter instead. The
  refusal then applies to pinned spans only. Both halves are stated in F1.6 so
  the two items cannot drift apart.
- **The layer budget counts spans**, and so does the `layer_count` reported in
  `PipelinePhaseResult`. `max_total_thickness_nm` stays a row sum — a sum is a
  sum either way, and it is already correct; this is stated so that a later
  reader does not "fix" it into a discrepancy.
- **`clamp_all`'s return value stops being discarded.** It grows into a small
  `ClampReport { spans_removed: Vec<String>, spans_capped: usize,
  rows_removed: usize }`, threaded into `PipelinePhaseResult` at all three
  pipeline sites and surfaced in the result dict. A removed span names its
  material and its thickness. The evaluator site
  ([evaluator.rs:373](rust/navette/src/smatrix/synthesis/evaluator.rs:373))
  accumulates into the phase rather than reporting per call — it fires dozens
  of times per cycle and a per-call message would be noise.
- **`inflate_design` and `round_to_qwot` exclude non-singleton-bulk spans**
  from their selection lists (N3). One filter on one `collect()` at
  [inflate.rs:135](rust/navette/src/smatrix/synthesis/inflate.rs:135). Without
  it F0.2 protects a span from three operations and hands it to a fourth.
- **`remove_thin_layers` takes the span exemption** (the rule about spans,
  alongside the existing `optimize` flag filter, not replacing it).

#### The licence — the exhaustive list of what may change

A run whose stack contains **no** multi-row span must be bit-identical, full
stop. For a run that does contain one, exactly these may differ, and nothing
else:

1. A graded span whose total is below the floor is removed **as a unit**
   instead of shedding rows (previously: partial or total silent deletion).
2. A graded span whose total exceeds the ceiling now **refuses at
   construction** (previously: accepted and silently uncapped).
3. `LAYER_BUDGET_REACHED` fires later, or not at all, on stacks that used to
   trip it on row count (previously: terminated before cycle 1).
4. `layer_count` in every phase result drops to the physical-layer count.
5. An `enable_inflate` run leaves graded spans alone (previously: distorted
   their per-sublayer ratios).
6. The result dict gains a clamp report that was never there.

Each of the six gets a twin. A seventh difference means the change is wrong.

#### Gates

- **Fingerprints unmoved on the no-span path**, both, named as in F0.1. This
  is the one that says the refactor did not leak.
- **Clamp twin (A2), measured, not hypothetical:** a background-pinned 5 nm
  graded film (δ = 0.1 → four 1.25 nm sublayers) under the default
  `clamp_min_nm = 2.0`. Before: the film vanishes. After: it survives whole.
  Assert film count **and** total thickness, not row count.
- **Deletion-report twin:** the same design with the floor raised above the
  span total → the span is removed, and the phase result names it. The
  assertion is on the report, because a silent correct deletion and a silent
  wrong one look identical from the outside, which is how case 1 survived
  this long.
- **Cap-refusal twin:** a pinned 1000 nm graded span with
  `clamp_max_nm = 300` → `NeedlePipeline::new` returns `Err`, message names
  the span material, `1000` and `300`. Same span at `clamp_max_nm = 1500` →
  constructs and runs.
- **Budget twin:** the measured case — one 1000 nm graded film at δ = 0.5,
  57 rows, `max_film_layers = 40`. Before: `LAYER_BUDGET_REACHED` on the
  pre-flight of cycle 1. After: the run proceeds, and `layer_count` reads 1.
- **Inflate twin (N3):** a pinned span through an `enable_inflate` cycle →
  per-row thicknesses unchanged, ratios preserved.
- **Thin-removal twin:** a graded film whose sublayers are below the floor
  **and whose rows are marked `optimize = true`** → all sublayers survive.
  The case the flag filter does not cover.
- `assert_spans_partition` after every clamp pass. The span-removal branch
  deletes a contiguous block and renumbers everything after it; that is the
  single most error-prone edit in the item.

**Risk.** L. It changes live behaviour on purpose, in the hottest loop in the
pipeline, and the six-item licence is the only thing standing between "fixed"
and "quietly different". Mitigated by the licence being enumerated *before*
the code, by a twin per line of it, and by F0.1 having already proved the
bookkeeping separately.

---

### F0.3 — `ThinLayerPolicy`: clamp up instead of removing (0.6.35)

**U1.** "Remove a layer that got too thin" and "set that layer to the minimum
thickness you can actually deposit" are two different operations with two
different purposes, and only the first exists.

#### What the tree has today

`clamp_all`'s low branch is unconditional removal
([structure.rs:381](rust/navette/src/smatrix/synthesis/structure.rs:381)),
and its own doc comment states the intent plainly: layers below `min_nm` are
removed, *not* clamped up, because the optimizer driving one to zero is the
optimizer saying it wants that layer gone. That is a coherent design — it is
how thickness optimization eliminates layers at all — but it means
`clamp_min_nm` **reads** like a manufacturing floor and **behaves** as an
elimination threshold, while `clamp_max_nm` at the other end is a genuine
manufacturing ceiling. The pair is asymmetric and the naming hides it.

The only clamp-up anywhere in the tree is `round_to_qwot`
([inflate.rs:186](rust/navette/src/smatrix/synthesis/inflate.rs:186)),
`(round_half_even(ratio) * step).max(step)` — a different quantity (quarter-wave
multiples, not nanometres), and it has no pipeline caller.

#### Change

A new `PipelineConfig::thin_layer_policy`, mirrored in the Python
`PipelineConfig` and in `_smatrix.pyi`:

| Variant | Floor behaviour during the run | Floor behaviour at the end | LM lower bound |
|---|---|---|---|
| `Remove` *(default)* | remove | remove | `0.0` |
| `ClampUpFinal` *(recommended)* | remove | **set to `clamp_min_nm`** | `0.0` |
| `ClampUpAlways` | **set to `clamp_min_nm`** | set to `clamp_min_nm` | **`clamp_min_nm`** |

`Remove` reproduces today bit for bit and is what the fingerprints hold.

#### Why the LM lower bound is in that table — the limit cycle

Clamping up **without** raising the LM lower bound builds an oscillator. LM
drives a film to 0.5 nm because that is where the merit gradient points; the
clamp puts it back to 2.0; the next `optimize_thicknesses` drives it to 0.5
again. The merit alternates between two values forever and the design never
moves. Worse, it does not look like a bug from the outside:
`StagnationDetector` (`stagnation_oscillation_ratio = 0.75`, window 5,
[config.rs:100](rust/navette/src/smatrix/synthesis/config.rs:100)) sees exactly
the signature it is built to catch and terminates the run with
`STAGNATION_OSCILLATION` — a true report of a false condition.

So `ClampUpAlways` must also change `lb` at
[evaluator.rs:341](rust/navette/src/smatrix/synthesis/evaluator.rs:341) from
`vec![0.0; n]` to `vec![clamp_min_nm; n]`. That coupling is the actual content
of this item; the config enum is the easy part.

#### And why raising the bound costs the search

With `lb = clamp_min_nm` nothing can reach zero, so **nothing can ever be
eliminated**. `remove_thin_layers` finds no candidates, cleanup does nothing,
and needle synthesis loses its rejection mechanism: a seed goes in at
`needle_seed_thickness_nm = 5.0`
([cycle.rs:47](rust/navette/src/smatrix/synthesis/cycle.rs:47)) and a *bad*
seed is rejected by the optimizer shrinking it through the 2.0 nm floor. Pin
the floor as a hard bound and the bad seed parks at 2.0 nm permanently. Layer
count then only ever grows, straight into `max_film_layers`.

**`ClampUpAlways` is therefore refused at `NeedlePipeline::new` when
`needles_per_cycle > 0`**, with a message that says which two settings
conflict and names `ClampUpFinal` as the thing the user probably wants. It is
the mode for re-optimizing a *fixed* architecture, and the docs say so in the
same breath as introducing it.

`ClampUpFinal` has none of this: the search runs exactly as it does today,
elimination and all, and only the final clamp pass
([pipeline.rs:273](rust/navette/src/smatrix/synthesis/pipeline.rs:273))
substitutes clamping for removal. A film that survived the whole run and
happens to sit at 0.8 nm comes out at 2.0 nm and is depositable. A film the
optimizer genuinely wanted gone is already gone by then. This is the variant
the documentation leads with.

**The reported merit must be recomputed after the final clamp-up**, not
before. A number that describes a stack other than the one returned is worse
than no number.

#### Interaction with F0.2 — stated here so the two cannot collide

On a **span**, clamping up is not a row operation: it is "scale the whole span
to `clamp_min_nm`", which is precisely F1.6's operation and does not exist
yet at 0.6.35. Until F1.6 lands, an under-thickness graded span is removed
whole with the F0.2 report under **every** policy, and the doc comment on
`ThinLayerPolicy` says that in one sentence. F1.6 then wires the span branch
in and removes the sentence.

#### Gates

- **Default is bit-exact:** both fingerprints unmoved with
  `thin_layer_policy` absent and with it explicitly `Remove`.
- **`ClampUpFinal` twin:** a run that ends with a 0.8 nm film → the returned
  stack has it at exactly `clamp_min_nm`, the film count is unchanged, and
  `final_mf` equals a fresh merit evaluation of the returned stack (this is
  the assertion that catches reporting the pre-clamp number).
- **`ClampUpAlways` bound twin:** `lb` is `clamp_min_nm`, and a design whose
  optimum wants a 0.5 nm film converges to 2.0 nm in **one** optimization
  rather than alternating. Assert the merit history is monotone over five
  cycles — that is the limit cycle's absence, positively stated.
- **Conflict-refusal twin:** `ClampUpAlways` + `needles_per_cycle = 3` →
  `NeedlePipeline::new` returns `Err`, message names both settings.
- **Span deferral twin:** a thin graded span under `ClampUpFinal` at 0.6.35 →
  removed whole, reported. The same twin is *inverted* at F1.6 and that
  inversion is listed in F1.6's licence.

**Risk.** M. The enum is trivial and the coupling is not: an implementation
that adds the policy and forgets the lower bound produces a run that
terminates early with a plausible-looking stagnation report, which is the
hardest class of bug to notice. The monotone-merit twin exists specifically
for that failure.

---

## 3. Phase A — gradient-index and inhomogeneous layers

Physics reference: `docs/plans/gradient_layers_plan.md` §D0–D4, §D9. That
document's vocabulary is locked and not restated here. What follows is the
sequencing and the gates.

### F1.1 — data model + `FixedSpan` expansion + homogenize path (0.6.36)

**Scope.** `GradientSpec`, the EMA model selector, `GradientMode::FixedSpan`,
`ProfileShape::Linear`, validation, the expansion branch, and the
non-background homogenize path (§D4.3) so the feature is safe from the first
commit.

**Where.**

- `structure/layer.rs`: `Layer` gains `gradient: Option<GradientSpec>`,
  default `None`. New module `structure/gradient.rs` for the spec types —
  `layer.rs` is 606 lines and already carries the validation surface.
- `structure/expansion.rs`: new branch beside the `inhomogen` branch at
  [expansion.rs:257](rust/navette/src/structure/expansion.rs:257), mutually
  exclusive with it by validation.
- `materials/ema.rs`: reused unchanged. No new mixing math.
- `driver.rs`: the two A4 plumbing lines (background predicate disjunct,
  `assemble_stack` field copy).

**AMENDED (N5) — the validation "Where" names the wrong file, and the "no new
gate" claim cannot hold for two of the refusals.**

Revision 1 said gradient issues join `property_issues` in
`structure/validation.rs`. `structure/validation.rs` is 82 lines and contains
only the `ValidationIssue` type and its `gate`; `property_issues` lives at
[layer.rs:155](rust/navette/src/structure/layer.rs:155). Correct the file.

More substantively: `property_issues(&self, label: &str)` takes **no
provider**. So of F1.1's four refusals, two can be `property_issues` findings
and two cannot:

| Refusal | Home |
|---|---|
| `gradient` + `inhomogen` together | `property_issues` — self-contained |
| `f_*` outside `[0, 1]` | `property_issues` — self-contained |
| `material_a == material_b` | `property_issues` — self-contained |
| **a material absent from the provider** | **cannot** — needs the provider |

The provider-existence refusal lands where the provider is: `expand`
(structure path) and `from_design` / `assemble_stack` (design path). That is a
second site, and pretending otherwise means it gets written twice or not at
all. State it: *one* rule surface for the self-contained numbers
(`property_issues`, covering the PyO3 constructor, the gating setters,
`set_properties`, `from_state` and `assemble_stack`), plus *one* resolution
check at the two provider doors, sharing a message helper so the wording is
identical.

**AMENDED (A9.1) — the setter count is five, not four.** Twelve setters are
exposed on the native `Layer`
([structure.rs:230–355](rust/navette-py/src/structure.rs:230)); **five** call
`gate_layer`: `thickness`, `inhomogen`, `inh_delta`, `roughness`,
`interface_thickness`. Four of them are numeric, which is what revision 1's
"four" meant and what the amendment's rewording assumed — but `inhomogen` is
a bool and gates too, which matters here precisely because `gradient` is the
next non-numeric field that must gate. Say "the five gating setters"; the
new `gradient` setter makes six.

**AMENDED (A5) — two nk spectra must reach the design path, and the decision
cannot wait for F1.5.**

The structure path is fine: `Layer.gradient` names `material_a`/`material_b`
and a `MaterialProvider` resolves both. The synthesis **design** path has no
materials library at all:

- `_film_dicts` ([pipeline.py:91–121](src/navette/synthesis/pipeline.py:91))
  evaluates **one** material per film into `nk` (`_eval_nk(mat, wl)`, `:117`),
  and refuses any key not in `_FILM_DEFAULTS`
  ([pipeline.py:85](src/navette/synthesis/pipeline.py:85)) — twice, once for
  global flags and once per film.
- `ArrayFilm` ([driver.rs:32](rust/navette/src/smatrix/synthesis/driver.rs:32))
  carries one `nk: Vec<Complex64>`.
- `assemble_stack` registers provider entries as film-name → film-nk only
  ([driver.rs:124](rust/navette/src/smatrix/synthesis/driver.rs:124)); the
  film's `material` **is** its name.

So a gradient design film's `material_b` has no nk anywhere in the design
path — and F1.1's own scope (expansion branch, homogenize path in
`from_design`) consumes it from the first commit that exercises a gradient
design film. **Resolved at F1.1** (§8.8):

- `ArrayFilm` gains `gradient: Option<GradientJson>`; the JSON carries
  `material_b` **and** `nk_b`, evaluated Python-side in `_film_dicts` at the
  door; `assemble_stack` injects `material_b → nk_b` into the provider
  entries when a gradient is present. This keeps the path self-contained —
  its existing convention is that every nk is film-supplied — and keeps
  §D2's same-provider refusal as the gate.
- `gradient` joins `_FILM_DEFAULTS` (default `None`) or the unknown-flag
  refusal rejects it before it reaches Rust.
- The design film keeps its **own name** (it is the contrast/needle identity,
  multi-env §2.3); `material_a`/`material_b` are provider keys only. Where
  `material_a` names the film itself, the nk is already registered and no
  second entry is needed.
- Refuse when `nk_b`/`material_b` is absent — never fall back to
  `material_a`'s nk. A half-mixture is worse than a refusal (§D2's words).

**AMENDED (A6, corrected) — reuse `MixRule`; do not declare a second enum.**

The amendment noted that Mori-Tanaka is in-tree and excluded from §D2's
four-variant `EmaModel`. Re-checking found the situation is broader and the
resolution cleaner. **Six** kernels are in-tree and Python-bound:
`bruggeman`, `maxwell_garnett`, `looyenga`, `lichtenecker`, `mori_tanaka`
**and `general_power_law` (Birchak)**
([ema.rs:23,37,52,70,87,136](rust/navette/src/materials/ema.rs:23);
[materials/__init__.py:232](src/navette/materials/__init__.py:232) dispatches
all six by name).

More to the point, **the enum already exists**: `MixRule`
([materials/mod.rs:54](rust/navette/src/materials/mod.rs:54)) has all six
variants and carries their parameters in the variant
(`Bruggeman { max_iter, tol }`, `MoriTanaka { l }`, `PowerLaw { alpha }`).
§D2's `EmaModel` would be a parallel, poorer copy of it — and note that even
§D2's "parameter-free four" is not parameter-free: Bruggeman already carries
two.

**Resolution: `GradientSpec.ema` is a `MixRule`.** All six work on day one,
the shape factor and exponent problems disappear (they live in the variant),
and there is one vocabulary for mixing in the crate instead of two.
Cost: `MixRule` derives only `Clone, Copy, Debug` today and needs
`Serialize`, `Deserialize` and `PartialEq` — a small, contained change, and
the serde representation is the schema-visible surface, so name the variants
and their fields deliberately. Default stays **Bruggeman** (§D9.1, locked).
This closes §8.9 rather than documenting an omission.

**The one open formula, resolved here.** §D3 step 2 leaves the sublayer-count
rule as "exact formula fixed at implementation". Fixing it now:

```
n_sub = clamp(ceil(thickness_nm / max_step_nm), 3, 64)
max_step_nm = min(20.0, lambda_min / (10.0 * n_max_re))
```

with `sublayers: Some(n)` overriding, clamped to `[2, 256]` with a warning
outside. Rationale: the step must be small against the local optical
wavelength or the staircase itself becomes the physics; `lambda_min / 10n` is
the conventional floor and the 20 nm cap keeps thick films from exploding the
row count. Pinned by a differential test over randomized thicknesses and
grids — integer boundaries must agree bit-exactly, the same technique that
pins `sub_layer_count`.

**Note the interaction with F0.1's clamp fix.** The `min 3` floor means any
gradient film thinner than 6 nm has sublayers below the default
`clamp_min_nm = 2.0`. Without F0.1 shipped first, such a film is erased by the
first clamp of the run. This is the concrete reason F0.1 is P0 and sequenced
before this item.

**AMENDED (N7) — the two sublayer-thickness conventions differ, on purpose.**
The legacy `inhomogen` branch divides uniformly:
`step_t = layer_thickness / f64::from(sub)`
([expansion.rs:273](rust/navette/src/structure/expansion.rs:273)), so `Σd`
is `thickness` only to float precision. §D3 step 3 and F1.1's gate specify
that the gradient branch has the **last sublayer absorb the remainder** so
`Σd ≡ thickness` exactly. Both are defensible; they are different rules in
adjacent branches of the same function. Say so in a comment at the branch, or
the next reader "fixes" one to match the other and moves the legacy
fingerprint.

**Gates.**

- Both fingerprints **unmoved** (gradient absent ⇒ `None` ⇒ untouched path),
  named per R1.
- Per-row nk vs a direct `ema_*` oracle on the same `(nk_a, nk_b, f_i)` —
  **bitwise**, same kernel, same inputs.
- `f_i` recomputed from the mode formula, EMA outputs compared to 1e-15.
- Thickness-independence: same `FixedSpan` spec at 50 nm and 500 nm →
  endpoint rows bitwise equal, counts differ.
- Sum of sublayer thicknesses `== thickness` exactly (N7).
- Inversion reverses row order; profile direction flips; tested both ways.
- Roughness lands on the first sublayer only (legacy convention).
- Refusals, each naming both sides: `gradient` + `inhomogen` together;
  `material_a == material_b`; `f_start == f_end` (N2 — a degenerate span
  merges away silently otherwise); a material absent from the provider (at
  the provider doors, N5); `f_*` outside `[0, 1]`.
- Homogenize warning names the mixture and `f_mid`.
- **R4 — homogenize cannot reuse the flag flip.** Today's path mutates a
  clone (`h.inhomogen = false`,
  [structure.rs:236](rust/navette/src/smatrix/synthesis/structure.rs:236))
  and re-expands uniform, which works because a base nk exists. **A gradient
  film has no base nk.** Homogenize must synthesise the row itself — EMA at
  `f_mid` over both provider spectra, which A5 guarantees are present. Stated
  in the "Where" so nobody reaches for the existing flip. Gate: the
  homogenized row is bitwise equal to a direct `ema_*` call at `f_mid`.
- `check_exposure.py` / `check_pyi_sync.py` clean.

**Docstring obligation (§D9.4).** `GradientSpec` must carry, at the call
site, the statement that `material_a` is **always the host** and `f` is the
volume fraction of B in A — Maxwell-Garnett and Mori-Tanaka are
host/inclusion-asymmetric, so A↔B is not the same physics and a silent swap
corrupts fits. This is a docstring requirement, and it is checked in review,
not by a tool. With A6 resolved toward `MixRule`, Mori-Tanaka is *available*,
which makes this obligation load-bearing rather than theoretical.

### F1.2 — gradient `RateCapped` mode (0.6.37)

`f(z) = f_start + rate * (z / ref_thickness)`, clamped to `[f_min, f_max]`.
`ref_thickness` default 100 nm, caps default `[0.0, 1.0]`, `rate` finite with
sign free.

**AMENDED (N2) — the saturation gate and `merge_adjacent` are the same fact.**
The proof that saturation works ("the tail rows are bitwise equal to pure
`material_b`") is also the proof that `merge_adjacent` will collapse them: it
keys on `material == && nk ==`, and all sublayers of one carrier share the
carrier's material name. `cleanup_design` runs it unconditionally, twice.
F0.1 must already handle intra-span merges for this item to be safe; this
item adds the twin that exercises it end to end.

**Gates.** 200 nm film, `f_start = 0.1`, `rate = 0.3` per 100 nm → ends at
0.7. Same spec at 400 nm → saturates at 1.0 from z = 300 nm, and the tail
rows are **bitwise equal to pure `material_b`** — that is the saturation
proof, not an approximate comparison. Negative rate saturates to pure
`material_a` at the other end.

**Merge twin (N2):** the saturated 400 nm film through a full
`cleanup_design` — the tail merges to one row, the span stays contiguous and
keeps its `logical`, `assert_spans_partition` holds, and the **simulated
spectra before and after the merge are bitwise equal** (the merge is optically
a no-op; if it is not, the bookkeeping is wrong).

Refusal when two slope-like keys are supplied for one layer, **A9.2: in
F1.2's ASCII form, not §D0's**:

```
gradient: 'rate' and 'delta' are the same slope - give one.
```

§D0's source version uses an em-dash. Implement the ASCII form (ground rule
7); do not copy §D0's punctuation, and do not take
[structure.rs:230](rust/navette/src/smatrix/synthesis/structure.rs:230)'s
existing em-dash as precedent. No silent precedence between the two keys.

**Risk.** M (was S) — the saturation behaviour that makes the feature
demonstrable is the same behaviour that exercises F0.1's least-tested
bookkeeping path.

### F1.3 — `InhMode::RateCapped` for single-material drift (0.6.38)

`delta_layer = min(rate * thickness / ref_thickness, cap)`, replacing the
constant `inh_delta` as the input to the **frozen** legacy arithmetic.

The combination order is unchanged and stays exactly as
[expansion.rs:259](rust/navette/src/structure/expansion.rs:259) has it:

```
delta_nominal = (delta_layer + group.inh_delta_summand) * 0.5
```

where `delta_layer` is `inh_delta` in `Fixed` mode and the rate formula in
`RateCapped`. In `RateCapped` only, `delta_nominal` is clamped to `[0, cap]`
**before** the stochastic error draw. The draw itself is not clamped —
clamping draws would bias Monte-Carlo statistics. Documented, and pinned by a
statistical test that the mean survives the clamp boundary.

`sub_layer_count()` consumes `delta_layer(thickness)`, a pure function of
thickness, so there is no circularity; evaluated once. Note that
`sub_layer_count` reads `self.inh_delta` directly today
([layer.rs:127](rust/navette/src/structure/layer.rs:127)) and has three other
callers — [expansion.rs:257](rust/navette/src/structure/expansion.rs:257),
[structure.rs:272](rust/navette/src/structure/structure.rs:272) and
[architect.rs:565](rust/navette/src/structure/architect.rs:565), the latter
two for row-count prediction. All four must see the same `delta_layer`, or the
predicted row count and the emitted row count diverge.

**Gates.** `Fixed` mode **bitwise** against the legacy oracle over the
existing randomized differential suite
(`test_differential.py:95`, `:119`) — this is the whole risk of the item and
the only acceptable evidence. `RateCapped` delta against hand-computed
`min(rate*t/t_ref, cap)` over a thickness sweep, with the saturation knee
exact. Double the thickness below the cap → double the delta. Row-count
prediction (`structure.rs`, `architect.rs`) agrees with emission in both
modes.

**Note.** `inhomogen` is **not deprecated**. Single-material drift is real
physics (oxidation gradients, nitrides) and coexists with gradients by the
§D2 rule.

### F1.6 — one thickness parameter per graded span (0.6.39)

**U2.** A graded layer is one physical layer, so its total thickness is one
number the optimizer should be allowed to move. Today it is either pinned
(`optimize = false, needle = false`, the background path) or homogenized away.
This item adds the third option for the profiles that scale exactly, and F1.7
adds it for the ones that do not.

#### Which profiles scale exactly, and why

A profile scales iff the index at a sublayer is a function of that sublayer's
**fractional** position in the layer, not of its absolute depth. Then
stretching the layer moves every sublayer's *thickness* and leaves every
sublayer's *index* alone.

| Mode | Fraction depends on | Scales freely | Item |
|---|---|---|---|
| `InhMode::Fixed` (ships today) | `i / (sub-1)` | **yes — verified** | F1.6 |
| `GradientMode::FixedSpan` (F1.1) | normalized depth between two endpoints | yes, by construction | F1.6 |
| `InhMode::RateCapped` (F1.3) | absolute `z` | no | F1.7 |
| `GradientMode::RateCapped` (F1.2) | absolute `z` | no | F1.7 |

The first row is not an assumption. [expansion.rs:265](rust/navette/src/structure/expansion.rs:265):

```rust
let mut factors: Vec<f64> = (0..sub)
    .map(|i| 1.0 - current_delta + 2.0 * current_delta * f64::from(i) / f64::from(sub - 1))
    .collect();
```

`factors[i]` is a function of `current_delta` and `i / (sub-1)` and of nothing
else. Thickness enters the branch exactly once, two lines later, as
`step_t = layer_thickness / f64::from(sub)`. So today's inhomogeneous profile
is already scale-free at a fixed row count, which is what makes this item
mostly plumbing rather than physics.

There is also precedent for scaling a whole layer: `Group::thick_factor` and
`thick_summand`, applied at
[expansion.rs:128](rust/navette/src/structure/expansion.rs:128) as
`layer.thickness * group.thick_factor + group.thick_summand` — a per-layer
multiplier that has always operated on the logical layer, before expansion.

#### Change

- **The LM parameter list stops being a row list.** `opt_indices: Vec<usize>`
  at [evaluator.rs:310](rust/navette/src/smatrix/synthesis/evaluator.rs:310)
  becomes `params: Vec<Param>`, where a `Param` is either `Row(usize)` or
  `Span { rows: Vec<usize>, fractions: Vec<f64> }`. A stack with no scalable
  span produces only `Row` params and the list is the same list it is today.
- **`fractions` are frozen when the parameter list is built** — `φ_r = d_r / D₀`
  over the span's bulk rows, captured once per `optimize_thicknesses` call.
  The residual closure at
  [evaluator.rs:327](rust/navette/src/smatrix/synthesis/evaluator.rs:327) then
  writes `set_thickness(r, φ_r · D)` for each row of the span. The interface
  slice is not in `rows` and is not scaled: it is an interface property, not
  part of the layer's thickness.
- **Which spans are scalable.** A graded carrier with `optimize = true`. That
  is the existing flag, given a meaning it did not have before — previously
  `optimize = true` on a graded film meant "homogenize me" — it falls out of
  the background set at
  [driver.rs:130](rust/navette/src/smatrix/synthesis/driver.rs:130) and is
  flattened at
  [structure.rs:236](rust/navette/src/smatrix/synthesis/structure.rs:236) with
  the profile dropped and a warning.
  **This changes the meaning of an existing flag combination**, so it is in
  F1.6's licence below and it is the one part of this item that is not
  additive. The homogenize warning stops firing for graded films that reach
  the new path; F1.1's `optimize = false, needle = false` pinned path is
  untouched.
- **Bounds are span quantities**, matching what F0.2 already made the clamp
  compare against: `lb = 0.0` (or `clamp_min_nm` under `ClampUpAlways`),
  `ub = clamp_max_nm`, both applied to `D`. The manufacturing ceiling that
  never fired on a graded film (F0.2 case 2) now fires correctly, as a bound
  rather than as a refusal — and F0.2's construction-time refusal narrows to
  **pinned** spans only, which is stated in both items.
- **F0.3's span branch lands here.** An under-thickness scalable span under
  `ClampUpFinal` is scaled to `clamp_min_nm` instead of removed, and F0.3's
  "removed whole under every policy" sentence is deleted from the doc comment.
- Needle host selection is **unchanged**: a scalable span is still not a host.
  Scaling a profile preserves it; splitting a foreign row into the middle of
  it does not.

#### The analytic Jacobian stays analytic, stays exact, and stays cheap

`∂/∂D = Σ_r φ_r · ∂/∂d_r` — with `φ_r` frozen and the nk frozen, this is
exact, not an approximation. Every `∂/∂d_r` is already what
`simulate_with_deposits` produces.

The cost is the surprise, and it is a good one. The deposits are documented at
[evaluator.rs:58](rust/navette/src/smatrix/synthesis/evaluator.rs:58) as "one
extra `StackFields` decomposition per point and polarization, plus **O(1) per
film**" — the whole Jacobian is a constant multiple of *one* simulate. So
adding a 57-row span to `par_films` adds 57 O(1) terms, not 57 simulates.
Central differences on the same parameter would cost two full simulates.
**Analytic is both exact and cheaper here, and there is no reason to decline
to finite differences** — which is worth saying explicitly, because the
obvious first guess is the opposite.

Mechanically: `assemble_jacobian`
([jacobian.rs:120](rust/navette/src/smatrix/synthesis/jacobian.rs:120)) gains
the parameter map and contracts `dep.column(channel, point)` (one entry per
*row*) into `out` (one entry per *parameter*) with the weights. A `Row` param
is a one-element group with weight `1.0`, and `1.0 * g` summed over one
element is bit-identical to `g` — that is the no-scalable-span bit-exactness,
and it is asserted, not assumed.

#### Licence — what a run may do differently

Only for a stack that has at least one graded film with `optimize = true`:

1. That film's thickness now moves during optimization (previously: it was
   homogenized to its base index with a warning).
2. The homogenize warning no longer fires for it.
3. Its total is now subject to `clamp_max_nm` as a bound (previously: F0.2
   refused it at construction).
4. Under `ClampUpFinal` an under-thickness one is scaled up rather than
   removed (inverts F0.3's span-deferral twin, deliberately).

A stack with no such film is bit-identical, including the Jacobian, to the
last bit.

#### Gates

- **The scale-invariance twin — this is the physics claim, and it is written
  before the feature.** Expand an `InhMode::Fixed` film at 100 nm and at
  200 nm with the row count pinned to the same value: the two nk row vectors
  are **bitwise equal** and the row thicknesses are exactly doubled. If this
  fails, F1.6 is not legal and nothing else in the item matters.
- **Jacobian twin, two ways:** the analytic span column against a central
  finite difference on `D` to 1e-6 relative; **and** the analytic span column
  against `Σ φ_r ·` (the analytic row columns) **bitwise**. The second is the
  one that catches a wrong weight.
- **Bit-exactness twin:** both fingerprints unmoved, plus a direct assertion
  that the assembled Jacobian of a no-span stack is bitwise identical before
  and after the parameter-map refactor.
- **Profile-preservation twin:** after an LM solve that moves `D` by 40 %,
  every `d_r / D` is unchanged to the last bit and every nk row vector is
  untouched. This is the assertion that the layer scaled rather than deformed.
- **Interface twin:** a scalable span with `interface = true` → the slice row's
  thickness is unchanged by the scale, and the span partition still holds.
- **Bound twin:** the F0.2 cap case (1000 nm graded, `clamp_max_nm = 300`),
  now marked `optimize = true` → constructs instead of refusing, and the
  optimizer never returns a `D` above 300.

**Risk.** M. The refactor reaches into the two hottest functions in the
optimizer, but every step of it has a bitwise oracle: the row-only path must
not move at all, and the span path has a closed-form check against the row
path it is built from. The real risk is the flag-meaning change (licence
item 1), which is why it is called out as non-additive rather than buried.

---

### F1.7 — profile refresh for the rate modes, at construction only (0.6.40)

**U3.** A rate-type profile is a function of absolute thickness — `f(z) =
f_start + rate · z / ref_thickness` — so scaling the span makes its stored nk
stale. Those spans need the profile recalculated.

**U4, and it governs the whole item.** That recalculation runs at
**construction points only**. Never inside an LM round, never inside the
Jacobian, never inside a needle scan or a cleanup trial. This is stated first
because it is the constraint the design is built around, not a tuning
decision taken afterwards.

#### Why the obvious implementation is the wrong one

Re-expanding inside the residual closure
([evaluator.rs:327](rust/navette/src/smatrix/synthesis/evaluator.rs:327)) is
what "recalculate after a scale" naively means, and it would re-run the
sublayer-count rule and the EMA kernel on **every residual evaluation and
every Jacobian call** — order 10³ times per `optimize_thicknesses`, which
itself runs seven times per macro cycle — to chase thickness excursions of a
fraction of a nanometre. Re-expansion is a construction operation and it stays
one.

**And the cost is not the strongest argument against it. Correctness is.**
`n_sub` comes from a `ceil()` ([layer.rs:126](rust/navette/src/structure/layer.rs:126)).
If the row count could change mid-solve, then:

- the residual vector changes length between iterations, which the LM backend
  is not built for; and
- the merit surface acquires a **step discontinuity** exactly at each `ceil`
  boundary. LM would see a jump in the residual that is not physics, it is
  discretization, and it would either take a garbage step or stall against it.

So the row count of a span is **frozen for the duration of one LM solve** —
a hard requirement, independent of performance. The user's constraint and the
numerics point the same way, which is the comfortable case.

#### The refresh points — the exhaustive list

`refresh_profiles(&mut stack)` is called at exactly three places:

1. **`DesignStack::from_design`**
   ([structure.rs:176](rust/navette/src/smatrix/synthesis/structure.rs:176)) —
   initial assembly. Already the natural expansion site.
2. **The top of each macro cycle**, immediately before the pre-flight budget
   check ([pipeline.rs:139](rust/navette/src/smatrix/synthesis/pipeline.rs:139)).
   One call per cycle, so the budget check and the needle scan both see a
   current profile.
3. **After the final clamp pass**
   ([pipeline.rs:273](rust/navette/src/smatrix/synthesis/pipeline.rs:273)) —
   so the returned stack and the reported `final_mf` describe the same object.

And nowhere else. In particular **not** in: the residual closure
([evaluator.rs:327](rust/navette/src/smatrix/synthesis/evaluator.rs:327)),
`DepositJacobian::fill`
([evaluator.rs:394](rust/navette/src/smatrix/synthesis/evaluator.rs:394)),
`build_scan_sites`, `cleanup_design`'s trial removals, or `inflate_design`.

#### What refresh does

For each rate-mode span, taking the span's current total `D`:

1. re-derive `n_sub` from the sublayer rule at the new `D`;
2. re-evaluate `f(z)` at the new absolute depths and re-run the EMA kernel;
3. rebuild the span's rows in place.

Row count may change **here and only here**, so `assert_spans_partition` runs
immediately after, and the F1.6 `fractions` are re-captured on the next
`optimize_thicknesses` call rather than carried across a refresh.

Fixed-mode spans are visited and left alone: for them refresh is a bitwise
no-op, because their profile was never a function of `D` in the first place.

#### What "stale within a cycle" actually costs — with a number

Inside one macro cycle the profile is the one built at the top of that cycle,
so LM minimizes the merit of a stack whose index profile is slightly behind
its thickness. The size of the error is bounded and small: `f` moves by
`rate · ΔD / ref_thickness`, so at the default `ref_thickness = 100 nm`, a
5 nm excursion moves the end fraction by `0.05 · rate` — for a typical
`rate = 0.3` that is a mixing-fraction error of 0.015 at the deepest sublayer
and less everywhere above it.

Two things make this an accepted approximation rather than a defect:

- the next cycle's refresh reconciles it, and macro cycles are short; and
- **every merit number shown to the user is evaluated after a refresh**, never
  on a stale profile — refresh point 3 exists for exactly that.

Documented in the user-facing docs as a stated property of rate-mode spans,
with the formula, not left for someone to discover.

#### The consequence nobody expects: the Jacobian stays analytic

With the profile frozen for the solve, `∂n/∂D` is zero *by construction of the
objective*. The stack LM is minimizing genuinely is the frozen-profile stack,
so F1.6's exact contraction `∂/∂D = Σ_r φ_r · ∂/∂d_r` is the exact derivative
of the actual objective — not an approximation of a live-profile derivative.

Had the profile been live, the derivative would carry an extra `∂n/∂d` term
that the deposit machinery does not produce, and rate-mode spans would have
had to decline to finite differences (the `Ok(None)` path that
`DepositJacobian` already uses for phase targets,
[evaluator.rs:406](rust/navette/src/smatrix/synthesis/evaluator.rs:406)) —
two full extra simulates per parameter per iteration. **U4 buys the analytic
Jacobian for rate modes.** It is worth recording that the constraint is a
speedup twice over, not a compromise.

#### Gates

- **The refresh-count twin — this is the test that enforces U4 mechanically.**
  Put a counter on the re-expansion entry point. Run a 3-macro-cycle pipeline
  with one rate-mode scalable span and assert the counter reads exactly
  **`1 + 3 + 1 = 5`**: the construction, the three cycle tops, the final pass.
  If anyone later wires refresh into the residual closure, the count jumps by
  three orders of magnitude and this test fails loudly. Written before the
  feature, and it is the item's most important test.
- **Row-count-frozen twin:** across one `optimize_thicknesses` call the
  residual vector length is constant on every iteration, including a run
  deliberately parameterized to cross a `ceil()` boundary in the sublayer
  rule mid-solve.
- **Refresh-correctness twin:** a rate-mode span scaled by hand from 100 nm to
  200 nm, then refreshed, is **bitwise equal** to the same span expanded from
  scratch at 200 nm. Refresh and construction must not be two different
  functions in disguise.
- **Fixed-mode no-op twin:** `refresh_profiles` on a stack of `Fixed` and
  `FixedSpan` spans is bitwise identity, over the randomized differential
  suite (`test_differential.py:95`, `:119`).
- **Drift-bound twin:** the merit recomputed immediately after a refresh
  differs from the pre-refresh merit by less than a documented tolerance for
  `rate` in the documented range. If it does not, the cycle is too long for
  that rate, and the docs say so rather than the code hiding it.
- **Reported-merit twin:** `final_mf` equals a fresh merit evaluation of the
  returned stack. Same assertion shape as F0.3's, for the same reason.
- Fingerprints unmoved: no existing stack has a rate-mode span.

**Risk.** M. The physics is a re-run of F1.1–F1.3 code and the plumbing is
three call sites. The risk is entirely that a fourth call site appears later —
someone fixes a staleness symptom by refreshing where the symptom showed up —
and the refresh-count twin is the guard against precisely that.

---

### F1.4 — schema v2 and a readable version range (0.6.41)

The correction from §1.2 above. **Sequenced before the Python surface**
because every state file on disk reads through this gate, and getting it
wrong is a data-loss-shaped bug rather than a feature defect.

**AMENDED (A3) — the gate has a Python half, and the range breaks three
existing tests, not one.**

Revision 1's change list named `version.rs` only. The full surface:

1. **`rust/navette/src/structure/version.rs`** — the constant, the range, the
   two distinct messages.
2. **`src/navette/structure/types.py:119`** — `SCHEMA_VERSION = 1`, a second
   literal. (The Python `check_schema_version` at `:123` is thin over the
   native gate — that half is fine; the constant is not.) The Python side
   also needs `MIN_READABLE_SCHEMA_VERSION` for the sync test to probe both
   endpoints.
3. **The policy comment, in *two* places.** `types.py:113–118` says "purely
   additive keys are safe without a bump (readers ignore unknown keys) — the
   fingerprint test in `test_roundtrip.py` enforces this decision on every
   key-set change." `test_roundtrip.py:70–76` says the same thing again from
   the other side. F1.4's fields *are* additive-with-defaults and this plan
   bumps anyway, because the newer-writer hazard justifies it. **Both
   comments are rewritten in the same commit**, or the code contradicts its
   own documented policy the week it lands.
4. **`validation/smoke/test_request_bits.py`** — and this is the part the
   amendment understated. The sync test does not read a constant; it
   *probes* the gate through `_accepted_version`
   ([test_request_bits.py:232](validation/smoke/test_request_bits.py:232)),
   which asserts `len(ok) == 1` — "the one version a native gate accepts."
   **A range gate fails that assertion by construction**, and it fails for
   the state gate at this item and for the program gate at F2.4. The helper
   becomes `_accepted_range(gate) -> (lo, hi)`; both tests move to asserting
   the pair. F1.4 writes the helper; F2.4 reuses it. Note the shared file:
   F1.4 and F2.4 both touch `test_request_bits.py` and must not conflict.
5. **`validation/regression/structure/test_roundtrip.py:51`** —
   `test_stale_schema_versions_refused` asserts that `SCHEMA_VERSION - 1` is
   **refused**. After the bump that is version 1, which must now be
   **accepted**. The assertion inverts. Change it deliberately: v1 accepted,
   v0 refused as stale, `SCHEMA_VERSION + 999` refused as newer-build, and
   an untagged state still refused as malformed.
6. **`test_state_fingerprint`** (`test_roundtrip.py:96`) asserts
   `FINGERPRINT["version"] == SCHEMA_VERSION`; the recorded key lists move to
   version 2 and `Layer` gains the new keys.

**Also (A3):** revision 1's gate sentence "`test_roundtrip.py` and
`test_program.py` extended rather than replaced" mixes schemas.
`test_program.py` (`validation/regression/config/`) exercises the *program*
schema. This item extends `test_roundtrip.py` only; the `test_program.py`
half moves to F2.4, where A1's bump lives.

**AMENDED (N4) — there are no serde defaults to add, because there is no
derive.**

`Layer`'s `Serialize` and `Deserialize` are **hand-written**
([layer.rs:335](rust/navette/src/structure/layer.rs:335) and
[layer.rs:356](rust/navette/src/structure/layer.rs:356)). Every field is read
by an explicit `map.get(k)` and every field is *required* — there is no
`#[derive(Deserialize)]` and therefore no `#[serde(default)]` to hang a
default on. §1.2's phrase "fields added in v2 carry serde defaults" describes
a mechanism this type does not have. Concretely:

- new fields are read as `map.get("gradient").map(...).transpose()?` with an
  explicit `None` fallback, matching the surrounding style;
- `s.serialize_map(Some(13))` ([layer.rs:338](rust/navette/src/structure/layer.rs:338))
  is a hard-coded count and must move, or become `None` (serde permits an
  unknown-length map);
- `GradientSpec` itself may derive serde normally — the hand-written impl is
  `Layer`'s alone.

**And a decision this forces (§8.10): emit `gradient` always, or only when
`Some`?**

- *Always* (as `null`): every v2 state file differs from its v1 equivalent
  even for stacks with no gradient; `test_state_fingerprint`'s flat key list
  grows unconditionally; every state file on disk changes shape at this
  release.
- *Only when `Some`*: a stack with no gradient serialises **byte-identically**
  at v1 and v2 apart from the version tag. `test_state_fingerprint` needs a
  conditional (two key lists: with and without a gradient).

**Recommend only-when-`Some`.** It is the plan's own bit-exactness ethic
applied to the state file, and it makes "v2 changes nothing unless you use
the feature" literally true rather than nearly true. Cost is one conditional
in one test.

**Also pin the nested key set.** `test_state_fingerprint` compares *top-level*
keys only. `gradient` is a nested object, so its internal key set is
unpinned by the existing mechanism — add a `FINGERPRINT["GradientSpec"]` list
and assert it when the field is present, or the new field has weaker
protection than every old one.

**R5 — sequence the fixture first.** Commit the v1 state fixture, and write
the refusal-side test (a v3 state naming the newer build) **against the
un-bumped gate**, before `version.rs` changes. The acceptance test then has
its oracle from the first line of the change rather than from the change
itself.

**Gates.** A v1 state fixture committed to the tree loads and expands
bit-identically. A `schema_version: 3` state refuses, with a message that
says the state was written by a newer build — not "stale". A
`schema_version: 0` state refuses as stale. Round-trip v2 → v2 preserves
every new field. `_accepted_range` asserts `(1, 2)` for the state gate.
`test_roundtrip.py` extended, not replaced. Both fingerprints unmoved.

### F1.5 — config rows and Python surface (0.6.42)

`design_config.rs::LayerRow` gains `gradient: Option<GradientJson>` next to
the existing `inhomogen` flag at
[design_config.rs:85](rust/navette/src/smatrix/synthesis/design_config.rs:85).
`ArrayFilm` in `driver.rs` gains the same — **already decided at F1.1 (A5)**;
this item mirrors the decision into `LayerRow` rather than making it.
`builders.py` passes through; Python `Layer` gains `gradient`, validated
natively — no Python pre-checks, no pydantic.

**Note (N9).** `LayerRow` carries `#[serde(deny_unknown_fields)]`
([design_config.rs:74](rust/navette/src/smatrix/synthesis/design_config.rs:74)),
so a v2 request meeting a v1 binary **refuses loudly** on this path. No range
gate is needed here and none should be added; the contrast with the state path
(§1.2) is worth one comment at the struct.

**Gates.** Python surface ≡ JSON path, hex-compared. Refusals arrive at
construction as `ValueError` with the native message, ASCII. The new
`gradient` setter gates (making six gating setters, A9.1). `.pyi` stubs
updated, `check_pyi_sync.py` clean. `_FILM_DEFAULTS` carries `gradient` and
the unknown-flag refusal still fires for genuinely unknown keys.

---

## 4. Phase B — multi-environment optimization

Design reference: `docs/plans/multi_environment_plan.md` §3–§6, including the
five rejected alternatives (§3.2) and the eight non-negotiable v1 limitations
(§6). Not restated. What follows is sequencing, gates, and the interactions
with Phase A.

**Why after gradients.** A design film is one row today. After Phase A it can
be a span of up to 64. The multi-environment routing table maps a shared
design parameter to its absolute position in each environment's assembly — if
that table is built row-keyed before gradients land, it has to be rebuilt
span-keyed afterwards. Built once, after, it is `design_slot → (env, span)`
and correct for both. That is the whole ordering argument, and it is why
F0.1 is a shared prerequisite rather than a gradient-only item.

### F2.1 — segment schema and compile (0.6.43)

Named design segments, per-environment ordered segment lists, name registry,
the exact-once rule, flag refusals inside fixed segments, unknown-environment
refusals, `None` → environment 0.

**Gates.** Old flat calls assemble **bitwise-identical** stacks (`assert_eq`
on `DesignStack` films **and spans**, post-F0.1). K=1 segmented ≡ flat. Each
refusal names the environment and the segment. Both fingerprints unmoved.

**AMENDED (A7) — the baseline needs a harness that does not exist.** Revision
1 said "Baseline recorded here: pre-segment eval timings from
`validation/benches/`". Inventory of that tree:

- `synthesis/` — **`bench_refold.py` only**: a refold-vs-static comparison
  with a simulate+fold vs LM cost breakdown. The closest proxy, and not a
  plain eval loop.
- `smatrix/` — `bench_backside_speed.py`, `bench_core_engine_scaling.py`.
- `structure/` — `bench_grid_assert.py`.

Nothing times assemble → simulate → merit of the design pipeline in
isolation, so as things stand there is nothing to record. **This item adds
`validation/benches/synthesis/bench_eval.py`** — one fixed design problem, a
timed assemble+simulate+merit loop, release-gated via
`_bench_common.require_release`
([_bench_common.py:75](validation/benches/_bench_common.py:75)) the way
`bench_refold.py` does at `:24`. Its output is the F2.2 baseline. This is the
plan's own logic ("a baseline taken after the change is not a baseline")
applied one step earlier: the baseline's *harness* must exist before the
baseline.

**R6 — record it as an artifact, not as prose.** The convention already
exists (`validation/benches/smatrix/results/*.json`). Write
`validation/benches/synthesis/results/eval_baseline.json` with the commit
hash and date, so F2.2's gate measures against a file rather than a number
someone typed into a plan.

### F2.2 — K assemblies, K solves, joint merit (0.6.44)

`residuals_multi(&[SimCurves])` routing each demand to `sims[d.env_idx]`;
per-`(env, key)` missing-curve penalties; residuals concatenate env-major,
insertion-minor. `MeritSpec::merit` and `residuals`
([merit.rs:736](rust/navette/src/smatrix/synthesis/merit.rs:736),
[merit.rs:755](rust/navette/src/smatrix/synthesis/merit.rs:755)) keep their
signatures; the multi variant is a thin wrapper, and the single-environment
call delegates through it with a one-element slice.

**No new merit formula, no new `NREQ_*`, no fold-arm change.** Any diff in
`residuals_into` means the item is over-scoped.

**The one branch.** The driver checks `K == 1` vs `K > 1` exactly once per
eval, outside every loop. K=1 calls the existing functions with existing
signatures, same call sequence op-for-op.

**AMENDED (A8) — the `merit.rs` citations are right and are not the change
surface.** The signatures at `merit.rs:736/755` are correctly quoted and
correctly kept, but the single-environment call sites that the multi variant
must route through live elsewhere, and sizing this item from its `merit.rs`
citation alone underestimates it:

- [evaluator.rs:282–285](rust/navette/src/smatrix/synthesis/evaluator.rs:282)
  — `DesignContext::evaluate_merit` calls `self.spec.merit(&sim, 1e6)`.
- [evaluator.rs:327-337](rust/navette/src/smatrix/synthesis/evaluator.rs:327)
  — the thickness-LM residual closure clones **one** base stack, sets
  thicknesses, simulates, calls `spec.residuals(&sim, out)`.
- [pipeline.rs:63–97](rust/navette/src/smatrix/synthesis/pipeline.rs:63) —
  `NeedlePipeline` owns exactly one `DesignStack` and one context;
  `run()` ([:122](rust/navette/src/smatrix/synthesis/pipeline.rs:122)) is the
  macro-loop, and it calls `clamp_all` three times.

K>1 therefore means: `NeedlePipeline` holds the shared design object plus
per-environment surroundings and rebuilds K stacks at run start and after
each insertion (multi-env §4.6's assembly-time rule); the LM residual closure
becomes span-aware across K stacks; `evaluate_merit` becomes the routing point
for `residuals_multi`. None of it touches `merit.rs` — the "no merit formula
change" contract holds. The "Where" is `pipeline.rs` (fields + run loop),
`evaluator.rs` (`DesignContext` impl + LM closure), `driver.rs` (K assembly +
the one K-branch). No re-sequencing; this is *why* F2.3 is effort-L.
Re-check the multi-env plan's "~350 Rust lines" estimate against this surface
before starting.

**Gates (§4.6, all three, as hard gates).**

- (a) **Bitwise twin:** K=1-segmented vs pre-segment binary, same seeds.
- (b) **Bench gate:** K=1-segmented vs pre-segment on `bench_eval.py`,
  threshold **< 1%** against the F2.1 baseline artifact (R6). Excess means
  new code leaked into the hot path.
- (c) **Scaling check:** K=2 with identical surroundings ≡ K=1 merit
  *exactly*; K=2 wall-time ≈ 2× single-solve ± concatenation noise.
- **Hand oracle:** one film behind two different cover sequences vs a numpy
  transfer-matrix computation, 1e-12.
- **Phase-A interaction (§1.4 + A4):** a graded film in a fixed segment
  expands with its full profile and is pinned — which requires the two A4
  plumbing lines to be present, so this twin is also their regression test.
  Twin: same stack expressed as a background-pinned flat design → bit-equal
  rows.

### F2.3 — needle and LM joint (0.6.45)

The hard item. Scan sites built per environment over the full stack, filtered
to design films; candidate locus translated back to (segment, intra-segment
position) so one insertion edits the single shared design object and
propagates everywhere by construction. Fold deposits route **by design-film
name** into shared buckets and sum across environments. Surroundings never
deposit — they are inadmissible films, which the existing flag mechanism
already handles.

**Post-Phase-A addition not in the source plan:** the locus translation is
span-aware. A needle host is a **singleton-bulk** span (F0.1 refuses the rest
— note the N1 correction: singleton-*bulk*, so an interface-carrying film is
still a legal host), so the translation maps span → span, and the alignment
assert becomes "the same design object produces the same span layout in every
environment" — a stronger and cheaper check than comparing row indices.

**Where (A8).** `pipeline.rs`, `evaluator.rs`, `driver.rs`, as F2.2.

**Gates.** K=1 fold **bitwise** vs pre-segment (against the R1 needle pin).
K=2-identical fold ≡ K=1. Finite-difference check with differing
surroundings: joint gradient == sum of per-environment analytic gradients,
1e-12. A two-environment U1 run improves both environments monotonically (the
anti-§3.2a proof — sequential per-env runs oscillate; this one must not).
Insertion lands in the design segment in every environment assembly:
positions differ, names match. An interface-carrying design film is still an
admissible host in every environment (the N1 twin, at K>1).

### F2.4 — Python surface and program sections (0.6.46)

`run_needle(layers | design={...}, ..., environments=[...])` at
[pipeline.py:188](src/navette/synthesis/pipeline.py:188); old flat-films plus
`ambient`/`substrate` kwargs keep working through a singleton default
environment on the bitwise path. Demand `environment=` tags pass through
`_dump`; the native `__post_init__` refuses unknown names. Program schema
gains `design:` and `environments:` sections, file-first refs resolved in the
live context like materials and groups.

**AMENDED (N8) — `layers` is a required positional.** `run_needle(layers,
targets, angles_deg, wavelengths, contrast, ...)` takes `layers` as the first
positional parameter, not a keyword with a default. `layers | design={...}`
is therefore a **signature change** (`layers` becomes optional, with a
refusal when neither or both are supplied), not a kwarg addition. Small, but
it is the difference between "add a parameter" and "change a public
signature", and the second needs the `.pyi` and a both/neither refusal test.

**AMENDED (A1, corrected) — the program envelope needs the §1.2 treatment,
but not for the reason the amendment gave.**

The amendment is right that a second equality gate exists and right that this
item must deal with it. It is wrong about the mechanism, and the difference
changes what the gates test.

*What is actually there* (note the path — it is `rust/navette/src/config.rs`,
not `smatrix/synthesis/config.rs`):

- [config.rs:26](rust/navette/src/config.rs:26) —
  `pub const PROGRAM_SCHEMA_VERSION: u32 = 1;`
- [config.rs:283](rust/navette/src/config.rs:283) —
  `Some(v) if v == PROGRAM_SCHEMA_VERSION as u64 => {}`, else
  `"program schema_version {other:?} unsupported"`. Same shape as
  `version.rs:16`.
- [config.rs:297-307](rust/navette/src/config.rs:297) — unknown **top-level**
  keys refuse. Pinned at [config.rs:722-743](rust/navette/src/config.rs:722).

*Where the amendment goes wrong (N9).* It reasons that `design:` and
`environments:` "are new top-level sections in a program document. Under the
v1 gate they refuse as unknown keys, so F2.4 *necessarily* bumps
`PROGRAM_SCHEMA_VERSION`." But for `kind == "program"` the only permitted
top-level keys are `schema_version`, `kind`, `name` and `sections`
([config.rs:295](rust/navette/src/config.rs:295)) — every payload, including
`materials`, `groups`, `structures` and `architect`, lives *inside* the
`sections` mapping. And `load_program_parts`
([config.rs:580–630](rust/navette/src/config.rs:580)) reads that mapping with
bare `sections.get("materials")`-style lookups: **there is no whitelist of
section names, and an unknown section is silently ignored.**

So a v1 binary handed a v2 program does not refuse `environments:` — it drops
it and optimizes a single-environment stack that looks entirely plausible.
Same conclusion as the amendment (bump, with a range), reached from the
opposite failure: **silent drop, not refusal.** Which is precisely §1.2's
failure class, one gate over — the amendment had the right instinct and the
wrong evidence.

*Two consequences the corrected reading adds:*

1. **State the shape explicitly.** If `design:`/`environments:` are placed at
   *top level* instead of inside `sections`, they refuse under the v1 gate
   and the amendment's reasoning becomes true. Either shape is defensible;
   inside `sections` matches every existing payload and is the recommendation.
   Write it down, because the two shapes have opposite failure modes.
2. **A section-name whitelist is worth adding in the same commit.** The
   silent-ignore behaviour is the bug that makes the bump necessary; leaving
   it in place means the next section added has the same problem. Collect
   unrecognised section names and refuse, the way top-level keys already are.

**Amendment (adopted, with the above).** F2.4 gains the F1.4 treatment in
both homes:

- Split the constant: `PROGRAM_SCHEMA_VERSION = 2` +
  `MIN_READABLE_PROGRAM_SCHEMA_VERSION = 1`; `gate_document` accepts the
  range and refuses below (stale) and above (newer build) with distinct
  messages. Both homes move together:
  [config.rs:26](rust/navette/src/config.rs:26) and
  [program.py:45](src/navette/config/program.py:45). Note
  `program.py`'s module docstring already claims the gate refuses "future"
  versions — it does, but with the wrong message; fix the wording with the
  constant.
- `validation/smoke/test_request_bits.py` — reuse F1.4's `_accepted_range`
  helper (A3, point 4); the program test asserts `(1, 2)`. **Same file as
  F1.4**: rebase carefully.
- The existing refusal test at
  [config.rs:726-731](rust/navette/src/config.rs:726)
  (`schema_version: 2` refuses) **flips meaning** after the bump: a v2
  program must be accepted and a v3 refused. Update it deliberately, not as
  drive-by fallout.
- `validation/regression/config/test_program.py` extends here, not at F1.4
  (A3).
- Cross-plan interaction row 4 gains: "the program envelope bump repeats the
  state-gate pattern — same range, same two-sided messages, same sync
  helper — but for the opposite unknown-key policy."

**Gates.** Surface ≡ JSON, hex-compared. Refusals at construction,
`ValueError`, named, ASCII. A v1 program fixture loads and assembles
bit-identically. A `schema_version: 3` program refuses naming the newer
build; `schema_version: 0` refuses as stale. Both endpoints asserted in
`test_request_bits.py`. An unknown section name refuses, naming it.
`run_needle` refuses both-`layers`-and-`design` and neither.
`check_pyi_sync.py` and `check_exposure.py` clean.

**Risk.** M (was S). **Effort.** L (was M). It carries a schema gate, a
public signature change and a new refusal class.

---

## 5. Phase C — documentation and release

### F3.1 — docs, examples, audit, release (0.6.47)

- `docs/spectralweave-target-kinds.md`: environment section, coverage matrix
  row, the ×K cost note.
- Worked example: the same AR stack bare and laminated, optimized jointly.
- Worked example: a rugate filter as a `FixedSpan` gradient, and the same
  filter as a discrete stack, with the row-count and merit comparison.
- Névot-Croce validity caveat cross-referenced from the gradient docs —
  gradient sublayers make interfaces numerous, and rtype 5 on each of 64
  sublayer boundaries is a fast way to leave the model's validity band.
- **Document the thin-gradient floor (A2/F0.1).** `clamp_min_nm` interacts
  with sublayer thickness: a gradient thinner than ~6 nm at default settings
  has sublayers at the clamp boundary. F0.1 stops them being deleted; the
  docs must still say that authoring a 3 nm gradient gives three 1 nm
  sublayers and that this is near the staircase-validity floor anyway.
- `tools/check_exposure.py` re-audit: new public functions either exported or
  allowlisted with a reason.
- Full battery, CI green, `main` fast-forward, tag.

---

## 6. Cross-plan interactions (the reason these are one plan)

| # | Interaction | Handled by |
|---|---|---|
| 1 | Multi-environment routing must address a *span*, not a row, once gradients exist | F0.1 first; F2.3 built span-aware from the start |
| 2 | Gradient sublayers are thin by construction; cleanup deletes thin films — **and `clamp_all` already deletes pinned ones today** | F0.1's clamp exemption (A2); the `optimize` filter covers `remove_thin_layers` only |
| 3 | A gradient in a fixed environment segment must keep its profile | Nearly free — the background rule does the physics; two plumbing lines at F1.1 (A4); pinned in F2.2 |
| 4 | Two schema surfaces bumping at once would be untraceable | State schema v2 at F1.4; program envelope v2 at F2.4; never in the same release. Both become ranges, with the same two-sided messages and the same `_accepted_range` helper — but for **opposite** unknown-key policies: state ignores, program-sections ignore, program top-level refuses (N9) |
| 5 | ×K solve cost gets worse when surroundings are graded (many rows) | Documented; the S-matrix embedding follow-up (§6.6) becomes more valuable, still not v1 |
| 6 | Needle refusal inside a span, and needle restriction to design segments, are the same predicate | One admissibility function, **three** callers: `build_scan_sites`, `insert_needle_seed` (R3), and F2.3's locus translation |
| 7 | **(new, N2)** A saturated `RateCapped` tail is bit-identical rows, which `merge_adjacent` collapses | F0.1 handles intra- and cross-span merges; F1.2 owns the end-to-end twin |
| 8 | **(new, N1)** An interface slice shares its carrier's span, so "non-singleton" is not "profiled" | The predicate is singleton-**bulk** everywhere: F0.1, F0.2, F1.6, F2.3 |
| 9 | **(new, U1)** Clamping a layer up to the minimum without raising the LM lower bound is a limit cycle that the stagnation detector reports as oscillation | F0.3 moves `lb` with the policy, and refuses `ClampUpAlways` alongside needle insertion — the two cannot both be on |
| 10 | **(new, U4)** Re-slicing a graded span is a construction operation; doing it in an inner loop is both 1000× the cost and a step discontinuity in the merit surface | F1.7 fixes three refresh points and pins the count with a counter twin; F1.6's frozen fractions and frozen nk are what make the LM solve well-posed |
| 11 | **(new, U2/U1)** "Too thin" means two different things once a span can be scaled: remove the layer, or scale it to the floor | F0.2 removes whole and reports; F0.3 adds the clamp-up policy but defers its span branch; F1.6 wires the span branch in and inverts F0.3's deferral twin |
| 12 | **(new, U2)** Marking a graded film `optimize = true` used to mean "homogenize me"; F1.6 gives it a second meaning | The flag-meaning change is in F1.6's licence, not implicit; F1.1's pinned path (`optimize = false, needle = false`) is untouched and still means "keep the profile, do not touch it" |

---

## 7. Risks

| Risk | Mitigation |
|---|---|
| Span bookkeeping drifts across five mutators | One named `assert_spans_partition` (R2) after every mutation and both constructors; fingerprints |
| **F0.1's span predicate silently re-classifies interface films and moves the needle fingerprint** | Singleton-**bulk**, not singleton (N1); the interface-film twin exists to fail if that is ever loosened |
| **`clamp_all` deletes a pinned thin graded span — today, not hypothetically** | A2 exemption in the removal branch; the 5 nm clamp twin asserts film count and total thickness |
| **A saturated `RateCapped` tail merges and the span shrinks under the bookkeeping** | N2: intra-span merge handled at F0.1, twinned end-to-end at F1.2, spectra asserted bitwise equal across the merge |
| `inflate_design` re-thicknesses a pinned profile per sublayer | N3: same exemption on its selection list; inflate twin (F0.2) |
| **The floor deletes a graded layer and the ceiling ignores it — both live, both measured on 0.6.32** | F0.2 makes both span quantities, with a six-item licence and a twin per item; `clamp_all`'s discarded return value becomes a report |
| **F0.2 is the one item allowed to change a number, so "it was already like that" stops being available as an excuse** | The licence is enumerated before the code; a seventh difference is a defect, not a discovery |
| **Clamp-up without a matching LM lower bound stalls the run and blames stagnation** | F0.3: `lb` moves with the policy; the monotone-merit twin is the positive statement that the limit cycle is absent |
| **`ClampUpAlways` silently disables layer elimination, so needle runs only ever grow** | Refused at `NeedlePipeline::new` when `needles_per_cycle > 0`, naming both settings and pointing at `ClampUpFinal` |
| **A refresh call appears in a fourth place because someone fixes a staleness symptom where it showed up** | F1.7's refresh-count twin asserts exactly five calls in a three-cycle run; a stray call moves it by orders of magnitude |
| **A mid-solve `ceil()` boundary changes the row count and puts a step in the merit surface** | F1.7 freezes the row count for the whole LM solve; the residual-length twin crosses a boundary deliberately |
| **F1.6's fraction weights are wrong and the Jacobian is plausibly wrong rather than obviously wrong** | Two oracles: central FD to 1e-6, and a *bitwise* comparison against the weighted sum of the row columns |
| `Fixed` inhomogeneous mode stops being bitwise | The existing randomized differential suite (`test_differential.py:95`, `:119`) is the gate, not review |
| Sublayer-count float boundaries diverge | Differential pin over randomized thicknesses (`sub_layer_count` precedent); all four callers see the same `delta_layer` |
| EMA per-draw cost surprises | Bench before/after; memo on `(nk_a, nk_b, f)` per draw held in reserve |
| Schema gate change breaks state loading | v1 fixture committed first (R5); both out-of-range directions tested |
| **The schema range breaks the sync helper and inverts two tests** | A3: `_accepted_range` replaces `_accepted_version`; `test_stale_schema_versions_refused` and the `config.rs` refusal test updated deliberately |
| **`Layer` has no derive-based serde, so "add a serde default" is not available** | N4: explicit `map.get` reads; the `serialize_map(13)` count moves; emit-only-when-`Some` decided at §8.10 |
| **A v1 binary silently drops a v2 program's `environments:`** | N9/A1: program range gate in both homes, plus a section-name whitelist in the same commit |
| The K>1 branch leaks into the hot path | Bench gate < 1% against the F2.1 **artifact** (A7 + R6), not a prose number |
| Fold routing across environments silently mis-sums | Alignment assert on span layout; FD check against summed per-env analytics |
| Users conflate `inh_delta` with `gradient` | Refusal messages cross-reference; docs §D0 |
| CI green locally but red remotely | `tools/check_toolchain.py` before every push; wait for CI before ff-ing `main` |
| **The plan's own claims drift from reality mid-series** | R7: every CORRECTIONS block names the A/N IDs it adopted |

---

## 8. Open decisions

Resolve at the item that first needs them; each is a one-liner.

1. **Default environment name** — explicit `"default"` vs `None`-only.
   Recommend named; explicit beats magic. *(F2.1)*
2. **Environment list home** — call-only, or also `TargetSet.environments`
   with JSON winning on clash. Recommend both. *(F2.1)*
3. **Per-environment `missing_penalty` scale** — recommend one floor for all.
   *(F2.2)*
4. **Exact-once relaxation** — allow an environment to omit a design segment?
   Recommend refuse in v1; omission is almost certainly a typo'd environment,
   and silent non-contribution would lie about jointness. Relax later behind
   a loud `contributes: false` marker if a real case appears. *(F2.1)*
5. **Fixed-film user names** — auto-only, or non-colliding custom names for
   debugging. Recommend auto-only v1; names that address nothing invite
   `contrast` entries that silently never split. *(F2.1)*
6. **Per-layer EMA host selection** — deferred by §D9.4 until a real
   asymmetric-mixture case lands. Do not design it now. *(none)*
7. **`ProfileShape` beyond `Linear`** — S-curve and exponential arrive as new
   enum variants, no schema break, each with its own twin batch. *(none)*
8. **(A5) `material_b`'s nk on the design path** — recommend `ArrayFilm` /
   `LayerRow` carry `gradient` with `material_b` **and** `nk_b`, evaluated
   Python-side in `_film_dicts`, injected into the provider entries by
   `assemble_stack`. Refuse when absent; never fall back to `material_a`.
   *(F1.1, mirrored at F1.5)*
9. **(A6, corrected) EMA selector type** — recommend `GradientSpec.ema` is
   the existing `MixRule`
   ([materials/mod.rs:54](rust/navette/src/materials/mod.rs:54)), not a new
   four-variant `EmaModel`. All six in-tree kernels work immediately, the
   parameterized ones carry their parameters in the variant, and the crate
   keeps one mixing vocabulary. Cost: add `Serialize`/`Deserialize`/
   `PartialEq` derives and name the serde representation deliberately, since
   it becomes schema-visible. *(F1.1)*
10. **(N4) `gradient` in the serialized state — always, or only when
    `Some`?** Recommend only-when-`Some`: a stack with no gradient then
    serialises byte-identically at v1 and v2 apart from the version tag,
    which makes "v2 changes nothing unless you use the feature" literally
    true. Cost: `test_state_fingerprint` needs two key lists. *(F1.4)*
11. **(N9) Program sections: inside `sections`, or top-level?** Recommend
    inside `sections`, matching every existing payload — and add the
    section-name whitelist in the same commit, because the silent-ignore
    behaviour is what makes the version bump necessary in the first place.
    *(F2.4)*
12. **(U2) How a span declares itself scalable.** Recommend reusing
    `optimize = true` on the graded carrier rather than adding a
    `scalable` flag: it is the flag that already means "the optimizer may
    move this layer's thickness", and a second flag would create four
    combinations where two are meaningful. Cost: the combination changes
    meaning (it used to trigger homogenization), which is why it sits in
    F1.6's licence. The alternative — an explicit `GradientSpec.scalable`
    — is more discoverable and adds a schema field; decide before F1.4,
    because F1.4 is the release that freezes the Phase A field set.
    *(F1.6, frozen by F1.4)*
13. **(U4) Should the per-cycle refresh be configurable?** Recommend no in
    v1: three fixed refresh points, no knob. A `refresh_profiles_each_cycle`
    toggle would let a user turn off the reconciliation that bounds the
    staleness, and the only thing they would buy is one re-expansion per
    macro cycle — which F1.7 already argues is negligible next to the
    seven `optimize_thicknesses` calls in the same cycle. Revisit only if a
    profiling run says otherwise. *(F1.7)*
14. **(U1) Does `clamp_min_nm` get renamed?** Recommend no. It genuinely
    does two jobs — elimination threshold under `Remove`, manufacturing
    floor under the clamp-up policies — but renaming a public config key
    breaks every caller to buy clarity that a doc comment can supply. The
    doc comment must state both jobs and which policy selects which.
    *(F0.3)*

---

## 9. Progress log

Updated as items land. Format: version, commit, what moved, what did not, and
which audit IDs the item's CORRECTIONS block adopted (R7).

| Version | Item | Commit | Amendments adopted | Status |
|---|---|---|---|---|
| 0.6.33 | F0.1 | — | A9.3, R1, R2, R3, N1, N2, **U5** (half) | not started |
| 0.6.34 | F0.2 | — | A2, N3, **U5** (half) | not started |
| 0.6.35 | F0.3 | — | **U1** | not started |
| 0.6.36 | F1.1 | — | A4, A5, A6, A9.1, R4, N5, N7 | not started |
| 0.6.37 | F1.2 | — | A9.2, N2 | not started |
| 0.6.38 | F1.3 | — | — | not started |
| 0.6.39 | F1.6 | — | **U2** | not started |
| 0.6.40 | F1.7 | — | **U3, U4** | not started |
| 0.6.41 | F1.4 | — | A3, R5, N4 | not started |
| 0.6.42 | F1.5 | — | A5 (mirror), N9 | not started |
| 0.6.43 | F2.1 | — | A7, R6 | not started |
| 0.6.44 | F2.2 | — | A4, A8 | not started |
| 0.6.45 | F2.3 | — | A8, N1 | not started |
| 0.6.46 | F2.4 | — | A1 (corrected), A3, N8, N9 | not started |
| 0.6.47 | F3.1 | — | A2 (docs), **U1/U2/U3** (docs) | not started |
