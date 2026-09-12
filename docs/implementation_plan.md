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

---

## 0. How to use this plan

### 0.1 Ground rules

1. **One item = one feature = one `0.0.1` version bump.** Seven version
   sites in five files move together
   (`pyproject.toml`, `Cargo.toml` ×2 hits, `src/navette/__about__.py`,
   `Cargo.lock` ×2 hits, `README.md`). Then: CHANGELOG entry, item marked
   **DONE** here with a CORRECTIONS block if reality differed from the
   design, master-table row struck, commit, push to `dev_feature`.
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
6. **Behaviour changes are announced, never silent.** Corrections at
   assembly (`from_design`) are warnings with the offending name in them;
   contradictions are refusals naming both sides and the alternative.
7. **Python-surface error messages stay ASCII** (cp1252 consoles).
   Clippy `-D warnings` rejects `///` doc comments on statements.

### 0.2 Priority / risk / effort

| Field | Meaning |
|---|---|
| **Priority** | P0 = blocks a later item; P1 = user-visible feature; P2 = plumbing the features need; P3 = docs |
| **Risk** | probability the change itself breaks something already shipped |
| **Effort** | S < 1 h · M = hours · L = 1–2 days · XL = multi-day |

### 0.3 Master item table

| ID | Item | Version | Priority | Risk | Effort | Source |
|---|---|---|---|---|---|---|
| F0.1 | Span provenance on `DesignStack` — protection by span, not by flag coincidence | 0.6.33 | **P0** | M (touches cleanup + needle host selection) | M | §D4.1–2, corrected §2 |
| F1.1 | Gradient data model + `FixedSpan` expansion + homogenize path | 0.6.34 | P1 | M (new expansion branch) | L | §D2–D3, §D4.3 |
| F1.2 | Gradient `RateCapped` mode — thickness-relative slope with caps | 0.6.35 | P1 | S | M | §D0(b), §D2 |
| F1.3 | `InhMode::RateCapped` — thickness-relative single-material drift | 0.6.36 | P1 | M (legacy path must stay bitwise) | M | §D0(b), §D2 |
| F1.4 | Schema v2 + a readable-version **range**, not a point | 0.6.37 | **P0** | M (every state file reads through this gate) | S–M | §D5 + correction §1.2 |
| F1.5 | `design_config` rows + Python `Layer.gradient` surface | 0.6.38 | P1 | S | M | §D5 |
| F2.1 | Environment segment schema + compile/validation | 0.6.39 | P2 | S | M | §4.1–4.2 |
| F2.2 | K assemblies, K solves, `residuals_multi` — joint merit | 0.6.40 | P1 | M (driver loop) | L | §4.3, §4.6 |
| F2.3 | Needle + LM joint: locus translation, name-routed fold sum | 0.6.41 | P1 | **L** (the hard one — fold routing) | L | §4.4 |
| F2.4 | Python `environments=` / `design=` surface + program sections | 0.6.42 | P1 | S | M | §4.5 |
| F3.1 | Docs, worked examples, exposure re-audit, release | 0.6.43 | P3 | S | M | §7-S5, §D7 |

Eleven items, `0.6.33 → 0.6.43`. Tag `v0.7.0` at F3.1 if the minor marker is
wanted; the version number itself stays on the `0.0.1` ladder throughout.

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

**Resolution (F1.4, promoted to P0 and sequenced before the Python surface):**
split the single constant into a readable **range**.

- `SCHEMA_VERSION = 2` — what this code writes.
- `MIN_READABLE_SCHEMA_VERSION = 1` — the oldest it accepts.
- `check_schema_version` accepts `MIN_READABLE ..= SCHEMA_VERSION`, refuses
  below (stale) and above (written by a newer build) with distinct messages.
  The above-range message is the important one: it is what stops a v1 binary
  from silently dropping a gradient.
- Fields added in v2 carry serde defaults that reproduce v1 behaviour exactly
  (`gradient: None`, `inh_mode: Fixed`, `shape: Linear`).
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
already computed; F0.1 is plumbing, not derivation. That is why it is cheap
enough to front-load.

### 1.4 The `optimize/needle = false` rule already routes gradients correctly

`§4.2` (multi-env) forces `optimize = false, needle = false` on every
fixed-segment film. `assemble_stack`
([driver.rs:130](rust/navette/src/smatrix/synthesis/driver.rs:130)) computes
the background set as `inhomogen && !optimize && !needle`, and `from_design`
expands background-pinned graded films **with** their profile.

Consequence, free: **a graded or gradient film inside a fixed environment
segment is automatically background-expanded and keeps its full profile.** A
graded cover glass in a surrounding stack does the physically right thing
without one line of new code. Stated here so nobody adds a special case for
it, and pinned by a twin in F2.2.

### 1.5 The two "unprotected span" holes in §D4.1–2 are already closed

The gradient plan presents thin-removal and needle insertion as live hazards
for profile spans. Against this tree they are not: `remove_thin_layers`
filters candidates on `l.optimize`, and needle host admissibility is
`layer.needle`. Background-pinned spans carry both as `false`, so both
operations already skip them.

The plan is not wrong about the danger, only about its timing — the
protection is a consequence of the span being wholly frozen, not a rule about
spans. F0.1 is re-scoped accordingly (§2): same code changes, honest
justification, and the one twin that actually distinguishes before from after
is an *optimizable* thin-sublayer span, not a pinned one.

---

## 2. Phase 0 — foundations

### F0.1 — span provenance on `DesignStack` (0.6.33)

**Not a bug today — a coupling.** Two pipeline operations would corrupt a
profile if they ever reached one:

1. `remove_thin_layers` ([cleanup.rs:40](rust/navette/src/smatrix/synthesis/cleanup.rs:40))
   deletes films below a thickness floor. Sublayers of a profile are thin *by
   construction* — a 200 nm graded film at 12 sublayers is twelve 16.7 nm
   rows.
2. `insert_needle_seed` ([structure.rs:296](rust/navette/src/smatrix/synthesis/structure.rs:296))
   splits a film in two and inserts a seed between the halves. Done inside a
   profile span, that orphans the rows on either side.

**Both are already safe, and it is worth being exact about why**, because the
source plan (§D4.1–2) presents them as open holes and they are not.
`remove_thin_layers` selects candidates with
`.filter(|(_, l)| l.optimize && l.d_nm < threshold)` — its own comment names
pinned graded sublayers as never-candidates. Needle host admissibility is
`layer.needle` ([needle_pass.rs:62](rust/navette/src/smatrix/synthesis/needle_pass.rs:62)).
Background-pinned graded spans carry `optimize = false, needle = false` on
every row, so both operations skip them today.

**What is actually wrong is the coupling.** A span is protected only as a
side effect of being *entirely frozen*. The protection is not a rule about
spans; it is a rule about flags that happens to cover spans. Two things break
it:

- **The obvious next feature.** A gradient whose *total* thickness is a free
  parameter — the first thing anyone asks for with a rugate filter — needs
  `optimize = true` on its rows. The moment that exists, cleanup eats the
  sublayers one at a time, each removal looking locally like an improvement.
- **F2.3.** Multi-environment fold routing maps a shared design parameter to
  its position in each of K assemblies. With profiles, a parameter owns a
  *span*, not a row. A flat `Vec<LayerSpec>` cannot express that, and nothing
  in the current structure can reconstruct it after `from_design` returns.

So this item makes the provenance explicit and the protection a stated rule
rather than a coincidence. It ships no user-visible behaviour change.

**Change.**

- `DesignStack` gains `spans: Vec<Span>` alongside `films`, populated in
  `from_design` from the vector it already receives and currently discards.
  Non-`from_design` constructors (`with_films`, used by tests) synthesise
  singleton spans — one row, one span — so the invariant "every row belongs
  to exactly one span" holds everywhere.
- Every mutation that changes row count maintains it:
  `insert_needle_seed`, `merge_adjacent`, `remove_film`, and
  `remove_thin_layers`. This is the actual work of the item; the data model
  is trivial and the bookkeeping is not.
- `remove_thin_layers` exempts rows whose span is non-singleton — as a rule
  about spans, *in addition to* the existing `optimize` filter, not
  replacing it. Same outcome today, different reason, and it survives a
  future optimizable span.
- Needle host selection refuses non-singleton spans, naming the span's
  material and row range — likewise alongside the `needle` flag.

**Gates.**

- Both fingerprints **unmoved**. No stack without a profile has a
  non-singleton span, so no existing run can take a different path.
- A `debug_assert` after every mutation: spans partition `0..films.len()`
  contiguously, in order, no gaps, no overlap. The assert is the proof that
  the bookkeeping is right — reasoning about four mutators is not.
- Thin-removal twin: a graded film whose sublayers are below the floor **and
  whose rows are marked `optimize = true`** → all sublayers survive. This is
  the twin that matters: it is the case the flag filter does not cover, and
  it fails before this item and passes after.
- Needle-refusal twin: message names material + row range.
- Existing `merge_adjacent` tests unchanged (merge is nk-keyed; distinct
  sublayer nk never collapse — verified by test, not by reasoning, per §D4).

**Risk.** M. Touches the needle host path, which is the most-exercised code
in the synthesis pipeline. Mitigated by the fingerprint gate: if either
fingerprint moves, the change is wrong, full stop.

---

## 3. Phase A — gradient-index and inhomogeneous layers

Physics reference: `docs/plans/gradient_layers_plan.md` §D0–D4, §D9. That
document's vocabulary is locked and not restated here. What follows is the
sequencing and the gates.

### F1.1 — data model + `FixedSpan` expansion + homogenize path (0.6.34)

**Scope.** `GradientSpec`, `EmaModel`, `GradientMode::FixedSpan`,
`ProfileShape::Linear`, validation, the expansion branch, and the
non-background homogenize path (§D4.3) so the feature is safe from the first
commit.

**Where.**

- `structure/layer.rs`: `Layer` gains `gradient: Option<GradientSpec>`,
  default `None`. New module `structure/gradient.rs` for the spec types —
  `layer.rs` is 606 lines and already carries the validation surface.
- `structure/expansion.rs`: new branch beside the `inhomogen` branch at
  [expansion.rs:258](rust/navette/src/structure/expansion.rs:258), mutually
  exclusive with it by validation.
- `structure/validation.rs`: gradient issues join `property_issues`, so the
  layer-construction gate from R3.5 covers them at every door — Python
  constructor, four setters, `set_properties`, `from_state`, and
  `assemble_stack`. **No new gate, no second rule.**
- `materials/ema.rs`: reused unchanged. No new mixing math. Gradient calls
  `bruggeman` / `maxwell_garnett` / `looyenga` / `lichtenecker` per sublayer.

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

**Gates.**

- Both fingerprints **unmoved** (gradient absent ⇒ `None` ⇒ untouched path).
- Per-row nk vs a direct `ema_*` oracle on the same `(nk_a, nk_b, f_i)` —
  **bitwise**, same kernel, same inputs.
- `f_i` recomputed from the mode formula, EMA outputs compared to 1e-15.
- Thickness-independence: same `FixedSpan` spec at 50 nm and 500 nm →
  endpoint rows bitwise equal, counts differ.
- Sum of sublayer thicknesses `== thickness` exactly (last sublayer absorbs
  the float remainder).
- Inversion reverses row order; profile direction flips; tested both ways.
- Roughness lands on the first sublayer only (legacy convention).
- Refusals, each naming both sides: `gradient` + `inhomogen` together;
  `material_a == material_b`; a material absent from the provider;
  `f_*` outside `[0, 1]`.
- Homogenize warning names the mixture and `f_mid`.
- `check_exposure.py` / `check_pyi_sync.py` clean.

**Docstring obligation (§D9.4).** `GradientSpec` must carry, at the call
site, the statement that `material_a` is **always the host** and `f` is the
volume fraction of B in A — Maxwell-Garnett and Mori-Tanaka are
host/inclusion-asymmetric, so A↔B is not the same physics and a silent swap
corrupts fits. This is a docstring requirement, and it is checked in review,
not by a tool.

### F1.2 — gradient `RateCapped` mode (0.6.35)

`f(z) = f_start + rate * (z / ref_thickness)`, clamped to `[f_min, f_max]`.
`ref_thickness` default 100 nm, caps default `[0.0, 1.0]`, `rate` finite with
sign free.

**Gates.** 200 nm film, `f_start = 0.1`, `rate = 0.3` per 100 nm → ends at
0.7. Same spec at 400 nm → saturates at 1.0 from z = 300 nm, and the tail
rows are **bitwise equal to pure `material_b`** — that is the saturation
proof, not an approximate comparison. Negative rate saturates to pure
`material_a` at the other end. Refusal when two slope-like keys are supplied
for one layer (`"gradient: 'rate' and 'delta' are the same slope - give
one."`, §D0) — no silent precedence.

### F1.3 — `InhMode::RateCapped` for single-material drift (0.6.36)

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
thickness, so there is no circularity; evaluated once.

**Gates.** `Fixed` mode **bitwise** against the legacy oracle over the
existing randomized differential suite — this is the whole risk of the item
and the only acceptable evidence. `RateCapped` delta against hand-computed
`min(rate*t/t_ref, cap)` over a thickness sweep, with the saturation knee
exact. Double the thickness below the cap → double the delta.

**Note.** `inhomogen` is **not deprecated**. Single-material drift is real
physics (oxidation gradients, nitrides) and coexists with gradients by the
§D2 rule.

### F1.4 — schema v2 and a readable version range (0.6.37)

The correction from §1.2 above. **Sequenced before the Python surface**
because every state file on disk reads through this gate, and getting it
wrong is a data-loss-shaped bug rather than a feature defect.

**Gates.** A v1 state fixture committed to the tree loads and expands
bit-identically. A `schema_version: 3` state refuses, with a message that
says the state was written by a newer build — not "stale". Round-trip
v2 → v2 preserves every new field. `test_roundtrip.py` and
`test_program.py` extended rather than replaced.

### F1.5 — config rows and Python surface (0.6.38)

`design_config.rs::LayerRow` gains `gradient: Option<GradientJson>` next to
the existing `inhomogen` flag at
[design_config.rs:85](rust/navette/src/smatrix/synthesis/design_config.rs:85).
`ArrayFilm` in `driver.rs` gains the same. `builders.py` passes through;
Python `Layer` gains `gradient`, validated natively — no Python pre-checks,
no pydantic.

**Gates.** Python surface ≡ JSON path, hex-compared. Refusals arrive at
construction as `ValueError` with the native message, ASCII. `.pyi` stubs
updated, `check_pyi_sync.py` clean.

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

### F2.1 — segment schema and compile (0.6.39)

Named design segments, per-environment ordered segment lists, name registry,
the exact-once rule, flag refusals inside fixed segments, unknown-environment
refusals, `None` → environment 0.

**Gates.** Old flat calls assemble **bitwise-identical** stacks (`assert_eq`
on `DesignStack` films). K=1 segmented ≡ flat. Each refusal names the
environment and the segment. Both fingerprints unmoved.

**Baseline recorded here**: pre-segment eval timings from
`validation/benches/`, which F2.2's < 1% gate measures against. Recording it
in the item *before* the one it gates is deliberate — a baseline taken after
the change is not a baseline.

### F2.2 — K assemblies, K solves, joint merit (0.6.40)

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

**Gates (§4.6, all three, as hard gates).**

- (a) **Bitwise twin:** K=1-segmented vs pre-segment binary, same seeds.
- (b) **Bench gate:** K=1-segmented vs pre-segment on the existing eval
  benchmarks, threshold **< 1%** against the F2.1 baseline. Excess means new
  code leaked into the hot path.
- (c) **Scaling check:** K=2 with identical surroundings ≡ K=1 merit
  *exactly*; K=2 wall-time ≈ 2× single-solve ± concatenation noise.
- **Hand oracle:** one film behind two different cover sequences vs a numpy
  transfer-matrix computation, 1e-12.
- **Phase-A interaction (§1.4):** a graded film in a fixed segment expands
  with its full profile and is pinned. Twin: same stack expressed as a
  background-pinned flat design → bit-equal rows.

### F2.3 — needle and LM joint (0.6.41)

The hard item. Scan sites built per environment over the full stack, filtered
to design films; candidate locus translated back to (segment, intra-segment
position) so one insertion edits the single shared design object and
propagates everywhere by construction. Fold deposits route **by design-film
name** into shared buckets and sum across environments. Surroundings never
deposit — they are inadmissible films, which the existing flag mechanism
already handles.

**Post-Phase-A addition not in the source plan:** the locus translation is
span-aware. A needle host is a singleton span (F0.1 refuses the rest), so the
translation maps span → span, and the alignment assert becomes "the same
design object produces the same span layout in every environment" — which is
a stronger and cheaper check than comparing row indices.

**Gates.** K=1 fold **bitwise** vs pre-segment. K=2-identical fold ≡ K=1.
Finite-difference check with differing surroundings: joint gradient == sum of
per-environment analytic gradients, 1e-12. A two-environment U1 run improves
both environments monotonically (the anti-§3.2a proof — sequential per-env
runs oscillate; this one must not). Insertion lands in the design segment in
every environment assembly: positions differ, names match.

### F2.4 — Python surface and program sections (0.6.42)

`run_needle(layers | design={...}, ..., environments=[...])` at
[pipeline.py:188](src/navette/synthesis/pipeline.py:188); old flat-films plus
`ambient`/`substrate` kwargs keep working through a singleton default
environment on the bitwise path. Demand `environment=` tags pass through
`_dump`; the native `__post_init__` refuses unknown names. Program schema
gains `design:` and `environments:` sections, file-first refs resolved in the
live context like materials and groups.

**Gates.** Surface ≡ JSON, hex-compared. Refusals at construction,
`ValueError`, named, ASCII. `check_pyi_sync.py` and `check_exposure.py`
clean.

---

## 5. Phase C — documentation and release

### F3.1 — docs, examples, audit, release (0.6.43)

- `docs/spectralweave-target-kinds.md`: environment section, coverage matrix
  row, the ×K cost note.
- Worked example: the same AR stack bare and laminated, optimized jointly.
- Worked example: a rugate filter as a `FixedSpan` gradient, and the same
  filter as a discrete stack, with the row-count and merit comparison.
- Névot-Croce validity caveat cross-referenced from the gradient docs —
  gradient sublayers make interfaces numerous, and rtype 5 on each of 64
  sublayer boundaries is a fast way to leave the model's validity band.
- `tools/check_exposure.py` re-audit: new public functions either exported or
  allowlisted with a reason.
- Full battery, CI green, `main` fast-forward, tag.

---

## 6. Cross-plan interactions (the reason these are one plan)

| # | Interaction | Handled by |
|---|---|---|
| 1 | Multi-environment routing must address a *span*, not a row, once gradients exist | F0.1 first; F2.3 built span-aware from the start |
| 2 | Gradient sublayers are thin by construction; cleanup deletes thin films | Already safe via the `optimize` filter; F0.1 restates it as a span rule so it survives an optimizable span |
| 3 | A gradient in a fixed environment segment must keep its profile | Free — the background rule already does it (§1.4); pinned in F2.2 |
| 4 | Two schema surfaces bumping at once would be untraceable | State schema v2 lands in F1.4; program schema gains sections in F2.4; never in the same release |
| 5 | ×K solve cost gets worse when surroundings are graded (many rows) | Documented; the S-matrix embedding follow-up (§6.6) becomes more valuable, still not v1 |
| 6 | Needle refusal inside a span, and needle restriction to design segments, are the same predicate | One admissibility function, two callers |

---

## 7. Risks

| Risk | Mitigation |
|---|---|
| Span bookkeeping drifts across four mutators | `debug_assert` partition check after every mutation; fingerprints |
| `Fixed` inhomogeneous mode stops being bitwise | The existing randomized differential suite is the gate, not review |
| Sublayer-count float boundaries diverge | Differential pin over randomized thicknesses (`sub_layer_count` precedent) |
| EMA per-draw cost surprises | Bench before/after; memo on `(nk_a, nk_b, f)` per draw held in reserve |
| Schema gate change breaks state loading | v1 fixture committed to the tree; both out-of-range directions tested |
| The K>1 branch leaks into the hot path | Bench gate < 1% with the baseline recorded one item earlier |
| Fold routing across environments silently mis-sums | Alignment assert on span layout; FD check against summed per-env analytics |
| Users conflate `inh_delta` with `gradient` | Refusal messages cross-reference; docs §D0 |
| CI green locally but red remotely | `tools/check_toolchain.py` before every push; wait for CI before ff-ing `main` |

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

---

## 9. Progress log

Updated as items land. Format: version, commit, what moved, what did not.

| Version | Item | Commit | Status |
|---|---|---|---|
| 0.6.33 | F0.1 | — | not started |
| 0.6.34 | F1.1 | — | not started |
| 0.6.35 | F1.2 | — | not started |
| 0.6.36 | F1.3 | — | not started |
| 0.6.37 | F1.4 | — | not started |
| 0.6.38 | F1.5 | — | not started |
| 0.6.39 | F2.1 | — | not started |
| 0.6.40 | F2.2 | — | not started |
| 0.6.41 | F2.3 | — | not started |
| 0.6.42 | F2.4 | — | not started |
| 0.6.43 | F3.1 | — | not started |
