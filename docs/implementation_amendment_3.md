# Implementation plan — third audit follow-up (amendment 3)

STATUS: worklist created from `docs/implementation_status.md` (the external
review of `docs/implementation_plan.md` revision 4 against the tree at
`31b8f6d`, v0.6.42, branch `dev_feature`). The review's verdict — Phase 0
and Phase A complete and matching the plan, ten of fifteen items landed,
every local gate green — is accepted as read. This file turns its eight
findings into sequenced, gated work that runs **before F2.1**, answers the
two decisions the review left open, and records the ritual changes the
review recommends adopting.

Findings here use the prefix **C\*** — disjoint from the plan's `A*`/`N*`
tags, the second audit's `B*` tags, and the plan's `U*`/`R*` blocks. Each
item carries a DONE marker that is filled in as the item lands, the same
way amendment 2's items were.

Predecessor files: `docs/implementation_plan_amendment.md` (first audit,
folded into the plan), `docs/implementation_plan_amendment_2.md` (second
audit, `B*` items — all landed).

**Revision 2 — this file re-checked against the tree.** Amendment 2's
findings were re-checked before adoption and three were corrected rather
than adopted verbatim; this file gets the same treatment. Every claim
below that names a line, a count or a mechanism was run against the tree
at `31b8f6d` (v0.6.42) before the work starts. Findings from that pass are
tagged **V\*** — disjoint from `A*`/`N*`/`B*`/`C*`/`U*`/`R*` — and are
folded in as **CORRECTED (V\*)** blocks under the items they touch. §0.1
is the ledger. Three of them change what an item *contains* (V1, V3, V4);
none changes the sequence or the versions.

---

## 0. What this amendment decides

The review's eight findings, and the disposition of each:

| Finding | One-line summary | Disposition |
|---|---|---|
| 1 | Needle pin fails on Linux CI; portability unexamined | **C1** — settle by probe, then re-specify per-platform |
| 2 | Interface-carrying film deleted under clamp-up where the same plain film is kept | **C2** — fix (not pin-as-deliberate); 0.6.43 |
| 3 | Two F0.2-introduced messages carry a 22-space run | **C3** — fix all five sites (**three** are this series', V1), pin exact messages, add a staged scanner tool (V2) |
| 4 | `evaluator.rs` doc comment says "until F1.7" post-F1.7 | **C3** — comment fix rides the sweep |
| 5 | F1.4 ladder row unstruck; §8 decisions 8–14 unstamped | **C3** — bookkeeping rides the sweep |
| 6 | `Layer.sub_layer_count` returns 1 for a gradient carrier, silently | **C3** — doc note now; sentinel decision logged for F3.1 |
| 7 | `uv.lock` neither tracked nor ignored | **C3** — track it (matches the repo's stated lockfile stance) |
| 8 | Phase A shipped with no user docs | **recorded** — F3.1 scope note; nothing ships before F3.1 without the gap |

### 0.1 Verification ledger (V\*)

Every claim in this file that names a line, a count or a mechanism, run
against the tree at `31b8f6d` before the work starts.

| ID | Claim as written | Verdict | Lands in |
|---|---|---|---|
| **V1** | Finding 3 is two F0.2 sites plus two pre-existing ones | **undercounts** — a **third** site comes from this series (`materials/mod.rs:107`, `010c0ef`/F1.1) and it is the most user-reachable of the three: the `MixRule` parse refusal, live on F1.5's Python surface | C3 §4.1 |
| **V2** | The scanner's allowlist "starts empty" | **wrong** — measured: 7 space-run sites (2 of them test asserts) and 15 non-ASCII literals on code lines, ~8 in production. The tool needs `#[cfg(test)]` scoping and comment awareness, or the allowlist *is* its product | C3 §4.1 |
| **V3** | C2's twin 4: an interface film under `ClampUpAlways` "converges to the floor" | **fails as worded** — the LM bound binds the BULK row at `clamp_min_nm` while C2's span floor binds `bulk + slice`; the film parks at `min + slice` | C2 §3 |
| **V4** | "A singleton-bulk span can never carry a `SpanRecipe`" | **holds except at t = 0** — a PINNED profiled film of zero thickness emits one uniform row and keeps its flag; harmless for C2's code, fatal if promoted to an assert | C2 §3 |
| **V5** | The digest "hashes `repr()`-rendered floats" | **imprecise** — `_canon` renders floats with `.hex()`; the conclusion is unchanged and slightly stronger | C1 §2 |
| **V6** | C1's probe compares like with like | **verified exact** — identical payload key list at `8a0ea5c`, `_strip_keys` a no-op at 0.6.32, `final_layer_count`'s F0.2 meaning change invisible on this design | C1 §2 |
| **V7** | `uv.lock` is "created by the `uv run` invocations the battery itself uses" | **wrong on fact** — `uv` appears nowhere in `ci.yml`, `README.md`, `pytest.ini` or `pyproject.toml`; no gate consumes the file | C3 §4.5 |
| **V8** | (new) Ground rule 7's own reality-check paragraph protects a message that no longer exists | new, low — `structure.rs:230` is now `ambient()`; both homogenize warnings are ASCII today | C3 §4.4 |
| **V9** | C2's arithmetic, its edge case, and the predicate swap | **verified exact** — `total < min_nm` ⇒ `min_nm − slice_total > bulk ≥ 0`; `is_singleton_bulk()` subsumes today's one-row case | C2 §3 |

The two decisions the review explicitly left open are answered here:

- **Finding 1 (settle how):** one throwaway-branch push of `8a0ea5c`
  answers both questions at once — whether a regression exists, and what
  the Linux digest at 0.6.32 is. §2 gives the procedure and both outcomes.
- **Finding 2 (fix or pin):** **fix.** The F0.3 deferral sentence that
  justified deleting non-singleton spans was written for *graded* spans,
  on the grounds that clamping a span up is F1.6's scale operation — and
  (a) F1.6 has since landed, so that justification now applies to
  scalable spans (which already scale, `structure.rs:1085`), and (b) for
  an interface-carrying plain film the clamp is not even a scale: it is
  `bulk = min_nm − t_slice`, well-defined, and B1 already established
  that the slice belongs to the authored thickness. Deleting the film is
  the opposite of U1, which is the policy's whole purpose. The narrow
  scope (floor side only) and the graded-span removal that stays are both
  argued in §3.

---

## 1. Sequenced items

| # | Item | Version | Depends on | Review finding |
|---|---|---|---|---|
| C1 | Settle the needle pin; re-specify it per-platform; ground rule 5 sentence | none (gate repair; test/docs/CI only) | — | 1 |
| C2 | Interface-span clamp-up: singleton-bulk spans clamp instead of deleting | 0.6.43 | C1 | 2 |
| C3 | The sweep: five message literals + staged scanner, stale comment, doc note, bookkeeping (incl. ground rule 7's stale paragraph), `uv.lock` | 0.6.44 | C2 | 3, 4, 5, 6, 7 |
| — | Phase B re-versioning (F2.1 → 0.6.45 … F3.1 → 0.6.49) | — | C3 | — |

Strict order. Nothing else on the ladder may claim "fingerprint unmoved"
until C1 has settled the gate — including C2, whose gates cite both
fingerprints.

---

## 2. C1 — the needle pin is not portable; settle it, then re-specify it

**Review finding 1 (blocking).** `validation/regression/synthesis/test_needle_pin.py:122`
compares a SHA-256 of hex-rendered floats against one constant recorded on
a **Windows release build** at 0.6.32. Linux CI red at the tip
([run 34749211804](https://github.com/opticsWolf/Navette/actions/runs/34749211804):
recorded `cf753a91…a63044d`, got `5c309813…f547da`), Windows green
locally and in CI. Whether the Linux digest *also* differed at 0.6.32 is
unknown, because the pin has never run on Linux — the whole series was
pushed as one batch, so CI has never evaluated any individual release.

Mechanics note (supports the re-specification). **CORRECTED (V5):** the
digest does not hash `repr()`-rendered floats — `_canon`
([test_needle_pin.py:71](validation/regression/synthesis/test_needle_pin.py:71))
renders every float with `.hex()`, and `repr()` is then applied to the
resulting structure of *strings*. The argument is the same and slightly
stronger: `float.hex()` is exact and bit-for-bit determined by the value
(shortest-roundtrip `repr` would also be value-stable, but hex leaves no
room to argue). So the digest moves exactly when the float **values**
move, and a platform difference in the digest is a platform difference in
the Rust LM trajectory's values (libm `exp`/`cos`/`powf` differ per
platform within their ulp contracts). The pin is measuring values, not
representation — which is what makes a per-platform digest map the right
repair rather than a tolerance.

### The probe (one push, two answers)

`8a0ea5c` is the commit that recorded the pin, on a 0.6.32 tree, and the
tree at that commit already contains the test (F0.1 added it). `ci.yml`
triggers on `push: branches: ["**"]`, so a branch push is enough — no
`workflow_dispatch` needed:

```
git push origin 8a0ea5c:refs/heads/probe/needle-pin-linux-0.6.32
```

The Linux job will fail, and its failure message **prints** the Linux
0.6.32 digest (`recorded …, got …` — the assert message at
`test_needle_pin.py:135`). Compare that value with the tip's Linux digest
from run 34749211804
(`5c30981320d71b9228a2dbc3b1942599a59a0bcdecd1c921c8b49cc8dbf547da`):

- **Outcome A — equal:** there is no regression. The pin is
  platform-divergent only, exactly as the review hypothesised. Proceed to
  the re-specification below, record both digests, delete the probe
  branch (`git push origin --delete probe/needle-pin-linux-0.6.32`).
- **Outcome B — different:** a genuine behaviour change somewhere in
  F0.1–F1.5 that Windows happens not to show. Bisect on Linux: throwaway
  branch pushes of intermediate commits (≈ ⌈log₂ 12⌉ = 4 runs), find the
  offending release, and fix it as a **versioned item before C2**. The
  re-specification then proceeds on top of the fix.

**VERIFIED (V6) — the two digests are comparable, which is the whole
premise.** Three things could have made the 0.6.32 Linux digest and the
tip's Linux digest measure different payloads, and none of them does:

1. `_digest`'s payload key list at `8a0ea5c` is identical to the tip's
   (`termination`, `final_mf`, `final_layer_count`,
   `final_total_thickness_nm`, `stagnation_detail`, `phases`, `stack`).
2. The tip's `_strip_keys(payload, {"clamp_report"})` is a **no-op** on a
   0.6.32 payload — that binary emits no `clamp_report` at all (F0.2
   added it) — so the stripped tip payload and the unstripped 0.6.32
   payload have the same shape.
3. `final_layer_count` changed meaning at F0.2 (rows → `spans().len()`),
   but the pin's design carries no profiled layer, so spans == rows and
   the number is the same on both sides.

`final_clamp_report` was never in the payload, so F0.2's other reporting
addition cannot reach the digest either. A difference between the two
Linux digests therefore means a trajectory difference and nothing else —
which is what makes Outcome A/B a real dichotomy rather than a hopeful
one.

### The re-specification (Outcome A assumed; Outcome B appends the fix first)

Replace the single constant with a platform-keyed map, keyed by
`platform.system()`:

```python
# Platform-keyed digests. The run's float values come from the Rust LM
# trajectory; libm exp/cos/powf differ per platform by design (each within
# its ulp contract), and Python's repr() is per-value stable, so the
# digest moves exactly when the trajectory's values move ON THAT PLATFORM.
# A digest recorded on one platform does not constrain another.
#
# Provenance per entry: (tree commit, runner image / OS, toolchain,
# build profile). Recorded at the 0.6.32 tree:
#   Windows: 8a0ea5c, windows-latest, MSVC, release (the original pin).
#   Linux:   8a0ea5c, <runner image>, GCC, release (CI probe run <id>).
RECORDED_DIGESTS = {
    "Windows": "cf753a910c7bf1d2e6c5ea5b9f5c696601ab675abaf497db28bac5340a63044d",
    "Linux":   "<from the probe run>",
}
```

with the test skipping — loudly, with the reason naming the procedure —
on a platform with no recorded digest:

```python
DIGEST = RECORDED_DIGESTS.get(platform.system())

@pytest.mark.skipif(
    DIGEST is None,
    reason=f"no recorded needle digest for {platform.system()} - "
           "record one at the current tree before gating it (amendment 3, C1)",
)
```

`dependency floors` (ubuntu) keeps checking the Linux digest: the payload
is Rust-computed, so a dependency-version shift that moved it would be a
real finding, not an environment quirk — that stays.

### CI hardening that rides this item

- Pin the `python` job's matrix to explicit runner images
  (`ubuntu-24.04`, and the current `windows-latest` equivalent) in the
  same commit. `ubuntu-latest`/`windows-latest` images rotate; a rotation
  that shifts libm would otherwise look like a fingerprint move with an
  unchanged tree. The provenance header then names exactly the image the
  digest was recorded on.
- **Ground rule 5 gains a sentence** (folded into the plan when this item
  lands): *"A recorded float digest is recorded per platform. 'Fingerprint
  unmoved' means every platform's recorded digest reproduces on that
  platform. A digest that moves with an unchanged tree is an environment
  move: re-record it with a provenance note in the commit message — that
  is a gate move, not a code move, and it never travels silently."*
- Adding a platform to the CI matrix requires recording its digest at the
  current tree first (the skipif message says so).

**Version:** none — test, docs and CI-file changes only (the R5 precedent:
test-only commits do not bump). **Commits:** one for the re-specification
(test + ci.yml), one for the plan's ground-rule sentence. Both pushed
immediately, so CI sees the repaired gate before C2 starts.

**Gates:** the pin test green on Windows locally; the Linux job green on
the push that carries the re-specification; every other gate uninvolved.

**DONE:** — done (2026-09-13). **Outcome A: no regression.** Probe run
[34770933572](https://github.com/opticsWolf/Navette/actions/runs/34770933572)
(`8a0ea5c` on branch `probe/needle-pin-linux-0.6.32`, deleted after): the
Linux 0.6.32 digest is
`5c30981320d71b9228a2dbc3b1942599a59a0bcdecd1c921c8b49cc8dbf547da` —
bit-identical to the tip's Linux digest from run 34749211804. The whole
series moved neither platform. Recorded Linux digest at 0.6.32: the value
above (provenance header in the test file). Runner images pinned:
`python` matrix → `ubuntu-24.04` / `windows-2025`; `dependency floors` →
`ubuntu-24.04`; the `rust` job stays on `ubuntu-latest` (Rust tests are
same-process comparisons — no recorded cross-machine constant). Ground
rule 5's cross-platform sentence folded into the plan
(`implementation_plan.md`, §0.1 rule 5). Local Windows pin green
(`2 passed`); the Linux job's verdict is read from the push that carries
this item (§7's push-per-item rule starts here).

---

## 3. C2 — an interface-carrying film is deleted where the same plain film is clamped

**Review finding 2 (measured).** `clamp_all_policy`'s clamp-up branch
selects with `sp.end - sp.start == 1`
([structure.rs:1107](rust/navette/src/smatrix/synthesis/structure.rs:1107))
— a **row count** — while the repo's own predicate for "one physical
layer" is `Span::is_singleton_bulk()`
([expansion.rs:54](rust/navette/src/structure/expansion.rs:54)), which N1
introduced precisely because an interface-carrying plain film is two rows
and still one layer. Measured: a 0.8 nm TiO₂ film with
`interface = true, interface_thickness = 0.2` under a 2.0 nm floor with
`clamp_up = true` is **removed whole** (total 100.8 → 100.0 nm), where
the identical film without an interface is clamped up to 2.0 nm and kept.
The deletion is reported (`spans_removed`), so it is not silent — but it
is the opposite of what the user asked for by selecting the policy whose
whole purpose (U1) is "set that layer to the minimum you can actually
deposit instead of deleting it."

### Why fix, and why the fix is narrow

- The F0.3 deferral sentence ("clamping a span up is F1.6's scale
  operation") was written for graded spans before F1.6 existed. It is now
  *live law* for scalable spans — the scalable branch at
  [structure.rs:1085](rust/navette/src/smatrix/synthesis/structure.rs:1085)
  scales a thin scalable span to the floor, fractions preserved — and
  *wrong* for the plain interface case, where the clamp is arithmetic,
  not a scale: the slice is interface-bounded and stays (B1's ownership
  rule), and the bulk takes the remainder.
- The remaining multi-row non-scalable spans **stay removed**, and that
  is deliberate: a plain multi-row span (several same-material films
  merged) has no principled way to choose which row grows, and a
  non-scalable profiled span clamped up would either mutate the profile's
  depth axis (the rate modes are defined on absolute depth — F1.2) or
  scale rows the user did not offer for scaling. Removal with a loud
  report is correct for those, and
  `f03_thin_graded_span_is_removed_whole_under_clamp_up_too` keeps
  pinning it.
- **Scoped out, deliberately:** the *ceiling* side keeps refusing a plain
  interface-carrying span above `clamp_max_nm` (the refusal at
  `structure.rs:1051` also keys on the row count and scalability). The
  refusal is loud, names both sides, and is safe; shrinking a plain film
  to the ceiling would also be well-defined, but changing the F0.2
  refusal's semantics is not what the review found and not what U1
  asks for. Recorded here as a known asymmetry; revisit only with a real
  user case.

### Design

The predicate at `structure.rs:1107` becomes the repo's own vocabulary:

```rust
if clamp_up_rows && sp.is_singleton_bulk() {
```

`is_singleton_bulk` already excludes graded and gradient spans (its doc
comment says so; the four-case twin pins it), so the recipe bookkeeping
(`new_recipes.push(old_recipes[sid].clone())`) rides unchanged whatever
the span carries.

**CORRECTED (V4) — "a singleton-bulk span can never carry a
`SpanRecipe`" is true in every case that matters and false at t = 0; do
not promote it to an assert.** Recipes attach on
`carrier.inhomogen || carrier.gradient.is_some()`
([structure.rs:574](rust/navette/src/smatrix/synthesis/structure.rs:574)),
read from the vector *after* the homogenize branch has cleared both flags
(`h.inhomogen = false` / `h.gradient = None`), so every homogenized film
is recipe-free and every profiled span is multi-row — the count rules
floor at 3 (gradient) and 2 (legacy) for t > 0. The hole is zero
thickness: a **pinned** profiled film at t = 0 never reaches the profiled
branch (`gradient_sub_layer_count`'s own comment: "zero-thickness layers
never reach the gradient branch"), emits one uniform row, and keeps its
flag — a singleton-bulk span with `Some(recipe)`. Harmless here: the
clone preserves it either way, and the next `refresh_profiles` re-emits
at the clamped-up total, which is the wanted behaviour. Stated so that
nobody writes `debug_assert!(recipe.is_none())` next to this branch.

The arithmetic targets the bulk row inside the copied slice:

```rust
let slice_rows = usize::from(sp.slice);
let slice_total: f64 = thin[..slice_rows].iter().map(|l| l.d_nm).sum();
thin[slice_rows].d_nm = min_nm - slice_total;
```

For the plain one-row case (`slice_rows == 0`) this is exactly today's
`thin[0].d_nm = min_nm` — the existing behaviour is a special case of the
new code, which is what keeps the change honest. Edge: `total < min_nm`
implies `min_nm − slice_total > bulk ≥ 0`, so the result is positive by
construction; the twin asserts it anyway.

The LM lower-bound coupling needs **no** widening: `lb = vec![lb_floor;
x0.len()]` (`evaluator.rs:361–366`) is a per-parameter floor and the
interface slice is not a parameter ("the interface slice is not in
`rows`" — the `Param` doc comment), so a plain interface film's bulk
parameter already gets the floor under `ClampUpAlways`.

**CORRECTED (V3) — but the two floors are then defined on different
quantities, and twin 4 below is where that shows.** C2 puts the span
floor on the slice-inclusive total (`bulk + slice ≥ min_nm`, B1's
ownership rule); the LM bound puts `bulk ≥ min_nm` on the bulk row alone,
because that is the parameter. The bound is therefore **conservative, not
exact**: under `ClampUpAlways` an interface-carrying film parks at
`min_nm + slice`, not at `min_nm`, and the clamp-up branch is then
unreachable for it through the LM path (it is still reachable through
cleanup and merge, which is why C2 is not a no-op under that policy).

Two ways to close it, and this item takes the first:

1. **Leave the bound conservative and say so.** A film that settles one
   slice-thickness above the floor is manufacturable by construction —
   the floor is a minimum, not a target — and the alternative touches
   `build_params`/bounds construction, which is F1.6 machinery under a
   fingerprint gate. The `ThinLayerPolicy` doc comment gains the
   sentence; twin 4 asserts the real number.
2. Make the bound `min_nm − slice` for interface spans. Rejected here:
   it buys exactness in a quantity nobody specified and widens C2 from a
   clamp fix into a bounds change.

Twin 4's wording is corrected accordingly below.

### Twins (the batch)

1. **The review's probe scenario, promoted to a twin:** a 0.8 nm plain
   TiO₂ film with `interface_thickness = 0.2` under a 2.0 nm floor,
   `ClampUpAlways` and `ClampUpFinal` variants → kept, `spans_removed`
   empty, slice stays bitwise 0.2, bulk == `min − slice`, and the total
   within the established post-scale tolerance (`|x − y| < 1e-9` —
   `2.0 − 0.2` need not round-trip through the sum bitwise).
2. **The one-row plain twin keeps passing unchanged** (the predicate
   subsumes it) — no edit, just proof it still holds.
3. **The graded-removal twin keeps passing unchanged** (non-scalable
   multi-row → removed whole, reported).
4. **The interface variant of the bounded-solve twin** (wording corrected
   per V3): a thin (10 nm) interface-carrying film under `ClampUpAlways`
   converges to the bound and is never deleted, mirroring the existing
   `clamp_up_always_moves_the_lm_bound_and_never_oscillates`. The
   assertion is `bulk == clamp_min_nm` — total `clamp_min_nm + slice` —
   **not** `total == clamp_min_nm`: the LM floor binds the parameter,
   which is the bulk row. Assert the merit history is monotone over the
   same five sweeps, so the twin still proves the absence of the limit
   cycle that F0.3's bound exists for.

### Docs and bookkeeping

- `ThinLayerPolicy`'s doc comment: "sets a surviving ONE-ROW span" →
  "sets a surviving singleton-bulk span (a plain film, with or without
  its interface slice); profiled multi-row spans are still removed, and
  scalable spans are scaled by F1.6's branch." Plus V3's sentence: under
  `ClampUpAlways` the LM bound binds the bulk row, so an
  interface-carrying film settles at `clamp_min_nm + t_slice`.
- The plan's F0.3 section gains a CORRECTIONS line pointing at C2 (the
  deferral sentence's justification retired by F1.6's landing; the
  plain-interface gap it masked fixed here).
- CHANGELOG `[0.6.43]` carries the behaviour change with the measured
  before/after table.
- The plan's ladder table gains the C2 row (marked when done) and the
  Phase B re-versioning of §4 below is stamped in the same docs commit.

**Gates:** cargo (workspace + the three feature variants), clippy
`-D warnings`, fmt, pytest (both fingerprints — Windows digest locally;
Linux per C1's outcome), the ten review harnesses, the four tools,
release build. **Version: 0.6.43** at all seven sites. **Pushed
immediately after its docs commit** — the first ladder release CI
evaluates individually.

**DONE:** — done (2026-09-13, feature commit pending). The predicate
swap + slice-aware arithmetic at `structure.rs`; the `clamp_all` doc
paragraph rewritten (its "does not exist yet" clause was stale too —
F1.6 landed the scale); `ThinLayerPolicy`'s deferral comment replaced by
the three-step landing story + V3's conservative-bound sentence. Twins:
`c2_interface_film_clamps_up_instead_of_being_deleted` (the review's
probe scenario, ClampUpFinal/Always/Remove variants — the Remove control
reports 2 rows removed) and `c2_interface_film_above_the_floor_is_
untouched`; the interface bound twin
`clamp_up_always_interface_film_binds_the_bulk_row_not_the_total`
(asserts `bulk == clamp_min_nm`, total `+ 0.5` slice untouched, merit
history monotone over five sweeps, never deleted; the Remove control
keeps the lead). The bound twin needed the lead's nk EXACTLY the
ambient's: optically invisible apart from a global phase, the
slice+bulk same-nk pair reduces to one L layer, and the plain twin's
descent argument carries over — with the demand's zero now at
bulk = −0.5, so the constrained optimum IS the bound. 546 lib tests,
workspace green.

---

## 4. C3 — the sweep: message rendering, stale comment, doc note, bookkeeping, `uv.lock`

One release (0.6.44), one theme: everything the review found that does
not change behaviour except message text. One feature commit (code) and
one docs commit (plan/CHANGELOG), as usual.

### 4.1 The five message literals (review finding 3 + V1's third series site + two pre-existing cousins)

Lost `\`-continuations left whitespace runs in wrapped literals:

| Site | Introduced | Defect |
|---|---|---|
| [structure.rs:1057](rust/navette/src/smatrix/synthesis/structure.rs:1057) | `e7a0e61` (F0.2) | 22-space run before "refusing" |
| [pipeline.rs:104](rust/navette/src/smatrix/synthesis/pipeline.rs:104) | `e7a0e61` (F0.2) | whitespace runs before "clamp_max_nm" and "rescaling" |
| **[materials/mod.rs:107](rust/navette/src/materials/mod.rs:107)** | **`010c0ef` (F1.1)** | **18-space run before "MaxwellGarnett" — V1** |
| [interpolate/mod.rs:210](rust/navette/src/interpolate/mod.rs:210) | pre-existing | whitespace run before "[lo" |
| [thick_opt.rs:323](rust/navette/src/smatrix/synthesis/thick_opt.rs:323) | pre-existing | whitespace run **and a non-ASCII `×`** (ground rule 7) |

**CORRECTED (V1) — three of these come from this series, not two, and the
new one is the most user-reachable of the three.** `MixRule::parse`'s
refusal is the message a user meets when they typo an EMA name on F1.5's
own Python surface. Measured at the tip:

```
>>> Layer(120.0, "TiO2", gradient={..., "ema": "Bruggman", ...})
ValueError: unknown mixing rule 'Bruggman' (one of Bruggeman,
                  MaxwellGarnett, Looyenga, Lichtenecker, MoriTanaka, PowerLaw)
```

It is also the one whose exact text is worth pinning **Python-side**: it
is the only one of the five reachable without constructing a
`DesignStack`, and F1.5's surface is what made it reachable.

The plan's "do not fix as drive-by" rule applied to the pre-existing
sites *when no message work was scheduled*; a dedicated message-hygiene
item is not a drive-by, and fixing all five in one commit is the honest
scope (the `×` → ASCII is the same ground rule the item exists to serve).
The three from this series matter most: two are the first messages a user
meets when licence item 2's refusal fires, and the third (V1) is one
typo away on every gradient a user authors.

- Fix all five literals (line-wrap the source, not the string).
- **Upgrade the existing twin**
  (`f02_new_refuses_pinned_span_above_the_ceiling`) from substring asserts
  to the **exact full message** — the messages are deterministic `format!`
  output over pinned inputs, so the strongest guard is free. Same for the
  pipeline twin, and a Python-side exact-message twin for V1's site.
- **New tool: `tools/check_message_whitespace.py`** — scans string
  literals under `rust/` for runs of ≥ 3 spaces and for non-ASCII
  characters; wired into the gate list.

**CORRECTED (V2) — the allowlist does not start empty, and the tool needs
real scoping before it can be blocking.** Measured over `rust/` at the
tip, with a literal-aware scan:

| Class | Sites | Of those, production |
|---|---|---|
| ≥ 3-space runs inside literals | **7** | 5 (the table above) |
| non-ASCII inside literals on code lines | **15** | ~8 |

What the raw counts mean for the tool:

- The two space-run sites the table does not list are **test assertion
  messages** — `evaluator.rs:1846` and `:1851`, the F1.7 staleness twin.
  They are not user-facing and their wrapping is nobody's problem, so the
  tool must skip `#[cfg(test)]` modules (and `rust/navette/tests/`) or
  carry them forever in an allowlist.
- A naive line scan also mis-reads `cycle.rs:228` — a quoted string
  inside a trailing `//` comment on a code line. Comment awareness is
  required, not optional.
- **Non-ASCII is a much bigger set than the single `×` this item names**,
  and a blanket ban is a bigger change than C3's scope: ~8 production
  literals carry `—`, `×`, `·`, `≥`, `λ`. Ground rule 7's stated reason
  is cp1252 consoles, and that sorts them: `—`, `×`, `·` **are**
  cp1252-encodable; `≥` ([config.rs:182](rust/navette/src/smatrix/synthesis/config.rs:182),
  today reachable only through `validated()`, which has no production
  caller) and `λ` ([inflate.rs:44](rust/navette/src/smatrix/synthesis/inflate.rs:44))
  are not. One more is worth naming because of where it lives:
  [navette-py/src/synthesis_pipeline.rs:924](rust/navette-py/src/synthesis_pipeline.rs:924)
  carries an em-dash on a live PyO3 error path.

So the tool ships in two stages, both inside C3:

1. **Space runs: blocking, allowlist empty**, tests excluded. That class
   is unambiguous — a run of three spaces inside a message is always a
   lost `\`.
2. **Non-ASCII: blocking for the cp1252-unencodable set** (which is what
   ground rule 7 actually protects), advisory-with-allowlist for the
   rest, and the allowlist is written out explicitly with the five or six
   entries it starts with rather than pretended away. Retiring those
   entries is F3.1's exposure re-audit, where the messages are being read
   anyway.

### 4.2 The stale F1.7 comment (review finding 4)

[evaluator.rs:425](rust/navette/src/smatrix/synthesis/evaluator.rs:425)
still says the RateCapped modes "stay homogenized until F1.7" — on the
`Param` enum, the first thing a reader of the parameter machinery meets.
F1.7 shipped: `span_is_scalable_rows` reads optimize flags only, and
`refresh_profiles` is what makes depth-dependent modes sound under
scaling. Rewrite to present tense, naming the reconciliation.

### 4.3 The `sub_layer_count` doc note (review finding 6)

`Layer::sub_layer_count()` ([layer.rs:144](rust/navette/src/structure/layer.rs:144))
is gated on `self.inhomogen`, so a gradient-only layer reports 1 while it
expands to ≥ 3 rows. The design constraint is real (the gradient rule
needs the wavelength grid and both endpoint spectra, which a bare `Layer`
does not have) but the getter is public API and silent about it. Doc
comment now: names the gradient exception and points at
`gradient_sub_layer_count` as the emitted count's source. Whether the
getter should *return* an honest sentinel (0? `None`? a documented 1?)
is a decision logged for F3.1's exposure re-audit — see §6.

### 4.4 Bookkeeping (review finding 5)

- Strike F1.4's master-table row
  ([implementation_plan.md:144](docs/implementation_plan.md:144)) — the
  only unstruck DONE row, a ground rule 1 miss.
- **NEW (V8) — retire ground rule 7's reality-check paragraph.** It reads
  "the tree does not currently keep this — the homogenize warning at
  `structure.rs:230` contains an em-dash … do not 'fix' the existing one
  as drive-by, and do not treat its punctuation as precedent." Both
  halves are now stale: `structure.rs:230` is `ambient()` (the file went
  853 → 3 060 lines), and the homogenize warnings that replaced it
  ([structure.rs:458](rust/navette/src/smatrix/synthesis/structure.rs:458),
  [:515](rust/navette/src/smatrix/synthesis/structure.rs:515)) are ASCII —
  F1.1 rewrote them. The exception the rule carved out no longer exists,
  which is exactly why §4.1's scanner can be blocking on the class it
  bans. Replace the paragraph with a pointer at the scanner and at the
  allowlist that records what is still outstanding.
- Stamp §8 decisions with RESOLVED blocks naming the commits where they
  resolved (decision 15 is the format precedent):

| Decision | Resolved at | How |
|---|---|---|
| 8 (A5) `material_b`'s nk on the design path | F1.1 (`010c0ef`), mirrored F1.5 (`3b08482`) | `GradientJson.nk_b`, injected by `assemble_stack`, refused when absent |
| 9 (A6) EMA selector | F1.1 (`010c0ef`) | `MixRule` as-is; serde representation became schema-visible and was pinned by the F1.5 nested-deny work |
| 10 (N4) `gradient` in state | F1.4 (`f956574`) | only-when-`Some`, as recommended; the fingerprint carries the nested key set |
| 12 (U2) how a span declares scalable | F1.6 (`86e4fed`), frozen by F1.4 (`f956574`) | `optimize = true`, no `scalable` field — the schema F1.4 froze has none |
| 13 (U4) refresh configurable? | F1.7 (`57b4602`) | no knob; three fixed points |
| 14 (U1) `clamp_min_nm` rename | F0.3 (`cb44c02`) | kept; both jobs documented (and C2 completes the second job's behaviour) |

  Decisions **1–5 and 11 stay open** — they are Phase B's to resolve at
  their items, exactly as the plan prescribes.

- Re-version the ladder's Phase B rows (see §5's table).

### 4.5 `uv.lock` (review finding 7)

Track it. `.gitignore`'s own header states the stance ("Source +
Cargo.lock + pyproject.toml stay tracked"; "(Keep Cargo.toml + Cargo.lock
— they pin reproducible builds)") and `uv.lock` is the Python
environment's analog. The third state (untracked, in every `git status`)
is noise in a workflow that reads `git status` fifteen times per item.

**CORRECTED (V7) — the stated reason is wrong, and the decision is
closer than it reads.** `uv` appears nowhere in `ci.yml`, `README.md`,
`pytest.ini` or `pyproject.toml`: the battery does not invoke it, and CI
never reads the file. It is local tooling (the dev interpreter is
uv-managed). So the `Cargo.lock` analogy is weaker than the sentence
above claims — `Cargo.lock` pins a build every gate consumes, `uv.lock`
pins an environment no gate consumes.

Both ends stay defensible and the item takes **track**, for one reason
that survives the correction: the battery's floors are declared in
`pyproject.toml` and CI's `dependency floors` job exists because "a floor
with no wheel is fiction" (R4.1) — a tracked `uv.lock` records what the
dev environment actually resolved to when a digest was recorded, which is
provenance C1 has just made load-bearing. If that argument is rejected,
the answer is `.gitignore`, not the status quo. Record whichever in the
commit message; do not leave the file in the third state.

**Gates:** the full battery (message twins now exact; the scanner
blocking on space runs with an empty allowlist for that class, staged
with its written-out entries for the cp1252-encodable non-ASCII set per
V2). **Version: 0.6.44** at all seven sites. The
message-text change is user-visible behaviour, which is why the sweep
gets a release rather than riding C2.

**DONE:** — (recorded here when the item lands, with the commit hashes).

---

## 5. Phase B re-versioning

The two follow-up releases consume 0.6.43 and 0.6.44. The ladder shifts:

| Item | Was | Becomes |
|---|---|---|
| C2 (interface clamp-up) | — | 0.6.43 |
| C3 (the sweep) | — | 0.6.44 |
| F2.1 environment segment schema + `bench_eval.py` | 0.6.43 | **0.6.45** |
| F2.2 K assemblies, `residuals_multi` | 0.6.44 | **0.6.46** |
| F2.3 needle + LM joint | 0.6.45 | **0.6.47** |
| F2.4 Python `environments=` surface | 0.6.46 | **0.6.48** |
| F3.1 docs, examples, release | 0.6.47 | **0.6.49** |

The plan's ladder table is corrected in C3's docs commit. Everything else
about the F2 items (their scope, their §8 decisions 1–5 and 11, their
gates) is unchanged.

## 6. F3.1 scope additions (recorded, implemented later)

- **The docs gap is now five releases deep** (review finding 8): no
  README section, example, or `docs/materials-*.md` line covers
  gradients, `ThinLayerPolicy`, scalable spans, or profile refresh. F3.1
  carries all of it; if anything ships from this branch before F3.1, the
  gap ships with it — which is a reason to keep the tail short, not to
  defer F3.1.
- The `sub_layer_count` sentinel decision (§4.3) — resolve or decline
  explicitly during the exposure re-audit.
- The two pre-existing message sites move **out** of F3.1's scope (C3
  fixes them); F3.1 keeps the `check_message_whitespace.py` tool in its
  gate list like the other tools.
- **NEW (V2):** F3.1 retires C3's non-ASCII allowlist — the ~6 remaining
  cp1252-encodable literals (`—`, `×`, `·` in `tables.rs`, `solver.rs`,
  `color_merit.rs`, `jacobian.rs`, `thick_opt.rs`,
  `navette-py/src/synthesis_pipeline.rs`). The exposure re-audit reads
  every one of those messages anyway; emptying the allowlist there costs
  one pass instead of a separate item, and the tool goes fully blocking
  at that release.

## 7. Ritual changes adopted from the review

1. **Push per item, from C1 on.** The series' twelve commits were pushed
   as one batch, so CI never evaluated any individual release — the pin
   defect (finding 1) lived unobserved through ten releases *because* of
   that. From C1 onward: each item's docs commit is immediately followed
   by `git push`, and the ladder does not advance until the pushed run is
   green. A red run now costs one item of reversion, not twelve.
2. **Ground rule 5's cross-platform sentence** (§2) enters the plan at
   C1, so "fingerprint unmoved" is defined for the whole CI matrix before
   any Phase B item needs the phrase.
3. **The probe pattern is legitimate gate work:** a throwaway branch push
   that makes CI *print* a value (here: a digest; previously nothing of
   the kind existed). `ci.yml`'s `push: branches: ["**"]` already
   supports it; record the branch in the item's DONE marker and delete it
   after.

---

## 8. Execution order (summary)

1. **C1 probe** → read the Linux job → outcome A or B.
2. **C1 re-specification** (platform-keyed digests, runner pin, ground
   rule 5 sentence) → push → green on both platforms.
3. **C2 (0.6.43)** interface clamp-up + twins + docs → push → green.
4. **C3 (0.6.44)** the sweep + scanner tool + bookkeeping + `uv.lock` →
   push → green.
5. **F2.1 (0.6.45)** — Phase B begins, one push per item, §8 decisions
   1–5 and 11 resolved at their items as the plan prescribes.
