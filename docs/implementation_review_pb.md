# Phase B review — the multi-environment series (F2.1–F3.1)

STATUS: review round 1, **refined** — every finding re-verified against
the tree, one citation corrected, M5 measured and closed, one new finding
(M7) added, M2/M3/M4 applied. §6 records what the refinement pass
changed. Scope: the six commits
`main..dev_multi` — `98ae05e` (F2.1, 0.7.1), `705b190` (F2.2, 0.7.2),
`948425f` (F2.3, 0.7.3), `90a09b3` (F2.4, 0.7.4), `e9484ca` (F3.1,
0.7.5), plus their three done-marking docs commits. The fork point is
`2bb57e2` (main at the review §8 closure), so this branch carries the
closed `nondispersive-bitwise` item from the PD series.

Reviewed against `docs/implementation_plan.md` §4–§8 and
`docs/plans/multi_environment_plan.md`. Method: full battery re-run,
every DONE section's claim checked against the code it describes, and
the new surfaces probed end-to-end from the installed 0.7.5 wheel.

## 0. Verdict summary

The series does what the plan says it does, and the things it promised
not to touch it demonstrably did not touch. Verified against the code,
not against the DONE sections: the flat path is delegated, not
re-assembled; `residuals_into` gained exactly one predicate; the K
branch is one `Option` read outside every loop; the K=1 merit is the
flat merit on the bit; the bench merit reproduces the F2.1 baseline
artifact bit-for-bit (`44dc154d6946ee40`, verified from the float, not
from the plan text); both fingerprints and the whole battery are green;
CI is green on every pushed commit.

Two real defects were found, both demonstrated end-to-end through the
public door on the installed wheel, and both are the same *class*: a
request that is accepted without a refusal and then answers a different
question than the one asked. **M1** — `per_film_flags` re-opens the
forced-false invariant on surroundings; **fixed at 0.7.6**
(`10d42dc`), see §3.1. **M7** — the merit spec and the
design are bound to each other by environment *position*, and position
is only re-checkable by count, so a roster written in a different order
in the two places is accepted and scores every demand against the wrong
surroundings. Neither is an arithmetic defect; both are refusal holes,
which is exactly why the twins did not catch them.

Everything else is bookkeeping the round-2 review of the PD series
already taught us to look for (a `PENDING` left in the progress table, a
gate tolerance that drifted from the plan text without a corrections
entry, forward-looking comments that outlived the item they point at).

| # | Finding | Sev | One line |
|---|---|---|---|
| M1 | `per_film_flags` can re-enable `optimize`/`needle` on fixed surroundings; under K > 1 the LM then optimizes environment 0's surroundings against environment-0-only residuals while every other environment keeps the compiled thickness — silently | **P1** → fixed (0.7.6, `10d42dc`) | §3.1 |
| M2 | F2.3's FD gate is 1e-4 relative in the test but the plan text says 1e-12, and no corrections entry recorded the reinterpretation — **recorded in this pass** | P2 → fixed | §3.2 |
| M3 | Progress-table rows 0.7.4/0.7.5 say commit `PENDING` (the same ledger/table split H2 fixed for 0.6.49) — **filled in this pass** | P2 → fixed | §3.3 |
| M4 | Four texts still pointed at F2.3 as future work; two were *false*, not merely stale — **reworded in this pass** | P3 → fixed | §3.4 |
| M5 | F2.2 gate (c)'s wall-time half (K=2 ≈ 2× single-solve) was pinned nowhere — **measured at this pass, and it passes** (1.99× / 2.14×) | P3 → closed | §3.5 |
| M6 | Pre-existing, faithfully mirrored: the needle seed's `nk` is looked up from the *host's* material while the seed's *name* is the winning sweep's partner — a cross-material site win pairs one name with another partner's `nk` | P3 (pre-existing) | §3.6 |
| M7 | `MeritSpec` keeps `n_envs` and discards the roster NAMES, so `run_environments` compares counts; a spec built with the same names in a different order is accepted and scores every demand against the wrong surroundings — silently | **P1** | §3.7 |

## 1. What was run (battery, this session, tree `e9484ca`)

| Gate | Claimed | Measured |
|---|---|---|
| pytest | 792 passed 1 skipped | 792 passed, 1 skipped, 8.9 s (re-run at the refinement pass: 792 / 1, 12.1 s) |
| cargo lib (base) | 585 | 585 passed, 0 failed (re-run after M4's reword: 585) |
| features | 590 / 594 / 599 | 590 (`opt-minpack-lm`) / 594 (`opt-argmin`) / 599 (all-features) |
| fmt | clean | clean (`cargo fmt --all --check`) |
| clippy `-D warnings` | workspace clean | clean (workspace, all targets) |
| five tools | clean | clean (`check_exposure` 229 pub fns / 110 allowlisted) |
| ten harnesses | all green | 10/10 green |
| fingerprints | unmoved | needle pin + roundtrip: 18 tests green; digests reproduce at 0.7.5 |
| PD twin harness | ALL OK | ALL OK (incl. the H1 half-space twin) |
| examples | all three run | both new examples run (`multi_environment_ar.py`, `rugate_gradient_vs_discrete.py`) + the pre-existing one |
| bench merit | bit-equal artifact | **bit-equal**, verified from the float: `62003.29065983792` → LE bytes `44dc154d6946ee40`, exactly the digest quoted in the plan at four places |
| CI | green per push | green on all six pushed `dev_multi` commits |

## 2. The implementation against the plan

### F2.1 — compile (0.7.1, `98ae05e`)

Verified in `rust/navette/src/smatrix/synthesis/environments.rs`:

- The absent-is-today rule is structural, not promised:
  `build_environments` on an empty roster **delegates** to
  `build_design` and wraps the result (`environments.rs:637–663`);
  there is no second flat assembler to drift from.
- `FixedLayerRow` is its own type with `optimize`/`needle` defaulting
  to `false` and refused when `true` (`environments.rs:58–101`,
  refusal at `:93` names the film) — corrections 2 and 3 hold.
- `EnvSegmentCfg` is two `Option`s plus the exactly-one check at
  `compile_env_rows` (`environments.rs:538–541`), so every refusal can
  name its environment and segment — correction 2 holds.
- The routing table is span-keyed and built from `Span::logical`
  (`environments.rs:774–795`) — correction 4 holds. *(Round 1 cited
  806–826 for this bullet, which is `resolve_env`'s range, quoted
  correctly two bullets down. Corrected.)*
- The three extra refusals (layers-with-environments,
  design-without-environments, environments-without-design) are all
  present (`environments.rs:633, 666–672`) — correction 5 holds.
- `resolve_env` (`:815`) answers `None → 0`, refuses unknown names with
  the known list, and is allowlisted with a rationale that names the
  reason (`tools/check_exposure.py:76` + comment block) — corrections
  6 and 7 hold, including F2.4's expiry repair (comment at
  `check_exposure.py:64–73`).
- `bench_eval.py` exists with the F2.1 baseline artifact and the
  module docstring states the eval-is-not-the-sum property —
  corrections 7 (A7) and 8 hold.

### F2.2 — K solves, joint merit (0.7.2, `705b190`)

- `merit_multi`/`residuals_multi` (`merit.rs:883`, `:932`) iterate
  env-major, key-minor; `merit`/`residuals` delegate through
  `std::slice::from_ref` and keep their signatures. At K = 1 the
  env loop has one element, so the accumulation order is the old one
  op for op — and the bitwise twin
  (`k1_joint_merit_is_the_flat_merit_bitwise`,
  `environments.rs:1007`) plus the bench artifact pin it.
- `residuals_into`'s only diff is the `&& t.env_idx == env` predicate
  (`merit.rs:1005`, color tail `:1220–1225`) — correction 1 holds, and
  the item's own "any diff means over-scoped" was corrected rather
  than quietly kept.
- The branch lives on `SmatrixContext.envs: Option<Arc<...>>`
  (`evaluator.rs:50–69`), read once per eval in `evaluate_merit`
  (`:322–336`), hoisted out of the LM closure as a bool (`:413–427`),
  and `run_environments` leaves it `None` at K = 1
  (`driver.rs:364`) — corrections 2 and the §4.6 single-branch
  contract hold in all four places that matter.
- `expand` copies design spans per row (`environments.rs:387–404`),
  surroundings untouched — correction 3 holds, and the alignment
  check (`:415–455`) is the fail-loud backstop.
- The `Ok(None)` Jacobian decline is explicit and documented
  (`evaluator.rs:617–627`) — correction 5 holds.
- The panic contract on a wrong-length `sims` slice is real
  (`merit.rs:970–990`, `check_sims`) and documented in `# Panics`.
- Gate (c)'s merit half is pinned bitwise:
  `two_identical_environments_double_the_merit` — doubling is exact in
  binary FP, so the env-major sum of two identical residual squares is
  2× on the bit. See M5 for the wall-time half.
- The hand oracle is the engine's own flat door bit for bit
  (`a_shared_film_behind_two_covers_solves_as_two_flat_coatings`) —
  correction 9 holds.

### F2.3 — the joint needle (0.7.3, `948425f`)

- `build_needle_targets` refuses a K-environment spec by name and
  delegates to `build_needle_targets_env` (`needle_pass.rs:237–269`),
  which filters demands by `env_idx` — the needle-side twin of the
  merit predicate.
- `SpectralInputs.folds` is one per environment, length
  `spec.n_envs()` (`pipeline.rs:426`, built at `:436`) — correction 2
  holds.
- Seed materials come from design rows of environment 0 only
  (`cycle.rs:140–152`) — correction 3 holds.
- The alignment assert is the bucket count
  (`cycle.rs:213–220`: every `(seed material, slot, step)` bucket
  filled exactly K times, named with the drift locus) — correction 4
  holds, and the deterministic sort before the min
  (`cycle.rs:224–229`) is the tie-break the identical-environments
  case needs.
- The elimination direction is closed: while `joint`, the floor is a
  hard LM bound (`evaluator.rs:446–460`) and the sweep clamps up
  instead of removing (`evaluator.rs:500–507`); `run_environments`
  refuses `enable_cleanup`, `enable_inflate` and
  `thin_layer_policy 'remove'` under K > 1 naming all three
  (`driver.rs:315–337`) — correction 1 and 8 hold, and the
  `a_joint_optimization_never_loses_a_span` twin pins the posture.
- `insert_seed` splits every environment's template at its own row
  (`environments.rs:291–361`), re-derives `intra` as
  position-within-segment (correction 7), and shifts later slots by
  two with the host's top portion keeping the span.
- The K=1/K=2 needle twin is structural with the reason stated
  (`a_joint_needle_over_identical_environments_lands_where_one_does`)
  and the N1 twin buries the interface film in *both* environments
  (corrections 5 and 6 hold).
- The FD gate exists (`the_joint_slope_is_the_sum_of_the_environments_slopes`)
  — see M2 for its tolerance.

### F2.4 — the Python surface (0.7.4, `90a09b3`)

- `run_needle(layers | design=...)`: `layers` optional, both/neither
  refused by name (`pipeline.py`, the `(layers is None) == (design is
  None)` gate); `environments` without `design` refused.
- `_environments_request` shapes the same `DesignRequest` a program
  document is, and `run_design_environments` (navette-py,
  `synthesis_pipeline.rs:591`) carries it to `build_environments` +
  `run_environments` — one compiler, one refusal set. The hex gate
  (`test_a_document_and_a_keyword_call_are_the_same_run`,
  `_hex` on `final_mf`) passes.
- The envelope: `PROGRAM_SCHEMA_VERSION = 2` +
  `MIN_READABLE_PROGRAM_SCHEMA_VERSION = 1` in both homes
  (`config.rs:36/:41`, `program.py:45/:54`), a range gate with
  two-sided messages (`config.rs:337–350`), the v1-refusal test
  flipped deliberately, and — the reason the bump exists — the
  `PROGRAM_SECTIONS` whitelist (`config.rs:50`, refused by name at
  `:664–671`). The N9 amendment's reading (silent drop, not refusal)
  is what the whitelist fixes.
- `program_to_dict` carries both sections as JSON text and
  `LoadedProgram` carries them on both `load_program` bodies
  (navette-py `structure.rs` + `program.py`) — correction 1 holds.
- `_load_program_native` quotes `PROGRAM_SCHEMA_VERSION` instead of
  the literal 1 (`program.py:300–311`) — correction 2 holds.
- `environment=""` refuses at construction (`target.py:294`), and an
  absent tag emits no key at all (verified: `_dump()` of an untagged
  target has no `environment` entry) — corrections 5 and the
  "keeps every pre-F2.1 target set meaningful" claim hold.
- Duplicate film names refuse where a name is an identity
  (`_environments_request.register`); a design row naming an unknown
  material refuses with the section in the message — correction 4
  holds.
- The F2.4 pytest figure correction (790 → 792, discrepancy recorded
  rather than overwritten) is exactly the right handling — gate
  numbers copied from a scrolling log are gate numbers that can be
  wrong.

### F3.1 — docs, examples, audit (0.7.5, `e9484ca`)

- `docs/graded-media.md` exists with the three facts (spans authored /
  rows solved with the 16 → 80 number; Névot-Croce at every sublayer
  boundary with σ ≪ thickness failing first; the ~6 nm thin-gradient
  floor). The README points at it.
- `docs/spectralweave-target-kinds.md` gains the Environments section
  and the coverage-matrix row.
- Both examples run and measure rather than assert (the AR example
  scores each design alone in each surrounding because the runs' own
  merits are not comparable — correction 2 holds; the rugate example
  states the triangular-not-sinusoidal fact in its own docstring —
  correction 1 holds; the scoring pass renames needle-inserted rows —
  correction 4 holds).
- Exposure re-audited: 229/110, clean.

## 3. Findings

### M1 (P1) — `per_film_flags` re-opens the forced-false invariant on surroundings

**What the plan promises.** §3.1/§4.2: surroundings are assembled with
`optimize/needle` **forced false**, and an explicit `true` *refuses* —
"silent freezing would lie about the optimum" and, more importantly for
K > 1, a free variable outside the shared design breaks the jointness
that is the whole point (multi_environment_plan.md §6 limitation 3: "a
film is either a shared design variable or fixed everywhere").

**The hole.** `FixedLayerRow::to_row` forces both flags false and
refuses an explicit `true` (`environments.rs:90–101`). But flag
application in `split_and_build_films` runs in the order *global flags
→ row → per-film override* (`design_config.rs:582–586`), and the
per-film override is keyed by **material code** — which is written down
for every surrounding row. So `per_film_flags = {"<code>":
{"optimize": true}}` lands **after** the forced-false row and
re-enables the flag. The `film_flags` global map does *not* have this
problem (the row wins there); only the per-film override does.

**Why it is P1 and not a nit.** The parameter list is built from the
stack's optimize-flagged spans with no design-slot filter
(`build_params`, `evaluator.rs:561–575`), and `expand` copies only
design slots' thicknesses across environments (`environments.rs:387–404`).
A re-flagged surrounding in environment 0's assembly therefore becomes
an LM parameter optimized against environment 0's residuals only, while
environments 1..K keep the compiled thickness — the exact
quiet-wrong-surroundings class the phase's refusals exist to prevent
(F2.2 correction 6's rationale, verbatim). The alignment check cannot
see it: span counts and materials are unchanged; only a thickness
moved. No twin covers it (none of the tests set `per_film_flags` on
the environment path).

**Demonstrated end-to-end, twice.** Round 1 drove
`run_design_environments` directly. The refinement pass reproduced it
through `run_needle(design=..., environments=[...])` — the keyword door,
no JSON by hand — on the installed 0.7.5 wheel: two environments, a
500 nm glass in environment 0's fixed segment and a 1000 nm one in
environment 1's, `enable_cleanup=False`, `enable_inflate=False`,
`thin_layer_policy='clamp_up_final'` (so the K > 1 refusals do not fire
first). Measured:

| call | environment 0's surrounding | shared film `cL` |
|---|---|---|
| no `per_film_flags` | 500.0000, `optimize=False` | 7.3137 |
| `per_film_flags={"glass": {"optimize": True}}` | **534.0298, `optimize=True`** | **2.0000** |
| `per_film_flags={"glass": {"needle": True}}` | 500.0000, **`needle=True`** | 7.3137 |
| `film_flags={"optimize": True}` (control) | 500.0000, `optimize=False` | 7.3137 |

Three things the table settles that round 1 asserted:

1. **No refusal fires**, and the surrounding moves — 500 → 534.03 nm.
2. **The design's own answer changes too.** `cL` lands at 7.31 nm
   unflagged and at the 2 nm clamp floor flagged. The extra degree of
   freedom is not an additive nuisance parameter: it competes with the
   shared films for environment 0's residuals, and the shared films
   lose. A reader who assumed "one stray optimized surrounding, one
   stray thickness" would be wrong about which numbers moved.
3. **The global `film_flags` map is genuinely safe** — round 1 claimed
   this from the code, and it is now measured. The order is
   `apply_flag_map(film_flags) → apply_row → per-film override`, so the
   row wins over the global map and loses to the per-film one. Only the
   per-film override is a hole, and the hole is the ORDER, not the
   mechanism.

Environment 1's glass stays at its compiled 1000.0 — `expand` copies
design slots only (`environments.rs:387–402`) — which is structural
rather than measured, and is what makes the divergence silent.

The same override also reaches `needle`: at K = 1 the run takes the
flat loop, where admissibility is the `needle` flag alone — so a
flagged surrounding can host a needle, splitting it in the one
assembly that exists. At K > 1 the joint sweep is protected
structurally (`slot_of_row` answers `None` for surrounding rows), so
the exposure there is the thickness LM only — but the compile then
carries a `needle: true` surrounding that its own schema refused at
the row level, which is the same lie told twice. One mechanism, two
doors: the JSON door and the keyword door both forward
`per_film_flags` (`_environments_request` forwards it verbatim,
`pipeline.py`).

**Recommendation.** Make the invariant structural rather than
order-dependent: carry the fixed-segment material codes (or a
per-row `from_fixed` bit) through `RowAssembly` and, after the
per-film override, refuse an override that sets `optimize` or
`needle` on a fixed-segment row — naming the environment, the segment
and the code, the way `to_row` does. Design-row overrides stay
legitimate (a design film's flags are addressable by name because the
name is the identity). Twins: (a) the compile refuses the
re-enable, naming all three; (b) the M1 probe request refuses at
compile instead of running. This is engine-behaviour repair and takes
its own bump (0.7.6) with the two commit-and-docs pushes.

**Fixed at 0.7.6 (`10d42dc`), as recommended, with one change of shape.**
The recommendation offered "the fixed-segment material codes (or a
per-row `from_fixed` bit)"; what shipped is neither, because a material
code can name several rows and a bit on `LayerRow` would travel further
than the invariant does. `compile_env_rows` returns the NAMES of the
rows it built from a fixed segment — films only, since a half-space row
there never reaches the per-film override — and `RowAssembly` carries
that set per environment, empty on the flat door. `split_and_build_films`
refuses before applying the override. Reported by the compile rather
than inferred downstream from the `{env}.fixed[{seg}][{i}]` naming
pattern, which would have worked today and is exactly the kind of
coupling that rots: a name is presentation, a forced flag is a
contract.

The refusal names the environment, the surrounding row, the flag, the
material code, and both ways out (put the flag on the design film you
meant, or move the film into a design segment). `film_flags` is left
unguarded on purpose — §3.1's own control measured it safe, because the
global map is applied *before* the row — and that exemption is now a
twin on both sides rather than a paragraph. Four Rust twins (the
refusal per flag, a benign override that still lands, a design film
that stays addressable, the flat door) and three Python ones (the
keyword door per flag, asserting it names `bare.fixed[0][0]`, and the
global-map control read off the returned stack). Lib tests 585 → 587.

*One narrowing worth recording.* Through the Python door the hole is
one row at a time: `_environments_request.register` refuses a duplicate
film name across the whole run, so every surrounding carries a unique
code and a `per_film_flags` entry can only ever name one row. Through
the raw JSON door two environments may share a material code for their
surroundings, and one override then re-flags both. The fix is the same
either way; the blast radius is not.

**Documented in the meantime.** `run_needle`'s docstring and the
environments section of `spectralweave-target-kinds.md` now name the
exception instead of promising "never optimized or split". A doc note is
not a fix and is not offered as one — it keeps the promise from being
false while the refusal is missing.

### M2 (P2 → fixed) The FD gate's tolerance drifted from the plan text, unrecorded

`docs/implementation_plan.md` §4.4 (F2.3 gates) says: "Finite-difference
check with differing surroundings: joint gradient == sum of
per-environment analytic gradients, **1e-12**." The implemented twin
(`the_joint_slope_is_the_sum_of_the_environments_slopes`,
`environments.rs:1591–1687`) inserts at delta = 1e-4 and asserts
`|measured − analytic| ≤ 1e-4 · max(|analytic|, 1e-12)`.

An FD measurement cannot reach 1e-12 — truncation error alone is
O(delta) — so the implementation's tolerance is the defensible one and
the plan's number was under-specified. But the reinterpretation is
recorded nowhere: the DONE section's gate list does not mention it,
and there is no corrections entry, while the sibling correction (F2.2
correction 8) did exactly this work for the K=2 run twin. A reader
comparing the item's gate list against the test reads a 1e-12 promise
backed by a 1e-4 check.

**Recommendation.** Bookkeeping, no bump: add the corrections entry
stating the reinterpretation (FD-vs-analytic-sum is truncation-bound;
1e-4 relative at delta = 1e-4 is the honest reading), and consider
whether the plan meant the analytic-sum identity
(Σ per-env analytic slopes vs the joint residual derivative, where
1e-12 *is* attainable) — if so, that second twin is worth having too.

**Applied** in this pass as F2.3 correction 9, including the second
half: the analytic-sum twin does not exist, 1e-12 *is* reachable for it
because no finite difference is involved, and it would pin the routing
rather than the physics. Recorded as a follow-up rather than folded into
the existing twin's tolerance — the two twins measure different things
and collapsing them would lose the one that can be tight.

### M3 (P2 → fixed) Progress-table rows 0.7.4/0.7.5 said `PENDING`

`docs/implementation_plan.md:3222–3223`. Every other row carries its
commit hash; these two do not, because the done-marking ritual
(feature commit → follow-up docs commit quoting the hash) ran for
F2.1–F2.3 (`af0eabe`, `6e8681b`, `495b4bb`) but not for F2.4 or F3.1 —
F2.4's DONE text rode the F3.1 commit and the F3.1 row cannot quote
its own hash. Same class as H2 in the PD round-2 review; round 2 fixed
0.6.49's instance and this one regressed.

**Recommendation.** One no-bump docs commit: fill both rows with
`90a09b3` and `e9484ca`, so the table's hash column is complete
through the tip.

**Applied** in this pass, plus a `PBR` row for this commit — the row
shape `9d5cbc3` established for a no-bump review commit (`—` in the
version column, the finding IDs in the notes).

### M4 (P3 → fixed) Texts that still pointed at F2.3 as future work

Round 1 called this "two reachable error texts". On re-reading it is one
reachable message and three comments, and **two of the four were not
merely stale but false**:

1. `environments.rs:419` (reachable) — `check_alignment` refused with
   "a structural move under K > 1 is F2.3's, not F2.2's". Stale: a
   caller hitting it today is looking at a drift the compile can no
   longer attribute to a pending item.
2. `environments.rs:378` — `expand`'s docstring listed "needle
   insertion, cleanup removal and inflate" as **refused by the driver**.
   *False since F2.3.* Needle insertion is not refused; it travels
   through `insert_seed`, which splits every environment's template at
   once. The driver refuses `enable_cleanup`, `enable_inflate` and
   `thin_layer_policy 'remove'` — three names, and needle insertion is
   not among them (`driver.rs:315–337`).
3. `environments.rs:18–19` (module header) — "Nothing here is reachable
   from Python yet." *False since F2.4.* `run_design_environments` is
   the door and `run_needle(design=..., environments=[...])` the Python
   spelling. Round 1 missed this one.
4. `environments.rs:1171` + the assertion at `:1196` — the twin's
   comment and its `e.contains("F2.3")` check, which pinned the stale
   wording in place.

`driver.rs:291` ("the compile has no inverse of `insert_seed` yet") was
checked and is **correct** present tense — there is no inverse, and the
"yet" is a real future rather than a shipped past. Left alone.

**Applied** in this pass: all four reworded, the assertion retargeted to
`insert_seed` (the word the new message actually contains), `cargo fmt`
clean, 585 lib tests still green, message hygiene OK. Message text is
observable, so it is noted here rather than filed as a comment-only
change — no version bump, nothing computed changes.

*Sibling, out of scope:* `structure.rs:133` says of `spans` that
"nothing observable of any existing run reads it yet", which F1.6/F1.7
made false when spans became first-class optimizer parameters
(`build_params`). Pre-existing, outside Phase B, recorded for whoever
touches F0.1's bookkeeping next.

### M5 (P3 → closed) Gate (c)'s wall-time half — measured at this pass

§4.6 (c) asks for two things: K=2-identical merit ≡ K=1 exactly (pinned
bitwise, `two_identical_environments_double_the_merit`) and "K=2
wall-time ≈ 2× single-solve ± concatenation noise". Round 1 found no
measurement of the second half anywhere — not the F2.2 DONE section, not
the CHANGELOG, not a test; the bench measures K=1 only, and §6
limitation 6 documents the ×K cost as inherent without measuring it.

**Measured now, and it passes.** 12 films, 61 wavelengths, 2 angles,
identical surroundings, LM pinned to one iteration, best of 9 whole
runs, two independent sessions:

| K | best run | vs K=1 |
|---|---|---|
| 1 (flat path) | 1.748 / 1.713 ms | 1.00 |
| 2 | 3.484 / 3.661 ms | **1.99 / 2.14** |
| 3 | 4.918 ms | 2.81 |
| 4 | 6.104 ms | 3.49 |
| 6 | 8.712 ms | 4.99 |

K=2 is 2× K=1 within the noise, and beyond K=2 each further environment
costs 1.31–1.35 ms against a 1.71–1.75 ms first one — about 0.74× a full
environment, because the per-run fixed cost (request JSON, compile,
spec) amortizes. Nothing superlinear; no concatenation penalty visible
at this granularity.

**Why the LM has to be pinned, which is the part worth keeping.** A run
left to converge answers a different question. The K=2 residual vector
is twice as long, so the LM takes a different number of iterations, and
K=1 additionally takes the flat path (`envs: if k == 1 { None }`) rather
than the segmented one. Measured that way the *same problem* reads
**11× K=1** (42.7 ms → 477.0 ms), and that number says nothing about the
branch — it is an iteration count wearing a wall clock. A whole-run
ratio cannot gate (c); pinning the evaluation count is what leaves
per-eval cost as the only variable. Recorded because the wrong
measurement is the easy one to take and it looks alarming.

**Disposition.** Recorded in the F2.2 corrections (entry 10) as
measured, with the pinning caveat. It is a session measurement, not a
test: a real gate wants a bench, and a bench wants `simulate_all` /
`merit_multi` on the Python surface (neither is exposed) or a Rust-side
harness. That is the follow-up; the gate itself is no longer
unanswered.

### M6 (P3, pre-existing — not a Phase B defect) Seed name vs seed nk pairing

In the flat needle cycle, the seed is built as `material = <the
winning sweep's partner name>` but `nk = contrast[<host material>].nk`
(`cycle.rs`, insertion step 4; the joint cycle mirrors it op for op at
`joint_cycle`). When the best site of one partner's sweep sits in a row
of a *different* material, the inserted seed pairs one partner's name
with the other partner's `nk`. Legacy port behavior
(`needle_synthesis.py` byte-comparability), faithfully mirrored, and
almost always invisible because a symmetric contrast map makes the
cross-material win rare. Recorded here so a future needle review does
not attribute it to Phase B, and because the joint path now has two
places to keep in sync if it is ever fixed.

### M7 (P1) The spec/design roster is bound by POSITION, and only the count is checked

**What the code promises.** `run_design_environments`' own doc comment:
"`spec` must carry the same number of environments as the request
compiles to, or the core refuses naming both counts: a demand scored
against the wrong surroundings is exactly the silent-wrong-answer this
whole phase exists to prevent"
(`navette-py/src/synthesis_pipeline.rs:580–583`). The target-kinds doc
says the same thing from the demand side: an unknown name refuses,
"because a typo has no shape error of its own".

**The hole.** Both statements are about *names*; the check is about
*counts*. `build_merit_spec(environments=[...])` resolves each demand's
tag to an index through `resolve_env` and then throws the roster away —
the compiled `MeritSpec` carries `n_envs: usize` and no names
(`merit.rs:607–633`). So the run door has nothing to compare but the
count (`driver.rs:313–321`). A roster written with the same names in a
different ORDER passes.

**Demonstrated end-to-end** (installed 0.7.5 wheel, `run_needle`, two
deliberately non-interchangeable environments — a 10 nm cover and a
2000 nm one — with asymmetric demands, R=0 on `bare` and R=0.25 on
`lam`):

| spec roster | request roster | result |
|---|---|---|
| `["bare", "lam"]` | `["bare", "lam"]` | `final_mf` 19876.13664037235 (`a442b7be0869d340`) |
| `["lam", "bare"]` | `["bare", "lam"]` | **completes, no refusal**, `final_mf` 8523.676912295265 (`36e50fa5d6a5c040`) |
| `["bare", "typo"]` | `["bare", "lam"]` | refuses — at `build_merit_spec`: `unknown environment "lam"` |
| `["bare"]` | `["bare", "lam"]` | refuses — same place, same reason |
| `TargetCollection` passed straight to `run_needle` | — | 19876.13664037235, bit-identical to row 1 |

The failure is precisely bounded. A **typo** refuses, early and well. A
**length mismatch** usually refuses for the same reason (a tag with
nowhere to resolve), and where it does not, the run door's count check
catches it — that check does have a job. A **permutation** has no shape
error at either gate: every name exists, the count matches, and the run
reports a merit for a pairing nobody described.

**Why P1 and not a nit.** Same class as M1 and the class the phase's
refusals exist to prevent; no unusual flag required; and the calling
style it bites is the one the docs recommend for reuse ("the roster must
match the one given to `build_merit_spec(environments=...)`" — an
instruction to write the list twice). The last row of the table is the
mitigation that already exists and was written down nowhere: handing
`run_needle` the `TargetCollection` makes the two rosters the same
object by construction (`pipeline.py`, `_run_needle_environments` builds
the spec from `request["environments"]`).

**Recommendation.** Carry the roster on `MeritSpec` — the names it was
compiled against, set by `compile_merit_spec`, where they are already in
scope — and have `run_environments` compare names position by position,
refusing with both lists when they differ. Additive: a spec compiled
without a roster (every pre-F2.4 spec, every K=1 spec) keeps today's
count check, so nothing existing refuses. Twins: (a) a permuted roster
refuses, naming both orders; (b) the matching roster still runs and is
bit-equal to the `TargetCollection` path. Same file set and the same
class of repair as M1, so this review put both in one bump; the
project's one-item-one-bump ritual wins instead and M7 takes **0.7.7**,
immediately after M1's 0.7.6 — shipping one refusal hole while leaving
its twin open would be an odd place to stop, but that is about the order
of the two bumps, not about merging them.

**Documented in the meantime.** The environments section of
`spectralweave-target-kinds.md` now states the positional binding and
names passing the collection as the safe spelling; `run_needle`'s
docstring says "in the same order" rather than "must match".

## 4. Probes and verification notes

1. **The bench artifact, from the float.** The plan quotes
   `44dc154d6946ee40` as the baseline merit's bit pattern four times
   (F2.2, F2.3, F2.4, F3.1 DONE sections). Verified directly:
   `62003.29065983792` packs (little-endian) to exactly that. The
   current tree's merit is bit-equal — the "K=1 branch costs nothing"
   gate's correctness half holds at 0.7.5 without trusting the plan.
2. **Baseline provenance is honest.** `eval_baseline.json` records
   commit `2bb57e2` and version `0.7.1`. `2bb57e2` is the `dev_multi`
   fork point (= `main`'s tip); the recording ran after the version
   bump was applied to the working tree but before the F2.1 feature
   commit existed. F2.1 compiles and does not evaluate, so the eval
   path at that tree IS the pre-segment path — the baseline is a valid
   pre-change measurement with an accurate commit label.
3. **Bench at 0.7.5 (this session, `--repeats 9` equivalent):**
   assemble 131.5 µs (baseline 126.2, +4.2%), simulate 37.5 (38.9),
   merit_only 3.3 (3.6), simulate_merit 40.8 (42.5), eval 299.8
   (292.0, +2.7%). Consistent with F2.4 correction 6's machine-state
   attribution; the merit bits are unchanged, which is the half that
   cannot drift for machine reasons.
4. **The M1 demonstration** (§3.1) was run through the public
   `run_design_environments` door on the installed wheel — no test
   harness, no internals. The refusal never fires; the surrounding
   moves. The refinement pass reproduced it through the *keyword* door
   (`run_needle(design=..., environments=[...], per_film_flags=...)`)
   and added the two controls in §3.1's table: the `needle` variant of
   the same override, and `film_flags` (the global map), which is safe.
4a. **The M7 demonstration** (§3.7) is the same shape: `run_needle` on
   the installed wheel, two environments whose surroundings differ by
   two orders of magnitude in thickness, asymmetric demands so a swap
   cannot be invisible. The permuted roster completes and returns a
   different merit; the `TargetCollection` path returns the correct one
   bit for bit.
4b. **Re-measured after this pass's edits** (the Rust message reword,
   the docstrings, the plan and CHANGELOG text): fmt clean; clippy
   `-D warnings` clean on the workspace and on each of `opt-minpack-lm`,
   `opt-argmin`, `--all-features`; 585 / 22 / 15 lib-parity-py and 15
   doc; feature variants 590 / 594 / 599; pytest 792 passed 1 skipped;
   five tools OK (`check_exposure` 229 / 110 unchanged — the reword
   added no pub fn); ten harnesses 10/10; both fingerprints 5 passed;
   bench merit recomputed from the float and **bit-equal**
   (`44dc154d6946ee40`); bench phases assemble 143.3 µs / simulate
   38.8 / merit_only 5.0 / simulate_merit 43.8 / eval 302.0, the same
   machine-state band §4.3 records.
5. **Fork-point note.** `2bb57e2` includes `main`'s `7649bc3`
   (review PD §8: the `nondispersive-bitwise` literals re-derived from
   a fresh v0.6.44 worktree build, all four bit-exact). The open item
   the 0.7.0 changelog carried is closed on the line this branch
   forks from; nothing here re-opens it.
6. **The single-branch contract, verified in code** at all four read
   sites: `evaluate_merit` (`evaluator.rs:322–336`), the LM residual
   closure (`:413–433`), `run_needle_cycles` (`cycle.rs:280`, read
   once outside the loop), and `run_environments`' `envs: if k == 1 {
   None }` (`driver.rs:364`). K = 1 segmented takes the flat path
   mechanically, which is why the bitwise twins have something to
   compare against.

## 5. Disposition

**Done in this pass** (one no-bump commit, the `138bdcb`/`9d5cbc3`
shape — nothing computed changes, seven version sites unmoved at
0.7.5):

- **M2** — F2.3 correction 9: the FD gate's real tolerance, why 1e-12 is
  unreachable for a finite difference, and the analytic-sum twin that
  *could* be tight, recorded as a follow-up.
- **M3** — progress rows 0.7.4 / 0.7.5 carry `90a09b3` and `e9484ca`;
  a `PBR` row added for this commit.
- **M4** — all four texts reworded, the twin's assertion retargeted.
- **M5** — measured and closed; F2.2 correction 10, with the pinning
  caveat and the 11× wrong-measurement recorded beside the 2×.
- **M1 / M7 documentation** — `run_needle`'s docstring and the
  environments section of `spectralweave-target-kinds.md` now state both
  gaps. Doc notes, not fixes.
- **M1 itself** — fixed at **0.7.6** (`10d42dc`): the assembler-side
  refuse in `split_and_build_films`, the fixed-row set reported by
  `compile_env_rows` and carried on `RowAssembly`, seven twins, and the
  two documentation paragraphs rewritten from "open gap" to "refuses".
  §3.1 records how the shipped shape differs from the recommendation.

**Open, for 0.7.7** (engine repair, feature commit + docs commit, push
per item):

- **M7** — the roster on `MeritSpec` + a name comparison at
  `run_environments`, additive so nothing existing refuses, + two twins.
  Same file set and the same class of repair as M1. The review put both
  in one bump; the project's one-item-one-bump ritual wins, so M7 takes
  its own.

**Left recorded, not fixed:**

- **M6** — pre-existing needle seed name/nk pairing. Fold it into future
  needle work rather than into this phase.
- `structure.rs:133`'s stale spans note (§3.4), same reasoning.

The battery, the fingerprints and the bit-exactness claims all verify at
0.7.5, and re-verify after this pass. Both P1s are refusal holes rather
than arithmetic defects: **nothing computed is wrong for any request the
compile accepts on the intended spelling** — which is exactly why they
slipped past every twin. Twins test what the code does with a request it
was designed for; a refusal hole is a request nobody thought to write
down.

## 6. What the refinement pass changed

Round 1's text is kept except where it was wrong; this section is the
record of the differences, so nothing was silently overwritten.

| Change | Why |
|---|---|
| **M7 added** (P1) | Not found in round 1. The roster/spec binding is positional and only the count is checked; demonstrated with a permuted roster that runs and gives a different answer. |
| **M5 closed** | Round 1 said "pinned nowhere" and offered measure-or-strike. Measured: K=2 = 1.99× / 2.14× K=1. The gate passes. |
| **M4 enlarged and applied** | Round 1 found 2 sites and called both stale. There are 4, and 2 are false rather than stale (needle insertion is not refused; the module is reachable from Python). |
| **M1 strengthened** | Reproduced through the keyword door, with the `film_flags` control and the `needle` variant measured rather than argued, and the shared film's own answer shown to move (7.31 → 2.00 nm). |
| **One citation corrected** | §2's routing-table bullet cited `environments.rs:806–826`, which is `resolve_env`. The routing table is `774–795`. |
| **M2, M3 applied** | Bookkeeping, as round 1 recommended. |
| Battery re-run | 792 pytest / 585 lib / five tools / the two examples, all green before and after this pass's edits. |

Round 1's verdict — "the series does what the plan says it does" — is
unchanged by any of it. Every arithmetic claim it verified verified
again; what the refinement pass found is one more hole of the kind it
had already identified as this phase's characteristic risk.
