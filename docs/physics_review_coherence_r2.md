# Physics review — coherent / incoherent propagation, round 2

STATUS: independent review round 2 of `docs/physics_review_coherence.md`
(C1-C7). **Applied so far: C8** (`47b2e18`, 0.7.8), **C1's refusal half**
(`f6f959d`, 0.7.9), and **C2 + C9** (`6380ca0`, 0.7.10). C3-C7 and C10 are
not.
This round re-verified every round-1 finding
against the code at 0.7.7 (`dev_phase_physics`, review round 1 committed as
`903fda7`), extended the scope to every door that reaches a solve, and ran
fresh probes on the rebuilt 0.7.7 wheel. Three new findings (C8-C10).

## 0. Verdict

**Round 1 is correct in all seven findings; two need sharpening (C2's
cross_T object, C4's coherent sign-erasure), one gains a sharper form
(C1: the honored path already exists and is exposed).** The photometric
core verdict of round 1 is not revisited — this round found nothing that
touches the cascade arithmetic itself.

New findings:

| # | Finding | Sev | § |
|---|---------|-----|---|
| C8 | The coherence flag is an unchecked `i32` with an `== 1` gate: any other non-zero value partitions the stack **without** the attenuation element — a flag typo silently erases the flagged layer's absorption | **P2** | §4.1 — **FIXED** `47b2e18` (0.7.8) |
| C9 | Mode A's `cross_T` is a third, undocumented object (product of per-block `t_p·t_s*` with no bounce series) — C2's disposition must cover it, and the record should state what it is | P3 | §4.2 — **RECORDED** `6380ca0` (0.7.10) |
| C10 | Guided-mode tooling (`smatrix/optimizer.rs`) is single-block by construction — correct for mode search, but silent; folds into the C6 doc sweep | P3 | §4.3 |

Round-1 findings, as verified this round:

| # | Round-1 claim | Verdict | § |
|---|---------------|---------|---|
| C1 | synthesis ignores `coherent: false` | **CONFIRMED**, sharpened: the honored path already ships | §3.1 — refusal **FIXED** `f6f959d` (0.7.9) |
| C2 | Mode A Stokes mixes objects | **CONFIRMED** + cross_T precision | §3.2 — **FIXED** `6380ca0` (0.7.10) |
| C3 | no thickness sanity for flagged layers | **CONFIRMED** | §3.3 |
| C4 | gain mangled differently per path | **CONFIRMED**, coherent path erases the sign entirely | §3.4 |
| C5 | half-space flags silently ignored | **CONFIRMED** (and the docstring invites it) | §3.5 |
| C6 | caveat on the wrong docstrings | **CONFIRMED** | §3.6 |
| C7 | no independent incoherent validation | **CONFIRMED** | §3.7 |

## 1. Method

Three tools, deliberately different in failure mode:

1. **Semantic code graph** (CodeGraph over the whole workspace, 6884
   nodes / 23545 edges) to enumerate every caller/callee relation touching
   coherence — then hand-verification of each edge, because the graph does
   not see through PyO3 `wrap_pyfunction!` doors (§2).
2. **Line-level reading** of every site the round-1 review cited, plus the
   sites it did not visit (needle pass production kernels, FFI doors,
   Structure bridge, expansion emission, serialization).
3. **Fresh numerical probes** on the installed 0.7.7 wheel, on stacks built
   independently of the round-1 probe stacks, so agreement is not a shared
   fixture.

## 2. Door inventory — every path that reaches a solve

Enumerated from the graph and closed by hand (the graph's one blind spot:
PyO3 `wrap_pyfunction!` exports are invisible to it — `solve_arrays`
reported "no callers" while being exported at `navette-py/src/
structure.rs:1662`; every FFI door below was confirmed by grep).

| Door | Where | Coherence flags | Mode | Notes |
|------|-------|-----------------|------|-------|
| `ScatterMatrix` (Python class) | `smatrix.py:463` -> `solver_new` (`smatrix.rs:159`) | honored | honored; validated Python-side via `CoherenceMode(int(...))` (`smatrix.py:580`) | the analysis door |
| `solve_structure` (bridge) | `structure/__init__.py:103` -> `solve_arrays_fn` | honored (Structure layers) | opt-in, default FRONT_BLOCK | flags flow through `get_solver_inputs` -> `SolverArrays.incoherent` |
| `solve_arrays_fn` (raw FFI) | `navette-py/src/structure.rs:1662` | honored (`&[bool]` -> i32) | **unvalidated i32** | thin over `solver::solve_arrays` (`solver.rs:1111`) |
| `SmatrixContext.simulate` (synthesis funnel) | `evaluator.rs:78/92/109` | **IGNORED** (C1) | n/a | one coherent block; 18 callers in driver/cycle/environments/evaluator |
| `needle_gradient` (Python needle door) | `needle.py:225` | honored | multiblock (Pmb), host refusal `needle_operator.rs:861` | flag gate `== 1` (C8) |
| synthesis needle pass | `needle_pass.rs:882/:902` | **IGNORED** (C1) | single block by construction | `p_coherent_*_from_fields` over `[0, nl-1)` |
| guided-mode tooling | `smatrix/optimizer.rs` (landscape + Nelder-Mead) | n/a (single block `0..len-1`) | n/a | coherent-stack concept; doc line only (C10) |

Consequence for the disposition: C1 is not "the feature is missing" — it
is "one of seven doors bypasses a flag the other six honor". A refusal in
the synthesis doors has a named, working alternative to point at.

## 3. C1-C7, verified line by line

### 3.1 C1 — the synthesis pipeline ignores `coherent: false` (P1) — CONFIRMED

> **Refusal half applied** at `f6f959d` (0.7.9): `DesignStack::from_parts`
> refuses a film with `coherent = false` and names the honored doors. The
> implementation half — wiring the synthesis funnel to
> `p_function_multiblock` — is still open. Everything below describes the
> state before that commit.

The full chain, each link verified in the tree at `ee1488b`:

* `src/navette/synthesis/pipeline.py:86` — `_FILM_DEFAULTS` ships
  `coherent=True` and accepts the flag on every film;
  `pipeline.py:228` writes it into the film dict.
* `rust/navette/src/smatrix/synthesis/design_config.rs` (`apply_flag_map`)
  parses `"coherent"` for the design door, the per-film override and the
  multi-environment `film_flags`/`per_film_flags` alike — one parser for
  all three spellings, so all three accept it.
* `rust/navette/src/smatrix/synthesis/structure.rs:1371`
  (`DesignStack::solver_arrays`) stores it:
  `incoherent_flags[slot] = i32::from(!layer.coherent)` — the flag is
  faithfully materialized into `SolverArrays`.
* `rust/navette/src/smatrix/synthesis/evaluator.rs:109` (`simulate_inner`)
  then reads `n_stack_cache`, `thicknesses`, `rough_vals`, `rough_types`
  from that same struct and calls
  `solve_coherent_block_fields_dual(0, nl-1, ...)` — `incoherent_flags`
  is built and **never read**. Confirmed via the graph:
  `solve_point`'s only production caller is `ScatterMatrix::solve`
  (`solver.rs:500`); no `synthesis/` symbol reaches it.
* `simulate` (the synthesis solve) has 18 callers across driver, cycle,
  environments and evaluator — the bypass covers the design door, the
  multi-environment door, the LM and the FD jacobian (all consume
  `SimCurves` from this one funnel).
* The needle half: `synthesis/needle_pass.rs:882` ("fully coherent block
  path") and `:902` ("Fully-coherent stacks: (0, nl-1)") — the scan runs
  the `p_coherent_*_from_fields` kernels over
  `build_stack_fields(0, nl-1)`, a single block. `p_function_multiblock`
  (`needle_operator.rs:1307`) exists with FD twins at `:2126/:2235/:2297/
  :2360` and a Python mirror (`test_physics_mirror.py`, "multiblock FD +
  single-block reduction") — validated, unwired.

Fresh probe (own stack, own numbers): a design
glass 50 um / L 110 nm, glass flagged `coherent: False` through
`stack_from_layers` + `SmatrixContext.simulate` + a flat-R merit:

```text
  merit flag=False  4469fe548c03c140  8711.09634380474
  merit flag=True   4469fe548c03c140  8711.09634380474   bit-identical
```

while the engine door on the same stack honors the flag (Rs at 500 nm
differs between the flagged and unflagged solve; round-1's merit-swing
table already showed the optimizer consequence).

**Sharpening (what round 1 undersold).** The machinery is not merely
"exists and is validated" — it is **already wired and Python-exposed on
the other door**:

* `src/navette/structure/__init__.py:103` (`solve_structure`) expands a
  `Structure` to `SolverArrays` — carrying each layer's coherent flag —
  and passes them to `solve_arrays_fn`, which honors them; it accepts
  `coherence_mode` in `solver_opts` (`:146`, default `FRONT_BLOCK`).
* `src/navette/smatrix/needle.py` (`needle_gradient`) passes
  `stack.incoherent_flags` (`:225`) and exposes the multiblock gradients
  `Pmb`/`P_MB*` plus the flagged-host refusal.

So the disposition is stronger than "implement or refuse": the honest
minimum refusal can name the working sibling (`solve_structure`) in its
message, and the implementation option is wiring, not new physics.

### 3.2 C2 — Mode A's Stokes vector is not a Stokes vector (P1) — CONFIRMED

Line evidence, all verified:

* `core_engine.rs:371` — `track_cross_channel = need_cross &&
  coherence_mode == MODE_B` (mode-gated, as round 1 said).
* `core_engine.rs:528-532` — Mode A returns
  `cross_r = rp0 * rs0.conj()` (first block) while `rs`/`rp` are the
  full intensity cascade.
* `solver.rs:994/:1003` — `s0r = rp + rs` (totals) vs
  `s2r = -2.0 * s.cross_r.re` (front block in Mode A): two different
  stacks in one Stokes vector, and `DOP_R` (`:1030`) divides them.
* `solver.rs:1022-1033` — the R6.2 clamp comment says "for a single
  coherent block"; the clamp is `.min(1.0)`, blind to a deficit.
* `refs/loom_matrix.py:621` — the quoted comment is verbatim in the tree.

Fresh probe (air / 95 nm H / 50 um flagged slab / 110 nm L / air, own
indices — different from round 1's, same mechanism):

```text
  angle    DOP_R (A, default)   DOP_R (B)     Delta_R (A)    Delta_R (B)
   0       0.763904485          1.000000000   +180.000       +180.000
   45      0.925557318          0.992686001   -163.057       -163.831
   70      0.858980933          0.884129668   -56.812        -35.088
```

Control confirmed as round 1 measured it: with no flag anywhere, A and B
agree bit-for-bit; `A/B intensities bit-equal: True` for `Rs`/`Tp` on the
flagged stack (the documented A/B photometric identity holds).

Two scope notes round 1 did not have:

* **The hazard is live on the Structure door too.** `solve_structure`
  defaults `coherence_mode=0` (`structure/__init__.py:146`), so the
  default DOP through the documented engine path is the Mode A object
  whenever a layer is flagged. Round 1 framed C2 as engine-door-only with
  the default mode; the door inventory (§2) shows the default is what
  `solve_structure` hands out unless the caller opts into Mode B.
* **Synthesis cannot request cross observables at all** (`merit.rs` has
  no DOP/DELTA/STOKES observables), so C2 is engine-door-only today —
  which is exactly why it can be fixed by refusal without touching
  parity.

**As applied** — `6380ca0`, 0.7.10. Refusal at both engine doors, on
`FRONT_BLOCK ∧ (mask & NEEDS_CROSS) ∧ any INTERIOR flag`. The second scope
note above is why the Rust half sits in `solver::solve_arrays`: that is
what `solve_structure` and the raw FFI both reach, so the Structure door is
covered by the same check. Interior-only keeps C2 from contradicting C5 —
a stack flagged on its half-spaces alone still answers, because the sweep
never consults those flags. The engine is untouched, and
`test_the_engine_still_computes_what_the_doors_refuse` fails first if the
check ever leaks inward.

### 3.3 C3 — no thickness sanity on flagged layers (P2) — CONFIRMED

`from_raw` (`solver.rs:392`) validates wavelengths, angles and index
shapes; thicknesses are checked `>= 0` and finite, nothing checks optical
thickness against the flag. Probes, own stack:

* d = 0 flagged interior layer: Rs = 0.109790 against the coherent
  0.193432 — the join breaks the p-s phase relation regardless of
  thickness, so a zero-thickness layer decoheres by 8.4e-2.
* Frustrated TIR, air gap flagged between two n = 1.52 half-spaces at
  45 deg (beyond the 41.2 deg critical angle): d = 200 nm gives
  `R = 1.000000, T = 0.000000` where the coherent stack gives
  `R = 0.763, T = 0.237`; at d = 2000 nm both agree at `R = 1`. The
  flagged answer is the thick limit, unconditionally — the model cannot
  see that the gap it was handed is thin.

Round 1's `tmm` comparison and disposition (warning at construction +
`L_C = lambda^2/dlambda` docstring; refuse flagged layers as optimization
parameters if C1 is ever implemented) stand as written.

### 3.4 C4 — gain media are silently mangled, differently per path (P2) — CONFIRMED

Line evidence: the coherent path flips decay by conjugation
(`coherent_block.rs:146` and `:345`, and the needle kernels mirror it at
`needle_operator.rs:248`); the incoherent path clamps
(`core_engine.rs:503` and `:665`, `spacer_tau` at `needle_operator.rs:1290`
same rule). Neither refuses. The Structure door refuses `Im(n) < 0`
(`structure/structure.rs:147-151`, "Nominal expansion produced k < 0");
`ScatterMatrix` and `solve_arrays` do not — one door guarded, two not,
exactly round 1's claim.

Fresh probe (2 um slab, |k| = 0.01, symmetric ambient):

```text
  coherent    k = -0.010   Rs = 0.006559395   Ts = 0.585491021   A = +0.407949585
  coherent    k = +0.010   Rs = 0.006559395   Ts = 0.585491021   A = +0.407949585
  incoherent  k = -0.010   Rs = 0.076953120   Ts = 0.923089546   A = -0.000042666
  incoherent  k = +0.010   Rs = 0.053518311   Ts = 0.557830426   A = +0.388651262
```

The sharpening: on this symmetric stack the coherent path gives
**bit-identical output for gain and loss** — the conjugation erases the
sign entirely, not merely the energy bookkeeping; the reported `A = +
0.408` for a gain layer is the loss-layer answer wearing its sign. In the
incoherent path the gain layer becomes the lossless slab (Ts = 0.923089546
is the exact lossless value) with a small negative absorption — the books
open by 4.3e-5. The disposition stands: refuse `Im(n) < 0` at the engine
doors with a shared explanation constant, the `AMBIENT_DROP_EXPLANATION`
pattern (`optics_core.rs:169`, pinned by
`test_both_doors_explain_it_the_same_way`), extended to cover the
needle-kernel flip sites by the same door-level refusal.

### 3.5 C5 — half-space flags are silently ignored (P3) — CONFIRMED

`core_engine.rs:394/:595` — the block sweep starts at
`current_idx + 1` and scans `while ni < idx_n && flags[ni] == 0`, and the
attenuation element is applied only under
`next_incoh < idx_n && flags[next_incoh] == 1` (`:493/:655`). Flags at
rows 0 and `idx_n` are never consulted. Fresh probe: rows 0, last, and
both flagged produce **bit-identical** Rs to the unflagged stack; the
interior flag moves it (0.394047181 -> 0.319475551).

Two things round 1 did not record:

* The Python constructor docstring (`smatrix.py:477-479`) documents the
  flag as "Non-zero where a layer breaks phase coherence (**thick
  substrate**)" — the docstring itself invites flagging the substrate
  row, which is the exact no-op this finding is about. The docstring is
  where the fix belongs, before or instead of a runtime warning.
* The model sentence exists next door: `solve_arrays` warns "first/last
  thickness is not 0 (ambient/substrate convention); the engine treats
  row 0/last as half-spaces" (`solver.rs:1130-1135`). The coherence
  sentence should be appended there and in the constructor docstring.

### 3.6 C6 — the front-block caveat is on the wrong docstrings (P3) — CONFIRMED

Verified in the tree: `ellipsometry` (`smatrix.py:635-636`), `stokes`
(`:647-648`) and `complex_amplitudes` (`:643-644`) carry no caveat;
`dispersion` (`:657-659`) and `differential_phase` (`:681-683`) carry it.
The Rust `OpticalState` doc (`core_engine.rs:294-296`) states the scope
honestly ("Complex amplitudes are the first coherent block (Modes A/B) or
the whole stack (Mode C)"); the Python convenience doors do not forward
it. Round 1's demonstration (|rs_c|^2 = 0.1517 against Rs = 0.1826 in one
output dict) is the user-visible consequence; the fix is a four-docstring
sweep.

### 3.7 C7 — no independent incoherent validation (P3) — CONFIRMED

Inventory check, whole `validation/` tree:

* The parity refs are `loom_matrix.py` — a port of the same block sweep
  (its own comment cites Katsidis & Siapkas; the `:621` sentence quoted
  by round 1 is verbatim in the tree).
* `core_engine.rs` has 3 tests, including
  `intensity_path_matches_full_path_bitwise` (`:877`) — full-walk vs
  intensity-walk self-consistency; `coherent_block.rs` has 0 tests;
  `needle_operator.rs`'s 23 tests are self-FD;
  `test_physics_mirror.py` translates Rust tests through the Python API
  "same stacks, same FD conventions, same hand values" — pinned
  self-consistency at both layers.
* No `inc_tmm`/`tmm` reference, no closed-form slab, no
  thickness-independence assert, no phase-average identity anywhere in
  `validation/` (the only "closed form" hits are colorimetry and Tauc).

Round 1's three proposed probes are indeed the missing twins. This round
adds one more candidate: the Stokes-vector phase average (the variant
that would have caught C2), which round 1 ran only for R/T.

## 4. New findings

### 4.1 C8 — the coherence flag is an unchecked `i32` with an `== 1` gate (P2)

The block sweep partitions at ANY non-zero flag
(`while ni < idx_n && flags[ni] == 0`, `core_engine.rs:394/:595`,
`needle_operator.rs:1265`), but the attenuation element is gated on
exactly 1 (`flags[next_incoh] == 1`, `core_engine.rs:493/:655`;
`partition_blocks` spacer at `:1269`; host refusal at `:861`). Nothing
validates the flag's range anywhere:

* `from_raw` (`solver.rs:392`) — no flag check;
* `solver_new` FFI (`smatrix.rs:159`, `incoherent_flags:
  PyReadonlyArray1<i32>`) — none;
* Python `_as_layer_array` (`smatrix.py:600`) — length check only;
* the bool-typed door (`solve_arrays` takes `&[bool]`) cannot produce a
  bad value; the i32-typed doors (`ScatterMatrix` via
  `_as_layer_array`, the needle door's flag array) can.

Probe — absorbing 2 um slab (k = 0.01), flag value swept:

```text
  flag =  0    Rs = 0.006559395   Ts = 0.585491021   (coherent, absorbed)
  flag =  1    Rs = 0.053518311   Ts = 0.557830426   (incoherent, tau applied)
  flag = +2    Rs = 0.076953120   Ts = 0.923089546   (split, tau SKIPPED)
  flag = -1    Rs = 0.076953120   Ts = 0.923089546   (same)
  flag = +7    Rs = 0.076953120   Ts = 0.923089546   (same)
```

With flag = 2 the stack is partitioned at the flagged layer but the
attenuation element is not applied: the slab's absorption (A ~ 0.36 in
the flag = 1 case) **vanishes from the energy books** — Ts comes out at
the exact lossless value. One typo away from a silently wrong spectrum,
on a lossless stack merely cosmetically wrong (tau = 1) and on an
absorbing one physically absurd.

~~The same asymmetry applies to `coherence_mode` at the raw FFI door
(`solve_arrays_fn` takes any `i32`; values outside {0, 1, 2} alias to
Mode A behavior via the `== MODE_B` / `== MODE_C` tests).~~ **This half
was wrong.** `Solver::validate` (`solver.rs:214`) range-checks the mode
and refuses with a named message; `from_wav_major_flat`, `from_raw` and
`new` all reach it, so `solve_arrays` and every PyO3 door inherit it, and
`validation_refuses` (`solver.rs:3392`) has pinned `mode = 5` since before
this review. Measured on the 0.7.7 wheel: `coherence_mode = 7` and `-1`
both raise `ValueError: coherence_mode must be 0 (front_block), 1
(coherency_matrix), or 2 (fully_coherent).` — and that message is
Rust-side, not the Python enum. The flag half stands exactly as written.

**Disposition (as stated).** Range-check both at `from_raw` and at the
two raw FFI doors, refusing with a shared message (the C4 shape). Small
commit, its own bump; it should land before C1's refusal because it is
strictly smaller and removes the worst silent wrongness.

**As applied** — `47b2e18`, 0.7.8. Canonicalization, not refusal, and
only for the flag:

* **Why not a refusal.** The only published statement of the contract is
  the `ScatterMatrix` docstring, "Non-zero where a layer breaks phase
  coherence", and the partition gate already honors it. Refusing a flag
  of 2 would narrow the published contract to match the buggy gate and
  would punish a caller who followed the documentation. Canonicalizing
  makes the documented contract true instead, and it is immune to a
  future gate picking either comparison, because downstream only ever
  sees 0 or 1.
* **Where.** `Solver::assemble` (`solver.rs:355`), the shared tail every
  constructor funnels through — so `new`, `from_wav_major_flat`,
  `from_raw`, and `solve_arrays` / `ScatterMatrix` / `core_engine` with
  them — plus the free `needle_gradient` (`solver.rs:1488`), whose flags
  arrive from the caller rather than from `self` and therefore never pass
  that door.
* **Pinned by** `any_nonzero_coherence_flag_is_the_same_flag`: 2, -1, 7,
  `i32::MIN` and `i32::MAX` bitwise equal to 1, the flagged answer
  asserted different from the coherent one (or the equality is vacuous),
  and tau asserted actually applied.
* **Measured after.** Every non-zero flag now gives
  `Rs = 0.054811350, Ts = 0.583945129, A = +0.361243521` — the flag-1
  answer, absorption intact — against the pre-fix
  `Ts = 0.923089546, A = -0.000042666`. Flags 0 and 1 are bit-for-bit
  unchanged, and all three fingerprints hold (bench merit
  `44dc154d6946ee40`).

### 4.2 C9 — Mode A's `cross_T` is a third, undocumented object (P3)

Round 1 (correctly) described Mode A's `cross_R` as the front block
alone. The transmitted channel is different code
(`core_engine.rs:441/:520`): in Mode A, `cross_t` accumulates
`t_p^block * conj(t_s^block)` **per block join** plus `tau` per flagged
layer — a product over blocks with no `1 - C_Ab*C_Bf` series, i.e. not
the front block and not the Mode B cascade either. The legacy-port
header (`test_core_engine_rigorous_ellipsometry.py`) names it: "a Mueller
cross-term product accumulated across blocks" — that phrase exists only
in the parity file, not in any user-facing doc.

Fresh probe, 550 nm, one flagged interior layer: cross_T = 0.887910062
(Mode A) vs 0.890209950 (Mode B) — different objects, both real.

**Disposition.** Fold into C2's fix: whatever Mode A does about
cross-channel observables (refuse/warn) must cover `CROSS_T` as well as
`CROSS_R`/DOP/DELTA — round 1's NEEDS_CROSS refusal already does, since
`REQ_CROSS_T` sits in `NEEDS_CROSS`. The record sentence belongs next to
the R6.2 clamp comment, whose "for a single coherent block" scope
condition round 1 already flagged.

**As applied** — `6380ca0`, 0.7.10. Recorded, not changed. `CROSS_T` is
refused with the other eleven bits, and
`test_mode_a_and_mode_b_cross_terms_are_different_objects` pins the two
values as distinct objects (alongside the A/B photometric identity, so the
difference cannot be read as a photometric disagreement). The record
sentence went where this section said it should, into the R6.2 clamp
comment, which now states its "for a single coherent block" scope and notes
that Mode A on a flagged stack produces a DEFICIT `.min(1.0)` cannot see.

### 4.3 C10 — guided-mode tooling is single-block by construction (P3)

`smatrix/optimizer.rs` (landscape scan, Nelder-Mead minimizer) solves
`solve_coherent_block_fields_inner(0, n-1)` — one block, flags not
consulted. Guided modes are a coherent-stack concept, so this is correct
by intent; it deserves the same one-line docstring as `field_profile`
(round 1's §5). Fold into the C6 docstring sweep; no code change.

## 5. Disposition — sequenced

Each item its own commit and version bump (docs-only entries excepted),
ordered by silent-wrongness per line of change:

1. **C8** — **DONE**, `47b2e18` (0.7.8). Shipped as canonicalization at
   `Solver::assemble` + the free `needle_gradient` rather than a refusal,
   and the flag only: the mode was already validated (§4.1). Closes the
   absorption erasure.
2. **C1 (refuse)** — **DONE**, `f6f959d` (0.7.9). The refusal sits in
   `DesignStack::from_parts`, the one internal constructor every stack
   passes through, rather than in each of the four sites where a false
   flag can enter (`design_config.rs:299/:318`, `driver.rs:115`,
   `environments.rs:102`) — one message, and no future route can slip
   past it. It names the film by index and material and names
   `solve_structure` / `ScatterMatrix` as the doors that do honor the
   flag. Ambient and substrate are not checked: they are the half-space
   rows where the engine ignores the flag too (C5), so a flag there is
   the same no-op on every door. The multiblock implementation (round 1
   option 2) stays a phase, unblocked by this.
3. **C2 (refuse/warn)** — **DONE**, `6380ca0` (0.7.10). Refuse, not warn.
   The deciding fact is that the legacy parity port enters through the raw
   `core_engine` pyfunction (`Solver::solve`), BELOW both doors, so a
   door-level refusal at `ScatterMatrix.compute` and `solver::solve_arrays`
   costs it nothing and the engine is untouched — a refusal inside
   `solve_point`/`resolve_plan` would have killed it. Predicate:
   `FRONT_BLOCK ∧ (mask & NEEDS_CROSS) ∧ any INTERIOR flag`, all twelve
   bits as one rule, interior-only so it cannot contradict C5.
   `NEEDS_CROSS` is exported and bound rather than re-declared Python-side.
   The R6.2 clamp comment, the `loom_matrix.py:621` comment and C9's
   cross_T record landed in the same commit, as did the
   `ellipsometry`/`stokes`/`complex_amplitudes` caveats (the start of C6).
   The battery passed without touching an existing test — the signal that
   the refusal is narrow.
4. **C4** — shared explanation constant + door refusal for `Im(n) < 0`.
5. **C3 + C5 + C6 + C10** — the doc/warning sweep: thickness warning at
   construction, half-space-flag sentence beside `solver.rs:1130`, four
   Python docstrings, one optimizer docstring line.
6. **C7** — land the three round-1 probes as
   `validation/review/incoherent_check.py` first (closed form,
   thickness independence, phase average); promote to regression twins
   alongside the C2 commit.

## 6. Recorded so a later pass does not re-derive them

* **The roundtrip pins the flag.** `coherent` is serialized
  (`structure/layer.rs:433`, deserialized `:482`), and the roundtrip
  regression carries a `coherent=False` layer
  (`test_roundtrip.py:31`) and compares the flag column (`:155`). Flag
  survival through save/load is not at risk.
* **Gradient emission preserves the flag.** Sublayer emission writes
  `layer.coherent` per emitted entry (`expansion.rs:519/:564/:582`,
  `:609/:615`) and expansion carries its own `incoherent` vec
  (`:70/:204`). When C1 is ever implemented, a flagged gradient parent
  would partition per emitted sublayer — that combination (hundreds of
  blocks from one profile) should be refused or explicitly pinned in the
  same commit, not discovered later.
* **The A/B photometric identity is pinned twice.** Documented in
  `docs/remediation_plan.md`, re-measured this round (bit-equal Rs/Tp on
  a flagged stack). Any C2 fix must keep it.
* **The duplicated block walk is contract-guarded.** `solve_point` and
  `solve_point_intensity` are deliberate copies with a written
  maintenance contract (`core_engine.rs:333-346`) and the bitwise test at
  `:877`. Any C1 implementation edit to the walk must be made twice; the
  comment block is the instruction.
* **The multiblock refusal is reachable and correct** — on the honored
  path. `locate_hosts_multiblock` refuses flagged hosts
  (`needle_operator.rs:861-865`); it is the C1-unreachable refusal round
  1 already noted. Its `== 1` gate joins C8's finding.
* **CodeGraph blind spot.** The semantic graph does not model PyO3
  exports as callers (it found `solve_arrays` "uncalled" while it is the
  body of an exported `#[pyfunction]`). For future reviews: enumerate
  FFI doors by grep over `wrap_pyfunction!`/`wrap_pyclass!` before
  trusting caller completeness. Everything in §2 was closed by hand.

## 7. Sources

As round 1 (§4 there): Byrnes arXiv:1603.02720 + `tmm`; Katsidis &
Siapkas Appl. Opt. 41, 3978 (2002); Troparevsky et al. Opt. Express 18,
24715 (2010). All probes this round ran on the installed 0.7.7 wheel
rebuilt from the working tree (`dev_phase_physics` @ `ee1488b` code
state, round-1 review docs commit on top); every number in §3/§4 is
reproducible from the stated stacks alone.
