# Navette Implementation Plan — differential phase (the PD series)

Successor-scoped plan, subordinate to `docs/implementation_plan.md`. It
covers exactly two defects in the differential-phase (`PDts`/`PDtp`)
surface, found by the status check of 2026-09-14:

- **(a)** differential phase is unreachable from the simulation surface —
  a simulation quantity that only the synthesis machinery can produce;
- **(b)** the reference index is a **scalar frozen at the centre
  wavelength**, so a dispersive incidence medium silently yields a wrong
  Δφ — and the same scalar is read at four more sites that must move with
  it.

Everything else about `PDts`/`PDtp` — the formula, the sign convention,
the native target machinery, the needle gain shift — was checked and is
correct. §1.4 says so explicitly so it does not get re-litigated.

---

## 0. How to use this plan

### 0.1 Item namespace

Items are **PD1–PD4**. Not `P1`: `P0`/`P1`/`P2`/`P3` are the *priority*
labels of `implementation_plan.md` §0.2 and a collision there would be
unreadable in the master table. `PD*` is verified free across `docs/`,
`rust/`, `src/` and `validation/` (only binary artefacts match `PD[0-9]`).

Namespaces already spent, for the next author: `A*`, `N*`, `U*`, `R*`,
`F*`, `§D` (main plan), `B*` (amendment 2), `C*` (amendment 3), `V*`
(amendment-3 verification ledger), `E*` (review C).

### 0.2 Ground rules

`implementation_plan.md` §0.1 applies **unchanged and in full** — one
item = one `0.0.1` bump across seven version sites in five files, N6 (no
blanket README `sed`), the release-build rule, the verification battery,
the toolchain-skew rule, bit-exactness, announced behaviour changes, and
ASCII Python-surface messages. Two clarifications this series needs:

- **Ground rule 3 now has five tools**, not four:
  `check_exposure.py`, `check_cie_sync.py`, `check_pyi_sync.py`,
  `check_toolchain.py`, `check_message_whitespace.py`. PD2 and PD4 both
  touch `.pyi` stubs, so `check_pyi_sync.py` is load-bearing here, not
  ceremonial.
- **Ground rule 7 / amendment 3 §7.1 — push per item.** Each item's docs
  commit is followed immediately by `git push`, and the ladder does not
  advance until that run is green.

### 0.3 Disposition

| Item | What it settles | Version | Priority | Risk | Effort |
|---|---|---|---|---|---|
| **PD1** | The reference index becomes wavelength-dependent in the core | 0.6.45 | **P0** | M (touches the merit inner loop and the needle gain shift) | M |
| **PD2** | The FFI and Python doors accept a per-λ reference; the scalar door gets a guard | 0.6.46 | P1 | S | S |
| **PD3** | GD/GDD over Δφ: the convention the fix disturbs, decided and pinned | 0.6.47 | P1 | S (a decision + twins, little code) | S |
| **PD4** | `PDts`/`PDtp` become first-class `compute()` observables | 0.6.48 | P1 | S | M |

**PD1 blocks PD4** and the ordering is not cosmetic: exposing differential
phase on the simulation surface while the reference index is frozen would
propagate (b) to a second, wider surface and give it a public API to be
compatible with. Fix the quantity, then publish it.

### 0.4 Version ladder, and the Phase B interaction

Phase B currently holds 0.6.45–0.6.49 (amendment 3 §5, itself already one
re-versioning). This plan takes 0.6.45–0.6.48 and pushes Phase B down by
four:

| Item | Was | Becomes |
|---|---|---|
| PD1 reference index per λ | — | **0.6.45** |
| PD2 array surface + scalar guard | — | **0.6.46** |
| PD3 GD/GDD convention | — | **0.6.47** |
| PD4 `compute()` observable | — | **0.6.48** |
| F2.1 environment segment schema | 0.6.45 | **0.6.49** |
| F2.2 K assemblies, `residuals_multi` | 0.6.46 | **0.6.50** |
| F2.3 needle + LM joint | 0.6.47 | **0.6.51** |
| F2.4 Python `environments=` surface | 0.6.48 | **0.6.52** |
| F3.1 docs, examples, release | 0.6.49 | **0.6.53** |

**Why before Phase B and not after.** (b) is a wrong answer, not a missing
feature, and it is live in a **published** release — 0.6.44 went to PyPI
and crates.io on 2026-09-14. Phase B has not started, so shifting it costs
one table edit; deferring PD costs five releases of shipping a silent
error. The series is also short and self-contained, so running it
contiguously avoids a second re-versioning later.

**The alternative, if the F2 ladder is not to be touched:** PD1–PD4 become
0.6.50–0.6.53 and Phase B keeps 0.6.45–0.6.49. That is the only other
coherent option — interleaving PD items between F2 items is not, because
PD3's GD/GDD decision and F2.3's joint needle both touch the phase fold.

**Superseded recommendation.** The status check proposed landing a cheap
interim guard first and the real fix later. That is dropped: the
instruction is to fix (b) properly, and an interim guard that PD1
obsoletes within one release is two bumps and two CI cycles of churn. The
*permanent* guard — the one that still has a job after PD1 — is the
scalar-door check, and it lives in PD2 where it belongs.

---

## 1. The findings this plan acts on

### 1.1 (b) — the reference index is frozen at the centre wavelength

The reference is defined (`docs/spectralweave-target-kinds.md` §Differential
phase) as

```
Δφ(λ) = arg t(λ) − passes · 2π · n_inc · D · cosθ_inc / λ
```

`reference_phase` itself is evaluated **per wavelength**
(`optics_core.rs:49`) — the formula is not the problem. The problem is its
`n_inc_re` argument, which is a scalar baked into the data structure:

```rust
// merit.rs — SimCurves
pub total_d: f64,
pub n_front_re: f64,   // one number for the whole band
pub n_back_re: f64,
```

filled at `evaluator.rs:209-214`:

```rust
let (total_d, n_front_re, n_back_re) = if self.spec.uses_differential() {
    let total_d: f64 = sa.thicknesses.iter().sum();
    let n_front_re = sa.n_stack_cache[(nw / 2) * nl * 2];
    let n_back_re = sa.n_stack_cache[(nw / 2) * nl * 2 + (nl - 1) * 2];
    (total_d, n_front_re, n_back_re)
} else {
    (0.0, 1.0, 1.0)
};
```

`n_stack_cache` is laid out per wavelength, two entries (re, im) per layer
(`needle_operator.rs:997`), so `(nw / 2) * nl * 2` is **layer 0's real
index at the centre wavelength** — applied across the entire band.

**The error is exactly**

```
err(λ) = 2π · D · cosθ · [ n(λ) − n(λ_centre) ] / λ
```

independent of the engine, and therefore checkable by hand:

| Incidence medium | Band (nm) | D | max \|err\| |
|---|---|---|---|
| Air | any | any | **exactly 0** |
| N-BK7 | 400–700 | 500 nm | 0.097 rad (5.6°) |
| N-BK7 | 400–700 | 3 µm | 0.581 rad (33.3°) |
| N-BK7 | 380–780 | 3 µm | 0.825 rad (47.2°) |
| N-SF11 | 400–700 | 3 µm | 2.537 rad (**145.3°**) |
| N-BK7 | 600–1000 | 10 µm | 0.578 rad (33.1°) |

The air row is why this has not bitten: the overwhelmingly common
incidence medium is non-dispersive, and there the approximation is exact.
It bites when light enters through the substrate (a coating on the inner
face), in immersion or cemented optics, and in prism coupling — none of
which are exotic.

**Status today: documented, unguarded.** Two comments call dispersive
ambients "pathological" (`evaluator.rs:206-207`,
`spectralweave-target-kinds.md:184-185`). There is **no runtime warning,
no refusal, and no test pinning the approximation**. A user with a
dispersive layer 0 gets a wrong number and no signal.

**The sites that read the scalar** — the exhaustive list PD1 must move:

| # | Site | What it does |
|---|---|---|
| 1 | `merit.rs` `SimCurves.n_front_re` / `n_back_re` | the fields themselves |
| 2 | `merit.rs:332-346` `impl Default for SimCurves` | defaults `1.0` / `1.0` |
| 3 | `evaluator.rs:209-225` | the fill (centre-λ pick) |
| 4 | `merit.rs:918-932` | merit inner loop: subtracts `reference_phase` |
| 5 | `needle_pass.rs:178-192` `sample_op_value` | op-point Δφ |
| 6 | `needle_pass.rs:327-335` | hoists `n_inc` per demand (front/back pick) |
| 7 | `needle_pass.rs:456-461` | pointwise `dM/dD` gain shift via `reference_wavenumber` |
| 8 | `needle_pass.rs:711-715` | integral-demand gain shift, same |
| 9 | `navette-py/src/synthesis_merit.rs:63-81` | `SimCurves` PyO3 ctor (`n_front`, `n_back`) |
| 10 | `navette-py/src/synthesis_merit.rs:244-255` | `reference_rotation` FFI (`n_inc: f64`) |
| 11 | `src/navette/synthesis/__init__.py:162` | `apply_reference_rotation` |
| 12 | `src/navette/synthesis/__init__.py` `sim_curves_from_arrays` | `n_front` / `n_back` |
| 13 | `src/navette/_smatrix.pyi:321-375` | the stubs for 9, 10 (gated by `check_pyi_sync.py`) |

Sites 1–8 are PD1. Sites 9–13 are PD2.

**One consistency question to settle during PD1, not before.**
`needle_pass.rs:184` reads `n_inc.unwrap_or(1.0)` where `merit.rs:919`
reads `sim.n_front_re` directly. If the `None` arm is genuinely reachable
with a differential demand, the two paths disagree (air reference vs stack
reference) and the fold/merit twins would not catch it, because both are
computed from the same `sim`. Resolve it explicitly: either prove `None`
is unreachable when `differential_passes.is_some()` and make it an
`expect` with that reason, or make both arms take the same branch. Record
the answer in PD1's CORRECTIONS block either way.

### 1.2 (a) — no `compute()` path

`Request` (`smatrix.py:84`, `IntFlag`; natively `u64` constants from
`core_engine.rs:53`) is a 49-bit observable mask. It carries

- `PHI_TS` / `PHI_TP` (bits 27, 28) — **absolute** transmitted phase,
- `TS_C` / `TP_C` (bits 31, 32) — complex transmitted amplitudes,

and **no differential-phase bit**. `compute()` returns exactly the keys
the request implies (`expected_keys`, `smatrix.py:201`), so Δφ is not
obtainable from the simulation surface at all. The Python docs call
`PDts`/`PDtp` "synthesis-only quantities"
(`src/navette/synthesis/__init__.py:63`).

A simulation user who wants Δφ today must import from
`navette.synthesis`, call `compute(TS_C | TP_C)`, call
`apply_reference_rotation(...)`, then `np.angle` — assembling a
simulation quantity out of a synthesis helper, and supplying `total_d`
and `n_inc` by hand from a stack the engine already holds.

Bits 49–63 of the `u64` mask are free.

### 1.3 The coupling PD3 exists for

`spectralweave-target-kinds.md` states that group delay and GDD over Δφ
"come for free (finite differences kill the reference anyway)". That is
true **only because the index is frozen**: a constant-`n` reference
contributes a constant group delay and *exactly zero* GDD. Make the
reference dispersive — which is precisely what PD1 does — and

- GD acquires a λ-dependent shift, and
- **GDD acquires a genuine contribution** it does not have today.

So PD1 does not merely correct a number; it changes what "GD/GDD over Δφ"
means. That is a convention decision with a user-visible consequence and
it gets its own item rather than riding silently inside PD1.

### 1.4 What is **not** wrong

Checked, correct, and out of scope — recorded so a later reader does not
re-open them:

- **The formula and the `passes` factor.** `reference_phase` is per-λ and
  hand-verified (`optics_core.rs:453-469`).
- **The sign convention.** `+kD`, pinned by
  `solver_propagation_sign_matches_reference`; self-consistent across
  phase demands, needle `P_PHI` and GD/GDD. Only textbook-imported target
  numbers need conjugating, and the doc says so.
- **The choice of "equivalent layer of incidence medium"** as the
  reference. It is a definition, it is documented, and it is not this
  plan's business to change.
- **The native target machinery.** `targets.rs:302` maps
  `PDts → (Ts, s, 1)` and `PDtp → (Tp, p, 1)`; `ingest_spectral` /
  `ingest_angular` force `mode = "phase"`; `MeritTarget.differential_passes`
  carries it; polarization mismatch and `phase=False` are refused. All
  correct.
- **The needle gain shift** `dM/dD = Σ −2·kz·w·(s−rt)`, and the fact that
  it is uniform in z and so never moves the needle site. Correct — but it
  reads the frozen `n`, so it is in PD1's site list.
- **The two equivalent evaluation paths** (native subtraction vs numpy
  rotation) agree to 1e-12, pinned by the `oracle` check of
  `validation/parity/synthesis/test_differential_phase.py`.

---

## 2. PD1 — the reference index becomes wavelength-dependent (0.6.45)

**DONE:** — done (2026-09-14, feature commit pending). The reference
indices are per-λ rows (`Arc<[f64]>`, length `wavelengths.len()` or 1,
broadcast) filled from the whole `n_stack_cache` column; the merit inner
loop samples the row with the complex rows' own two-pointer machinery,
the needle op-point and both gain-shift sites interpolate via
`interp_n_row` at the demand wavelength. The plan's length-1 broadcast,
length-refusal (`reference_length_issue`, naming both numbers) and
 licence items 1–4 landed as written.

CORRECTIONS:

1. §7 decision 1 (PD before Phase B) is **resolved by instruction**:
   implementing the plan adopts its own §0.4 recommendation (0.6.45–
   0.6.48; Phase B slides to 0.6.49–0.6.53).
2. §7 decision 2 (`n_inc.unwrap_or(1.0)` reachability): **`None` is
   unreachable with a differential demand** — both `sample_op_value`
   call sites match `(Some(sim), Some(rows))` before calling and `n_inc`
   derives from that same `current_sim`, so the air-reference arm could
   never fire. Replaced by `expect` with the proof in its message (the
   first arm of the plan's §1.1 resolution menu).
3. The broadcast row needs an explicit branch in every reader (the
   merit loop's aligned path would index `row[offset + i]` out of
   bounds on a length-1 row): `sample_n` returns `row[0]` for the
   broadcast shape, and `interp_n_row` does the same before delegating
   to `interp_clamped`. Not in the plan's site list because the plan's
   sketch assumed a single shape — the broadcast rule made it two.
4. The engine fill produces full-length rows even for a constant
   medium (the plan's fill snippet); length 1 exists only through the
   ctor/default. The bit-exactness argument covers both (a constant
   row interpolates to its one value exactly), verified bitwise.
5. The engine-fill test pinned the NEW shape as well: the
   non-differential branch now yields `[1.0]` (one element) and the
   differential branch length-`nw` columns — the pre-PD1 test asserted
   scalars and was updated with the per-λ expectation.

Gates: fmt, clippy -D warnings, 550 lib tests (564 with `lm`/all
features), pytest 766 passed, five tools, ten harnesses, both
fingerprints (per-platform needle digests reproduce; CI evaluates the
Linux side). Release wheel 0.6.45 built and installed.

### What the tree has today

Sites 1–8 of §1.1's table. One `f64` per side, picked at `nw / 2`.

### Change

**The fields become per-λ rows, shaped like every other row in the
struct:**

```rust
pub struct SimCurves {
    pub angles: Arc<[f64]>,
    pub wavelengths: Arc<[f64]>,
    pub total_d: f64,
    /// Real incidence index per wavelength (differential-phase reference).
    /// Length is either `wavelengths.len()` or 1 (a non-dispersive medium
    /// broadcasts). Default `[1.0]`.
    pub n_front_re: Arc<[f64]>,
    /// Real exit index per wavelength; same shape rule. Default `[1.0]`.
    pub n_back_re: Arc<[f64]>,
    ...
}
```

**Length 1 broadcasts.** This is the load-bearing design decision of the
item and it buys three things at once: `Default` stays `[1.0]` so
`total_d = 0` still reproduces absolute phase bit-for-bit; the air case
allocates one `f64` and takes the same branch it takes today; and PD2's
scalar FFI door keeps working without a shim. A length that is neither 1
nor `n_wavs` is a refusal at construction, naming both numbers.

**The fill** (`evaluator.rs:209-225`) collects the whole column instead of
one element:

```rust
let n_front_re: Arc<[f64]> =
    (0..nw).map(|w| sa.n_stack_cache[w * nl * 2]).collect();
let n_back_re: Arc<[f64]> =
    (0..nw).map(|w| sa.n_stack_cache[w * nl * 2 + (nl - 1) * 2]).collect();
```

`total_d` is unchanged — the ambient and substrate carry zero thickness,
so the plain sum is still the coating thickness.

**The readers** (sites 4–8) gain an index. The demand loop already walks
wavelengths, so this is an index into a slice, not a search. The
front/back pick stays exactly where it is (`key.curve.is_back()`), it just
selects a slice instead of a scalar.

For the two gain-shift sites (7, 8), `reference_wavenumber(twl_i, n, …)`
takes the index **at the demand wavelength `twl_i`**, not at the sim grid
index — the demand grid and the sim grid need not coincide, which is why
the existing code interpolates the rows. Use the same
`interp_clamped` the rows use, against `sim.wavelengths`. Getting this
wrong is the one way to make the `dM/dD` twin fail while merit still
looks right, so the FD twin below is the gate that matters.

### The licence — the exhaustive list of what may change

1. Δφ for a spec with a differential demand **and** a wavelength-varying
   incidence (or exit) index. By `err(λ)` of §1.1, up to 145° in the
   N-SF11 case.
2. Merit, residuals and LM trajectories for such a spec — a consequence
   of 1, not a separate change.
3. The needle gain shift `phi_gain_shift` for such a spec, and with it the
   *predicted* gain bookkeeping. The needle **site** must not move: the
   shift is uniform in z and stays uniform when `kz` varies with λ,
   because the uniformity is in z, not in λ.
4. `SimCurves`'s field types and its constructor's arity (an API change,
   crate-internal until PD2 exposes it).

**An eighth difference is a defect.** Nothing else may change — in
particular nothing at all may change for a non-dispersive incidence
medium, which is the bit-exactness argument below.

### Bit-exactness

Ground rule 5. **PD1 may move neither fingerprint.** Verified reachable as
a gate: neither the random-stack differential test nor the needle-run pin
carries a differential-phase demand (`test_needle_pin.py`'s "phases" are
pipeline phases, not optical phase), so both must reproduce exactly.

Stronger, and the gate to actually lean on: **for a non-dispersive
ambient the whole differential path must be bitwise unchanged.** With
`n_front_re = [n]` the arithmetic reduces to today's
`reference_phase(λ, n, …)` with the identical `f64` operand in the
identical order, so "bitwise" is the honest requirement, not "within
tolerance". A twin that asserts `==` on the merit value for an air-ambient
PD spec, before and after, is the cheapest proof and it is exact.

### Twins

Extending `validation/parity/synthesis/test_differential_phase.py` (367
lines, nine checks today — `hand`, `oracle`, `numpy-2x2`, `dM/dD`, `fold`,
`zero-D`, `errors`, `fold-equiv`, `angular`):

1. **`dispersive-hand`** — a three-point band with a hand-built dispersive
   `n_front`, Δφ against hand arithmetic per point. The one twin that
   would have caught the defect.
2. **`nondispersive-bitwise`** — air ambient, PD spec, merit and residuals
   `==` (not `approx`) the values the same spec produced at 0.6.44.
   Recorded as literals in the test, so it keeps working after PD2 and
   PD4 move the surface.
3. **`dM/dD` extended** — the existing FD-slope check re-run with a
   dispersive index, proving the gain shift interpolates the index at the
   demand wavelength rather than the sim index.
4. **`needle-site-invariance`** — same design, dispersive ambient, needle
   site `argmax` identical with and without the differential demand.
   Pins licence item 3.
5. **`length-refusal`** — `n_front_re` of length 2 against a 5-wavelength
   grid refuses, naming both lengths, ASCII (ground rule 7).
6. Rust unit: `reference_phase` unchanged (it is — this item does not
   touch `optics_core.rs`), plus a `SimCurves` construction test for the
   broadcast rule.

### Gates

The full battery of `implementation_plan.md` §0.1 rule 3, five tools
included, plus:

- both fingerprints reproduce (ground rule 5, per platform per C1);
- `nondispersive-bitwise` passes as `==`;
- the nine existing differential-phase checks still pass unchanged;
- `cargo clippy --workspace --all-targets -- -D warnings` — the
  `Arc<[f64]>` change will surface `needless_range_loop` or
  `clippy::type_complexity` if the indexing is written carelessly; fix,
  do not `allow`.

---

## 3. PD2 — the array surface at the FFI and Python doors (0.6.46)

### Change

Sites 9–13 of §1.1. Each door that takes a scalar index learns to take an
array, keeping the scalar as the broadcast case:

- **`SimCurves` PyO3 ctor** (`synthesis_merit.rs:63-81`): `n_front` /
  `n_back` accept `float | FloatArray`. A float becomes length 1.
- **`reference_rotation`** (`synthesis_merit.rs:244-255`): `n_inc`
  accepts `float | FloatArray`, length 1 or `len(wavelengths)`.
- **`apply_reference_rotation`** and **`sim_curves_from_arrays`**
  (`src/navette/synthesis/__init__.py`): same, passed straight through —
  these stay thin, the validation stays native.
- **`_smatrix.pyi:321-375`**: stubs updated. `check_pyi_sync.py` is
  blocking, so this is not optional and not deferrable.

### The guard — the part that outlives PD1

After PD1 the native path is correct by construction, because it reads the
stack. The doors above are the remaining way to supply a *scalar* index
for a stack that is actually dispersive, and a user assembling
`SimCurves` by hand has no engine to catch them.

So: when a scalar (or length-1) index is supplied **and** a differential
demand is active, the door emits a warning naming the quantity, the
supplied value, and the remedy — pass the per-λ array. Ground rule 6 —
announced, never silent. Ground rule 7 — ASCII.

It is a warning, not a refusal: air is a legitimate scalar and by far the
common case, and the door cannot know the stack it was not given. A
refusal would break every correct existing caller.

### Twins

- Scalar and length-1 array produce **bitwise identical** results at every
  door.
- A dispersive array through `apply_reference_rotation` equals the native
  differential demand to 1e-12 — the existing `oracle` check, re-run with
  a per-λ index, which is what keeps the numpy path a real oracle instead
  of a co-drifting copy.
- The warning fires once, carries the remedy, and is ASCII.
- Length mismatch refuses at the door with both lengths named.

### Gates

Full battery. `check_pyi_sync.py` and `check_exposure.py` are the two that
will actually catch mistakes here; run them before the Rust suites, they
are seconds.

---

## 4. PD3 — GD/GDD over Δφ: decide the convention and pin it (0.6.47)

### The problem

§1.3. Before PD1, the differential reference contributes a constant group
delay and zero GDD, which is what makes the current documentation true.
After PD1 it contributes a λ-dependent group delay and a non-zero GDD.
Users computing GD/GDD over Δφ get different numbers, from a change whose
stated purpose was a *different* quantity.

### The three options

1. **Report GD/GDD over the corrected Δφ.** Physically honest: the
   reference layer really is dispersive and really does have a GDD. Users
   comparing against a measured group delay through an immersion medium
   get the right answer. Changes numbers for dispersive ambients.
2. **Report GD/GDD over `arg t` (absolute), documenting that the
   differential reference is not applied to dispersion orders.** Keeps
   today's numbers; makes Δφ and its own derivatives inconsistent, which
   is a trap.
3. **Report both, under distinct keys.** Most information, widest surface.

### Decision

**Option 1.** The series exists to make the reference correct; carving out
the dispersion orders would reintroduce, one level up, exactly the
inconsistency being removed — and option 2's inconsistency is the kind
that is discovered years later by someone chasing a factor nobody can
explain. Option 3 is deferred, not refused: if a real case wants the
absolute-phase orders alongside, it is an additive key and can be added
without revisiting this decision.

### Change

- `spectralweave-target-kinds.md`: replace the "finite differences kill
  the reference anyway" sentence. It was true, it is about to stop being
  true, and it must not survive as a stale claim — this is the same class
  as amendment 3 §4.2's stale `Param` comment.
- CHANGELOG: state the numeric consequence for dispersive ambients.
- A twin that pins the decision so it cannot silently regress.

### Twins

- GD over Δφ with a dispersive reference equals the analytic
  `d(ref)/dω` correction applied to GD over `arg t` — hand-computed for a
  two-point case.
- GDD over Δφ is **non-zero** for a dispersive reference and **zero** for
  a non-dispersive one. Two asserts, and they are the whole decision.
- Non-dispersive: GD/GDD bitwise unchanged from 0.6.46.

### Gates

Full battery. This item is mostly a decision, a doc edit and three
asserts; the risk is that it gets skipped because it looks like docs. It
is not docs — it is the item that keeps PD1 from quietly changing a
second quantity.

---

## 5. PD4 — `PDts`/`PDtp` as first-class `compute()` observables (0.6.48)

### Why this is small, once PD1 has landed

The objection to a `compute()`-level differential phase used to be that
`compute()` has no `total_d` and no `n_inc`. It does — it just never
looked: `ScatterMatrix` holds the stack, so the coating thickness is the
sum of the interior thicknesses (ambient and substrate carry zero, the
same rule `evaluator.rs` relies on) and the incidence index is layer 0's
real part **per wavelength**, which is exactly the column PD1 taught the
evaluator to collect. No new user plumbing; the same two derivations, at a
second call site.

### Change

- Two request bits, `PD_TS = 1 << 49` and `PD_TP = 1 << 50` (bits 49–63
  of the `u64` mask are free), with `REQ_PD_TS` / `REQ_PD_TP` in
  `core_engine.rs` and the `Request` IntFlag in `smatrix.py`.
- Emitted keys `PDts` / `PDtp`, wired into `expected_keys`.
- Derivation: the complex forward-t rows the engine already computes for
  `TS_C`/`TP_C`, minus `reference_phase(λ, n_front_re[λ], θ, D, 1)`.
- A convenience view alongside `complex_amplitudes()` and `dispersion()`:
  `differential_phase(*, s_pol=True, p_pol=True)`.
- `.pyi` stubs; `check_pyi_sync.py` blocking.

### Two decisions this item must state, not discover

1. **Wrapped or unwrapped.** `arg()` wraps to `(−π, π]`; `reference_phase`
   does not (`evaluator.rs:773` notes the 2.5π case). Decision: **emit the
   wrapped principal value**, consistent with `PHI_TS`/`PHI_TP`, and say
   in the docstring that unwrapping across the grid is the caller's job.
   A user-facing observable that silently differs in convention from its
   absolute-phase neighbour is worse than one that needs `np.unwrap`.
2. **Incoherent stacks.** `dispersion()` already carries the "physically
   meaningful only for coherent stacks" caveat. Δφ inherits it verbatim —
   copy the caveat, do not paraphrase it.

### Twins

- `compute(PD_TS)` equals `apply_reference_rotation` on
  `compute(TS_C)` followed by `np.angle`, to 1e-12 — the same oracle
  relation §1.4 already relies on, now spanning the new surface.
- `compute(PD_TS)` equals the native merit's op-point Δφ for the same
  stack and grid, to 1e-12 — this is what proves the two surfaces agree
  rather than merely both existing.
- Dispersive ambient: both of the above still hold (they would not, with
  the frozen index — which is the concrete reason PD4 is sequenced last).
- `expected_keys` round-trip for the two new bits, and the existing
  request-bit sync test still green.

### Gates

Full battery, plus `validation/smoke/test_request_bits.py` — the request
mask has a dedicated smoke test and two new bits must appear in it.

---

## 6. Execution order

Per amendment 3 §7.1, each item is: code commit → battery green → seven
version sites → CHANGELOG → this file's item marked **DONE** with a
CORRECTIONS block if reality differed → docs commit → **push** → wait for
green CI → next item.

1. **PD1** (0.6.45) — the core fix. Blocks everything else.
2. **PD2** (0.6.46) — the doors and the surviving guard.
3. **PD3** (0.6.47) — the GD/GDD convention PD1 disturbed.
4. **PD4** (0.6.48) — the simulation surface, last by design.
5. Phase B resumes at **0.6.49** (F2.1), per §0.4's ladder.

The master-item table of `implementation_plan.md` §0.3 gains four rows and
the Phase B rows are re-versioned; both edits ride PD1's docs commit. Per
review C/E1, **count the cells** — every master-table row is exactly seven
cells, and the last two table edits in this repo joined two rows into one.

## 7. Open decisions

| # | Decision | Status |
|---|---|---|
| 1 | PD before Phase B (0.6.45–0.6.48) vs after (0.6.50–0.6.53) | **Resolved — before, adopted by instruction** (implementing this plan takes its §0.4 recommendation). |
| 2 | `n_inc.unwrap_or(1.0)` at `needle_pass.rs:184` — reachable with a differential demand? | **Resolved during PD1:** `None` is unreachable when a differential demand is live (`sample_op_value` runs only under `(Some(sim), Some(rows))`, and `n_inc` derives from the same sim); now an `expect` with the reason (PD1's CORRECTIONS 2). |
| 3 | GD/GDD over the corrected Δφ (option 1) | **Decided** — PD3 §Decision. Option 3 deferred until a real case. |
| 4 | Wrapped principal value for `PDts`/`PDtp` keys | **Decided** — PD4. |
| 5 | Scalar index at the FFI doors: warn, not refuse | **Decided** — PD2. |
| 6 | Back-side differential labels (`PDrs`, `PDtbs`, …) | **Out of scope.** `n_back_re` is carried and now correct, but no label maps to it; adding labels is a separate item with its own `passes = 2` question. |
