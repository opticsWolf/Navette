# Implementation review — the PD series (0.6.45–0.6.48)

Review of `82998a7`, `12ed85a`, `8d400c6`, `fab63d1` (plus their four docs
commits) against `docs/implementation_plan_pd.md`. Reviewed at
`54af967`, branch `dev_feature`.

**Verdict: the series does what the plan said it would.** Both defects
are fixed, the licence held, the bit-exactness demand held, and every
twin the plan named exists and passes. Five findings follow; **G1 is the
one that should not ship** — it is a false-positive warning in the exact
item whose entire deliverable is that warning.

**Nothing is released yet.** `main` is still at `a095f60` (0.6.44) and
the only tag is `v0.6.44`. All four items live on `dev_feature`. Every
finding below is therefore fixable before the release tag, at no cost to
users.

---

## 0. Finding namespace

Findings are **G1–G5**. Namespaces already spent: `A*`, `N*`, `U*`,
`R*`, `F*`, `§D` (main plan), `B*` (amd 2), `C*` (amd 3), `V*` (amd-3
ledger), `E*` (review C), `PD*` (the PD plan). `G*` is verified free
across `docs/`. `P*` remains unusable (P0–P3 are priority labels).

---

## 1. What was verified, and how

Everything in this section was re-derived independently — the oracle
arithmetic was written from `spectralweave-target-kinds.md`'s definition,
not imported from the repo's own twins — and then measured against the
built 0.6.48 extension.

### 1.1 Gates (all green, all re-run)

| Gate | Result | Plan's claim |
|---|---|---|
| `cargo test --workspace` | 552 lib + 22 parity + 15 py | 552 ✓ |
| `cargo test -p navette --features opt-minpack-lm` | 557 | — |
| `cargo test -p navette --features opt-argmin` | 561 | — |
| `cargo test -p navette --all-features` | 566 | 566 ✓ |
| `pytest` | 771 passed, 1 skipped | 771 ✓ |
| `cargo fmt --all --check` | clean | ✓ |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean | ✓ |
| Five `tools/check_*.py` | all PASS | ✓ |
| Ten `validation/review/*.py` | all PASS | ✓ |
| Both fingerprints (needle pin, roundtrip) | 18 passed | ✓ |
| Seven version sites | all `0.6.48` | ✓ |

### 1.2 PD1 — the fix reaches every site

All eight core sites of the plan's §1.1 table moved. The readers sample
at the **demand** wavelength, which is the part that could have been got
wrong silently:

- `merit.rs:940` `sample_n` — shares the two-pointer with `sample_c`.
  Sound: `sample_c` runs first in the `Phase` arm and advances the index
  to the bracket for `t_wl[i]`; `sample_n`'s `while` then finds it
  already reached. `offset` is a **wavelength**-grid offset
  (`sim_wl.partition_point`), not an angle offset — the complex row is
  pre-sliced per angle at `merit.rs:850` — so `row[offset + i]` is
  in-bounds for a length-`n_wav` reference row.
- `needle_pass.rs:794` `interp_n_row` — length-1 short-circuits, longer
  rows go to `interp_clamped` against `sim.wavelengths`.
- Both gain-shift sites (`:484`, `:741`) interpolate at `twl_i`.

Measured: with a frozen centre-λ index vs the per-λ row, the merit
differs on **every** grid regime — aligned, offset, coarse, and clamped
past both ends. The row is genuinely reaching the arithmetic, not
sitting unread.

The `n_inc.unwrap_or(1.0)` → `expect` change (plan §7 decision 2) is
sound. Both `sample_op_value` call sites (`:402`, `:650`) sit inside
`(Some(sim), Some(rows))` arms and `n_inc = current_sim.map(...)`, so
`current_sim.is_some()` ⟹ `n_inc.is_some()`. The `None` arm really was
unreachable.

### 1.3 Bit-exactness — proven, by a different route than the plan's

The plan demanded that a non-dispersive ambient be **bitwise** unchanged
and offered `nondispersive-bitwise` (recorded 0.6.44 literals) as the
proof. I could not reproduce those literals through my own oracle — my
`oracle_tf` differs from the test's at ULP level, which moves the
*target* by ~1e-15 rad and so moves the residual. That makes the literal
check unverifiable from outside without rebuilding the 0.6.44 tree.

So I proved the same claim a different way. Pre-PD1 the reference was one
`f64` used at every point; post-PD1 a **length-1** row takes the `row[0]`
branch, which *is* that scalar arithmetic unchanged. If a length-1 row
and a constant **full-length** row agree bitwise, the full-row path is
bitwise the pre-PD1 path:

| Target grid | merit | residuals | |
|---|---|---|---|
| aligned (target == sim grid) | `==` | `==` | bitwise |
| non-aligned, offset +7.3 nm | `==` | `==` | bitwise |
| non-aligned, coarse 3-point | `==` | `==` | bitwise |
| clamped below and above the grid | `==` | `==` | bitwise |

This is **stronger** than the literal check in one respect — it exercises
the interpolating and clamping branches, which the literal check's
aligned grid never reaches — and weaker in another: it cannot rule out
both paths having moved together. Taken with the two untouched
fingerprints, the licence holds.

### 1.4 PD2 — the doors

Measured at `reference_rotation`: a per-λ array passes silently; a
length-2 array against a 3-point grid refuses with
`reference_rotation: n_inc has length 2, but the grid has 3 wavelengths -
supply length 1 (constant medium) or 3` — both numbers named, ASCII
(ground rule 7). The scalar guard fires, carries the remedy, and is
ASCII. See **G1** for what is wrong with *which* scalar it reports.

### 1.5 PD3 — a decision, honestly scoped

No production code changed, and the item says so. The stale
"finite differences kill the reference anyway" sentence is gone from
`spectralweave-target-kinds.md:176-177`, replaced with the analytic
corrections. The engine's own `GD`/`GDD` keys are over **absolute**
phase and were not touched — so the item's practical content is the
decision, the doc, and three asserts. That is exactly what the plan
specified, and the CORRECTIONS block is candid about the two places
reality was messier (the GDD *correction* is what vanishes for a
constant ambient, not the coating's own GDD; the hand-exact construction
needed film index == substrate index to kill the Fabry–Pérot
denominator). Both are the kind of trap worth having written down.

### 1.6 PD4 — measured against the definition, not against itself

Independent oracle, dispersive ambient, 4 angles × 5 wavelengths, both
polarizations:

| Check | Result |
|---|---|
| `PDts`/`PDtp` vs `arg(t) − 2πn(λ)D cosθ/λ` | max err **1.8e-15** |
| same, single angle | max err **8.9e-16** |
| what the frozen index would have given | off by **0.226 rad** |
| wrapped principal value in (−π, π] | in range |
| `D` excludes the ambient (d = 0 vs d = 999 nm) | **identical** |
| `D = 0` → differential ≡ absolute | **bitwise** |
| merit op point vs `compute(PDts)`, dispersive | **0.0e+00, bitwise** |
| `compute(PD_TS)` alone emits | `['PDts']` only — no `ts_c` leak |

The last row is CORRECTIONS 1 working: a PD-only request does not leak
the amplitude buffers. The op-point agreement at exactly zero is the
result of deriving both through the same wrap expression, as designed.

The coherent-stacks caveat is copied **verbatim** from `dispersion()`,
as the plan required, and `differential_phase()` has the same posture as
its neighbour — docstring only, no runtime gate. Confirmed by running
`compute(PD_TS)` on a stack with a 50 µm incoherent layer: it returns
finite numbers, silently, with `D` including that layer. That is the
caveat doing its job rather than a defect, and it is recorded here so it
is not rediscovered as one.

---

## 2. Findings

### G1 — the PD2 guard warns about a side no demand can read (**P0, fix before release**)

`warn_scalar_reference` (`synthesis_merit.rs:354`) gates on
`spec.uses_differential()` — a **spec-wide** boolean — and then reports
*every* reference row of length 1:

```rust
if !spec.uses_differential() { return Ok(()); }
let nf = sim.n_front_re.as_ref();
let nb = sim.n_back_re.as_ref();
if nf.len() != 1 && nb.len() != 1 { return Ok(()); }
let sides = [("n_front", nf), ("n_back", nb)]
    .into_iter()
    .filter(|(_, row)| row.len() == 1)
```

But `n_back_re` is read only under `key.curve.is_back()`
(`merit.rs:996`, `needle_pass.rs:354`), and `targets.rs:302`
`differential()` maps only `"PDts" => ("Ts", "s", 1.0)` and
`"PDtp" => ("Tp", "p", 1.0)` — **both front**. No label reaches the back
reference. The plan says so itself (§7 decision 6: "`n_back_re` is
carried and now correct, but no label maps to it").

So the back half of the warning is unconditionally a false positive.
Measured, front-only `PDts` demand, hand-built `SimCurves`:

| `n_front` | `n_back` | warning |
|---|---|---|
| scalar | default | `(n_front=1, n_back=1)` — half of it spurious |
| **per-λ array** | default | **`(n_back=1)` — entirely spurious** |
| per-λ array | per-λ array | silent |

The middle row is the problem. A user who does exactly what the guard
asks — supplies the per-λ front reference — is still warned, about an
index nothing reads, and the only way to silence it is to hand over a
full-length row for a quantity no label can consume. `n_back_re`
defaults to `[1.0]`, so this fires on every hand-assembled `SimCurves`
with a differential demand, forever.

This matters more than its size suggests. PD2's whole deliverable is
this warning; the plan called it "the part that outlives PD1". A guard
that cries wolf at correct callers is the one failure mode that
guarantees it gets ignored or suppressed — and then it is not a guard.

**Fix.** Report only the side(s) actually demanded. `MeritSpec` already
holds what is needed: walk `targets()` for `differential_passes.is_some()`
and collect `keys()[t.key_idx].curve.is_back()`. Front-only demands then
report `n_front` alone and stay silent on a default `n_back`. The
engine-fill path is unaffected (it produces full-length rows on both
sides and never warns either way).

Add a twin: front-only demand + per-λ `n_front` + default `n_back`
emits **zero** warnings. That is the case that is wrong today and no
existing check covers it.

### G2 — two `total_d` derivations, one claiming to be the other (**P1**)

`solver.rs:723` (PD4):

```rust
let total_d: f64 = self.thicknesses[1..self.n_layers.saturating_sub(1)]
    .iter().sum();
```

`evaluator.rs:212` (PD1, and pre-existing):

```rust
let total_d: f64 = sa.thicknesses.iter().sum();
```

PD4's comment says its `D` is "the sum of the interior thicknesses
(ambient and substrate carry zero — **the same rule the synthesis
evaluator relies on**)". It is not the same rule; it is a different
expression that happens to agree whenever the half-spaces are zero.

They can disagree. `DesignStack(ambient, substrate, films)` on the
Python surface (`synthesis_pipeline.rs:598`) takes arbitrary
`PyLayerSpec`s; `PyLayerSpec::new` rejects only negative and non-finite
thicknesses, and `with_films` (`structure.rs:152`) passes the ambient
through untouched. So a `LayerSpec("amb", nk, 100.0)` puts 100 nm at
slot 0, `solver_arrays()` writes it to `thicknesses[0]`, and the
evaluator's `D` gains 100 nm that the engine ignores for propagation.

The *pipeline* path is safe — `driver.rs:63` `fixed_half` hardcodes
`d_nm: 0.0` — so this is the direct-`DesignStack` door only, and the
evaluator half predates PD1. What PD4 newly created is the **second
surface** and the untrue equivalence claim. Verified that
`compute(PDts)` is invariant to an ambient thickness of 999 nm (PD4 is
the correct one); the merit path would not be, so the two would disagree
on precisely the stack PD4's own twin asserts they agree on.

**Fix.** One line: make `evaluator.rs:212` use `[1..nl - 1]` too. That
makes both expressions identical and the comment true. Optionally also
zero or refuse a non-zero half-space thickness in `with_films` — the
engine already warns about it at `solver.rs:1131`, but that warning does
not reach the `DesignStack` door.

### G3 — the master item table marks only PD1 done (**P2, bookkeeping**)

`implementation_plan.md:164` carries
`| ~~PD1~~ **DONE (0.6.45)** | … |`, matching the convention every other
finished item uses (F1.2…F1.5, C1, C2, C3). Lines 165–167 still read
plain `| PD2 |`, `| PD3 |`, `| PD4 |` — while the version ledger at
`:2754-2757` marks all four **done**. The two tables contradict each
other.

Cause: PD2's, PD3's and PD4's docs commits (`72cc061`, `ed3cd4a`,
`54af967`) each touched `implementation_plan.md` for exactly one line —
the ledger row — and never returned to §0.3. PD1's (`55ecdd1`) did both.

Trivial to fix, three lines. Worth naming because it is the same class
review C's E-series existed to close, and because §0.3 is the table a
reader consults first.

### G4 — `merit()` and `residuals()` now panic (**P2, API posture**)

`MeritSpec::merit` (`merit.rs:768`) and `::residuals` (`:790`) call
`panic!` on a malformed reference-row length. Both are `pub` on a crate
published to crates.io.

The choice is defensible and the comment explains it: the FFI ctor
already refuses at construction (`synthesis_merit.rs:130`), so this is
unreachable from Python; `merit` returns `f64` with no error channel,
and `residuals` returns `Result<(), CurveId>` whose error type cannot
carry a message. Failing closed beats mis-sampling silently.

Two things to tidy rather than redesign:

1. It is not in the public doc comments. A Rust caller hand-building
   `SimCurves` — the exact caller the guard is for — learns about it by
   crashing. One `# Panics` line on each.
2. `curve_sensitivity` (`:1272`) has no such guard, which reads as an
   oversight and is not: the reference is **additive and independent of
   the curve**, so it drops out of `d(residual)/d(curve)` and
   `curve_sensitivity_into` never reads the rows (`reference_phase`
   appears exactly twice in `merit.rs`, at `:576` in a doc comment and
   `:1000` in `residuals_into`). Worth one line saying so, or the next
   reader adds a guard that does nothing.

### G5 — PD4 is absent from the README (**P2, before the tag**)

`README.md` contains **zero** mentions of differential phase, `PDts`, or
`differential_phase`. PD4 adds a first-class `compute()` observable, two
request bits and a convenience view alongside `complex_amplitudes()` and
`dispersion()` — both of which the README does describe.

This is not a deviation: the plan's §5 never asked for a README edit. It
is flagged because `a095f60`, the commit immediately before this series,
existed specifically to close the gap between the README's feature
surface and the code — after twelve releases in which the README got
nothing but its version bump. PD4 reopens that gap by one item, on the
same branch, three commits later.

Ground rule N6 applies: targeted edit, never a blanket `sed`.

---

## 3. Incidental (pre-existing, outside this series)

`CHANGELOG.md:2934` has an unescaped `|n|` inside a two-column table
cell, splitting the row into four cells:

```
| refractive index with `|n|` past `sqrt(DBL_MAX)` | `n**2` overflowed inside the solve; NaN out |
```

GFM drops the surplus cells, so the row renders wrong today. From
`9cc1349` (R3.1), long before the PD series — reported only because the
cell-count check that review C's E1 recommended was run across all four
touched documents and this is the one hit. `implementation_plan.md` (14
tables), `implementation_plan_pd.md` (5) and `implementation_amendment_3.md`
(7) are all consistent.

Note for whoever builds that check for real under F3.1: it must split on
**unescaped** pipes only. A naive `split('|')` false-positives on
`max \|err\|`, which appears in the PD plan's own error table.

---

## 4. Disposition

| # | Finding | Priority | Effort | Where |
|---|---|---|---|---|
| **G1** | PD2 guard warns about `n_back`, which no demand reads | **P0** | S | `synthesis_merit.rs:354` + one twin |
| **G2** | Two `total_d` derivations; the comment claims they are one | P1 | S | `evaluator.rs:212` |
| **G3** | Master table §0.3 marks only PD1 done | P2 | S | `implementation_plan.md:165-167` |
| **G4** | `merit`/`residuals` panic, undocumented; sensitivity asymmetry unexplained | P2 | S | `merit.rs:763`, `:789`, `:1272` |
| **G5** | PD4 absent from the README | P2 | S | `README.md` |

**Recommended shape.** G1 is a behaviour fix and wants a version bump of
its own — `0.6.49`, which pushes Phase B to `0.6.50–0.6.54`. G2 is a
one-line correctness alignment and can ride it. G3, G4 and G5 are
bookkeeping and docs: one commit, no bump, the same shape as `138bdcb`.

If the preference is not to move the F2 ladder a third time, G1 and G2
can instead fold into the release commit that merges `dev_feature` to
`main` — nothing is published, so there is no user-visible history to
preserve, and 0.6.48 would ship correct. That is a call for the owner;
the ladder cost is the only difference.

---

## 5. What this review did not verify

Stated plainly so the gaps are known:

- **The `nondispersive-bitwise` literals were not reproduced.** §1.3
  explains why and what was proven instead. Confirming them needs a
  rebuild of the 0.6.44 tree, which would clobber the installed 0.6.48
  extension. The claim is independently supported, but not by the route
  the test takes.
- **The fingerprints were checked on Windows only.** Per C1 the needle
  digests are per-platform; the Linux side is CI's to evaluate, and CI
  has not run on `54af967` at the time of review.
- **No long-horizon synthesis run.** The licence's item 2 (LM
  trajectories for dispersive specs) was verified by construction — it
  follows from item 1 — not by running an optimization to convergence
  and comparing designs.

---

## 6. Applied

Every finding was re-verified against the tree before fixing; all five
stood as written (G1's three-row table reproduced exactly, including the
entirely-spurious middle row). Applied in two commits on `dev_feature`:

| # | Action | Where |
|---|---|---|
| **G1** | `MeritSpec::demanded_reference_sides()` (walk targets, collect `key.curve.is_back()` per differential demand); `warn_scalar_reference` reports only demanded sides. Front-only + per-λ front + default back is now **silent**; a scalar front is reported alone. Twins: `guard-sides` in the PD2 door twin + Rust `demanded_reference_sides_tracks_the_labels` (front/back/none, incl. the hypothetical back curve). | `32adee2` (0.6.49) |
| **G2** | `evaluator.rs` `total_d` now sums `[1..nl-1]`, matching the engine's PD keys; the comment claiming the two expressions were one rule is now true. Measured: merit op point vs `compute(PDts)` agree **bitwise** at half-space thicknesses 999/777 nm (they would have disagreed pre-fix); bitwise identical for legal stacks. | `32adee2` (0.6.49) |
| **G4** | `# Panics` sections on `MeritSpec::merit` and `::residuals`; the deliberate-no-guard note on `curve_sensitivity` (the reference is additive and independent of the curve, so it drops out of d(residual)/d(curve)). | `32adee2` (docs-only in substance; they ride the bump commit) |
| **G3** | §0.3 rows for PD2/PD3/PD4 struck (`~~PD2~~ **DONE (0.6.46)**` etc.), matching PD1's convention; every row verifies at 7 unescaped-pipe cells. | this commit |
| **G5** | §4 gains a **Differential Phase Observables** bullet (`PDts`/`PDtp` keys, the view, the per-λ reference, the wrap, the shared derivation, the coherent caveat) — targeted edit, not a blanket sweep. | this commit |
| incidental | `CHANGELOG.md`'s `\|n\|` cell now escapes its pipe; the row renders at 2 cells again. | this commit |

The review's recommended shape was followed: G1 took its own bump
(`0.6.49`, G2 riding it), Phase B slides to `0.6.50–0.6.54` (§0.4's table
restamped: Was = amendment 3's assignment, Becomes = today's), and
G3/G4/G5 + the incidental landed as one no-bump bookkeeping commit — the
`138bdcb` shape.

CI: the review was written before CI ran on `54af967`; that run
(34898484657) has since completed green on all four jobs, which closes
the Linux-fingerprint gap above. What remains open is the first bullet:
the `nondispersive-bitwise` literals are still not reproduced from
outside (rebuilding the 0.6.44 tree would clobber the installed
extension; §1.3's route carries the claim independently).

CORRECTIONS:

1. G2's optional half — "zero or refuse a non-zero half-space thickness
   in `with_films`" — was **not** done. With the interior-sum alignment
   the two surfaces agree on any stack, so the disagreement the option
   guarded against no longer exists; the engine's own warning at
   `solver.rs:1131` still fires. A door-level refusal remains future
   work if a stack with a thick half-space is ever judged worth
   refusing rather than warning about.
2. The disposition's "one commit" for the bookkeeping items became two
   files' worth of Rust doc comments landing in the *bump* commit
   instead — they are code-adjacent (fmt/clippy apply to them), and
   splitting them out would have made the docs commit compile-unclean
   for no reader benefit. Recorded because the shape differs from the
   review's table.
3. The review's gate table reports feature-specific lib-test counts
   (557 / 561 / 566); post-application the counts are 553 lib / 567
   all-features (the new `demanded_reference_sides` test). All green;
   noted so the next reader does not treat the drift as a regression.

---

## 7. Second round — review of the applied fixes

Reviewed at `b86998f`, branch `dev_feature`, against a **fresh
`maturin develop --release` build** of that tree (the extension in the
working copy predated the last `merit.rs` edit, so every behavioural
claim below is measured against a rebuild, not against the diff).

**Verdict: all five findings are genuinely fixed, and one of the two
fixes was not pinned.** Findings are **H1–H4**; `H*` is verified free
across `docs/`.

### 7.1 What the fixes actually do

**G1.** Re-measured on the rebuilt extension, the three-row table that
made the case now reads:

| `n_front` | `n_back` | before | after |
|---|---|---|---|
| scalar | default | `(n_front=1, n_back=1)` | `(n_front=1)` — demanded side alone |
| per-λ array | default | **`(n_back=1)` — entirely spurious** | **silent** |
| per-λ array | per-λ array | silent | silent |

`demanded_reference_sides()` mirrors its consumers exactly: the
row-selection branches at `merit.rs:1032` (`residuals_into`) and
`needle_pass.rs:353` (the gain-shift hoist) use the same
`key.curve.is_back()` predicate the new method walks. Pinned twice —
the `guard-sides` checks in the PD2 door twin and the Rust
`demanded_reference_sides_tracks_the_labels`.

**G2.** Merit op point vs `compute(PDts)`, measured across four
half-space configurations on the rebuilt extension:

| ambient `d` | substrate `d` | merit op point vs `compute(PDts)` |
|---|---|---|
| 0 | 0 | bitwise equal |
| 999 | 0 | bitwise equal |
| 999 | 777 | bitwise equal |
| 0 | 777 | bitwise equal |

Pre-fix the reference would have been displaced by up to **2.47 rad**
at 999 nm and **2.21 rad** at 1776 nm. The slice is also safe at the
boundary: `n_layers = ambient + films + substrate ≥ 2` by construction
(`solver_arrays` chains `once(ambient) + films + once(substrate)`), so
`[1..nl-1]` is at worst empty, never inverted.

**G4.** The doc comment's claim was checked, not taken: the only live
read of `n_front_re`/`n_back_re` in `merit.rs` is `residuals_into` at
`:1032–1035`. `curve_sensitivity_into` genuinely never touches them, so
the absent guard is correct.

**G3, G5, incidental.** Present and correct; all five markdown files
verify cell-consistent, zero raw `||` in any table.

**Gates, all re-run at `b86998f`:** 553 lib / 22 parity / 15 py; 567
all-features; `cargo fmt --all --check`; `cargo clippy --workspace
--all-targets -- -D warnings`; 771 passed + 1 skipped; five
`tools/check_*.py`; ten `validation/review/*.py`; both bit-exactness
fingerprints; seven version sites at `0.6.49`. Every number matches
what `32adee2` claimed.

### 7.2 H1 — the G2 fix was pinned by nothing (**P1**)

Reverting `evaluator.rs:217` to the plain `.iter().sum()`, rebuilding,
and running the full battery left **everything green**: 553 Rust lib
tests, the differential-phase harness (`ALL OK`), 771 pytest. Nothing
in the tree distinguished the two expressions.

The reason is structural: every test stack sets both half-spaces to
zero, which is exactly the regime where the two sums agree. `ar_stack`
builds `LayerSpec::constant("air", 1.0, 0.0, 0.0, nw)`; all four
`DesignStack(...)` constructions in the PD harness pass `0.0`. The
`simulate_fills_pd_metadata_and_complex_t` assertion `total_d == d`
holds under both expressions.

`32adee2`'s message lists the 999/777 measurement under `Twins:`
alongside two real twins, which reads as three. The CHANGELOG was
honest about it (`Measured:`, not a twin), so the durable record never
claimed more than was done — but the effect was the same: a P1 fix in a
repository whose ritual is *one item, one twin* had no twin, and any
later refactor could have undone it silently.

**Applied.** `compute-observable` (the PD4 twin) gains the half-space
case, inside the polarization loop so both `PDts` and `PDtp` are
covered. Verified in both directions: it passes on the fixed tree
(`max|d| = 0.00e+00`, bitwise) and **fails on the reverted tree by
5.03 rad, exit 1**. A twin that was not watched fail is not a twin.

### 7.3 H2 — §0.3 had no row for the 0.6.49 rung (**P2**)

The progress ledger recorded `0.6.49` at `implementation_plan.md:2758`;
the master item table jumped `0.6.48 → 0.6.50`. That is the same
ledger-vs-table split G3 was about, one commit later and in mirror
image.

**Applied.** §0.3 gains `| ~~PDR~~ **DONE (0.6.49)** | ... |`, and the
ledger's ID column now reads `PDR` to match. All 14 tables still verify
cell-consistent.

### 7.4 H3 — the harness docstring numbered 15 strategies for 14 functions (**P3**)

The `guard-sides` entry took the next free number (`15`) but was filed
between `12` and `13`, and had no section of its own: its two checks
live inside item 12's `test_dispersive_rotation_door`.

**Applied.** Renumbered `12a` (same function as 12), which fixes the
ordering and the count together; item 14's entry gains the half-space
twin H1 added.

### 7.5 H4 — §6 CORRECTIONS 1 names a warning that does not fire on this path (**P3**)

CORRECTIONS 1 argues G2's optional half (refusing a non-zero half-space
in `with_films`) is unnecessary partly because "the engine's own warning
at `solver.rs:1131` still fires". It does not fire on the path G2 was
about. `solve_arrays` — the only function carrying that warning — is
reached from `navette-py/src/structure.rs:1644`, the `ScatterMatrix`
door, and from its own tests. The synthesis evaluator's `simulate`
never calls it. Measured: building a `DesignStack` with 999/777 nm
half-spaces and running it through `SmatrixContext.simulate` emits
**zero** warnings.

The conclusion CORRECTIONS 1 reaches is unaffected — with the interior
sum, the two surfaces now agree on any stack, so there is nothing left
to warn about on the merit path, and H1's twin pins that. Only the
supporting clause was wrong. Recorded here rather than edited into §6,
which is an append-only record of what was believed at the time.

### 7.6 Not a regression

§0.3's item count has been stale since the C rows landed: at `a095f60`
it read "Fifteen items" against 18 rows, and PD1 correctly added four
to an already-low number, giving "Nineteen" against 22. Inherited, not
introduced — same class as the CHANGELOG `|n|` cell. Corrected while
editing the same sentence: the count is now derived from the table
(**twenty-three** rows, including `C1`, which carries no version).

### 7.7 Disposition

| # | Finding | Priority | Action |
|---|---|---|---|
| **H1** | The G2 fix was pinned by nothing | **P1** | Half-space twin added to `compute-observable`; watched fail at 5.03 rad |
| **H2** | §0.3 had no row for 0.6.49 | P2 | `PDR` row added; ledger ID matched |
| **H3** | 15 numbered strategies, 14 functions | P3 | `guard-sides` renumbered `12a`; item 14 updated |
| **H4** | CORRECTIONS 1 names a warning off this path | P3 | Recorded here; the conclusion stands, the clause does not |

**No version bump.** Nothing here changes shipped behaviour: one added
regression test, three documentation corrections, and one CHANGELOG
sentence brought in line with the tree it describes. This is the
`138bdcb` / `b86998f` shape — a bookkeeping commit that corrects the
existing section in place rather than opening a new one.

CORRECTIONS:

1. §5's first bullet stands: the `nondispersive-bitwise` literals are
   still not reproduced from outside, and this round did not attempt it
   again. §1.3's independent route continues to carry the claim.
2. This round verified the *fixes*, not the PD series a second time.
   §1's measurements were not re-derived from scratch; the gates were
   re-run and the four G-finding sites re-measured.

---

## 8. Closure — the `nondispersive-bitwise` literals reproduce

Measured after the `v0.7.0` tag, at `eb977bf`. This section closes §5's
first bullet, which §6's CI note and §7's CORRECTIONS 1 both left
standing as the one open gap in this review.

### 8.1 Why it stayed open, and what changed

The obstacle was never the arithmetic — it was the build. Reproducing
check 11 needs a running `_smatrix` from the `0.6.44` tree, and
rebuilding in place would have clobbered the installed extension the
rest of the battery was being measured against. §1.3 therefore proved
the same claim structurally (length-1 row ≡ pre-PD1 scalar; length-1 and
constant-full-length agree bitwise) and the literal route was left
unwalked.

What changed is that the extension under test is now a *released* one.
A second worktree at `v0.6.44` (`a095f60`, detached) with its own
`.venv` builds and installs without touching the live tree at all —
confirmed clean throughout (`## main...origin/main`, no stale `.pyd` in
the worktree, so the PD2 shadowing trap is structurally absent).

### 8.2 The recipe

`maturin develop --release` into the worktree's own venv, then a dump
script that copies check 11's construction **verbatim** — including
`oracle_tf`. That last part is the whole reason the check was
unverifiable from outside: an independently written oracle differs at
ULP level, which moves the embedded *target* by ~1e-15 rad and so moves
the residual. The recipe must reuse the construction, not re-derive it.

The scalar call site is signature-compatible across the hop: at
`0.6.44`, `sim_curves_from_arrays` already carried `total_d`,
`n_front`, `n_back` as `f64` scalars (`#[pyo3(signature = (angles,
wavelengths, total_d=0.0, n_front=1.0, n_back=1.0))]`); PD2 widened
them to arrays without changing what a scalar call means. The same call
runs on both trees.

### 8.3 The comparison

Both recordings were compared separately. They are not the same
measurement twice: the door is fed the oracle's `t`, the engine its own
solver's, so the two residual **vectors differ from each other** — they
agree only on merit. Check 11 pins two independent recordings, and both
had to be re-derived.

| Recording | Pinned at 0.7.0 | Measured at 0.6.44 | |
|---|---|---|---|
| door merit | `7.888609052210118e-29` | `0x1.9000000000000p-94` | `==` |
| door residuals | `[0.0, 0.0, 8.881784197001252e-15]` | `[0x0.0p+0, 0x0.0p+0, 0x1.4000000000000p-47]` | `==` |
| engine merit | `7.888609052210118e-29` | `0x1.9000000000000p-94` | `==` |
| engine residuals | `[0.0, -8.881784197001252e-15, 0.0]` | `[0x0.0p+0, -0x1.4000000000000p-47, 0x0.0p+0]` | `==` |

Hex is shown on the measured side because `repr` round-trips are not
proof of bit-equality; the comparison itself was machine-made (`==` on
the floats, not on their strings) and reported four `MATCH` lines and
`VERDICT: all four literals reproduce from a 0.6.44 build`.
`ctx.evaluate_merit(st) == m2` also holds at `0.6.44`, as check 11
asserts.

### 8.4 What this does and does not change

It changes no code and no conclusion. §1.3's structural route was
already sufficient, and this is a second, independent route to the same
place — the one the test itself takes. What it removes is the residue:
the `[0.7.0]` entry's "A non-dispersive medium is unaffected, bit for
bit" now rests on a proof *and* on a recording that has been re-derived
from the tree it names, rather than on a proof plus a same-hand
recording nobody had reproduced.

It also confirms the provenance line in check 11's own docstring —
*"Literals recorded at 0.6.44 (a095f60..138bdcb tree)"* — which until
now was an unverified claim about where four magic numbers came from.

**No version bump.** Docs only: this section and one CHANGELOG
follow-up sentence. The worktree was removed after measurement.

CORRECTIONS:

1. §5's first bullet is **closed**, and with it the "what remains open"
   clause of §6's CI note and §7's CORRECTIONS 1. Both are left in
   place — they were true when written, and this document is an
   append-only record of what was believed at the time.
2. The closure is retrospective. It was run after `v0.7.0` was tagged
   and published, so it confirms a release rather than gating one. Had
   the literals *not* reproduced, this would have been a finding about
   a shipped release, not a docs commit.
