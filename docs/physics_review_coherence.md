# Physics review — coherent / incoherent propagation

STATUS: review round 1, nothing applied. Scope: where Navette keeps optical
phase, where it destroys it, and whether both are right. Read against
`smatrix/core_engine.rs`, `smatrix/coherent_block.rs`,
`smatrix/optics_core.rs`, `smatrix/needle_operator.rs`,
`smatrix/synthesis/evaluator.rs` and the Python doors in
`src/navette/smatrix/smatrix.py`.

Method: the formalism was read out of the code, then checked three ways that
do not share an error mode.

1. **Against an independent implementation** — Byrnes' `tmm` (the package
   written alongside arXiv:1603.02720), installed in a scratch venv, driven
   over six stacks × 2 polarizations × 3 angles × 3 wavelengths.
2. **Against the definition of incoherence itself** — an incoherent layer is
   one whose round-trip phase is uniform over 2π, so the incoherent answer
   must equal the coherent answer averaged over one full phase period of
   that layer. This uses no second implementation at all.
3. **Against closed form** — the textbook two-surface slab.

All three were run on the installed 0.7.7 wheel.

## 0. Verdict

**The photometric core is correct, and correct for the right reasons.** The
block partition, the attenuation element and the intensity cascade are the
standard Katsidis–Siapkas / Byrnes formalism, implemented in a form that is
numerically better conditioned than the reference. R and T agree with `tmm`
to 1.6e-15 worst case across every mixed stack probed, and the incoherent
answer reproduces the phase-averaged coherent answer to 8e-17 — that is the
definition of incoherence, satisfied to machine precision.

**Two defects sit on top of that core, both in the layers above it.** Neither
is an arithmetic error in the cascade; both are places where a correct engine
is either bypassed or read with the wrong object.

| # | Finding | Sev | § |
|---|---------|-----|---|
| C1 | The synthesis pipeline ignores `coherent: false` entirely — it solves every design as one coherent block | **P1** | §3.1 — refusal **FIXED** `f6f959d` (0.7.9) |
| C2 | Mode A's Stokes vector mixes a front-block cross term with total intensities; Mode A is the default | **P1** | §3.2 — **FIXED** `6380ca0` (0.7.10) |
| C3 | An incoherent flag carries no thickness sanity check — `d = 0` decoheres | P2 | §3.3 — **FIXED** `e982e50` (0.7.12) |
| C4 | Gain media are silently mangled, differently in the coherent and incoherent paths | P2 | §3.4 — **FIXED** `5b59aed` (0.7.11) |
| C5 | Coherence flags on the half-spaces are silently ignored | P3 | §3.5 — **FIXED** `e982e50` (0.7.12) |
| C6 | `ellipsometry()` / `stokes()` / `complex_amplitudes()` carry no front-block caveat, while `dispersion()` does | P3 | §3.6 — **FIXED** `6380ca0` + `e982e50` (0.7.12) |
| C7 | The incoherent cascade has no independent validation — parity is against a port of itself | P3 | §3.7 — **FIXED** `84513dd` (0.7.13) |

## 1. Where phase is kept, and where it is destroyed

The engine walks the stack once (`core_engine.rs::solve_point`), cutting it
into maximal runs of coherent layers:

```text
  block 0        τ        block 1        τ        block 2
[0 .. i₁]  →  layer i₁  →  [i₁ .. i₂]  →  layer i₂  →  [i₂ .. N]
 complex       real          complex       real         complex
```

* **Phase is kept inside a block.** `coherent_block.rs` accumulates the
  complex Redheffer star product over interfaces, applying the propagation
  factor `e^{iβ}` only for layers strictly interior to the block
  (`if idx + 1 < end_idx`). That guard is what makes the partition exact:
  the layer that terminates a block is a half-space or the incoherent layer
  itself, and its propagation belongs to the next stage, not to this one.
* **Phase is destroyed at a block boundary.** The block's four complex
  amplitudes are squared into `(R_f, T_b, T_f, R_b)` and joined with the
  *real* star product `redheffer_product_real_inner`. Squaring before
  joining is the phase scramble; the geometric series
  `1/(1 − R_A,back·R_B,front)` that survives is the incoherent
  multiple-reflection sum.
* **The incoherent layer itself contributes attenuation only** — a two-port
  with `r = 0` both ways and `t = τ = exp(−2·Im β)`. It carries no interfaces
  of its own, because both of them already belong to the blocks on either
  side. This is right, and it is the part implementations most often get
  wrong by double-counting an interface.

The blocks that bound an incoherent layer are the load-bearing detail. Block
`k` runs `[i_k .. i_{k+1}]` and block `k+1` runs `[i_{k+1} .. i_{k+2}]`, so the
flagged layer's front interface is the last interface of one block and its
back interface is the first of the next. No interface is used twice and none
is dropped. `needle_operator.rs::partition_blocks` reproduces the same sweep
for the gradient path, and the comment there ("a flagged medium doubles as
the entrance half-space of the next block") names the invariant correctly.

### Verified against the reference implementation

`tmm`'s `inc_tmm` builds the same object a different way — coherent
sub-stacks, then a 2×2 intensity transfer matrix per incoherent layer. Its
attenuation is

```python
P_list[i] = exp(-4 * pi * d * (n * cos(th)).imag / lam)
```

which is `exp(−2·(2π d/λ)·Im(n cosθ))` — identical to Navette's `τ`.
Worst-case disagreement over the probe set:

| case (mixed coherent/incoherent) | max &#124;ΔR&#124; | max &#124;ΔT&#124; |
|---|---|---|
| bare thick slab | 1.2e-16 | 7.8e-16 |
| AR coating on a thick slab | 5.6e-17 | 1.2e-15 |
| absorbing film + absorbing slab | 5.0e-16 | 0.0 |
| two incoherent slabs, coherent film between | 3.5e-18 | 1.6e-15 |
| two incoherent layers in contact | 1.5e-17 | 4.4e-16 |
| all coherent, via both code paths | 8.3e-17 | 1.2e-15 |

s and p, 0°/30°/60°, 450/550/700 nm, absorbing layers included.

Two places where Navette is **better** than the reference, worth keeping on
the record:

* **The cascade never divides by a transmittance.** `tmm` builds
  `L = [[1/P, 0], [0, P]] · […] / T` and must clamp `P ≥ 1e-30` to keep the
  division alive in an opaque layer. Navette's intensity star product
  divides only by `1 − R·R`, which is bounded away from zero whenever the
  reflectances are, so an opaque incoherent layer needs no clamp and loses no
  digits. The probe with a 500 µm absorbing slab returned `T = 0.0` exactly
  rather than a denormal.
* **The p-polarization energy factor falls out rather than being special-cased.**
  Byrnes needs `Re(n_f cos θ_f)/Re(n_i cos θ_i)` for s but
  `Re(n_f cos*θ_f)/Re(n_i cos*θ_i)` for p — the conjugate is easy to miss in
  an absorbing medium. Navette uses the tilted admittance `η_s = n cosθ`,
  `η_p = n/cosθ` and one formula, `T = |t|²·Re(η_last)/Re(η_first)`. Working
  the algebra through, `|t_η|² = |t_E|²·|cos θ_f|²/|cos θ_i|²` and
  `Re(η_p) = Re(n cos*θ)/|cosθ|²`, so the `|cosθ|²` cancel and the tilted form
  *is* Byrnes' p-pol expression, conjugate included. The oblique absorbing
  p-pol cases above confirm it numerically to 1.6e-15.

### Verified against the definition of incoherence

The stronger check needs no second implementation. Coherence is a property of
the source and the geometry, not of the layer: a layer is incoherent when its
round-trip phase varies by more than 2π across the source bandwidth, the beam
divergence or the layer's own thickness non-uniformity, so that `2β` is
uniformly distributed. The incoherent answer is then *by definition* the
coherent answer averaged over one phase period.

Stack: air / 95 nm H / 50 µm slab / 110 nm L / air, the slab flagged
incoherent, compared against the fully coherent stack averaged over 2048
thicknesses spanning one period:

```text
   0 deg  Rs: incoherent = 0.172728629   averaged = 0.172728629   diff = 8.3e-17
          Tp: incoherent = 0.827271371   averaged = 0.827271371   diff = 1.0e-15
  45 deg  Rs: incoherent = 0.335147051   averaged = 0.335147051   diff = 4.4e-16
          Tp: incoherent = 0.904442825   averaged = 0.904442825   diff = 3.3e-16
```

For an absorbing slab the naive version of this test disagrees at 6e-4, but
that is the test's artifact, not the engine's: sweeping `d` also sweeps `τ`.
Holding `k·d` invariant while the phase turns drops the residual to 4e-7 —
the remainder being the slab's own Fresnel coefficients moving as `k` is
rescaled. The formalism is correct; the test was measuring itself.

Closed form, for completeness — air / n = 1.5 slab / air, where
`R = 2R₁/(1+R₁)` and `T = (1−R₁)/(1+R₁)`:

```text
  FRONT_BLOCK       Rs = 0.076923076923   Ts = 0.923076923077
  COHERENCY_MATRIX  Rs = 0.076923076923   Ts = 0.923076923077
  closed form       Rs = 0.076923076923   Ts = 0.923076923077
```

and `R` is bit-identical across slab thicknesses 10 µm → 3.7 mm, which is the
thickness-independence a lossless incoherent layer must have.

Energy conservation over mixed lossless stacks is exact (`A = ±0.000000000000`
at 0°, 45°, 75°).

## 2. What the three modes actually are

`MODE_A` (`FRONT_BLOCK`, **the default**), `MODE_B` (`COHERENCY_MATRIX`) and
`MODE_C` (`FULLY_COHERENT`).

For **intensities the A/B distinction does not exist** — `solve_point` gates
only `track_cross_channel` on the mode, so `Rs/Rp/Ts/Tp` are the same
computation. `docs/remediation_plan.md` already says this and asserts it.
Mode C ignores the flags and solves one block.

The modes differ only in the **p–s coherency channel** `C = E_p·E_s*`, which
feeds Δ, DOP, S2/S3 and retardance:

* **Mode B** cascades `C` through the incoherent joins with
  `redheffer_product_cross_inner`, whose denominator `1 − C_Ab·C_Bf` is the
  cross-channel geometric series under the rule that *different bounce orders
  are mutually incoherent*. Through the flagged layer itself `C` picks up
  `t_p·t_s* = |e^{iβ}|² = τ`, correctly real, because propagation does not
  distinguish polarization.
* **Mode A** drops the cascade and reports `C_r = r_p⁰·conj(r_s⁰)` — the front
  block alone — while `S0`/`S1` still come from the full intensity cascade.

Mode A is a deliberate bit-for-bit port of the legacy numba kernel, and
`docs/remediation_plan.md` §"CORRECTIONS (0.6.11)" says so. That is a good
reason for the mode to exist. It is not a reason for it to be the default,
and §3.2 is why.

## 3. Findings

### 3.1 C1 — the synthesis pipeline ignores `coherent: false` (P1)

`SmatrixContext::simulate_inner` (`synthesis/evaluator.rs:159`) calls

```rust
solve_coherent_block_fields_dual(start /* = 0 */, end /* = nl-1 */, …)
```

one coherent block over the whole stack, and its own docstring says so:
"Fully-coherent path: block = [0, nl−1) with the substrate as half-space."
`sa.incoherent_flags` is built by `structure.rs:1371` and then never read.
Grepping the whole `synthesis/` module for `solve_point`, `solve_arrays` or
`coherence_mode` returns nothing: the design path never reaches the
coherence-aware engine at all. The needle pass is the same story — it calls
`needle_operator::p_function`, never `p_function_multiblock`, although the
multiblock operator exists, is finite-difference validated, and correctly
refuses to seed a needle inside a flagged spacer.

Meanwhile `coherent` is a **documented, accepted per-film flag** all the way
down: `_FILM_DEFAULTS` in `pipeline.py:86`, `_film_dicts` writes it into the
film dict, `design_config.rs:299` parses it, `LayerSpec.coherent` exposes it.
A user can set it, read it back, and never learn it did nothing.

Measured, on a stack air / 95 nm H / 50 µm glass / 110 nm L / air with the
glass flagged `coherent: False`:

```text
  slab flagged coherent=False          slab flagged coherent=True
  d = 50000.000  merit = 2051.332381309369    merit = 2051.332381309369
  d = 50090.461  merit =  255.295479689127    merit =  255.295479689127
  d = 50137.000  merit = 1368.488113331268    merit = 1368.488113331268
  d = 71234.000  merit =  586.724534898677    merit =  586.724534898677
```

The merit is bit-identical with the flag on and off, and it swings by a
factor of eight as the thickness of a **lossless incoherent** layer moves —
a parameter the physics the user requested says is unobservable. The engine
door, on the same stack, holds `Rs` at `0.172728629140` across that sweep and
moves it by 3.6e-2 only when the layer is genuinely coherent.

So a design carrying a thick substrate layer is optimized against
interference fringes that Navette itself says are not there. The optimizer
will happily tune a 50 µm layer to a quarter wave.

**Disposition.** Two honest options, and they are not equivalent.

1. **Refuse.** `stack_from_layers` / `design_from_program` reject
   `coherent: false` with a message naming the limitation and pointing at
   `ScatterMatrix` for analysis of such stacks. Small, safe, and it converts
   a silent wrong answer into a loud one. This is the minimum.
2. **Implement.** Route `simulate_inner` through the block partition it
   already has in `partition_blocks`, and the needle pass through
   `p_function_multiblock`. The machinery exists and is validated; what is
   missing is the wiring and the thickness-gradient deposits through the
   intensity cascade (the adjoint for that is already derived in
   `needle_operator.rs`, "Incoherent (Mode A/B) block cascade adjoint").

Option 1 should land regardless, and first — option 2 is a phase, not a fix.

### 3.2 C2 — Mode A's Stokes vector is not a Stokes vector (P1)

In Mode A, `S0 = R_p + R_s` counts every incoherent echo, while
`S2 = −2·Re(C_r)` and `S3 = −2·Im(C_r)` are built from the front block alone.
Those are two different physical objects, and the ratio between them is
reported as a degree of polarization.

At **normal incidence on an isotropic stack** the reflected light is fully
polarized — s and p are the same problem. The true answer is `DOP = 1`
exactly. Measured on air / 95 nm TiO₂ / 1 mm incoherent glass / air at 550 nm:

```text
 angle      DOP_R (A)   |rs_c||rp_c|·2/S0      DOP_R (B)   Stokes-average
     0    0.831085970          0.831085970    1.000000000    1.000000000
    45    0.926036997          0.739943962    0.995633865    0.995633865
    70    0.886808788          0.327138226    0.887044148    0.887044148
```

The second column pins the mechanism: at normal incidence `S1 = 0`, so Mode
A's DOP collapses exactly to `2|r_s⁰||r_p⁰|/(R_s+R_p)`, i.e. the front-block
reflectance over the total reflectance — `0.151729360/0.182567587 =
0.831085970`, to the last digit. It is not an approximation error, it is a
ratio of two different stacks.

Δ is worse, because it is not clamped:

```text
 angle    Delta_R (A)    Delta_R (B)      gap
     0    +180.000000    +180.000000     0.000 deg
    45    -166.451530    -166.882551     0.431 deg
    70     -75.727155     -53.503196    22.224 deg
```

**Mode B is the correct one, and that is measured rather than assumed.**
Stokes parameters add incoherently, so the true incoherent limit is the
*Stokes vector* of the fully coherent stack averaged over one phase period of
the flagged layer. That average reproduces Mode B to every printed digit at
every angle (column 4 above, 4096 samples), and Mode A at none of them. The
control matters too: with no incoherent layer anywhere, A and B agree
bit-for-bit and both give `DOP = 1.000000000000`.

Two claims in the tree should be corrected along with the behaviour.

* `validation/parity/smatrix/refs/loom_matrix.py:621` — "S₂, S₃ from first
  coherent block only (incoherent echoes lose phase)". Echoes do not lose
  *their own* p–s phase relation; they lose phase **relative to each other**.
  That is exactly what `1 − C_Ab·C_Bf` encodes, and dropping the echoes from
  S2/S3 while keeping them in S0 is the error the sentence licenses.
* `solver.rs` R6.2, on the DOP clamp — "For reflection the algebra is exact
  (`s1r² + s2r² + s3r² = s0r²` for a single coherent block) and the excess is
  pure round-off, measured at 1.0000000000000004". True for a single block.
  With an incoherent join in Mode A the residual is not round-off and it is
  not an excess: it is a 17 % deficit, and the clamp cannot see it because it
  only clamps from above.

This is textbook-recognised physics, not a Navette peculiarity: backside
reflection from a thick substrate depolarizes, substrates go incoherent above
roughly 50–100 µm in the UV–VIS–NIR, and the standard treatment is exactly a
coherency/Mueller object rather than a field cross term.

**Disposition.** Do not change Mode A's numbers — the parity port depends on
them. Instead:

1. Make `NEEDS_CROSS` observables in Mode A **refuse** (or warn loudly) when
   the stack has at least one interior incoherent flag, naming Mode B. The
   engine already resolves `need_cross` in `resolve_plan`; the flag scan is
   one pass.
2. Alternatively auto-promote to the Mode B cross channel when cross
   observables are requested, keeping Mode A as an explicit opt-in for
   legacy reproduction. This changes default output, so it is a bump.
3. Either way, change the Python default to `COHERENCY_MATRIX`, or document
   in `CoherenceMode` that `FRONT_BLOCK`'s ellipsometric outputs are a legacy
   front-surface convention and not the stack's.

**As applied** — `6380ca0`, 0.7.10. Option 1, as a refusal rather than a
warning, at the two doors (`ScatterMatrix.compute` and
`solver::solve_arrays`) with the engine untouched; option 2 was rejected as
the house anti-pattern (silently answering a different question than the one
asked), and option 3's default flip was rejected because `FRONT_BLOCK`
matches the legacy port by a remediation-plan decision and flipping it would
re-litigate that for stored results — with the refusal in place the default
is no longer dangerous, only loud. Option 3's second half shipped: the
`CoherenceMode` docstring now says it. What made refusal cheap is that the
parity port enters through the raw `core_engine` pyfunction, below both
doors — this section's "the engine already resolves `need_cross` in
`resolve_plan`" would have put the check inside the engine and killed it.
See r2 §5 item 3.

### 3.3 C3 — an incoherent flag carries no thickness sanity check (P2)

Coherence is destroyed by path-length spread, so the flag is only meaningful
when the layer is thick against the coherence length `L_C ≈ λ²/Δλ`. Navette
checks nothing:

```text
  d =   0.0 nm   incoherent Rs = 0.040040419   coherent Rs = 0.042579995
  d =   1.0 nm   incoherent Rs = 0.040040419   coherent Rs = 0.042578510
  d =  10.0 nm   incoherent Rs = 0.040040419   coherent Rs = 0.042432867
```

A layer of **zero thickness** decoheres the stack and shifts `Rs` by 2.5e-3.
Frustrated total internal reflection disappears entirely — a 200 nm air gap
flagged incoherent returns `R = 1.000000, T = 0.000000` beyond the critical
angle at every gap thickness, where the coherent answer at 45° is
`R = 0.716, T = 0.284`. That is the right answer for a *thick* gap and a
wrong one for the gap that was actually described.

`tmm` has the same hole, so this is a shared limitation of the formalism's
usual implementations rather than a Navette-specific error. But Navette is a
design tool: if C1 is ever fixed by implementation rather than refusal, an
optimizer will be free to move a flagged layer's thickness, and nothing
currently stops it from walking one to zero while the model keeps answering
in the incoherent limit.

**Disposition.** A warning at stack construction when a flagged layer's
optical thickness is below a few wavelengths, plus a sentence in
`incoherent_flags`' docstring giving `L_C = λ²/Δλ` as the criterion and
50–100 µm as the practical substrate threshold. If C1 is implemented, a
flagged layer should additionally be refused as an optimization parameter
unless it is absorbing.

**As applied** — `e982e50`, 0.7.12. Both halves, at both doors. "A few
wavelengths" became five, but the message does not lean on that number: it
quotes the source bandwidth the layer would need,
`Δλ > λ²/(2nd)`, so the threshold is derived and the caller can check it
against their own source. The optimizer clause stays parked with C1's
implementation phase, which has not happened.

### 3.4 C4 — gain media are silently mangled, and differently per path (P2)

`coherent_block.rs` forces decay by conjugating the propagation phase
(`if beta.im < 0.0 { beta = conj }`); the incoherent path clamps instead
(`if beta_imag < 0.0 { 0.0 }`). Neither refuses. Measured on a 2 µm layer at
`k = −0.01`:

```text
  coherent    k = -0.010   Rs = 0.014003623  Ts = 0.608577324  A = +0.377419053
  coherent    k = +0.010   Rs = 0.014029268  Ts = 0.609691819  A = +0.376278913
  incoherent  k = -0.010   Rs = 0.076953120  Ts = 0.923089546  A = -0.000042666
  incoherent  k = +0.010   Rs = 0.054811350  Ts = 0.583945129  A = +0.361243521
```

A gain layer reports **positive absorption** in the coherent path — the
propagation was flipped to loss while the interfaces still saw gain, so the
result is not any physical system, merely close to the absorbing one. In the
incoherent path the same layer becomes transparent with a small negative
absorption, i.e. the energy books do not close.

`structure.rs:149` already refuses `k < 0` ("Nominal expansion produced
k < 0"), so one door is guarded and the other is not — the same "both doors
must explain the same rule" shape as `AMBIENT_DROP_EXPLANATION`, which is the
model to copy.

**Disposition.** Refuse `Im(n) < 0` at the `ScatterMatrix` / `solve_arrays`
door with a shared explanation constant, exactly as the absorbing-ambient
rule is handled.

**As applied** — `5b59aed`, 0.7.11. Refused, with the shared constant
(`GAIN_MEDIUM_EXPLANATION`) this section asked for, but one level deeper
than "at the door": `Solver::assemble`, which both named doors and the raw
FFI pass through, so the rule is written once rather than three times. The
Structure door named here as the guarded one now carries the same text —
its "check provider data" said nothing about what the other paths would
have done with the number. The free `needle_gradient` bypasses the
constructor and carries its own copy, over the needle material as well as
the host stack. See r2 §5 item 4 for why this could go deeper than C2's.

### 3.5 C5 — half-space coherence flags are silently ignored (P3)

`solve_point` scans `ni < idx_n` and applies `τ` only when
`next_incoh < idx_n`, so `incoherent_flags[0]` and `incoherent_flags[last]`
are never consulted. Flags on rows 0 and last change nothing:

```text
  flags=[0,0,0]  Rs = 0.038375631357     flags=[0,0,1]  Rs = 0.038375631357
  flags=[1,0,0]  Rs = 0.038375631357     flags=[1,0,1]  Rs = 0.038375631357
```

This is **correct** — a semi-infinite half-space has no second surface and no
finite path to decohere — and it matches `tmm`, which requires its first and
last entries to be `'i'`. It is also silent, and "flag the substrate
incoherent" is the first thing a user reaching for this feature will try.
`solve_arrays` already warns when row 0 / last carries a non-zero thickness;
the coherence flag deserves the same sentence, saying that a thick substrate
is modelled as an interior layer plus an exit half-space.

**As applied** — `e982e50`, 0.7.12. The sentence sits beside that one, at
both doors, and says exactly that. The constructor docstring also lost the
"thick substrate" example, which was inviting the mistake in the first
place. The measurement above is now a test rather than a table
(`test_a_half_space_flag_really_is_bit_identical`), with an interior
control so it cannot pass vacuously.

### 3.6 C6 — the front-block caveat is on the wrong docstrings (P3)

`dispersion()` and `differential_phase()` both say "Physically meaningful only
for coherent stacks (`FULLY_COHERENT` mode or a stack with no incoherent
boundaries)". Good. But `rs_c`/`ts_c`/`phi_*` have the same property — the
Rust `OpticalState` doc says "Complex amplitudes are the first coherent block
(Modes A/B) or the whole stack (Mode C)" — and the Python
`complex_amplitudes()` says only "Complex r/t coefficients". `ellipsometry()`
and `stokes()` say nothing at all, and per §3.2 they are the two that need it
most.

Demonstrated: with an incoherent slab present, `|rs_c|² = 0.151729` against
`Rs = 0.182568`. Both numbers are in the same output dict, and nothing says
they describe different stacks.

### 3.7 C7 — the incoherent cascade has no independent validation (P3)

The parity suite compares the engine against
`validation/parity/smatrix/refs/loom_matrix.py`, which is a port of the same
algorithm with the same block sweep (its own comments cite Katsidis &
Siapkas). That pins the port, not the physics. The Rust unit tests exercise
the incoherent path for bitwise agreement between the full and
intensity-only walks (`intensity_path_matches_full_path_bitwise`), which is
again self-consistency. Searching the tree, there is no closed-form or
external check of the incoherent result anywhere.

The three probes in this review are cheap and are the missing twins:

* the two-surface slab against `2R₁/(1+R₁)` and `(1−R₁)/(1+R₁)`;
* thickness-independence of a lossless flagged layer, asserted on `to_bits()`;
* the phase-average identity of §1 — the only one that tests the *definition*
  rather than an agreement, and the one that would have caught C2 had it been
  written for the Stokes vector as well as for R/T.

**As applied** — `84513dd` (0.7.13). All three landed as
`validation/review/incoherent_check.py` (the eleventh review harness) with
the Stokes-vector variant as a fourth part, plus
`validation/smoke/test_incoherent_physics.py` (17 tests) as the regression
twins. Measured: the closed form to 1e-12 in both modes; bit-identical `Rs`
from 10 µm to 3.7 mm; the phase average to 8.3e-17 (`Rs`, 0°), 1.0e-15
(`Tp`, 0°), 4.4e-16 / 3.3e-16 at 45°, and 7.4e-07 / 9.0e-08 on an absorbing
slab with `k*d` held invariant. Each part carries a control that fails it.

The Stokes part did what round 1 predicted it would. `COHERENCY_MATRIX`
reproduces the averaged vector to ≤6.7e-16 on all four components — mode B
verified against the definition rather than against a port — while
`FRONT_BLOCK`, reached past its 0.7.10 refusal through the raw engine,
misses it by 2.26e-02 (`S2_R`) and 3.91e-03 (`S3_R`) *while its `S0_R` and
`S1_R` pass the same average to 1e-16*. That split is C2 measured from the
definition; every earlier piece of C2 evidence was one implementation
disagreeing with another.

Two traps in the identity itself are now pinned as tests, because both were
live mistakes while writing it and both look like engine bugs. Sweeping `d`
to turn the phase also sweeps `tau = exp(-2*Im(beta))`, so the naive
absorbing average is taken over a stack whose absorption moves underneath it
(3.4e-04, measuring itself); the test asserts both that holding `k*d` fixed
works and that the naive version fails. And `S0..S3` are bilinear in the
fields while `DOP` is a nonlinear function of them: the first draft averaged
`DOP` and reported a 5.8e-03 "disagreement" that was entirely its own, since
this stack is non-depolarizing and every coherent sample has `DOP = 1`
exactly, while the averaged Stokes vector has `DOP = 0.994186`. That gap is
the physics — partial depolarization is what incoherent superposition
produces — and mode B reproduces it to 1e-16.

The twins average over 64 phase samples rather than the harness's 2048: the
coherent answer is analytic and periodic in the round-trip phase, so an
equispaced Riemann sum converges geometrically (8 samples → 1.7e-07, 16 →
5e-15). The whole file runs in 0.18 s, and the convergence *rate* is itself
asserted, so a future discontinuity in the sweep cannot quietly turn these
identities into approximations.

## 4. Sources

* S. J. Byrnes, *Multilayer optical calculations*, arXiv:1603.02720, and the
  `tmm` package written alongside it — used as the reference implementation.
* C. C. Katsidis and D. I. Siapkas, *General transfer-matrix method for
  optical multilayer systems with coherent, partially coherent, and
  incoherent interference*, Appl. Opt. **41**, 3978 (2002) — the formalism the
  block cascade implements; already cited by `loom_matrix.py`.
* M. C. Troparevsky et al., *Transfer-matrix formalism for the calculation of
  optical response in multilayer systems: from coherent to incoherent
  interference*, Opt. Express **18**, 24715 (2010) — the coherent→incoherent
  transition by random-phase averaging, and partial coherence as the
  intermediate case Navette does not model.
* Coherence length `L_C = λ²/Δλ`; substrates depolarize above ~50–100 µm in
  the UV–VIS–NIR (Mueller-matrix ellipsometry practice) — the basis for §3.3
  and the physical justification for §3.2.

## 5. Not findings

Recorded so a later pass does not re-derive them.

* **Partial coherence is out of scope and that is fine.** Navette's flag is
  binary; Katsidis–Siapkas and Troparevsky both describe a partially coherent
  middle ground. Nothing claims to implement it.
* **Roughness is not decoherence.** The `w_function` / Névot–Croce factors
  damp specular amplitudes inside a coherent block; they are a separate axis
  from the coherence flag and are documented with their own validity budget.
* **The needle already refuses a flagged host.** `locate_hosts_multiblock`
  errors with "host layer {j} is incoherent-flagged", which is the right
  refusal — a coherent needle inside a layer declared incoherent is not a
  thing. That this refusal is currently unreachable from the synthesis path
  is C1, not a separate finding.
* **`field_profile` takes no coherence flags** because an eigenmode field
  profile is a coherent-stack concept. Worth a docstring line, not a finding.
