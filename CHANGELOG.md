# Changelog

All notable changes to Navette are recorded here. Work items reference
`docs/remediation_plan.md` (Rx.y) and `docs/code_review.md` (§).

## [0.6.18] — Hygiene: the archive, the headers, and two docs that lied (R6.4, part 1)

No behaviour change; the solve fingerprint is unmoved at `30d96909…3c6c`. Four
of R6.4's seven items — the ones that are about the repository rather than
about the physics.

### Removed

- **`attic/`** — 26 retired files, import-isolated and never referenced by live
  code, but responsible for roughly 300 of the code review's 1 528 dead-code
  findings and most of its scaffolding-comment hits. Deleted; git history keeps
  them. `docs/plans/` and every `validation/**/refs/` stay — those are the
  parity oracles, not archive.

### Added

- **SPDX headers** on 213 sources that lacked them (`// SPDX-License-Identifier:
  LGPL-3.0-or-later`, `#` for Python), placed after a shebang and PEP 263 coding
  line where present so nothing moves ahead of a module docstring. 8 files
  already carried one. `COPYING` and `COPYING.LESSER` were already in the tree.

### Documented

- **`validation/README.md` rewritten.** Its tree predated `review/`,
  `parity/synthesis`, `parity/color`, `regression/**` and nine of the twelve
  smoke suites, and still pointed at `crates/navette-materials`, a path that has
  not existed since the Rust consolidation. Now carries the real tree, the
  per-directory test counts (287 smoke / 327 regression / 56 parity / 22
  goldens), the distinction between what `pytest validation` collects and the
  `review/` harnesses it does not, and the numba-thread-pool timing
  contamination that the two whole-engine parity scripts call `cooldown()` for.
  The "the rewrite is behind the kernel it replaced" note is now history with a
  date on it: that was true at 0.64–0.77× and R5.3 fixed it.
- **`docs/plans/bug_fix_plan.md` closed.** BUG-A, BUG-B and BUG-C had unticked
  checkboxes and, directly underneath, their own notes recording the fix — the
  tracker said open while the tree said fixed, which is precisely the audit-trail
  rot the review filed. Ticked against named passing tests
  (`test_asymmetric_interface_position_and_pair`,
  `test_inverted_roughness_follows_plane`,
  `test_partial_inversion_run_edge_clean`, all in
  `validation/regression/structure/test_inversion.py`), with a header saying the
  page is a record rather than a work item and that `expander.py` has been
  `rust/navette/src/structure/expansion.rs` since the port.

## [0.6.17] — The small physics nits, and which of them were real (R6.2)

Four items the review filed as minor. Two were: a missing clamp and a stale
docstring. One was the opposite of what it looked like — the `+ 0.0` the review
called "a no-op … harmless" is load-bearing, and removing it would have swung
Delta_R by 2π at normal incidence. The fourth needed a decision and got one.

### Fixed

- **`DOP_R` is clamped to ≤ 1**, as `DOP_T` already was. The asymmetry was the
  review's finding; the excess on a physical stack is pure round-off
  (`1.0000000000000004`, measured), but a degree of polarization above 1 is not
  a number a caller can do anything with. **This changes output bits** — see
  Changed.
- **`needle_slopes` docstring**: τ̂ is `iβ′(1+r₁₂²)/(1−r₁₂²)`, not `2iβ′/(1−r₁₂²)`.
  The code, the module header and the star-product test always had the right
  form; only the one-line summary was stale.

### Added

- **Sellmeier domain guard.** `n² = 1 + Σ Bᵢλ²/(λ²−Cᵢ)` is only a function
  between its poles: on a pole it is infinite, and inside the first resonance
  n² < 0 and `sqrt` returns NaN — which then travels through every layer of the
  stack and arrives as a NaN spectrum with nothing left pointing at the
  material that produced it. `sellmeier_domain_check` now refuses that grid with
  a message naming the wavelength and where the fit's resonances actually are.
  Both user-facing entry points run it (the PyO3 binding and the spec
  evaluator); the kernels themselves are untouched, so the check costs one scan
  for a non-finite value and the arithmetic that works out *which* pole only
  runs once something is already wrong.
  - The domain edge is not the pole. For BK7 the first resonance is at 77.5 nm,
    but n² stays positive down to ~70.7 nm (the other terms hold it up), so
    50 nm evaluates to a finite, physically meaningless n = 0.47. The guard
    catches non-finite values; it does not pretend to know where a fit stops
    being *trustworthy*. That is what the coefficient set's paper is for.

### Changed

- **Solve fingerprint moved**, for the first time in this series:
  `a99e8383…f154` → `30d96909…3c6c`. The delta is the `DOP_R` clamp and nothing
  else — verified by recomputing the unclamped ratio from the `S0_R`…`S3_R`
  channels (which stay unclamped) and asserting the emitted channel equals
  `min(raw, 1)` at every one of the fingerprint's 170 points. 20 of those
  exceeded 1: 18 by one ulp, and 2 by a lot (6.1 and 51.1) — those two come
  from the fingerprint harness randomizing the *incident* medium to an
  absorbing index, where `R = |r|²` is not an energy ratio at all and `Rs`
  reaches 86 and `Rp` goes negative. With a real ambient the same stack's
  worst `DOP_R` is 0.99997. Noted as an input-validation gap, not fixed here.

### Documented

- The `+ 0.0` on `s2r`/`s3r`/`s2t`/`s3t` is a negative-zero flush, not dead
  weight: `-2.0 * 0.0` is `-0.0`, and `atan2(-0.0, negative)` answers −π where
  `atan2(+0.0, negative)` answers +π. An isotropic stack at normal incidence
  has `cross_r` exactly real, so this is the ordinary case, not a corner. The
  numba reference carries the same flush at the same four places, with a
  comment; the Rust side now has one too, plus a test pinning Delta_R = +π.

## [0.6.16] — The two argmin backends, and what they are worth (R4.4c-argmin)

The last row left open under R4.4: `OptimizerBackend::ArgminGaussNewton` and
`ArgminTrustRegion`, behind a new `opt-argmin` cargo feature. Both are
reference points rather than candidates, and the interesting part of this
release is the measurement, not the plumbing.

### Added

- **`opt-argmin`** (off by default): `LmConfig(optimizer="argmin_gauss_newton")`
  and `"argmin_trust_region"`. Both unbounded, so both run on the same
  `IntervalMap` reparametrization `minpack_lm` uses — the map and the chain
  rule on the analytic Jacobian are shared, not re-derived. The feature reuses
  `nalgebra` through `argmin-math`'s `nalgebra_v0_34` backend, so enabling it
  alongside `opt-minpack-lm` pulls one linear-algebra crate, not two.
- **`OptimizerBackend::runs_to_max_iterations()`** — true for the argmin trust
  region alone. argmin's `TrustRegion::terminate` returns `NotTerminated`
  unconditionally: it has no convergence test, so `max_iterations` is its only
  exit and `MaxIterations` is its normal outcome, not a failure. A caller
  reading the termination reason needs that declared rather than inferred.
- **Part E of `validation/review/lm_check.py`**, skipped unless the wheel was
  built with the feature: the two backends on the same three refold starts the
  rest of the harness uses, against scipy.
- Ten cargo tests across both backends (interior optimum, box containment,
  cost convention, error propagation, analytic-Jacobian precedence, and the
  two failure modes below), plus name round-tripping and feature-hint coverage
  now driven from `ALL_BACKENDS` rather than a hand-kept list.

### Known — measured, on `validation/review/lm_check.py`'s refold starts

| backend | two films | two films, far | three films |
|---|---|---|---|
| `trf` | 8 evals | 21 evals | 32 evals |
| `argmin_trust_region` | 218 evals | 224 evals | 393 evals |
| `argmin_gauss_newton` | refuses | refuses | 201 evals, cost 2690.8 vs 2670.0 |

- **`argmin_trust_region` is correct and expensive.** It reaches the same
  optimum as scipy on all three starts (agreement 1e-13 to 1e-9) and pays
  7–27× `trf`'s evaluation count to do it, because it cannot stop early.
- **`argmin_gauss_newton` is not usable on this problem class.** The undamped
  step is `(JᵀJ)⁻¹Jᵀr`, so a singular `JᵀJ` makes it undefined — and a film
  driven toward zero thickness, or a bound the interior map has flattened, is
  enough to produce one. It refuses two of the three starts; on the third it
  runs to the iteration cap without converging. The refusal now carries a
  message that says which of those happened and names a damped backend to use
  instead, rather than surfacing argmin's bare "Non-invertible matrix".
- Neither backend changes anything on a default build: `available_optimizers()`
  still returns `["builtin", "trf"]`, and naming one the wheel lacks is still a
  `ValueError` with the rebuild command.

## [0.6.15] — The derive pass is no longer the serial half (R5.2)

R5.3's scaling bench made the gate for this item concrete: after the
construction fix, a rigorous twelve-channel request at 20 000 points spent
~1.1 ms of its 2.4 ms in one serial loop — the per-point derive, where every
requested channel's algebra and transcendentals (`sqrt`, `atan`, `atan2`,
`arg`) are evaluated. The solve itself was already on the rayon pool; this pass
was not.

### Changed

- **The derive pass splits across the rayon pool.** The point range and *every
  live destination buffer* are halved together, down to a 512-point leaf, so
  each output index is still written exactly once, by one task, with the
  expression it had when the pass was serial. There is no reduction and no
  shared accumulator: bit-identity is structural, not a tolerance.
- The 45 channel buffers and their halving are generated from one list, so a
  new channel cannot be added to the pass and forgotten by the split.

Measured with `bench_core_engine_scaling.py`, 6 layers, 1 angle, median of 25
(0.6.14 → 0.6.15):

| points | photometry | | rigorous | |
|---:|---:|---|---:|---|
| 500 | 0.183 → 0.177 ms | 1.03× | 0.231 → 0.231 ms | — |
| 2 000 | 0.364 → 0.344 ms | 1.06× | 0.440 → 0.401 ms | 1.10× |
| 5 000 | 0.592 → 0.540 ms | 1.10× | 0.769 → 0.666 ms | 1.15× |
| 20 000 | 1.539 → 1.499 ms | 1.03× | 2.387 → 1.684 ms | 1.42× |
| 60 000 | 4.229 → 3.700 ms | 1.14× | 6.977 → 4.022 ms | 1.73× |

The gain tracks how many channels the request asks for, which is what the pass
costs: a four-channel photometry request has almost nothing to divide, a
twelve-channel ellipsometric one has 1.7× at 60 000 points. Against the numba
kernel the rigorous request now runs at ~1.1–1.8× at 20 000–60 000 points,
where 0.6.13 had it at ~0.4–0.5×.

### Added

- Four solver tests: the split pass is bit-identical to a single serial call on
  a grid that is not a multiple of the leaf; the split reaches every index
  (buffers start as NaN, so a range nobody visits is caught directly rather
  than by comparing two equally-gapped runs); a grid below the leaf still
  derives; and absent channels — 41 of the 45 on a photometry request — survive
  the halving at every level.

### Known

- **The leaf size is a scheduling knob, not a numerical one.** Swept at 256 /
  512 / 2048 / 8192: the grid sizes that matter land within ~3% of each other,
  the one real effect being that a leaf above the grid means no split at all
  (2 000 points cost 0.428 ms at 2048 and 0.386 ms at 512). 512 is the smallest
  leaf that still buys something measurable.
- Small grids remain dispatch-bound and behind the numba kernel; see 0.6.14.
- Bit-identity held: the differential fingerprint is unchanged, all ten review
  harnesses exit 0, both parity oracles pass.

## [0.6.14] — `core_engine` was paying for its own construction (R5.3)

R2.4a timed the Rust engine at 0.5–0.65× the numba kernel it replaced and the
plan put the blame on the emit path — a `Vec<OpticalState>` walked once per
requested channel on the way out. The profile says otherwise. At 20 000 points
the emit is 0.01–0.25 ms of a 4.3–4.8 ms call; **`Solver::new` is 2.2–2.5 ms of
it**, and 1.77–2.24 ms of *that* is building the index cache as a
`Vec<Vec<Complex64>>` — roughly 40 000 heap allocations per call, on a path
whose callers construct a fresh solver for every evaluation.

### Changed

- **The index cache is one flat wav-major `Vec<Complex64>`** (stride
  `n_layers`) instead of a `Vec` of per-wavelength `Vec`s, and so is its
  reciprocal. Same numbers, same order, two allocations instead of 40 000.
- **`Solver::from_wav_major_flat`** — a second constructor that takes the
  interleaved `[re, im]` buffer the Python boundary already holds. `core_engine`
  was transposing that into layer-major so `Solver::new` could transpose it
  straight back; both passes are gone. `new` keeps its layer-major signature and
  its behaviour.
- **`flat_cache` is built on first use** (`OnceLock`), not in the constructor.
  Only `needle_gradient` reads it; every ordinary solve was paying for it.
- **The emit takes its `Solution` by value.** `PyArray::from_vec` moves, so
  borrowing meant cloning every channel — one extra full copy of the result.
- **Both parity oracles and the new scaling bench pause between engines.** See
  Known: the 0.5–0.65× figure was partly a measurement artefact.

### Added

- `validation/benches/smatrix/bench_core_engine_scaling.py` — the ratio across
  the grid (500 → 60 000 points) at two request breadths, which is what the
  single-size speed line in the parity scripts could not show.
- Four solver tests (the two constructors agree bit-for-bit on a dispersive
  stack; a mis-sized flat cache is refused; the flat constructor validates what
  `new` validates; the serial and parallel reciprocal paths agree across the
  threshold) and two Python tests that a wrong-length index cache raises instead
  of panicking through the binding.

Measured with `bench_core_engine_scaling.py`, 6 layers, 1 angle, median of 25:

| points | photometry before → after | | rigorous before → after | |
|---:|---:|---|---:|---|
| 500 | 0.254 → 0.183 ms | 1.39× | 0.324 → 0.230 ms | 1.41× |
| 2 000 | 0.598 → 0.364 ms | 1.64× | 0.718 → 0.440 ms | 1.63× |
| 5 000 | 1.233 → 0.592 ms | 2.08× | 1.474 → 0.769 ms | 1.92× |
| 20 000 | 4.149 → 1.539 ms | 2.70× | 5.023 → 2.387 ms | 2.10× |
| 60 000 | 12.63 → 4.229 ms | 2.99× | 15.74 → 6.977 ms | 2.26× |

Against the numba kernel the same bench now reports roughly 0.3–0.5× at 500
points, 0.6–0.8× at 5 000, and 0.8–1.2× at 20 000–60 000 — the rewrite passes
the kernel it replaced on large grids and is still behind on small ones, where
~150 µs of rayon dispatch is most of the call.

### Known

- **numba's default threading layer contaminated the original ratio.** Its
  workers spin-wait after a call returns, so timing the two engines back to back
  charges whichever runs second. The same Rust call at 20 000 points measures
  1.52 ms on its own, 2.00 ms immediately after a block of numba calls, and
  1.58 ms after a one-second pause. Forcing numba onto the non-spinning
  `workqueue` layer also removes the effect, but makes numba itself ~70% slower
  here, so both benches take the pause instead and leave each engine in its
  native configuration.
- **Small grids are still dispatch-bound.** Below ~5 000 points the parallel map
  costs more to start than the work it distributes, but serializing them is
  worse: a `SOLVE_PAR_THRESHOLD` experiment made 2 000 points 5× slower
  (0.35 → 1.75 ms), because serial is ~875 ns/point against ~80 ns/point
  effective under rayon. Not taken.
- **What is left is the derive loop.** On a rigorous 20 000-point request it is
  ~1.1 ms of a 2.4 ms call, and it is serial. That is plan item R5.2, whose gate
  ("only if the scaling bench shows the serial fraction matters") this bench now
  satisfies.
- Bit-identity held throughout: the differential fingerprint is unchanged, and
  all ten review harnesses and both parity oracles pass.

## [0.6.13] — What the batch unweave was actually spending (R5.1)

`unweave_collection` ran at 0.22–0.48× the Python reference and the plan blamed
per-key hash-map rebuilds, `String` clones and a whole-collection clone. None of
those were there any more. What it spends is a `memcpy`, and the reason the
reference does not is that the reference **does not own its data**.

### Changed

- **The batch unweave copies only the span the plan reaches**, not the whole
  input curve, and shares that one copy across a key's fragments. Frames rarely
  tile the entire curve handed to them; the part no fragment reads was being
  copied per key.
- **One collection-wide lock per key instead of one per fragment.** The key→
  frames map was taking its write lock once per (key, frame) pair — 16 384 times
  for a 512-key call, all of it on one lock. `OpticalCollection::
  map_frames_to_key` takes it once per key; the singular form now delegates to
  it.
- **`SpectralDataFrame::set_data` hashes the key once**, not twice
  (`contains_key` then `insert` → `insert` and read what it returns).
- **Above 8192 fragments the per-key loop runs on the rayon pool.** Gated on the
  fragment count rather than the byte volume, because the bytes are the part
  that does *not* parallelise — see Known.

Measured on the bench grid, Rust-side `unweave_collection`, median of 7:

| grid × frames × keys | before | after | |
|---|---:|---:|---|
| 1 000 × 32 × 128 | 0.296 ms | 0.234 ms | 1.27× |
| 10 000 × 64 × 256 | 2.585 ms | 1.346 ms | 1.92× |
| 10 000 × 256 × 64 | 1.024 ms | 0.446 ms | 2.29× |
| 100 000 × 512 × 32 | 2.957 ms | 2.057 ms | 1.44× |
| 100 000 × 32 × 512 | 41.29 ms | 35.92 ms | 1.15× |
| 100 000 × 16 × 32 | 1.848 ms | 1.8 ms | — |
| 500 000 × 32 × 16 | 3.711 ms | 3.7 ms | — |

### Fixed

- **A curve that is not one value per wavelength no longer panics in a batch.**
  The distribution plan addresses the source by position, so a short curve was
  indexed straight out of its buffer — a Rust panic surfacing as
  `pyo3_runtime.PanicException`, which does not even derive from `Exception`, so
  an `except Exception:` around the call did not catch it. The Python wrapper
  checked this for a single `unweave` and never for `unweave_batch`. It is a
  `ValueError` naming the key now, on both paths.
- **A rejected batch writes nothing.** Everything that depends only on the plan
  and the grid is checked before the first key is distributed, so a grid that
  misses part of a frame no longer rejects the call after an arbitrary prefix of
  it has already been applied — which, once the loop can run on a pool, would
  have meant a *scheduling-dependent* prefix.

### Added

- **`batch_equals_per_key`** in the spectral bench's correctness section, on both
  engines: a 40-frame × 256-key batch must write exactly what 256 single
  `unweave` calls write. The batch path resolves one plan, shares one
  materialised source per key and may run pooled; nothing in the per-op timings
  would have noticed the two drifting.
- 13 tests for `opticalweaver`, which had none of its own (it was covered only
  from Python), and four for the wrapper's batch guards.

### Known

- **`unweave_collection` is still ~0.28× the Python reference at 512-key scale,
  and that gap is not an optimization target.** The Python reference stores
  *views into the caller's numpy buffer*; the Rust engine owns its fragments.
  Mutating the caller's array after the call changes the stored curve on the
  Python engine and does not on the Rust one. Owning it means copying it: 512
  keys × 100 000 points is 400 MB in and 400 MB out, and 36 ms of that is
  22 GB/s of DRAM traffic — this machine's ceiling, not a code path. Measured
  1 / 4 / 32 rayon threads: 39.7 / 35.1 / 36.5 ms. The plan's stated goal for
  this item, "≥ 1× vs numba at 512-key scale", is reachable only by adopting the
  aliasing, which would make a fragment change under a caller who edited their
  own array afterwards. Recorded in `docs/remediation_plan.md`; not done.

## [0.6.12] — The method the docs have been naming (R4.6)

`thick_opt.rs` has said since the rewrite that it replaces
`scipy.optimize.least_squares(method="trf")`. What it implements is a
Levenberg-Marquardt that keeps `lb ≤ x ≤ ub` by vetoing and clamping the step
it already solved — bounds applied after the fact, to a step computed as
though they were not there. `synthesis::trf` is the method that sentence was
naming: trust-region reflective (Branch-Coleman-Li), where the bounds enter
the subproblem.

Hand-rolled, no new dependency, no cargo feature — so `optimizer="trf"` works
on a standard wheel. The built-in stays the default.

### Added

- **`LmConfig(optimizer="trf")`** — trust-region reflective. The trust region
  is reshaped every iteration by the Coleman-Li scaling `D = diag(√v)`, `v`
  being the distance to the bound the anti-gradient points at, and each step
  is the best of three candidates: the trust-region step cut back to the first
  bound it hits, that step **reflected** off the bound, and the constrained
  Cauchy step.
- **`synthesis::trf::trust_region_reflective`** — same signature and same
  `LmResult` as the built-in LM, so the two are interchangeable behind
  `run_optimizer`. The subproblem is Moré's, solved from the augmented QR
  R4.4b already built rather than from an SVD; the Newton recurrence on the
  secular equation is identical term for term.
- **`lm_check.py` part D** — the scipy comparison this repository could not
  make before. Parts A–C compare two *different* algorithms and can only ask
  for the same optimum; D is the same algorithm on both sides, so it asks for
  the same answer: costs to 1e-9 relative, thicknesses to 1e-4 nm, including
  on a merit whose optimum sits **on** the clamp. This also closes R4.4d's
  deferred "parametrize over every enabled backend" row.
- **`OptimizerBackend::feature()` is public** — the Python binding was
  carrying its own copy of the backend→cargo-feature mapping for its rebuild
  hint, which is the kind of pair that drifts.

### Known

- **TRF never returns a thickness exactly on a bound.** Its iterates must stay
  strictly interior — that is what keeps the Coleman-Li scaling
  differentiable — so where the built-in returns `50.0` it returns
  `49.999999999986834`, at an identical merit. The removal sweep compares
  against `clamp_min`, so a film driven to the bound is still removed, but a
  caller testing `x == ub` will be disappointed. This is the reason TRF is not
  made the default here.
- **`lambda_init`, `lambda_up`, `lambda_down`, `damping` and
  `gtol_scale_invariant` do nothing on this backend.** They are
  Levenberg-Marquardt settings; the trust-region radius plays their role and
  is not user-settable, exactly as in scipy.
- **The plan's algorithm sketch for this item described a different method**
  (a frozen MINPACK column-norm `D`, and "2-D subspace minimization" — which
  is scipy's sparse `tr_solver`, not its reflections). Both corrected in
  `docs/remediation_plan.md`; following either would have cost the scipy
  oracle that makes the item checkable.

## [0.6.11] — The two whole-engine parity oracles are back (R2.4a)

`parity/smatrix/test_core_engine_photometry_only.py` and
`test_core_engine_rigorous_ellipsometry.py` were the only tests that compared
the *entire* engine — inputs in, thirteen channels out — against the numba
kernel it replaced. Since R2.4 (0.5.8) they had been skipping, because they
loaded a `navette_matrix` extension that does not exist in this repository.
They now run against `navette._smatrix.core_engine`, and the parity layer has
no NEEDS PORT rows left: 13/13 channels, agreement ~1e-14 to ~1e-16.

### Changed

- **Both scripts ported to the request-driven `core_engine`.** The legacy
  kernels took `calc_s, calc_p` and a `debug_flag` and always returned a
  fixed-width tuple; the successor takes a request mask and emits only what
  was asked for. Every element of both tuples is still compared.
- **The coherence mode is `FRONT_BLOCK`, and the mapping is asserted.** The
  legacy ellipsometry kernel took reflection S₂/S₃ from the first coherent
  block and transmission S₂/S₃ from a Mueller cross-term product across
  blocks — exactly `core_engine.rs`'s `track_cross_channel == false` path.
  Half the cases flag an incoherent layer so the modes actually differ, and
  the script requires `COHERENCY_MATRIX` to *disagree* with the reference:
  without that, "mode A is the legacy treatment" would be an untested claim
  that happens to hold because nothing separates the modes.
- **A skipped channel is asserted absent, not compared against zeros.** The
  legacy kernel zero-filled the polarization it was told to skip. Comparing
  the missing channel against those zeros was tried deliberately: it passes,
  and would keep passing if the mask were ignored entirely. The port asserts
  the key is missing, and separately that `Rs` from an s-only request is
  bit-identical to `Rs` from the full one.
- **The 13th value moved rather than disappeared.** `conservation_err` is no
  longer an engine channel; it is `solver_energy_conservation` (R1.2), and
  the port reconstructs and compares it from the four intensities.

### Known

- **The rewrite is slower than the kernel it replaced** on this workload —
  0.5–0.65× across 500 to 60 000 grid points — and roughly half the cost
  looks to be outside the physics. Recorded as **R5.3 (P2)** rather than
  fixed here: a performance change does not belong in a parity port.
- **`Delta` is compared on the circle as a guard, not a finding.** The last
  ellipsometry case is built to straddle ±π, and the script asserts it really
  does — but both engines pick the same side everywhere, so a plain
  difference would pass this file today. Removing the circle handling was run
  as a deliberate break and was *not* caught; it is recorded in the plan that
  way rather than dressed up.

## [0.6.10] — Which solver runs is now a choice, and a named one (R4.4c, work item B)

The thickness optimizer's residual system was always solver-agnostic — a
closure over `x` and a `JacobianSource` over `x` — but only one solver was
ever reachable, and the plan's claim of scipy parity had no second
implementation to be parity *with* (review §3.6, §18.2). `synthesis::optimizer`
is now the single place that decides who gets handed that pair, and the
argmin ecosystem's reference LM is a `--features` flag away.

**Nothing changes for a standard build.** The default backend is the built-in
LM, the optional crates are default-off, and the whole cargo and Python suites
pass unchanged.

### Added

- **`LmConfig(optimizer=...)`** — `"builtin"` (default) or `"minpack_lm"`. The
  latter is the `levenberg-marquardt` crate (rust-cv, MINPACK `lmdif`-derived,
  MIT), available when the wheel was built `--features opt-minpack-lm`.
- **`navette._smatrix.available_optimizers()`** — what *this* wheel can run.
  Backends are cargo features, so the answer is a property of the build, and
  a caller should be able to ask rather than discover it by failing.
- **`synthesis::optimizer`** — `OptimizerBackend`, `OptimizerResult`,
  `run_optimizer`, and `IntervalMap`. A new solver is an arm here, not a
  second copy of the call site; R4.6's TRF backend takes the same seam.
- **`optimize_thicknesses_report`'s dict gained `"backend"`** — which solver
  produced the answer, beside which Jacobian path it took.

### Changed

- **`LmConfig` gained `backend`.** The plan's separate `OptimizerConfig`
  "mapping 1:1 onto `LmConfig`" would have been `LmConfig` plus one field, and
  two structs that must agree field for field are a synchronization bug
  waiting to be written. Every other knob already carries over — ftol, xtol,
  gtol and max_evals are MINPACK's own names.

### Notes

- **Bounds are not a shared contract, and this is the reason the built-in
  stays the default.** Navette's LM vetoes and clamps the step it just solved,
  so an optimum may sit *exactly* on a bound — a film driven to zero thickness
  is a real answer the synthesis loop then removes. Every reference
  implementation in the ecosystem is unbounded, so those backends run on
  `x = mid + half·tanh(u)` and their optima are *strictly inside* the box,
  with the gradient vanishing as a bound is approached. A contract
  difference, not a bug.
- **A missing backend is refused, never substituted.** Naming one the wheel
  was not built with raises `ValueError` carrying the rebuild command. A
  silent fall back to the built-in would make every "compared against the
  reference LM" claim a comparison with ourselves.
- **The plan's ndarray/nalgebra friction does not arise.** argmin-math 0.5
  ships a `nalgebra_0_34` backend and `levenberg-marquardt` 0.15 is built on
  nalgebra 0.34, so the two share one linear-algebra crate and navette's
  ndarray 0.17 pin is untouched — no adapter, no per-iteration conversion, no
  version bump. (The plan's alternative, argmin-math's default `Vec<f64>`
  backend, would not have worked: it has no `ArgminInv`, which GaussNewton
  requires.)
- **`stepbound` is not exposed.** It is a MINPACK-only knob; putting it on the
  shared config would add a setting that does nothing on the default backend.
  The crate's default stands until a real workload argues otherwise.
- Eight deliberate breaks, all caught: the interval map made linear, the chain
  rule dropped, the chain rule applied along rows instead of columns, the
  result left in `u`-space, MINPACK's ½‖r‖² not converted to Navette's ‖r‖²,
  an unavailable backend falling back to the built-in, differenced Jacobians
  counted as analytic, and the `atanh` clamp removed.

## [0.6.9] — The thickness optimizer stops guessing its own Jacobian (R4.5, increment ii)

`build_jacobian` was central differences: **2n full residual evaluations per
LM iteration**, each one a complete solver sweep over the angle × wavelength
grid, and each column carrying the difference noise that matters most exactly
where it hurts — small steps near an optimum (review §3.6). The Jacobian is
now assembled analytically:

```text
    J[i,k] = Σ_terms  ∂r_i/∂(curve value) · ∂(curve value)/∂d_k
```

with the left factor from 0.6.8's `MeritSpec::curve_sensitivity` and the right
from the same solver sweep that produces the curves. Analytic is the default;
`jacobian="fd"` restores the differences.

**No pinned optimum moved** — the whole cargo and Python suites pass unchanged,
so the §8 golden protocol had nothing to triage — and the S-matrix engine's
bit-exactness fingerprint is unchanged.

### Changed

- **`LmConfig.jacobian`** (`"analytic"` default, `"fd"`) chooses the path. The
  analytic one applies when every residual row is covered by the sensitivity
  chain; a spec with a phase target or a color demand declines *per run* and
  is differenced, so the setting is a preference, not a promise.
- **The thickness derivative is the needle operator with the needle material
  set to the host's own index.** Growing film *j* by δ is inserting a slab of
  *n_j* inside layer *j*: r₁₂ vanishes, ρ̂ with it, and τ̂ = iβ_j is the bare
  propagation slope. The existing dual-number composition `U ⊗ N ⊗ L` then
  differentiates the whole block through the same Redheffer star product the
  forward solver uses — every multiple-reflection path included, no
  hand-expanded algebra, and no way for the derivative to drift from the value
  it differentiates. Intensities follow the solver's own `finalize`: R = |r_f|²
  and T = |t_back|²·f_back, with the boundary-admittance factor constant under
  a film thickness and riding outside the derivative.

### Added

- **`SmatrixContext::simulate_with_deposits`** — the simulate that also reads
  off ∂(Rs, Rp, Ts, Tp)/∂(thickness) per grid point. The curves come back
  bit-identical to `simulate`'s, which is a test, not a hope.
- **`synthesis::jacobian`** — `CurveDeposits` and `assemble_jacobian`, the seam
  where the two halves meet. J is materialized m×n with per-destination writes
  and a fixed accumulation order per row: no cross-point reduction, nothing
  whose order depends on scheduling (§13).
- **`SmatrixContext.optimize_thicknesses_report(stack)`** → `(merit, report)`,
  with `iterations`, `evals`, `cost`, `termination`, `gain_ratio` and
  **`analytic_jacobians`** — how many of the run's Jacobians came from the
  analytic chain. Which path ran is now something a caller can *ask*, rather
  than infer from a timing; the bench and the fallback tests assert on it.
- **`levenberg_marquardt_with`** and the `JacobianSource` trait. A source may
  supply a Jacobian, decline (`Ok(None)` → difference this iteration), or fail
  (`Err` → abort). Declining and failing are deliberately different: a source
  that believed it had a Jacobian and was wrong should not hide behind a
  slower run.
- **Cargo tests** — the per-element analytic-vs-FD cross-check the plan names
  as required (deposits against `simulate` differenced in thickness space, and
  the assembled J against `residuals` differenced the same way: worst relative
  deviation ~1e-9), the coverage test at the seam the optimizer uses, the
  absorption two-channel row, and five driver-level tests of the hook.
- **`validation/smoke/test_analytic_jacobian.py`** — the end-to-end half, CI's
  subset: the two modes reach the same optimum, the report says which ran, a
  phase demand falls back and says so.
- **`bench_refold.py`** grew a Jacobian section. At ten films, iteration for
  iteration: **1.7× wall-clock, 16× fewer residual evaluations**.

### Notes

- **The speedup is not 2n, and the plan's estimate needs correcting.** Building
  the deposits costs one `StackFields` decomposition per point *and*
  polarization, so an analytic iteration is roughly three solver sweeps rather
  than one — against 2n + 1 for the differenced one. That is 1.7× at n = 10,
  not 20×. What does scale as promised is the residual evaluation count (16×
  fewer at n = 10, and it keeps growing with the stack) and, more importantly,
  the noise: the columns are exact.
- **`validation/review/lm_check.py`'s A2 stationarity check was testing the
  iteration cap.** At the default 200 iterations, both solvers are still
  crawling on the far-start cases — a run that stopped on `MaxIterations` has
  not claimed a stationary point, so asserting one of it says nothing about
  the optimizer. The harness now gives those cases room to converge and checks
  the termination criterion explicitly. (The analytic path is what exposed
  this: it takes a different path through the same flat valley and was a few
  parts in 10⁵ behind at iteration 200. With room, both converge to the same
  cost and both are fixed points.)
- Eight deliberate breaks, all caught: the T deposit losing its flux factor,
  differentiating the forward amplitude instead of the backward one, the film
  index not offset past the ambient, the factor of two in d|z|²/dd, the needle
  taking the ambient's index instead of the host's, a transposed J, a row
  keeping only its last term, and the coverage gate removed.

## [0.6.8] — Every residual row can say what it depends on (R4.5, increment i)

The merit half of the analytic Jacobian. `MeritSpec` can now report, for each
residual row it produces, which simulated curve values that row reads and with
what derivative — the factor `J[i,k] = Σ_terms ∂r_i/∂(curve value) ·
∂(curve value)/∂θ_k` needs on the left. The right-hand factor, the per-point
deposits from the solver sweep, is increment ii; nothing in the optimizer uses
this yet, so no synthesis result moves.

Splitting R4.5 this way is deliberate: the merit half is verifiable on its own
against a finite difference of `residuals()` **in curve space**, with no solver
in the loop, so when the two halves are joined a disagreement has only one
place left to be.

### Added

- **`MeritSpec::curve_sensitivity(&sim) -> MeritSensitivity`** — `rows[i]` is
  the list of `CurveTerm { curve, angle_row, wavelength, d_residual }` for
  residual `i`, in `residuals()` order, index for index. Covered: pointwise
  and integral-mean targets, over intensity and absorption curves, for every
  constraint kind and every real transform, on aligned and interpolated target
  grids.
- **`MeritSensitivity::uncovered` / `is_complete()`** — rows whose dependence
  this pass does not carry: phase targets (`arg()` of a complex row, plus a
  differential reference that depends on the stack's total thickness directly)
  and color demands (a spectrum-wide integral through the CIE chain, whose
  derivative `build_needle_targets` already emits in a different shape). They
  occupy their rows and are *listed*, never returned as empty-and-covered — a
  caller that finds one among its active demands must fall back to a
  finite-difference Jacobian, and `is_complete()` is how it asks.
- **Ten cargo tests**, led by `sensitivity_is_a_finite_difference_of_residuals`
  — five kinds × two transforms, central-differenced through `residuals()`
  itself. That is the anti-drift guard: this pass walks the target grids a
  second time rather than threading a sink through `residuals_into`, which is
  the bit-exactness-critical path, and the cross-check is what keeps the two
  walks in step.

### Fixed

- **`n_residuals()` over-counted integral targets.** It reported one component
  per target point; `residuals()` pushes exactly **one** row per integral
  frame, because the constraint kind applies to the mean, not to each point.
  A spec with one nine-point integral frame advertised 9 components and
  produced 1. The count is now what the vector actually holds, and its
  docstring says the remaining caveat out loud: a frame whose grid misses the
  simulated grid is skipped by `residuals_into` and cannot be counted from the
  spec alone.

### Notes

- The kinks are the constraint, not an approximation: `a`/`b` rows are exactly
  flat on their satisfied side and `r` rows exactly flat inside the band, so
  an inactive constraint comes back with *no terms at all* rather than terms
  that happen to be zero. `Log` is likewise flat below its 1e-12 clamp — the
  residual genuinely stops responding there, and so does a finite difference.
- Bit-exactness fingerprint unchanged; `bench_refold` unchanged (one LM
  optimize 2.1 ms — this increment is not yet on any hot path).

## [0.6.7] — The LM stops squaring its own condition number (R4.4b)

The bounded Levenberg-Marquardt behind every thickness optimization solved the
**normal equations** — `(JᵀJ + λ·diag(JᵀJ))δ = −Jᵀr` — which squares the
condition number of `J`. Thin-film stacks with correlated layers are exactly
where `JᵀJ` goes singular, and the λ-floor was what kept bailing the solve out
(review §3.6). The damping semantics are unchanged; how the step is obtained
is not.

Work item A of R4.4. Backend selection and the optional argmin-ecosystem
solvers (R4.4c) are a separate item; the built-in stays the only backend.

### Changed

- **The damped step is solved by QR.** The step is now the least-squares
  solution of the augmented system `[J; √λ·D] δ ≈ [−r; 0]`, whose condition
  number is the square root of the normal-equation matrix's. It is the *same*
  step — `RᵀR = JᵀJ` exactly — obtained from the square roots. Following
  MINPACK's `qrsolv` structure, `J` is factored once per iteration (m·n²) and
  each λ trial then costs a 2n×n QR (n³), so the per-iteration work does not
  grow; building `JᵀJ` is gone entirely, since only its diagonal was ever
  needed. Unpivoted Householder: the damped system is full rank for every
  λ > 0 with a floored `D`, so pivoting would add rank diagnostics, not
  solvability. The normal-equation solve survives as the fallback when the
  factorization degenerates, and as the oracle the QR is cross-checked
  against.
- **Damping is Nielsen's gain ratio**, not a fixed ×5 / ÷3 ladder: ρ = actual
  over predicted reduction; on acceptance λ ← λ·max(⅓, 1−(2ρ−1)³) and ν ← 2,
  on rejection λ ← λ·ν and ν ← 2ν. The old ladder is kept as
  `LmDamping::Fixed` (`damping="fixed"` from Python) so the two can be
  compared on the same problem.
- **The predicted reduction is computed for the step actually taken.** The
  bound veto and clamp rewrite δ *after* it is solved. Scoring it with the
  full LM step's prediction overstates what the model promised whenever a
  bound is active, so ρ comes out too small — over-damping, and premature
  ftol exits right at the boundary. This is the one place a naive MINPACK port
  goes wrong, and it is pinned by a test that fails with ρ = 0.049 instead of
  1.0 on an exactly-linear clamped problem.
- **ftol needs both reductions** (MINPACK semantics): actual *and* predicted
  relative reduction below `ftol`. A step clipped by a bound can deliver very
  little while the model still sees plenty of room — small progress is not the
  same fact as no progress left.
- **gtol gained the scale-invariant form**, `cos∠(J·e_j, r) ≤ gtol`, alongside
  the existing ‖Jᵀr‖∞ test (`gtol_scale_invariant`, default on;
  `gtol_scale_invariant=False` from Python restores the old behaviour alone).
  ‖Jᵀr‖∞ answers a different question after a parameter is rescaled — nm
  versus µm — and the cosine answers the same one. They are different notions
  of stationarity and can exit at different points on flat valleys.

### Added

- **`LmResult.gain_ratio`** — ρ of the last accepted step, against the
  prediction for the step actually taken. NaN when nothing was accepted. A
  diagnostic, not a control: it tells a synthesis stalling at its bounds apart
  from one that is finished.
- **`validation/review/lm_check.py`** — the scipy parity §18.2 says the plan
  docs have been claiming and nobody had reproduced. Runs
  `SmatrixContext.optimize_thicknesses` and
  `scipy.optimize.least_squares(method="trf")` over the *same* residual system
  and bounds, plus the cargo tests' pinned optima recomputed by scipy.
  `validation/smoke/test_lm_parity.py` is the subset CI runs.
- **Fourteen cargo tests** in `thick_opt.rs`, including the two the plan names
  as the guards that make the port trustworthy: the QR-vs-normal-equations
  step cross-check (rel ≤ 1e-8 where both are valid, and a Vandermonde case
  where the QR is measurably better), and the clipped-step prediction.

### Notes

- Two corrections to R4.4d, both found by running it:
  - **Cost parity from a far start is not a well-posed criterion.**
    Reflectance against thickness is oscillatory, so two local solvers started
    far from an optimum legitimately land in different basins. The harness
    compares them where the basin is unambiguous, and asserts the basin-free
    properties (improvement, stationarity) on the far starts.
  - **The QR's advantage is conditional, and the harness says so.** At the λ
    the solver starts from, Marquardt damping regularizes `JᵀJ` enough that
    the two formulations are indistinguishable. The QR matters where λ has
    decayed towards nothing near a good optimum — precisely where the last
    digits are decided.
- No behaviour changed in the S-matrix engine: the bit-exactness fingerprint
  is unchanged. `bench_refold` reports one LM optimize at 2.0 ms against the
  2.1 ms baseline.

## [0.6.6] — The clippy gate is blocking (R2.3a)

CI ran `cargo clippy` from 0.6.0 onward, but with `continue-on-error: true`
against 267 findings. An advisory gate over a tree that cannot go green is a
gate nobody reads. The tree is now clean and the step fails the build.

### Changed

- **`cargo clippy --workspace --all-targets -- -D warnings` is a blocking CI
  step.** 267 findings → 0. Nothing in the diff changes behaviour: the
  bit-exactness fingerprint (12 seeded random stacks × every `Request` bit,
  needle gradients, an eigenmode landscape, hashed over the raw f64 bytes)
  is byte-identical before and after, and was re-checked after each pass.
- **108 `excessive_precision` findings came from one generated file.**
  `rust/navette/src/color/matrices.rs` is emitted by
  `validation/benches/color/gen/gen_matrices.py`, which formatted with
  `f"{x:.17e}"` — enough digits to name every double, and ~5 more than needed
  for most of them. The generator now uses `repr(float(x))`, which is the
  shortest string that round-trips. All 108 literals were verified bit-identical
  through `struct.pack("<d", ...)` against the previous tree. Fixing only the
  output would have reintroduced the warnings on the next regeneration.
- **Remaining lint decisions are crate-level `allow`s with written rationale**,
  in both `lib.rs` files, rather than scattered per-site suppressions:
  - `neg_cmp_op_on_partial_ord` — all 15 sites are `!(x > 0.0)`-shaped
    validators. The negation is what rejects NaN; clippy's suggested
    `x <= 0.0` would wave NaN through every one of them. Load-bearing.
  - `too_many_arguments` — the wide kernel signatures are what R6.1 / R6.3
    exist to restructure. Silenced so the gate can go blocking now.
  - `needless_range_loop` — flat angle-major index arithmetic
    (`k = a * num_wavs + w`), where the variable indexes several arrays at
    different strides.
  - `type_complexity`, `should_implement_trait` — plan/cache tuples, and
    inherent `from_str` parsers returning the crate's own `String` error.

### Notes

- `cargo fmt --check` stays advisory (R2.3b): making it blocking requires a
  tree-wide reformat that rewrites `git blame` for the whole crate, and the
  plan asks for a style decision before that lands.

## [0.6.5] — The one shipped example did not run (R4.3)

`examples/spectralweave_example.py` failed on a fresh clone, at line 22 of 24.
Running it was how three separate defects came to light.

### Fixed

- **The example is rewritten against the documented API** and now runs end to
  end. It used the raw native `OpticalWeaver` with a key tuple; the documented
  path is `SimulationWeaver` + `OpticalFragment`. It also handed `unweave` a
  `linspace` over the same *range* as the stored frame, which is not the same
  thing as the same *grid* — it shared exactly one wavelength out of a hundred.
  The example now weaves three fragments on three grids, takes the continuous
  curve, unweaves an edited target back onto it, asserts the round trip, and
  shows the rejection path deliberately.
- **`unweave`'s error explains the contract instead of counting.** The failure
  surfaced as `Length mismatch` raised four levels down in
  `SpectralDataFrame::set_data`, which knows nothing about grids or targets.
  `unweave` and `unweave_collection` now check coverage where the context
  exists and say how many of the frame's wavelengths the supplied grid carried,
  which frame (by span), that matching is exact to 1e-12 rather than
  interpolated, and where to get a grid that works. `set_data`'s own message
  now names the frame and both counts.
- **`OpticalFragment` is hashable, so `unweave_batch` can be called at all.**
  Its signature is `dict[OpticalFragment, np.ndarray]`, but `frozen=True` on a
  dataclass holding two numpy arrays synthesises a `__hash__` over those
  arrays, which raises `unhashable type: numpy.ndarray` — and a `__eq__` that
  raises `truth value of an array is ambiguous`. `eq=False` gives identity
  semantics for both, which is the honest answer for a container of mutable
  buffers. Found by writing the example's last line.

### Added

- **`SimulationWeaver.unweave` rejects a non-`OpticalFragment` template** with
  a `TypeError` naming what it got and how to build one. The backend's key
  tuple is an implementation detail; passing one used to surface as
  `AttributeError: 'tuple' object has no attribute '_rust_key'`. A curve whose
  value count does not match its wavelength count is now caught here too,
  rather than becoming a distribution-plan error about something else.
- **`validation/smoke/test_examples.py`** — every `examples/*.py` runs in a
  subprocess (repo root as cwd, `PYTHONIOENCODING=utf-8`, `MPLBACKEND=Agg`,
  120 s timeout) and must exit 0 *and* print something. A subprocess keeps
  pytest's process state and any module-level side effects out of the suite.
  An empty `examples/` fails the guard rather than passing vacuously.
- **`validation/smoke/test_spectralweave_wrapper.py`** — 10 tests over the
  round trip, the coverage message's contents, the guards, and
  `OpticalFragment` as a dict key. Includes the positive control that matters:
  a grid covering *one* of two frames exactly is still accepted, so the new
  check cannot start over-rejecting partial updates.

### Evidence

Teeth proven by dropping two probe files into `examples/`: one exiting 3
(caught by the returncode assert) and one that runs silently (caught by the
output assert). 656 passed + 2 skipped; cargo 376; zero `cargo check`
warnings; all nine review harnesses exit 0.

### Known gaps

- `examples/` holds exactly one file. The runner is written to cover whatever
  lands there, but there is no second example to prove the parametrization
  against — the probe files above stood in for one.

## [0.6.4] — Color demands were optimized against zero in Python (R4.2)

The documented Python needle flow is `build_needle_targets` →
`needle_gradient`. Every pointwise demand came through it correctly. A
**color** demand — Lab/DE2000, White, Yellow, dominant wavelength — came
through as a gradient of exactly zero, with no warning and no error.

A color demand integrates the whole spectrum into one residual, so it has no
per-point target to fold into the `r`/`t` (target, weight) pairs. The native
fold emits its analytic `dF/dcurve` per solver point instead
(`NeedleTargets.grad_r`/`grad_t`), and `needle_pass.rs` adds those into the
same accumulator as the pointwise terms. But the **Python** dict from
`build_needle_targets` dropped both arrays, and `needle_gradient` had no
argument that could have accepted them — so the folded `r` pair was all
zeros, the engine faithfully returned zero, and a caller assembling their own
needle cycle in Python optimized a color demand against nothing.

Not a crash and not a wrong shape: a silent wrong-optimizer. `run_design` was
never affected — it folds and deposits entirely inside Rust.

### Added

- **`build_needle_targets` returns `"grads_r"` / `"grads_t"`** — flat `na*nw`
  f64 arrays, the same angle-major layout as the target/weight pairs. They are
  *not* a (target, weight) pair: each entry is the chain-rule factor
  `g = dF/dcurve` for one solver point, with the demand's weight, its current
  residual and the U-curve half already folded in. All-zero for a spec with no
  color demands, so passing them unconditionally is safe.
- **`needle_gradient` accepts `grads_r` / `grads_t`** (Python keyword, PyO3
  method, and the free `needle_engine` entry), deposited into `P` and `P_T`
  respectively via the existing `p_coherent_grad_r_from_fields` /
  `p_coherent_grad_t_from_fields` kernels — the same kernels, with the same
  zero-skip, that `needle_pass.rs:889-896` calls. The native internal path is
  untouched; its cargo tests still pin the deposit semantics.
- **A non-zero bucket handed to a channel that is not being computed is an
  error**, naming the missing bit. Dropping it silently would be this same bug
  wearing a new hat. An all-zero array is not a demand and passes, because the
  fold hands both arrays through unfiltered. Non-finite entries are rejected by
  index, matching the R3.1/R3.2 convention.
- `validation/review/color_grad_python.py` (new, 24 checks, exits 1 on drift)
  and `validation/smoke/test_color_needle_python.py` (new, 13 tests — the part
  CI runs).

### Changed

- `navette.synthesis`'s module docstring — the documented flow — now shows the
  color branch with a worked two-channel call, and says outright that omitting
  the two kwargs loses the color contribution.
- The Part C comment in `validation/review/color_merit_check.py` no longer
  describes the Python dict as omitting the buckets; it now says why that
  harness deliberately keeps hand-assembling the chain rule (it is the oracle
  the new harness is checked against, so it must not consume the engine's
  own answer).

### Evidence

`grads_r` vs a central-difference `dF/dR` on the sim row: max relative
deviation **7.1e-10** across all 31 points. The deposit against the pointwise
kernel driven to the same scalar (`g·P_ref/(2R)`, the two differ only in one
factor): **2.3e-15 … 3.1e-15** — 1e-12 is a real bound here, not a rounded
1e-6. End to end, `dF/d(thickness)` through the Python path against a
thickness FD of the merit, for all three layers: **1.6e-8 … 1.7e-7**; against
`color_merit_check.py`'s hand-assembled chain rule: **8.3e-16 … 9.4e-16**.
Superposition `P(pointwise + color) = P(pointwise) + P(color)` is bit-exact.

Teeth proven by disabling the R deposit and rebuilding: 2 of the 13 smoke
tests and 10 of the 24 harness checks fail; restored and re-verified.

### Known gaps

- Only the front R and T channels carry color buckets — that is the native
  fold's own v1 scope (`CurveId::Ru`/`Tu` take the ÷2 U-curve half; the back
  siblings and absorption fold nothing). The Python path now mirrors exactly
  what the native path computes, no more.
- The deposit rides `P`/`P_T`, so a color demand requires those bits. That is
  a real constraint, now stated in an error rather than discovered by a
  gradient that quietly reads zero.

## [0.6.3] — Every dependency floor was fiction (R4.1)

`requires-python = ">=3.12"`, and **not one** of the declared floors has a
cp312 wheel:

| declared | first cp312 wheel | now |
|---|---|---|
| `numpy>=1.22.0` | 1.26.0 | `numpy>=2.0` |
| `scipy>=1.8.0` | 1.11.2 | `scipy>=1.13.0` |
| `pyyaml>=6.0` | 6.0.1 | `pyyaml>=6.0.1` |
| `numba>=0.56.0` (extra) | 0.60.0 | `numba>=0.61.0` |

A floor is a promise that the package works with at least that version. These
could not be installed at all on the Python the project declares.

### Fixed

- **`numpy>=2.0`** — the floor R4.1 named. Two independent floors meet in this
  package and they are not the same thing: the extension is `abi3-py312` (one
  wheel for 3.12 and every later 3.x), while the `numpy` Rust crate (0.28)
  compiles against the numpy **2** C-API and makes no abi3-style promise of
  its own. A wheel built that way cannot load against numpy 1.x whatever the
  Python ABI says. Both floors are now documented next to each other in
  `pyproject.toml`.
- **`scipy>=1.13.0`** — the first scipy that both supports numpy 2 and ships
  cp312 wheels. Beyond R4.1's stated scope; found while making the gate
  blocking, and the same one-line defect.
- **`pyyaml>=6.0.1`**, **`numba>=0.61.0`** (extra) — same.

### Changed

- **`ci.yml`'s floor job is blocking, and no longer checks numpy by name.**
  It reads *every* `>=` floor out of `pyproject.toml`, pins them all at once
  with `--only-binary=:all:`, runs `pip check`, and then runs the full suite
  against that combination on the **lowest** supported Python. Checking numpy
  alone would have left three broken floors behind a green check — which is
  exactly how they survived. Renamed `numpy-floor` -> `dependency-floors`.
- The job now runs the whole suite rather than importing. A floor that imports
  but breaks a kernel is still a lie, and pinning the numba extra means the
  parity layer compares against its reference there instead of skipping.

### Verified

Locally, before pushing, on a clean Python 3.12.13 venv: the four floors
install together as wheels, `pip check` is clean, the PEP 517 build produces a
**release** extension, and `pytest validation` gives **631 passed, 2 skipped**
at `numpy 2.0.0 / scipy 1.13.0 / pyyaml 6.0.1 / numba 0.61.0` — identical to
the development environment.

## [0.6.2] — Eigenmode search is boxed, and says when it found nothing (R3.3)

`char_func` is `|1/r(n_eff)|^2`, so a pole drives it to zero — but so does
letting `|n_eff|` run away, and the search was unbounded. Seeded on a stack
with no s-polarized mode, it reached `n_eff = -2.0e8 + 8.4e6j` and reported a
characteristic value of `1.3e-217`: twelve orders *better* than the real
surface-plasmon pole it never found, for a trial index two hundred million
times the largest index in the stack. Nothing said anything.

### Fixed

- **The minimizer is confined to a physical box**, `3 x max |n|` over the
  stack. Every guided mode satisfies `min|n| <= |n_eff| <= max|n|`, so 3x is
  generous even for leaky and substrate-side modes. Candidates are *projected*
  into the box rather than rejected, so Nelder-Mead stays well-defined and
  slides along the boundary instead of reflecting off an invisible wall.
- **`refine_mode` refuses to call something a mode when the characteristic
  value says otherwise**: new `max_residual` kwarg, default `1e-6`, raising
  `ValueError` with the settled `n_eff`, the value, and what to do instead.
  Pass `max_residual=None` for the old unconditional behaviour. The margin is
  not tight — on the pinned SPP a real pole reaches `1e-17`, eleven orders
  below the threshold, while the boxed non-mode settles around `1e-3`.
- `char_func` now reports a non-finite `r` as `1e30` rather than propagating
  NaN into `find_minima`, whose comparisons are undefined on NaN.

### Changed

- **The box is in `char_func_xy`, not `char_func`.** The plan put it in the
  shared `char_func`, which the landscape scanner also calls — that would
  paint `1e30` across any part of a user-requested scan range lying outside
  the box, corrupting a diagnostic the caller explicitly asked for at the
  range they asked for it. The minimizer picks its own points and needs walls;
  the scanner is already bounded by its caller. A test pins the distinction.
- `find_eigenmodes` is untouched and does **not** go through the residual
  contract: its seeds come from a bounded landscape scan, which is what makes
  them seeds. On this stack it correctly returns `[]` for s-polarization — the
  documented workflow was never the broken one.

### Added

- `validation/smoke/test_eigenmodes.py` — 19 tests. Nothing pinned the
  eigenmode path before, which is why the runaway survived a full review
  cycle. Covers the SPP from three seeds, the field profile peaking at the
  metal/glass interface (the physics check a characteristic value cannot
  give you), the runaway from five seeds including two already outside the
  box, the flat s-polarized landscape, the threshold behaviour, and the
  landscape *not* being walled off. Removing the box fails exactly seven of
  them, including `DID NOT RAISE` on the original runaway seed.

### Known gaps

- The plan pinned the SPP at `n_eff = 1.0459458 + 0.0015949j, val < 1e-10`.
  That is not reproducible from the stack as described (the wavelength is not
  stated), and the nearby feature this code converges to — `1.0471183 + 0j`,
  val `1.8e-6` — is a shallow resonance, not a pole: `Im(n_eff) -> 0` on a
  lossy metal, and it sits below the substrate index. The genuine pole on this
  geometry is the glass-side plasmon at `1.7139417 + 0.0226069j`, val `1e-17`,
  and that is what is pinned. See the plan's R3.3 notes.

## [0.6.1] — Needle depth is range-checked in release builds (R3.2)

`needle_slopes4_ddz` guarded its host-layer invariant with `debug_assert!`,
which is compiled out of every release build — the one everybody runs. And
`locate_depth_in` does not reject an out-of-range depth either: it falls
through to the last layer of the range and returns a depth past that layer's
thickness. So `z = 1000` on a 400 nm stack, or `z = -10`, returned a gradient
for a needle that is not where the caller put it, with no error.

### Fixed

- **`z` outside the eligible span now raises `ValueError`**, naming the index,
  the value and the span. The two paths have different spans *within the same
  call* — the coherent kernels are confined to `[start_idx, end_idx]`, the
  multiblock cascade walks every non-ambient layer — so each is checked only
  when the request actually reaches it. A multiblock request is not held to
  the block's narrower bound; a call using both must satisfy both.
- **A non-finite `needle_n_per_wav` entry now raises**, naming the index. It
  previously made the whole gradient NaN with nothing to point at. (Beyond
  R3.2's stated scope, but the same check in the same place.)

Both endpoints stay legal and a `1e-9` tolerance absorbs round-off: a caller
writing `np.linspace(0, sum(d), k)` lands on the bottom endpoint with whatever
error the sum accumulated, and rejecting that would make the obvious way to
write the call fail intermittently.

### Changed

- The check lives in `solver::needle_gradient`, not in the two PyO3 entry
  points. The plan put it binding-side; one implementation in the core covers
  both bindings *and* Rust callers, and it returns `Result<_, String>` which
  the bindings already map to `PyValueError` — so the user-visible error is
  identical either way.
- The kernel keeps its `debug_assert!`s as Rust-caller invariants, with
  messages that now say which contract was skipped. A release `assert!` in an
  O(1) hot kernel would trade a silent wrong answer for a panic, which is
  worse on a library path.
- `validation/review/garbage_in.py`: the three needle rows flip from
  `silent-clean` / `NaN-in-output` to `raises`. The fourth (`n = 0.3`) stays
  permissive — an index below 1 is a real metallic index, not garbage.

### Added

- 13 more rows in `validation/smoke/test_input_validation.py`, including that
  `end_idx` narrows the span, that a multiblock request keeps the wide one,
  and that both endpoints plus the tolerance band are accepted.

## [0.6.0] — `ScatterMatrix` rejects malformed input (R3.1)

**Behaviour change.** The constructor now validates its inputs and raises
`ValueError`. Inputs that used to return plausible numbers for a stack nobody
asked for are errors; the message names the offending index and value.

### Fixed — these used to be silent

| input | old behaviour |
|---|---|
| negative thickness | same numbers as **deleting the layer** |
| NaN / inf thickness | same numbers as **deleting the layer** |
| NaN or inf refractive index | every output NaN, nothing naming the layer |
| refractive index with `|n|` past `sqrt(DBL_MAX)` | `n**2` overflowed inside the solve; NaN out |
| NaN wavelength | NaN output |
| wavelength <= 0 | `k = 2*pi/lambda` divides by zero |
| duplicated wavelength | NaN `GD`/`GDD`/`TOD`/`FOD` (the kernels divide by the grid spacing) |
| descending wavelength grid | worked, but silently a second convention |
| angle > 90 deg (e.g. 120 deg) | **aliased onto its mirror** (60 deg) — only `sin(theta)` reaches the engine |
| negative or NaN angle | aliased / NaN |

Descending grids are **rejected, not sorted**: silently reordering would
desynchronize the grid from the caller's own wavelength-indexed arrays, and
nothing would say so. Sorting is the caller's call.

### Unchanged — deliberately still accepted

Grazing incidence (90 deg, `R = 1`), a single-point wavelength grid, zero
thicknesses, a 1e9 nm layer (no upper cap — verified safe), and metallic
indices with `n < 1`. A validation layer that rejected any of these would
have broken the library to fix a bug; `test_input_validation.py` holds that
line with 9 positive controls.

The native `Solver` underneath stays permissive. It is also the optimizer's
inner loop and the Rust test surface, where a re-check per call is pure
overhead — so the checks live at the Python user surface, run once per
construction, and have no opt-out flag.

### Added

- `validation/smoke/test_input_validation.py` — 30 tests: 17 rejected inputs
  (each asserting the index and value appear in the message), 9 accepted
  ones, and four positive controls, including that a rejected construction
  does not mutate the caller's arrays and that the native `Solver` is never
  constructed for a stack that is about to be rejected.

### Changed

- `validation/review/garbage_in.py` recorded the old silent behaviour as
  prose. It now carries an expected verdict per case and **exits 1** when one
  drifts (it previously computed an `OK` flag it never used and always
  exited 0). The four `needle_gradient` rows are still permissive and are
  marked `[R3.2]` rather than quietly listed.
- `validation/parity/smatrix/test_physics_mirror.py` — the Nelder-Mead
  thickness search is unbounded and its minimum sits on the `d = 0` boundary,
  so the simplex reached for negative thickness. That was silently accepted
  before, meaning the optimizer was steered by the merit of a *different*
  stack; it now gets a sloped barrier at zero.

### Known gaps

- `needle_gradient` is unguarded: a NaN needle index still gives NaN, and
  `z` outside the stack or negative is accepted silently. That is R3.2.
- `roughness_values` and `incoherent_flags` are not validated (negative sigma
  is still accepted). Not in R3.1's scope; no silent-wrong-answer path is
  known for them.

## [0.5.9] — Request-bit and schema sync guards (R2.5)

Three constant tables are written once per language with nothing tying them
together: the 49 `REQ_*` solver bits, the 18 `NREQ_*` needle bits, and two
schema versions. Writing the guard found one of them already broken.

### Fixed

- **`PROGRAM_SCHEMA_VERSION` was decorative.** `config.rs`'s `gate_document`
  matched the envelope version against the literal `Some(1)` while quoting the
  constant in its rejection message, so bumping the constant would have made
  the gate reject exactly the version it claimed to read, and accept the one
  it claimed was stale. It compares against the constant now. Demonstrated:
  with the literal in place, setting the constant to 2 left the gate accepting
  1 and the new sync test green; with the fix, the same edit fails the test
  (`python PROGRAM_SCHEMA_VERSION=1, rust accepts 2`).

### Added

- `validation/smoke/test_request_bits.py` — 65 tests, guarding each table from
  both ends:
  - **Source-level**, parsing `pub const REQ_*` / `NREQ_*` out of
    `core_engine.rs` and `needle_engine.rs` and comparing name-for-name and
    bit-for-bit against `Request` / `NeedleRequest`. Catches a constant added
    on one side only — which the plan's semantic probe alone cannot see, since
    a bit Python never mirrored is a bit Python can never request. Skips when
    the Rust tree is absent (installed wheel).
  - **Behavioural**, driving all 49 single bits through the engine and
    asserting the returned channel keys exactly, plus the six convenience
    bundles and an all-bits-at-once call. Catches a renumbering both source
    tables agree on. Proven by swapping `RS`/`RP`: 3 failures.
  - **Density/uniqueness** on both flags, and completeness of the bit ->
    output-key map (two bits claiming one key is otherwise undetectable).
  - **Schema**, probing the two native gates at runtime rather than parsing
    them — `check_schema_version` and `load_document` are thin shims over
    Rust, so the version Rust accepts is observable from Python, wheel or
    checkout. Untagged states and untagged documents must be refused.

### Known gaps

- The `REQ_*` constants are still not exported by the extension, so the
  source-level half is the only exact check and it needs the Rust tree. The
  behavioural half covers the wheel case.

## [0.5.8] — Parity tests were passing without comparing anything (R2.4)

`validation/conftest.py` had ignored the whole `parity/` tree, so nobody had
run those files under pytest. Collecting them was supposed to be a one-line
change. It uncovered two defects that made five of the seven smatrix parity
scripts worthless as tests.

### Fixed

- **The numba reference never loaded, and the scripts called that a pass.**
  Five parity scripts imported `loom_matrix` from `validation/parity/loom/`,
  a directory that does not exist; the reference actually lives in
  `validation/parity/smatrix/refs/`. The `except ImportError` branch then
  scored every comparison `"PASS (rust-only)"` — a green verdict for a run
  that compared the Rust result against nothing. The import now points at
  `refs/`, and a genuinely missing reference **skips** instead of passing.
- **Failure exited 0.** The same scripts printed `OUTPUT_STATUS FAIL` and
  returned success, so no caller — CI, a shell loop, a human — could gate on
  them. Each now ends in a shared `report()` verdict, a `test_parity()`
  assertion for pytest, and `sys.exit(1)` for direct invocation. Verified by
  injecting a 1e-3 error into `test_w_function.py`: the script exits 1 and
  pytest reports FAILED.
- **`test_needle_t_a_phi.py` collected zero tests** — it had a `main()` and
  no test function. Given a pytest wrapper that asserts `main() == 0`.
- **Two `core_engine_*` scripts crashed collection** with a module-level
  `sys.exit(1)`, which pytest reports as INTERNALERROR. They now skip with
  the real reason (see Known gaps).

### Added

- `validation/parity/_parity.py` — shared preamble for the parity layer:
  `console_utf8()` (cp1252 consoles choked on `μ`), `require_reference()`
  (skip under pytest, exit 0 when run directly), and `report()`.

### Changed

- `validation/conftest.py` collects `parity/` and ignores only what is not a
  test: `benches`, `refs/*` (the reference implementations themselves),
  `gen_*.py` fixture generators, and `golden_mirror.py`.
- Suite is now **504 passed, 2 skipped** (was 450). The plan estimated the
  directory held 11 pytest-style tests; it holds 17 files and 54 tests.

### Known gaps

- **R2.4a**: `test_core_engine_photometry_only.py` and
  `test_core_engine_rigorous_ellipsometry.py` load a standalone
  `navette_matrix` extension from `validation/parity/target/release/`. That
  crate does not exist anywhere in this repository, so the plan's remedy
  ("document its build step") is not possible — the two files need porting
  onto `navette._smatrix.core_engine` and its request mask. They are the only
  whole-engine parity oracles; every other parity file covers one kernel.
- The `edge_nc` case in `test_solve_coherent_block_fields.py` is
  self-referential: `refs/loom_matrix.py` was edited in lockstep with the
  Rust change it is supposed to check independently.

## [0.5.7] — CI fixes found by running it (R2.3)

The first live run of `ci.yml` was green on all three blocking jobs (rust,
python/windows, python/ubuntu) and surfaced two defects in the workflow
itself. Inspection would not have caught either.

### Fixed

- **The numpy-floor job failed for the wrong reason.** With no cp312 wheel,
  pip fell back to building numpy 1.22.0 from sdist and died with
  `BackendUnavailable: Cannot import 'setuptools.build_meta'` — a build-env
  error that buries the actual finding. It now installs with
  `--only-binary=:all:` and, on failure, prints the finding as a
  `::error::`: the declared floor has no wheel for the declared
  `requires-python`. (Still advisory; R4.1 fixes the floor itself.)
- **`actions/checkout@v4` and `actions/setup-python@v5` are Node 20**, which
  the runners now force onto Node 24 with a deprecation warning. Bumped to
  `@v5` / `@v6` (both `node24`) in `ci.yml` **and** `release.yml`.

### Known gaps

- `release.yml` still pins `actions/upload-artifact@v4` and
  `download-artifact@v4`. Left alone deliberately: that workflow only runs on
  a tag, so a major bump there cannot be verified before it matters, and it
  is the workflow that publishes to PyPI.

## [0.5.6] — push/PR CI gate (R2.3, §7)

### Added

- **`.github/workflows/ci.yml`** — runs on every push and pull request.
  Until now the only workflow was `release.yml`, and it only ran on a tag:
  a broken example, unsynced request bits and a tree full of live compiler
  warnings all survived the 0.5.0 release with nothing to catch them.

  Blocking: `cargo test --workspace`; `cargo check --workspace --all-targets`
  under `RUSTFLAGS=-D warnings` (possible as of 0.5.5); `pytest validation` on
  Windows **and** Linux; `tools/check_exposure.py` — the README has claimed
  this is "enforced in CI" for some time and it now is; `tools/check_cie_sync.py`;
  and an assertion that the PEP 517 install really is a release build, which
  exercises R2.1's `build_profile()` probe end to end.

  Advisory (`continue-on-error`), each with the reason recorded in the file:
  `cargo clippy` (248 findings), `cargo fmt --check` (981 files differ from
  rustfmt defaults), and the numpy-floor job (see below). Making any of them
  blocking means landing a large mechanical diff first; doing that inside the
  commit that introduces the gate would bury the gate. A CI that cannot go
  green teaches people to ignore CI.

### Fixed

- **`release.yml`'s `check-versions` could not catch the version skew that
  actually happened.** It compared the tag against `pyproject.toml`,
  `Cargo.toml` and `__about__.py` — but not `Cargo.lock` (which lagged at
  0.5.1 while `Cargo.toml` said 0.5.2, see 0.5.3) nor the internal `navette`
  path-dependency requirement in `Cargo.toml`. All six are now checked, each
  mismatch reported individually as a `::error::` rather than one opaque
  "version mismatch". Verified both ways against the working tree.

- `tools/check_exposure.py`: allowlisted `max_disp_order`. It had counted as
  exposed only because `navette-py/src/smatrix.rs` carried an **unused `use`**
  of it; deleting that dead import in 0.5.5 surfaced it. A bare `use`
  satisfies this lint, so an unused-import cleanup can legitimately turn it
  red — worth knowing before assuming a red run means a lost binding.

### Known gaps

- The **declared numpy floor is fiction**: `pyproject.toml` says
  `numpy>=1.22.0` while `requires-python` is `>=3.12`, and numpy only gained
  cp312 wheels at 1.26. The CI job that would catch this is present but
  advisory; **R4.1 must raise the floor and delete its `continue-on-error`**.

## [0.5.5] — Zero compiler warnings (R2.3 prerequisite, §7)

Prerequisite for the push/PR CI gate: `-D warnings` cannot be enabled while
the tree emits any. All 25 are gone; the build is warning-free.

The §7 inventory listed 13. `cargo check --workspace --all-targets` finds 25 —
the extra 12 are four `unused_mut`, one never-read enum field, one
non-snake-case binding, two private-interface warnings, and four pyo3
deprecations. As with R1.1's site list, the inventory was a lead, not a census.

### Fixed

- **Two dead computations, not just dead names.** `solver::find_minima` copied
  the whole landscape into `land_vec` and bound `land` to it; neither was ever
  read (the median is taken from a second copy). The copy is gone — one fewer
  full-grid allocation per eigenmode scan. `IntGap::EdgeMean`'s first field
  (the op mean) was never read either: the edge form uses `d_bar`, re-derived
  from `opm` at the consumption site. The variant now carries only the
  bandwidth.
- Unused imports (`NeedleTargets`, `cplx`, `PI`, `max_disp_order`,
  `rayon::prelude`, `ArraySeed`, two glob imports), unused bindings (`land`,
  `m` ×3, `wavelengths`, `num_angles`, `ok`), and four redundant `mut`.
  `NeedleTargets` is used only by `cycle.rs`'s tests, so it moved into the
  test module rather than being deleted.
- `PyFilmInput` was private while appearing in the `pub(crate)` signatures of
  `assemble_design` and `run_design`; it is now `pub(crate)` to match.
- pyo3 0.28 deprecation: `#[pyclass]` types deriving `Clone` get an automatic
  `FromPyObject` that is becoming opt-in. The five `config_type!` classes and
  `PyLayerSpec` now declare `from_py_object` explicitly, which **preserves
  today's behavior** rather than silently losing the conversion at the next
  pyo3 bump.

Behavior-neutral by construction, and verified: eigenmode scan / coarse-minima
/ landscape checksums are **bit-identical** across the change (three
resolutions × three median factors, built before and after and diffed), and
450 pytest + 376 cargo tests pass.

## [0.5.4] — Build-profile probe + bench guard, UTF-8 benches (R2.1, R2.2, §5.1.1, §17)

### Added

- **`navette.build_profile()`** — reports the Cargo profile the installed
  extension was compiled with, `"release"` or `"debug"`. A debug build imports
  and computes correctly but runs several times slower, so nothing detects it
  by eye; this makes it a one-line check.
- **`validation/benches/_bench_common.py`** — shared bench preamble.
  `require_release()` exits rather than time a debug build; `setup_bench()`
  reconfigures stdout/stderr to UTF-8 and puts the repo's `src/` on
  `sys.path`; `bench_provenance()` returns the profile/version/platform block
  now stamped into bench results JSON.
- All **seven** benches (`smatrix/bench_backside_speed.py`,
  `structure/bench_grid_assert.py`, `synthesis/bench_refold.py`,
  `spectralweave/navette_{spectral,target}_bench.py`,
  `interpolate/1dinterpol_test_bench.py`, `color/bench_validate_color.py`)
  call the preamble before importing `navette`.
- `validation/smoke/test_build_profile.py` (12 tests) — asserts the gate
  refuses `"debug"` and an unbuilt extension, that the package-level accessor
  matches the native symbol, and that **no bench can be added without
  `require_release()`**. It deliberately does *not* assert the installed
  profile: testing against a debug build is fine, *timing* one is not.

### Fixed

- **Benches died with `UnicodeEncodeError` on a cp1252 console** (§17) — the
  interpolate, target and color benches print `->` arrows and emoji and needed
  `PYTHONIOENCODING=utf-8` to run at all. `setup_bench()` fixes this at the
  source; verified by running every bench with the variable unset.
- Benches that did `sys.path.insert(0, "src")` only worked when launched from
  the repo root; they now resolve `src/` from `__file__`.

### Changed

- **Every `maturin develop` build instruction now says `--release`** — the
  README, `navette/__init__.py`, and the sixteen `ImportError` build hints in
  the wrapper modules. The README carries an explicit warning; following any
  of the old hints produced exactly the debug build of §5.1.1.
- Committed bench results (`validation/benches/smatrix/results/*.json`)
  predate the gate and are of **unknown profile** — `bench_backside_speed.py`'s
  `A_legacy` case measures ~0.06 ms on a release build against the 0.50 ms
  recorded there, so those files are almost certainly debug-build numbers and
  must not be compared against new runs. Noted in `validation/README.md`.
- `tools/check_exposure.py`: allowlisted `nevot_croce_factors` (internal, via
  every solve/field/needle path). The lint had been **failing since 0.5.3** —
  it is not wired into any workflow yet, which is R2.3.

## [0.5.3] — Névot-Croce: remaining two sites + single source of truth (R1.1, §3.2)

### Fixed

- **Two further Névot-Croce sites still applied the reflection factor to
  transmission.** R1.1 (0.5.1) fixed the two sites the remediation plan listed;
  the type-5 branch actually exists in **four** places. Still broken were:

  | Site | Reached from |
  |------|--------------|
  | `solver::field_prof` | `ScatterMatrix.field_profile()` |
  | `needle_operator::interface_matrix` | `needle_gradient()` / needle synthesis |

  The needle case mattered most: synthesis merit is evaluated through
  `coherent_block` (fixed in 0.5.1) while its needle sensitivity ran through
  `needle_operator` (not fixed), so the P-function stopped being a derivative
  of the merit it optimizes. Measured on a 21-layer σ = 4 nm stack, the needle
  gradient moved by **20.9 %**; the field profile by 0.08 %.

### Changed

- **`optics_core::nevot_croce_factors()` is now the single source of truth**
  for the type-5 factors. All four interface builders call it instead of
  open-coding the exponentials, which is what allowed a partial fix to land.
  Bit-identical output on every path (verified by checksum: roughness type 0
  is unchanged everywhere, and `coherent_block` type 5 is unchanged).
- `RoughnessType` and `ScatterMatrix.energy_conservation()` docstrings now
  state the observable consequence of Névot-Croce's perturbative unitarity:
  **`T` can exceed 1 and the residual `A = 1 - R - T` can go negative** outside
  `|kz·σ| ≪ 1`; neither is clamped, and `energy_conservation()` reports a
  magnitude, so an energy *gain* is indistinguishable from absorption.
  Measured: n 1 → 4.28, σ = 20 nm, 550 nm gives T = 1.0768, A = −0.2347.
- `Cargo.lock` was one version behind (it recorded 0.5.1 while `Cargo.toml`
  said 0.5.2). `release.yml`'s `check-versions` greps `pyproject.toml`,
  `Cargo.toml` and `__about__.py` but never the lock, so this passed CI
  silently while dirtying the tree on any `cargo build`.

### Added

- `validation/smoke/test_rtype5_cross_path.py` (15 tests) — makes a partial
  fix impossible to repeat, three independent ways:
  - **index-matched invisibility**: a type-5 interface between two identical
    media must be bit-identical to no interface at all (r = 0 and the
    transmission factor collapses to exactly 1), checked on `compute()`,
    `field_profile()` and `needle_gradient()` at zero tolerance;
  - **needle ↔ solver agreement**: `2·P_T` must equal a central difference of
    `ΣT²` taken through the ordinary solver, to 1e-6, on a rough stack;
  - **source guard**: every `rtype == 5` branch must call the shared helper,
    and the reflection exponential may be written out in `optics_core.rs` only.

  Verified to have teeth: reintroducing the old factor fails 11 of the 15,
  while the pre-existing 423-test suite passes on that same broken build.

## [0.5.2] — energy_conservation() accepts 2-D input (R1.2, §15)

### Fixed

- **`ScatterMatrix.energy_conservation()` no longer raises `TypeError`.** The
  wrapper feeds 2-D `[n_angles, n_wavs]` arrays to the native
  `solver_energy_conservation`, which previously accepted only 1-D input, so
  the advertised API failed on every call. The binding now accepts 2-D, checks
  shape agreement, and reshapes the elementwise residual back to 2-D; the
  single-angle path still squeezes to 1-D.

### Added

- `validation/smoke/test_energy_conservation.py`: lossless (~0 residual),
  absorbing (residual == absorptance, monotone in thickness), 2-D shape, and
  1-D squeeze coverage.

### Changed

- `test_package_version` validates semver format instead of pinning a literal,
  so it survives per-item version bumps (release workflow still enforces sync).

## [0.5.1] — Névot-Croce transmission factor (R1.1, §3.2)

### Fixed

- **Energy conservation in the Névot-Croce roughness path (roughness type 5).**
  The type-5 branch previously applied the reflection Debye-Waller factor
  `f = exp(-2·kz1·kz2·σ²)` to the transmission amplitudes as well, silently
  destroying transmitted energy at rough interfaces. Transmission now takes the
  canonical wavevector-difference factor `exp(+((kz1-kz2)·σ)²/2)` (Névot & Croce
  1980; de Boer, Phys. Rev. B 44, 498 (1991); Stearns, J. Appl. Phys. 65, 491
  (1989)), so `R + T → 1` for lossless interfaces.

  Measured single interface, n = 1.0 → 1.5001, λ = 550 nm, normal incidence
  (no absorption — energy must be conserved):

  | σ (nm) | old R+T (buggy) | new R+T (fixed) |
  |--------|-----------------|-----------------|
  | 0      | 1.000000        | 1.000000        |
  | 3      | 0.992977        | 1.000001        |
  | 6      | 0.972202        | 1.000016        |
  | 10     | 0.924678        | 1.000125        |

  **This changes numerical results for every stack using roughness type 5**
  (XRR is the main consumer). Fixed at both `coherent_block.rs` call sites; the
  in-repo loom parity mirror was patched in lockstep. **Incomplete** — the
  remediation plan's site inventory said "two verbatim sites" but there are
  four; `solver::field_prof` and `needle_operator::interface_matrix` were
  missed and are fixed in 0.5.3.

### Added

- `validation/smoke/test_rtype5_energy.py`: energy-conservation sweep
  (rtype × σ × contrast × angle), a regression guard against re-damping
  transmission, and an independent closed-form Névot-Croce literature anchor
  checking engine R and T individually to machine precision.
