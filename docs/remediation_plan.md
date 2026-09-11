# Navette Remediation Plan

Derived from `docs/code_review.md` (24 sections, current at commit `b567d26`).
This document is the executable counterpart to the review: every finding that
requires a change is turned into a work item with a precise fix design, a
validation strategy built on the **existing** test suites plus the **review
harnesses** (`validation/review/`), and the option of additional tests where
coverage is thin.

---

## 0. How to use this plan

### 0.1 Ground rules (apply to every item)

1. **Release build only.** Every build/test/bench action uses
   `maturin develop --release` (see R2.2). Before any timing claim or any
   fix that touches hot code, verify the installed extension is not the
   debug build: `python -c "import navette; print(navette.build_profile())"`
   must print `release` (R2.1, landed 0.5.4). The benches enforce this
   themselves and exit on a debug build.
2. **Suite commands.**
   - Python: `.venv/Scripts/python.exe -m pytest validation` (656 passed +
     2 skipped at 0.6.5, parity included; `PYTHONIOENCODING=utf-8` no longer
     needed after R2.2/R2.4)
   - Rust: `cargo test --workspace` (376 tests today)
   - Parity: the pytest-style parity tests after R2.4 makes them collectable
   - Review harnesses: `for f in fd_step1 fd_rchannel color_merit_check color_grad_python garbage_in weaver_race tauc_check kk_validate kk_conv; do python validation/review/$f.py; done` (all must print `ALL OK` unless the item explicitly changes documented behavior)
3. **One item = one commit** (or a small series), message referencing the
   item ID and the review section. Never mix behavior changes with refactors.
4. **Behavior changes are releases.** Fixes that alter numerical results
   (R1.1) or reject previously accepted input (R3.x) ship together in a minor
   version bump (0.6.0) with a changelog entry per item. Pre-1.0, breaking is
   acceptable; silent is not.
5. **Bit-exactness is a contract.** Anything that must not change numbers
   (refactors R6.x) is validated by the random-stack differential
   bit-identity test plus the review harnesses. Any new parallel **reduction**
   must keep per-destination-indexed writes (§13 determinism note).

### 0.2 Priority, impact, risk, effort

| Field | Meaning |
|---|---|
| **Priority** | P0 = wrong numbers or crashes shipped to users; P1 = silent garbage in edge cases / missing guard rails; P2 = process, tooling, perf; P3 = hygiene, docs, optional refactors |
| **Impact** | who is affected and how badly (physics correctness > API crashes > robustness > speed > docs) |
| **Risk** | probability the fix itself breaks something (golden churn, API breakage, perf regression) and the mitigation |
| **Effort** | S < 1 h · M = hours · L = 1–2 days · XL = multi-day |

### 0.3 Master item table

| ID | Item | Priority | Impact | Risk | Effort | Review § |
|---|---|---|---|---|---|---|
| R1.1 | Névot-Croce (rtype 5) transmission factor | **P0** | physics correctness | M (golden/parity churn) | S–M | §3.2 |
| R1.2 | `energy_conservation()` TypeError | **P0** | advertised API broken | S | S | §15 |
| R2.1 | `build_profile()` probe + bench guard | P0-enabler | benchmarking hygiene | S | S | §5.1.1 |
| R2.2 | README `--release` + bench UTF-8 | P0-enabler | docs/benches | S | S | §5.1.1, §17 |
| R2.3 | push/PR CI workflow | **P0-enabler** | gates everything after it | M (fix live warnings first) | M | §7 |
| ~~R2.3a~~ | ~~clippy clean → blocking gate~~ | — | **DONE (0.6.6)** | — | — | §7 |
| R2.3b | rustfmt adoption → blocking gate — **DEFERRED**, decision pending | P3 | 981 files; style decision first — put to the maintainer at 0.6.12, answer: not now | S (blame churn) | S/L review | §7 |
| R2.4 | parity tests collected; `sys.exit` → skip | P1 | test suite honesty | M (env dependency) | M | §15, §6.3 |
| ~~R2.4a~~ | port `test_core_engine_*` onto `core_engine` — **DONE (0.6.11)** | P1 | only whole-engine parity oracles — both restored, 13/13 channels, ~1e-14 | M | M | §15 |
| R2.5 | request-bit + schema sync tests | P1 | prevents silent corruption | S | S | §4.3, §9.3 |
| R3.1 | `ScatterMatrix` input validation | P1 | silent-garbage class closed | M (behavior change) | M | §21.2, §19.3 |
| R3.2 | needle z-range: debug-assert → error | P1 | release-build garbage | S | S | §21.1 |
| R3.3 | eigenmode `char_func`/`refine_mode` bound | P1 | silent n_eff = −1.7e8 | S–M | M | §16 |
| R4.1 | numpy floor → `>=2.0` | P1 | broken installs | S | S | §7 |
| R4.2 | color gradients on Python needle path — DONE (0.6.4) | P2 | documented flow incomplete | M (binding change) | M | §20.2 |
| R4.3 | fix + test all `examples/` — DONE (0.6.5) | P1 | shipped example broken | S | S | §6.2 |
| ~~R4.4~~ | Optimizer backends: hardened built-in LM + optional ecosystem solvers — **DONE** | P2 | R4.4b (0.6.7), R4.4c (0.6.10), R4.4c-argmin (0.6.16). Five backends; the two argmin ones are baselines that measured badly and are documented as such (see corrections) | M (pinned optima may shift) | M–L | §3.6, §18.2 |
| ~~R4.5~~ | ~~Analytic Jacobian for the refold optimizer (deposit chain, FD fallback)~~ | P2 | **DONE (0.6.8 merit rows, 0.6.9 deposits + J)** | M (fold-kink semantics; ordered accumulation) | M | §3.6, §19.1, §20.2 |
| ~~R4.6~~ | TRF backend (trust-region-reflective) — **DONE (0.6.12)** | P2 | correct boundary behavior; direct scipy parity; retires the clamp-prediction caveat — `optimizer="trf"`, 10/10 breaks caught | M–L (largest algorithmic lift; scipy as oracle) | L | §3.6 |
| ~~R5.1~~ | `unweave_collection` batch optimization — **DONE (0.6.13)** | P2 | 1.15×–2.29× where the cost is per-fragment; the rest is DRAM bandwidth, and the reference wins by aliasing the caller's buffer (see corrections) | M | L | §5.3 |
| ~~R5.2~~ | parallelize serial derive loop — **DONE (0.6.15)** | P3 | 1.42× at 20k and 1.73× at 60k on a twelve-channel request; ~1.03× on a four-channel one, which is the whole story (see corrections) | M (bit-identity) | M–L | §5.4 |
| ~~R5.3~~ | `core_engine` — **DONE (0.6.14)**, and it was not the emit path | P2 | 1.4×–3.0× across the grid; the cost was `Solver::new` (see corrections). Now ~0.8–1.2× numba at 20k–60k points, still ~0.3–0.5× at 500 | M (must stay bit-identical) | M–L | §5, R2.4a |
| R6.1 | `needle_gradient` refactor | P3 | cyclomatic 96, 7× copy-paste | M (must stay bit-exact) | XL | §4.1 |
| R6.2 | small physics nits batch | P3 | DOP_R clamp, docstring, `+0.0` | S | S | §3.3, §19.3 |
| R6.3 | solver triplication (optional) | P3 | maintenance | M (perf-sensitive) | L | §4.2 |
| R6.4 | docs/hygiene batch | P3 | audit-trail rot | S | M | §6.3, §24.2, §18.4 |
| R6.5 | `.pyi` stubs for `_smatrix`/`_spectralweave` | P3 | IDE/mypy coverage | S | M | §24.1 |
| R6.6 | Rename `color/func_NN.rs` → descriptive module names | P3 | readability; the mod.rs doc comment is currently the only decoder ring | S (internal paths only) | S–M | §11 |

Recommended execution order = the phase order below. R1.1 and R1.2 come
first despite the CI item being "enabling": the physics bug is the single
highest-value change in the repo and its validation does not depend on CI.

---

## 1. Phase 1 — Physics correctness (P0)

### R1.1 Névot-Croce (roughness type 5) transmission factor

**Review:** §3.2 (measured R+T = 0.924678 at σ = 10 nm, n: 1.0 → 1.5001, λ = 550 nm).

**What.** The type-5 branch applies the reflection Debye-Waller factor
`f = exp(−2·kz1·kz2·σ²)` to the transmission amplitudes as well. Canonical
Névot-Croce corrects **r** with `exp(−2 kz1 kz2 σ²)` and **t** with
`exp(±(kz1 − kz2)²σ²/2)` (low-contrast → 1), so transmission must not take
`f`.

**Where.** ~~Two verbatim sites:~~ **CORRECTION (0.5.3): four sites, not two.**
This inventory was wrong and the 0.5.1 fix was consequently partial — see the
0.5.3 CHANGELOG entry. Treat every `**Where.**` block in this plan as a lead to
verify, not a complete list: `grep` for the construct before calling an item done.

- `rust/navette/src/smatrix/coherent_block.rs:128-131` (single-pol `solve_pol_specialized`) — fixed 0.5.1
- `rust/navette/src/smatrix/coherent_block.rs:319-323` (dual `solve_coherent_block_fields_dual`, the `rg_t` tuple member) — fixed 0.5.1
- `rust/navette/src/smatrix/solver.rs:1687` (`field_prof`, behind `ScatterMatrix.field_profile()`) — **missed**, fixed 0.5.3
- `rust/navette/src/smatrix/needle_operator.rs:145` (`interface_matrix`, behind `needle_gradient()` and needle synthesis) — **missed**, fixed 0.5.3

All four now delegate to `optics_core::nevot_croce_factors()`; a source-level
guard in `validation/smoke/test_rtype5_cross_path.py` fails if a fifth appears.

Current code (both sites):
```rust
} else if rtype == 5 {
    let kz1 = two_pi_lam * n_curr * cos_curr;
    let kz2 = two_pi_lam * n_next * cos_next;
    let f = (-2.0 * kz1 * kz2 * sigma * sigma).exp();
    (r12 * f, r21 * f, t12 * f, t21 * f)   // site 1; site 2: (f, f, f)
```

**Fix design.**

```rust
} else if rtype == 5 {
    let kz1 = two_pi_lam * n_curr * cos_curr;
    let kz2 = two_pi_lam * n_next * cos_next;
    let f  = (-2.0 * kz1 * kz2 * sigma * sigma).exp();
    // Névot-Croce transmission factor: exp(+((kz1-kz2)σ)²/2) per amplitude.
    // ≥ 1 (compensates the r damping), → 1 in the low-contrast limit.
    let d = (kz1 - kz2) * sigma;
    let ga = (d * d * 0.5).exp();
    (r12 * f, r21 * f, t12 * ga, t21 * ga)   // site 2: (f, f, ga)
}
```

Sign choice is settled by arithmetic, not taste. For the review's repro
interface: `f_amp = e^(−0.039144) = 0.96161`, `ga_amp = e^(+0.0016318) =
1.001633`. Intensities: `R = 0.03700`, `T = 0.96304` → **R+T = 1.000041**
(residual 4.1e-5 — NC is a perturbative factor model). The `−` sign gives
R+T = 0.99395 (6e-3 off): the `+` sign is the canonical, energy-compensating
form. Encode this expectation in the test.

**Impact.** Every rtype-5 stack currently loses up to 7.5 % of transmitted
energy at a single interface (worse with σ, better with contrast). XRR — an
advertised use case — is the main consumer. All shipped numbers for rtype-5
stacks change.

**Risk & mitigation.**
- *Parity/golden churn*: the loom reference has the same bug, so mirror and
  parity tests pinning rtype-5 numbers will fail **by design**. Protocol:
  1. Apply the Rust fix.
  2. `cargo test --workspace` → triage every failure into "pins the old bug"
     (update pinned values with a comment referencing R1.1) vs "unexpected"
     (investigate; stop).
  3. Patch the in-repo mirror `validation/parity/smatrix/refs/loom_matrix.py`
     with the identical corrected factor (comment: *matches corrected NC
     physics; intentionally diverges from historical pins*). Python↔Rust
     parity stays meaningful (implementation parity); *physical* correctness
     is guarded by the new independent test below, replacing the old oracle.
  4. Grep Python-side goldens for rtype-5 stacks (materials goldens are
     unaffected — no roughness; spectralweave goldens: verify). Never
     regenerate goldens *from* the buggy reference (the §18.2 trap).
- *Users' results change*: ship in 0.6.0 with a changelog entry showing the
  before/after repro table.

**Validation.**
- **Existing:** full `cargo test` + `pytest validation` after the churn
  protocol; review harness `fd_rchannel.py` (roughness-free, must stay
  bit-exact); differential bit-identity test (same fix both sides of the
  mirror, still bit-identical).
- **New — required:** `validation/smoke/test_rtype5_energy.py`
  (pytest, auto-collected):
  - sweep `rtype ∈ {0..5} × σ ∈ {0,3,6,10 nm} × contrast ∈ {1.0→1.5001, 1.0→2.35, 1.0→Si(4.28)} × θ ∈ {0°,30°,60°}`,
    assert `max|1−R−T| < tol` with `tol` calibrated on the *fixed* build
    (start 5e-4 at near-matched, ≤ 2e-3 at high contrast), tightened with 2×
    headroom and documented in the docstring;
  - assert the pre-fix failure mode is covered: the test must fail by ≥ 10×
    tol if the old `t12*f` line is restored (guard against regression);
  - assert the worked example above numerically (R+T = 1.00004 ± 1e-4);
  - assert type-4 (Gaussian w-function path) unchanged (R+T ≈ 0.9948 at
    σ=10 — documented graded-profile approximation, §3.2 note).
- **New — optional:** a one-paragraph doc note on `RoughnessType` separating
  "graded profile (types 1–4, per-side factors)" from "NC decorrelation
  (type 5, correlated)" — folds into R6.4 or ships with this item.

**Definition of done.** Fix at both sites; mirror patched; new energy test
green on the fixed build and red on the old one; full suites green; changelog.

**Effort.** S–M (the edit is minutes; the churn triage is the work).

---

### R1.2 `ScatterMatrix.energy_conservation()` — guaranteed TypeError

**Review:** §15 (wrapper `src/navette/smatrix/smatrix.py:393-401`; binding
`rust/navette-py/src/smatrix.rs:368-381` 1-D-only; `_energy_conservation` is
the native fn imported at `smatrix.py:49`).

**What.** The wrapper hardcodes `squeeze=False` (2-D `[n_angles, n_wavs]`
arrays) and feeds them to a binding that only accepts `PyReadonlyArray1<f64>`
→ `TypeError` on **every** call. Zero test coverage (§15).

**Fix design (binding side, preferred).** Accept 2-D in
`solver_energy_conservation`: take `PyReadonlyArray2<f64>` for rs/rp/ts/tp,
assert equal shapes, flatten (contiguous) through the existing elementwise
engine fn, reshape the result to `[na, nw]`, return a 2-D `PyArray1`→
`PyArray2`. The wrapper then works unchanged and 2-D semantics are preserved
(its `_squeeze` at the end gives 1-D for 1×N). Add a debug assert on shape
equality with a clear message.

**Impact.** The README advertises the energy-conservation feature; it cannot
be called at all today.

**Risk.** S — additive binding change; no other caller (grep to confirm).

**Validation.**
- **Existing:** none exist for this path (that's the bug's habitat).
- **New — required:** `validation/smoke/test_energy_conservation.py`:
  - lossless stack (3-layer oxide) → `max(cons) < 1e-12`, both polarizations;
  - absorbing stack (k = 0.05 film) → `cons ≈ A` from `compute()` and ∈ [0,1],
    monotone in thickness;
  - 2-D input: shape `[n_angles, n_wavs]` preserved; 1-D path (single angle
    via `squeeze=True` request) still returns 1-D;
  - the R1.1 sweep file reuses this helper for its invariant (shared util or
    duplication — prefer duplicating 5 lines over new coupling).

**Effort.** S.

---

## 2. Phase 2 — Guard rails: build hygiene, CI, test honesty (P0-enablers)

### R2.1 Build-profile probe + bench guard — DONE (0.5.4)

**Review:** §5.1.1 (the venv shipped a debug build; every timing was wrong).

**Fix design.**
1. New pyfunction in `rust/navette-py/src/lib.rs` (or `smatrix.rs` module
   root): `#[pyfunction] fn build_profile() -> &'static str { if cfg!(debug_assertions) { "debug" } else { "release" } }`;
   register on the top-level `navette` module. Export the same string in the
   `__about__`/docs note.
2. New `validation/benches/_bench_common.py`:
   ```python
   import sys
   def setup_bench():
       for s in (sys.stdout, sys.stderr):        # cp1252-safe benches (R2.2)
           try: s.reconfigure(encoding="utf-8", errors="replace")
           except AttributeError: pass
   def require_release(profile_fn):
       p = profile_fn()
       if p != "release":
           sys.exit(f"REFUSING TO BENCHMARK: extension built with '{p}' "
                    f"profile. Run: maturin develop --release (§5.1.1).")
   ```
   Import + call at the top of every bench
   (`smatrix/bench_backside_speed.py`, `synthesis/bench_refold.py`,
   `spectralweave/navette_{spectral,target}_bench.py`,
   `structure/bench_grid_assert.py`, `interpolate/1dinterpol_test_bench.py`,
   `color/bench_validate_color.py`).
3. Stamp provenance into committed results JSON
   (`benches/smatrix/results/*.json`): add `"build_profile"` key on the next
   regeneration; note the historical files are of unknown profile.

**Impact.** Makes the §5.1.1 situation structurally impossible; benches
self-identify their build.

**Risk.** S — one new exported symbol; benches hard-fail (intended).

**Validation.** Existing: benches run green under `--quick` after rebuild.
New: trivial — assert `build_profile() == "release"` on the dev venv *after*
R2.2's README rebuild; do **not** add a pytest assert on the profile (would
fail for every contributor with a debug build — the check belongs at bench
time, not test time).

**Effort.** S.

### R2.2 README + UTF-8 benches — DONE (0.5.4)

- `README.md:93`: `maturin develop` → `maturin develop --release`, with a
  one-line warning that plain `maturin develop` produces a debug build that
  invalidates every benchmark (§5.1.1).
  **CORRECTION (0.5.4): seventeen sites, not one.** The same instruction is
  repeated in `navette/__init__.py` and in every wrapper's `ImportError`
  build hint (`color/`, `interpolate/`, `materials/`, `spectralweave/`,
  `structure/materials.py`, `smatrix/{needle,smatrix}.py`, and the six
  `_*.py` shims). A user following any of them landed on the debug build, so
  fixing only the README would have left the trap in place.
- Bench scripts: route prints through `_bench_common.setup_bench()` (or add
  `sys.stdout.reconfigure(...)` inline) — fixes the cp1252
  `UnicodeEncodeError` in `1dinterpol_test_bench.py` and the target/color
  benches (§17). Either keep the emoji (now safe) or strip them; prefer keep.
- Re-run all benches on the release build and refresh §5.2's committed JSON
  provenance (optional, ties to R6.4's README perf re-audit).

**Validation.** Existing: full bench suite runs (UTF-8 guard verified by
running one bench in a `cmd.exe` console without `PYTHONIOENCODING`).
**Effort.** S.

### R2.3 push/PR CI workflow — DONE (0.5.6; warning cleanup 0.5.5)

**Review:** §7 (only tag-triggered `release.yml` exists; README claims
exposure lint is "enforced in CI" — it is not).

**Fix design.** New `.github/workflows/ci.yml` on push + PR:
- job `rust`: `cargo fmt --check`; `cargo clippy --workspace --all-targets -- -D warnings`;
  `cargo test --workspace`.
  - **Prerequisite:** fix the live warnings first (§7 inventory: unused
    imports `NeedleTargets`, `cplx`, `PI`, `max_disp_order`,
    `rayon::prelude`, `ArraySeed`; unused vars `land`, `m`×3, `wavelengths`,
    `num_angles`, `ok`). Do this as its own commit *before* enabling
    `-D warnings`.
    **CORRECTION (0.5.5): 25 warnings, not 13.** The extra 12 are four
    `unused_mut`, one never-read enum field (`IntGap::EdgeMean.0`), one
    non-snake-case binding, two private-interface warnings (`PyFilmInput`)
    and four pyo3 `FromPyObject` deprecations. Done in 0.5.5, and two of
    them were dead *computations*, not dead names.
  - **What actually shipped as blocking (0.5.6):** `cargo test`,
    rustc `-D warnings`, `pytest validation` (Windows + Linux),
    `check_exposure.py`, `check_cie_sync.py`, and a `build_profile() ==
    "release"` assertion. **`cargo clippy -- -D warnings` and
    `cargo fmt --check` are ADVISORY**, because the tree has 248 clippy
    findings and 981 files that differ from rustfmt defaults — neither was
    in the plan's estimate, and both are large mechanical diffs that must
    not ride along with the gate itself. Two follow-up items below.
- job `python` (windows-latest + ubuntu-latest): install Python 3.12 +
  numpy≥2 → `maturin develop --release` → `pytest validation` (post-R2.4
  scope) → `python tools/check_exposure.py` → `python tools/check_cie_sync.py`.
- job `numpy-floor`: create a fresh venv, `pip install -e .` with `numpy>=2`
  pinned, `python -c "import navette; print(navette.__version__)"` (guards R4.1).
- `release.yml` gains `needs: [ci]` (or the equivalent gate) so publishes can
  only follow a green tree.

**Impact.** Closes the gap that let a broken example (§6.2), unsync'd bits
(§4.3), and live warnings survive a 0.5.0 release.

**Risk.** M — the parity tests (R2.4) may need the standalone
`navette_matrix.pyd`; if it cannot be built in CI cheaply, keep those
specific tests skipping loudly (`importorskip`) and say so in the job summary.
Windows runner minutes are cheap at this project's size.

**Validation.** Existing: CI green = suites green. New: the workflow itself;
test by opening a PR with a deliberate failure (e.g. revert one warning fix)
and watching it gate.

**Effort.** M.

### R2.3a Clippy clean → make the clippy gate blocking — DONE (0.6.6)

**Discovered by R2.3 (0.5.6).** `cargo clippy --workspace --all-targets`
reports **248** findings (217 in `navette`, 31 in `navette-py`), of which
clippy offers to auto-fix ~138. The CI job runs clippy with
`continue-on-error: true` today.

**Fix design.** Triage in passes, each its own commit, bit-exactness checked
per ground rule 5: (1) `cargo clippy --fix` for the mechanical ones, reviewed
hunk by hunk — auto-fixes in numeric code are *not* automatically
behavior-neutral; (2) hand-fix or `#[allow(...)]`-with-rationale the rest;
(3) delete `continue-on-error` from the clippy step in `ci.yml`.

**Risk.** M — touches hot physics code. Nothing here is a bug fix, so any
numeric change is a regression by definition.
**Effort.** M/L.

**CORRECTIONS / NOTES (0.6.6).**

* **267 findings, not 248.** The count in this plan was taken at 0.5.6; six
  releases of new code landed since. Final state: 0.
* **108 of them — 40% of the whole backlog — came from one generated file.**
  `color/matrices.rs` is emitted by `validation/benches/color/gen/gen_matrices.py`
  with `f"{x:.17e}"`, which names every double but writes ~5 more digits than
  most of them need. Fixing the *output* alone would have reintroduced all 108
  on the next regeneration and quietly un-greened a blocking gate, so the
  generator was changed too (`repr(float(x))` — shortest round-tripping form).
  All 108 literals verified bit-identical via `struct.pack("<d", ...)` against
  the pre-change tree.
* **`neg_cmp_op_on_partial_ord` must never be blanket-rewritten.** All 15 sites
  are validators shaped `!(x > 0.0)`. The negation is what rejects NaN;
  clippy's `x <= 0.0` suggestion would wave NaN through every one of them.
  Allowed crate-wide with that written in, next to the attribute.
* **Remaining allowances are crate-level, not per-site**, so they read as a
  short list of decisions in `lib.rs` rather than as noise at 40 call sites:
  `too_many_arguments` (owned by R6.1 / R6.3), `needless_range_loop` (flat
  angle-major index arithmetic), `type_complexity`, `should_implement_trait`.
* **Ground rule 5 held throughout.** The bit-exactness fingerprint
  `a99e83835e8b44ae3416a715ee94410b6123f0ca62aa4610ef8e37b3ab00f154` (12 seeded
  random stacks × every `Request` bit + needle gradients + an eigenmode
  landscape, hashed over raw f64 bytes) was unchanged after every pass, and
  `cargo clippy --fix` output was reviewed hunk by hunk rather than trusted.
* **The `--fix` pass is *not* a separate commit per pass as the plan proposed.**
  Splitting mechanical auto-fixes from their hand-reviewed corrections would
  have left intermediate commits that do not build or do not pass the
  fingerprint. It ships as one reviewed change with the gate flip.

### R2.3b rustfmt adoption → make the fmt gate blocking — DEFERRED (decision pending, asked 0.6.12)

**Discovered by R2.3 (0.5.6).** `cargo fmt --all --check` reports diffs in
**981** files — i.e. the tree has never been rustfmt-formatted and the
codebase's own style differs from rustfmt defaults in places (notably the
2-space-indented modules under `synthesis/`).

**Decision needed before doing anything:** adopt rustfmt defaults (one
tree-wide reformat commit, which rewrites `git blame` for the whole crate —
mitigate with `.git-blame-ignore-revs`), or add a `rustfmt.toml` that encodes
the existing house style and reformat to *that*. Do not start until this is
chosen; the two produce very different diffs.

**Effort.** S to do, L to review — hence its own item.

### R2.4 Parity tests collected; `sys.exit(1)` → skip — DONE (0.5.8)

**Review:** §15 (`validation/conftest.py` `collect_ignore = ["parity", "benches"]`
hides 11 pytest-style parity tests — the strongest oracle in the repo; two
NEEDS-PORT files call `sys.exit(1)` at import → pytest INTERNALERROR).

**Fix design.**
1. In `validation/parity/smatrix/test_core_engine_photometry_only.py` and
   `test_core_engine_rigorous_ellipsometry.py`: replace the import-time
   `sys.exit(1)` with `pytest.skip("NEEDS PORT: …", allow_module_level=True)`
   (keep the NEEDS PORT status in `validation/README.md`, now actually tracked).
2. `conftest.py`: replace the blanket ignore with
   `collect_ignore_glob = ["benches/*", "parity/refs/*", "parity/**/refs/*"]`
   so `test_*.py` files under `parity/` collect while numba reference scripts
   and `refs/` mirrors stay script-only.
3. Files that import the standalone `navette_matrix` parity extension: guard
   with `navette_matrix = pytest.importorskip("navette_matrix")` at module
   top so they **skip with a reason** (not error) when the pyd is absent, and
   document its build step (`maturin build --release` of the parity crate →
   `validation/parity/target/release/`) in `validation/README.md`.
4. Verify the collected count: 355 + 11 (+N skipped) ≈ 366+.

**CORRECTIONS (0.5.8).** Three things the plan did not anticipate, all found
by actually collecting the directory:

* **The parity comparisons were not running.** The numba reference was
  imported from `parity/loom/`, which does not exist; the `except ImportError`
  branch then scored every comparison `"PASS (rust-only)"`. Five "parity"
  scripts had been passing while comparing nothing. Fixed: the reference loads
  from `smatrix/refs/loom_matrix.py`, and a missing reference **skips**.
* **Failure exited 0.** The scripts printed `OUTPUT_STATUS FAIL` and returned
  success, so no caller could gate on them. Each now ends in
  `report()` + a `test_parity()` assertion + `sys.exit(1)` for direct runs.
  Proven by injecting a 1e-3 error: `rc=1` and pytest FAILED.
* **The standalone `navette_matrix` crate does not exist in this repo**, so
  step 3's "document its build step" is not possible — there is no such
  Cargo.toml. The two files skip with the real reason (see R2.4a).

Collected count went 450 → **504 passed + 2 skipped**, i.e. 54 real tests
that the blanket ignore had been hiding, not the 11 estimated.

**Impact.** The documented run stops silently excluding the best tests.

**Risk.** M — some parity tests may be slow (numba import ~seconds) or
platform-quirky; if a test is genuinely script-shaped (top-level asserts,
no test functions), leave it ignored and list it explicitly in conftest with
a comment (whitelist style beats glob style when mixed).

**Validation.** Existing: the 11 parity tests themselves. New: none needed;
update `validation/README.md`'s inventory (§6.3 rot) in the same commit.

**Effort.** M.

### R2.4a Port `test_core_engine_*` onto the request-driven `core_engine` — DONE (0.6.11)

**Discovered by R2.4 (0.5.8).** `parity/smatrix/test_core_engine_photometry_only.py`
and `test_core_engine_rigorous_ellipsometry.py` load a standalone
`navette_matrix` extension from `validation/parity/target/release/`. **That
crate does not exist anywhere in the repository** — there is no Cargo.toml for
it — so this is a port, not a build step. Both currently skip with that
reason and are visible in the run.

**Fix design.** The other five smatrix parity scripts were ported by swapping
the loader for `import navette._smatrix as rust_mod`. These two need more:
the two legacy numba kernels were merged into a single
`navette._smatrix.core_engine` driven by a request mask, so the port must map
the old positional signature onto that mask and compare the correct output
slices. The numba reference for both still exists in
`parity/smatrix/refs/loom_matrix.py`, so the oracle is available.

**Impact.** These are the only whole-engine parity oracles — every other
parity file covers a single kernel. Restoring them closes the largest gap in
the parity layer.

**Risk.** M — a wrong mask mapping produces a false FAIL, which is the safe
direction, but debugging it needs care with the ellipsometry channel order.
**Effort.** M.

**CORRECTIONS / NOTES (0.6.11) — both scripts ported; the parity layer has no
NEEDS PORT rows left.** Every channel of both legacy kernels is compared
against the numba reference, agreement ~1e-14 to ~1e-16 throughout, and both
files are collected by pytest rather than skipped.

* **The mode mapping is `FRONT_BLOCK`, and the reason is in the Rust, not in
  a guess.** The legacy ellipsometry kernel takes reflection S₂/S₃ from the
  *first coherent block's* field amplitudes and transmission S₂/S₃ from a
  Mueller cross-term product accumulated across blocks. In `core_engine.rs`
  that is exactly the `track_cross_channel == false` path (`cross_r =
  rp₀·conj(rs₀)`, `cross_t = cross_t_acc`), which is mode A.
  `COHERENCY_MATRIX` (mode B) cascades the complex p-s channel through the
  incoherent echoes and is a *later, different* treatment. Half the cases flag
  an incoherent layer so the two modes actually differ, and the script asserts
  that mode B disagrees with the reference — otherwise "mode A is the legacy
  treatment" would be an untested claim that happens to hold because nothing
  separates the modes.
* **For photometry the A/B distinction does not exist.** The coherency channel
  is only tracked when something requests it, and pure intensities never do,
  so `FRONT_BLOCK` and `COHERENCY_MATRIX` are the same computation there. The
  photometry port says so and asserts it, and asserts the distinction it *can*
  make: `FULLY_COHERENT` ignores the flags — identical to the other two when
  every layer is coherent, different as soon as one is not.
* **`calc_s`/`calc_p` became mask bits, and the old comparison would have
  passed for the wrong reason.** The legacy kernel returned four arrays always,
  zero-filling the polarization it was told to skip; the mask simply omits the
  key. Comparing the missing channel against those zeros was tried
  deliberately: **it passes**, and would keep passing if the mask were ignored
  entirely. The port asserts the key is *absent* instead, and separately that
  `Rs` from an s-only request is bit-identical to `Rs` from the full one — the
  fast path is a fast path, not different physics.
* **The 13th return value moved rather than disappeared.** `conservation_err`
  is no longer an engine channel; it is `solver_energy_conservation` (R1.2).
  The port reconstructs it from the four intensities and compares it, so the
  tuple is covered end to end and nothing is dropped on the grounds that it
  was refactored.
* **`incoherent_flags` changed length.** The binding requires one flag per
  layer; the numba kernel only ever reads the first n−1. The old test
  generated n−1 and would have been rejected outright. One array of length n
  now feeds both.
* **The comparison inputs are built once, in both layouts.** numba wants a
  complex (n_wavs, n_layers) cache, the binding wants the same numbers
  interleaved re/im and flattened wav-major. Deriving both from one source is
  what keeps a layout bug from reading as a physics difference — swapping the
  interleave was one of the deliberate breaks, and it moves `Rp` by 4e5.
* **The angle comparison is a guard, not a finding, and is labelled as one.**
  Δ is compared on the circle and the last ellipsometry case is built to sit
  on ±π (transparent, low contrast, near normal incidence — the script asserts
  it really does straddle the cut). The two engines nonetheless pick the same
  side everywhere, so a plain difference would pass this file today. Removing
  the circle handling was run as a deliberate break and was **not** caught;
  that is recorded here rather than dressed up.
* **Six deliberate breaks run, five caught**: the wrong coherence mode (Δ_R
  off by 2.8 rad), the re/im interleave swapped, `debug_flag=0` on the legacy
  call, the single-polarization identity check given the wrong mask, and the
  request mask ignored. The sixth — comparing unrequested channels against the
  legacy zeros — passes, which is the point being made above.
* **A performance finding, recorded as R5.3 rather than acted on here.** The
  request-driven engine is **slower than the numba kernel it replaced** on
  this workload: 0.5–0.65× across 500 to 60 000 grid points. A first
  measurement puts roughly half the cost outside the physics — an `Rs`-only
  solve of 20 000 points takes 0.98 ms where the full twelve-channel emit
  takes 2.08 ms, against numba's 1.94 ms for all thirteen. Mixing a
  performance change into a parity port would violate §0.1 rule 3, so it is
  its own item. **Superseded by R5.3 (0.6.14):** the 0.5–0.65× was measured
  with numba's threads still spinning, and the cost was `Solver::new`, not the
  emit. See R5.3's corrections.

### R2.5 Request-bit and schema sync tests — DONE (0.5.9)

**Review:** §4.3 (49 `REQ_*` bits hand-copied in `smatrix.py:85+`, 18
`NREQ_*` in `needle.py`; native exports the NREQ constants at
`navette-py/src/smatrix.rs:702+`; no test asserts sync), §9.3
(`PROGRAM_SCHEMA_VERSION` duplicated in `config/program.py` and `config.rs`).

**Fix design.** New `validation/smoke/test_request_bits.py`:
- **NREQ (exact):** for every exported native constant
  (`navette.NREQ_P`, …), assert `getattr(needle_module, "NREQ_" + name_tail)
  == native_value`. If the Python side builds the IntFlag from the constants
  instead of copying them (preferred long-term), this test still guards the
  module surface.
- **REQ (semantic probe):** the 49 solver bits are not exported, so probe
  semantically: for each `Request.*` member, build a single-bit request,
  `compute()` on a small fixed stack, and assert the returned dict contains
  **exactly** the expected channel key(s) for that bit (map defined once,
  49 rows). This catches renumbering, duplicates, and omissions in one test.
- **Schema:** assert `config.program.PROGRAM_SCHEMA_VERSION == 1` and that
  loading a v1 document through the native config path succeeds (a bump of
  one side without the other fails loudly here).

**Impact.** One pytest file retires an entire silent-corruption class.

**Risk.** S — the 49-row map is hand-written once; the test itself is the
documentation of the bit map.

**Validation.** Existing: none (that's the point). New: as above.

**Effort.** S.

**CORRECTIONS (0.5.9).**

* **The `NREQ_*` constants are not hand-copied.** `needle.py` imports all 18
  from `navette._smatrix` and binds them into the IntFlag, so their *values*
  cannot drift. What can drift is *membership*: a constant added in Rust and
  exported but never surfaced in `NeedleRequest` is unreachable from Python
  and nothing notices. The test checks all three tables (Rust source ->
  extension export -> IntFlag) for membership as well as value.
* **There are two schema constants across the boundary, not one.**
  `PROGRAM_SCHEMA_VERSION` (`config.rs` vs `config/program.py`) *and*
  `SCHEMA_VERSION` (`structure/version.rs` vs `structure/types.py`). Both are
  covered. Neither needed source parsing: `check_schema_version` and
  `load_document` are thin shims over the native gates, so the version Rust
  accepts is discoverable at runtime by probing.
* **The semantic probe is necessary but not sufficient.** A `REQ_*` bit that
  exists in Rust and was never mirrored in Python cannot be requested, so no
  behavioural test can reach it. The file therefore also parses the Rust
  constants and compares both directions exactly; that half skips on an
  installed wheel, the behavioural half always runs.
* **Found a live bug while writing the guard.** `gate_document` matched the
  envelope version against a literal `Some(1)` while quoting
  `PROGRAM_SCHEMA_VERSION` in its error message — the constant was
  decorative, and a bump would have rejected the version the message claimed
  to read. Fixed to compare against the constant, and proven both ways
  (constant = 2 with the literal: gate still accepts 1, test green; constant
  = 2 with the fix: test fails).
* Size: the plan said "one pytest file", risk S — accurate, but it is 65
  tests, not a handful, because the 49 bits and 6 bundles are parametrized.

---

## 3. Phase 3 — Input validation & robustness (P1)

### R3.1 `ScatterMatrix` construction validation — DONE (0.6.0)

**Review:** §21.1/§21.2 + §19.3 (dup wavelengths → NaN GD/GDD/TOD/FOD;
negative/NaN thickness silently deletes the layer; θ > 90° aliases to its
mirror; NaN indices propagate — **zero** construction-time checks; verified
in `validation/review/garbage_in.py`).

**Fix design.** Module-level `_validate_stack(n, d, wavelengths, angles)`
called from `ScatterMatrix.__init__` (the Python user surface; the Rust core
stays permissive for internal/test callers):
- `n`: every entry finite (real and imaginary) — `np.isfinite` on both parts;
  error names the layer index.
- `d`: finite, `≥ 0` (0 legal — ambient rows). Reject NaN and negative;
  **no upper cap** (1e9 nm verified safe, §21).
- `wavelengths`: finite, **strictly increasing** (unique). Descending grids
  are mathematically handled (sign-invariance verified, §21) but one rule is
  better than two; the error message tells the caller to sort ascending.
  This also protects the dispersion channels (÷Δω) and matches the
  `UniInterpolator` and color `check_grid` precedent.
- `angles`: finite, `0° ≤ θ ≤ 90°` (reject the 120°→60° silent alias).
- No opt-out flag: silent physics changes are worse than loud errors; the
  changelog lists every newly-rejected input (0.6.0).

**Impact.** Converts a family of silent-wrong-numbers paths into
actionable errors; costs microseconds at construction (not hot).

**Risk.** M — behavior change. Existing tests constructing odd inputs
(descending grids? negative-thickness probes?) must be audited: run the full
suite, fix any test that *relies* on silent clamping (they would now error —
that is the intended diff, but verify each is a test-harness artifact, not a
legitimate use case). Review harness `garbage_in.py` must be updated to
expect the new errors (it currently *documents* the silent behavior).

**Validation.**
- **Existing:** full pytest + cargo; `validation/review/garbage_in.py`
  (updated expectations); `validation/review/fd_step1.py` (its stacks are
  well-formed, must stay `ALL OK`).
- **New — required:** `validation/smoke/test_input_validation.py`,
  parametrized: dup λ / NaN λ / descending λ / NaN d / negative d / NaN n /
  inf n / θ = 120° / NaN θ — each raises `ValueError` whose message names
  index and value; plus positive controls (descending-reject ≠
  sort-mangle: assert no silent reordering happens anywhere).

**Effort.** M.

**NOTES (0.6.0).**

* **Two checks beyond the plan.** `wavelengths > 0` (`k = 2*pi/lambda`
  divides by zero at 0) and an upper bound on `|n|`: the review harness's
  "inf index" case actually passes `1e308`, which *is* finite, so the finite
  check let it through and `n**2` then overflowed inside the solve. The bound
  is `sqrt(DBL_MAX)` — the machine's, not a taste judgement: exactly where
  the engine's own arithmetic dies.
* **One test relied on the silent clamping**, as the plan anticipated:
  `test_physics_mirror.py::test_pd_optimizer_recovery` runs an unbounded
  Nelder-Mead over thickness whose minimum sits on the `d = 0` boundary, so
  the simplex reached for `d = -4`. It was a harness artifact and worse than
  it looked — the optimizer had been steered by the merit of a stack with the
  layer deleted. Fixed with a sloped barrier at zero, not by relaxing the
  check.
* **`garbage_in.py` had the parity scripts' defect too**: it computed an `OK`
  flag nothing ever set and always exited 0. It now asserts a verdict per
  case and exits 1 on drift (verified by flipping one expectation).
* Out of scope and still unvalidated: `roughness_values` (negative sigma) and
  `incoherent_flags`. No silent-wrong-answer path is known for either.

### R3.2 Needle z-range: debug-assert → real error — DONE (0.6.1)

**Review:** §21.1 (`rust/navette/src/smatrix/needle_operator.rs:383-387`
`debug_assert!(xi <= ds[j])` / `debug_assert!(xi >= -1e-9 && xi <= fields.ds[j] + 1e-9)`;
release builds silently compute garbage for out-of-range z).

**Fix design.**
1. Binding-side validation (clean error, no panic): in the
   `needle_gradient` PyO3 entry (and the eigenmode/`field_profile` entries
   that take z arrays), check every z ∈ [0, Σd_total] within `1e-9` and
   raise `PyValueError` naming the bound. Do it before `py.detach`.
2. Keep the internal `debug_assert!`s as invariants for Rust callers, but
   upgrade the message to `expect`-style wording; **do not** put a release
   `assert!` in the hot loop (panic-in-release is worse than the binding
   check for a user-facing path — the binding check is the contract).

**Impact.** Removes the last debug-only user-input check found by the audit.

**Risk.** S — additive validation; internal callers unaffected.

**Validation.** Existing: `validation/review/fd_step1.py` Part B (legal z,
bit-exact). New: extend `validation/smoke/test_input_validation.py` with
z-range rows (z = 1000 on a 400 nm stack; z = −1) asserting the precise
`ValueError`.

**Effort.** S.

**NOTES (0.6.1).**

* **The check went in the core, not the bindings.** The plan said
  binding-side, before `py.detach`. One implementation in
  `solver::needle_gradient` covers both PyO3 entry points *and* Rust callers;
  it returns `Result<_, String>`, which the bindings already map to
  `PyValueError`, so the user-visible error is identical. No panic is
  introduced either way.
* **There is no single span.** The coherent kernels are confined to
  `[start_idx, end_idx]`; the multiblock cascade calls
  `locate_depth_in(thicknesses, 0, nl - 1, z)` — the whole stack. A call can
  request both. Checking one bound for both would reject legitimate
  multiblock depths; each span is checked only when the request reaches that
  path.
* **Out-of-range z was not merely unasserted, it was silently relocated.**
  `locate_depth_in`'s `j == end_idx - 1` arm returns the last layer for any z
  past the end, with xi > d_j. The gradient came back finite and plausible.
* Also closed here (beyond the plan's scope, same check, same place): a
  non-finite `needle_n_per_wav` entry, which used to NaN the whole gradient.

### R3.3 Eigenmode `char_func` upper bound + `refine_mode` no-mode contract — DONE (0.6.2)

**Review:** §16 (`optimizer.rs` `char_func` guards the lower edge only;
unbounded Nelder–Mead returned n_eff = −1.7e8 + 7.4e6j with val ≈ 1e-217 for
a seed where no s-pol mode exists).

**Fix design.**
1. `char_func` (`optimizer.rs:56`): add the upper guard symmetric to the
   lower one. Compute a physical box from the stack it already receives:
   `n_bound = 3.0 * max(|n_i| for all layers incl. ambient)` (O(layers), no
   signature change needed); return `1e30` when `|n_eff| > n_bound`,
   `!n_eff.is_finite()`, or the propagation exponent would overflow (guard
   `Im(β)·d` magnitude before `exp`, mirroring the existing branch-safety
   pattern). This kills the `|1/r|² → 0` attractor: out-of-box values are
   maximally *bad*, so the simplex cannot sprint into overflow.
2. `refine_mode`: clamp simplex candidates into `[−n_bound, +n_bound]²`
   (projection, not rejection — keeps Nelder–Mead well-defined), and after
   convergence raise `ValueError("no eigenmode near seed … (best residual …)")
   from the Python wrapper when `val > 1e-6` (configurable threshold kwarg,
   default strict). The 1e-6 threshold is calibrated against the verified
   SPP: a true pole polishes to val ≈ 4e-19 (§16), six orders below.
3. `find_eigenmodes`/`landscape` unchanged (bounded grid already protects
   them; the clamp is inert inside the box).

**Impact.** The primary documented workflow (landscape → refine) is already
protected; this closes the manual-seed garbage path with a confident-looking
1e-217 "success".

**Risk.** S–M — touches the shared `char_func`; the SPP scenario (verified
§16, not currently pinned in the suite) must stay bit-stable. The box is
3× the max index — generous for guided modes (all physical n_eff ∈
[min|n|, max|n|]); note the bound in the docstring so multi-mode users with
high-index substrates understand the rejection region.

**Validation.**
- **Existing:** none pin eigenmodes today.
- **New — required:** `validation/smoke/test_eigenmodes.py`:
  - pin the verified SPP: 50 nm metal film (ε = −12+1j) between air and
    glass, p-pol; `refine_mode` from seeds (1.45, 0.2) and (1.3, 0.1) both
    → n_eff = 1.0459458 + 0.0015949j (tol 1e-6), val < 1e-10;
  - `field_profile` of the converged mode peaks at the film/substrate
    interface (z = 50 nm) — the textbook LR-SPP check;
  - **runaway regression:** s-pol same stack, seed in the box but no mode
    present → must raise `ValueError` (or return val > 1e-6), and
    `|n_eff| ≤ n_bound` for *any* intermediate the minimizer touched
    (assert via the returned value; optionally expose the best-seen point);
  - landscape s-pol stays flat (min ≈ 0.99).

**Effort.** M.

**CORRECTIONS (0.6.2).**

* **The pinned SPP numbers are not reproducible.** The plan asks for
  `n_eff = 1.0459458 + 0.0015949j` with `val < 1e-10` on "50 nm metal film
  (eps = -12+1j) between air and glass, p-pol" — but does not state the
  wavelength. At 632.8 nm this geometry has *two* landscape minima. The one
  near 1.046 refines to `1.0471183 + 0j` with `val = 1.8e-6`: a shallow
  resonance, not a pole — `Im(n_eff) -> 0` on a lossy metal, and it lies below
  the substrate index. Worse, `1.8e-6` is *above* the threshold the plan
  proposes, so pinning it would have made the plan's own example fail the
  plan's own default. The genuine pole is the glass-side plasmon at
  `1.7139417 + 0.0226069j`, `val = 1e-17` — eleven orders below the threshold,
  and its field profile peaks at the metal/glass interface as a bound surface
  mode must. That is what the test pins.
* **The box went into `char_func_xy`, not `char_func`.** `char_func` is shared
  with the landscape scanner, over a range the *user* chooses. Guarding it
  there would paint `1e30` across any part of a requested scan lying outside
  the box — corrupting a diagnostic the caller explicitly asked for. The
  minimizer picks its own points and is the thing that needs walls.
  `test_the_landscape_is_not_walled_off_outside_the_box` pins the distinction.
* **Step 3 was already true and is now tested**: `find_eigenmodes` returns
  `[]` for s-polarization on this stack. The documented workflow never
  produced the garbage; only the manual seed did.
* The s-polarized landscape minimum is `1.0000000000000004`, not "≈ 0.99".
* Boxing alone would not have been enough. Confined, the non-mode settles on
  the boundary with `val ~ 1e-3` — the residual contract is what turns that
  into an error instead of a number.

---

## 4. Phase 4 — API & packaging completion (P1/P2)

### R4.1 numpy floor → `numpy>=2.0` — DONE (0.6.3)

**Note (0.5.6):** `ci.yml`'s `numpy-floor` job is `continue-on-error: true`
and **expected to fail** — the declared floor `numpy>=1.22.0` cannot install
on `requires-python = ">=3.12"` (cp312 wheels start at numpy 1.26). This item
must delete that `continue-on-error` once the floor is raised, or the job
stays decorative.

**Review:** §7 (wheels built on the numpy-2 C-API cannot load against
numpy 1.x; `pyproject.toml:30` declares `numpy>=1.22.0`).

**Fix.** One line: `"numpy>=2.0"`. Also re-state the abi3 caveat in the
comment: the `numpy` Rust crate (0.28) makes no abi3 promise; the wheel is
py312-abi3 *and* numpy-2-API — document both floors next to each other.

**Validation.** Existing: R2.3's `numpy-floor` CI job. New: none beyond it
(testing against numpy 1.x verifies the *rejection*, not needed).

**Risk/Effort.** S.

**CORRECTIONS (0.6.3).**

* **It was not one floor, it was all four.** `requires-python` is `>=3.12`
  and no declared floor has a cp312 wheel: numpy 1.22 (first cp312: 1.26),
  scipy 1.8 (1.11.2), PyYAML 6.0 (6.0.1), numba 0.56 (0.60.0). Raised to
  `numpy>=2.0`, `scipy>=1.13.0` (first with both numpy-2 support and cp312),
  `pyyaml>=6.0.1`, `numba>=0.61.0`. Fixing numpy alone and *then* making the
  gate blocking would have shipped three known-broken floors behind a green
  check.
* **The CI job no longer names numpy.** It reads every `>=` floor from
  `pyproject.toml`, pins them together, runs `pip check`, and runs the full
  suite on the lowest supported Python. A by-name check is what let the other
  three hide.
* `continue-on-error` deleted, as this item required; job renamed
  `numpy-floor` -> `dependency-floors`.
* Verified locally on Python 3.12.13 before pushing: floors install as wheels,
  `pip check` clean, PEP 517 build is a release build, `pytest validation` =
  631 passed / 2 skipped — identical to the dev environment.

### R4.2 Color gradients on the Python needle path — DONE (0.6.4)

**Review:** §20.2 (Rust fold computes `grad_r`/`grad_t` deposit buckets;
`build_needle_targets`' Python dict omits them; `needle_gradient` has no
grad inputs → the documented Python fold→`needle_gradient` flow yields a
**silent zero** color gradient. Production `run_design` is unaffected).

**Fix design.**
1. Binding (`rust/navette-py/src/synthesis_merit.rs`): the
   `build_needle_targets` result dict gains `"grads_r"` / `"grads_t"` keys
   (flat `na*nw` f64 arrays, same layout as targets/weights).
2. Binding (`rust/navette-py/src/smatrix.rs` needle entry): accept optional
   `grads_r`/`grads_t` (same shape rules as targets) and add the deposit
   terms into the fold accumulation exactly as the native pipeline does
   (`needle_pass.rs:119-120, 436` is the reference implementation to mirror —
   including the documented ÷2 U-curve rule).
3. Python (`src/navette/smatrix/needle.py`): plumb the optional kwargs
   through `NeedleRequest`. Additive only — existing callers unaffected.
4. Docstring of `synthesis/__init__` (the documented flow) updated to show
   the color branch.

**Impact.** Users assembling custom needle cycles in Python with color
demands currently optimize against pointwise targets only, with no error —
a silent wrong-optimizer, not a crash.

**Risk.** M — binding signature grows; keep kwargs optional and ordered;
the native internal path must remain byte-identical (cargo tests pin the
deposit semantics — they must not move).

**Validation.**
- **Existing:** `validation/review/color_merit_check.py` (its Part already
  hand-assembles the color chain rule and matches the native pipeline to
  1.6e-7 — that hand-assembly is the oracle); cargo synthesis tests
  (deposit semantics pinned).
- **New — required:** `validation/review/color_grad_python.py`: drive the
  *Python* path (fold dict → `needle_gradient` with grads) for the same
  3-layer color demand and assert agreement with (a) the hand-assembled
  chain rule ≤ 1e-6 and (b) the native pipeline's gradient ≤ 1e-12. Extend
  `color_merit_check.py` rather than duplicating its setup if convenient —
  new file preferred to keep harnesses single-purpose.

**Effort.** M.

**CORRECTIONS/NOTES (0.6.4).**

* **There is no "fold accumulation" in `needle_gradient` to add into.** Step 2
  says to "add the deposit terms into the fold accumulation exactly as the
  native pipeline does". `needle_pass.rs` accumulates every quantity into ONE
  `acc` per depth; `solver::needle_gradient` instead returns a separate map
  per channel (`P_s`, `P_T_s`, …) and leaves the summing to the caller. The
  deposits therefore go into `P` (from `grads_r`) and `P_T` (from `grads_t`)
  respectively — which reproduces the native total exactly once the caller
  sums the channels it asked for, and is the only placement that keeps the
  existing output contract.
* **The binding change is in the core, not only in `navette-py`.** Step 2
  names `rust/navette-py/src/smatrix.rs`, but the needle entry there is a
  thin forwarder; the arguments had to be added to
  `navette::smatrix::solver::needle_gradient` (free fn *and* the `Solver`
  method) and threaded through **both** PyO3 entries — the `Solver.needle_gradient`
  method, which is what `src/navette/smatrix/needle.py` actually calls, and
  the free `needle_engine` function. The plan mentions only one.
* **A non-zero bucket with no matching request bit is now an error.** The plan
  did not say what happens when a caller passes `grads_r` without `NREQ_P`.
  Dropping it would reproduce this item's own bug — a color demand optimized
  against nothing, silently — so it raises, naming the missing bit. An
  all-zero array passes: the fold hands both arrays through unfiltered and
  callers must not have to filter them.
* **The 1e-12 oracle in the Validation block is not the native pipeline.** No
  native entry point returns a raw `P` profile to Python (`NeedlePipeline.run`
  runs a whole synthesis), so "agreement with the native pipeline's gradient
  ≤ 1e-12" is not directly measurable. Two tighter equalities stand in, and
  both hold at 1e-15: the deposit equals the *pointwise* kernel driven to the
  same scalar (`g·P_ref/(2R)` — the two kernels differ in exactly one factor,
  so this pins the kernel, the channel and the depth row), and the end-to-end
  chain matches `color_merit_check.py`'s hand-assembled rule.
* **`color_merit_check.py` Part C stays hand-assembled on purpose.** It is the
  oracle, so it must not consume the buckets it is meant to check. Only its
  stale "the PYTHON fold dict omits them" comment changed.
* A pytest file was added alongside the required review harness
  (`validation/smoke/test_color_needle_python.py`): the review harnesses are
  not run by CI, and a fix guarded only by something nobody runs is the R4.1
  failure mode again.

### R4.3 Fix and continuously execute `examples/` — DONE (0.6.5)

**Review:** §6.2 (`examples/spectralweave_example.py:22` fails on a fresh
clone — API drift; opaque `Length mismatch` from deep inside native).

**Fix design.**
1. Rewrite the example against the current `OpticalWeaver` API (weave →
   fragment → unweave); run it end-to-end before committing.
2. New `validation/smoke/test_examples.py`: for each `examples/*.py`, spawn
   a subprocess (`sys.executable`, cwd = repo root, env with
   `PYTHONIOENCODING=utf-8`, timeout 120 s), assert returncode 0 and name
   the failing example in the assert message. Subprocess keeps pytest's
   process state and matplotlib/global side effects isolated.
3. While in the wrapper (`OpticalWeaver.unweave`): reject a non-`OpticalFragment`
   template argument with a precise `TypeError` instead of the deep native
   `Length mismatch` (small DX fix, same commit).

**Impact.** Examples stop being dead weight; the API-drift class gets caught.

**Risk.** S — one example fixed, one runner added.

**Validation.** Existing: none existed (§6.2). New: the runner.

**Effort.** S.

**CORRECTIONS/NOTES (0.6.5).**

* **The example's bug was not only API drift.** §6.2 reads as a rename
  problem; the deeper fault is that `unweave` distributes by **exact**
  wavelength match (1e-12) and the example passed a `linspace` over the same
  *range* as the stored frame. Same interval, different grid, one shared point
  out of a hundred. Rewriting against `SimulationWeaver`/`OpticalFragment` was
  necessary but would not by itself have made the script run.
* **The DX fix is in the wrong place as specified.** Step 3 says
  "`OpticalWeaver.unweave`: reject a non-`OpticalFragment` template argument".
  That guard belongs on the wrapper (`SimulationWeaver.unweave` — the class is
  not called `OpticalWeaver`; that name is the *native* type) and it was added
  there. But a bad template never produced `Length mismatch`: it produced
  `AttributeError: 'tuple' object has no attribute '_rust_key'`. The opaque
  `Length mismatch` §6.2 complains about comes from
  `SpectralDataFrame::set_data`, four levels down, and was fixed separately by
  checking coverage in `unweave`/`unweave_collection` where the grid context
  still exists.
* **`unweave_batch` could not be called at all.** Its documented signature is
  `dict[OpticalFragment, np.ndarray]`, but `OpticalFragment` is
  `@dataclass(frozen=True)` over two numpy arrays, so the synthesised
  `__hash__` raises `unhashable type: numpy.ndarray` (and `__eq__` raises the
  ambiguous-truth-value error). `eq=False` fixes both. Not in this item's
  scope as written; found by writing the example's last line, and left broken
  it would have meant shipping an example that quietly avoids a dead API.
* **The runner needs an empty-glob guard.** `examples/` holds exactly one
  file; a parametrization over a glob that matches nothing passes silently,
  so `test_the_examples_directory_is_not_empty` is there to make that visible.
  The runner also asserts the example printed *something* — a script that
  exits 0 having shown the reader nothing has not demonstrated anything.

---

### R4.4 Optimizer backends: hardened built-in LM + optional argmin-ecosystem solvers

**Review:** §3.6 (the LM solves the **normal equations (JᵀJ)** — squaring the condition number; thin-film stacks with correlated layers are exactly where JᵀJ goes singular; the λ-floor bails it out today), §18.2 ("plan docs claim scipy parity; this review never reproduced it").

**Current state** (`rust/navette/src/smatrix/synthesis/thick_opt.rs`, 675 lines): Marquardt damping `(JᵀJ + λ·diag(JᵀJ))δ = −Jᵀr` with the diagonal floored at 1e-14 (`:206-210`), fixed λ factors ×5 / ÷3 (`lambda_up/down`), strict-decrease acceptance (`new_cost < cost`), normal-equation solve via `solve_symmetric` (`:373`, Gaussian elimination), bound veto+clamp (`:238-246`), central-difference rayon Jacobian (h = ∛ε·max(|x|,1)), terminations gtol‖Jᵀr‖∞ / xtol / ftol(actual-only) / MaxIterations / Stalled.

#### R4.4a Reference implementations surveyed (Sept 2025)

| Crate | What it is | Relevant facts |
|---|---|---|
| **argmin 0.11.0** (argmin-rs, MIT/Apache-2.0) | solver *framework*: `Solver<O>` trait, `IterState`, observers/checkpointing, `argmin-math 0.5` linalg backends (nalgebra / ndarray) | Ships GaussNewton (Operator+Jacobian, step width γ, cost-difference tolerance), TrustRegion (dogleg/Cauchy/Steinhaug), quasi-Newton, Nelder-Mead, PSO, SA. **No LM solver anymore** — LM was removed from the argmin core; the ecosystem's reference LM lives in the dedicated crate below. None of its local solvers support box bounds natively. |
| **levenberg-marquardt 0.15.0** (rust-cv, MIT, nalgebra 0.34) | the argmin-ecosystem's reference LM — a MINPACK-`lmdif`-derived implementation | `LeastSquaresProblem` trait (`set_params`/`params`/`residuals`/`jacobian`) + `differentiate_numerically` checker; hyperparameters ftol / xtol / gtol / stepbound; **step solved by column-pivoted QR on J** (`lm.rs:279 PivotedQR::new(jacobian)`) — *not* normal equations; termination on actual **and** predicted relative reduction ≤ ftol (MINPACK semantics); scale-invariant gtol = cos∠(Jeᵢ, r) ≤ gtol; eval-count termination. Unbounded. |

Three consequences drive the design: (1) the condition-squaring critique applies to the current built-in — the reference avoids it via pivoted QR; (2) "add argmin" buys the framework plus extra algorithms, **not** LM — the LM reference is the `levenberg-marquardt` crate; (3) both references are unbounded — Navette's bounds contract must be preserved on the default path and wrapped explicitly on alternative backends.

#### R4.4b Work item A — harden the built-in LM (default backend, zero new dependencies) — DONE (0.6.7)

1. **QR step solve.** Replace `solve_symmetric(&a, &neg(&jtr), &mut delta)` (thick_opt.rs:373) with a QR solve of the damped augmented system `[J; √λ·D]δ ≈ [−r; 0]` (D = the same diag-scaled Marquardt weights, keeping the current 1e-14 flooring semantics). Two scoping facts: (a) the damped augmentation is **always full-rank for λ > 0** (`[J;√λD]ᵀ[J;√λD] = JᵀJ + λD² ≻ 0`), so pivoting is a *stability/diagnostic* feature, not a mathematical necessity — a plain unpivoted Householder QR (~80 lines) already eliminates the condition-squaring problem, while column-pivoted QR (~200 lines, permutation bookkeeping, dropped-column → zero-step-component mapping) additionally detects and degrades gracefully through degenerate columns (pivot-driven column drop) instead of relying on the current λ-escalation retry after solve failure — which remains as fallback either way. Scope decision: start unpivoted, add pivoting only if the rank-diagnostics are actually needed; (b) whatever the flavor, add the **QR-vs-legacy cross-check cargo test**: on well-conditioned problems the new step solve must reproduce the old normal-equations step to high precision (rel ≤ 1e-8 on δ) — the cheapest strong guard against linear-algebra implementation bugs.
2. **Gain-ratio damping update** (MINPACK, replaces fixed ×5/÷3): ρ = actual/predicted reduction; accept when ρ > 0; on acceptance λ ← λ·max(1/3, 1−(2ρ−1)³), ν ← 2; on rejection λ ← λ·ν with ν ← 2ν. Fewer residual evaluations near the optimum, better behavior on correlated layers.
3. **Clipped-step prediction (the one place a naive MINPACK port goes wrong).** The bound veto+clamp modifies the step *after* it is solved, but the MINPACK predicted reduction is computed for the **full** LM step — using the unclipped prediction systematically overestimates it, so ρ is underestimated whenever a bound is active (over-damping, premature ftol exits near boundaries). The predicted reduction MUST be recomputed for the clipped step from the linear model on the surviving components (−gᵀδ − ½δᵀJᵀJδ restricted to the non-vetoed entries); pin it with a dedicated unit test (bound-active case asserting the clipped prediction < unclipped prediction and the accepted step still reduces cost).
4. **MINPACK termination semantics.** ftol: require **both** actual and predicted reductions ≤ ftol (currently actual-only), with "predicted" understood per item 3. gtol: add the scale-invariant angle criterion (cos∠(Jeᵢ, r) ≤ gtol) alongside the existing ‖Jᵀr‖∞ test (config flag `gtol_scale_invariant: bool`, default true for new configs; the scale-dependent check stays available). Note the two criteria are *different notions of stationarity*: the angle test can exit at a different point than the norm test on flat, noise-dominated valleys — expected, but document it.
5. **Eval accounting can drift either way.** Gain-ratio damping typically reduces evals near the optimum, but the ν-doubling rejection path can burn more evals than the fixed ×5 ladder on adversarial steps. Net effect is empirical — the `bench_refold` gate is the arbiter, not intuition.
6. **Keep unchanged:** the bound veto+clamp contract (documented as Navette's bounds semantics — neither reference has bounds), the FD rayon Jacobian, and the `LmTermination` variants (map new conditions onto existing values where possible to avoid API churn; extend only if a new reason is genuinely distinct).

**CORRECTIONS / NOTES (0.6.7).** Item A shipped as written -- QR step solve
(unpivoted, as the scoping decision allowed), gain-ratio damping,
clipped-step prediction, MINPACK ftol, scale-invariant gtol, bounds contract
and FD Jacobian untouched. What the plan did not say:

* **The QR replaces the JᵀJ build outright, not just the solve.** Only
  diag(JᵀJ) is still needed (Marquardt scaling, and the column norms the
  scale-invariant gtol divides by), and that is one m·n pass. Building the
  full n×n product every iteration is gone. `J` is factored once per iteration
  and each λ trial is a 2n×n QR -- MINPACK's `qrsolv` structure. Factoring the
  augmented (m+n)×n system per λ trial, which is the obvious reading of item
  1, would have been O(m·n²) *per trial* and a real regression.
* **The QR's advantage is conditional, and the test says which condition.**
  At the λ the solver starts from, the Marquardt term regularizes JᵀJ enough
  that the two formulations agree to the last digits; the difference appears
  when λ has decayed towards nothing near a good optimum (a Vandermonde J at
  n = 14, λ = 1e-16: the QR step's objective is 0.8 % lower). The claim the
  cross-check pins is therefore "never worse, and better where the damping
  stops covering", not "always better".
* **Item 3's prescription is not quite right, and the right rule is simpler.**
  The plan says to recompute the prediction "on the surviving components",
  which handles the veto but not the *clamp* -- a step from an interior point
  that overshoots a bound is shortened, not zeroed. The prediction is computed
  for `trial − x`, which is both cases at once.
* **A unit test of the prediction was not enough to pin the wiring.** Swapping
  the clipped step for the unclipped one in the driver left all of it green.
  `LmResult` therefore gained `gain_ratio` (a standard LM diagnostic, useful
  in its own right for telling a bound-stalled synthesis from a finished one),
  and a driver-level test asserts ρ = 1 on an exactly-linear clamped problem.
  That test fails with ρ = 0.049 under the swap.
* **Item 4's ftol change could not be pinned end to end either.** No problem
  tried separated the two rules by a run-level observable -- which rule fires
  depends on the damping path, so an end-to-end assertion would have pinned a
  coincidence. The predicate is extracted (`ftol_converged`) and tested
  directly instead.
* **R4.4d's cost-parity criterion is not well posed as written.** Reflectance
  against thickness is oscillatory: from a distant start the engine and scipy
  land in different local minima, and on one of the three cases tried the
  engine's was the *better* one. `lm_check.py` compares them where the basin
  is unambiguous (anchor on scipy's answer, perturb, start both there) and
  asserts the basin-free properties on the far starts.
* **Part of R4.4d is structurally impossible and should stay that way.**
  Running "the bounded problems from the cargo tests through both" needs a
  Python entry point for the LM over an arbitrary residual closure; the FD
  Jacobian is rayon-parallel, so a Python callback would take the GIL inside
  every worker. The harness checks scipy against the *pinned constants* in
  those tests instead -- which is what the item is actually for: the pins were
  written by hand and the solver agrees with them by construction.
* **`LmConfig` gained two fields** (`damping`, `gtol_scale_invariant`), both
  with the new behaviour as the default, and both exposed on `LmConfig` in
  Python. `lambda_up` still drives the error ladder in both damping modes;
  `lambda_down` applies to `Fixed` and to the accept path when the model
  predicted no reduction.
* **No pinned optimum moved.** All 11 pre-existing `thick_opt` cargo tests and
  the whole Python suite passed unchanged, so the §8 golden protocol had
  nothing to triage. `bench_refold`: one LM optimize 2.0 ms against the 2.1 ms
  baseline.

#### R4.4c Work item B — backend selection + optional ecosystem solvers (feature-gated) — DONE (0.6.10)

1. New `synthesis/optimizer.rs`: `enum OptimizerBackend { BuiltinLm, MinpackLm, Trf /* R4.6 */, ArgminGaussNewton, ArgminTrustRegion }`; `OptimizerConfig { backend, ftol, xtol, gtol, max_iterations, max_evals, stepbound, … }` mapping 1:1 onto `LmConfig` for the default; one entry point `run_optimizer(problem, x0, bounds, cfg) -> OptimizerResult` with the existing result shape (x / cost / iterations / evals / termination). Residual system stays the injected closure (`MeritSpec::residuals`) — it is already solver-agnostic.
2. Cargo features, **default OFF** (zero new dependencies for standard builds): `opt-minpack-lm = ["dep:levenberg-marquardt"]`, `opt-argmin = ["dep:argmin", "dep:argmin-math"]` (pin `argmin = "0.11"`, `levenberg-marquardt = "0.15"`). Python surface: `optimizer: str = "builtin"` kwarg on the synthesis entry points; absent feature → `PyValueError` with a rebuild hint (same pattern as the module ImportError fallbacks, §13). `design_config.rs` gains the backend field (`deny_unknown_fields` envelope updated in lockstep — the schema-version rule of §9.3 applies).
3. **Bounds contract per backend:** `MinpackLm`/`Argmin*` are unbounded — wrap with an interior logit reparametrization u = atanh(2(x−lb−ε)/(ub−lb−2ε)), x = lb+ε+(ub−lb−2ε)·(tanh(u)+1)/2. The transform is elementwise and the LM Jacobian is FD per parameter, so gradients pass through it exactly. Document explicitly that boundary semantics on these backends (interior parametrization, optimum strictly inside) **differ** from the built-in veto+clamp (optimum may sit on the bound) — this is a contract difference, not a bug; the built-in stays the default for bounded problems.
4. **argmin integration friction, recorded up front:** argmin-math 0.5 ships backends for nalgebra and ndarray ≤ 0.16; navette pins **ndarray 0.17** (§7) → either write the small manual adapter for our residual types (~100 lines; argmin explicitly supports user-supplied implementations) or convert J to nalgebra matrices per iteration (m×n copy per iteration — negligible against the TMM residual cost). Decide at implementation time; do not bump the ndarray pin for it.

**CORRECTIONS / NOTES (0.6.10) — item B shipped; `MinpackLm` is the one
backend it ships with.** Items 1, 2 and 3 are done: `synthesis/optimizer.rs`
with `OptimizerBackend` / `OptimizerResult` / `run_optimizer`, a default-off
cargo feature per backend, the Python `optimizer=` surface with a rebuild hint
in place of a silent fallback, and the interior reparametrization with its own
tests. What the plan did not say, or said differently:

* **The `OptimizerConfig` / `LmConfig` split is not worth having.** A struct
  "mapping 1:1 onto `LmConfig` for the default" *is* `LmConfig` plus a
  `backend` field — every other knob it names (ftol, xtol, gtol,
  max_iterations, max_evals) already carries over verbatim, because they are
  MINPACK's own names. Two structs that must agree field for field are a
  synchronization bug waiting to be written, so `backend` went on `LmConfig`
  and `run_optimizer` takes it directly.
* **`stepbound` is deliberately not exposed.** It is a MINPACK-only knob, and
  a setting that does nothing on the default backend is a support question,
  not a feature. The crate default stands until a workload argues otherwise.
* **The backend field does not belong in `design_config.rs`.** That module is
  the *design* request — structure, materials, film flags — and carries no
  optimizer settings at all, so §9.3's schema-version rule does not come into
  play. The backend travels on `LmConfig`, which is already threaded through
  `SmatrixContext`, the pipeline and the driver.
* **Item 4's ndarray/nalgebra friction does not arise.** argmin-math 0.5 ships
  a `nalgebra_0_34` backend and `levenberg-marquardt` 0.15 is built on
  nalgebra 0.34, so both crates share one linear-algebra dependency and the
  ndarray 0.17 pin is untouched: no adapter, no per-iteration conversion, no
  bump. Worth recording that the *other* escape the plan hints at does not
  work — argmin-math's default `Vec<f64>` backend has no `ArgminInv`, which
  `GaussNewton` requires, so a Vec-typed argmin problem will not compile.
* **`ArgminGaussNewton` / `ArgminTrustRegion` are not declared yet.** An enum
  variant that nothing can ever select is worse than no variant: it appears on
  the Python surface as a name that is accepted and then fails. They arrive
  with their implementation. `MinpackLm` is declared unconditionally *because*
  it has one — the feature gates the implementation, so a build without it can
  still name the backend and be told how to get it.
* **A missing backend is refused, never substituted.** This is the load-bearing
  half of item 2 and has its own test on both build configurations: a silent
  fall back to the built-in would make every "compared against the reference
  LM" claim a comparison with ourselves.
* **The interior map needs no ε.** The plan's
  `u = atanh(2(x−lb−ε)/(ub−lb−2ε))` is also missing a `− 1` (its argument runs
  over [0, 2], outside `atanh`'s domain). The map shipped is
  `x = mid + half·tanh(u)` with the `atanh` argument clamped to
  `±(1 − 1e-14)`, which is the same guard with one fewer parameter and keeps
  the map exactly symmetric: an `x0` sitting on a bound maps to a large finite
  `u`, and comes back inside the bound by `half·1e-14`.
* **The contract difference is documented in three places** (module docs,
  `LmConfig`'s Python docstring, README) because it is the reason the built-in
  stays the default: the unbounded backends converge to a boundary optimum
  only in the limit, and the gradient vanishes as they approach it — whereas
  the built-in's veto+clamp lets a film land *on* zero, which is how the
  synthesis loop learns to remove it.
* **The MINPACK adapter has three conversions worth naming**, each with a
  test: the crate reports ½‖r‖² where Navette's cost is ‖r‖²; it reports
  evaluations but not iterations, and since MINPACK builds the Jacobian
  exactly once per outer iteration, counting builds *is* the iteration count;
  and it has no error channel at all (`residuals()` returns `Option`), so a
  failure is stashed and re-raised — otherwise "the merit spec is missing a
  curve" arrives as "the solver gave up".
* **The crate does not difference on its own.** Its `differentiate_numerically`
  is a checker, not a fallback, so the adapter calls Navette's own
  `build_jacobian` when there is no analytic source or the source declines —
  which keeps the FD step rule identical across backends.
* **Eight deliberate breaks, all caught**: the interval map made linear, the
  chain rule dropped, the chain rule applied along rows instead of columns,
  the result left in `u`-space, the ½‖r‖² conversion removed, an unavailable
  backend falling back to the built-in, differenced Jacobians counted as
  analytic, and the `atanh` clamp removed.
* **What R4.6 inherits.** `OptimizerBackend::Trf` is a new arm in
  `run_optimizer` and nothing else: the `JacobianSource` seam, the result
  shape, the Python name plumbing and the availability/refusal machinery are
  all in place. TRF has bounds of its own, so it will set
  `bounds_are_native()` and skip `IntervalMap` entirely.
* **`R4.4c-argmin` — DONE (0.6.16); see the corrections block below for what
  the two backends turned out to be worth.** As specified:
  `OptimizerBackend::ArgminGaussNewton` and `ArgminTrustRegion` behind
  `opt-argmin = ["dep:argmin", "dep:argmin-math"]`, with
  `argmin-math`'s `nalgebra_0_34` backend (shared with `MinpackLm`, so no new
  linear-algebra crate). Both are unbounded and take the same `IntervalMap`.
  Note `TrustRegion` is a general-purpose minimizer, not a least-squares one:
  it wants gradient `Jᵀr` and Hessian, so it needs the Gauss-Newton
  approximation `JᵀJ` built explicitly — a different problem shape from the
  other backends, and the reason it is a separate increment rather than a
  second arm in the same one.
* **R4.4d's remaining row.** "Parametrize `lm_check.py` over every enabled
  backend" is not done here: on a default build there is exactly one enabled
  backend, so the harness would parametrize over a single row and assert
  nothing new. It belongs with the first CI job that builds the feature —
  noted under R4.4d rather than silently dropped. **Closed by R4.6 (0.6.12)**:
  `trf` needs no cargo feature, so `lm_check.py` part D is that second row and
  it runs on every build.

**CORRECTIONS / NOTES (0.6.16) — `R4.4c-argmin` shipped, and the answer it
produced is that neither backend should be used.** Both are in, behind
`opt-argmin`, with the plumbing the item specified. What the item did not
specify was a measurement, and the measurement is the result.

* **Measured on the three refold starts `lm_check.py` already uses** (residual
  evaluations to reach scipy's optimum):

  | backend | two films | two films, far | three films |
  |---|---|---|---|
  | `trf` | 8 | 21 | 32 |
  | `argmin_trust_region` | 218 | 224 | 393 |
  | `argmin_gauss_newton` | refuses | refuses | 201, cost 2690.8 vs 2670.0 |

* **`argmin_trust_region` is correct and 7–27× more expensive.** It finds the
  same optimum as scipy every time (1e-13 to 1e-9), but argmin's
  `TrustRegion::terminate` returns `NotTerminated` unconditionally — the solver
  has **no convergence test at all**, so `max_iterations` is its only exit. The
  item did not anticipate this and there is no way to fix it from the adapter;
  it is declared instead, as `OptimizerBackend::runs_to_max_iterations()`, so a
  caller reading `termination` is not misled by a `MaxIterations` that is the
  normal outcome.
* **`argmin_gauss_newton` cannot run this problem class.** The item described
  it as "unbounded and takes the same `IntervalMap`", which is true and beside
  the point: the step is `(JᵀJ)⁻¹Jᵀr` with no damping, so a singular `JᵀJ`
  makes it *undefined*, and a film driven toward zero thickness — the single
  most ordinary thing that happens in a refold — produces exactly that. It
  refuses two of the three starts outright and fails to converge on the third.
  This is the algorithm behaving as the algorithm does; damping is what every
  practical method adds to it. Shipped with a message that says so and names a
  backend that copes, rather than surfacing argmin's bare "Non-invertible
  matrix" to a user who will go looking for a broken stack.
* **The item's premise about `TrustRegion` was right and its consequence
  understated.** "It wants gradient `Jᵀr` and Hessian, so it needs the
  Gauss-Newton approximation `JᵀJ` built explicitly — a different problem
  shape." True; the adapter implements `CostFunction`/`Gradient`/`Hessian` on
  ½‖r‖², `Jᵀr`, `JᵀJ`. What follows from it is the cost: `gradient` and
  `hessian` each rebuild `J`, so one trust-region iteration differences the
  Jacobian twice where a least-squares solver does it once. That is part of the
  7–27×, and it is structural — argmin's trait split gives the adapter nowhere
  to cache it that is not a lie about when the parameter changed.
* **One dependency-plumbing trap, recorded because it fails silently.**
  `argmin-math`'s nalgebra support is the feature `nalgebra_v0_34`.
  `nalgebra_0_34` — the name the item used, and the one cargo generates for the
  optional dependency — enables the *crate* and none of the trait impls, so the
  build fails much later with a wall of unsatisfied `ArgminSub`/`ArgminInv`
  bounds that say nothing about the cause.
* **R4.4d's last row is now closed for real.** "Parametrize `lm_check.py` over
  every enabled backend" was closed by R4.6 on the grounds that `trf` needs no
  feature; part E now does it for the feature-gated ones too, skipping with the
  rebuild command when they are absent. It asserts what the backends *do*, not
  that they are good — the trust region's evaluation count and the
  Gauss-Newton refusal are both pinned, so a future argmin that fixes either is
  noticed rather than assumed.
* **Recommendation, recorded so the next reader does not have to re-measure.**
  `builtin` and `trf` are the backends to use. `minpack_lm` is a useful
  independent LM to check them against. The two argmin backends earn their
  place only as a demonstration of what damping and a stopping test are for,
  and the feature stays off by default.
* **No behaviour change on a default build.** `available_backends()` still
  returns `["builtin", "trf"]`; 687 pytest, 440/445/449/454 cargo tests across
  the four feature combinations, clippy clean on all four, fingerprint
  unchanged.

#### R4.4d Validation

- **Existing:** the 11 `thick_opt` cargo tests (linear exact recovery, exponential fit, bounds corner, mixed interior/bound variables, residual error propagation, max-iterations, x0-outside-bounds, invalid inputs) must stay green or be updated **deliberately** with rationale; `bench_refold` (LM = 2.1 ms baseline — no-regression gate; note the bench's self-reported "no insertions" degeneracy remains item §5.4, unaffected); synthesis regression tests (pinned optima — see risk).
- **New — required:** `validation/review/lm_check.py` — the §18.2 scipy parity, finally executed. Run the bounded problems from the cargo tests **and** a thin-film refold case through both `scipy.optimize.least_squares(method="trf")` and the engine optimizer; assert final-cost agreement (rel ≤ 1e-6), optimum agreement within stated tolerance, and plausible termination reasons. Parametrize the harness over every enabled backend.
- **New — required:** cargo tests comparing backends on identical problems (optimum/cost agreement within tolerance; iteration counts reported, not pinned). The analytic-vs-FD Jacobian question is initially covered by `validation/review/fd_step1.py` / `fd_rchannel.py` for the gradient chain; once R4.5 lands, the per-element cross-check becomes a first-class required cargo test (see R4.5 validation).
- **New — required (work item A internals):** the QR-vs-legacy step cross-check and the clipped-step prediction unit test (R4.4b items 1 and 3) — these two are the guards that make the MINPACK port trustworthy without a new dependency.
- **Golden protocol (§8):** termination-semantics changes may shift pinned synthesis optima (same-or-better cost expected); triage per §8 with one hand-check per file before updating pins.

#### R4.4e Risk, effort, sequencing

- **Risk M:** (a) termination-semantics change can shift regression-pinned optima → §8 protocol; (b) a new backend could regress runtime → bench gate; built-in stays default until an alternative demonstrably wins on real refold workloads; (c) reparametrized backends change boundary behavior → documented as a distinct contract; (d) MIT-licensed crates are LGPL-compatible — note in Cargo.toml comment; (e) gain-ratio damping changes the iteration *path* (not the optimum): iteration counts, ftol/xtol exit points, and everything pinned downstream of them (synthesis regression optima, `bench_refold` timings) can move within tolerance — mandatory §8 triage with one hand-check per changed pin; (f) eval-count variance (#5 above) means bench medians should be compared over a fixed workload before/after, not single runs.
- **Effort:** item A = M (self-contained, precisely specified by MINPACK); item B = M–L (feature plumbing, adapters, Python surface).
- **Sequencing:** A before B (the default must be worthy of being the default). Independent of R6.1 (different file) but both touch synthesis — do not land simultaneously. Place after R4.2/R4.3 in Phase 4.

---

### R4.5 Analytic Jacobian for the refold optimizer — assemble J from the verified deposit chain — DONE (0.6.8 + 0.6.9)

**Review:** §3.6 (`build_jacobian` is central-difference: **2·n full residual evaluations per LM iteration**, thick_opt.rs:305 called at :160), §19.1 (the analytic chain is independently verified: fold evaluator 2e-16, end-to-end gradient 2e-8–7e-7 vs solver FD, dispersion ladder bit-exact), §20.2 (the deposit machinery exists in the Rust fold).

**What.** Each target row's residual depends on one curve value (or an integral aggregate over points), whose per-parameter derivative the needle engine computes **in the same solver sweep that produces the curves** (needle slopes; the dispersion sensitivity chain dφ/dgd/dgdd/dtod/dfod was verified bit-exact in §19.1 Part C). The LM can therefore get exact Jacobian columns for one sweep instead of 2n — and without the FD noise floor that pollutes J exactly where it matters (small steps near the optimum; the solver docstring's own TOD-noise warning).

**Fix design.**
1. Expose **per-point deposits** for the residual assembly: d(curve_c[j])/dθ_k per point/channel, already produced by the sweep — no new solver work, reuse the pass.
2. Expose **fold-chain weights per residual row**: J[i,k] = (d residual_i / d curve value) · (deposit), where the first factor is the transform derivative × rscale (linear/log/phase — verified in `fd_step1` Part A) and, for integral-mean rows, the aggregation weights over points. This is a refactor of the existing fold accumulation (which currently collapses directly into the r-weighted merit gradient) into row-wise form.
3. **Materialize J (m×n)** — sizes are small (m ~ 10³–10⁴ rows, n ≲ 50 layers) — and hand it to the optimizer. R4.4b's QR step solve consumes J directly, so work item A slots in unchanged.
4. `OptimizerConfig.jacobian = "analytic" | "fd"`, analytic the default when deposit coverage exists for every active channel (R/T/A/φ + the dispersion channels all have verified chains); FD remains the fallback for uncovered channels **and** the cross-check oracle.
5. **Determinism:** J assembly is per-(row,point) writes into a fixed-address buffer — the per-destination pattern §13 certifies. Prefer materializing J over any JᵀJ-style accumulation; if an accumulation is ever added, it must keep an ordered reduction (§13's contract).

**Impact.** 2n solver sweeps → 1 per LM iteration (at n = 30 layers, ~98 % of the Jacobian's share of the 2.1 ms refold LM gone); exact J instead of noisy J; benefits every backend (R4.4b/c and R4.6 all consume the same J).

**Risk & mitigation.** M: (a) sub-gradient semantics at fold kinks (kinds a/b/r clamps) must match the FD conventions the native pipeline already uses — same source of truth, but document it at the assembly site; (b) per-element analytic-vs-FD cross-check becomes a **first-class required cargo test**: rel ≤ 1e-6 per element away from kinks on a representative spec; (c) a coverage test asserting every channel activated by a representative spec has deposit support (else the fallback must trip visibly); (d) bench gate.

**Validation.** Existing: `fd_step1` Parts A/B (fold + end-to-end), `bench_refold` (expect a large refold-cost drop; update the bench to assert which jacobian path ran). New: the per-element J cross-check and the coverage test above.

**Effort.** M. **Sequencing:** after R4.4b (feeds it J directly), before R4.6 (TRF consumes the same J).

---

**CORRECTIONS / NOTES (0.6.8).** Split into two increments; **increment i
(fix-design item 2, the merit half) is DONE (0.6.8)**. Increment ii — item 1
(per-point deposits), item 3 (J assembly), item 4 (`jacobian =
"analytic" | "fd"`) and the bench gate — is open.

* **The split is the point, not an accident of scheduling.** Item 2 is the
  only half that can be verified without a solver in the loop: differencing
  `residuals()` in *curve* space gives an exact oracle for ∂r/∂(curve value),
  and `sensitivity_is_a_finite_difference_of_residuals` uses it over five
  kinds × two transforms. When the deposits arrive, a disagreement with the
  FD Jacobian has one place left to be.
* **Item 2 says "a refactor of the existing fold accumulation into row-wise
  form"; that is not what was done, on purpose.** The row-wise pass walks the
  target grids a *second* time rather than threading a sink through
  `residuals_into`. That function is bit-exactness-critical — pinned point for
  point against the Python original — and threading a sink through it to save
  one traversal would put every residual in the engine at risk to buy a
  Jacobian. The duplication is the safer trade; the FD cross-check is what
  keeps the two walks from drifting, and in debug builds
  `curve_sensitivity` additionally asserts its row count against
  `residuals()` itself.
* **"Each target row's residual depends on one curve value" is not true of two
  cases the engine supports.** Absorption rows read **two** curves
  (A = 1 − R − T, both companions entering with −1), and interpolated target
  grids read **two adjacent** simulated points with weights (1−f, f). A row is
  therefore a *list* of terms, not one term. A one-term design would have
  looked right on every aligned intensity spec — which is most of them.
* **Phase and color rows are reported, not silently zeroed.** Neither chain is
  carried here (phase because `arg()` plus the differential reference depends
  on total thickness directly; color because `build_needle_targets` already
  emits that derivative per solver point). They occupy their rows and are
  listed in `uncovered`, so item 4's "FD remains the fallback for uncovered
  channels" has something concrete to test: `is_complete()`.
* **`n_residuals()` was wrong for integral frames** — it counted one component
  per point where `residuals()` pushes one per *frame*. Found by the row-count
  assertion this pass needs; fixed, with the residual grid-overlap caveat now
  documented (a frame that misses the simulated grid contributes nothing, and
  that cannot be known from the spec alone).
* **Sub-gradient semantics at the kinks (risk (a)) are settled and tested.**
  `a`/`b` are exactly flat on their satisfied side, `r` exactly flat inside
  the band, `Log` exactly flat below its 1e-12 clamp — in every case the
  derivative of the arm `kind_residual` itself takes, so the two agree on
  which side of the boundary a point is on. An inactive constraint yields no
  terms at all rather than zero-valued terms, so J assembly does no work for
  it.
* **Seven deliberate breaks, all caught** (dropped 1/n on the integral mean,
  `Log` derivative replaced by the norm factor, absorption companion sign
  flipped, sample reconstruction collapsed to the left grid point,
  interpolation weights swapped, dead-zone derivative made live, grid-miss
  frame emitting empty rows). The first draft of the interpolation test did
  *not* catch two of them: under `Linear`/`Exact` the derivative is the same
  number wherever the row is evaluated, so that pair pins the interpolation
  weights and nothing else. The test now also runs `Log`, which pins the point
  the derivative is taken *at*.

**CORRECTIONS / NOTES (0.6.9) — increment ii, and R4.5 is DONE.** Items 1, 3,
4 and 5 shipped: per-point deposits, a materialized m×n J, `jacobian =
"analytic" | "fd"` with analytic the default, and the determinism contract.
Both required tests exist (the per-element cross-check and the coverage test),
and the bench gate reports which path ran.

* **The thickness derivative did not need new solver work — it is the needle
  operator with the needle material set to the host's own index.** Growing
  film *j* by δ is inserting a slab of *n_j* inside layer *j*: r₁₂ = 0, so ρ̂
  vanishes and τ̂ = iβ_j is the bare propagation slope. `needle_slopes4_ddz`
  then differentiates the whole block through the same dual-number
  composition, which is why the first run of the cross-check agreed with
  finite differences to ~1e-9 without a single correction. The plan's item 1
  says "already produced by the sweep — no new solver work, reuse the pass";
  that is right about the *machinery* and wrong about the cost (below).
* **The impact estimate is too optimistic and should be restated.** "2n solver
  sweeps → 1 per LM iteration" is not what happens: the deposits need a
  `StackFields` decomposition per point *and* polarization, so an analytic
  iteration is about **three** sweeps against **2n + 1** for the differenced
  one. Measured at ten films, iteration for iteration: **1.7× wall-clock and
  16× fewer residual evaluations**. The evaluation ratio keeps growing with n;
  the wall-clock ratio approaches ~(2n+1)/3. The claim that does hold in full
  is the second one — the columns are exact, with no difference noise near the
  optimum.
* **A residual row reads a LIST of curve values, which is what made the split
  worth it** (see the 0.6.8 note). The assembly therefore accumulates per row
  and the absorption case (two channels) is its own test.
* **Item 4's "analytic the default when deposit coverage exists" is decided
  per run, not per build.** `JacobianSource::fill` answers `Ok(None)` for a
  spec with uncovered rows and the driver differences that iteration —
  visibly, via `LmResult::analytic_jacobians`. Declining and *failing* are
  separate: a source that errors aborts the run, because a Jacobian that
  cannot be built is a different fact from one that does not apply.
* **Which path ran is now observable rather than inferred.**
  `LmResult::analytic_jacobians` and the new
  `SmatrixContext::optimize_thicknesses_report` (bound as
  `optimize_thicknesses_report`) exist because the plan's "update the bench to
  assert which jacobian path ran" is not possible otherwise — a timing cannot
  tell a fast analytic run from a silent fallback. The report also carries
  `termination` and `gain_ratio`, which is what the correction below needed.
* **No pinned optimum moved.** The whole cargo suite and the whole Python
  suite pass unchanged with analytic as the default, so risk (a)'s §8 triage
  had nothing to triage. The bit-exactness fingerprint is unchanged:
  `simulate_with_deposits` shares the sweep with `simulate` and returns the
  same bits, which is itself a test.
* **A correction to `validation/review/lm_check.py`, found by this change.**
  Part A2's "the engine's answer is a fixed point of itself" was, on the
  far-start cases, testing the iteration cap: at the default 200 iterations
  both solvers are still crawling, and a run that stopped on
  `MAX_ITERATIONS_REACHED` has not claimed a stationary point. The analytic
  path exposed it by taking a different route through the same flat valley
  (a few parts in 10⁵ behind at iteration 200). The harness now gives those
  cases room and checks the termination criterion explicitly; with room, both
  paths converge to the same cost and both are fixed points.
* **Eight deliberate breaks, all caught**: the T deposit losing its flux
  factor, differentiating the forward amplitude instead of the backward one,
  the film index not offset past the ambient, the factor of two in d|z|²/dd,
  the needle taking the ambient's index instead of the host's, a transposed J,
  a row keeping only its last term, and the coverage gate removed.
* **What R4.6 inherits.** TRF consumes the same `J` — `JacobianSource` is the
  seam, and a TRF backend gets the analytic Jacobian for free by taking the
  same argument.

### R4.6 TRF backend — trust-region-reflective, the bounded-LS reference method — DONE (0.6.12)

**Review/context.** `thick_opt.rs`'s own docstring says it replaces scipy `least_squares(method="trf")` — the Rust rewrite dropped to a simpler clamped LM, and the clamp is exactly the weak point R4.4b item 3 patches. TRF (Branch–Coleman–Li 1999) is the reference bounded-LS method: bounds enter the trust-region subproblem (scaled space + reflection steps), iterates stay strictly interior, and boundary optima are handled by construction — the existing corner/boundary cargo cases are precisely where clamp-LM is weakest.

**Fix design.**
1. New `OptimizerBackend::Trf` in the R4.4c enum. Hand-rolled like the rest of the module — scipy's `_lsq/trf.py` is the semantic reference, no new dependency.
2. Structure: D = column-norm scaling (frozen, MINPACK-style); **generalized Cauchy point** in the scaled box for the initial step; 2-D subspace minimization with reflections at active bounds; ρ-based radius update (same gain-ratio machinery as R4.4b item 2); the subproblem's least-squares solves **reuse R4.4b's QR infrastructure** — this is why item A is not wasted work when TRF arrives.
3. Termination: the same ftol/xtol/gtol MINPACK semantics — which is why scipy is a *direct oracle* here (same algorithm family): extend `lm_check.py` with `method="trf"` rows and assert the tightest agreement in the whole plan (same algorithm vs same algorithm).
4. **Bounds contract:** replaces veto+clamp entirely on this backend — the R4.4b item-3 clipped-prediction caveat *does not exist* for TRF (the predicted reduction is computed for the actual feasible step). When TRF is made the default for bounded problems, that caveat retires rather than being patched forever.
5. Requires J — analytic (R4.5) preferred, FD fallback; same problem interface as the other backends.

**Impact.** Correct boundary behavior, direct scipy parity, kills the clamp-prediction subtlety by construction.

**Risk & mitigation.** The largest algorithmic lift in the R4.4 family (Cauchy point + reflection logic are fiddly): mitigated by (a) scipy-as-oracle in `lm_check.py`, (b) reuse of QR + gain-ratio pieces, (c) backend-equivalence tests, (d) bench gate; shipped behind the backend enum — built-in stays default until TRF demonstrably wins on real refold workloads.

**Effort.** L. **Sequencing:** after R4.4b (QR + gain-ratio) and R4.5 (analytic J).

**CORRECTIONS / NOTES (0.6.12) — shipped, and the plan's algorithm sketch was
not TRF.** `synthesis/trf.rs` implements `scipy.optimize._lsq.trf.trf_bounds`
with `tr_solver="exact"`, reachable as `LmConfig(optimizer="trf")`. Every
piece of the fix design landed except the two in item 2 that described a
different method:

* **The scaling matrix is Coleman-Li, not "column-norm, frozen,
  MINPACK-style".** Those are two different matrices doing two different
  jobs, and the plan conflated them. `D = diag(√v)` where `v` is the
  *distance to the bound the anti-gradient points at* — recomputed every
  iteration, because the whole method is that the trust region changes shape
  as the iterate approaches a bound. MINPACK's column-norm diagonal is
  scipy's `x_scale`, an orthogonal user-facing knob whose default is 1, and
  it is fixed at 1 here because `LmConfig` has no field for it. Freezing a
  column-norm `D` and calling it TRF would have produced a clipped
  trust-region method with none of the boundary behaviour the item is for —
  and it would have cost the scipy oracle, which is the item's own
  validation plan.
* **"2-D subspace minimization with reflections" is two unrelated things.**
  The 2-D subspace (indefinite dogleg) is scipy's `tr_solver="lsmr"`, for
  sparse Jacobians of millions of rows; the dense default is
  `tr_solver="exact"`, the Moré subproblem the item's own point 2 already
  asks for. Reflections are not part of either solver — they are in
  `select_step`, which scores three candidates (the trust-region step cut
  back to the first bound it hits, that step **reflected** off the bound, and
  the constrained Cauchy step) and takes the best. Thickness Jacobians are
  dense with one column per film, so `"exact"` is both the right choice and
  the one the oracle runs.
* **The QR reuse the item hoped for is real, and replaces the SVD.** scipy
  solves the subproblem from an SVD of the augmented Jacobian. Moré's
  algorithm actually needs only `p(α)` and the `‖q‖` with `Rαᵀq = p`, and the
  augmented QR from R4.4b supplies both from one factorization — the Newton
  recurrence on the secular equation is then identical term for term. The
  single place the two part company is the rank test: scipy compares σ_min to
  σ_max, this compares `min|Rᵢᵢ|` to `max|Rᵢᵢ|`, which brackets it. A
  borderline call changes which branch seeds `α`, not the answer. So R4.4b
  was not wasted work, for the reason the plan gave.
* **Point 4 holds, and it is visible.** The clipped-prediction caveat does
  not exist here: the predicted reduction is computed for the feasible step
  that was actually taken. `lm_check.py` D2 runs a merit whose optimum is
  *on* the clamp; both bounded backends find it and agree to 1e-11 nm.
* **There is one contract difference, and it goes the other way.** TRF's
  iterates must stay *strictly* interior — that is what keeps `v`
  differentiable — so it returns a thickness one ULP short of the bound where
  the built-in returns the bound itself (49.999999999986834 against 50.0, a
  13-femtometre gap at an identical merit). The removal sweep compares
  against `clamp_min`, so a film still gets removed; a caller testing
  `x == ub` would be surprised. Documented in the module, the Python
  docstring and the README, and asserted in both test layers. It is also why
  making TRF the *default* is not automatic, and is not done here.
* **`x_scale` and `stepbound` stay unexposed** for the same reason
  `stepbound` did in R4.4c: a knob that does nothing on the default backend
  is a support question, not a feature.
* **Ten deliberate breaks, all ten caught** — Coleman-Li removed, the `dv`
  sign flipped so `C` goes negative, the reflection candidate disabled, the
  `theta` step-back removed, `φ'` sign-flipped in the secular equation, the
  step left in hat space, the half-cost conversion dropped, the radius never
  shrinking, the gain-ratio guard dropped from the ftol test, and the
  curvature diagonal `C` dropped from the augmented system. Worth recording
  *which* layer caught them: eight fell to the cargo tests, but the
  **reflection candidate and the `theta` step-back were caught only by
  `lm_check.py`** — and those two are exactly what separates TRF from a
  trust-region method that clips. Without the scipy oracle this item would
  have shipped with its defining feature untested.
* **R4.4d's remaining row is now done.** "Parametrize `lm_check.py` over every
  enabled backend" was deferred in 0.6.10 because a default build had one
  backend and the harness would have asserted nothing new. `trf` needs no
  feature, so part D is that parametrization and it runs everywhere.
* **The comparison is tolerance-level, not exact, and the reason is
  recorded.** Costs agree to 1e-9 relative; thicknesses to 1e-4 nm. The two
  solvers difference the merit differently — Navette hands TRF the analytic
  Jacobian (R4.5) and falls back to central differences, scipy's default is a
  forward difference — so the iterates diverge slightly inside the same ftol
  basin. 1e-4 nm is a tenth of a picometre.
* **What this does not change.** `builtin` remains the default and no pinned
  optimum moved: the whole suite passes unchanged, and the bit-exactness
  fingerprint is untouched (TRF is opt-in and no existing call site selects
  it).

---

## 5. Phase 5 — Performance (P2/P3)

### R5.1 `unweave_collection` batch path — DONE (0.6.13)

**Review:** §5.3/§5.2 (only genuine perf regression on the release build:
0.22× @ 100k pts × 512 keys; 0.48× @ 500k × 16 keys; diagnosis: per-key
`AHashMap` rebuild + `String` clones + `frames.read().clone()`).

**Where.** `rust/navette-py/src/spectralweave_optical.rs:435+` (pointer-capture
entry) and the weaver engine's collection path it calls.

**Fix design.**
1. Hoist the key map: the entry already collects `items: Vec<((f64,String,String), ptr, len)>`
   once — build the `AHashMap<OpticalKey, …>` **once per call**, keyed by an
   owned-but-cheap key (intern the `(String, String)` pair into an `Arc<(str,str)>`
   or use indices into a `Vec` of extracted keys instead of cloning per key).
2. Frame access: replace `frames.read().clone()` (clones the whole collection
   per operation) with a single `frames.read()` guard held across the
   per-key work, cloning only the *needed* frame `Arc`s (or borrowing them
   under the guard — the guard's lifetime spans the detached closure).
3. Results materialization: write per-key outputs into one preallocated
   arena keyed by index rather than per-key `Vec` churn.
4. Keep the pointer-capture pattern intact (§13 flags it fragile but
   correct); while editing, add the missing "lifetime anchor" comment on
   `py_arrays` (§13 — fold into this item).

**Impact.** Brings the last regressing data-plane op to ≥ 1× vs numba at
512-key scale; removes per-call allocation churn for all callers.

**Risk.** M — concurrency-sensitive file (§22 verified the locking holds
*with the current structure*; re-run the race harness after). Borrow-checker
friction around the read guard inside `py.detach` is the likely time sink.

**Validation.**
- **Existing:** `validation/review/weaver_race.py` (must stay: 0 exceptions,
  0 torn writes, all keys bit-exact); spectralweave parity tests;
  `navette_spectral_bench --quick` (the measurement itself, before/after
  JSON with `build_profile` stamp per R2.1).
- **New — optional:** extend the spectral bench's parity section to assert
  `unweave_collection` output == per-key `unweave` output on the same data
  (cheap, kills the whole class of batch-vs-single drift).

**Effort.** L.

**CORRECTIONS / NOTES (0.6.13) — shipped, and the diagnosis in this item was
wrong twice over.** The batch path is 1.15×–2.29× faster where the cost is
per-fragment work, unchanged where it is bandwidth, and it no longer panics on a
short curve. The Impact line is not reachable and the reason is a semantic
difference, not a missing optimization.

* **The named causes were not the cost.** "Per-key `AHashMap` rebuild +
  `String` clones + `frames.read().clone()`": the binding already builds the
  map once per call (fix-design point 1 was done before this item was picked
  up), `OpticalKey` holds `Arc<str>` so a clone is two atomic increments, the
  plan is resolved once, and `unweave_collection` never called
  `frames_snapshot` at all. Measured, the per-key cost is memcpy bandwidth:
  10.4, 7.1 and 15.8 GB/s across three configurations before the change, which
  is what a copy of that size costs on this machine and nothing like what a
  hash-map rebuild costs.
* **The Python reference wins because it does not own its data.** `full_data[s:e]`
  on a numpy array is a *view*; the Rust engine copies, because a fragment
  outlives the borrowed buffer the call was handed. A probe settles it — mutate
  the caller's array after `unweave_collection` returns and the stored curve
  changes on the Python engine (first value becomes 99.0) and does not on the
  Rust one (0.0). So the **Impact line — "brings the last regressing data-plane
  op to >= 1x vs numba at 512-key scale" — is not reachable.** 512 keys x 100k
  points is 400 MB copied in and 400 MB out; the 36 ms that takes is 22 GB/s of
  DRAM traffic, this machine's ceiling. Matching the reference means adopting
  its aliasing, which would make a stored fragment change under a caller who
  edited their own array afterwards. That is a correctness hazard, not an
  optimization, and it is not done.
* **Fix-design point 2 does not apply and point 3 was already true.** There is
  no `frames.read().clone()` on this path to replace with a held guard; the
  plan's `Arc<SpectralDataFrame>` handles are what the loop walks, and they are
  cloned once when the plan is built and cached. Results are not materialized
  per key either — fragments are written straight into the frames.
* **What was actually worth doing**, all of it found by measuring rather than
  from the item:
  1. **Copy only the span the plan reaches.** Frames rarely tile the whole
     curve handed to them, and the rest of it was being copied once per key.
  2. **One collection-wide lock per key, not per fragment.** `map_frame_to_key`
     took the `key_map` write lock once per (key, frame) pair — 16 384
     acquisitions of one lock for a 512-key call. `map_frames_to_key` takes it
     once per key; the singular form delegates to it, so there is one
     implementation.
  3. **`set_data` hashes the key once**, not twice.
  4. **Rayon across keys above 8192 fragments.** Gated on the fragment count,
     not the byte volume, because the bytes are the part that does not
     parallelise: the same 512-key call measures 39.7 / 35.1 / 36.5 ms on
     1 / 4 / 32 threads. What scales is the per-fragment work, and across the
     bench grid the split is clean — every configuration at 16 384 fragments
     gains (1.06× to 2.11×) and every one at 4096 or fewer loses, by up to 2.2×
     where the hand-off costs more than the whole call.
* **Two correctness defects, neither in the item.** A curve that is not one
  value per wavelength was indexed straight out of its buffer: a Rust panic,
  surfacing as `pyo3_runtime.PanicException`, which does not derive from
  `Exception` — so `except Exception:` around `unweave_batch` did not catch it.
  The Python wrapper checked this for a single `unweave` and never for a batch.
  And the coverage check ran *inside* the per-key loop, so a bad grid rejected
  the call after an arbitrary prefix of the batch had been written — which,
  once the loop can run on a pool, would have become a scheduling-dependent
  prefix. Both checks are now hoisted ahead of the first write and the batch is
  all-or-nothing.
* **Determinism (§13).** Each key is handled by one thread start to finish, so
  the frame order recorded under a key is still the plan's — which is what
  `get_converted` hands back. The keys arrive in an `AHashMap`'s iteration
  order, which was already not the caller's, so the loop never ordered anything
  observable. Asserted directly (`every_key_records_its_frames_in_plan_order`)
  rather than argued. `weaver_race.py` re-run: 0 exceptions, 0 torn writes, all
  192 keys bit-exact.
* **Point 4 of the fix design (the `py_arrays` lifetime-anchor comment) is not
  folded in here.** The binding is untouched by this item; the comment belongs
  with §13's review of the pointer-capture pattern and is left with R6.x.
* **The optional validation item is done**, and on both engines:
  `batch_equals_per_key` in the spectral bench asserts a 40-frame × 256-key
  batch writes exactly what 256 single `unweave` calls write. Sized to cross
  the parallel threshold, so the pooled path is the one under test.
  `opticalweaver.rs` also gains the 13 tests it had none of — it was covered
  only from Python.

### R5.2 Parallelize the serial derive loop — DONE (0.6.15)

**Review:** §5.4 (`Solver::solve`'s per-point derive pass — atan/sqrt/arg —
is serial; the next Amdahl bottleneck at 10⁵+ points).

**Fix design.** Rayon with **per-thread output buffers indexed by
destination** (the pattern §13 certifies as bit-deterministic; no parallel
*reduction*). Scope the first cut to the pure per-element trig chain; leave
anything with cross-point accumulation serial.

**Risk.** M — bit-identity is contractual: the random-stack differential
test must stay **bit-identical**; any last-ULP drift means the parallel
structure is wrong (a reduction crept in) — treat as failure, not tolerance.

**Validation.**
- **Existing:** bit-identity test; `bench_backside_speed` before/after
  (target: measurable win at 10⁵+ points, zero change at 120 points).
- **New — optional:** a scaling sweep extension of `bench_backside_speed`
  (40 → 10⁶ points, report µs/point vs size) — also fills §5.4's "no solver
  scaling bench" gap.

**Effort.** M–L. **Priority note:** P3 — do it after R5.1 and only if the
scaling bench shows the serial fraction matters at real workloads.

**CORRECTIONS / NOTES (0.6.15) — shipped; the item was right, and its gate was
the useful part.** The derive pass now splits across the rayon pool. 1.42× at
20 000 points and 1.73× at 60 000 on a twelve-channel request; ~1.03–1.14× on
a four-channel one.

* **The priority note's gate was met by R5.3's bench, not by argument.** "Only
  if the scaling bench shows the serial fraction matters at real workloads" —
  it does, but only after R5.3: with `Solver::new` still costing 2.2–2.5 ms of a
  4.3–4.8 ms call, the derive loop was ~23% of the call and this item would
  have been worth ~1.1×. Once construction was fixed the same 1.1 ms became
  ~46% of a 2.4 ms call. **The sequencing in the plan (R5.3 before R5.2) was
  right for a reason the plan did not know**: R5.3 is what made R5.2 worth
  doing.
* **The fix design is taken as written.** Per-destination indexed writes, no
  reduction. The form is a recursive halving of the point range *and every live
  buffer together* (`rayon::join` down to a 512-point leaf) rather than
  per-thread scratch buffers — same guarantee, no merge step: each output index
  is written once, by one task, with the expression the serial pass used.
  Bit-identity is therefore structural and holds at any leaf size, which is
  asserted directly (`the_split_derive_is_bit_identical_to_the_serial_one`)
  rather than argued.
* **"Scope the first cut to the pure per-element trig chain" was not needed.**
  Nothing in the pass crosses points: every channel is a function of one
  `OpticalState`. The cross-point work is the dispersion pass, which is a
  separate stage and untouched. So the whole loop parallelises, not a subset.
* **Two alternatives measured and rejected.** One parallel pass *per channel*
  (`buf.par_iter_mut().zip(states)`) is simpler but walks the 3.8 MB state
  array once per requested channel — 12 passes for a rigorous request, which is
  more traffic than the single serial pass it replaces. Emitting inside the
  existing solve map (R5.3's fix sketch) would make that map's write pattern
  channel-dependent, which is the §13 property worth keeping. The halving
  keeps one pass over the states and one writer per index.
* **The gain is a function of the request, not the grid.** The pass costs what
  the caller asked for: four trivial channels (photometry) have almost nothing
  to divide and gain ~1.03–1.14×, twelve channels with `sqrt`/`atan`/`atan2`/
  `arg` gain 1.42× at 20 000 points and 1.73× at 60 000. A benchmark that only
  ever measured a photometry request would have concluded this item was not
  worth doing.
* **The optional validation item is not done and is now redundant.** The plan
  offered "a scaling sweep extension of `bench_backside_speed` (40 → 10⁶
  points) — also fills §5.4's 'no solver scaling bench' gap." R5.3 filled that
  gap with `bench_core_engine_scaling.py`, which sweeps the grid at two request
  breadths — and the request breadth is the axis this item actually moves, which
  a point-count sweep alone would have missed.
* **The leaf size is a scheduling knob, not a numerical one.** Swept at 256 /
  512 / 2048 / 8192: within ~3% of each other at the sizes that matter, the one
  real effect being that a leaf above the grid means no split (2 000 points:
  0.428 ms at 2048, 0.386 ms at 512). 512 chosen as the smallest leaf that
  still buys something measurable.
* **Bit-identity (§0.1 rule 5) held.** Fingerprint unchanged
  (`a99e8383…f154`), ten review harnesses exit 0, both parity oracles pass,
  687 pytest / 476 cargo tests green, clippy clean with and without
  `opt-minpack-lm`.

---

### R5.3 `core_engine`'s emit path — DONE (0.6.14) — the rewrite is behind the kernel it replaced

**Discovered by R2.4a (0.6.11).** With the whole-engine parity oracles finally
running, the two engines can be timed on the same inputs, and the answer is
not the one the rewrite assumed: `navette._smatrix.core_engine` is **0.5–0.65×
the speed** of the numba `loom_matrix` kernel it replaced, at every grid size
measured (500 to 60 000 points, 6 layers, release build, both parallel —
numba `prange`, Rust rayon).

**Where the time goes — one measurement, not a diagnosis.** At 20 000 points:
an `Rs`-only solve (`Level::Intensities`, one emitted array) takes 0.98 ms; the
full twelve-channel request takes 2.08 ms; numba computes all thirteen in
1.94 ms. So roughly half the Rust time is spent after the physics. The obvious
suspect is the shape of the pipeline — `solve` materializes a
`Vec<OpticalState>` over every grid point, then walks it once per requested
channel into its own `Vec<f64>`, then again into a numpy array — where numba
writes each channel directly into its output array inside the point loop. That
is a hypothesis; profile before rewriting.

**Fix sketch.** Emit into the destination buffers inside the parallel map
(per-destination indexed writes, §13 determinism — no reduction), so the
intermediate `Vec<OpticalState>` and the per-channel walk both disappear.
`PyArray::from_vec` already moves rather than copies, so the last hop is
probably not the problem.

**Risk.** M — this is the hot path of the whole package and every golden in
the suite runs through it; it must stay **bit-identical** (§0.1 rule 5: the
random-stack differential test plus the review harnesses), and the parallel
write pattern must stay per-destination-indexed.

**Validation.** Existing: the two R2.4a parity scripts are now the oracle and
print the speedup ratio; the differential bit-identity test; the smatrix
benches. New: a grid-scaling row in the smatrix bench so the ratio is tracked
rather than rediscovered.

**Effort.** M–L. **Priority note:** P2 — this is the package's advertised
reason for existing, so it outranks R5.2; sequence it after R5.1.

**CORRECTIONS / NOTES (0.6.14) — shipped, and the item pointed at the wrong
half of the call.** `core_engine` is 1.4×–3.0× faster across the grid. The item
named the emit path; the emit path was never the cost. Two of the numbers this
item was built on were also measurement artefacts.

* **The item's own hypothesis, quoted so the correction is unambiguous:** "the
  obvious suspect is the shape of the pipeline — `solve` materializes a
  `Vec<OpticalState>` over every grid point, then walks it once per requested
  channel into its own `Vec<f64>`, then again into a numpy array … That is a
  hypothesis; profile before rewriting." Profiled, at 20 000 points, 6 layers,
  before the change: `Solver::new` 2.2–2.5 ms, binding transpose 0.14–0.45 ms,
  the parallel `solve` map 0.9–1.1 ms, the derive loop 0.3 ms (photometry) to
  1.1 ms (rigorous), dispersion + pack ~0.02 ms, `solution_to_dict` 0.01–0.25 ms
  — 4.3–4.8 ms total. **Construction was ~55% of the call and the emit was
  under 5%.** The pipeline also does not walk the states once per channel; it
  walks them once, filling every requested channel in the same pass.
* **What it actually was: 40 000 heap allocations per call.** `n_cache` and
  `inv_n_cache` were `Vec<Vec<Complex64>>` — one inner `Vec` per wavelength,
  twice — and 1.77–2.24 ms of the 2.2–2.5 ms constructor was building them.
  They are now two flat wav-major `Vec<Complex64>`s of stride `n_layers`, read
  through `layer_n(w)` / `layer_inv_n(w)`. Same values, same order, bit-identical
  output.
* **The cache was transposed twice to arrive where it started.** The binding
  holds a wav-major interleaved `[re, im]` buffer; `Solver::new` wants
  layer-major complex and immediately re-transposes to wav-major internally. So
  `core_engine` built a whole intermediate `Vec` for nothing.
  `Solver::from_wav_major_flat` takes the caller's layout directly; `new` keeps
  its signature for the callers that genuinely have layer-major data. A
  side-effect worth having: the length is now a checked precondition with a
  message, where the transpose loop was a raw index and a mis-sized cache was an
  out-of-bounds panic — which through PyO3 is a `PanicException`, not a
  `ValueError`.
* **`flat_cache` was eager and almost nobody read it.** Only `needle_gradient`
  does; every ordinary solve was building it. Now a `OnceLock`.
* **The emit was still worth one line.** `solution_to_dict` borrowed its
  `Solution` and cloned every channel into `PyArray::from_vec`, which moves — so
  the whole result was copied once on the way out. Taking the `Solution` by
  value removes the copy. Worth 0.01–0.25 ms, i.e. the item's entire named
  target.
* **The 0.5–0.65× ratio was measured through numba's spinning thread pool.**
  numba's default threading layer leaves its workers spin-waiting after a call
  returns, so timing the two engines back to back charges whichever runs second.
  The same Rust call at 20 000 points: **1.52 ms on its own, 2.00 ms immediately
  after a block of numba calls, 1.58 ms after a one-second pause.** R2.4a's
  speed sections timed numba then Rust, with nothing in between. Both parity
  scripts and the new bench now pause between engines. `NUMBA_THREADING_LAYER=
  workqueue` removes the effect too, but makes numba itself ~70% slower here
  (1.004 → 1.721 ms), which is not a fair comparison either — a pause leaves
  both engines in their native configuration.
* **The fix sketch is not taken, and should not be.** Emitting into destination
  buffers inside the parallel map would remove ~0.02 ms of packing and the
  `Vec<OpticalState>` it walks, at the cost of making the parallel map's write
  pattern channel-dependent — real §13 risk for a target that the profile says
  is ~0.5% of the call. The `Vec<OpticalState>` stays.
* **One experiment failed and is recorded so it is not retried.** Small grids
  are dispatch-bound (~150 µs of rayon start-up is most of a 500-point call),
  so a `SOLVE_PAR_THRESHOLD = 2048` was tried to run them serially. It made
  2 000 points **5× slower** (0.35 → 1.75 ms) and 500 points 1.9× slower
  (0.192 → 0.358 ms): serial costs ~875 ns/point against ~80 ns/point effective
  under rayon, so the pool wins even at 500 points despite its start-up. Not
  taken.
* **Where the ratio stands now**, from the new bench (photometry / rigorous,
  against the numba kernels, noisy because numba's own timings swing ~2× run to
  run): ~0.3–0.5× at 500 points, ~0.6–0.8× at 5 000, ~0.8–1.2× at 20 000 and
  60 000. The rewrite now passes the kernel it replaced on large grids. It is
  still behind on small ones, and that is dispatch cost, not physics.
* **What is left, quantified.** The serial derive loop is now the biggest
  remaining block — ~1.1 ms of a 2.4 ms rigorous 20 000-point call. **That is
  R5.2, and its stated gate ("only if the scaling bench shows the serial
  fraction matters at real workloads") is now met**; the master row is annotated.
  Two smaller ones, measured but not done: folding `inv_n_cache` into
  `solve_point` (~0.3 ms, removes a second pass over the cache) and letting
  `n_cache` alias the caller's buffer via a `bytemuck` cast (~0.17 ms, needs a
  lifetime parameter on `Solver` and would tie the solver to the caller's array
  — the same ownership trade R5.1 declined).
* **The optional validation item is done.**
  `validation/benches/smatrix/bench_core_engine_scaling.py` sweeps 500 → 60 000
  points at two request breadths, so the ratio is tracked rather than
  rediscovered — which is the whole reason this item existed: the single 500-
  point measurement in the parity scripts could not show a cost that scaled with
  the grid.
* **Bit-identity (§0.1 rule 5) held.** The random-stack differential fingerprint
  is unchanged (`a99e8383…f154`), all ten review harnesses exit 0, both parity
  oracles pass, 687 pytest / 472 cargo tests green, clippy clean with and
  without `opt-minpack-lm`.

---

## 6. Phase 6 — Maintainability, docs, hygiene (P3)

### R6.1 `needle_gradient` refactor (channel-demand struct)

**Review:** §4.1 (31 params, 7× copy-pasted target/weight plumbing,
cyclomatic 96 — the #1 smell; adding an 8th demand type touches ~6 places).

**Fix design.** `struct ChannelDemand<'a> { target: Option<&'a [f64]>,
weight: Option<&'a [f64]> }` + one accessor; the seven `load_pair` copies
become an iteration over `[ChannelDemand; 7]`; the in-loop `want_*` ladders
collapse to loops over an enum; `PointOut` fields likewise. At the PyO3
boundary, introduce a `NeedleRequest` builder (R4.2's new optional grads
slots in cleanly here — sequence R4.2 first so the refactor absorbs it once,
not twice).

**Risk.** M — must remain **bit-exact**: differential bit-identity test +
all three fd harnesses are the gate; any numeric drift = refactor bug.
Hot-loop allocation profile must not regress (bench_backside as the gate).

**Validation.** Existing: `fd_step1.py`, `fd_rchannel.py`,
`color_merit_check.py`, `color_grad_python.py` (R4.2), needle parity tests,
bit-identity test, `bench_backside_speed`. New: none (this refactor is what
the existing harnesses were built to protect).

**Effort.** XL.

### R6.2 Small physics nits batch (one commit)

- **DOP_R clamp** (§3.3): clamp reflected DOP to ≤ 1 symmetric with DOP_T;
  add one assert-row to the ellipsometry smoke coverage.
- **`+0.0` removals** (§3.3, twice) — cosmetic.
- **`needle_slopes` docstring** (§19.3): τ̂ is `iβ′(1+r₁₂²)/(1−r₁₂²)` — the
  code is right, the doc is stale; one line.
- **Sellmeier pole guard** (§3.5): decide error-vs-documented-NaN (see
  §8.3 optional tests); minimum: doc note stating λ² = Cᵢ ⇒ inf/NaN.

**Validation.** Existing suites; `garbage_in.py` unchanged behavior rows.

**Effort.** S.

### R6.3 Solver triplication (optional)

**Review:** §4.2 — const-generic `LEVEL` consolidation of
`solve_point`/`solve_point_intensity`, single-pol delegating to dual.
Only with benches before/after (per-call overhead is the risk); abort if
`bench_backside_speed` regresses > 2 % or bit-identity breaks. Low priority.

**Effort.** L.

### R6.4 Docs & hygiene batch

From §6.3, §18.4, §24.2 (each small; group into 2–3 commits):
1. `validation/README.md` — refresh the inventory (post-R2.4 counts,
   `regression/**`, `parity/synthesis`, `parity/color`, goldens scope,
   the navette_matrix build step, NEEDS-PORT status).
2. `docs/plans/bug_fix_plan.md` — close or re-scope: BUG-1/2 fixed via the
   Rust port; 🔴 BUG-A now pinned by `test_asymmetric_interface_position_and_pair`
   (passing) — mark RESOLVED with the test reference, or re-scope to the
   Rust module if genuinely open (§6.3 says the plan cannot currently tell).
3. `attic/` — `git rm -r` (history preserves it); re-run CodeRadar smells
   after to confirm the noise drop; keep `docs/plans/` + `validation/**/refs/`
   (parity oracles).
4. Naming drift — decide "Loom" vs "Navette" for docs headers
   (`materials-architecture.md` etc.); minimum: one line in each doc mapping
   the internal name to the public one; full rename only if it stays a sed
   (no code identifiers involved — verified §24.2).
5. Doc-line batch (each verified by this review): Bradford matrix
   row-vector convention (§11); interpolation `extrap` semantics vs scipy +
   `extrap='error'` NaN footgun (§10 — decide: raise `PyValueError` from the
   binding when the flag is set, or document; prefer raising, one-line fix);
   KK 80 eV ceiling (§3.5); coverage-rule asymmetry for color demands
   (§20.3); weaver same-key concurrent-write semantics (§22); seed=None RNG
   asymmetry comment (§9.1); **`SpectralTarget` class docstring: enumerate
   the full spectral-label vocabulary including the differential-phase
   labels `PDts`/`PDtp`** (verified reachable end-to-end and error-guarded
   during the R4.x review pass; today the class docstring says only
   "radians (phase)" + a pointer to `docs/spectralweave-target-kinds.md`,
   while the vocabulary lives in the `synthesis` module docstring — and
   note the transmission-only scope: `reference_phase` supports `passes=2`
   for a reflection round trip, but no reflection differential label
   exists).
6. SPDX headers — add `// SPDX-License-Identifier: LGPL-3.0-or-later` to
   sources missing them (mechanical; `.py` + `.rs`).
7. README performance table — re-audit against §5.2 release numbers (some
   claims may now be understated); stamp build profile.

**Effort.** M total.

### R6.5 `.pyi` stubs for `_smatrix` / `_spectralweave`

**Review:** §24.1 (`_materials.pyi` verified zero-drift; the other modules
are unstubbed). Build stubs from the PyO3 registration lists the same way
the drift audit did; add a CI check (extend the §24.1 extraction into
`tools/check_pyi_sync.py`, wire into R2.3's ci.yml) so drift cannot
recur — the same pattern as `check_cie_sync.py`.

**Effort.** M.

---

### R6.6 Rename `color/func_NN.rs` → descriptive module names

**Review:** §11 (the color catalog is organized as `func_01`–`func_16` after the original port task numbering).

**What.** `rust/navette/src/color/func_01.rs` … `func_16.rs` carry opaque names; the only map is the doc-comment list in `mod.rs:11-26` (func_01 xyY … func_16 CIEDE2000). The sibling files (`common.rs`, `tables.rs`, `matrices.rs`, `metrics.rs`, `composites.rs`, `golden.rs`, `parity.rs`) are already descriptive — the item covers exactly the 16 `func_NN` files. The Python surface is unaffected (bindings call functions, not module paths).

**Proposed mapping** (from the `mod.rs` doc, which becomes the migration record):

| Old | New | Content |
|---|---|---|
| `func_01` | `xyy` | XYZ ↔ xyY |
| `func_02` | `lch` | Lab ↔ LCh |
| `func_03` | `luv` | XYZ ↔ CIELUV |
| `func_04` | `oklab_xyz` | XYZ ↔ Oklab (direct XYZ matrices) |
| `func_05` | `oklab_srgb` | sRGB ↔ Oklab (legacy sRGB matrices) |
| `func_06` | `uvw1964` | CIE 1964 U\*V\*W\* |
| `func_07` | `ucs1960` | CIE 1960 UCS & chromaticity |
| `func_08` | `bradford` | Bradford chromatic adaptation |
| `func_09` | `delta_e_76` | ΔE 76 |
| `func_10` | `delta_e_94` | ΔE 94 |
| `func_11` | `delta_e_cmc` | ΔE CMC(l:c) |
| `func_12` | `din99` | DIN99 |
| `func_13` | `spectral_srgb` | spectral pipeline (SPD × CMF → sRGB) |
| `func_14` | `photometry` | photometry engine |
| `func_15` | `shapes` | shape handling & broadcasting |
| `func_16` | `delta_e_2000` | CIEDE2000 |

**Mechanics.**
1. `git mv` each file (rename detection keeps history; verify with `git log --follow` on one file afterwards).
2. Update `mod.rs`: the 16 `pub mod func_NN;` declarations, the doc-comment list (keep a one-line `func_NN → new_name` mapping table at the bottom of the module docstring for traceability with the port records), and any intra-module `crate::color::func_NN::` paths.
3. Update all external path references — grep-verified sites: `rust/navette/src/smatrix/synthesis/color_merit.rs`, `rust/navette-py/src/color.rs`, plus intra-module uses from `common.rs`/`composites.rs`/`metrics.rs`/`parity.rs`; re-run the grep after the rename to confirm zero remaining `func_` module-path hits (`smatrix/optics_core.rs` also matched the grep — verify whether it is a real path or a coincidental doc string before renaming).
4. Check `tools/check_exposure.py` (its allowlist may key on module paths) and the CodeRadar index (`codegraph_reindex` after the change) — both must reflect the new names.
5. **Do not touch** `docs/plans/color/func_NN_*.md` — they are historical port-task records under their original names; instead add the same old→new mapping row to `docs/plans/color/README.md` so the records stay findable. The review doc's §11 catalog keeps the func_NN labels as historical text.

**Risk.** S — internal module paths only, no public API change (verify: after the rename, `check_exposure.py` passes unchanged and the Python suite is green untouched). Mechanical but wide; single-purpose commit, nothing else mixed in.

**Validation.** Existing: `cargo test --workspace` (376 — the color units live in `#[cfg(test)]` inside the renamed files and move with them), `pytest validation` (355 — Python surface must be bit-identical), `check_exposure.py`, color benches untouched. New: none required — this is a pure rename; the verification is the unchanged-suite run plus the zero-remaining-`func_` grep.

**Effort.** S–M.

---

## 7. Test strategy — master view

### 7.1 Existing suite inventory (what guards what)

| Suite | Size | Guards |
|---|---|---|
| `pytest validation` (smoke + goldens + regression) | 355 | imports, goldens, structure round-trips, mirror semantics |
| `cargo test --workspace` | 376 | engine units, materials parity (22), bindings (15), synthesis deposit semantics |
| pytest-style parity tests (`validation/parity/**`) | 13 | Python-vs-Rust semantic parity — **hidden until R2.4**; the last two NEEDS PORT rows closed by R2.4a (0.6.11) |
| review harnesses (`validation/review/`) | 8 scripts | independent physics: fold eval, FD gradients, color kernel, garbage-in, weaver race, TL golden, KK machinery |
| benches (`validation/benches/`) | 7 | timing + numeric parity sections (release-only after R2.1) |

### 7.2 New tests required by this plan (committed with their item)

| Path | Item | One-line purpose |
|---|---|---|
| `validation/smoke/test_rtype5_energy.py` | R1.1 | roughness energy-conservation sweep + NC worked example |
| `validation/smoke/test_energy_conservation.py` | R1.2 | 1-D/2-D/lossless/absorbing residual helper |
| `validation/smoke/test_input_validation.py` | R3.1, R3.2 | construction-time rejection matrix + needle z-range |
| `validation/smoke/test_request_bits.py` | R2.5 | 49 REQ semantic probe + 18 NREQ exact + schema |
| `validation/smoke/test_eigenmodes.py` | R3.3 | SPP pin, field profile, runaway rejection |
| `validation/smoke/test_examples.py` | R4.3 | executes every `examples/*.py` |
| `validation/review/color_grad_python.py` | R4.2 | Python-path color gradient vs native + hand chain rule |
| `validation/benches/_bench_common.py` | R2.1/2.2 | release-assert + UTF-8 setup shared by benches |
| parity: R-row in `test_needle_t_a_phi.py` `check_fd` | §19.2 | permanent R-channel FD coverage with a distinct needle material (one-line addition, do it with R1.1's commit or R2.4) |
| `validation/review/lm_check.py` | R4.4 | engine optimizer vs `scipy.optimize.least_squares(trf)` parity (closes §18.2's never-reproduced claim), parametrized over enabled backends |
| cargo: backend-equivalence tests | R4.4 | same problem through `BuiltinLm`/`MinpackLm`/`Trf`/`Argmin*` — optimum/cost agreement within tolerance |
| cargo: analytic-vs-FD per-element J cross-check | R4.5 | deposit-assembled J vs FD J on a representative spec (rel ≤ 1e-6 away from fold kinks) + channel-deposit coverage test |
| `lm_check.py` TRF rows | R4.6 | TRF vs `scipy.optimize.least_squares(method="trf")` — same algorithm family, the tightest parity in the plan |

### 7.3 Optional additional tests (recommended, not blocking)

These close the §18.3/§6.4 leftovers; write them opportunistically, one per
PR, after Phases 1–3:

1. **Property test — reciprocity:** random lossless stacks, assert
   `R_forward(n↔substrate)` identities and `T_p/s` reciprocity at multiple
   angles (§6.4's suggestion; cheap, catches branch-sign regressions).
2. **Monte-Carlo statistics:** seeded `apply_errors` vs the Python expander
   distribution (mean/σ per row, fixed seed) — §18.3's never-run claim.
3. **Pickle round-trip:** `UniInterpolator.__reduce__` (§18.3).
4. **Multiblock + hybrid A/B/C hand-computed cases:** single thin block ×
   known interference, one thick spacer (analytic incoherent limit) — the
   last parity-only physics (§18.2) without an independent check.
5. **FH/Sprague interpolation vs scipy** (`FloaterHormann`, `sprague`) —
   same protocol as §10's pchip/makima (§18.1).
6. **Coverage tooling:** `pytest --cov` config + `cargo llvm-cov` in CI
   (§6.4) — measure before optimizing anything in Phase 5.
7. **`get_weaved` key-probe API** (§22 minor polish) + its test.

---

## 8. Golden & parity regeneration protocol (applies to R1.1 and any future physics change)

1. **Never** regenerate goldens from a reference before fixing the
   reference's own bug (the §3.2 lesson — goldens generated from a buggy
   reference are bug archives, not oracles).
2. Fix engine + in-repo mirror (loom) **in the same commit** where mirror
   tests exist; keep the divergence comment.
3. Triage in this order: (a) which tests pin the changed numbers; (b) for
   each, verify the delta matches the model change analytically (one
   hand-check per file — for R1.1, the worked R+T = 1.000041 example); (c)
   update pins with `// updated: <item id> (<review §>)` comments; (d) only
   then regenerate any `.npy` goldens, as a **separate commit** with a
   manifest of which files changed and why.
4. The independent-physics test (§7.2 rows) must exist **before** pin
   updates — it is the justification for the new numbers.

---

## 9. Sequencing & dependencies

```
Phase 1 (R1.1, R1.2)          physics fixes; independent tests first
   └─ R1.1 needs the R-row parity addition (§7.2) and mirror patch
Phase 2 (R2.1→R2.5)           guards; R2.3 requires the warning-fix commit
   └─ R2.3 (CI) becomes the gate for everything after; R2.4 widens its scope
Phase 3 (R3.1→R3.3)           behavior changes (0.6.0); R3.1 before R3.2
   └─ R3.1's audit may touch tests R2.4 just un-hid — run both
Phase 4 (R4.1→R4.6)           R4.2 before R6.1 (refactor absorbs the new slots);
                              R4.4b (QR + gain-ratio) before R4.4c (backends);
                              R4.5 (analytic J) after R4.4b, before R4.6 (TRF consumes J)
Phase 5 (R5.1, R5.3, R5.2)    perf; all three DONE. R5.3 before R5.2 turned
                              out to be load-bearing: until construction was
                              fixed the derive loop was only ~23% of the call
Phase 6 (R6.1→R6.5)           R6.1 last among code items (XL, bit-exact gate)
```

Parallelizable: R4.1/R4.3 with anything; R5.1 with R6.x docs; R6.4 with any
phase (docs only).

**Ship gate:** Phases 1–4 together = release 0.6.0 (physics fix + behavior
changes + new tests + CI). Phases 5–6 land on `dev` and ride the next
release.

## 10. Risk register

| Risk | Items | Mitigation |
|---|---|---|
| Golden/parity churn from R1.1 | R1.1 | §8 protocol; mirror fixed in-lockstep; hand-checked expectation |
| API breaks (validation raising) break user code | R3.1 | 0.6.0 changelog with exact rejected-input list; error messages name indices |
| CI red from pre-existing debt | R2.3 | warning-fix commit precedes `-D warnings`; parity skips stay loud, not silent |
| Perf regression from R6.1/R6.3 | R6.x | benches before/after; abort criteria (>2%) written into the item |
| Bit-identity drift from parallelization | R5.2 | §13's per-destination pattern; differential test is the hard gate |
| Review harnesses encode soon-obsolete behavior | R3.1 | update `garbage_in.py` expectations in the same commit as the validation change |
| Termination-semantics change shifts pinned synthesis optima | R4.4 | §8 protocol; cost must be same-or-better; hand-checked expectation per changed pin; built-in stays default until alternatives win on real workloads |
| Analytic-J assembly introduces a cross-point reduction that breaks bit-identity | R4.5 | materialize J per-destination (no reduction); any accumulation keeps ordered reduction per §13; per-element cross-check against FD |

## 11. Explicitly deferred (reviewed, no change required now)

- QR-based LM solve instead of normal equations (§3.6 — scipy parity
  anchors current behavior; add the planned comment only).
- criterion benches + CI perf gate (§5.4) — after R5.2 establishes the
  scaling baseline.
- sdist contents / macOS universal2 audit (§18.4).
- `tools/check_cie_sync.py` source-logic audit (§18.4 — it passes).
- Remaining §18.1 unread-code deep reads (navette-py/structure.rs bindings,
  pipeline stages) — schedule as review work, not fixes.
