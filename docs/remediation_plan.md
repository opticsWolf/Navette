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
   - Python: `.venv/Scripts/python.exe -m pytest validation` (450 tests at
     0.5.4; `PYTHONIOENCODING=utf-8` no longer needed for the benches after R2.2)
   - Rust: `cargo test --workspace` (376 tests today)
   - Parity: the pytest-style parity tests after R2.4 makes them collectable
   - Review harnesses: `for f in fd_step1 fd_rchannel color_merit_check garbage_in weaver_race tauc_check kk_validate kk_conv; do python validation/review/$f.py; done` (all must print `ALL OK` unless the item explicitly changes documented behavior)
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
| R2.3a | clippy clean → blocking gate | P2 | 248 findings in hot code | M (numeric drift) | M–L | §7 |
| R2.3b | rustfmt adoption → blocking gate | P3 | 981 files; style decision first | S (blame churn) | S/L review | §7 |
| R2.4 | parity tests collected; `sys.exit` → skip | P1 | test suite honesty | M (env dependency) | M | §15, §6.3 |
| R2.5 | request-bit + schema sync tests | P1 | prevents silent corruption | S | S | §4.3, §9.3 |
| R3.1 | `ScatterMatrix` input validation | P1 | silent-garbage class closed | M (behavior change) | M | §21.2, §19.3 |
| R3.2 | needle z-range: debug-assert → error | P1 | release-build garbage | S | S | §21.1 |
| R3.3 | eigenmode `char_func`/`refine_mode` bound | P1 | silent n_eff = −1.7e8 | S–M | M | §16 |
| R4.1 | numpy floor → `>=2.0` | P1 | broken installs | S | S | §7 |
| R4.2 | color gradients on Python needle path | P2 | documented flow incomplete | M (binding change) | M | §20.2 |
| R4.3 | fix + test all `examples/` | P1 | shipped example broken | S | S | §6.2 |
| R4.4 | Optimizer backends: hardened built-in LM + optional argmin-ecosystem solvers | P2 | synthesis robustness on ill-conditioned stacks | M (pinned optima may shift) | M–L | §3.6, §18.2 |
| R4.5 | Analytic Jacobian for the refold optimizer (deposit chain, FD fallback) | P2 | 2n→1 solver sweeps per LM iteration; removes the FD noise floor; benefits every backend | M (fold-kink semantics; ordered accumulation) | M | §3.6, §19.1, §20.2 |
| R4.6 | TRF backend (trust-region-reflective — the bounded-LS reference method) | P2 | correct boundary behavior; direct scipy parity; retires the clamp-prediction caveat | M–L (largest algorithmic lift; scipy as oracle) | L | §3.6 |
| R5.1 | `unweave_collection` batch optimization | P2 | 0.22–0.48× at scale | M | L | §5.3 |
| R5.2 | parallelize serial derive loop | P3 | next Amdahl bottleneck | M (bit-identity) | M–L | §5.4 |
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

### R2.3a Clippy clean → make the clippy gate blocking

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

### R2.3b rustfmt adoption → make the fmt gate blocking

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

### R2.4 Parity tests collected; `sys.exit(1)` → skip

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

**Impact.** The documented run stops silently excluding the best tests.

**Risk.** M — some parity tests may be slow (numba import ~seconds) or
platform-quirky; if a test is genuinely script-shaped (top-level asserts,
no test functions), leave it ignored and list it explicitly in conftest with
a comment (whitelist style beats glob style when mixed).

**Validation.** Existing: the 11 parity tests themselves. New: none needed;
update `validation/README.md`'s inventory (§6.3 rot) in the same commit.

**Effort.** M.

### R2.5 Request-bit and schema sync tests

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

---

## 3. Phase 3 — Input validation & robustness (P1)

### R3.1 `ScatterMatrix` construction validation

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

### R3.2 Needle z-range: debug-assert → real error

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

### R3.3 Eigenmode `char_func` upper bound + `refine_mode` no-mode contract

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

---

## 4. Phase 4 — API & packaging completion (P1/P2)

### R4.1 numpy floor → `numpy>=2.0`

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

### R4.2 Color gradients on the Python needle path

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

### R4.3 Fix and continuously execute `examples/`

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

#### R4.4b Work item A — harden the built-in LM (default backend, zero new dependencies)

1. **QR step solve.** Replace `solve_symmetric(&a, &neg(&jtr), &mut delta)` (thick_opt.rs:373) with a QR solve of the damped augmented system `[J; √λ·D]δ ≈ [−r; 0]` (D = the same diag-scaled Marquardt weights, keeping the current 1e-14 flooring semantics). Two scoping facts: (a) the damped augmentation is **always full-rank for λ > 0** (`[J;√λD]ᵀ[J;√λD] = JᵀJ + λD² ≻ 0`), so pivoting is a *stability/diagnostic* feature, not a mathematical necessity — a plain unpivoted Householder QR (~80 lines) already eliminates the condition-squaring problem, while column-pivoted QR (~200 lines, permutation bookkeeping, dropped-column → zero-step-component mapping) additionally detects and degrades gracefully through degenerate columns (pivot-driven column drop) instead of relying on the current λ-escalation retry after solve failure — which remains as fallback either way. Scope decision: start unpivoted, add pivoting only if the rank-diagnostics are actually needed; (b) whatever the flavor, add the **QR-vs-legacy cross-check cargo test**: on well-conditioned problems the new step solve must reproduce the old normal-equations step to high precision (rel ≤ 1e-8 on δ) — the cheapest strong guard against linear-algebra implementation bugs.
2. **Gain-ratio damping update** (MINPACK, replaces fixed ×5/÷3): ρ = actual/predicted reduction; accept when ρ > 0; on acceptance λ ← λ·max(1/3, 1−(2ρ−1)³), ν ← 2; on rejection λ ← λ·ν with ν ← 2ν. Fewer residual evaluations near the optimum, better behavior on correlated layers.
3. **Clipped-step prediction (the one place a naive MINPACK port goes wrong).** The bound veto+clamp modifies the step *after* it is solved, but the MINPACK predicted reduction is computed for the **full** LM step — using the unclipped prediction systematically overestimates it, so ρ is underestimated whenever a bound is active (over-damping, premature ftol exits near boundaries). The predicted reduction MUST be recomputed for the clipped step from the linear model on the surviving components (−gᵀδ − ½δᵀJᵀJδ restricted to the non-vetoed entries); pin it with a dedicated unit test (bound-active case asserting the clipped prediction < unclipped prediction and the accepted step still reduces cost).
4. **MINPACK termination semantics.** ftol: require **both** actual and predicted reductions ≤ ftol (currently actual-only), with "predicted" understood per item 3. gtol: add the scale-invariant angle criterion (cos∠(Jeᵢ, r) ≤ gtol) alongside the existing ‖Jᵀr‖∞ test (config flag `gtol_scale_invariant: bool`, default true for new configs; the scale-dependent check stays available). Note the two criteria are *different notions of stationarity*: the angle test can exit at a different point than the norm test on flat, noise-dominated valleys — expected, but document it.
5. **Eval accounting can drift either way.** Gain-ratio damping typically reduces evals near the optimum, but the ν-doubling rejection path can burn more evals than the fixed ×5 ladder on adversarial steps. Net effect is empirical — the `bench_refold` gate is the arbiter, not intuition.
6. **Keep unchanged:** the bound veto+clamp contract (documented as Navette's bounds semantics — neither reference has bounds), the FD rayon Jacobian, and the `LmTermination` variants (map new conditions onto existing values where possible to avoid API churn; extend only if a new reason is genuinely distinct).

#### R4.4c Work item B — backend selection + optional ecosystem solvers (feature-gated)

1. New `synthesis/optimizer.rs`: `enum OptimizerBackend { BuiltinLm, MinpackLm, Trf /* R4.6 */, ArgminGaussNewton, ArgminTrustRegion }`; `OptimizerConfig { backend, ftol, xtol, gtol, max_iterations, max_evals, stepbound, … }` mapping 1:1 onto `LmConfig` for the default; one entry point `run_optimizer(problem, x0, bounds, cfg) -> OptimizerResult` with the existing result shape (x / cost / iterations / evals / termination). Residual system stays the injected closure (`MeritSpec::residuals`) — it is already solver-agnostic.
2. Cargo features, **default OFF** (zero new dependencies for standard builds): `opt-minpack-lm = ["dep:levenberg-marquardt"]`, `opt-argmin = ["dep:argmin", "dep:argmin-math"]` (pin `argmin = "0.11"`, `levenberg-marquardt = "0.15"`). Python surface: `optimizer: str = "builtin"` kwarg on the synthesis entry points; absent feature → `PyValueError` with a rebuild hint (same pattern as the module ImportError fallbacks, §13). `design_config.rs` gains the backend field (`deny_unknown_fields` envelope updated in lockstep — the schema-version rule of §9.3 applies).
3. **Bounds contract per backend:** `MinpackLm`/`Argmin*` are unbounded — wrap with an interior logit reparametrization u = atanh(2(x−lb−ε)/(ub−lb−2ε)), x = lb+ε+(ub−lb−2ε)·(tanh(u)+1)/2. The transform is elementwise and the LM Jacobian is FD per parameter, so gradients pass through it exactly. Document explicitly that boundary semantics on these backends (interior parametrization, optimum strictly inside) **differ** from the built-in veto+clamp (optimum may sit on the bound) — this is a contract difference, not a bug; the built-in stays the default for bounded problems.
4. **argmin integration friction, recorded up front:** argmin-math 0.5 ships backends for nalgebra and ndarray ≤ 0.16; navette pins **ndarray 0.17** (§7) → either write the small manual adapter for our residual types (~100 lines; argmin explicitly supports user-supplied implementations) or convert J to nalgebra matrices per iteration (m×n copy per iteration — negligible against the TMM residual cost). Decide at implementation time; do not bump the ndarray pin for it.

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

### R4.5 Analytic Jacobian for the refold optimizer — assemble J from the verified deposit chain

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

### R4.6 TRF backend — trust-region-reflective, the bounded-LS reference method

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

---

## 5. Phase 5 — Performance (P2/P3)

### R5.1 `unweave_collection` batch path

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

### R5.2 Parallelize the serial derive loop

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
| pytest-style parity tests (`validation/parity/**`) | 11 (+2 NEEDS PORT) | Python-vs-Rust semantic parity — **hidden until R2.4** |
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
Phase 5 (R5.1, R5.2)          perf; R5.2 gated by the optional scaling bench
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
