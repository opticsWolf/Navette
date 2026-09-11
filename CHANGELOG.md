# Changelog

All notable changes to Navette are recorded here. Work items reference
`docs/remediation_plan.md` (Rx.y) and `docs/code_review.md` (§).

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
