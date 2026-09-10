# Changelog

All notable changes to Navette are recorded here. Work items reference
`docs/remediation_plan.md` (Rx.y) and `docs/code_review.md` (§).

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
