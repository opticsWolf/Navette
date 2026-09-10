# validation/

All test, parity, benchmark, golden and reference material in one place.
`pytest validation` runs everything except the benches — **including
`parity/`**, which was blanket-ignored until R2.4 (0.5.8). The parity files
are dual-shaped: they still run as standalone scripts, and each also exposes
a pytest entry point that asserts the verdict.

```
validation/
├── conftest.py               # collects everything but benches/ (+ refs/, gen_*)
├── smoke/
│   ├── test_navette_imports.py        # import surface
│   ├── test_build_profile.py          # release-build gate (R2.1)
│   ├── test_rtype5_cross_path.py      # Nevot-Croce, all four sites (R1.1)
│   ├── test_rtype5_energy.py          # roughness energy balance
│   └── test_energy_conservation.py
├── goldens/
│   └── spectralweave/
│       └── test_golden.py             # pytest: pinned numeric goldens (22 tests)
├── regression/               # differential + bit-identity suites
│   └── color/ config/ structure/ synthesis/
├── parity/                   # numba/NumPy-vs-Rust parity (collected AND runnable)
│   ├── _parity.py            # console_utf8 / require_reference / report
│   ├── smatrix/
│   │   ├── test_w_function.py                    # PASS, ~0.4x (µs kernel)
│   │   ├── test_redheffer_product_real.py        # PASS, ~1.7x
│   │   ├── test_redheffer_product_complex_field.py # PASS, ~1.7x
│   │   ├── test_solve_coherent_block_fields.py   # PASS, ~1.1x
│   │   ├── test_core_engine_photometry_only.py   # SKIPS: NEEDS PORT (R2.4a)
│   │   ├── test_core_engine_rigorous_ellipsometry.py # SKIPS: NEEDS PORT (R2.4a)
│   │   └── refs/
│   │       └── loom_matrix.py         # numba reference (all kernels)
│   └── materials/
│       ├── gen_goldens.py             # regenerates goldens/ via NumPy
│       └── goldens/*.npy              # read by crates/navette-materials/tests/parity.rs (22 tests)
└── benches/                  # timing scripts (run directly)
    ├── _bench_common.py      # UTF-8 console + release-build gate (import first)
    ├── smatrix/
    │   ├── bench_backside_speed.py    # request-mask throughput
    │   └── results/*.json             # committed timings (see note below)
    ├── structure/
    │   └── bench_grid_assert.py       # provider-grid assertion cost
    ├── synthesis/
    │   └── bench_refold.py            # per-cycle needle-target re-fold cost
    ├── spectralweave/
    │   ├── navette_spectral_bench.py  # [--quick] Python-vs-Rust weave/unweave
    │   ├── navette_target_bench.py    # [--quick] ingest + merit timings
    │   └── refs/
    │       └── loom_spectraldata.py   # pure-Python reference weaver
    ├── interpolate/
    │   ├── 1dinterpol_test_bench.py   # SciPy/loom/Rust interpolation shootout
    │   └── refs/
    │       └── loom_unispline.py      # numba reference interpolator
    └── color/
        ├── bench_validate_color.py    # vs colour-science + loom (needs `colour` pkg)
        ├── refs/
        │   └── loom_colorengine.py    # numba reference color engine
        └── gen/
            ├── gen_golden.py / gen_matrices.py  # regenerate golden.rs
            └── golden.rs                  # copy of crates/navette-color/src/golden.rs
```

## Commands

```powershell
# pytest suites (smoke + spectral goldens)
pytest validation

# Rust suites (incl. materials parity against parity/materials/goldens/)
cargo test --workspace

# parity — collected by `pytest validation`; also runnable directly
# (from the repo root, .venv active; needs `numba` for the reference)
python validation/parity/smatrix/test_w_function.py
# ... etc.

# benches (require a --release extension; they exit on a debug build)
python validation/benches/smatrix/bench_backside_speed.py
python validation/benches/structure/bench_grid_assert.py
python validation/benches/synthesis/bench_refold.py
python validation/benches/spectralweave/navette_spectral_bench.py --quick
python validation/benches/spectralweave/navette_target_bench.py --quick
python validation/benches/interpolate/1dinterpol_test_bench.py
python validation/benches/color/bench_validate_color.py

# regenerate materials goldens (then: cargo test -p navette-materials)
python validation/parity/materials/gen_goldens.py
```

## Notes

- `parity/smatrix/test_core_engine_*.py` predate the request-driven
  `core_engine` API (`navette._smatrix.core_engine`) and need porting. They
  now **skip with that reason** instead of calling `sys.exit(1)` at import,
  which used to turn any collection of `parity/` into a pytest
  INTERNALERROR. They load a standalone `navette_matrix` extension from
  `parity/target/release/` — a crate that does not exist anywhere in this
  repo, so building it is not the fix; porting onto `core_engine` is
  (**R2.4a**). Everything else in `parity/` passes against the release build.
- **The parity comparisons were not running at all before 0.5.8.** The numba
  reference was imported from `parity/loom/`, a directory that does not
  exist; the import failed, and the scripts responded by scoring every
  comparison `"PASS (rust-only)"` — passing while comparing nothing. They
  also printed `OUTPUT_STATUS FAIL` and exited **0**, so a failure was
  invisible to any caller. The reference now loads from `smatrix/refs/`
  (requires `numba`), a missing reference is a **skip**, and a failure exits
  non-zero and fails under pytest. Verified by injecting a 1e-3 error into
  the Rust side: script `rc=1`, pytest FAILED.
- Legacy pre-unification scripts that tested removed modules (`smatrix`,
  `navette_interpolator`, `request_flags`) were deleted as unportable;
  git history retains them.
- Speedups quoted above are illustrative (Windows, release LTO build);
  re-run on your hardware.

## Benchmarking: build profile

Every bench imports `_bench_common` before `navette` and calls
`require_release()`, which exits unless `navette.build_profile() == "release"`.

This exists because the dev venv once shipped a plain `maturin develop`
(debug) extension: it imported and computed correctly, so nothing looked
wrong, but it is several times slower than the release build and every timing
taken against it was worthless. Concretely, `bench_backside_speed.py`'s
`A_legacy` case measures **~0.06 ms** on a release build against the
**0.50 ms** recorded in `smatrix/results/backside_speed_*.json` — those
committed files predate the gate and are of **unknown profile**; do not
compare new runs against them. Results written from now on carry a
`"_provenance"` block (profile, version, Python, platform).

`_bench_common.setup_bench()` also reconfigures stdout/stderr to UTF-8 (the
benches print `->` arrows and emoji, which raised `UnicodeEncodeError` on a
cp1252 console) and puts the repo's `src/` on `sys.path`, so a bench no longer
has to be launched from the repo root.
