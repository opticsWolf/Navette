# validation/

All test, parity, benchmark, golden and reference material in one place.
`pytest validation` runs the collected suites; parity/bench scripts are
standalone (they execute workloads on import) and must be run explicitly.

```
validation/
├── conftest.py               # collects smoke/ + goldens/, ignores parity/ + benches/
├── smoke/
│   └── test_navette_imports.py        # pytest: import surface (6 tests)
├── goldens/
│   └── spectralweave/
│       └── test_golden.py             # pytest: pinned numeric goldens (22 tests)
├── parity/                   # numba/NumPy-vs-Rust parity scripts (run directly)
│   ├── smatrix/
│   │   ├── test_w_function.py                    # PASS, ~0.4x (µs kernel)
│   │   ├── test_redheffer_product_real.py        # PASS, ~1.7x
│   │   ├── test_redheffer_product_complex_field.py # PASS, ~1.7x
│   │   ├── test_solve_coherent_block_fields.py   # PASS, ~1.1x
│   │   ├── test_core_engine_photometry_only.py   # NEEDS PORT (old entry point)
│   │   ├── test_core_engine_rigorous_ellipsometry.py # NEEDS PORT
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

# parity scripts (from the repo root, .venv active)
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
  `core_engine` API (`navette._smatrix.core_engine`) and need porting;
  everything else in `parity/` passes against the release build.
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
