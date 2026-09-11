# validation/

All test, parity, benchmark, golden and reference material in one place.
`pytest validation` runs everything except `benches/` and `review/` — **including
`parity/`**, which was blanket-ignored until R2.4 (0.5.8). The parity files are
dual-shaped: they still run as standalone scripts, and each also exposes a
pytest entry point that asserts the verdict.

692 tests collected at 0.6.18: 287 smoke, 327 regression, 56 parity, 22 goldens.

```
validation/
├── conftest.py               # collects everything but benches/, review/ (+ refs/, gen_*)
├── smoke/                    # 287 tests — API surface, guards, one-behaviour checks
│   ├── test_navette_imports.py        # import surface
│   ├── test_build_profile.py          # release-build gate (R2.1)
│   ├── test_request_bits.py           # request-mask ↔ schema sync (R2.5)
│   ├── test_input_validation.py       # ScatterMatrix + needle + core_engine
│   │                                  #   refusals (R3.1/R3.2/R5.3), Sellmeier
│   │                                  #   domain + DOP_R clamp (R6.2)
│   ├── test_rtype5_cross_path.py      # Névot-Croce, all four sites (R1.1)
│   ├── test_rtype5_energy.py          # roughness energy balance
│   ├── test_energy_conservation.py    # 2-D wrapper (R1.2)
│   ├── test_eigenmodes.py             # char_func / refine_mode bounds (R3.3)
│   ├── test_analytic_jacobian.py      # analytic-vs-FD J, coverage fallback (R4.5)
│   ├── test_optimizer_backend.py      # backend names, availability, refusals (R4.4c)
│   ├── test_lm_parity.py              # built-in LM vs scipy
│   ├── test_color_needle_python.py    # colour gradients on the Python needle path (R4.2)
│   ├── test_spectralweave_wrapper.py
│   └── test_examples.py               # every script in examples/ runs (R4.3)
├── goldens/
│   └── spectralweave/test_golden.py   # pinned numeric goldens (22 tests)
├── regression/               # 327 tests — differential + bit-identity suites
│   ├── color/                #   CIE table sync
│   ├── config/               #   program / pipeline config round-trips
│   ├── structure/            #   CRUD, providers, inversion, bit-identity
│   └── synthesis/            #   colour merit + targets, graded pipeline
├── parity/                   # 56 tests — reference-vs-Rust (collected AND runnable)
│   ├── _parity.py            #   console_utf8 / require_reference / report
│   ├── smatrix/              #   numba reference, refs/loom_matrix.py
│   │   ├── test_w_function.py
│   │   ├── test_redheffer_product_real.py
│   │   ├── test_redheffer_product_complex_field.py
│   │   ├── test_solve_coherent_block_fields.py
│   │   ├── test_physics_mirror.py
│   │   ├── test_backside_solver.py
│   │   ├── test_needle_t_a_phi.py
│   │   ├── test_core_engine_photometry_only.py        # whole-engine (R2.4a)
│   │   └── test_core_engine_rigorous_ellipsometry.py  # whole-engine (R2.4a)
│   ├── synthesis/            #   merit bridge/mirror, needle refold, weighting,
│   │                         #   differential phase, pipeline
│   ├── color/                #   golden mirror of color/golden.rs
│   └── materials/
│       ├── gen_goldens.py    #   regenerates goldens/ via NumPy
│       └── goldens/*.npy     #   read by rust/navette/tests/materials_parity.rs
├── review/                   # standalone harnesses — NOT collected by pytest;
│   │                         # each exits 0/1 and prints its own verdict
│   ├── fd_step1.py           #   fold + end-to-end gradient vs FD
│   ├── fd_rchannel.py        #   R-channel needle slopes vs FD
│   ├── lm_check.py           #   optimizer parts A–E (E needs opt-argmin)
│   ├── color_merit_check.py / color_grad_python.py
│   ├── kk_validate.py / kk_conv.py    #   Kramers–Kronig
│   ├── tauc_check.py
│   ├── garbage_in.py         #   input-abuse survey
│   └── weaver_race.py        #   concurrent weaver writes
└── benches/                  # timing scripts (run directly)
    ├── _bench_common.py      #   UTF-8 console + release-build gate (import first)
    ├── smatrix/
    │   ├── bench_backside_speed.py       # request-mask throughput
    │   ├── bench_core_engine_scaling.py  # 500 → 60 000 points vs numba (R5.3)
    │   └── results/*.json                # committed timings (see note below)
    ├── structure/bench_grid_assert.py
    ├── synthesis/bench_refold.py         # per-cycle needle-target re-fold cost
    ├── spectralweave/                    # [--quick] weave/unweave + target ingest
    ├── interpolate/1dinterpol_test_bench.py
    └── color/
        ├── bench_validate_color.py       # vs colour-science + loom (needs `colour`)
        └── gen/gen_golden.py, gen_matrices.py, golden.rs
```

Every `refs/` directory holds a reference implementation (numba or pure Python)
that the parity tests measure against. They are inputs, not tests, and
`conftest.py` excludes them from collection along with `gen_*` generators.

## Commands

```powershell
# everything pytest collects (smoke + regression + parity + goldens)
pytest validation

# Rust suites (incl. materials parity against parity/materials/goldens/)
cargo test --workspace

# parity — collected by `pytest validation`; also runnable directly
# (from the repo root, .venv active; needs `numba` for the reference)
python validation/parity/smatrix/test_w_function.py

# review harnesses — not collected; run one, or all of them
python validation/review/lm_check.py

# benches (require a --release extension; they exit on a debug build)
python validation/benches/smatrix/bench_backside_speed.py
python validation/benches/smatrix/bench_core_engine_scaling.py
python validation/benches/synthesis/bench_refold.py
python validation/benches/spectralweave/navette_spectral_bench.py --quick

# regenerate materials goldens (then: cargo test -p navette)
python validation/parity/materials/gen_goldens.py
```

## Notes

- `parity/smatrix/test_core_engine_*.py` were ported onto the request-driven
  `core_engine` in **R2.4a (0.6.11)**; they had been skipping since 0.5.8
  because they loaded a standalone `navette_matrix` extension from
  `parity/target/release/`, a crate that does not exist anywhere in this repo.
  Both now compare every channel of their legacy 4- and 13-tuples against the
  numba reference (~1e-14 to ~1e-16), including `conservation_err`, which is
  reconstructed from `solver_energy_conservation` (R1.2). No NEEDS-PORT rows
  are left.
- **These two are also the only whole-engine speed comparison.** They ran
  0.64–0.77× the numba kernel until **R5.3 (0.6.14)** rebuilt `Solver::new`;
  `bench_core_engine_scaling.py` now measures ~0.8–1.2× at 20 000–60 000 points
  and still ~0.3–0.5× at 500, where per-call setup dominates. The single-point
  ratios these two scripts print (~1.0× photometry, ~0.45× rigorous) are a
  small-grid measurement — read the scaling bench, not these, for the shape.
- **Timing here is contaminated by whatever ran before it.** numba's thread
  pool keeps spinning after its block returns: the same Rust call at 20 000
  points measures 1.52 ms alone, 2.00 ms immediately after a numba block, and
  1.58 ms after a one-second pause. Both parity scripts call `cooldown()`
  before each timing loop for that reason. `NUMBA_THREADING_LAYER=workqueue`
  also removes it but makes numba ~70 % slower, which would flatter the Rust
  side instead.
- **The parity comparisons were not running at all before 0.5.8.** The numba
  reference was imported from `parity/loom/`, a directory that does not exist;
  the import failed, and the scripts responded by scoring every comparison
  `"PASS (rust-only)"` — passing while comparing nothing. They also printed
  `OUTPUT_STATUS FAIL` and exited **0**, so a failure was invisible to any
  caller. The reference now loads from `smatrix/refs/` (requires `numba`), a
  missing reference is a **skip**, and a failure exits non-zero and fails under
  pytest. Verified by injecting a 1e-3 error into the Rust side: script `rc=1`,
  pytest FAILED.
- Legacy pre-unification scripts that tested removed modules (`smatrix`,
  `navette_interpolator`, `request_flags`) were deleted as unportable; git
  history retains them.
- Speedups quoted here are illustrative (Windows, release LTO build); re-run on
  your hardware.

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
