# Basin optimizer backends

Build with `maturin develop --release --features opt-basin`, then select
`LmConfig(optimizer="basin_lm_qr")` or `LmConfig(optimizer="basin_trf")`. Both
use Basin 1.15's `Vec<f64>` and `DenseMatrix` implementations, with default
features disabled. The default Navette backend remains `builtin`.

`basin_lm_qr` uses `LevenbergMarquardtQr`, including pivoting and Nielsen
damping. It receives the existing residual and Jacobian through `IntervalMap`.
Analytic Jacobian columns include the map's chain rule; fallback central
differences operate in the transformed coordinates. Like MINPACK, this backend
can stop near a bound because the transformed gradient becomes small. Starting
exactly at a bound can therefore prevent it from finding an interior optimum.

`basin_trf` uses the full `TrustRegionReflective`, not Basin's earlier `Trf`. It
uses Coleman-Li scaling, a rank-aware dense SVD subproblem, and native box
constraints. Free coordinates stay strictly interior; equal bounds fix
coordinates. All-fixed problems need one residual evaluation and no Jacobian.
Fallback differences use Basin's central `BoundedFiniteDiff`, with feasible
one-sided stencils near bounds and zero columns for fixed coordinates. Navette's
finite-bound input contract is retained.

## Configuration mapping

  | Navette setting                       | QR LM                                                                                               | Full TRF                                                                          |
  | ---                                   | ---                                                                                                 | ---                                                                               |
  | `ftol`                                | Both actual and predicted reduction relative to the current cost, with Basin's gain-ratio safeguard | Absolute observed cost change relative to the previous cost                       |
  | `xtol`                                | Attempted step norm relative to the current transformed parameter norm                              | Change between published iterates relative to the current physical parameter norm |
  | `gtol`                                | Infinity norm of `Jᵀr` in transformed coordinates                                                   | Infinity norm of the Coleman-Li scaled gradient over free coordinates             |
  | `gtol_scale_invariant`                | Also enables the Jacobian-column/residual orthogonality test at `gtol`                              | Unused; the native gradient test already accounts for bounds                      |
  | `lambda_init`                         | Initial relative damping `tau`                                                                      | Unused                                                                            |
  | `damping`, `lambda_up`, `lambda_down` | Unused; Basin controls Nielsen damping                                                              | Unused                                                                            |
  | `max_iterations`                      | Maximum fully completed Basin iterations                                                            | Same                                                                              |
  | `max_evals`                           | Residual-call budget checked after initialization and between iterations                            | Same                                                                              |

The tolerances are not identical to Navette's built-in LM or SciPy's TRF. In
particular, Basin's observed TRF cost test has no predicted-reduction or
gain-ratio condition, and its relative step test has no additive `xtol` term.
LM's tolerances apply after the interval transform. No extra absolute cost or
step checks are enabled. Rust callers may set tolerances to zero for exact-zero
tests; Python retains its existing strictly positive validation.

The adapter counts residual calls, including its finite-difference probes. Work
performed internally by a supplied `JacobianSource` remains opaque, as with the
other adapters. Basin's own residual counter does not include the adapter's
probes, so a separate shared counter drives the evaluation budget.
Initialization or the final iteration can exceed `max_evals`; it is not a hard
callback cap. `iterations` excludes an iteration that terminates internally
before Basin marks it complete.

## Results and errors

Costs are `sum(r_i**2)`, twice Basin's internal objective. `gain_ratio` is `NaN`
because Basin does not expose it. `analytic_jacobians` counts successful
requests served by `JacobianSource`, including hybrid sources.

The adapter reads Basin's native convergence diagnostics from
`run_with_solver()` and maps them onto Navette's `LmTermination`:

| Basin outcome | `LmTermination` |
| --- | --- |
| absolute, scaled, or orthogonality gradient test | `Gradient` |
| relative model reduction; explicit observed cost check | `Cost` |
| relative trial step or trust radius; explicit observed step check | `Step` |
| no free parameters (all coordinates fixed by equal bounds) | `NoFreeParameters` |
| converged, native test not recognized by this adapter | `Converged` |
| `MaxIter` | `MaxIterations` |
| `MaxEvaluations` (the adapter's own `max_evals` check) | `MaxEvaluations` |
| `NumericalNoProgress` | `Stalled` |

If several native tests pass, Navette prefers `Gradient`, then `Cost`, then
`Step`, independently of the order Basin reports them. `SolverFailed`, invalid
inputs, inconsistent callback dimensions, and callback errors return an error
with the backend name. Non-finite trial residuals remain available to Basin's
rejection logic.

`basin_lm_qr` is additionally subject to the check every interval-mapped
backend gets in `run_optimizer`: after a convergence reason, each coordinate
within `1e-8` of the box width of a bound is probed once, `√ε` of the width
back inside. If the cost falls there, the bound is not a constrained optimum
and the reason becomes `StuckAtBound`. The probes count towards `evals`; `x`
and `cost` are returned unchanged.

`LmTermination` is `#[non_exhaustive]`, and `is_converged()` / `is_budget()`
give the groups. Selecting a Basin backend without `opt-basin` returns the
rebuild hint; it never falls back.

## Reproducing the comparisons

```sh
cargo test -p navette --features opt-basin
maturin develop --release --features opt-basin,opt-minpack-lm,opt-argmin
python -m pytest validation/smoke/test_optimizer_backend.py
python validation/review/lm_check.py
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python validation/benches/synthesis/bench_optimizer_backends.py \
  --output validation/benches/synthesis/results/basin_backends.json
```

Rust tests cover interior, lower-bound, upper-bound, and corner solutions;
boundary starts; fixed coordinates; zero columns; singular `JᵀJ`; analytic and
differenced Jacobians; budgets; and errors. Python checks exercise actual
thin-film merits and the removal sweep. SciPy comparisons use matching starts
near established optima; distant starts require feasibility and improvement,
since different local minima are legitimate. Rank-deficient tests compare
identifiable combinations rather than arbitrary parameter representatives. For
thickness comparisons, retained films are compared directly at `1e-4 nm`
tolerance. If the removal sweep changes the film count, sub-tolerance films
count as zero on both sides. The physical removal sweep is tested separately
at the usual `2 nm` threshold.

The benchmark warms each case once, times 21 fresh-stack runs, and records a
five-number summary of the wall time (min, quartiles, median, max, in ms)
alongside solver work and build provenance, one JSON line per row. `--raw`
also keeps every sample, for local investigation. Timing is informational and
does not gate CI. The argmin Gauss-Newton failure
demonstration remains in the review harness.

## Measured results and limitations

The [recorded run](../validation/benches/synthesis/results/basin_backends.json)
used Basin 1.15.0, Rust 1.97.0, one thread, and 21 timed runs per row on a
shared cloud VM, so absolute times are slower and noisier than on a desktop;
compare backends within a row, not against other machines. All near-optimum
and active-bound optical comparisons passed at `rtol=1e-9, atol=1e-12` for
cost and `1e-4 nm` for thickness, in both analytic and finite-difference
modes. Times include fresh context and stack construction; they measure
end-to-end optimization rather than the linear solve alone.

Each cell is median wall time in ms, the interquartile range in brackets, and
residual calls. QR LM's calls include the inward probe `run_optimizer` adds
for each coordinate it leaves saturated at a bound.

| Case, analytic Jacobian | Built-in TRF | Basin full TRF | Basin QR LM |
| --- | --- | --- | --- |
| Two films, near optimum | 2.02 [2.00–2.04] / 25 | 2.67 [2.57–2.75] / 28 | 3.79 [3.38–3.99] / 38 |
| Three films, near optimum | 2.92 [2.85–3.02] / 32 | 3.06 [3.02–3.13] / 33 | 1.72 [1.64–1.80] / 18 |
| Upper thickness bound | 0.516 [0.509–0.559] / 7 | 0.581 [0.563–0.602] / 7 | 0.239 [0.231–0.256] / 3 |
| Lower thickness bound | 0.386 [0.378–0.416] / 5 | 0.566 [0.553–0.624] / 5 | 0.418 [0.397–0.438] / 3 |
| One-film interior control | 1.54 [1.42–1.68] / 21 | 1.23 [1.17–1.27] / 21 | 0.997 [0.953–1.02] / 17 |

From a distant `10 nm` start in the one-film interior problem, Basin QR LM
saturates at the `1000 nm` upper bound with cost `249.17`. Moving inward lowers
the cost there, so the bound is not a constrained optimum: the run is reported
as `StuckAtBound`, at the price of one extra residual evaluation. Full TRF
reaches about `103.17 nm` and cost `54.90`. The built-in LM and MINPACK also
reach different local solutions from that start. From `40 nm`, all four reach
the same optimum. The benchmark retains both starts.

Basin 1.15 fixes the full TRF failure at a floating-point cost floor. The Rust
regression test uses `r(x) = [x[0] - 0.3, x[1] - 0.7, 2]`, starting at `[0, 1]`
in the unit box. With default tolerances, both analytic and finite-difference
Jacobians stop about `3e-9` from the optimum because subsequent cost reductions
round to zero. Basin now reports `NumericalNoProgress`, which Navette maps to
`Stalled`, preserving the accurate parameters and cost. This stop does not imply
that a convergence tolerance passed.
