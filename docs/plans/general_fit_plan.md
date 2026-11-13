# General fit engine: thin films, groups and dispersion models together

Status: **reviewed** (2026-11-07); decisions in §0. Refreshed
2026-11-13 against 0.9.0 (synthesis now runs through incoherent plates and
scores ratio, Ψ, retardance, GD…FOD and Stokes/DOP targets; §1, §5.4).
Target release: **0.10.0** (phase 0 renames public names, a breaking
change). Work branch: `dev_fit`.

| Phase | State |
|---|---|
| 0 rename | done (`dev_fit`) |
| 1 parameters, `StackForward`, one stage | done (`dev_fit`) |
| 2 dispersion fitting | done (`dev_fit`) |
| 3 strategies, environments, ties, priors, converter | done (`dev_fit`) |
| 4 tournament, generated strategies, truth draws, benchmark | done (`dev_fit`) |
| 5 Python door, program `fit` section (schema v4), docs | done (`dev_fit`) |
| 6 fit Jacobian (block skipping, analytic thickness and nk columns), `global` stage on basin | done (`dev_fit`) |
| 7 | not started |

Phase 1 settled §3.4.1's open item: the fit owns a `Send + Sync` mirror
(`FitModel`: the architect's traversal with each entry's structure row,
the merged groups, the material entries by value), and an evaluation
clones it, writes the parameters, and expands. A layer address names a
*structure* row, so a shared or repeated structure is one parameter.
Recovery is calibrated: over 40 noise seeds, 118 of 120 fitted values
of the three-parameter test land within 2σ, mean χ²_red 0.999.

Phase 2 kept dispersion fitting inside the same engine: `NkData` rows
follow the spectra in one residual vector, so a joint fit needs nothing
extra. Choices made here:

- AIC/BIC are the noise-scale-estimated Gaussian forms
  (`m·ln(χ²/m) + k·penalty`), so `table_to_model`, whose table has no σ,
  ranks independently of the σ it assumes.
- `fit_dispersion` frees every number by default, except a `Table`'s data
  and a mixer's constituents (an air inclusion's `n` would only trade
  against the fraction).
- A trial point the forward model refuses becomes a rejected step
  instead of a failed fit.
- Bruggeman's Newton stop (|Δε| ≤ 1e-9) leaves an error of order 1e-18,
  so the FD derivative in the fraction matches the closed form to 1e-6
  relative for steps from 1e-3 to 1e-6, and to 1e-5 (round-off) down to
  1e-8; the fit's own step is about 2e-8.

Calibration over 200 seeds: Cauchy 564/600 and Tauc-Lorentz (ε data)
950/1000 within 2σ, mean χ²_red 1.000. The Sellmeier test holds `C2`,
because the two UV poles trade off almost exactly.

Phase 3 kept every stage an ordinary `fit` of a copy of the problem:
parameters outside the stage's selection are held at the model's current
values, so a strategy adds no second solver and a stage's report is a
plain `FitReport`. Choices made here:

- Environments are named chains over one flat entry list; a structure
  row keeps one address however many samples use it, fixed rows
  (ambient, substrate, a cover) are held and not addressable, materials
  and groups are shared.
- A tie's follower is never a parameter, has one driver, and a driver
  is never itself a follower. Priors are rows after the data, reported
  as the `prior` dataset.
- The converter presets are data (`fit/conventions.json`) with no vendor
  names; σ defaults to unit weight and is set per call. tan Ψ/cos Δ
  returns Δ in [0, π].

Deferred from phase 3: declaring parameters by selector (enumerating a
model's scalars for `{"select": ...}`), the `OptMask` default selection,
and a per-environment `coherence_mode` (one mode for the problem today).

Phase 4 (`fit::explore`) drew its truths in **parameter space**, not
through `expand` with errors on: `expand` perturbs every traversal entry,
each layer's index included, and most of that has no address, so a
recovery score would have measured the model's blind spots instead of
the strategy. A free parameter is drawn from its group's law (a layer
channel as `expand` draws it; a group factor from the law's relative
part, a summand from its absolute part); everything else holds. Choices
made here:

- `StrategyGen::Groupings` takes `variants` (one strategy each, a list of
  blocks of sets) rather than one `blocks` list.
- Rankings are deterministic: a benchmark ranks by failures, mean
  recovery error (relative per parameter, RMS), evaluations, name; a
  tournament by χ², evaluations, name. Wall times are reported, never
  ranked on.
- Budget (open question 2): orderings and cumulative orders run in full
  up to 24 (`DEFAULT_ORDERING_CAP`) and are sampled with seed 0 above
  that, with a note; seeds are the caller's choice.

`cargo run --release -p navette --example fit_tournament` ranks seven
strategies on four truths of a 20-layer chirped edge filter (22 free
parameters) in about 4 s and checks a second run gives the same bits.
Freeing the group factors first recovers the truth best (mean relative
error 2e-3); freeing everything at once lands in a wrong minimum on one
of the four truths.

In phase 5 the request became a program document. A program with a `fit`
section (schema v4) is what Python sends to Rust, what `FitProblem.save`
writes, and what `fit::fit_program` runs, so a request and its saved fit
use one format. Choices made here:

- Program structures load as films only. A measured sample (ambient,
  substrate, stack) is therefore declared under `fit.environments`, and
  data with no `env` goes to the first declared environment.
- A `select` declaration expands over `FitModel::addresses` when the
  request is read. Addresses declared explicitly elsewhere in the list
  are skipped, so an explicit declaration wins.
- The doors (`fit::doors`) take JSON text and return JSON, with NaN and
  Inf written as null. Python only wraps them in dataclasses, and their
  shapes are tested in Rust.
- Progress is a `Watcher` with three events: stage start, χ² improved,
  and stage end. When its callback errs, the watcher refuses further
  trial points and the fit ends with that error. Python raises it as a
  `ValueError`.
- `serde_json` `float_roundtrip` is on, so a saved fit reloads
  bit-exactly. Without it, the last ulp moved.

`examples/fit_coating.py` and `examples/fit_optical_constants.py` run in
the battery. The second shows that a six-term Sellmeier fitted over
250–2000 nm reproduces n to 5e-5 RMS, while B1 alone has an error of
±0.19. Assert on the curve, not on the coefficients.

Phase 6 gave the fit its own Jacobian source (`fit::jacobian`); the old
path stays as `lsq.jacobian = "fd"`. What it does and what it found:

- **Block skipping.** A column recomputes only the environments and nk
  datasets its parameter (or the followers it drives) can reach; reach is
  read from the expansion, so gradients and anything `expand` asks the
  provider for count. Differenced entries are bit-identical to the old
  Jacobian, unreached rows exactly zero (a test checks the bits).
- **Thickness columns, analytically.** Not via `simulate_with_deposits`
  or `SmatrixContext` (§3.7's two options): a column is classified per
  environment by expanding at u±h without solving. If only slot
  thicknesses moved, the column is `Σ ∂X/∂d_slot · ∂d_slot/∂u`, the first
  factor read from `needle_operator`'s `thickness_slopes_ddz` on one
  dual field build per point for all slots, the second the (linear)
  expansion differenced. This covers films, `thick_factor`,
  `thick_summand`, graded and sliced layers with no per-kind code. It
  needs a coherent interior and R/T/A channels; Psi/Delta and phases fall
  back to block differences. All thickness columns of an environment cost
  about three evaluations together.
- **Dispersion columns** on nk data are closed form for Konstant, Cauchy,
  Sellmeier and Drude top-level numbers. Urbach and Lorentz were left:
  Lorentz oscillators live in arrays (pointer addresses), and nk-only fits
  are microseconds per evaluation already.
- **Index perturbation** (index columns on spectra) was not done. They
  are the largest remaining cost: in `coat20` the two Cauchy `A` columns
  take four of the Jacobian's ~8 evaluations. The needle function's depth
  integral over a layer is the route, for a later phase.
- **Exact columns beat the oracle where the oracle is noisy.** Central
  differences with the default scale step a summand that starts at 0 by
  6e-9 nm, and an `A_avg` of a lossless stack read 1e-8 of rounding noise
  where the true slope is 0. Tests compare analytic and FD at 1e-5 of each
  column's scale, with parameters whose FD step is sane.
- `examples/fit_perf.rs` runs every case both ways. Best of three on
  this machine (32 threads): coat20 (22 free, 453 residuals) 53 → 18 ms,
  2.9×; coat20 with witness samples 54 → 25 ms, 2.2×; coat20 with an nk
  table 37 → 20 ms, 1.9×; same iterations and χ² in each. Single-threaded,
  a coat20 Jacobian is 9.9 evaluations instead of 44. Fits on nk data
  alone (0.2 ms) are unchanged within 10 %. The Sellmeier nk case still
  hits the iteration cap with exact derivatives: it is ill-conditioned
  (B/C trade-off), not an FD-accuracy problem.
- **`global` stage** (`fit::global`, `opt-basin`): basin's `De`,
  `DeInject` (Nelder-Mead polish) or `BoundedCmaEs` over the stage's box
  scaled to [0, 1]ⁿ, χ² compiled once, then the stage's LSQ from the best
  point. `tests/fit_global.rs`: on the phase-4 coating from six random
  starts in a 0.5–1.8× thickness, ±10 % index box, all three methods end
  at or below χ²(truth); plain LSQ does too there. On the same design
  four times thicker, read through one curve, LSQ is caught on all six
  starts (χ² thousands of times χ²(truth)) and DE on none.
- Found on the way: a strategy's final report re-applied the declared
  start values over the fitted model when a parameter had an explicit
  `value`. Fixed; phase-4 problems never set one, so it had not shown.

## 0. Decisions (review of 2026-11-07)

| # | Question | Decision | Where |
|---|---|---|---|
| 1 | Rename policy | **Hard rename, no backwards compatibility**: no aliases, no re-exports, no refusal shims for old names | §2.2 |
| 2 | First door | **Rust first**; Python binds the finished Rust API | §3.10, phase order |
| 3 | Ellipsometry / instrument conventions | the engine works **in the solver's conventions only**; loaded instrument data passes through a **converter** into those conventions | §3.5.5 |
| 4 | Point-wise / spline nk models | in the plan as a **deferred model option** | §4.6, phase 7 |
| 5 | Global search | **both**: basin's DE / CMA-ES behind `opt-basin` first, as the reference; then a hand-rolled DE in the default build, checked against it | §3.8.3, phases 6-7 |
| 6 | Strategy comparison | tournament as planned, plus **parameter sets** (named groups of parameters) whose order and grouping the tournament can rearrange | §3.3.5, §3.8.1 |

This plan covers three requests:

- **(a)** Give today's optimizer path a name that says what it does
  (thickness refinement inside needle synthesis), so it cannot be confused
  with the general engine.
- **(b)** A general fit engine that varies any scalar in a model (film
  thicknesses, grading, roughness, group factors, material dispersion
  parameters, EMA fractions) and runs them in **ordered stages**, so that
  fitting strategies for reverse engineering can be compared.
- **(c)** Fitting dispersion models (including EMA mixtures) to dispersion
  data (n, k / ε₁, ε₂ tables), which has no path today.

The short version: the solvers are already general, and most of the
modelling machinery can be reused unchanged. What is missing is a
**parameter layer**, a **forward model built on the general solver** (not on
the synthesis one), **measurement residuals** (incl. ellipsometry), a
**strategy runner**, and **fit statistics**.

---

## 1. Where things stand

### 1.1 What the current optimizer reaches

| Piece | File | What it does | General? |
|---|---|---|---|
| LM solver | `smatrix/synthesis/thick_opt.rs` | bounded LM over a residual closure; FD or analytic Jacobian; `LmConfig`, `JacobianSource`, `LmTermination` | **yes**, despite the name |
| TRF solver | `smatrix/synthesis/trf.rs` | bounded trust-region reflective | **yes** |
| Backend dispatch | `smatrix/synthesis/optimizer.rs` | `OptimizerBackend` (7 backends), `run_optimizer`, `IntervalMap`, stuck-at-bound probe | **yes** |
| Basin adapters | `smatrix/synthesis/optimizer/basin_backends.rs` | Basin LM-QR / TRF | **yes** |
| Parameter list | `evaluator.rs` `Param`, `build_params`, `apply_params` | `Row(i)` thickness or `Span` total thickness | thickness only |
| Bounds | `evaluator.rs:466-488` | one global `[0 or clamp_min, clamp_max]` | thickness only |
| Analytic Jacobian | `jacobian.rs`, `evaluator.rs` `DepositJacobian`, `simulate_with_deposits` | ∂(curve)/∂d via the needle operator, through the incoherent cascade since W6 | thickness only |
| Forward model | `SmatrixContext::simulate` | coherent blocks joined by the engine's intensity cascade (fixed incoherent plates, W6); R/T/A per polarization, front and back, complex r/t rows (phase only where no plate intervenes) | design stacks only: `LayerSpec` rows, no groups, materials or grading parameters as inputs |
| Merit | `merit.rs` `MeritSpec` | design targets: kinds e/a/b/r/c, transforms, point weights, colour, phase, differential phase, retardance, s/p ratios, Ψ, GD…FOD, Stokes/DOP for an input state (W8), multi-environment | 18 `CurveId` rows plus quantities derived from them in the merit |
| Thread pool sizing | `evaluator.rs` `in_lm_pool`, `lm_threads_for` | per-run rayon pool sized by work | **yes** |
| Post-run sweep | `DesignStack::clamp_all_policy` | removes/clamps thin films | thickness only |

So **only thicknesses** are reachable (per film, or one total per graded
span). The scoring side is now broad (every W8 target, plates, several
environments); the parameter side is not.

### 1.2 Declared but unreachable

- `Group.optimization_mask: [i32; 7]` with `OptMask { Thickness, N, K,
  Roughness, InhDelta, Interface, Material }` (`structure/enums.rs:114`).
  Only `Thickness` is read, by `Structure/Architect::optimization_entries`;
  synthesis ignores the mask completely and reads `LayerSpec.optimize`.
- Layer values that are inputs only: `inh_delta`, `InhMode::RateCapped
  {rate, ref_thickness, cap}`, `roughness`, `interface_thickness`,
  `GradientSpec` (`f_start`, `f_end`, `rate`, `f_min`, `f_max`).
- Group values: `thick_factor`, `thick_summand`, `n_factor`, `k_factor`,
  `inh_delta_summand`, `roughness_summand`, `interface_summand`.
- Material parameters of all 23 `MaterialSpec` models
  (`structure/specs.rs:20`), including EMA `fraction` and the nested
  `host` / `inclusion` specs.
- `materials/mod.rs` already names its `Dispersion` trait and `Model` enum
  "the seam for a future Rust-side fitting loop", and
  `docs/materials-architecture.md` §(b) / phase 6 defers that loop "until
  measurements justify it". This plan is that loop.

### 1.3 What a measurement fit needs that synthesis does not have

| Need | Synthesis path (0.9.0) | General path (`smatrix/solver.rs`) |
|---|---|---|
| Thick incoherent substrate (back-side reflection) | fixed plates only, intensity cascade (W6) | yes: `incoherent_flags`, `coherence_mode` |
| Ellipsometry Psi/Delta | Ψ and Δ (`PsiR/T`, `DeltaR/T`) as design targets, front, coherent for Δ | yes: `Psi_R/T`, `Delta_R/T` (`REQ_PSI_*`, `REQ_DELTA_*`) |
| Stokes, DOP, retardance, amplitudes, GD/GDD | front elements as design targets (W8); none through a plate when they read phase | yes (45 output maps), every coherence mode |
| Partial coherence modes (`COHERENCY_MATRIX`) | no | yes |
| Graded films, interfaces, roughness, groups | via `from_design` spans | via `structure::expansion::expand` |
| Monte-Carlo error draws | no | `ExpandOptions { apply_errors, seed }` |

**Decision: the fit engine's stack forward model is
`expand` → `solve_arrays`, not `SmatrixContext`.** The synthesis path stays
the design engine. W6-W8 narrowed the channel gap, not the reason for the
decision: a fit varies groups, grading, materials and EMA fractions, which
only `expand` turns into solver rows, and it needs every coherence mode and
the back elements, which the synthesis evaluator does not carry.

---

## 2. (a) Naming: split "solver" from "thickness refinement"

Three things are called "optimizer" today, and none of them is the general
engine:

1. the least-squares solvers and their dispatch (general, but live under
   `synthesis` with thickness names),
2. the thickness parameterisation inside synthesis (the specific path),
3. `smatrix::optimizer`, which is the **eigenmode** Nelder-Mead landscape
   search and has nothing to do with either.

### 2.1 New names

| Today | New | Why |
|---|---|---|
| `smatrix::synthesis::thick_opt` (LM, `LmConfig`, `JacobianSource`, `NoJacobian`, `LmResult`, `LmTermination`, `levenberg_marquardt[_with]`) | `crate::lsq::lm` | generic solver; shared by both paths |
| `smatrix::synthesis::trf` | `crate::lsq::trf` | generic |
| `smatrix::synthesis::optimizer` (`OptimizerBackend`, `OptimizerResult`, `run_optimizer`, `IntervalMap`, `available_backends`) | `crate::lsq::backends` with `LsqBackend`, `LsqResult`, `solve_least_squares`, `IntervalMap`, `available_backends` | "optimizer" becomes free for nothing; the backend is a least-squares backend |
| `optimizer/basin_backends.rs` | `crate::lsq::basin` | generic |
| `evaluator.rs` `in_lm_pool`, `lm_threads_for`, `LM_THREADS`, `LM_BASE_WORK` | `crate::lsq::pool` | both paths need it |
| `evaluator.rs` `Param`, `build_params`, `apply_params`, `DepositJacobian` | `smatrix::synthesis::thickness_refine` with `ThicknessParam`, `thickness_params`, `apply_thickness_params`, `ThicknessJacobian` | the specific path gets a specific name |
| `smatrix::synthesis::jacobian` | `smatrix::synthesis::thickness_jacobian` | it is ∂/∂d only |
| `SmatrixContext::optimize_thicknesses[_report]` | `refine_thicknesses[_report]` | matches the module; "optimize" stays for the pipeline as a whole |
| `smatrix::optimizer` (eigenmode NM) | `smatrix::eigenmodes` | removes the third "optimizer" |
| Python `LmConfig(optimizer=...)` | `LmConfig(backend=...)` | the keyword names a least-squares backend; `optimizer=` is simply an unknown keyword afterwards (decision 1) |
| Python `available_optimizers()` | `available_lsq_backends()` | same |
| Docs: `docs/basin-backends.md`, pipeline docs | "thickness refinement" for the stage, "LSQ backend" for the solver | |

`LmConfig` keeps its name: it is solver settings and both paths take it.

### 2.2 How

- **Hard rename (decision 1)**: every caller is updated in the same change;
  no alias, re-export, deprecation path or old-name refusal message is
  kept, in Rust or Python. Two commits for reviewability: (1) file moves
  with call sites updated, (2) identifier renames. Neither leaves a
  compatibility shim behind.
- Size: 15 Rust files reference the moved modules (`thick_opt` 25 hits,
  `synthesis::optimizer` 8, `super::optimizer` 6, `trf` 3, `jacobian` 4);
  8 Python/doc files use `optimizer=` or `available_optimizers`.
- `tools/check_exposure.py` allowlist entries move with the functions.
- Gate: bench fingerprint `46dc154d6946ee40` and both needle digests
  unchanged (pure rename); full battery green.

---

## 3. (b) The general fit engine

### 3.1 Module layout

```text
rust/navette/src/
  lsq/                    (moved, §2)   solvers + dispatch + pool
    global/de.rs          (phase 7)     hand-rolled DE, default build
  fit/
    mod.rs                FitProblem, fit(), FitReport
    params.rs             ParamAddr, Parameter, ties, transforms
    sets.rs               ParamSet, built-in sets, selector resolution
    convert.rs            InstrumentConvention, to_solver / from_solver, presets
    model.rs              FitModel: owns architect/structures, groups, material specs
    forward/
      stack.rs            StackForward: apply → expand → solve_arrays → channel maps
      dispersion.rs       DispersionForward: MaterialSpec.evaluate on data grids   (§4)
    data.rs               Dataset kinds, residual blocks, sigma/weights, masks
    residuals.rs          assembles all blocks into one r(x); Jacobian columns
    strategy.rs           Stage, Selector, Strategy, run_strategy
    explore.rs            multi-start, scans, strategy tournament, StrategyGen
    stats.rs              χ², covariance, correlations, AIC/BIC, uncertainties
    config.rs             FitRequest (serde), program section "fit"
rust/navette-py/src/fit.rs       Python door
src/navette/fit/                 Python wrapper, report helpers
```

`fit` depends on `lsq`, `structure`, `materials`, `smatrix::solver`,
`interpolate`, and optionally on `smatrix::synthesis::merit` (for
design-style targets, §3.5.3). Nothing in synthesis depends on `fit`.

### 3.2 The model the engine fits

`FitModel` is the container one evaluation reads:

- **Materials**: `HashMap<code, MaterialSpec>` (all 23 models; EMA specs
  nest their constituents). Tables stay `Entry::Array` / `Table` specs.
- **Structure**: an `Architect` (blocks, inversion, shared structures) or a
  single `Structure`, with its `Group`s.
- **Measurement configurations**: per dataset, the ambient / substrate /
  back-side context, angles, polarisation, `coherence_mode`, which output
  maps to request.

It is built from what already loads:

| Source | Reuse |
|---|---|
| Program file | `config::load_program_json` → `LoadedProgram { materials: SpecProvider, groups, structures, architect }` — **directly usable** as the model |
| Structure state | `Architect::from_state` / `Structure::from_state` (schema v3 migration included) |
| Design request | `DesignRequest` + `eval_library` / `spec_from_def` for the library |
| Python objects | `Navette_Architect`, `Navette_Structure`, `Group`, material specs via existing pyo3 classes |

### 3.3 Parameters

#### 3.3.1 Addresses

One `ParamAddr` names one scalar in the model. Proposed textual form (used
by configs, reports and selectors):

| Kind | Address | Written into |
|---|---|---|
| Thickness | `layer[main:3].thickness` or `layer[main:#TiO2_2].thickness` | `Layer.thickness` |
| Grading | `layer[...].inh_delta`, `layer[...].inh_rate`, `layer[...].inh_cap` | `Layer.inh_delta`, `InhMode::RateCapped` |
| Roughness / interface | `layer[...].roughness`, `layer[...].interface_thickness` | `Layer` |
| Gradient | `layer[...].gradient.f_start`, `.f_end`, `.rate`, `.f_min`, `.f_max` | `GradientSpec.mode` |
| Group | `group[TiO2].n_factor`, `.k_factor`, `.thick_factor`, `.thick_summand`, `.inh_delta_summand`, `.roughness_summand`, `.interface_summand` | `Group` |
| Material | `material[TiO2].A`, `material[TiO2].osc[1][0]`, `material[porous].fraction`, `material[porous].host.B1` | JSON pointer into `MaterialSpec.params`, recursive through `host` / `inclusion` / `bottom` / `top` |
| Instrument (phase 7) | `data[ell1].angle_offset`, `.wavelength_shift`, `.scale`, `.offset`, `.bandwidth` | the dataset's measurement config |

Layer addresses resolve through `Architect::map_global` / per-structure
indices; shared structures are one address (editing a shared `Structure`
reaches every block that points at it, which is the existing `Rc` sharing).

Material addresses reuse `MaterialSpec`'s plain-JSON params: the address is
a JSON pointer, written with `serde_json::Value::pointer_mut`. **No
per-model code is needed to make a model fittable.** Validation that the
pointer lands on a number happens once, when the problem is built.

#### 3.3.2 Parameter record

```rust
pub struct Parameter {
    pub addr: ParamAddr,
    pub value: f64,          // start, then current
    pub lb: f64, pub ub: f64,
    pub scale: f64,          // typical size; the solver sees (x - x0)/scale
    pub transform: Transform,// Linear | Log (positive params: alpha0, Γ, f, ...)
    pub fixed: bool,
    pub tags: Vec<String>,   // for selectors: "dispersion", "thickness", "H-layers", ...
}
```

- **Normalisation lives here, not in the solvers.** The FD step in
  `lsq::lm` is `h = ∛ε·max(|x|, 1)`, which is wrong by orders of
  magnitude for an Urbach `alpha0` of 1e3 or a Sellmeier `C1` of 0.01. The
  fit layer hands the solver normalised coordinates `u = (t(x) − t(x0)) /
  scale`, so every solver, backend and FD stencil works unchanged.
- **Default bounds and scales per kind**, from the existing gates:
  thickness `[0, clamp_max]`; `inh_delta` `(-2, 2)` (`Layer::validate`);
  `cap` `[0, 1]` (`InhMode::issues`); group factors `> 0`; EMA fraction
  `[0, 1]`; gradient fractions `[0, 1]` (`GradientSpec::issues`);
  material defaults from a small per-model table (positivity of
  oscillator widths and strengths, `eu > 0`, `lambda_g > 0`).
- **Validity after every write**: the existing `Layer::validate`,
  `Group::validate` and `GradientSpec::issues` run on the touched objects
  once per problem build; per evaluation only the cheap numeric checks are
  repeated (a negative k is a residual failure, not a panic).

#### 3.3.3 Ties

`Tie { driver: ParamAddr, targets: Vec<(ParamAddr, f64 /*factor*/, f64 /*offset*/)> }`

- Same thickness for all L layers; same material in two structures; the
  substrate of sample A equals that of sample B.
- Groups already express "all layers of material X scale together"
  (`thick_factor`, `n_factor`), so **prefer group parameters over ties**
  for that case; ties are for what groups cannot say.

#### 3.3.4 Selectors and the optimisation mask

Selectors pick the free set of a stage: by kind (`kind=thickness`), by tag,
by glob on the address (`material[TiO2].*`, `group[*].n_factor`), or
explicitly. `OptMask` becomes meaningful as a **default selector**: a
group's mask channel `N` set to 1 adds `group[g].n_factor`, `Thickness` adds
its layers' thicknesses, `Roughness` adds `roughness_summand`, `InhDelta`
adds `inh_delta_summand`, `Interface` adds `interface_summand`. `Material`
(discrete material swap) stays out of scope; the needle pool already covers
discrete material choice in synthesis.

#### 3.3.5 Parameter sets (decision 6)

A **parameter set** is a named, reusable group of parameters: the unit a
strategy orders, and the unit the tournament rearranges.

```rust
pub struct ParamSet {
    pub name: String,            // "dispersion_H", "thicknesses", "group_factors", ...
    pub members: Vec<Selector>,  // resolved once, at problem build, to addresses
}
```

- Built-in sets come from the model's structure, so a user gets sensible
  units without writing selectors: one per material (`material[TiO2]`),
  one per group (`group[H]`: its factors and summands), `thicknesses`,
  `grading`, `roughness`, `interfaces`, and one per parameter kind.
- User sets are any selector union (`"H_all" = material[TiO2].* +
  group[H].*`), so sets can be **regrouped**: split a material's
  oscillators into two sets, or merge a group's factors with its
  material's dispersion.
- A stage's `free` list may name sets (`set:dispersion_H`) as well as
  selectors. Stages written against sets stay valid when the stack
  changes, because sets resolve against the model.
- Sets may overlap; a stage frees the union. The report lists, per stage,
  which sets were free and the resolved addresses.

### 3.4 Forward model: `StackForward`

One evaluation, for one parameter vector:

1. **Write** parameters into a working copy of `FitModel` (apply ties).
2. **Materials**: re-evaluate only materials whose parameters changed.
   `SpecProvider` already memoises by name and has `invalidate(Some(name))`
   (`structure/providers.rs:264`) — exactly the dirty tracking needed.
3. **Expand**: `Architect::solver_inputs(provider, wavelengths, opts)` →
   `SolverArrays` + spans, deterministic (`ExpandOptions::deterministic()`).
   Groups, grading, gradients, interfaces, roughness and inversion all come
   for free; the schema-v3 sign is respected.
4. **Solve**: `smatrix::solver::solve_arrays(indices, thicknesses,
   incoherent, rough_types, rough_vals, wavelengths, angles, requested,
   coherence_mode)` → `Solution { f64maps, c64maps, dispmaps }`.
5. **Pick** the channels each dataset reads.

**Batching**: datasets that share (model variant, wavelengths, angles,
coherence mode) are solved once with their request bits OR-ed.
**Grids**: each solve runs on the union of its datasets' wavelengths, or on
the measured wavelengths directly when they coincide; otherwise the
simulated curve is resampled with `interpolate::UniInterpolator`
(`pchip`/`linear`), never extrapolated.

#### 3.4.1 Adaptations needed

| Item | Today | Change |
|---|---|---|
| `SpecProvider` cache | `RefCell` → not `Sync` | the residual closure must be `Sync` (column-parallel FD). Give each evaluation its own provider (cheap: clone entries, keep a per-worker cache), or move the cache to `Mutex`/per-thread. **Prefer per-evaluation working copies**; they also make the problem embarrassingly parallel across strategies. |
| `Architect` | `Rc<RefCell<Structure>>` → not `Send` | same answer: the fit owns a `Send` description (states / `Vec<(Layer, bool)>` traversal + groups) and builds the expansion input per evaluation, or a fit-local `Send` mirror of the traversal. Decide in phase 1 by measuring: `iter_entries()` + `merged_groups()` once per problem, then only parameter writes per evaluation. |
| Sample/measurement context | ambient/substrate are layers of the structure | each dataset may override ambient, substrate and back-side (several samples, one film recipe) — modelled as **environments**, mirroring the synthesis `EnvironmentCfg` vocabulary (§3.6). |
| Request bits | Python view helpers (`rt_request`, `ellipsometry_request`, …) | reuse `solver.rs` request helpers to map channel names → bits. |

### 3.5 Data and residuals

All blocks produce `(sim − target)/σ` rows; the solver minimises the sum of
squares.

#### 3.5.1 Measured spectra (new)

```rust
pub struct SpectrumData {
    pub id: String,
    pub env: String,             // which sample / measurement configuration
    pub channel: String,         // any Solution map key: Rs, Tp, R_avg, Psi_R, Delta_R, S1_R, ...
    pub angle: Radians,
    pub wavelengths: Vec<f64>,
    pub values: Vec<f64>,
    pub sigma: Vec<f64>,         // per point; or a scalar
    pub weight: f64,
    pub window: Option<(f64, f64)>, // nm; stages can override
}
```

- **Delta** residuals wrap: `wrap(sim − meas)` into (−180°, 180°].
- **Unit conventions are the solver's** (decision 3). `SpectrumData`
  always holds solver-convention values; anything an instrument wrote is
  converted on load (§3.5.5), never inside the residual.
- Optional instrument response: absolute scale/offset per dataset
  (fittable, phase 7), bandwidth averaging (phase 7).

#### 3.5.2 Dispersion data (new, §4)

`NkData { material, wavelengths, n, k, sigma_n, sigma_k }` and the ε₁/ε₂,
n-only, k-only variants.

#### 3.5.3 Design targets (reused)

A `MeritSpec` (compiled from a `TargetSet` by the existing
`compile_merit_spec` path) can be one residual block: adapt `Solution` →
`SimCurves` (keyed by `CurveId`, row-major `[angle, λ]`) and call
`MeritSpec::residuals`. This gives the general engine **design
optimisation of non-thickness parameters** for free: tune `inh_delta`, an
EMA fraction or a group factor against spectral, colour and phase targets
with every constraint kind (e/a/b/r/c), weight and count normalisation the
synthesis merit has, including the W8 types (s/p ratios, Ψ, retardance,
GD…FOD, Stokes/DOP for an input state). Multi-environment specs map env
names one to one.

#### 3.5.4 Priors and regularisation (new)

- Prior per parameter: `(p − p0)/σ_p`. Stabilises reverse engineering when
  thickness and index are correlated.
- Smoothness for point-wise or spline materials (phase 7): second
  differences of n(λ), k(λ).
- These are ordinary residual rows, so every solver handles them.

#### 3.5.5 Instrument data converter (decision 3)

Measured data enters through `fit::convert`, which maps an instrument's
conventions onto the solver's and records what it did:

```rust
pub struct InstrumentConvention {
    pub wavelength_unit: WavelengthUnit,  // nm | um | A | eV | cm-1  -> nm
    pub angle_unit: AngleUnit,            // deg | rad                -> solver angle (Radians)
    pub value_scale: ValueScale,          // fraction | percent       -> fraction
    pub psi_delta: PsiDeltaConvention,    // units, Delta range ([0,360) / (-180,180]),
                                          // sign (Nebraska vs. exp(+iwt) / exp(-iwt) convention)
    pub ellipsometric_form: EllForm,      // Psi/Delta | tan(Psi), cos(Delta) | N, C, S -> Psi/Delta
    pub sigma: SigmaSource,               // column | constant | relative | none (-> unit weight)
}

pub fn to_solver(raw: RawDataset, conv: &InstrumentConvention)
    -> Result<SpectrumData, String>;
```

- **Presets** for the conventions common instruments and file types use,
  kept as data (`conventions/*.json`), plus a fully explicit form; no
  vendor file parsers in the core (reading a CSV/text file is the Python
  door's job, or a later reader module).
- **Pinned by round-trip tests**: solver output → `from_solver` (the
  inverse) → `to_solver` reproduces the solver output for every preset,
  and a known film's Psi/Delta in each convention converts to the same
  solver values.
- σ is propagated through the conversion (e.g. N, C, S → Psi/Delta uses
  the Jacobian of the map).
- The convention used is stored in the `FitReport`, so a result says how
  its data was interpreted.

### 3.6 Several samples, one model

Reverse engineering usually combines measurements: R and T of the coating,
the bare substrate, a witness single layer, Psi/Delta at two angles. The
engine fits **one parameter vector against all of them**:

- an environment = one sample configuration (ambient, substrate,
  back-side, which structure / blocks), reusing the naming and validation
  ideas of `smatrix/synthesis/environments.rs` (`EnvironmentCfg`,
  segment references);
- every dataset names its environment; materials and group parameters are
  shared across environments by construction (one `FitModel`);
- residual ordering is env-major, like `MeritSpec::residuals_multi`.

### 3.7 Jacobian

| Columns | Phase 1 | Later |
|---|---|---|
| any parameter | central differences in normalised coordinates, column-parallel (existing `lsq::lm` FD path; pool sizing from `lsq::pool`) | — |
| dispersion parameters on nk data | FD | analytic ∂n/∂p for closed forms (Cauchy, Sellmeier, Urbach, Lorentz, Drude) — trivial derivatives, phase 6 |
| thickness on coherent R/T/phase datasets | FD | reuse the needle-operator deposits (`simulate_with_deposits`) by running those datasets through `SmatrixContext` when they are fully coherent; or port the deposit read-off onto `solve_arrays` (phase 6) |
| index of one layer (→ material params, group n/k factors) | FD | first-order perturbation: ∂r/∂ε of a layer is the depth integral of the needle function over that layer; chain with ∂ε/∂p. Research item, phase 6+ |

**Sparsity worth exploiting early**: a material parameter only changes the
datasets whose environments contain that material; a dispersion-data block
depends only on its material's parameters. Phase 1 evaluates all blocks
per column (simple, correct); phase 6 skips unaffected blocks per column.

### 3.8 Strategies (staged fitting)

A strategy is an ordered list of stages. Each stage says **what is free,
what is scored, and how to solve**.

```rust
pub struct Stage {
    pub name: String,
    pub free: Vec<Selector>,           // everything else is fixed at its current value
    pub datasets: Vec<DatasetSel>,     // which blocks are active, with weight overrides
    pub window: Option<(f64, f64)>,    // wavelength window override
    pub kind: StageKind,
    pub lsq: LmConfig,                 // backend, tolerances, budgets (existing struct)
    pub repeat: Repeat,                // Once | UntilNoGain { rel, max }
    pub accept: Accept,                // Always | IfBetter (else roll back to the stage's start)
}

pub enum StageKind {
    Lsq,                                           // local least squares (any lsq backend)
    Scan { param: Selector, points: usize },       // 1-D grid, keep the best, then Lsq
    MultiStart { starts: usize, sampler: Sampler, seed: u64 }, // LHS / uniform in bounds, Lsq each
    Global { method: GlobalMethod, budget: usize, seed: u64 }, // §3.8.3
}
```

Examples a user can now express directly:

1. **Materials first**: free `material[*].*`, datasets = witness-film nk
   data and single-layer T; then free thicknesses; then everything.
2. **Thickness first**: free `kind=thickness` on R/T in a window where the
   materials are transparent; then `group[*].n_factor`, `k_factor`; then
   dispersion parameters; then release all for a polish.
3. **Group-factor model of deposition errors**: free `group[*].thick_factor`
   and `n_factor` only (the classic "systematic error" reverse-engineering
   model), compare χ² against strategy 2.
4. **Escape a local minimum**: `Scan` on the thickest layer, then `Lsq`.

Every stage starts from the previous stage's end state. Stages are pure
data, so a strategy round-trips through JSON/YAML.

#### 3.8.1 Comparing strategies

`explore::tournament(problem, strategies, seeds)` runs several strategies
(and/or seeds) on the same problem and ranks them by final χ², reduced χ²,
parameter uncertainty, run time and evaluations. Runs are independent, so
they parallelise across strategies; `lsq::pool` already yields to an outer
rayon worker (`in_lm_pool` runs inline on a worker), so nesting is safe.
Results are deterministic for a fixed seed (ordered collection, the rule
the solver already follows).

**Rearranging parameter sets (decision 6).** Besides hand-written
strategies, the tournament can generate strategies from parameter sets
(§3.3.5):

```rust
pub enum StrategyGen {
    Explicit(Vec<Strategy>),
    /// Every order of the given sets (or a seeded sample when n! is large),
    /// one Lsq stage per set, then a final stage with all sets free.
    Orderings { sets: Vec<String>, sample: Option<(usize, u64)>, polish: bool },
    /// Fixed order of "blocks", each block a list of sets freed together:
    /// [[dispersion_H, dispersion_L], [thicknesses], [group_factors]].
    /// Variants regroup the same sets into different blocks.
    Groupings { blocks: Vec<Vec<Vec<String>>>, polish: bool },
    /// Cumulative: free set 1, then 1+2, then 1+2+3 ... in each order.
    Cumulative { sets: Vec<String>, sample: Option<(usize, u64)> },
}
```

- Each generated strategy is named after its arrangement
  (`dispersion_H > thicknesses > group_H | polish`), so the ranking table
  reads directly as "which order works".
- Generated and explicit strategies run in the same tournament and are
  scored the same way (χ², recovery error on synthetic truths, cost).
- `n!` grows fast: above a limit (default 24 orderings) the generator
  samples orderings with a seed, and says so in the report.

#### 3.8.2 Benchmarking strategies against a known truth (reuse)

The structure layer can already draw **realistic deposition errors**:
`Group` error params per channel and `ExpandOptions { apply_errors: true,
seed }`. That gives a ground-truth generator for free:

1. take a design, draw K error realisations (seeds),
2. simulate their "measurements" with the forward model, add noise,
3. run every strategy on every realisation,
4. score recovery error of the true parameters, not just χ².

This answers "which strategy works best for this kind of coating" with
numbers, which is the point of (b).

#### 3.8.3 Global stage (decision 5)

Two implementations of one `GlobalMethod` stage, in this order:

1. **Reference, behind `opt-basin`** (phase 6). basin 1.15, already a
   dependency of that feature, provides:
   - `De`: differential evolution over a box, SciPy-named mutation rules
     (`Rand1`, `Best1`, `Rand2`, `Best2`, `RandToBest1`, `CurrentToBest1`),
     binomial / exponential crossover, dithering, seeded, deterministic
     batch evaluation under `parallel`;
   - `DeInject`: memetic DE, top-k candidates refined by a local solver
     every generation (the natural "DE then LM" hybrid);
   - `BoundedCmaEs` (and `CmaInject`), `GlobalBestPso`, `BasinHopping`,
     `Direct`, `SimulatedAnnealing` as further references.
   The stage maps the normalised fit box straight onto basin's box and
   uses the χ² of the residual vector as the cost.
2. **Hand-rolled DE in the default build** (phase 7, `lsq::global::de`):
   `DE/rand/1/bin` and `DE/best/1/bin`, dithering, reinit-per-coordinate
   bound repair (basin's choice), synchronous generations evaluated in
   parallel with ordered collection, an optional LM polish of the best
   member. No new dependency (uses the crate's existing RNG).
   - **Checked against the reference**: same benchmark functions and the
     phase-4 fit problems, compared on success rate and evaluations to a
     given χ² over many seeds (the trajectories will differ; the
     statistics must not be worse).
   - The basin variants stay available under `opt-basin` for comparison,
     the way the LSQ backends are today.

### 3.9 Statistics and reporting

Per stage and at the end (`FitReport`):

- χ², reduced χ² (`χ²/(m − n_free)`), RMS per dataset;
- covariance `C = (JᵀJ)⁻¹ · χ²_red` at the solution (QR of J, reusing the
  LM's factorisation code), standard errors, **correlation matrix** (the
  thickness ↔ index degeneracy is the main failure mode of reverse
  engineering, and it shows here first);
- parameters at bounds (and whether the bound is active, reusing
  `flag_stuck_at_bound`'s probe idea);
- `LmTermination` per stage (existing), evaluations, Jacobian kind;
- AIC/BIC for model comparison (Cauchy vs Sellmeier vs Tauc-Lorentz);
- the final model as **state** (`Architect::to_state`, material specs as
  `MaterialDef`), so a fit result loads anywhere a design loads.

### 3.10 Configuration and doors

#### 3.10.1 Rust request

```jsonc
{
  "model": { /* program sections: materials, groups, structures, architect */ },
  "environments": [ {"name": "coating", "ambient": "air", "substrate": "BK7",
                     "backside": "air", "coherence_mode": 1}, ... ],
  "data":   [ {"id": "R8", "env": "coating", "channel": "Rs", "angle_deg": 8,
               "wavelengths": [...], "values": [...], "sigma": 0.002} ],
  "parameters": [ {"addr": "material[TiO2].A", "lb": 1.8, "ub": 2.6, "scale": 0.05},
                  {"select": "kind=thickness", "lb": 0, "ub": 400} ],
  "ties": [ ... ],
  "strategy": { "stages": [ ... ] }
}
```

- `FitRequest` with `#[serde(deny_unknown_fields)]`, angle keys through the
  existing `resolve_angle_keys` convention (`angle_deg` / `angle_rad`).
- Program files gain a `fit` section: add `"fit"` to `PROGRAM_SECTIONS`
  (`config.rs:43`) and bump `PROGRAM_SCHEMA_VERSION` (4), which the section
  whitelist makes a loud refusal on older builds — the mechanism F2.4 built
  for exactly this.

#### 3.10.2 Python

`navette.fit`: `Parameter`, `Dataset`, `Stage`, `Strategy`, `FitProblem`,
`fit(problem, strategy, progress=callback) -> FitReport`, `tournament(...)`.
Data in / data out like `synthesis_pipeline.rs`; progress callback per stage
and per iteration like `run_design`'s callback. Plot helpers stay in Python.

---

## 4. (c) Dispersion-model fitting

### 4.1 What exists

- 23 models behind `MaterialSpec::evaluate` (`structure/specs.rs`):
  Konstant, Table, Cauchy, CauchyUrbach, Sellmeier, SellmeierUrbach,
  Lorentz, Drude, DrudeLorentz, CodyLorentz, Forouhi-Bloomer (4 variants),
  TaucLorentz, UBF, and the EMA mixers Bruggeman, MaxwellGarnett,
  Looyenga, Lichtenecker, MoriTanaka, PowerLaw, plus Roughness.
- EMA specs **nest** (`host`, `inclusion`, `bottom`, `top` are full specs),
  so a porous TiO2 with a Sellmeier host and air inclusion is one tree.
- Kernels in `materials/` are pure, rayon-parallel and bit-reproducible.
- `materials/kk.rs` (KK transform, energy grid) exists as a scaffold.
- Python docs show fitting done in scipy (`docs/materials-dispersion-usage.md`).

### 4.2 What is missing

- a parameter layer over `MaterialSpec` (solved by §3.3: JSON pointers),
- data containers for nk / ε data with σ,
- a forward model on arbitrary data grids,
- model-specific default bounds and scales,
- statistics and model comparison,
- analytic derivatives (performance only),
- physical-consistency checks (KK, k ≥ 0, no poles inside the fitted range
  for Sellmeier `C` terms).

### 4.3 Design: same engine, second forward model

Dispersion fitting is **not** a separate optimizer. It is the general
engine with `DispersionForward` instead of (or next to) `StackForward`:

```rust
pub struct DispersionForward;   // evaluates MaterialSpec on each NkData grid
```

- Residuals: `(n_sim − n)/σ_n`, `(k_sim − k)/σ_k`, or ε₁/ε₂ forms.
- Parameters: `material[...]` addresses only; the same `Parameter`,
  bounds, transforms, ties, stages and reports.
- Cost: microseconds per evaluation, so FD Jacobians are fine; multi-start
  is cheap and should be the default stage for oscillator models.

Because both forward models feed the same residual vector, **joint fits**
come for free: fit a film's Sellmeier parameters to its measured T
spectrum *and* to a literature nk table (as a weighted prior block), or fit
an EMA fraction against both a witness nk table and the coating's R.

### 4.4 EMA specifics

- Fittable: `fraction` (`[0, 1]`), `L` (Mori-Tanaka depolarisation),
  `alpha` (power law), and any parameter of the nested constituents.
- Bruggeman uses Newton per point with `max_iter`/`tol` fixed in
  `specs.rs` (100 / 1e-9). Its result is smooth in `fraction` except where
  Newton stops early; FD steps must be larger than `tol`'s effect.
  Phase 2 adds a test that the FD Jacobian of a Bruggeman fraction is stable
  across step sizes.
- Graded EMA layers (`GradientSpec` with `ema`) are stack parameters
  (`gradient.f_start`, `.f_end`), already addressable in §3.3.

### 4.5 Convenience doors

- `fit_dispersion(spec, data, strategy=None) -> (fitted_spec, report)` in
  Rust and Python: builds a one-material `FitProblem`, default strategy
  multi-start + LSQ.
- `table_to_model(table_material, candidates=[...])`: fit each candidate
  model to a `Table` material and rank by AIC/BIC; returns the best spec
  and the ranking.
- Library write-back: a fitted spec replaces the library entry
  (`SpecProvider::upsert`) and serialises as a `MaterialDef`.

### 4.6 Deferred model option: point-wise and spline nk (decision 4)

For unknown materials, where no closed-form model fits, a free-form model:

- **`PointWise`**: n and k (or ε₁, ε₂) at every data wavelength are
  parameters. Only meaningful with several independent measurements per
  wavelength (R + T, several angles, Psi/Delta) and a smoothness prior.
- **`Spline`**: n, k on a coarse knot grid, interpolated with
  `interpolate::UniInterpolator` (pchip / makima); far fewer parameters,
  smooth by construction.
- **`BSplineEps`**: ε₂ as a B-spline with ε₁ from Kramers-Kronig (on
  `materials/kk.rs`, which needs finishing first), so the model is
  KK-consistent by construction.

What it needs from the engine: vector-valued material parameters
(`material[X].knots.n[7]`, addressable as JSON pointers already), the
smoothness residual block (§3.5.4), k ≥ 0 bounds per knot, and the
per-column block skipping of phase 6 (hundreds of parameters).
**Deferred**: not in phases 1-6; it becomes a new `MaterialSpec` model in
phase 7 and stays out of the default strategies until then.

---

## 5. Reuse map

### 5.1 Directly usable, unchanged

| Component | File | Used for |
|---|---|---|
| LM, TRF, all LSQ backends, `IntervalMap`, stuck-at-bound probe | `lsq/*` (after §2) | every `Lsq` stage |
| `LmConfig`, `LmTermination`, `LsqResult` | `lsq/lm.rs` | stage solver settings, reports |
| `in_lm_pool` / `lm_threads_for` | `lsq/pool.rs` | pool per fit run; nesting under tournaments |
| `expand`, `ExpandOptions`, `SolverArrays`, `Span` | `structure/expansion.rs` | stack forward model; truth generation with error draws |
| `Architect`, `Structure`, `Layer`, `Group`, `GradientSpec`, `InhMode` | `structure/*` | the fitted model; state round trip; validation gates |
| `MaterialSpec::evaluate`, all kernels, EMA | `structure/specs.rs`, `materials/*` | material forward model |
| `SpecProvider` memo + `invalidate` | `structure/providers.rs` | dirty tracking of materials (per evaluation copy, §3.4.1) |
| `solve_arrays`, `Solver`, request helpers, `Solution` maps | `smatrix/solver.rs` | all measurement channels incl. incoherent and ellipsometry |
| `UniInterpolator` | `interpolate/mod.rs` | resampling simulated curves onto measured grids |
| Program loader, version gate, section whitelist | `config.rs` | model input; `fit` section |
| Schema versioning + migration | `structure/version.rs` | fit results as states |
| `Radians`, angle-key resolution | `smatrix/angle.rs`, `targets.rs` | dataset angles |

### 5.2 Usable with an adapter

| Component | Adapter |
|---|---|
| `MeritSpec` / `TargetSet` / `compile_merit_spec` | `Solution` → `SimCurves` conversion so design targets become a residual block (§3.5.3) |
| `environments.rs` vocabulary (`EnvironmentCfg`, rosters) | fit environments reuse names and validation rules; not the `CompiledEnvironments` routing (that is design-segment specific) |
| `optimization_mask` / `OptMask` | default selectors (§3.3.4) |
| Thickness deposits (`simulate_with_deposits`, `thickness_jacobian`) | analytic thickness columns for coherent R/T/phase datasets (phase 6) |
| LM's QR factorisation | covariance at the solution (§3.9) |
| `pyo3` patterns in `navette-py/src/synthesis_pipeline.rs` | template for `navette-py/src/fit.rs` |
| Golden forward baseline (`validation/regression/forward`) | its cases double as synthetic fit problems (known truth) |

### 5.3 New

- `fit::params` (addresses, normalisation, ties, selectors),
- `fit::data` (spectra, nk data, priors, Delta wrapping, σ handling),
- `fit::forward::{stack, dispersion}`,
- `fit::strategy`, `fit::explore` (multi-start, scans, tournament;
  differential evolution in phase 6),
- `fit::stats` (covariance, correlation, AIC/BIC),
- `FitRequest`, program `fit` section, Python door.

### 5.4 Deliberately not reused

- `SmatrixContext` / `DesignStack` as the fit forward model: its inputs
  are design rows (`LayerSpec`), not the model a fit varies (groups,
  grading, materials, EMA); it has fixed plates but no coherence modes;
  and its mutation operators (needles, merge, clamp sweeps) are design
  operations, not fit operations. Its *merit* is reused (§3.5.3), and
  since W8 that brings every design target type with it.
- The thin-film clamp/remove sweep: a fit must not delete layers; bounds do
  the job. (A structural "drop a layer and refit" step can be a later
  explore stage if wanted.)

---

## 6. Phases

Each phase ends with the full battery green and the bench fingerprint
`46dc154d6946ee40` unchanged (the fit engine never touches the forward
model or the synthesis path).

| Phase | Scope | Done when |
|---|---|---|
| **0** | §2 rename: `lsq` module, `thickness_refine`, `thickness_jacobian`, `eigenmodes`, Python `backend=` / `available_lsq_backends()` | pure moves; fingerprint and needle digests unchanged; CHANGELOG "Changed (breaking)" under 0.10.0 |
| **1** | `fit::params` (thickness, grading, roughness, interface, group, material JSON-pointer addresses; bounds, scale, log transform), `StackForward` on `expand` + `solve_arrays`, `SpectrumData` incl. Psi/Delta, single `Lsq` stage, FD Jacobian, χ² + covariance report | recovery tests: a 3-layer stack with known truth + noise is recovered within 2σ for thickness, `n_factor`, and a Cauchy `A`; Psi/Delta fit of a SiO2 film on Si |
| **2** | `DispersionForward`, `NkData` (n/k, ε₁/ε₂, n-only), EMA nested addresses, default bounds per model, `fit_dispersion`, `table_to_model`, AIC/BIC | Sellmeier/Cauchy/Tauc-Lorentz recover their own synthetic data; Bruggeman fraction recovered; FD-step stability test for Bruggeman |
| **3** | `Strategy` / `Stage` / selectors / parameter sets / `Repeat` / `Accept`, `Scan`, `MultiStart`, ties, priors, environments (several samples), instrument converter with presets | the four example strategies in §3.8 run from JSON; multi-sample fit (coating + bare substrate) shares the substrate parameters; converter round trips pass for every preset |
| **4** | strategy tournament incl. `StrategyGen` (orderings, groupings, cumulative), truth generation from group error draws (§3.8.2), recovery scoring | a benchmark script ranks explicit and generated strategies on K seeds of a 20-layer coating; deterministic per seed |
| **5** | Python door (`navette.fit`), program `fit` section (program schema v4), docs (`docs/fitting.md`), examples | examples run in the battery; `.pyi` stubs; `check_exposure` clean |
| **6** | performance: per-column block skipping, analytic dispersion derivatives, analytic thickness columns on coherent datasets, index-perturbation research; `Global` stage on basin (`De`, `DeInject`, `BoundedCmaEs`) behind `opt-basin` | measured speed-ups recorded in a bench like `bench_lm_breakdown.py`; basin DE solves the phase-4 problems from random starts |
| **7** | hand-rolled DE in the default build, checked against the basin reference; instrument parameters (angle offset, λ shift, scale/offset, bandwidth); deferred model option: point-wise / spline / KK-B-spline nk (§4.6) | DE statistics not worse than basin's over many seeds; each other item with a recovery test |

Phases 1 and 2 are independent after phase 0 and can run in parallel
(different files: `fit::forward::stack` vs `fit::forward::dispersion`);
phase 3 needs both.

---

## 7. Tests and validation

- **Recovery tests** (the core contract): truth → forward → noise → fit →
  parameters within kσ; per parameter kind and per channel kind.
- **FD Jacobian checks**: normalised-coordinate FD vs a fine-step oracle,
  per parameter kind; analytic columns (later) vs FD (rel. dev. ~1e-8, as
  the thickness Jacobian tests do).
- **No-op invariants**: fitting with every parameter fixed returns the
  model unchanged, bit for bit; a fit started at the truth with no noise
  stays there.
- **Determinism**: same request + seed → bit-identical report, across
  thread counts (`NAVETTE_LM_THREADS`), as the synthesis path guarantees.
- **State round trip**: fitted model → state → load → forward equals the
  fit's final forward bit for bit.
- **Parity**: a small problem against `scipy.optimize.least_squares` (the
  `test_lm_parity.py` pattern).
- **Guard rails**: invalid addresses, pointers to non-numbers, bounds
  outside the layer/group gates, Delta wrapping at ±180°, σ ≤ 0 — each
  refused with a message naming the address or dataset.

---

## 8. Risks

| Risk | Mitigation |
|---|---|
| `Rc`/`RefCell` model types are not `Send`/`Sync` | per-evaluation working copies; decide representation in phase 1 by measurement (§3.4.1) |
| FD step sizes across parameters spanning 1e-3 … 1e4 | normalised coordinates with per-parameter `scale`, log transform for positive scale parameters |
| Thickness–index degeneracy | correlation matrix in every report; priors; staged strategies; multi-sample fits |
| Bruggeman Newton tolerance makes FD noisy | step-size stability test; larger relative step for EMA fractions |
| Cost of expand + solve per column for big stacks | dataset batching, block skipping, analytic columns (phase 6) |
| Two forward models drifting (synthesis vs fit) | none needed: the fit uses the general solver, which the forward golden already pins |

## 9. Open questions

The five questions of the first draft are settled (§0). Still open, to
decide when the phase gets there:

1. **Converter presets** (phase 3): which instruments / file conventions
   to ship presets for first.
2. ~~**Generated-strategy budget** (phase 4)~~: settled in phase 4, 24
   orders before sampling, seeds left to the caller (see the status notes).
