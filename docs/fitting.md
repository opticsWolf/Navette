# Fitting

`navette.fit` fits any scalar of a model to measured data: film
thicknesses, group factors (index, extinction, thickness), and the numbers
of any dispersion model. A spectrum fit and a dispersion fit use one engine,
and both can run in a single residual vector. The engine is Rust
(`navette::fit`). The Python module builds the request and reads the reply.

Two runnable examples:

- `examples/fit_coating.py` reverse-engineers a 6-layer coating. It uses a
  staged strategy, a tournament, a benchmark and progress, and shows how to
  stop a fit.
- `examples/fit_optical_constants.py` covers a dispersion model for a
  measured n table, and a film's thickness and Cauchy index taken from an
  instrument file in percent.

## 1. A problem

```python
from navette.fit import Dataset, Environment, FitProblem, Parameter, fit

problem = FitProblem(
    model=sections,       # a program's sections: materials, groups, structures, ...
    environments=[Environment("sample", ambient="air", substrate="glass",
                              stack=[{"structure": "ar"}])],
    data=[Dataset("R", "R_avg", wl, measured, sigma=0.002, angle_deg=8.0)],
    parameters=[Parameter(select="kind=thickness"),
                Parameter("group[H].n_factor", lb=0.95, ub=1.05)],
)
report = fit(problem, progress=print)
print(report.summary())
```

| Piece | Meaning |
|---|---|
| `model` | A program's `sections`. The model is the program's architect, or its single structure, plus its material library and groups. With no structures it is a dispersion-only fit. |
| `Environment` | A measured sample: a chain of `stack` items between an optional `ambient` and `substrate` half-space. A stack item is `{"structure": s, "inverted": False, "repeat": 1}` or `{"layer": material, "thickness": d}`. Program structures are films only, so a sample is normally declared here. |
| `Dataset` | A spectrum in solver convention: `channel`, `wavelengths` in nm, `values`, and `sigma` (one number or one per point). Optional keys are `angle_deg`, `env`, `weight` and `window=(lo, hi)`. `env` defaults to the first declared environment. |
| `NkDataset` | Optical constants of one `material`. `columns` maps `n`, `k`, `eps1` or `eps2` to `(values, sigma)`. |
| `Parameter` | One scalar (`addr`) or a selector (`select`). Optional keys are `value`, `lb`, `ub`, `scale`, `log`, `fixed`, `tags` and `prior=(p0, sigma)`. |
| `Tie` | A follower that moves with its driver: `target = factor * driver + offset`. A follower is never a parameter itself. |
| `sets` | Named unions of selectors, the units a strategy orders. |
| `lsq` | Solver settings: `backend` (`"lm"`, `"trf"`), `max_iterations`, `max_evals`, `ftol`, `xtol`, `gtol`. |

Channels are the synthesis channel names:

- reflection and transmission: `Rs`, `Rp`, `Ts`, `Tp`, `R_avg`, `T_avg`
- absorption: `A_s`, `A_p`, `A_avg`
- ellipsometry: `Psi_R`, `Psi_T`, `Delta_R`, `Delta_T`
- polarization: `DOP_R`, `DOP_T`, `Diattenuation_R`, `Diattenuation_T`, `Retardance_R`, `Retardance_T`
- phase: `phi_rs` and the rest of the `phi_*` channels, and `PDts`/`PDtp`

An unknown name lists the known ones.

### Addresses and selectors

An address names one scalar:

| Address | Scalar |
|---|---|
| `layer[S:i].thickness` | the thickness of film row `i` (from 0) of structure `S`. A shared or repeated structure is one parameter. |
| `group[G].n_factor`, `k_factor`, `thick_factor`, `thick_summand` | a group's factors |
| `material[M].A` | a top-level number of material `M`'s model |

A selector picks several:

| Selector | Picks |
|---|---|
| `*` | everything |
| `kind=thickness` | by kind: `thickness`, `grading`, `roughness`, `interface`, `gradient`, `group`, `material` |
| `tag=H` | parameters carrying the tag |
| `set:name` | the members of a set |
| anything else | a glob on the address: `material[TiO2].*`, `group[*].n_factor` |

A `select` declaration expands over every address the model has. It skips
addresses that are declared explicitly elsewhere in the list, so an
explicit declaration overrides what a selector would give.

The built-in sets need no declaration:

- one per kind: `thicknesses`, `grading`, `roughness`, `interfaces`,
  `gradients`, `groups`, `materials`;
- `material[X]` and `group[G]` for every material and group a parameter
  names.

## 2. The report

`fit` and `evaluate` return a `FitReport` with these fields:

- the fit quality: `chi2`, `chi2_reduced`, `aic`, `bic`;
- the size of the problem: `n_residuals`, `n_free`;
- one `FittedParameter` per parameter (`start`, `value`, `std_error`,
  bounds, `at_bound`);
- `covariance` and `correlation`;
- how the solve ended: `termination`, `iterations`, `evals`, `backend`;
- the `curves`: measured, simulated and sigma per spectrum and per nk
  column.

Use `report.parameter(addr)` for one parameter, `report.values()` for all
values, and `report.summary()` for a printable table. `plot_fit(report)`
draws the curves (it needs matplotlib).

`evaluate(problem)` returns the same report at the start values without
moving anything. It is the model's prediction for your data.

Strongly correlated parameters show up as large error bars and
off-diagonal correlations near ±1. The fitted curve can still be right
when the individual numbers are not; the Sellmeier part of
`fit_optical_constants.py` shows this.

## 3. Strategies

A fit runs in one stage named `fit` with everything free, or by a
`Strategy`. A strategy is a list of `Stage`s, and each stage frees the
selectors in `free` and holds everything else at its current value:

```python
from navette.fit import Stage, Strategy

strategy = Strategy("thick first", [
    Stage("thicknesses", free=["kind=thickness"], window=(500, 900)),
    Stage("all", free=["*"], accept="if_better"),
])
report = fit(problem, strategy)
```

| Stage key | Meaning |
|---|---|
| `datasets`, `weights`, `window` | Restrict or reweight the data for this stage. |
| `kind` | `"lsq"` (the default). `{"scan": {"param": sel, "points": n}}` tries `n` values of one parameter across its bounds, then runs LSQ. `{"multi_start": {"starts": n, "sampler": "lhs", "seed": 0}}` runs LSQ from `n` sampled points and keeps the best. |
| `lsq` | Solver settings for this stage. |
| `repeat` | `"once"`, or `{"until_no_gain": {"rel": 1e-3, "max": 5}}`. |
| `accept` | `"always"` or `"if_better"`. |

`report.stages` holds each stage's χ², runs, evals and whether it was
accepted.

### Progress and stopping

`fit(problem, progress=cb)` calls `cb(event)` with one of three events:

- `{"event": "stage_start", "stage", "name", "free"}` when a stage begins;
- `{"event": "improved", "stage", "evals", "chi2"}` each time χ² drops
  within a stage;
- `{"event": "stage_end", "stage", "name", "chi2", "accepted"}` when a
  stage ends.

The callback may be called from several threads during a multi-start
stage. If it raises, the fit stops and comes back as a `ValueError` that
names the exception. This is how a GUI cancels a fit.

### Comparing strategies

Generators describe families of strategies:

- `orderings(sets)`: one stage per set, in every order, then a polish stage.
- `groupings(variants)`: each variant is a list of blocks freed together.
- `cumulative(sets)`: free the first set, then the first two, and so on.

`generate([...])` expands them and returns `(strategies, notes)`. A family
with more orders than `cap` (default 24) is sampled with a seed, and the
notes say so. Pass `sample=(count, seed)` to choose the sample yourself.

- `tournament(problem, strategies)` runs every strategy on the problem's
  data and ranks them by χ². It returns `(entries, notes)`.
- `benchmark(problem, strategies, seeds)` checks whether strategies find
  the right answer. For each seed it calls `draw_truth`, which draws a
  "true" model from the program's error laws (group `error_mask`, std
  devs), simulates its data, fits from the nominal start, and scores how
  well the truth is recovered. Its `summary()` ranks the strategies.

## 4. Instrument data

Instruments write data in their own conventions: units, percent, Δ range
and sign, and tan Ψ/cos Δ forms. `to_solver(raw, convention)` converts one
instrument record into `Dataset`s. `from_solver(id, data, convention)`
converts them back, so the round trip gives the same file.

A convention is a preset name from `conventions()`, such as `"nm_deg"`,
`"nm_deg_percent"`, `"ev_deg"` or `"nm_deg_tan_cos"`. It can also be a
dict, typically a preset with σ set:

```python
conv = {**conventions()["nm_deg_percent"], "sigma": {"constant": 0.2}}
data = to_solver({"id": "R", "quantity": {"spectrum": "R_avg"}, "angle": 8.0,
                  "abscissa": wl, "columns": [r_percent]}, conv)
```

σ is given in the instrument's units and is converted with the values. It
defaults to unit weight.

## 5. Dispersion models

- `fit_dispersion(spec, nk_data, parameters=None, lsq=None)` fits one
  material spec to `NkDataset`s. It returns `(fitted_spec, report)`.
  Without `parameters`, every number of the spec is free, except a
  `Table`'s data and a mixer's constituents.
- `table_to_model(table, candidates)` fits each candidate model to a
  `Table` spec and ranks them by BIC, best first. BIC here uses the
  Gaussian form that estimates the noise scale, so the ranking does not
  depend on an assumed σ.

To fit a film's dispersion together with its thickness, put the model's
numbers (`material[M].*`) and the thickness in one `FitProblem`.
`NkDataset`s and spectra can share that problem.

## 6. Saved fits: the program `fit` section

A `FitProblem` is a program document (schema version 4) with a `fit`
section. `problem.save(path)` and `FitProblem.load(path)` write and read
it, and `problem.to_document()` gives the dict. A program loader also
accepts the section, and `LoadedProgram.fit` holds it as a dict. The
section's keys mirror the dataclasses:

```json
{"schema_version": 4, "kind": "program", "sections": {
  "materials": [...], "groups": [...], "structures": [...],
  "fit": {
    "environments": [{"name": "sample", "ambient": "air", "substrate": "glass",
                      "stack": [{"structure": "ar"}]}],
    "data": [{"id": "R", "channel": "R_avg", "angle_deg": 8,
              "wavelengths": [...], "values": [...], "sigma": 0.002}],
    "nk_data": [{"id": "e", "material": "TiO2", "wavelengths": [...],
                 "columns": [{"quantity": "n", "values": [...], "sigma": 0.01}]}],
    "parameters": [{"select": "kind=thickness"},
                   {"addr": "group[H].n_factor", "lb": 0.95, "ub": 1.05}],
    "ties": [{"driver": "layer[a:1].thickness",
              "targets": [{"target": "layer[b:1].thickness", "factor": 1, "offset": 0}]}],
    "sets": [{"name": "H", "members": ["tag=H"]}],
    "strategy": {"name": "thick first", "stages": [...]},
    "coherence_mode": 0,
    "lsq": {"backend": "trf"}
  }}}
```

Floats are written and read back bit-exactly, so a saved fit reproduces
its χ².

## 7. From Rust

The engine has these entry points:

- `navette::fit::fit`, `run_strategy` and `fit_program`, which take a
  program JSON with a `fit` section;
- the explorers in `navette::fit::explore`;
- the JSON doors in `navette::fit::doors`, which are what `navette.fit`
  calls.

`fit_watched` and `run_strategy_watched` take a `Watcher`, a progress
callback that can stop the fit.
