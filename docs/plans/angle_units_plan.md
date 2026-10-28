# Angle units — one unit inside, explicit units at the door

Status: planned, not built. This replaces the `angles_in_radians` bool
(253d8f6). Backwards compatibility is not a goal: names, keyword
arguments and target JSON keys change.

## A1. What it changes

Every calculation in Rust and Python works in radians. A degree value
exists only in two places: a user typing one at a door, and a value
printed for a person (reports, warnings, phase dicts, plot labels).
Nothing in between stores, matches or passes degrees.

The user can still give angles in either unit, but has to name the unit
each time. There is no default unit, so an angle can no longer be read
in the wrong one without anyone noticing.

Colour hue (CIE ΔE formulas, DIN99, `color/common.rs::DEG2RAD`) is out
of scope. It is a colour-space angle, it is defined in degrees by the
standards, and it never meets an incidence angle.

## A2. Rust: a `Radians` newtype

```rust
#[repr(transparent)]
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Radians(f64);

impl Radians {
    pub fn new(rad: f64) -> Self;           // range-checked, see A4
    pub fn from_degrees(deg: f64) -> Self;  // deg.to_radians(), range-checked
    pub fn get(self) -> f64;
    pub fn sin(self) -> f64;
    pub fn cos(self) -> f64;
    pub fn to_degrees(self) -> f64;         // display only
}
```

It lives in a new `smatrix/angle.rs` and is re-exported from the crate
root. There is no `From<f64>` and no `Deref`: the only ways in are the
two named constructors, so a bare `f64` cannot become an angle by
accident. `repr(transparent)` makes `&[Radians]` the same memory as
`&[f64]`; the type costs nothing at run time.

What takes `Radians` (or `&[Radians]`) instead of a degree `f64`:

| Site | Today | After |
|---|---|---|
| `merit.rs` `MeritKey.angle` | degrees, f64 | `Radians` |
| `merit.rs` `SimCurves.angles` | degree row labels | `Radians` labels |
| `merit.rs` `reference_rotation` / `reference_rotation_radians` | two doors | one, `angle: Radians` |
| `targets.rs` angle fields (:36, :68, :161) | degrees | `Radians` |
| spectralweave `OpticalKey` | generic f64 slot | stays f64 (see below) |
| `optics_core.rs` `reference_wavenumber(.., theta_inc_deg)` | degrees | `theta_inc: Radians` |
| `pipeline.rs` `SpectralInputs.angles_deg`, `from_spec`, `from_spec_radians` | degrees + radian twin | `angles: Vec<Radians>`, one `from_spec` |
| `driver.rs` `run_design(.., angles_in_radians)`, `run_environments_angles`, `fold_and_sines` | bool | `&[Radians]`, one `run_environments` |
| `solver.rs:464`, `:1198` `angles` + `radians: bool` | bool | `&[Radians]` |
| `evaluator.rs` row labels (`asin().to_degrees()`) | degrees | `Radians::new(asin)` |
| `core_engine.rs:830` (asin → degrees) | degrees | `Radians`, degrees only if printed |

`OpticalKey` is spectralweave's generic curve key (its number field is
named `wavelength` and keys any 1-D axis). The synthesis compile only
borrows it as a dictionary key, so it stays an `f64`: `targets.rs` stores
the radians' bits there and reads them back with `Radians::new`. The
Python TargetWeaver doors convert at the binding in commit 2, like every
other door.

The solver kernels keep `sin_theta: &[f64]`. They do not see angles at
all, only the precomputed sines, so the hot loops do not change.

Angle matching (`angle_row`, the needle fold's nearest match in
`build_needle_targets_env`) compares `Radians` with `Radians`. It is a
plain nearest match (argmin, no tolerance), so it works the same in
either unit, and it stays as it is.

## A3. Python: `*_deg` or `*_rad`, exactly one

Every door that takes an angle has two keyword-only arguments, one per
unit. Exactly one must be given; both or neither is a `ValueError` that
names the two arguments.

- Grids: `angles_deg=` / `angles_rad=`.
- Single angles (targets, merit keys, reference rotation):
  `angle_deg=` / `angle_rad=`.

The doors:

- `ScatterMatrix`, `solve_structure`
- `SmatrixContext`, `NeedlePipeline`, `run_design`, `run_design_environments`
- `run_needle` (flat and environments)
- `reference_rotation`, `apply_reference_rotation`
- `SimCurves`, `sim_curves_from_arrays`
- `build_needle_targets`
- `MeritSpec.add_key`
- `SpectralTarget`, `AngularTarget`, `ColorTarget`
- spectralweave `TargetWeaver.add_spectral_target` / `add_angular_target`

The conversion happens once, in Rust, at the binding: the PyO3 layer
builds `Radians` from whichever argument came in and passes only
`Radians` on. Python-side wrappers forward the keyword as they got it;
they do not convert. That leaves one conversion site per unit in the
whole code base.

Python dataclasses (`SpectralTarget` etc.) store one field, `angle_rad`
(or `angles_rad`), and expose `angle_deg` as a read-only property for
display. `__post_init__` takes whichever unit was given.

Target JSON uses the key `"angle_deg"` or `"angle_rad"` (and
`"angles_deg"` / `"angles_rad"`). A file with the old `"angle"` key is
rejected with a message that names the new keys. It is not read as
degrees.

Removed: the `angles_in_radians` bool everywhere (including
`ScatterMatrix` and `solve_structure`), and the helpers `deg_to_sin`,
`angles_to_sin`, `angles_as_degrees`, `_angles_to_degrees`,
`fold_and_sines`, `from_spec_radians`, `run_environments_angles` and
`reference_rotation_radians`.

Returned dicts that carry angles (phase dicts, run reports) carry
radians under `*_rad` keys. Where a person reads them (report text,
warnings), the text shows degrees.

## A4. Range check

`Radians::new` and `Radians::from_degrees` reject anything outside
[0, π/2] (or [0°, 90°]) and any NaN, with an error that gives the value
and the unit it was read in. A grid `[0, 30, 60]` passed as `angles_rad`
then fails at once (30 rad is out of range) instead of solving a
nonsense problem. That is the main point of the check: the most likely
mistake, degrees given as radians, cannot pass silently.

90° itself is allowed, because grazing incidence is a valid input to
the solver today. If some caller turns out to rely on angles above 90°
or below 0, the survey in A7 finds it before the check goes in.

## A5. Bit-exactness: the one decision to make

The code base converts degrees two ways today:

- `a.to_radians()`, which is `a * (PI / 180)`: in `solver.rs`, `driver.rs`,
  `merit.rs::reference_rotation` and `optics_core.rs`.
- `a * PI / 180`, which is `(a * PI) / 180`: in `pipeline.rs::from_spec`
  and the bindings' `deg_to_sin`.

The two differ in the last bit for 2517 of the 9001 angles from 0.00° to
90.00° in 0.01° steps. On the grid the forward golden baseline pins (0,
15, 30, 45, 60, 75, 89°), they differ only at 89°. Even there the sine
comes out the same, and only the cosine differs.

`Radians::from_degrees` has one formula, `to_radians`. My expectations,
to be confirmed in the build:

- **Forward golden baseline:** stays bit-exact. Its synthesis path feeds
  only sines into the solver, and those match at every pinned angle.
- **Needle pin:** stays bit-exact. It runs at normal incidence, where
  both formulas give 0.
- **Elsewhere:** a run at, say, 37.3° through the synthesis door can move
  by one ulp in its angle and so in its result. No test pins such a run
  bit for bit.

If either pin does move, the build stops and I report the diff. Then
the choice is the user's: re-pin with the reason in the commit, or
add a second constructor that keeps `(a * PI) / 180` for the old sites.
My recommendation is to re-pin. Two constructors would bring back the
very inconsistency this refactor removes.

## A6. Performance

The newtype is zero-cost (`repr(transparent)`, `Copy`, inlined
accessors). Each solve converts its angle grid once, which costs
O(n_angles) and is negligible next to the O(n_angles · n_wl · n_layers)
solve. The hot loops keep their precomputed sines and are not touched.
Angle matching compares two f64s, as it does today.

To check this, I run `validation/benches/` before the change, on
253d8f6, and after it, and put the numbers in the commit message. A
slowdown outside run-to-run noise blocks the commit.

## A7. Callers to update

Counted on 253d8f6 (`.py`, `.md`, `.ipynb`, outside `.venv` and `target/`):

| Door | Uses | Files |
|---|---|---|
| `ScatterMatrix(` | 117 | 38 |
| `SpectralTarget(` | 87 | 21 |
| `run_needle(` | 71 | 23 |
| `sim_curves_from_arrays(` | 71 | 10 |
| `add_key(` | 70 | 11 |
| `build_needle_targets(` | 41 | 10 |
| `ColorTarget(` | 29 | 2 |
| `SmatrixContext(` | 24 | 13 |
| `solve_structure(` | 14 | 9 |
| `apply_reference_rotation(` | 12 | 4 |
| `NeedlePipeline(` | 9 | 5 |
| `AngularTarget(` | 8 | 6 |
| `angles_in_radians` | 54 | 8 |

Not every use passes an angle. Calls that leave the angle out get the
normal-incidence default. On that default, see the open question in A9.

Rust: 93 `MeritKey { .. }` literals, most of them tests with
`angle: 0.0`. They become `angle: Radians::ZERO` (a `const`).

Identifiers with "deg" in their name (`angles_deg` 80, `angle_deg` 16,
`theta_inc_deg`, `TH_DEG` 36 in validation, `brewster_deg`, …) are
renamed to `*_rad` where they hold a working value. Where they hold
what the user typed or what gets printed, they stay `*_deg`.

Also updated: `src/navette/_smatrix.pyi`, `_spectralweave.pyi`,
`tools/check_pyi_sync.py` (it must pass), the docs under `docs/` that
show calls, `docs/forward-golden.md`, and the notebooks/examples.

The mechanical rewrite is done by a scratchpad script that refuses any
anchor that does not match exactly once, followed by a grep for
leftover `angles_in_radians`, `deg_to_sin`, bare `angle=` keyword uses
and `"angle":` JSON keys. The grep must come back empty.

## A8. Order of work and commits

Each commit needs the user's approval and leaves the full battery green.

1. **Rust core:** `Radians`, every site in the A2 table, the Rust tests
   and the 93 literals. The bindings convert at their edge for now,
   behind the old Python signatures, so commit 1 alone changes no Python
   behaviour. Checks: the cargo tests, the golden dry run, and the
   needle pin.
2. **Python boundary:** the `*_deg`/`*_rad` doors, the range check, the
   new target JSON keys, the removals from A3, and the stubs plus
   `check_pyi_sync.py`. `test_angle_units.py` is rewritten: both units
   agree, both or neither raises, out of range raises, and the old
   keyword and the old JSON key are rejected.
3. **Callers:** validation, docs, examples, notebooks. Then the benches,
   before and after, and the full battery.

The battery is the release one: `cargo test --release -p navette`, then
`maturin develop --release`, then `pytest validation -q`. On top of it
come the golden dry run ("OK: baseline reproduced") and the needle pin
fingerprint.

## A9. Open questions for the user

- **Normal-incidence default.** Today an omitted angle means 0°. Should
  an omitted angle stay valid and mean 0 (0 is the same in both units,
  so it cannot be misread), or should every door require an angle? The
  plan assumes 0 stays a default, because it is unambiguous.
- **A one-ulp re-pin** (A5), if either fingerprint moves.
