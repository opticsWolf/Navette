# navette::color — design notes

Rust rewrite (parity port) of the **Unified Color Engine**
(`loom_colorengine.py`, "Gold Master", LGPL-3.0, opticsWolf).

> **Naming.** The project is **Navette**; *Loom* was its internal name while
> this page was written, and the page keeps it below. The mapping is
> mechanical: `loom` → `navette` for the project, `loom-core`/`loom._core` →
> the `navette` engine crate, `loom-py` → `navette-py` (which builds
> `navette._navette`). Where a name below belongs to the **numba reference
> implementation** — `loom_colorengine.py`, `loom_matrix.py`, the `refs/`
> oracles — it is not drift and is left alone: those files really are called
> that, and the parity tests import them by name.

The module shipped as `navette::color` (`rust/navette/src/color/`), bound as
`navette._color`; `loom_color` was the working name on this page.

## Module map

The port-task numbering `func_01`–`func_16` became descriptive filenames in
0.6.22 (R6.6). The **Was** column is the decoder ring: the port records on this
page and its siblings (`func_NN_*.md`) are dated documents and keep their
original names, as does the `func_01–16` catalog in §11 of the code review.

| Module            | Was      | Conversion / metric                       |
|-------------------|----------|-------------------------------------------|
| `xyy`             | func_01  | XYZ ↔ xyY                                   |
| `lch`             | func_02  | Lab ↔ LCh                                   |
| `luv`             | func_03  | XYZ ↔ CIELUV                                |
| `oklab_xyz`       | func_04  | XYZ ↔ Oklab (direct-XYZ matrices)           |
| `oklab_srgb`      | func_05  | sRGB ↔ Oklab (legacy sRGB-baked matrices)   |
| `uvw1964`         | func_06  | CIE 1964 U*V*W*                             |
| `ucs1960`         | func_07  | CIE 1960 UCS & chromaticity                 |
| `bradford`        | func_08  | Bradford chromatic adaptation               |
| `delta_e_76`      | func_09  | Delta E 76                                  |
| `delta_e_94`      | func_10  | Delta E 94 (`De94Params`)                   |
| `delta_e_cmc`     | func_11  | Delta E CMC(l:c)                            |
| `din99`           | func_12  | DIN99                                       |
| `spectral_srgb`   | func_13  | spectral pipeline (SPD → sRGB)              |
| `photometry`      | func_14  | photometry engine                           |
| `shapes`          | func_15  | shape handling & broadcasting (`map_pairs`) |
| `delta_e_2000`    | func_16  | CIEDE2000                                   |

The single-digit `func_0`–`func_5` that appear in `smatrix/` are the *numba
reference implementation's* module names — a different numbering, imported by
name by the parity tests, and untouched by this rename.

The foundational sRGB↔XYZ and XYZ↔Lab conversions plus transfer functions live
in `common.rs`; generated matrices in `matrices.rs`.

## Build & test
```bash
cargo build
cargo test     # 55 unit + 15 parity + 15 doc tests, all pass
```
Parity vs the Python reference: worst-case deviation anywhere is **2.5e-14**
(`luv`, CIELUV), otherwise mostly exactly 0 — machine precision.

### Optional `parallel` feature (rayon)
`rayon-core` 1.13 needs `rustc >= 1.80`. Built/tested on rustc/cargo **1.75**
(edition 2021) with the default sequential build. To use `--features parallel`
on 1.75, pin an older rayon (`cargo update -p rayon-core --precise <older>`),
otherwise use rustc >= 1.80.

## Reproducing the golden vectors
```bash
cd refgen
python gen_matrices.py   # -> ../src/matrices.rs
python gen_golden.py     # -> ../tests/golden.rs
```
`numba.py` is a no-op shim so the engine runs in pure NumPy. Drop a copy of
`loom_colorengine.py` into `refgen/` first.

## Parity convention notes
- Matrices stored natural; `mat3_mul_vec(M,v)[j] = Σ M[j][i]·v[i]`.
- `calc_transform_matrix` (Bradford) returns a **row-vector** matrix:
  `adapted = white @ M`, i.e. `out[j] = Σ_i white[i]·M[i][j]`.
- Oklab uses sign-preserving `signed_pow` (not `cbrt`); `sign(0)=0` per NumPy.
- `LAB_EPSILON = (6/29)³` (cubic/linear switch), `LAB_DELTA = 6/29`.
- Delta-E 94 / CMC weights use the reference sample `lab1`; `dH²` clamped ≥ 0.

## License
LGPL-3.0-or-later, matching the upstream engine.
