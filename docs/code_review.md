# Navette — Deep-Dive Code Review

**Project:** Navette — 1D optical engine (scattering-matrix based) for thin-film simulation
**Stack:** Rust core (`rust/navette` engine crate + `rust/navette-py` PyO3 bindings) with a thin Python package (`src/navette`), maturin build, pytest validation suite
**Scope of review:** ~36,500 lines of Rust, ~18,800 lines of Python (package + validation), 286 files
**Method:** CodeRadar static analysis (smells, scaffolding/AI-slop detection, dead code, clone detection) + manual physics review of the S-matrix engine, materials models, synthesis pipeline + review of benchmarks and test coverage

**Legend:** 🔴 critical · 🟠 major · 🟡 minor · 🔵 observation/nit

---

## 1. Executive Summary

The full summary is in **§8** (written after the first pass) with a post-expansion addendum in **§15** (covering structure/persistence, interpolation, color, materials dispatch, targets, the binding layer, and edge-case numerics). Read §8 first, then §15.

---

## 2. Code Smells (CodeRadar static analysis)

CodeRadar indexed 216 files / 316 classes / 2,882 functions / 5,705 call edges and reported **2,251 findings**: 29 critical, 1,833 high, 305 medium, 84 info.

| Rule | Count | Severity mix | Verdict |
|---|---|---|---|
| dead-code | 1,528 | mostly High | **mostly false positives** (PyO3 entry points + `attic/` archive), a handful real |
| long-parameter-list | 286 | High/Med | real; worst offender has **31 parameters** |
| long-method | 178 | High/Med | real; worst `solve_coherent_block_fields` LOC=198 |
| deep-nesting | 105 | Med | real; worst `eps2_at` nesting=7, `calculate_merit` nesting=9 |
| high-cyclomatic-complexity | 89 | Critical/High | real; see hotspot table |
| data-class | 29 | Info | mostly config/result carriers, acceptable |
| excessive-returns | 17 | Med | validation-heavy code, mostly justified |
| brain-method / too-many-fields | 19 | Critical | real design pressure in structure/synthesis types |

### 2.1 Complexity hotspots (cyclomatic ≥ 30)

| Entity | File | Cyclomatic | Notes |
|---|---|---|---|
| `needle_gradient` | `rust/navette/src/smatrix/solver.rs:698` | **96** | 31 parameters, 7 copy-pasted target/weight blocks (see §4.1) |
| `build_needle_targets` | synthesis | 67 | target-kind dispatch ladder |
| `fold_integral_demand` | synthesis | 40 | envelope/integral folding |
| `Layer::set_properties` | `rust/navette/src/structure/layer.rs` | 43 | setter doing type-dispatch + validation |
| `Group::set_properties` | `rust/navette/src/structure/group.rs` | 40 | twin of the above |
| `check_color_demand` | `synthesis/color_merit.rs` | 35 | 13 return statements |
| `levenberg_marquardt` | `synthesis/thick_opt.rs` | 34 | optimizer core loop |
| `residuals_into` | `synthesis/merit.rs` | 32 | merit assembly |
| `MeritSpec` (brain class) | `synthesis/merit.rs` | 32 max / WMC 65 | accumulation + eval + serialization in one type |
| `Solver` (brain class) | `smatrix/solver.rs` | 31 max / WMC 70 | solve + derivatives + eigen + fields |
| `Layer` (brain class) | `structure/layer.rs` | 43 max / WMC 56 | |
| `Group` (brain class) | `structure/group.rs` | 40 max / WMC 62 | |
| `expand` | `structure/expansion.rs` | 31 | graded-layer expansion |
| `solve` | `smatrix/solver.rs` | 31 | |
| `new` (Parser?) | misc | 31, 21 returns | |
| `evaluate` ×2 | synthesis | 30/31, 22 returns | |
| `MaterialSpec` (brain) | `structure/specs.rs` | 30 max | |
| `ColorDemand` (brain) | `synthesis/color_merit.rs` | 31 max | |
| `ClampedNeedleSynthesizer` (brain) | `attic/needle_pipeline.py` | 18 max | archived |

**Assessment:** the smell profile is characteristic of a **performance-first numeric kernel written to be bit-identical with a legacy reference**: hand-monomorphized s/p solvers, macro-generated request plumbing, and 7-way copy-paste in the needle gradient. The physics core itself (`optics_core.rs`, `coherent_block.rs`) is clean; the debt is concentrated in (a) the needle/synthesis gradient machinery, (b) the `structure` CRUD types (`Layer`/`Group` setters), and (c) validation ladders in `config`.

### 2.2 Dead code — interpretation required (1,528 findings)

CodeRadar marks 1,528 entities unreachable. Two systematic false-positive patterns must be excluded before acting:

1. **PyO3-boundary blindness.** Every `#[pymethod]`/`#[pyfunction]` is an entry point from Python, but the analyzer cannot see the interpreter. Examples flagged as dead that are **exposed and used from Python**: `Solver.solve`, `Solver.new`, `PySolver.needle_gradient`, `PyTargetWeaver.add_spectral_target` / `add_angular_target` / `export_entries`, `UniInterpolator.new/evaluate`, `_color`, `_smatrix`, `ColorDemand.new`, `PyNeedlePipeline.new`, `PyDesignStack.from_design`. `tools/check_exposure.py` (bidirectional exposure lint, run in CI) is the correct authority for what is required at the boundary — **not** the dead-code report.
2. **Archived code.** `attic/` (~28 files: `needle_pipeline.py`, `structure-monolith.py`, `materials-rust-wrappers/*`, …) is deliberately retired; ~300 findings land there. Either exclude `attic/` from indexing (`.coderadar` config) or delete the directory — it currently pollutes every static-analysis run.

Dead-code findings that look **real** (verify with `affected` before removal):

| Entity | Location | Why it looks real |
|---|---|---|
| `MaterialSpec::evaluate` | `rust/navette/src/structure/specs.rs` | ~156 lines; superseded by the native weaver path |
| `Structure::bake_materials` / `bake_films` | `rust/navette/src/structure/structure.rs` | replaced by `DesignStack` flow |
| `Layer::set_properties` / `Group::set_properties` | `structure/layer.rs`, `group.rs` | 40+ cyclomatic setters with no caller in the crate |
| `MeritSpec::add_target` | `synthesis/merit.rs` (~71 lines) | targets appear to be added via `compile_merit_spec` now |
| `SmatrixContext::optimize_thicknesses` | `synthesis/evaluator.rs` (~50 lines) | superseded by `thick_opt.rs` LM |
| `Group::from_state`, `Architect::from_state`, `Layer::deserialize` | `structure/` | serialization round-trip leftovers |
| `load_program` (`src/navette/config/program.py`), `save_material_library` (`config/loader.py`) | Python | config API moved on |
| `TargetCollection.build_weaver` (`spectralweave/target.py`) | Python | superseded by native `PyTargetWeaver` |
| `marked_stacks_and_films_blocks`, `norm_factor_for` | `structure/architect.rs`, `spectralweave/targetweaver.rs` | retired helpers |

🟡 **Minor:** `rust/navette/src/smatrix/solver.rs` keeps a `flat_cache` (re/f64 interleaved) built in `Solver::new` that no `solve*` path reads — a candidate for deletion unless the binding layer consumes it (grep `flat_cache` before removing).

### 2.3 Long parameter lists — worst offenders

| Function | Params | Comment |
|---|---|---|
| `needle_gradient` (free fn) | **31** | 7 target/weight `Option<&[f64]>` pairs → should be `struct NeedleRequest` |
| `Solver::new` / `from_raw` | 18 / 10 | borderline acceptable for a numeric ctor |
| `redheffer_product_real_inner` / `..._complex_field` | 8 | structural (S4 tuple unpack); fine |
| `p_coherent_from_fields` | 11 | gradient plumbing |
| `scan_box` | 12 | eigenvalue search box config |
| `solve_coherent_block_fields_inner` | 10 | hot loop; keep as-is for perf |
| `sellmeier_urbach_nk` | 10 | 3×(B,C) + 3 tail params — a `SellmeierUrbachParams` struct exists conceptually in `config` |
| `delta_e_2000` | 7 | CIE ΔE2000 signature is standard; leave |
| `__init__` (Python holders) | 9–13 | config holders; use dataclasses with fields |

🟡 The 7 `Option<&[f64]>` target/weight pairs in `needle_gradient` repeat a scalar-or-vector broadcast (`load_pair`) **seven times verbatim** — the single most mechanical AI-slop pattern in the crate (see §4.1).

### 2.4 Scaffolding / AI-slop signals (CodeRadar `find_scaffolding`)

79 findings:

- **27 placeholder bodies** — 21 are `.pyi` type-stub signatures (`src/navette/_materials.pyi`: `...` bodies). **False positives**: that is what stubs are. Real ones: `src/navette/structure/materials.py:41-42` (`get_nk`, `contains` raise-not-implemented style stubs — verify these are documented protocol stubs, not missing features).
- **52 phase/comment markers** — 41 in `attic/` and `docs/plans/*` (planning documents, fine). In live code: none material. No `TODO`/`FIXME`/`HACK` debt in the Rust core — unusually clean.
- **No temp-file naming, no hardcoded secrets detected.**

---

## 3. Physics Review (manual deep-dive)

The physics core was reviewed line-by-line: `optics_core.rs`, `coherent_block.rs`, `core_engine.rs`, `solver.rs`, `needle_operator.rs`, `thick_opt.rs`, all of `materials/`, and the color `func_01–16` modules. Overall verdict: **the engine is physically sound and unusually well-documented** (conventions, sign choices, and known noise amplifications are written down at each site), with **one confirmed energy-conservation bug in the Névot-Croce roughness path** (§3.2) and several model choices that deserve explicit documentation or verification.

### 3.1 S-matrix engine — correct

* **Admittance formulation**: s → `y = n·cosθ`, p → `y = n/cosθ` with a `1e-12` magnitude floor on `cosθ` (prevents the p-pol singularity at grazing). Amplitudes `t = 2y₁/(y₁+y₂)`, `r = (y₁−y₂)/(y₁+y₂)`; intensity renormalization `T = |t|²·Re(y_exit)/Re(y_inc)`. This is the standard wave-admittance (Macleod) treatment; the p-pol T-factor is correct because `y_p = n/cosθ` is the true wave admittance.
* **Branch safety**: `cosθ = √(1−(n₀sinθ/n)²)` with a forced `Im ≥ 0` root, and propagation phase `β` sign-flipped to keep `Im(β) ≥ 0` — a consistent decay convention throughout (`e^{i(kz−ωt)}` style). The flip appears at every site that could pick the growing branch (coherent solvers, incoherent `tau`, needle operators). ✔
* **Redheffer star products** (three variants in `optics_core.rs`): complex-field (amplitudes), real (intensities), and cross-coherency `C = p·s̄`. All three are structurally the same `A⋆B` with denominator `1 − r_A^back·r_B^front`. Regularization is sensible everywhere: complex-field preserves phase at `|denom| < 1e-100`; the real product drops the feedback term (`inv = 0`) when `|1 − R·R| < ε` (probabilities, so `denom ∈ (0,1]` — no negative-denominator pathology possible); the cross product mirrors the real one. 🔵 nit: the complex-field regularization adds `1e-300` to the phase denominator (`LOG_MIN·phase + 1e-300`) — harmless but unexplained next to the other two.
* **Hybrid coherence (Modes A/B/C)**: incoherent spacers propagate intensity only (`τ = e^{−2·Im(β)}`, clamped at `β_imag ≥ 0` so numerical noise can't manufacture gain); the cross-coherency channel through incoherent spacers decays by the same `τ`, with a documented geometric-series argument for why the first-order product (not the full complex star) is correct for the cross channel. The docstring reasoning is sound and the `Delta_T off-by-π` bug from the original port was genuinely fixed (cross-term now uses `t_fwd`, verified in code and in `docs/plans/smatrix/REVIEW.md` §2.1).
* **Energy conservation**: `A = 1 − R − T` is computed from totals including all incoherent echoes; a dedicated `energy_conservation()` residual helper exists and is tested. Verified empirically below for roughness (where it currently fails — §3.2) and clean for the standard paths.

### 3.2 🔴 CONFIRMED BUG — Névot-Croce (roughness type 5) violates energy conservation on transmission

`coherent_block.rs`, `solve_pol_specialized` / `solve_coherent_block_fields_dual`, `rtype == 5` branch applies the **same** Névot-Croce factor to the transmission coefficients as to reflection:

```rust
} else if rtype == 5 {
    let f = (-2.0 * kz1 * kz2 * sigma * sigma).exp();
    (r12 * f, r21 * f, t12 * f, t21 * f)   // ← t12, t21 must NOT get f
}
```

Canonical Névot-Croce corrects **r** by `exp(−2·kz1·kz2·σ²)` but **t** by the wavevector-difference factor `exp(±(kz1−kz2)²σ²/2)` — for a low-contrast interface (`kz1 ≈ kz2`) the transmission factor must approach **1**, while the code applies `exp(−2kz²σ²) ≪ 1`.

Measured repro (single near-matched interface, n = 1.0 → 1.5001, λ = 550 nm, normal incidence, **no absorption — there is nowhere for energy to go**):

| σ (nm) | type 5 R+T | type 4 (Gaussian) R+T |
|---|---|---|
| 0 | 1.000000 | 1.000000 |
| 3 | 0.992977 | 0.999530 |
| 6 | 0.972202 | 0.998128 |
| 10 | **0.924678** | 0.994837 |

The engine silently destroys 7.5 % of the energy at a *single almost-invisible interface*. The general w-function path (types 1–4) handles transmission correctly (`ga = W((kz1−kz2)σ)`) — only type 5 is affected. Severity grows with σ and with n-contrast; XRR stacks (the advertised use case) have large contrast, where the loss is smaller but never zero. Likely inherited from the numba reference (parity tests pin the behavior), which is why the mirror tests do not catch it. **Fix:** `t12·ga, t21·ga` with `ga = exp(−(kz1−kz2)²σ²/2)` (or `+`, after fixing the convention against the reference); add an energy-conservation regression test sweeping `rtype × σ × contrast × angle`. The incoherent-`tau` clamp (`beta_imag ≥ 0`) shows the team knows this pattern — the roughness path just never got the same treatment.

🟡 Related model choice worth documenting: types 1–4 apply **per-side** characteristic factors to reflection (`r12·W(2kz1σ)`, `r21·W(2kz2σ)`), which is the "graded transition layer" interpretation (measured above: R+T ≈ 0.995 at σ=10 nm — behaves like a smooth profile, only ~0.5 % loss). This differs from the *correlated* Debye-Waller form `exp(−2kz1kz2σ²)` that NC uses for reflection; the two models can disagree substantially at large σ. A short doc note in `RoughnessType` distinguishing "graded profile (1–4)" vs "NC decorrelation (5)" would prevent misuse.

### 3.3 Ellipsometry / Stokes — correct, two nits

* `tanΨ = √(Rp/Rs)`, `Ψ→90°, Δ→0` below the signal floors (`RS_FLOOR = 1e-12`, `TS_FLOOR = 1e-20` — the asymmetric floors are deliberate and documented: T can be vanishingly small in absorbing stacks). Δ extracted as `atan2(S3, S2)` on the mode-resolved cross-coherency, i.e. `arg(rp·rs̄)` plus the convention π-shift; reflected S2/S3 carry a minus sign, transmitted do not — the standard reflected-basis handedness flip. The "Azzam & Bashara" claim in the README is anchored by parity with the numba reference (`test_core_engine_rigorous_ellipsometry.py`, currently marked NEEDS PORT — see §6.3) rather than by a textbook cross-check; the port should be finished so the claim has a live guard.
* 🟡 `DOP_T` is clamped to ≤ 1.0 but `DOP_R` is not — floating-point round-off can push the reflected DOP epsilon-above 1. Asymmetric guard; clamp both or neither.
* 🔵 `s2r = -2.0 * s.cross_r.re + 0.0;` — the `+ 0.0` is a no-op (leftover from mirroring Python float semantics). Harmless; appears twice.

### 3.4 Dispersion post-pass — correct with documented caveats

Phase unwrap + repeated `grad_nonuniform` (bit-faithful `numpy.gradient` semantics, verified formula) for GD/GDD/TOD/FOD; ω = 2πc/λ with `C_NM_PER_FS = 299.792458`. The docstring honestly flags that TOD/FOD amplify solver noise and repeated-gradient error, recommending spline-fit validation — good practice. The same derivative chain is reused by the needle operator's spectral sensitivity (`phase_dispersion_sensitivity`), keeping the two consistent. 🟡 The dispersion helpers assume the wavelength grid is **sorted ascending per angle row**; nothing validates this (an unsorted grid silently produces garbage derivatives). One `assert`-grade check in `Solver::new` would close it.

### 3.5 Materials — faithful and mostly clean

* **Sellmeier (+Urbach)**: standard 3-term form, k=0. 🟡 No pole guard: λ² = Cᵢ → `n = inf`; n² < 0 below resonance → `sqrt(neg) = NaN` that propagates silently into the stack. Faithful to the reference, but a domain check (or at least a doc note) is warranted since the Python side validates almost nothing else this strictly.
* **Drude / Drude–Lorentz**: ε = ε∞ − ωₚ²/(E²+iγE) + Σ fⱼE₀²/((E₀²−E²) − iEΓⱼ); the `E + 1e-12` guard and the (documented!) quirk that the Lorentz damping term also uses the guarded E are exact-fidelity choices — good that the deviation is written down.
* **EMA**: Lichtenecker, Looyenga, Birchak power law (α→0 correctly limits to Lichtenecker), Maxwell-Garnett (spheres, correct standard formula), Mori–Tanaka (ellipsoids, correct MT reduction), Wiener bounds (harmonic/parallel, correct denominators), and the 50:50 Looyenga roughness interface. All verified against the standard forms. ✔
* **Kramers–Kronig (FFT Hilbert)**: odd extension, unnormalized-realfft correction, Nyquist/DC handling to match `np.fft.irfft` bit-for-bit, 8192-point 0.01–80 eV grid behind a `OnceLock`. Careful work. 🟡 The 80 eV ceiling is a silent modeling limit for UBF/Cody-Lorentz users in the EUV/X-ray range; the interpolation clamps (rather than extrapolates) beyond the grid — worth a line in the docs.

### 3.6 Synthesis (needle + LM) — numerically reasonable, structure-heavy

* The needle gradient convention (`½`-merit form `2w(R−Rₜ)·Re{r̄ ∂r/∂δ}`, full factor for phase, flux-corrected T/A) is documented precisely in `needle.py` and implemented consistently in `needle_operator.rs`. The differential-phase reference-phase machinery (`reference_phase`, pinned by a unit test on the propagation sign) keeps absolute and differential demands mutually consistent — the docstring's warning that textbook-imported targets need conjugating is exactly the kind of note that prevents silent mis-targeting.
* **Bounded Levenberg–Marquardt** (`thick_opt.rs`): central-difference Jacobian with `h = cbrt(ε)·max(|x|,1)` (standard optimal step), Marquardt diagonal damping, λ-escalation (≤40 tries), bound-aware projection, and five distinct termination reasons. 🟡 It solves the **normal equations (JᵀJ)** — squaring the condition number. Thin-film stacks with correlated layers (the classic case the needle method creates) are exactly where JᵀJ goes singular; the λ-floor of `1e-14·diag` bails it out, but a QR-based solve would be sturdier. scipy parity exists (per plan docs) so current behavior is anchored — fine for now, worth a comment in the code.
* 🔵 `build_jacobian` re-allocates two residual buffers per column per iteration (parallel tasks own them, so it's correct, just allocator-chatty). Negligible vs the solve cost.

### 3.7 Color — spot-checked, correct

`delta_e_2000_single` (func_16) verified term-by-term against Sharma's CIEDE2000 reference implementation: G/a′/C′/h′, hue-difference wrapping, T, RT with Δθ = 30°exp(−((h̄′−275)/25)²), and the SL/SC/SH weighting all match, with a sensible `c1·c2 > 1e-12` degeneracy guard. CIE tables are CI-verified against sources (`tools/check_cie_sync.py` — passes). The color module carries its own golden matrix (`golden.rs`) generated from `colour-science`.

---

## 4. AI-Slop / Duplication Patterns (manual + CodeRadar)

The codebase is **not** typical AI slop — there is no placeholder meat, no contradictory comments, no invented APIs; scaffolding detection is essentially clean outside `attic/` and `.pyi` stubs. What *is* present is a specific, mechanical copy-paste fingerprint:

### 4.1 🟠 `needle_gradient` — 31 parameters, 7× copy-pasted target/weight plumbing

`rust/navette/src/smatrix/solver.rs:698` (cyclomatic **96**, the #1 smell in the project):

* Seven `Option<&[f64]>` target/weight pairs (`targets_r/weights_r` … `targets_ab/weights_ab`), each handled by a **verbatim repetition** of the same 6-line `load_pair` + `target_x_of`/`weight_x_of` closure block — seven copies of ~8 lines that differ only by the suffix. A `struct ChannelDemand<'a> { target: Option<&'a [f64]>, weight: Option<&'a [f64]> }` plus one accessor removes ~50 lines and seven chances to mistype a suffix.
* The same seven-fold pattern then reappears twice inside the point loop (coherent `want_*` ladder, multiblock ladder) and once more in `PointOut` (14 fields of `[Option<Vec<f64>>; 2]`). Adding an 8th demand type today means touching ~6 places.
* 31-parameter signature: at the PyO3 boundary this is already a builder-shaped call; introduce `NeedleRequest` and the cyclomatic count drops by roughly the number of `want_` flags (each becomes a loop over an enum).

### 4.2 🟡 Solver triplication (deliberate, but undocumented as a maintenance contract)

Three near-identical interface sweeps exist: `solve_pol_specialized::<S>` (single-pol), `solve_coherent_block_fields_dual` (shared roughness/phase), and the intensity-only variant inlined in `core_engine::solve_point_intensity` (~60 lines duplicating `solve_point`). The docstrings claim byte-identical parity with the reference, which explains — but does not excuse — the duplication: `solve_point` and `solve_point_intensity` differ only in whether complex amplitudes are captured, and could be one function with a const-generic `LEVEL` (the dead branches already constant-fold, as the dual-solver docstring itself argues). Same for single-pol vs dual: the dual docstring says "per-polarization arithmetic is byte-identical", so the single-pol entry could delegate to `dual` with one pol off unless profiling says otherwise.

### 4.3 🟡 Request-bit duplication Python↔Rust with no sync guard

`src/navette/smatrix/smatrix.py:85+` hand-copies all 49 `REQ_*` bit positions into a Python `IntFlag` (and `needle.py` mirrors the `NREQ_*` set), with only a comment: "bit positions MUST stay in sync with core_engine.rs". Unlike the CIE tables (which have `tools/check_cie_sync.py` in CI), **no test asserts `Request.RS.value == native.REQ_RS`** — even though the native module already exports the constants (`m.add("NREQ_P", …)`). Inserting a bit mid-list in Rust silently corrupts every Python request. One pytest looping over the native constants closes the gap; the same hook can replace the hand-copy entirely (build the IntFlag from the exported constants).

### 4.4 🔵 Macro plumbing in `Solver::solve` — good pattern, keep

The `f64buf!`/`put!`/`keep_f64!` macro family that wires 40+ observable buffers is the *right* way to do this; it reads well and keeps the derive loop uniform. Mentioned only to contrast with §4.1: where the codebase reaches for macros it stays clean — the needle path predates that discipline.

### 4.5 CodeRadar scaffolding detail

27 placeholder bodies: 21 are `.pyi` stubs (false positives by design), 2 are `attic/`, 4 are `src/navette/structure/materials.py:41-42` protocol-stub raises — documented as such in the docstring above them. 52 comment markers: all in `attic/`, `docs/plans/`, or release notes. No secrets, no temp files. The AI-slop surface of the live code is the §4.1/§4.3 duplication only.

---

## 5. Benchmarks & Performance

### 5.1 Inventory

Benches live in `validation/benches/` (Python scripts, run explicitly): `smatrix/bench_backside_speed.py`, `synthesis/bench_refold.py`, `spectralweave/navette_{spectral,target}_bench.py`, `structure/bench_grid_assert.py`, `interpolate/1dinterpol_test_bench.py`, `color/bench_validate_color.py`. Results JSON is committed (`benches/smatrix/results/backside_speed_{before,after}.json`). There are **no Rust-native `criterion` benches** and **no CI performance gate** — every timing depends on someone remembering to run a script.

### 5.1.1 🔴 CRITICAL BENCHMARKING-HYGIENE FINDING — the venv ships a debug build

The first pass of this review measured "Rust slower than Python" across the board and the numbers were consistent and reproducible — because **the installed extension was an unoptimized debug build**. Verified: `src/navette/_navette.pyd` (22.8 MB, timestamps matching `target/debug/_navette.dll`) while the optimized release build (4.5 MB, `opt-level=3, lto=true, codegen-units=1`, Sep 6) sat unused in `target/release/`. The README's standard setup instruction is `maturin develop` — which defaults to **debug** — and every bench then compares that build against fully-JIT-optimized numba reference code. The interpolated throughputs and spectralweave "0.04×–0.64× speedups" in the first version of this review were artifacts of this.

Rebuilding with `maturin develop --release` (1 m 28 s) and re-running every bench produced a completely different picture (§5.2). **Consequences:**

* All committed benchmark results and any conclusions drawn from `maturin develop`-built environments are unreliable for performance claims.
* The README setup instructions must say `maturin develop --release` (and the bench scripts ideally assert a non-debug build, e.g. via a build-feature marker or `#[cfg(debug_assertions)]` guard exposed to Python).
* The CI absence (§7) is how this state persisted: nothing distinguishes "works" from "works fast".

### 5.2 Measured results (this review, Windows 32-core, **release build** after rebuild)

| Bench | Result | Verdict |
|---|---|---|
| `bench_backside_speed` (40 λ × 3 θ, 5 layers, full mask incl. complex amps + dispersion) | full mask **0.164 ms** median (was 0.56 ms debug) ≈ 1.4 µs/point all-observables | excellent |
| `bench_refold` (synthesis) | one LM optimize **2.1 ms**; re-fold 0.19 ms → **9.1 %** overhead; self-verifying OK flags | good — **but the bench prints "no insertions in either mode — bench is degenerate"**: the needle *insertion* path — the expensive part of synthesis — is never exercised |
| `navette_spectral_bench` (weaver, Python vs Rust) | Rust **faster everywhere except one path**: set_data 1.8–5.2×, get_weaved 1.5–7.8×, unweave_cached 2.8–4.2×, unweave_collection 4.1–8.0× at small/mid grids; **only `unweave_collection` regresses at scale** (0.22× at 100k pts × 512 keys; 0.48× at 500k × 16) | the §5.3 "Rust slower" claim is retracted — it was the debug build; the batch path at high key counts is the one real remaining issue |
| `1dinterpol_test_bench` (interp, Python vs Rust) | release: pchip 1M pts = **1.22 ms** (1.2 ns/pt) vs numba loom 1.59 ms — **Rust faster**; accuracy identical (≈1e-14 vs analytic) both sides | retracted regression — was 17× with debug |
| `bench_validate_color` (color) | **"Exact" parity vs colour-science** for ΔE76/94/CMC/DIN99/2000 including edge rows (black/white/near-black); Navette 23–36× faster than loom on ΔE batches | exemplary |
| `bench_grid_assert` (structure) | 0.9 µs / assert | fine |
| `navette_target_bench` (TargetWeaver, rust-only) | 200 targets: ingest 4.6 ms, merit 3.3 ms, merit_interp 6.3 ms | reasonable; no cross-engine comparison exists |

### 5.3 Re-evaluated: the weaver data-plane (was "up to 25× slower")

With a release build the weaver beats the numba reference at every operation except `unweave_collection` at large grids with many keys (measured 0.22–0.48×: at 100k pts × 512 keys, 42.6 ms vs 9.5 ms). The regression is consistent with per-key overhead in the batch entry: the binding rebuilds an `AHashMap<OpticalKey, Vec<f64>>` from a Python dict on every call (String clone per key), `frames.read().clone()` clones the whole frame collection per operation, and results are materialized per key. Recommendations: (1) the README's setup section must say `maturin develop --release` — without this the first-pass conclusion of this review would have stood and shipped; (2) make `unweave_collection` the documented batch API and optimize its boundary (borrow the Python dict arrays zero-copy like `color.rs` already does); (3) keep the numeric parity checks in the bench — they are the reason the correction here is trustworthy.

### 5.4 Gaps

* No scaling bench for the **core solver** (layers × wavelengths × angles) — the one number a thin-film engine lives or dies by. The backside bench fixes 120 points; a sweep to 10⁵–10⁶ points would show the rayon scaling curve and the serial derive-loop fraction (§3.3 note: the per-point derive pass in `Solver::solve` is a **serial** loop — `atan`/`sqrt`/`arg` across millions of points will eventually dominate; it is the obvious next Amdahl bottleneck).
* No memory benchmark (the README claims "minimized cache misses" / memory efficiency).
* The synthesis bench never triggers needle insertions (its own output says so) — refactor to a workload that actually inserts/merges layers, otherwise it validates only the cheap path.
* Consider `criterion` benches in the workspace (cargo bench) so the Rust-only core can be gated without Python.
* **Bench scripts should refuse to run against a debug build** — a `debug_assertions` probe or file-size sanity check would have made the §5.1.1 situation impossible. The interpolate bench additionally crashes on a stock Windows console (emoji + cp1252, §15); it should run with UTF-8-safe output.

---

## 6. Testing & Validation

### 6.1 Status: green and fast

* `pytest validation`: **355 passed, 0 failed in 1.14 s** (this review, Windows).
* `cargo test --workspace`: **376 passed, 0 failed** (339 engine + 22 materials-parity + 15 bindings).
* Strategy is unusually strong: **parity-vs-reference** (numba `loom_*` mirrors for solver/interpolate/color/spectralweave, `.npy` goldens for materials, `golden.rs` for color vs `colour-science`), **mirror tests** (Python surface vs native semantics), **random-stack differential bit-identity** (`test_random_stacks_bit_identical`), round-trip/persistence suites, and an exposure lint (`check_exposure.py`: 208 pub fns → 94 allowlisted internals, passes). The known-hard scenarios are pinned: inversion mirrors, shared-block aliasing, Redheffer field cases, needle T/A/φ, differential phase, merit weighting.

### 6.2 🟠 Broken shipped example

`examples/spectralweave_example.py` **fails on a fresh clone** at line 22: `weaver.unweave(key, full_wl, full_data)` → `ValueError: Length mismatch`. The Python wrapper (`OpticalWeaver.unweave`) now requires an `OpticalFragment` template, not a bare key tuple; the example predates that API change. Worse, the failure mode is an opaque `Length mismatch` from deep inside the native layer rather than an argument-type error — the wrapper should reject a non-`OpticalFragment` with a clear message. Fix the example, and add a smoke test that executes every file in `examples/` (they are currently run by no one and nothing).

### 6.3 🟡 Stale documentation of the test suite (audit trail rot)

* `validation/README.md` predates the current tree: it omits `validation/regression/**` (structure/config/synthesis/color — ~20 test files), `parity/synthesis/**` (6 files), `parity/color/**`, `goldens/` scope, and still marks `test_core_engine_photometry_only.py` and `test_core_engine_rigorous_ellipsometry.py` as "NEEDS PORT" without tracking them anywhere.
* `docs/plans/bug_fix_plan.md`: BUG-1/BUG-2 marked fixed via the (since-deleted) Python `expander.py` path; the open 🔴 BUG-A (inverted interface slice) references files that no longer exist — the behavior is now pinned by `test_asymmetric_interface_position_and_pair` (passing) in the Rust port, but the plan no longer says so. Close or re-scope the plan against the Rust modules; an open 🔴 that is actually fixed (or actually still open — the plan can't tell you which) is worse than none.

### 6.4 🟡 Coverage gaps

* **No energy-conservation test across roughness** (the §3.2 bug is invisible to the suite: parity tests compare against a reference with the same bug; nothing asserts the physical invariant). One parametrized test (`rtype × σ × contrast × angle`, assert `max|1−R−T| < tol`) would have caught it — the helper `energy_conservation` already exists and is even exported.
* No test executes `examples/` (see §6.2).
* No pytest coverage measurement configured (`pytest.ini` has no `--cov`), no coverage badge/threshold; for a numerics project, line coverage of the *validation* suite vs the Rust core would be more meaningful — the Rust side has no `tarpaulin`/`llvm-cov` setup either.
* No property/fuzz tests on the solver invariants (e.g., reciprocal stack `R_forward(n↔substrate)` identities, `T_p/s` reciprocity, random stacks at grazing angles with roughness + incoherent spacers — the differential test is a great start but fixes the mask to photometry).

---

## 7. Build, CI & Tooling

* **🟠 No CI on push/PR.** The only workflow is `release.yml` (tag-triggered): CIE-sync + version-match checks, 3-OS wheel build, PyPI/crates.io publish. Nothing runs `pytest validation`, `cargo test`, `check_exposure.py`, `clippy`, or `fmt` on ordinary commits — yet the README asserts the exposure lint is "enforced … in CI". All quality gates currently live in the developers' fingers. A 20-line `ci.yml` (pytest + cargo test + exposure + clippy `-D warnings`) closes the gap; the release workflow can then `needs: ci`.
* **🟠 Packaging: `numpy>=1.22.0` floor is wrong for the shipped wheels.** The extension is built with the `numpy` 0.28 crate (numpy-2 C-API) under `abi3-py312`; such a binary cannot load against numpy 1.x (clean `ImportError` at import time, but still a broken install for anyone satisfying the declared floor with e.g. numpy 1.26). The floor must be `numpy>=2.0` (or the claim "abi3" must be re-verified against the numpy crate's actual guarantees — its docs/Cargo.toml make no abi3 promise; the empirical build here works only because the runtime is numpy 2.5).
* **🟠 The standard dev setup ships a debug build** (`maturin develop` per README, no `--release`): every benchmark run against the documented setup compares unoptimized code against optimized numba (§5.1.1). One word in the README fixes it.
* **🟡 Compiler warnings are live** (observed during `cargo test`): unused imports (`NeedleTargets`, `cplx`, `PI`, `max_disp_order`, `rayon::prelude`, `ArraySeed`), unused variables (`land`, `m`×3, `wavelengths`, `num_angles`, `ok`) in `navette` and `navette-py`. No clippy gate exists; `-D warnings` in CI would have kept the tree clean.
* **unsafe** is confined to `navette-py/src/color.rs` zero-copy helpers — each block is guarded (contiguity via `as_slice()`, shape pre-checked, `[f64;3]` transmute is layout-valid) and commented. Sound. Consider `debug_assert` on pointer alignment for paranoia.
* Dependency hygiene is good (workspace deps, pinned `ndarray017` with a rationale comment, lean dev-deps). Release automation (trusted publishing, version sync across pyproject/Cargo/`__about__.py`) is exemplary.
* `attic/` (28 retired files, ~5k lines) ships in the repo and pollutes every static-analysis run (≈300 of the 1,528 dead-code findings, all scaffolding comment hits). Delete it or move it out-of-tree (git history preserves it); the `docs/plans/` and `validation/**/refs/` numba references should stay (they are the parity oracle).

---

## 8. Executive Summary

*(Sections 9–24 below extend this review: structure/persistence, interpolation, color, materials dispatch, targets, the FFI layer, edge-case numerics, the first-expansion addendum, the eigenmode tools, the second-expansion addendum, coverage gaps and their six-step closure (§19–§24 — which correct one §8 item, see §5.1.1/§17).)*

Navette is a **serious, physically-literate codebase**: the S-matrix core is correct (admittance formulation, branch-safe complex arithmetic, three properly regularized Redheffer products, hybrid coherence with a fixed cross-term bug documented in its own review history), the materials models are faithful ports verified against goldens, CIEDE2000 is term-exact, and the validation culture (376 cargo tests + 355 pytest green via the documented invocation, numba-parity mirrors, bit-identical differential tests, exposure lint) is far above average. Documentation of *conventions* — the thing that usually rots first in numerics — is consistently excellent.

The debt is concentrated and identifiable:

1. **🔴 Physics:** Névot-Croce roughness (type 5) applies the reflection factor to transmission → energy non-conservation (measured: R+T = 0.925 at a near-matched interface, σ=10 nm). One-line fix + regression test.
2. **Performance & benchmarking:** the first-pass claim "spectralweave 4–25× slower than Python" was **retracted** — it was an artifact of benchmarking a debug build (§5.1.1). In release the weaver and interpolation beat the numba reference nearly everywhere. The remaining true performance work: the `unweave_collection` batch path at high key counts (0.22–0.48×), the serial derive loop in `Solver::solve` (next Amdahl bottleneck), and — most importantly — the benchmarking/CI hygiene that allowed a debug build to pose as the reference setup (§5.1.1, §7).
3. **🟠 Process:** no push/PR CI (tests, exposure, clippy all ungated despite a README claim), which is exactly how a broken example (§6.2), 49 unsync'd request bits (§4.3), and live compiler warnings survived a 0.5.0 release.
4. **🟠 Maintainability:** the needle-gradient complex (31-param signature, 7× copy-paste demand plumbing, cyclomatic 96) and the structure CRUD brain-methods (`Layer::set_properties` 43, `Group::set_properties` 40) are the two refactor targets with real payoff; the rest of the smell report is either PyO3/attic false positives or justified hot-loop style.
5. **🟡 Hygiene:** stale `validation/README.md` and `bug_fix_plan.md` (audit trail contradicts the tree), dead exports (`MaterialSpec::evaluate`, `bake_materials`, `set_properties`…), `attic/` in-tree, serial derive loop in `Solver::solve` (next bottleneck), unclamped `DOP_R`, unguarded Sellmeier poles, 80 eV KK ceiling undocumented.

**Recommended order of attack:** (1) fix rtype-5 transmission factor + add the roughness energy-conservation test; (2) add a minimal push-CI workflow and a request-bit sync test (both are <30 lines); (3) fix/execute examples in CI; (4) profile and repair the spectralweave boundary; (5) refactor `needle_gradient` around a `NeedleRequest` struct and an enum of demand channels; (6) delete `attic/`, refresh `validation/README.md`, close out `bug_fix_plan.md`.

---

## 9. Structure, Architect, Persistence & Config

### 9.1 Two-phase expansion (`structure/expansion.rs`) — careful, bit-pinned port ✔

Phase 1 resolves bulk (group scaling, no RNG), phase 2 emits rows in traversal order with error draws in a **documented, preserved stream order** (thick → nk → rough → iface → inhg) — the detail that makes seeded Monte-Carlo runs reproducible and Python-identical. The tricky physics is all handled: interface slices are carved from the owner layer with a buffered-row rescale, mixed at f = 0.5 via Looyenga, roughness follows the plane (first row of a block), inverted-run planes inherit the traversal predecessor's flags; inhomogeneous layers subdivide into `z·f`-scaled rows; thickness floors at 0 and interface carve clamps to the available bulk. Deterministic output is pinned bit-for-bit against the Python expander (the random-stack differential test exercises it heavily). One 🔵 nit: `seed=None, apply_errors=true` silently uses thread RNG (non-reproducible) while `seed=None, errors off` builds a seeded(0) RNG that is never used — the asymmetry is correct but deserves a one-line comment.

### 9.2 Architect & structure model ✔

`architect.rs` mirrors the Python `Navette_Architect` including the **aliasing semantics**: blocks hold `Rc<RefCell<Structure>>` so edits through one block propagate to all references, and `clone_structure` exists to break aliasing — faithful to the Python object model rather than "improving" it into subtle behavior drift. Global-index addressing (`map_global`/`map_solver`/`layer_at_global`/`insert_at_global`/`split_at_global`) handles repeat counts and inverted blocks. The Python `Navette_Structure` wrapper is thin and Pythonic (`__add__`, `__iter__`, `bake_films`). Validation is warning-based (`gate_validation` → `UserWarning`) — fine, but the random-stack tests emit dozens of clamping warnings that flood pytest output; consider asserting on them or filtering in the test config.

### 9.3 Config documents ✔

Rust-first versioned envelopes (`config.rs`): `#[serde(deny_unknown_fields)]` on all sections, single canonical `PROGRAM_SCHEMA_VERSION = 1` gate, legacy flat-file fallback, prefix namespacing for multi-load collision-free merging, dependency-order assembly (materials → groups → structures → architect). JSON-only native (YAML stays Python-side for authoring) — a sensible split. 🟡 `PROGRAM_SCHEMA_VERSION = 1` is **duplicated** in `config/program.py` and `config.rs` with no sync check (same pattern as the Request bits, §4.3, though lower risk since bumping one without the other fails loudly on load).

## 10. Interpolation Module — scipy-exact ✔ with two API warts

`interpolate/mod.rs` (847 lines) implements pchip, makima, sprague (≥6 pts), Floater–Hormann (degree d), and linear, with `extrap` modes (linear/clamp/error), sorted fast-paths vs general kernels, and rayon batching (8192-point split threshold, 1024 chunk minimum). Verified numerically against scipy:

| Check | Result |
|---|---|
| PCHIP in-range vs `PchipInterpolator`, random walk, 500 pts | **max diff 7.1e-15** (rel 2e-16) |
| PCHIP interior slope formula | identical to scipy's weighted harmonic mean (`(w1+w2)·δ₁δ₂/(w1·δ2+w2·δ1)`) |
| PCHIP endpoint | identical to scipy `_edge_case` (sign guard + 3× clamp) |
| MAKIMA vs `Akima1DInterpolator(method='makima')` | **max diff 7.1e-15** |
| Monotonicity preservation on monotone data | ✔ (min diff ≥ 0 over 2000 pts) |
| Oklab/Lab/Bradford (§11) | ✔ |

Findings:

* 🟡 **Extrapolation semantics differ from scipy** and bite hard: `extrap='linear'` extends linearly from the end secants, while scipy's PCHIP extends the cubic. Porting a scipy workflow that queries outside the knot range shows enormous divergence (measured rel diff ≈ 11× on a random walk) — behavior is documented, but a note in the module docstring calling out the scipy difference would save every porter an hour.
* 🟡 **`extrap='error'` is silent NaN from Python**: the Rust kernel returns NaN + an internal flag, but the binding never raises — `eval` at an out-of-range point quietly returns `[nan]`. Either raise `PyValueError` or surface the flag; "error" mode that errors not is a footgun.
* 🔵 Derivatives for makima/FH/sprague fall back to central differences (`h = 1e-6·span`) rather than analytic forms — accurate enough for GD-style post-processing, but worth a doc line (pchip's derivative IS analytic).
* 🔵 The lib-side `as_slice().unwrap()` calls are infallible (self-built contiguous arrays) — fine, but `expect("internal arrays are contiguous")` would self-document.

## 11. Color Module — all 16 function families verified ✔

Catalog: func_01 xyY, 02 LCh, 03 Luv, 04/05 Oklab (XYZ and sRGB pipelines), 06 CIE-1964 UVW*, 07 UCS, 08 Bradford adaptation, 09–12 ΔE76/94/CMC/DIN99, 13 spectral→sRGB (SPD × CMF integration), 14 photometry (photopic/scotopic/mesopic), 15 broadcasting, 16 ΔE2000. Spot-checks via the Python API:

| Check | Measured | Expected |
|---|---|---|
| `XYZ_to_Oklab(D65 white)` | (1.0, −2.2e-5, −1.23e-4) | (1, 0, 0) ✔ |
| Oklab round-trip | 5.6e-16 | — ✔ |
| `XYZ_to_Lab(D65 white)` | (100, 0, 3e-6) | (100, 0, 0) ✔ |
| ΔE2000, Sharma test pair 1 | **2.04246** | 2.0425 ✔ |
| ΔE76, same pair | 4.00106 | 4.0012 ✔ |
| Bradford D65→A via `chromatic_adaptation_VonKries` | matches first-principles matrix·xyz to **1e-9** | ✔ |
| Photometry constants | Kmₚ = 683.002, Kmₛ = 1700.05 | standard ✔ |

* 🔵 **Convention trap**: `calc_transform_matrix` returns the matrix in **row-vector convention** (the transpose of the textbook column-convention Bradford matrix; applied internally via `vec_mul_mat3(v, M)`). Internally consistent and correct — but anyone extracting the matrix and applying `M @ xyz` in numpy (as this reviewer first did) gets a silently wrong answer. One doc line on the binding fixes it.
* 🔵 Mesopic is a caller-supplied linear blend (`m`), not the CIE 191 MOVE model — the docstring is upfront; fine for the film-colorimetry use case.
* The CIE tables are CI-verified against source data (`tools/check_cie_sync.py`) and the module carries a golden matrix generated from `colour-science` — the strongest parity story in the repo.

## 12. Materials Dispatch & Synthesis Targets ✔

`materials/__init__.py::evaluate` is a flat if-elif over ~20 native models (Konstant, Table, Cauchy[+Urbach], Sellmeier[+Urbach], Lorentz, Drude[+Lorentz], CodyLorentz, Forouhi-Bloomer ×3, Tauc-Lorentz, UBF, six EMA composites, Roughness) with a tidy `req()` helper producing precise missing-param errors, recursive `_eval_nested` for EMA host/inclusion, and an explicit, helpful refusal of non-linear Table interpolation ("resample the table offline") rather than silent misbehavior. Cyclomatic count is high but the shape is table-like and readable — leave it. `spectralweave/target.py` defines the target DSL: `SpectralTarget`/`AngularTarget`/`ColorTarget` dataclasses with kinds `e`(exact)/`a`(above)/`b`(below)/`r`(range box)/`c`(soft center), per-point `band`, relative tolerances, `weight` × `normalize_count` (sampling-density fairness — the 200-pt vs 10-pt weighting question, solved), and `integral` mode. The merit-shape documentation (each kind's contribution formula) is precise enough to reimplement from the docstring — the standard this codebase hits when it is at its best.

## 13. FFI / Binding-Layer Audit ✔ (with two fragility flags)

* **GIL discipline is good**: 48 `py.detach` sites across the bindings — the solver, materials, spectralweave, and synthesis long-running calls all release the GIL; the synthesis driver detaches for the whole pipeline and **re-attaches (`Python::attach`) only per user callback**, converting a Python exception into the documented `USER_ABORT` termination. Textbook PyO3.
* **Panic hygiene is good**: the ~350 raw `unwrap/expect/panic` hits are overwhelmingly in `#[cfg(test)]` modules (merit.rs 1008+, needle_pass 930+, color_merit 723+, solver 1811+, needle_operator 1385+ are test-line offsets). Library code holds ~17 across the hot files, mostly infallible internal invariants. PyO3 converts residual panics to `PanicException` (no UB), but the handful of lib unwraps should still carry `expect("invariant: …")` messages.
* **Float determinism under rayon holds**: parallelism in the solver/merit paths is per-element maps or per-thread output buffers indexed by destination (deterministic write addresses); the only float reductions found are sequential (`iter().sum()`), so run-to-run bit-identity — which the differential tests pin — is not at risk from scheduling order. Worth keeping as an invariant note: any future parallel *reduction* (e.g. a parallel merit sum) would introduce last-ULP nondeterminism into a suite that asserts bit-identity.
* 🟡 **`unweave_collection` pointer dance is correct but fragile** (`spectralweave_optical.rs:435`): numpy buffer addresses are captured as `usize`, the GIL is released, and slices are reconstructed via `from_raw_parts` — kept valid only because the `py_arrays` Vec (never read again) pins the Python objects across the detached section. Nothing marks that Vec as load-bearing; a refactor that "cleans it up" produces use-after-free. Add a comment naming it as the lifetime anchor (and consider `PyReadonlyArray::as_slice` + holding `Bound` refs instead of raw addresses).
* 🔵 The `maturin develop` ImportError fallbacks (materials/interpolate/color/spectralweave wrappers) all give the same helpful build hint — nice consistent touch.

## 14. Edge-Case Numerics (experiments run for this review)

| Experiment | Result |
|---|---|
| TIR (glass 1.5 → air, 43°–89°) | Rs = Rp = 1.000000 exactly, Ts = Tp = 0; Ψ → 45°, Δ continuous through the critical angle ✔ |
| Near-Brewster angles | Rs/Rp split consistent with Fresnel ✔ |
| Absorbing slab (k = 0.05), d = 0 → 100 µm | A = 1−R−T **monotone increasing**, saturates at the incoherent limit (0.9596); A ∈ [0,1] ✔ (d = 0 shows −0.0000 round-off — unclamped, cosmetic) |
| Interpolation vs scipy (in-range) | 7e-15 ✔ (§10) |
| Colorimetry vs references | all pass (§11) |
| Roughness type 5 energy | **fails** — R+T = 0.9247 (§3.2) ✗ |

---

## 15. Addendum to the Executive Summary (post-expansion findings)

The expansion pass (structure/architect/config, interpolation, full color catalog, materials dispatch, target DSL, binding layer, edge-case numerics) **raises the overall grade**: every module reviewed after the first report turned out to be correct, scipy/colour-science-exact where a reference exists, and honestly documented. New findings, merged into the priority list:

* **🟠 `ScatterMatrix.energy_conservation()` can never work**: it hardcodes `squeeze=False` (always returns 2-D `[n_angles, n_wavs]`) but the native `solver_energy_conservation` binding only accepts 1-D arrays → guaranteed `TypeError` on every call. The method is advertised by the README's energy-conservation feature and has zero test coverage. Fix: `squeeze=True` in the wrapper (or accept 2-D in the binding).
* **🟠 Parity tests are silently excluded from the documented run**: `validation/conftest.py` sets `collect_ignore = ["parity", "benches"]`, so `pytest validation` runs 355 tests while 11 pytest-style parity tests (the strongest oracle in the repo) never execute. Worse, the two un-ported files (`test_core_engine_photometry_only.py`, `test_core_engine_rigorous_ellipsometry.py` — still marked NEEDS PORT) call `sys.exit(1)` at module import when targeted directly, crashing pytest with an INTERNALERROR instead of skipping. Replace `sys.exit` with `pytest.skip`, narrow the ignore list to the true script files, and include the live parity tests in the default run.
* **🟡** interpolation `extrap='error'` silently returns NaN (§10); extrapolation semantics vs scipy deserve a doc warning (§10).
* **🟡** `unweave_collection` lifetime-pin fragility (§13); duplicated `PROGRAM_SCHEMA_VERSION` (§9.3).
* **🔵** row-vector convention on `calc_transform_matrix` (§11); clamping-warning flood in tests (§9.2); lib unwraps without invariant messages (§13).

Everything else examined in this pass is verified correct: bit-parity expansion with disciplined RNG streams, faithful aliasing semantics, scipy-exact interpolation, textbook-correct colorimetry, correct TIR/absorption energetics, and a binding layer whose GIL and panic discipline most production PyO3 projects don't match. The project's weaknesses remain exactly where §8 located them — the rtype-5 roughness factor, the spectralweave performance regression, CI/test-gating gaps — and are joined only by small, easily-fixed API warts.


---

## 16. Eigenmode Tools (landscape / refine / field profile) — one robustness bug

The guided-mode feature set (`ScatterMatrix.eigenmode_landscape`, `refine_mode`, `find_eigenmodes`, `field_profile`; Rust: `solver.rs` + `optimizer.rs`) scans `|1/r(n_eff)|²` over a complex effective-index box, finds grid local minima (true median threshold — the docstring even records *why* median was chosen over mean), polishes seeds with an unbounded 2-D Nelder–Mead, and extracts `|E(z)|` profiles. Verified numerically against physics:

| Check | Result |
|---|---|
| p-pol SPP of a 50 nm metal film (ε = −12+1j) between air and glass | landscape dip val 1.2e-4; `refine_mode` from seeds (1.45, 0.2) and (1.3, 0.1) both converge to **1.0459458 + 0.0015949j** (LR-SPP), characteristic value ≈ 4e-19 ✔ |
| `field_profile` of the converged mode | |E| peaks **exactly at the film/substrate interface** (z = 50 nm) — textbook LR-SPP ✔ |
| s-pol (no mode in box) | landscape correctly flat (min ≈ 0.99) — **but `refine_mode` diverges**: returned n_eff = −1.7e8 + 7.4e6j with val = 1.7e-217 ✗ |

The runaway mechanism: `char_func` guards the **lower** edge (`|r| < 1e-15 → 1e30`, protecting Brewster zeros from becoming fake poles) but has no upper guard — at large |n_eff| the propagation exponent overflows, `|r| → ∞`, and `|1/r|² → 0`, which the unbounded Nelder–Mead reads as a *perfect* eigenmode and sprints toward. Numerically the returned `val ≈ 1e-217` is indistinguishable from success. The primary workflow is partially protected (find_eigenmodes seeds come from the bounded grid), but a manual `refine_mode` call with a poor guess returns silent garbage with a confident-looking value. **Fix:** clamp/reject simplex points outside a physical box (e.g. `|n| ≤ 3·max|n_layer|`), and treat non-finite/overflowed `char_func` values as `+inf` for the minimizer. Asymmetric with the existing lower guard — same pattern as §3.3's DOP_R asymmetry and §3.2's roughness t-factor: one-sided invariants are this codebase's recurring failure shape.

---

## 17. Final Addendum (second expansion) — corrections & new findings

The second expansion pass (needle engine, eigenmode tools, remaining synthesis internals, packaging, benches re-run in release) **corrects the single biggest claim of the first report and adds one genuine physics-adjacent bug**:

* **🔴→✔ RETRACTED**: "the spectralweave data-plane is 4–25× slower than Python" — an artifact of benchmarking a **debug** build against optimized numba (§5.1.1). With the documented-once-only `maturin develop --release`, Rust wins 1.4–8× across the weaver (except `unweave_collection` at high key counts, 0.22–0.48×) and 1.0–1.3× on interpolation (pchip 1.22 ms vs loom 1.59 ms at 1M pts). The README setup instruction (and every committed bench conclusion made from debug builds) is the finding, not the engine.
* **🟠 NEW** `refine_mode` returns silent garbage (n_eff ≈ −1.7e8, val ≈ 1e-217) when no pole exists in range — unbounded Nelder–Mead + one-sided guard in `char_func` (§16).
* **🟠 NEW** packaging: `numpy>=1.22.0` floor incompatible with the numpy-2-built abi3 wheels (§7).
* **🟡 NEW** `1dinterpol_test_bench.py` crashes on stock Windows consoles (emoji in print + cp1252 `UnicodeEncodeError`) — and the target/color benches carry the same emoji prints; guard with `sys.stdout.reconfigure(encoding="utf-8")` or drop the emoji. (The bench runs fine with `PYTHONIOENCODING=utf-8`.)
* ✔ Verified good in this pass: needle engine request-bit architecture (18 NREQ channels, single-sweep per point), median-threshold minima search, field-profile extraction, color bench exact parity incl. edge rows, release-build performance across the board (backside full mask 0.164 ms ≈ 1.4 µs/point all-observables).

**Revised recommended order of attack:** (1) fix rtype-5 transmission factor + roughness energy test (§3.2 — unchanged, the physics bug stands); (2) fix `energy_conservation()` wrapper (§15) and bound `refine_mode` (§16) — both one-liners; (3) README → `maturin develop --release`, bench scripts assert non-debug + UTF-8-safe output; (4) minimal push-CI + request-bit sync test + narrow the parity `collect_ignore` (§15); (5) numpy>=2.0 floor; (6) fix/execute examples in CI; (7) optimize `unweave_collection` batch path; (8) refactor `needle_gradient` (§4.1); (9) delete `attic/`, refresh stale docs.

---

## 18. Remaining Coverage Gaps — what this review did NOT fully examine

Honest inventory after two expansion passes. Approximately **14,000 of 36,500 Rust lines** and **~2,100 Python lines** were never read line-by-line (they were assessed only through smells, module docs, benchmarks, and test outcomes).

### 18.1 Unread code (by expected impact)

| Area | Files (LOC) | Status |
|---|---|---|
| **Needle math core** | `needle_operator.rs` (2,479) | Header/docstrings read; the gradient algebra (P(z), `dr/ddz`, four-point slopes, dispersion chains) never term-verified — only loom-parity-pinned |
| **Largest binding** | `navette-py/structure.rs` (2,231) | Smells assessed (brain methods), binding logic (py↔native conversion, defaults, validation) not read |
| **Synthesis core** | `needle_pass.rs` (2,128), `merit.rs` (1,931), `color_merit.rs` (1,349), `targets.rs` (1,083), `structure.rs` (694), `design_config.rs` (617), `evaluator.rs` (545) ≈ 8,300 LOC | Architecture + smells + bench smoke only. The target-kind fold algebra (e/a/b/r/c → clamp/quadratic forms), integral demands, and **color merit (illuminant/observer handling — the largest unread feature)** never term-checked |
| **Pipeline stages** | `inflate.rs` (426), `stagnation.rs` (338), `cleanup.rs` (320), `driver.rs` (274), `context.rs`/`config.rs` | Headers only; the documented "clamp before AND after re-optimize" ordering unverified |
| **Structure model** | `weaver.rs` (448), `group.rs` (581), `providers.rs` (415), `layer.rs` (382), `enums.rs`, `validation.rs`, `version.rs` ≈ 2,300 LOC | The brain-method smell counts came from here; the actual semantics (group error draws, provider resolution, structure weaving/merge) not read |
| **Bindings** | `materials.rs` (414), `spectralweave_target.rs` (426), `synthesis_merit.rs` (314), `config.rs` (166), `lib.rs` | Not read |
| **Interpolation kernels** | `eval_fh`/`calc_fh_weights`, `eval_sprague_*`, `robust` mode | Only pchip/makima were scipy-verified; FH weights, Sprague 5th-order, and the `robust` guarded mode are parity-pinned only |
| **Python surface** | `needle.py` (237, docstring only), `architect.py` (289), `config/*` (≈612), `target.py` (383, header only), `optical.py` (133, partial), `data/__init__` (31) | Thin wrappers by all evidence, but unverified |

### 18.2 Physics verified only by parity (inherited-bug risk)

The rtype-5 finding (§3.2) proved a reference can be wrong; everything below shares that exposure:

* **Needle gradients**: no independent finite-difference check of P(z), the T/A/φ channels, or the dispersion sensitivity chain (dφ/dgd/dgdd/dtod/dfod) was performed by this review.
* **LM optimizer**: the plan docs claim scipy parity; this review never reproduced it.
* **Material kernels not read or independently verified**: Cody-Lorentz, the three Forouhi–Bloomer variants, Tauc-Lorentz, UBF — goldens exist but are generated *from the reference* (same trap as §3.2).
* **Color formulas not term-verified**: ΔE94, CMC, DIN99, CIE-64 UVW*, func_13 SPD→sRGB (only ΔE76/2000, Oklab, Lab, sRGB-white, Bradford were independently checked).
* **Roughness types 1–4**: spot-checked at normal incidence on a single interface only — no angular/multi-interface sweep against an analytic graded-layer reference.
* **Multiblock cascade (P_MB channels)** and **hybrid Modes A/B/C**: never exercised against hand-computed cases.
* **Eigenmodes**: one SPP scenario verified; no multi-mode stacks, no tolerance/robustness sweep of `median_factor`.

### 18.3 Runtime dynamics never tested

* **`OpticalWeaver` thread-safety** — the docstring advertises concurrent setups; no race test was run (Arc/parking_lot structure looks right, but untested).
* LRU `cache_size` eviction semantics; `n_stack_cache` invalidation on mutation.
* **Garbage-in audit**: NaN/inf indices, negative/absurd thicknesses, sinθ > 1, zero-point masks, empty wavelength grids — behavior unspecified and untested (the eigenmode §16 finding shows how these become silent garbage).
* Monte-Carlo error-path statistics (`apply_errors` — the header asserts statistical agreement with the Python expander; never run).
* Pickle round-trip of `UniInterpolator` (`__reduce__` exists, untested).

### 18.4 Process/docs/packaging gaps (discovered while compiling this section)

* **License is present but under GNU naming**: `COPYING` (GPL-3.0) + `COPYING.LESSER` at the repo root, with `license = "LGPL-3.0-or-later"` declared consistently in both `Cargo.toml` and `pyproject.toml`. Only nit: not all source files carry the SPDX header (e.g. `smatrix.py`, `solver.rs` have none) — fine legally, inconsistent stylistically.
* `docs/` holds ~30 more documents never inventoried (`materials-architecture.md`, `materials-dispersion-usage.md`, `materials-implementation.md`, per-function color plans, synthesis plans) — doc-vs-code drift unaudited.
* `tools/check_cie_sync.py` was run (passes) but not read — its source-comparison logic is unverified.
* `.pyi` stubs never cross-checked against the actual binding signatures (drift audit missing).
* The README performance table was written against (unknown-profile) builds and not re-audited after the release rebuild — with the corrected numbers (§5.2) some claims may now be *understated*.
* Committed `backside_speed_*.json` results: build-profile provenance unknown.
* release.yml: macOS wheel architecture (universal2/arm64) not examined; sdist contents not inspected.
* The `attic/` tree is import-isolated (verified: no live imports), but its 28 files were never diffed against their live replacements.
* Git history/branch hygiene not examined.

### 18.5 If the review continues, in order of expected information gain

*(Worked through in §19–§24; status updated in place.)*

1. ~~**Term-verify `needle_operator.rs` + `merit.rs`/`needle_pass.rs` fold algebra with independent finite differences**~~ → **DONE, all green (§19)**.
2. **`color_merit.rs`** — the largest unread feature (1,349 lines), and color synthesis is a headline use case. → §20.
3. **Garbage-in audit** (NaN/inf/sinθ>1) across `ScatterMatrix`, `UniInterpolator`, `OpticalWeaver` — cheap, and proven valuable by §16. → §21.
4. **Threaded race test of `OpticalWeaver`** + LRU eviction test. → §22.
5. `color_merit`-adjacent: regenerate one material golden from first principles (break the reference-oracle dependence for one model). → §23.
6. `.pyi` drift audit + docs inventory. → §24.

---

## 19. Independent Verification of the Synthesis Gradient Chain ✔ (one parity gap found)

Step 1 of §18.5, executed. Method: everything below was verified against **finite differences through the real solver and hand-rolled Python reimplementations** — no loom, no parity reference. Harness: `validation/review/fd_step1.py` + `fd_rchannel.py` (kept for reproduction). A 5-layer roughness-free stack (air | TiO₂ 120 | SiO₂ 200 | Ta₂O₅ 80 | glass, 31 λ × 2 angles) and a mixed spec (kinds e / a / r with band, weights 0.9–1.7, linear + log + phase transforms) were used.

### 19.1 Results — all verified

| Check | Result |
|---|---|
| **Fold evaluator** (Part A): my independent Python reimplementation of the merit fold (kind_residual per e/a/b/r/c, transforms linear/log/phase with residual wrapping, rscale = √(weight/count_norm), integral-mean branch) vs native `MeritSpec.merit` | **2e-16** — the evaluator is exactly what its docstrings claim |
| **End-to-end gradient** (Part B): analytic dF/dθₖ from `build_needle_targets` fold + `needle_gradient` (Σ 2·(P_R+P_T+P_A) + P_φ) vs central-difference FD of the native merit through the real solver, for each of the 3 film layers | **rel 2e-8 – 7e-7** across all layers — the fold→needle chain (including band-activated per-point weights: only 33/62 active on the 'r' target) is coherent end to end |
| **Dispersion ladder** (Part C): `dphi` vs physical FD of the phase (unwrap-free Im{conj(a)·∂a/∂δ}/\|a\|²) | clean **O(h²) convergence** to the engine: 1.2e-1 (h=8nm) → 3.9e-4 (h=0.5nm), rate 4.1; both s- and p-pol verified |
| **Spectral chain**: dgd/dgdd/dtod/dfod vs independent `np.gradient` in ω applied to the FD-verified dphi | **0.0 (bit-exact)** at every order — the engine's `grad_nonuniform` is a plain 2nd-order central difference and correctly chained |
| **P(z) constancy** (needle material = host material ⇒ inserting host = growing the layer) | flat to **8e-16** inside a layer |
| **R-channel slope with a distinct needle material** — *never FD-verified anywhere in the repo* (see §19.2) | O(h²) convergence confirmed: 5.9e-3 (h=2nm) → 2.2e-5 (h=0.125nm), rate 4.0 ✔ |

The needle-slope algebra itself (ρ̂ = −2iβ′r₁₂/(1−r₁₂²), τ̂ = iβ′(1+r₁₂²)/(1−r₁₂²)) was additionally re-derived analytically and the dual-number composition (`CDual`/`star_dual`) checked term-by-term against a hand-written Redheffer star — exact. **The synthesis gradient chain is correct, independently of the reference.**

### 19.2 Finding: the parity suite never FD-checks the R channel

`test_needle_t_a_phi.py`'s `check_fd` exercises only **T** and **A** (and φ via wiring-to-zero tests); the pure-reflectance needle channel (`NeedleRequest.P`) appears only in (a) wiring tests where P ≈ 0 by construction and (b) the loom mirror — which shares bugs with the reference by design. My §19.1 Part B inadvertently inherited the same gap (host-material needles make P ≡ 0, so the aggregate matched through the T channel alone). The supplementary check (`fd_rchannel.py`) closes it: P is correct — but a one-line addition to the parity suite (an R-row in `check_fd` with a distinct needle material) would make this permanent. **Recommend adding.**

### 19.3 Incidental discoveries

* **🟠 Negative layer thicknesses are silently clamped** (feeds §21): `ScatterMatrix([…], […, −5.0, …])` is accepted at construction and the solver's `d > 1e-12` guard makes the layer vanish — `r(δ=−5nm slab)` is bit-identical to the no-slab stack, **no warning**. Any user typo (sign error in a thickness array) silently changes physics. Recommend validating `thickness >= 0` at construction (or warning).
* **Sub-nm positive thicknesses are fine** (δ = 0.001nm is honored) — only the negative side is clamped.
* **🟡 Stale docstring**: `needle_slopes` documents "τ̂ = 2iβ′/(1−r₁₂²)" but the code (correctly) computes τ̂ = iβ′(1+r₁₂²)/(1−r₁₂²) — verified analytically. One-line doc fix.
* **🟡 Harness lesson (for future benches)**: FD across δ must use positive-only stencils — the clamping in the first bullet turns any ±h stencil into a half-slope, silently.

---

## 20. `color_merit.rs` — Color Synthesis Merit ✔ (verified exact; one Python-API gap)

Step 2 of §18.5: the full 1,349-line color-demand kernel, read and then verified numerically against independent reimplementations (`validation/review/color_merit_check.py`). Setup: D65/1931-2° demand on the solver's Rₛ row (450–750 nm, 31 pts, 30°), weight 1.7.

### 20.1 Verified (all exact)

| Check | Result |
|---|---|
| **Lab + ΔE2000 merit value** vs my own chain (rectangular XYZ integral with k-normalization, native-grid white, my Lab, my Sharma-ΔE2000) | **2.7e-16** — the demand's own white correctly integrates on the illuminant's **native 1 nm grid** (300–830), not the sim grid; the sample integral uses the covered sim subset with forward-difference Δλ |
| **End-to-end color gradient** dF/dθₖ through the real solver (my chain rule: my gᵢ = d(residual)/dRᵢ from FD, times dRᵢ/dθ from the needle R channel) | **1.6e-7 / 2.8e-8 / 1.6e-8** for the three layers — the deposit mathematics (analytic dXYZ/dR = E·cmf·k·Δλ, exact-linear, chained with a 3-pt FD of the XYZ→objective map) is sound |
| **CIE whiteness pair** (W = 100Y + 800(xₙ−x) + 1700(yₙ−y), Tw = 1000(xₙ−x) − 650(yₙ−y)) vs closed form | exact |
| **E313 yellowness** (100(CxX − CzZ)/Y, explicit coefficients — schema defaults D65/10°) | exact |
| **Dominant wavelength + purity** vs my own locus-ray intersection (forward purity 1/t, complementary branch negative, achromatic (0,0) kink) | exact (dl = 567.5475, dp = 0.5969) |
| **XyY / Oklab / Luv Channels** (Oklab: Bradford-adapt to D65 then the §11-verified map; LCh hue wrapped ±180° before scaling) | 1.5e-16 / 2.9e-10 (my Bradford uses a numeric inverse — expected) / 7.8e-15 |

Design choices verified as sensible: ΔE refused for non-Lab quantities (Oklab-ΔE would double-count; equal-tol Channels is the algebraic twin); sRGB deliberately **unclipped** (linear extension keeps gradients honest); `wavelength_range` windows the sample integral only (white/locus stay full-range, keeping numbers comparable); partial table overlap narrows silently, empty overlap errors — all documented in-module.

### 20.2 🟡 Finding: color gradients are invisible on the Python needle path

The Rust fold computes dedicated color buckets (`grad_r`/`grad_t` — chain-rule deposits per covered sim point) but the **Python dict returned by `build_needle_targets` omits them** (only `targets`/`weights` per channel surface), and `needle_gradient` has no grad-bucket inputs. Consequences:

* The production path is unaffected — `run_design` consumes the buckets internally (cargo tests pin the deposit semantics, including the U-curve ÷2 rule).
* But the **documented Python flow** (fold → `needle_gradient`, shown in the `synthesis.__init__` docstring) silently yields **zero color gradient**: a user assembling custom needle cycles in Python with color demands optimizes against pointwise targets only, with no error. Recommend exposing `grad_r`/`grad_t` in the fold dict and matching `grads_r`/`grads_t` inputs on the native `needle_gradient` binding.
* My Part C check assembles the chain rule by hand to close this gap for verification purposes; the numbers confirm the hidden deposits are correct.

### 20.3 Read-only observations

* Input validation in `ColorDemand::new` is thorough (reference-shape gates per quantity, compat matrix, grid monotonicity/finiteness, weight ≥ 0, window lo<hi, degenerate-locus/illuminant refusals) — consistent with the rest of the codebase.
* The objective-map gradient step is a **central FD** (h = 1e-6·(1+\|XYZ\|)) rather than analytic — appropriate for the non-smooth maps (DomWl ray intersection, hue wrap) where an analytic derivative would be wrong at kinks; cost is 6 tiny map evals per demand.
* Coverage-rule asymmetry (partial overlap narrows silently vs pointwise targets which interpolate) is documented — a porter trap worth a README line.

---

## 21. Garbage-in Audit — NaN/inf/sinθ>1 (§18.5 step 3)

Systematic pass across `ScatterMatrix`, `needle_gradient`, `UniInterpolator` (`validation/review/garbage_in.py`). Results sorted by severity:

### 21.1 Findings

* **🟠 Duplicate wavelengths poison the dispersion channels silently.** A grid `[550, 550, …]` (a merge/concat bug away in any user code) produces **NaN in GD/GDD/TOD/FOD** with no warning — the spectral differentiation divides by Δω = 0. The color-demand path validates strictly-increasing grids (`check_grid`), but `ScatterMatrix` construction does not. One-line fix: validate strictly increasing wavelengths at construction (and thicknesses — see next).
* **🟠 Negative or NaN thickness silently deletes the layer.** `d > 1e-12` is false for both, so the layer vanishes from the propagation sweep — bit-identical results to removing it (verified), no warning (§19.3). A NaN thickness is indistinguishable from intent.
* **🟡 Angles > 90° silently alias to their mirror.** cos is reconstructed as √(1−sin²θ) (sign discarded), so R(120°) ≡ R(60°) (verified equal to 1e-10) — physically self-consistent (R+T = 1 still holds at the aliased angle) but almost certainly not what a caller passing 120° meant. Construction-time validation should reject θ ∉ [0°, 90°].
* **🟡 `needle_gradient` z out of range is debug-only.** `debug_assert!(xi <= ds[j])` in `needle_slopes4_ddz` means release builds silently evaluate the needle slope with ξ outside the host layer (z=1000 on a 400nm stack, or negative z — both "silent-clean" garbage). The assert must be a release `assert!`/`Err` — it's a pure user-input check, no hot-loop cost.
* **🟡 NaN/inf indices and NaN needle indices propagate NaN** with no construction-time validation — standard GIGO behavior (acceptable), but combined with the above, `ScatterMatrix` has **no input validation at all**: everything (NaN indices, NaN/duplicate/descending wavelengths, negative/NaN thickness, out-of-range angles) is accepted silently. All of it is checkable in O(n) at construction.
* ✔ **Handled well**: huge thickness (1e9 nm — no overflow, finite output); single-wavelength grids with dispersion requests (finite zeros, no division error); descending wavelength grids for dispersion (the central-difference quotient is sign-invariant under reversal — algebraically correct, verified finite and consistent); `UniInterpolator` rejects duplicate/unsorted knots with a clean `ValueError` at construction and propagates NaN inputs as NaN outputs.
* ✔ TIR-adjacent needle inputs (needle index below the critical floor) compute evanescent-branch sensitivities without NaN — undocumented but numerically stable.

### 21.2 Recommendation

A single `_validate_stack(n, d, wavelengths, angles)` at `ScatterMatrix` construction (finite n/d, strictly increasing unique wavelengths, 0 ≤ d, 0° ≤ θ ≤ 90°) plus promoting the two needle `debug_assert!`s would close every silent-garbage path found. Cost: microseconds per call; the construction path is not hot.

---

## 22. OpticalWeaver Concurrency & LRU Eviction ✔ (§18.5 step 4)

The docstring advertises "optimized for concurrent simulation loops" — tested (`validation/review/weaver_race.py`):

* **LRU plan-cache eviction**: `cache_size=2` with 4 tiling frames and 10 keys (deliberate thrash) — every key reassembles **bit-exactly** (bad = 0/10). Evicted plans rebuild correctly; data frames are independent of the plan cache.
* **8-thread race**: 320 mixed operations (`unweave`, `unweave_collection`, `get_weaved`, `get_weaved_collections`, `invalidate_cache`) against one weaver with cache_size=2 — **no exceptions, no torn writes**: all 192 single-writer keys reassemble exactly to their writer's curve (0 mismatches). The `RwLock`-protected frame map + `Arc`-shared frames + per-key fragment writes hold under contention.
* One behavioral note: `get_weaved` on a not-yet-written key raises `ValueError: Key not found.` — correct strictness, but a "contains_key"-style probe would help polling consumers (minor API polish).
* Harness note: concurrent writers to the *same* key are not transactional across the 4 frame writes (per-frame locking only) — last-writer wins per fragment, so same-key concurrent writes can interleave. The intended usage (one key per simulation step) never hits this, but it belongs in the concurrency docs.

---

## 23. Material Kernel Golden from First Principles: Tauc–Lorentz ✔ (§18.5 step 5)

The reference-oracle dependence is broken for one KK-family model. `TaucLorentz` (Jellison–Modine 1996) was verified entirely from first principles (`validation/review/tauc_check.py`, `kk_validate.py`, `kk_conv.py`) — no loom, no goldens:

| Check | Result |
|---|---|
| **ε₂ closed form** (A·E₀·C·(E−Eg)²/(((E²−E₀²)²+C²E²)·E), 2 oscillators, 10 energies incl. below-gap) | engine == my independent formula to **1e-12** (scaled); below-gap ε₂ ≡ 0 exact |
| **ε₁ via Kramers–Kronig** — my independent **pair-sampled principal-value quadrature** (ξ = E±(k+½)h pairing cancels the PV singular part exactly; no exclusion window), validated against an **analytic Lorentz KK pair** first | agreement **0.07–1.1 %** at all 10 energies (worst at the sharp E₀=4 eV resonance); my quadrature itself converges (h-refinement stable to 9e-5) |
| **Residual ~1 %** near resonances does **not** close with quadrature refinement on either side | consistent with the module docstring: ε₁ is deliberately the **FFT-KK on the 0.01–80 eV grid** ("for consistency" across KK models), whose near-resonance accuracy is ~1 % — a documented design tradeoff, not an error. The closed-form Jellison–Modine ε₁ is explicitly noted as not implemented. |
| n̂ = √ε, below-gap k | exact |

Methodological note: my first two PV schemes (fixed exclusion window, naive pair masking that dropped the + branch for u > E) were themselves wrong by 20–55 % — the analytic-Lorentz validation step is what caught *my* harness, which is exactly why the oracle chain (analytic case → quadrature → engine) matters.

**Code-level read of the sibling kernels** (`cody_lorentz.rs`, `ubf.rs`, `forouhi_bloomer.rs` — first read of these files):
* All KK-family models share the same `kk::kk_fft` path (documented in each module header) — so this first-principles validation of the KK machinery transfers to Cody-Lorentz, UBF, and FB at the transform level; their closed-form ε₂ generators remain parity-pinned only.
* Cody-Lorentz ε₂: band region above Eₜ with the Cody factor (E−Eg)²/((E−Eg)²+Ep²), C⁰-matched Urbach tail below Eₜ (amplitude summed over oscillators at Eₜ) — the Ferlauto-style construction, with overflow-safe softplus in UBF (x>50 / x<−50 guards) and exact γ∈{0.5,1,2} fast paths. No red flags spotted; the guard style is consistent with the codebase's one-sided-invariant pattern (worth the same FD spot-check as §19 someday, low priority).
* The `forouhi_bloomer.rs` FB term (25–67) implements the standard FB interband ε₂ = (E−Eg)²/(E²·((E²−C)²+B²E²)) shape via `fb_term` with the metal variant dropping the (E−Eg)² numerator — matches the Ferouhi-Bloomer 1986 forms at read-through.

**Updated coverage statement**: the earlier §18.2 worry ("goldens are generated from the reference — same trap as rtype-5") is now resolved for the KK transform itself: the transform is verified against independent physics, and it is shared by all four remaining models. What remains parity-only: the closed-form ε₂ shapes of Cody-Lorentz/FB/UBF (low risk — each is ~15 lines of algebra against a textbook formula).

---

## 24. `.pyi` Drift Audit & Docs Inventory ✔ (§18.5 step 6)

### 24.1 Stub drift: none

The repo ships exactly one stub file, `_materials.pyi`. Extracting the actual PyO3 registration list (`wrap_pyfunction!` in `navette-py/src/materials.rs`) and diffing against the stub's `def`s: **22 = 22, exact match** — every binding function is stubbed, no extras, no signature drift detected (param-name spot-checks align). 

Gap of a different kind: `_materials` is the **only** stubbed module — the entire `_smatrix` surface (49 request bits, `needle_gradient`, `ScatterMatrix`, eigenmode methods, synthesis pipeline, `SimCurves`/`MeritSpec`) and `_spectralweave` (weaver) are unstubbed, so IDE support and mypy coverage stop at materials. Low-risk (the Python wrappers carry their own type hints) but inconsistent.

### 24.2 Docs inventory (39 files under `docs/`)

| Group | Files | State |
|---|---|---|
| Materials | `materials-architecture.md`, `materials-dispersion-usage.md`, `materials-implementation.md` | authoritative; headers still say "Loom" (the project's internal name — appears across several docs; cosmetic naming drift vs the public "Navette") |
| Color port records | 16 × `plans/color/func_NN_*.md` + `plans/color/README.md` | per-function implementation task records — historical, consistent with the shipped code as verified in §11 |
| Synthesis/smatrix plans | `smatrix-SYNTHESIS_PLAN.md`, `plans/smatrix/REVIEW.md`, 6 × `unit_task_func_N.md`, `color_targets_optionB_plan.md`, `color_targets_proposal.md`, `needle_gradients_plan.md`, `gradient_layers_plan.md`, `multi_environment_plan.md`, `spectralweave-optimization-notes.md` | design/plan records; `plans/smatrix/REVIEW.md` documents previously-fixed bugs (consistent with §3 cross-references) |
| Structure | `plans/structure_plan.md`, `plans/rust_first_plan.md`, `plans/rust_port_plan.md`, `plans/bug_fix_plan.md`, `plans/exposure_audit.md` | `structure_plan.md` is explicitly marked **SUPERSEDED** (good practice); `bug_fix_plan.md`'s two items are stale-as-undone — the checkbox style without completion markers invites the exact misreading caught in §9 (plan docs read as open bugs) |
| Public contract | `spectralweave-target-kinds.md` | the fold/convention contract — verified accurate by §19/§20 numerics |

### 24.3 Status after steps 1–6 — updated review stance

Every §18.5 item is now closed: the synthesis gradient chain (§19), color merit (§20), garbage-in behavior (§21), weaver concurrency (§22), the first-principles material golden (§23), and stubs/docs (§24) all came back green except the specific, small findings already logged (silent clamps/aliasing §21, Python color-gradient gap §20.2, debug-only needle asserts §21.1, stale `needle_slopes` docstring §19.3). The two confirmed physics/API bugs from the first pass (rtype-5 transmission factor §3.2; `energy_conservation()` wrapper §15) remain the only defects that change numerical results, and the remediation list in §17 stands.
