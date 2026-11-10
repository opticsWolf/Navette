# Synthesis targets: polarization, curve coverage, weighting, coherence

Status: **in progress.** W1, W2 and W3 are implemented (2026-11-10);
W4-W8 are open. Drafted 2026-11-09.

This plan follows a review of the needle / thin-film synthesis path
(`smatrix/synthesis/*`) against five reported issues:

| # | Report | Verdict |
|---|---|---|
| a | Target weighting is scalar only; no weighting profile over wavelength or angle | **Partly true**: per-point tolerances already act as inverse weights, but they are coupled to band semantics and lost in integral targets |
| b | The evaluator produces only front amplitudes; the fold does not map polarized targets to their own polarization | **True, both parts** |
| c | Many target types are not usable in synthesis | **True**: unpolarized, back-side, front-reflection phase, colour on unpolarized curves, log transform in the needle pass |
| d | Materials with the same index are deduplicated | **True, by design**; rationale sound, several sharp edges |
| e | Incoherent layers are refused in synthesis | **True, by design** (a refusal that replaced a silent wrong answer); the engine already has the operator to lift it |

Workstreams W1-W8 below fix or close each item. W1-W3 are bugs and come
first. Everything is implemented in Rust; Python only wraps the Rust
doors (no target, weight or validation logic in Python).

---

## 1. Findings (verified)

### 1.1 The needle fold merges polarizations (b)

`build_needle_targets_env` (`needle_pass.rs:359-381`) chooses a bucket by
**quantity only**:

```text
Rs | Rp | Ru -> r      Ts | Tp | Tu -> t      As | Ap | Au -> a
RBs| RBp| RBu-> rb     TBs| TBp| TBu-> tb     ABs| ABp| ABu-> ab
phase: Rs|Rp -> phi[0], Ts|Tp -> phi[2], RBs|RBp -> phi[3], TBs|TBp -> phi[1]
colour: Rs|Rp|Ru -> grad_r, Ts|Tp|Tu -> grad_t
```

`needle_pass_scan` then feeds **every** bucket through **every** enabled
polarization branch (`needle_pass.rs:1072-1224`), and production always
enables both (`cycle.rs:782`, `cycle.rs:1036`, `environments.rs:3726`).

The fold's "exact merge" (`Σ w_e (R − t_e) = W (R − t_eff)`,
`needle_pass.rs:22-28`) is exact only when the entries constrain **the same
variable**. Rs and Rp are different variables. Consequences:

| Target set | What the needle step optimizes | What the merit optimizes |
|---|---|---|
| Rs → 1, Rp → 0 at 45° (PBS) | Rs → 0.5 **and** Rp → 0.5 | Rs → 1, Rp → 0 |
| Rs only, at 45° | Rs → t **and** Rp → t | Rs → t |
| Ru → t at 45° | Σ_pol w (R_pol − t)² (per-branch pull) | w ((Rs+Rp)/2 − t)² |
| Rs phase and Rp phase | one blended phase target on both branches | two independent phases |
| colour on Rs | dF/dRs applied to the s **and** p branch | s only |

The thickness refinement (LM) uses the correct merit, so final numbers are
consistent. The damage is to **needle site choice and predicted gain**:
wrong insertion sites, wasted cycles, premature "no negative P" stops, and
convergence checks fed by a wrong P. At normal incidence (Rs = Rp) the
sites agree but P magnitudes are off by factors of 2 (Ru: 4w vs 2w).

The existing tests (`dual_pol_equals_sum_of_branches`,
`u_curve_half_is_sum_of_single_branch_halves`,
`grad_buckets_add_across_branches_bitwise`) check internal consistency
(both branches add), not agreement with the merit's derivative, so they
pass with the bug.

The same design is in the engine door: `solver.rs:1721` (`needle_gradient`,
`navette.smatrix.needle`) takes one target/weight pair per quantity and
runs it through `calc_s`/`calc_p`. The Python fold dict
(`synthesis_merit.rs:446-468`: keys `r, t, a, rb, tb, ab, grads_r, grads_t,
phi0..3`) carries the same layout.

### 1.2 The evaluator fills four intensity curves (b, c)

`SmatrixContext::simulate_inner` (`evaluator.rs:203-229`) stores Rs, Rp,
Ts, Tp and, when phase is used, the complex forward t_s, t_p. Nothing else.
Probe (2-layer stack, 45°, one target each, `SmatrixContext` +
`optimize_thicknesses`):

| Target | `evaluate_merit` | `optimize_thicknesses` |
|---|---|---|
| Rs, Ts, As, Ts-phase | correct | runs |
| Ru, Tu, Au | **1e6** (missing-curve penalty, silent) | Err "missing curve Ru" |
| RBs, TBp (back) | **1e6** | Err |
| Rs-phase (front reflection phase) | **1e6** | Err |

Notes:

- The dual solver **already computes** what is missing. `BlockResult` is
  `(rf, tb, tf, rb, R, Tb, T, Rb)` per polarization; `simulate_inner` keeps
  `.4` (R), `.6` (T) and `.2` (t_f), and drops `.0` (complex r_f), `.3`
  (r_b), `.1` (t_b), `.5` (T_b) and `.7` (R_b).
- The needle pass already has operators for every back quantity
  (`p_coherent_{rb,tb,ab}_from_fields`), so back targets fold fine and
  then fail in the merit.
- `SpectralInputs::from_spec` (`pipeline.rs:435`) does no capability
  check, so the failure surfaces late (first LM step) or silently (merit).
- `merit()` returns a 1e6 penalty per missing key group: a target set
  made only of Ru demands reports a flat merit that never moves.

### 1.3 Other target gaps (c)

- `normalization_mode = "log"` is accepted by the merit and LM, refused by
  the needle fold ("needle pass requires linear-normalized targets").
- Colour demands on back curves are refused at compile (documented v1
  limit); on Ru/Tu they compile and then hit 1.2.
- Solver channels without a `CurveId` (Psi/Delta, DOP, retardance, GD/GDD,
  Stokes) cannot be synthesis targets at all.

### 1.3b Full target-type matrix

Probe of 2026-11-09: every label × polarization × phase flag, plus kinds,
integral, angular, log and colour, through `build_merit_spec` →
`SmatrixContext.evaluate_merit` → `optimize_thicknesses` →
`build_needle_targets` (2-layer stack, 45°). "Penalty" = merit exactly
1e6, the silent missing-curve value.

**Synthesis vocabulary (what a target can name today)**

| Target | Compile | Merit | LM | Needle fold | Status |
|---|---|---|---|---|---|
| Rs, Rp | ok | ok | ok | ok* | **usable** (*W1 bug at oblique incidence) |
| Ts, Tp | ok | ok | ok | ok* | **usable** (*W1) |
| As, Ap (1 − R − T) | ok | ok | ok | ok* | **usable** (*W1) |
| Ts, Tp phase (absolute) | ok | ok | ok | ok* | **usable** (*W1: s and p phase share `phi[2]`) |
| PDts, PDtp (differential phase) | ok | ok | ok | ok | **usable** |
| kinds e/a/b/r/c, integral, angular targets, `normalize_count` | ok | ok | ok | ok* | **usable** (*W1) |
| colour on Rs, Rp, Ts, Tp | ok | ok | ok | ok* | **usable** (*W1: grad bucket shared by s and p) |
| Ru, Tu, Au | ok | **penalty** | Err | folds (wrong form) | **easy** (W2: derive from s/p) + W1 gradient form |
| RB, TB, AB × s/p/u (back side) | ok | **penalty** | Err | ok | **easy** (W2: dual solver already computes `.7`/`.5`; needle operators exist); refuse on absorbing substrate |
| Rs, Rp phase (front reflection) | ok | **penalty** | Err | ok | **easy** (W2: keep complex r_f `.0`; needle `phi[0]` exists) |
| RBs/RBp/TBs/TBp phase | ok | **penalty** | Err | ok | **easy** (W2: keep `.3`/`.1`; `phi[3]`/`phi[1]` exist) |
| colour on Ru, Tu | ok | **penalty** | Err | ok | **easy** after W2 |
| `normalization_mode="log"` | ok | ok | ok | **Err** | **easy** after W1 (gradient form, W5) |
| colour on back curves (RB/TB) | **refused** (v1) | — | — | — | **easy** after W2 (lift the compile refusal; needs the back rows) |
| colour on absorptance (A) | refused | — | — | — | out of scope (colour of absorbed light is not a measured colour) |

**Engine channels with no synthesis target yet** (`core_engine.rs` REQ_*
bits; needle gradients in `needle_engine.rs` NREQ_* bits)

| Channel (engine) | Forward data in synthesis | Needle derivative | Status |
|---|---|---|---|
| R_avg, T_avg, A_avg | = Ru/Tu/Au | via Ru | covered by Ru/Tu/Au above |
| complex amplitudes rs, rp, ts, tp (+ back) | W2 keeps them | `needle_slopes4_ddz` per channel | building block, not a target |
| φ_RS, φ_RP, φ_TS, φ_TP, φ_RB*, φ_TB* | = phase targets above | phi channels | covered above |
| PD on reflection (PDrs, PDrp; `passes = 2`) | needs r_f (W2) | phi + gain shift (`reference_phase` handles 2 passes) | **deferred** (decision 2026-11-09): reference-plane question for the round trip stays open, see §5 |
| Delta_R, Delta_T | from the retardance curve | as retardance | **in W8a**: constant offset of retardance (`Delta_R = Ret_R − π`, `Delta_T = Ret_T`, wrapped), two extra labels |
| Psi_R, Psi_T | from rp/rs (W2) | gradient form on complex amplitudes (new small operator `p_coherent_grad_c_from_fields`: Re(conj(g)·∂a)) | **moderate**: blocked by W1 + W2 + new operator |
| Rs/Rp ratio, extinction ratio, diattenuation | from intensities | gradient form on R/T (W5) | **moderate**: blocked by W1/W5 |
| retardance Ret_R, Ret_T (φp − φs) | Ret_T from the t_f rows kept today; Ret_R needs r_f (W2) | gradient form on two phases (+g on p, −g on s) | **in scope, W8a**: blocked by W1 (phi buckets are per S-matrix element and feed both branches today) and W2 (Ret_R) |
| DOP, Stokes S0-S3 (for a given input state) | from complex rows | gradient form, complex | **moderate**: needs an input-polarization field on the target; blocked by W1/W2/W5 |
| GD, GDD, TOD, FOD (R/T, s/p) | **not computed** in the synthesis evaluator (engine `dispmaps`) | **exists** in the engine (`phase_dispersion_sensitivity`, NREQ_DGD…DFOD) | **harder**: port the dispersion forward into `simulate_inner`, add curve ids, map the fold to Q rows; LM via FD first |
| cross_R, cross_T (cross-polarized) | — | — | **not meaningful** for isotropic stacks (identically zero); blocked by the absence of anisotropic media in synthesis |
| any channel through an incoherent layer / real substrate back face | — | multiblock operators exist (`p_function_multiblock`) for R/T/A/RB/TB/AB | **blocked by W6**; phase through a spacer is physically undefined and stays refused |
| angle-cone or bandwidth-averaged R/T | — | linear in the pointwise terms | **moderate**: needs point weights (W4) and a quadrature on the solver grid |

### 1.4 Weighting (a)

- Scalar `weight` per target: residual × √(weight / count).
- Per-point `tolerances`: residual = d / tol_i; the needle fold uses
  nf² / tol_i². So **tol_i = tol₀ / √w_i is a working per-point weight
  today**, for Exact/Above/Below pointwise targets, on spectral (over λ)
  and angular (over θ) targets alike.
- Where it breaks:
  - Range (`r`) uses tol as the fallback band half-width; CenterBand (`c`)
    uses it for the outside slope. Changing tol to weight a point changes
    the constraint.
  - Integral targets compute mean(d) / mean(tol): per-point weighting is
    averaged away.
  - Colour demands have a `wavelength_range` window but no profile.
  - No set-level profile (e.g. "weight 500-600 nm ×3 for every target").

### 1.5 Material deduplication (d)

`dedup_pool` (`cycle.rs:587`) keeps the first label of each bit-identical
index array. Throughout synthesis, **a material is its nk array**: insert
rules resolve hosts with `entry_of(nk)` (`cycle.rs:211`), adjacency and
boundary checks compare nk (`cycle.rs:363-373`, `426`, `461-470`), and
merges join adjacent films of equal nk (`structure.rs:1221`).

Rationale (sound): the needle operator sees only the index. Two entries of
one index cost two identical scans and tie in the argmin; a needle of a
neighbour's index is that neighbour grown thicker, which LM covers.

Risks found:

1. A dropped label vanishes: insert rules or edge lists naming it fail
   with "not a pool material; the pool holds [...]" although the user
   listed it.
2. Two labels of one index share one insert-rule row; per-label rules are
   impossible.
3. Adjacent films of equal nk merge; the second film's name is lost
   (optically right, bookkeeping loss).
4. Bitwise equality: a table and a model of the same material, or two
   resamplings, are "different" and both get swept.
5. Group `n_factor`/`k_factor` change a layer's nk, so a scaled H layer is
   a different material from the H pool entry (a needle of H next to it is
   allowed), and `pool_from_layers` can emit two entries with the same
   label, which `label_index` then refuses as ambiguous.
6. Two labels with one raw index can carry **different groups** and
   **different insertion rules**. Rules are looked up by the host's nk
   (`entry_of`), so the second label silently gets the first label's
   rules; needles join no group; and the default pool (`pool_from_layers`,
   nk after group factors) and the explicit pool (`build_pool`, library nk
   before them) sweep different indices for the same material. W7 is
   therefore a redesign (identity by label), not hygiene.

### 1.6 Coherence (e)

`DesignStack` construction refuses `coherent = false`
(`structure.rs:334`), because `simulate_inner` solves `[0, nl−1)` as one
coherent block and the flag used to be ignored silently. The substrate is
always a half-space, so a real plate's back surface cannot be part of a
synthesis. The engine already has:

- `partition_blocks`, the incoherent cascade (`star_real`, `cascade_step`,
  `block_intensities`, `block_flux_factors`);
- `p_function_multiblock` / `p_multiblock_point`
  (`needle_operator.rs:1086`, `1491`): needle gradient through incoherent
  spacers, adjoint over the real-star cascade, FD-validated; per call one
  `PmbQuantity` (R, T, A, TB, RB, AB);
- the engine's general solver for multi-block forward simulation.

Synthesis calls none of them.

---

## 2. Workstreams

### W1 — Polarization-correct needle fold (bug, b)

**Goal:** the needle P(z) equals dM/dδ of the merit the LM minimises, for
every demand type, at any angle.

**Done.** Where the implementation differs from the design below:

- The intensity kernels return w·(X − t)·∂X, half the slope; the phase
  kernel returns the full slope. So an s or p intensity pair folds at
  2w on its branch, and a phase pair at w.
- An unpolarized demand stays a pair, not the gradient form: on each
  branch at weight w, with target `X_b − (X_u − rt)`, so that
  `X_b − t_b = X_u − rt` and both branches carry the shared u residual.
  Without an operating point it folds `(rt, w)` on both, which is exact
  where Rs = Rp. `required_rows` asks for the s/p rows of every
  unpolarized demand, so the operating point is there whenever the fold
  has a simulation.
- Colour on an s or p curve deposits 2g on its branch; on a u curve, g on
  both branches.
- `phi_gain_shift` lives in each `PolFold`; `NeedleTargets::phi_gain_shift()`
  sums it. Two older bugs surfaced: a differential-phase demand handed the
  kernel its target without the reference phase, and its gain shift used
  the unwrapped residual. Both are fixed.
- The engine door added the shift with the wrong sign and once per point;
  it now adds `shift / n_points` to each point.
- Oracle: `evaluator.rs` `needle_oracle` checks P + shift against a
  Richardson finite difference of the merit after a real needle insertion,
  for 26 demand sets at 0°, 30° and 60° (≤ 1e-5 relative).

**Fold layout.** Buckets per (quantity, polarization branch):

```rust
pub struct PolFold {            // one per branch: s, p
    pub r: Pair, pub t: Pair, pub a: Pair,
    pub rb: Pair, pub tb: Pair, pub ab: Pair,
    pub phi: [Pair; 4],
    pub grad_r: Vec<f64>, pub grad_t: Vec<f64>,
    pub grad_rb: Vec<f64>, pub grad_tb: Vec<f64>,   // for W2/W5 forms
}
pub struct NeedleTargets {
    pub s: PolFold,
    pub p: PolFold,
    pub phi_gain_shift: [[f64; 4]; 2],
}
```

The scan runs branch s with `fold.s` only and branch p with `fold.p` only.
A branch whose buckets are all zero is skipped (an s-only target set costs
one field build per point instead of two: a speed-up).

**Routing.**

| Demand | Form | Branch |
|---|---|---|
| Rs, Ts, As, RBs, TBs, ABs | (target, weight) pair | s |
| Rp, Tp, … | pair | p |
| phase on an s/p curve | phi pair on its channel | its own branch |
| Ru, Tu, Au, RBu, TBu, ABu | **gradient form** (below) | s and p |
| colour on s/p curve | grad bucket | its own branch |
| colour on Ru/Tu | grad bucket, g/2 | s and p |

**Gradient form for unpolarized demands.** The merit term is
w·(½(Rs+Rp) − t)², so dM/dδ = w·(Ru − t)·(∂Rs + ∂Rp). This is not a
(target, weight) pair per branch: the residual couples the branches. It is
exactly the colour mechanism: deposit the chain factor
g = w·(Ru_op − t) (incl. kind activation, tolerance, count) into
`grad_r` of **both** branches; `p_coherent_grad_r_from_fields` multiplies
by ∂R_pol. Exact at the operating point, which the fold always has when
`refold_per_cycle` is on (the default).

**First fold without a simulation** (`from_spec` passes `None`): fall back
to the per-branch pair form with weight w/2 per branch, which equals the
exact form when Rs = Rp (normal incidence) and is the documented
conservative fold elsewhere. Better: `from_spec` gets the initial
simulation (one solve, cheap) and the static form disappears. Proposal:
do the latter; keep the pair fallback only for mock contexts.

**Engine door and Python.** The same change in `solver.rs`'s
`needle_gradient` and the Python fold dict (`s`/`p` sub-dicts). Hard
change, no compatibility (project policy); update `navette.smatrix.needle`,
`.pyi`, `test_needlefold_mirror.py`, `test_weighting.py`,
`test_differential_phase.py`, `bench_refold.py`.

**Tests.**

- New oracle test per demand type: P(z) from the fold vs a central FD of
  the **merit** after inserting a needle of thickness ±δ at z, at 0°, 30°
  and 60°, for Rs, Rp, Ru, Ts, Tp, Tu, As, Au, RBs, TBp, ABu, Rs-phase,
  Tp-phase, PDts, colour on Rs and Ru, every kind (e/a/b/r/c) and
  integral. Tolerance ~1e-6 relative.
- PBS case: Rs → 1, Rp → 0 at 45°: fold P matches merit FD (fails today).
- Replace the three consistency tests with the oracle (they stay as
  internal checks, reworded).

**Effect on pins.** Needle runs will change: sites at oblique angles, P
magnitudes everywhere (normal-incidence Ru halves; s-only sets lose the p
branch). The bench fingerprint `46dc154d6946ee40` and the needle digests
**will move** and must be re-baselined in a dedicated commit, with a
before/after table of final merits on the bench problems as the evidence
that the change is an improvement (or at worst neutral at 0°).

### W2 — Evaluator curve completeness (bug, b/c)

**Goal:** every `CurveId` the target compiler accepts is produced by
`simulate_inner`, or refused up front (W3).

1. **Unpolarized**: Ru = (Rs+Rp)/2, Tu = (Ts+Tp)/2 into `curves`; Au/As/Ap
   are already derived from companions. Free.
2. **Back side**: R_b = `.7`, T_b = `.5` into `back`. Caveat already noted
   in `evaluator.rs:131-134`: with an **absorbing substrate** as the
   incidence medium, back incidence from a semi-infinite absorber is
   ill-defined and `.5` divides by Re(y_last) ≈ 0. Refuse back demands
   when the substrate has k > threshold at any demand wavelength (named
   error), rather than produce T_b = 84.
3. **Complex amplitudes**: r_f (`.0`) for front-reflection phase; r_b
   (`.3`), t_b (`.1`) for back phase, into `cplx` / `cplx_back`, gated by
   `uses_phase()` as today.
4. **Deposits (analytic LM Jacobian)**: `deposit_row` covers ∂(Rs, Rp, Ts,
   Tp) and ∂t_f. Extend:
   - Ru/Tu: linear combination of existing rows (free);
   - R_b, T_b, r_f, r_b, t_b: `thickness_slopes_ddz` returns [dr, dt] for
     the forward channels; add the back channels (same field
     decomposition; the slopes for s_left/s_right back elements), with the
     FD-vs-analytic test the forward deposits have;
   - until a channel has a deposit, the LM falls back to FD Jacobian for
     specs that use it (existing `JacobianMode::Fd` path), stated in the
     report.
5. **Cost**: intensity-only specs stay bit-identical (new rows only when a
   demand needs them; gate per curve with `spec.uses_curve(id)`).

Tests: probe table of 1.2 becomes a test (every CurveId × phase flag:
merit finite, LM runs, merit vs `ScatterMatrix` oracle at 1e-12).

### W3 — Capability check and loud failures (bug, c)

1. `DesignContext::provides() -> CurveSet` and
   `MeritSpec::required_curves() -> CurveSet` (incl. complex/back/phase
   needs). `SmatrixContext` construction, `SpectralInputs::from_spec`,
   `run_design`, `run_environments` and the Python doors refuse a spec
   whose requirements exceed what the context provides, naming the demand
   (`spectral[3]: RBs at 45°`) and the reason.
2. Missing-curve **penalty** in `merit()`: keep it only for genuine runtime
   failures (a solve error mid-run), never for a curve the context never
   produces. After W3 that case cannot reach `merit()`.
3. Needle-pass refusals (log transform, unsupported colour) move to the
   same up-front check, so a run never starts that cannot finish.

### W4 — Per-point weight profiles (feature, a)

**Schema.** Optional `weights: Vec<f64>` on `SpectralTarget` and
`AngularTarget` (same length as values, finite, ≥ 0; default all ones),
mirrored in the Python dataclasses and `.pyi`. `MeritTarget` gains
`point_weights: Arc<[f64]>` (empty = all ones, so legacy specs stay
bit-identical).

**Semantics.** Point weight multiplies the scalar weight:
residual_i × √(weight · w_i / count). Tolerances keep their meaning
(scale and band); weights carry emphasis. One rule for all kinds.

| Place | Change |
|---|---|
| `residuals_into` (`merit.rs:1229-1298`) | `rscale_i = rscale · √w_i` pointwise |
| integral targets | weighted mean: R = Σw_i d_i / Σw_i over Σw_i tol_i / Σw_i |
| derivative passes (`merit.rs:1698-1717`, `1904-1917`) | same factors |
| needle fold (pointwise and `fold_integral_demand`) | `w *= interp(w_i)` at emission |
| TargetWeaver ingestion / `register_metadata` | carry the array; angular targets split it per angle like tolerances |
| colour demands | **later** (W4b): a spectral weight profile inside the colour integral changes what the colour *is*; out of scope, keep `wavelength_range` |

**Profiles live in Rust.** All profile logic is a Rust module
(`synthesis::weights`); Python only wraps it, like every other door. A
target carries either explicit `weights` or a `weight_profile`, never both:

```rust
#[derive(Deserialize, Serialize)]
#[serde(tag = "kind", deny_unknown_fields)]
pub enum WeightProfile {
    Constant { value: f64 },
    Ramp { lo: f64, hi: f64, w_lo: f64, w_hi: f64 },          // linear, clamped outside
    Gaussian { center: f64, fwhm: f64, peak: f64, base: f64 },
    Bands { bands: Vec<(f64, f64, f64)>, outside: f64 },        // (lo, hi, w), overlaps refused
    Table { x: Vec<f64>, w: Vec<f64>, interp: Interp },         // pchip | linear (crate::interpolate)
    Product(Vec<WeightProfile>),                                // emphasis × window
}
impl WeightProfile {
    pub fn evaluate(&self, x: &[f64]) -> Result<Vec<f64>, String>; // x = λ (nm) or θ (rad)
}
```

- The axis is implied by the target: wavelength for spectral targets,
  angle for angular targets (angles in radians internally; the JSON door
  accepts `*_deg` / `*_rad` keys through the existing angle-key resolver).
- **Named set-level profiles**: `TargetSet.weight_profiles: {name →
  WeightProfile}`, referenced as `"weight_profile": "visible_emphasis"`,
  resolved at compile into point weights. Saves repeating arrays across
  many targets and keeps one definition per emphasis.
- Validation in Rust (finite, ≥ 0, not all zero, monotone `Table.x`,
  `fwhm > 0`, `lo < hi`), with messages naming the target and profile.
- Python: `WeightProfile` pyclass wrapping the enum (constructors
  `WeightProfile.ramp(...)`, `.gaussian(...)`, …, `.evaluate(x)`), and the
  `weight_profile=` keyword on the target dataclasses that just forwards
  the JSON. No profile maths in Python.

**Tests.** Weights all ones = bit-identical to today; weight w_i ≡ tol
rescaling equivalence for Exact pointwise; integral weighted mean vs a
hand computation; fold vs merit FD with non-uniform weights (W1 oracle).

**Docs.** Add a section to `docs/spectralweave-target-kinds.md`:
tolerance (scale/band) vs weight (emphasis), and the old
tol₀/√w trick as superseded.

### W5 — Transform coverage in the needle pass (feature, c)

The gradient form of W1 generalises: any residual r(R) folds as
g = w·r(R_op)·r'(R_op) into the branch's grad bucket. That covers:

- `SimTransform::Log` (r' = 1/(R ln 10)·nf),
- unpolarized demands (W1),
- colour (already),
- later: ratio targets (Rs/Rp), extinction ratio, polarization contrast.

Activation (Above/Below/Range/CenterBand) uses the same op-point logic as
today. Without an op point the demand folds nothing for non-linear forms
(documented), which W1's "initial simulation in `from_spec`" makes moot.

### W6 — Incoherent synthesis (feature, e)

Phased, larger:

1. **Model**: allow `coherent = false` films, and an **exit medium**
   behind the substrate (`DesignStack { ambient, films, substrate,
   exit: Option<LayerSpec> }`, the substrate then a finite incoherent layer
   of given thickness). This is the real use case: a coated plate with its
   back surface (optionally coated too).
2. **Forward**: `simulate_inner` partitions with `partition_blocks`; solves
   each coherent block with the dual solver; cascades intensities with the
   real-star product; phases only for demands inside a single block
   (front reflection of the first block; refuse phase demands that cross a
   spacer, by name).
3. **Needle**: map each fold bucket to `PmbQuantity` and call
   `p_function_multiblock` per (quantity, branch); host mask excludes
   spacers; sites only inside coherent blocks.
4. **LM**: FD Jacobian first; analytic deposits through the cascade
   (adjoint weights W[b][c] already derived for the needle) as a follow-up.
5. **Environments (F2)**: a backside segment is an environment segment;
   check the joint door's routing handles an incoherent substrate.
6. **Pins**: fully coherent stacks must stay bit-identical (the single
   block path is unchanged when no flag is set).

Tests: forward vs `solve_structure`/`ScatterMatrix` on stacks with an
incoherent substrate and a coated back side; needle P vs merit FD
(oracle of W1) through one spacer; LM recovery.

### W7 — Material identity: label, not index (redesign, d)

**Problem with "a material is its nk array".** Two materials with
bit-identical raw indices and different labels can legitimately differ in
everything except the index:

- **(a) group factors**: `H` in group `G1` (n_factor 1.00) and `H2` in
  group `G2` (n_factor 1.03, thick_factor 0.97). Their *raw* index is
  equal, their *deposited* index and thickness errors are not.
- **(b) insertion rules**: host `H` accepts needles {L}, host `H2` accepts
  {M}. Today the host row is found by `entry_of(host.nk)`, so `H2` layers
  silently get `H`'s row (or, if `H2` was dropped by dedup, the rule
  naming it fails as "not a pool material").
- **(c) needle membership**: a needle of `H2` should join `H2`'s group and
  carry its factors. Today `new_needle` copies the pool entry's nk and
  joins no group. Worse, the two pool doors disagree: `pool_from_layers`
  takes the layer's nk **after** group factors, `build_pool` takes the
  library nk **before** them, so the same material can be swept with two
  different indices depending on how the pool was given.

So deduplicating **identity** by nk is wrong in general. Deduplicating
**computation** by nk is still right: two entries with the same effective
index produce the same P(z) profile.

**Proposal: identity by label, computation shared by index.**

1. **Pool entry identity** = pool label (library code), plus its group
   (if any). `MaterialPool` entries are never dropped. `dedup_pool` is
   replaced by a scan-sharing map: entries whose **effective** needle nk
   (after their group's n/k factors) is bit-identical share one P(z)
   sweep; the sweep result is reused, and each entry applies its **own**
   admissible-site mask (its own insertion rules, edge lists, adjacency).
   Same cost as today when rules agree; correct results when they differ.
2. **Films carry their pool label** explicitly (a `pool_label:
   Option<Arc<str>>` on `LayerSpec`, set for design layers from their
   material code and for needles from the entry), instead of being
   recognised by nk or by parsing `"H#3"`. Insertion rules resolve the
   host by `pool_label`; nk matching remains only as the fallback for
   films without a label (the `"*"` rule).
3. **Groups for needles**: a needle of entry `m` joins `m`'s group (from
   the design's group table; ambiguous if the label is in two groups →
   named error) and gets the group's n/k factors applied to its index.
   Thickness factors/summands of the group apply as for any member
   (decide: needles are *design* thicknesses, so thick_factor applies only
   in error draws, never to the optimized value — to confirm in review).
4. **One pool convention**: both pool doors evaluate the effective index
   the same way (library nk → group factors), so explicit and default
   pools agree.
5. **Adjacency / "needle of its own index"**: a needle is excluded next to
   a film of the **same label and group** (it would be that film grown
   thicker). Same nk but different label/group is allowed (it is a
   different deposited material). Bit-identical nk with no label: current
   behaviour.
6. **Merges** of adjacent films: only same label **and** same group **and**
   same nk merge. Today's "same nk at a layer boundary" merge would fuse
   `H` and `H2` films and lose a group membership: that rule becomes
   label-aware. The merge report names absorbed films.
7. **Near-equality** (optional, `pool_share_tol`, default 0 = bitwise):
   only for sharing the **scan**, never for identity.
8. **Report**: the run report lists shared sweeps ("H and H2 share one
   index scan; rules applied separately") and every needle's label and
   group.

**Risks.** Changes needle naming/merge behaviour for designs that relied
on nk-merging of differently named films (rare; named in CHANGELOG);
`LayerSpec` gains a field (Python `LayerSpec` signature, `.pyi`, state
round-trip); group-factor application to needles changes results for
group-scaled designs (intended: today they are inconsistent between the
two pool doors). Pins: unchanged for designs without groups or duplicate
indices; the bench problems have neither (verify in review).

**Tests.** Two labels, one raw index, different groups → two distinct
needles with different effective indices; two labels, one index,
different insertion rules → each host gets its own rule; explicit vs
default pool give identical runs; shared-scan cost equals one sweep
(timing counter); merge of `H`|`H2` refused, `H`|`H` merged.

---

## 3. Order and dependencies

| Step | Content | Depends on | Pins |
|---|---|---|---|
| 1 | W3 capability check (refuse what cannot be scored) | — | unchanged |
| 2 | W2 unpolarized + back + complex curves, deposits for Ru/Tu | W3 (to lift refusals as curves land) | unchanged (new rows gated) |
| 3 | W1 per-polarization fold, gradient form, initial sim in `from_spec`, engine door, Python dict | W2 (needs Ru/back rows at the op point) | **re-baseline** fingerprint and needle digests, with evidence |
| 4 | W2 back-channel deposits (analytic LM) | 2 | unchanged |
| 5 | W4 point weights (+ helpers, docs) | 3 (fold emission point) | unchanged at all-ones |
| 6 | W5 log and general transforms in the fold | 3 | unchanged |
| 7 | W7 label identity, shared scans, needle groups, label-aware merges | 3 (both touch the scan loop in `cycle.rs`) | unchanged without groups/duplicate indices |
| 4a | W8a retardance Ret_R/Ret_T (+ Delta_R/Delta_T labels) | 2, 3 (not W5: own gradient-form phase operator) | unchanged |
| 8 | W8 further target types, in this order: Rs/Rp ratio, extinction, diattenuation → Psi (+ complex gradient operator) → DOP/Stokes → GD…FOD. PDrs/PDrp deferred (§5) | 2, 3, 6 | unchanged |
| 9 | W6 incoherent synthesis | 2, 3 | unchanged for coherent stacks |

Steps 1-3 are the bug fixes and should land before any new feature on the
needle path. W7 follows W1 because both change the per-entry scan loop.

### W8 — New target types (feature, c)

From the matrix in §1.3b. Each new type is: a `CurveId` (or a derived
curve id with its formula), a row in W3's capability set, a forward value
from the complex rows W2 keeps, and a gradient-form fold (W5). Per type
the W1 oracle test (fold P vs merit FD) is the acceptance gate. GD…FOD
additionally need the dispersion forward in `simulate_inner` (port of the
engine's `dispmaps`) and a fold onto the engine's existing Q-row needle
sensitivities.

### W8a — Retardance targets (feature, c; decision 2026-11-09)

Neither retardance nor ellipsometric Δ is a synthesis target today (no
`CurveId`, no fold). The engine defines both (`solver.rs`, `core_engine.rs`):

- `Ret_R = arg(cross_r)`, `cross_r = r_p·conj(r_s)` → φ_rp − φ_rs, wrapped
  to (−π, π].
- `Ret_T = arg(cross_t)`, `cross_t = t_p·conj(t_s)·√(f_s f_p)`; the flux
  factor is real positive, so `Ret_T = φ_tp − φ_ts`.
- `Delta_R = atan2(−2 Im cross_r, −2 Re cross_r) = Ret_R − π` (wrapped);
  `Delta_T = Ret_T` (docs/forward-golden.md). Synthesis follows the
  engine's numbers exactly; the oracle test compares against them.

**Design.**

1. **Curves.** Derived curve ids `RetR`, `RetT` (labels `RetR`, `RetT`,
   `DeltaR`, `DeltaT`; polarization field not used, refused if given as
   anything but `u`/none). Delta labels compile to the same curve with a
   constant offset (−π for R, 0 for T) applied to the target values before
   wrapping, so one forward and one fold serve all four labels.
2. **Forward.** From the complex rows: `arg(p·conj(s))`, not
   `arg(p) − arg(s)`, so the wrap happens once. Ret_T needs only the t_f
   rows kept today; Ret_R needs W2.3's r_f. Degenerate points follow the
   engine exactly: `Retardance_R/T` have no floor (arg of a zero cross
   term is 0), `Delta_R` is 0 where R_s < `RS_FLOOR` (1e-12). The fold
   contributes zero gradient wherever |a|² < 1e-20 on either branch (the
   existing phase-operator rule).
3. **Merit.** Phase-style residual: `wrap(Ret − t)/tol`, same kinds and
   `integral` handling as phase demands today. `norm_factor` must be 1
   (raw radians), same rule as phase.
4. **LM.** FD Jacobian first; an analytic deposit (`Im(conj(a)·∂a)/|a|²`
   on p minus the same on s, both already in the phase deposit path)
   follows behind an FD-vs-analytic test.
5. **Needle fold.** Gradient form: `g = 2·w·wrap(Ret_op − t)` deposited
   as **+g** into the p-branch phase-gradient bucket of the element and
   **−g** into the s-branch one. New buckets `phi_grad: [[f64; 4]; 2]`
   (branch × element) next to W1's per-branch `phi` pairs, and a new
   operator `p_coherent_grad_phi_from_fields(…, channel, grad, …)` =
   `grad · Im(conj(a)·∂a)/|a|²`. The existing
   `p_coherent_phi_from_fields` becomes a thin call of it with
   `grad = 2·w·wrap(φ − t)` (bitwise-identical, pinned by the existing
   phase tests).
6. **Capability (W3).** `RetR`/`DeltaR` require complex r_f s and p;
   `RetT`/`DeltaT` require complex t_f s and p. Refuse on incoherent
   stacks (phase through a spacer, §1.6).
7. **Python.** `SpectralTarget(…, spectrum="RetR"|"RetT"|"DeltaR"|"DeltaT")`
   and the same on `AngularTarget`; the wrapper only forwards labels.

**Why W1 first.** The phase buckets are indexed by S-matrix element
(0 = r_front, 1 = t_back, 2 = t_fwd, 3 = r_back) and both polarization
branches consume them, so +g on p and −g on s would cancel in a shared
bucket. Retardance is the cleanest demonstration that W1 is needed.

**Tests.** Forward equals the engine's `Retardance_R/T` and `Delta_R/T`
bitwise on a random stack at 0°, 45°, 70°; fold P vs merit FD (W1
oracle) for RetR, RetT, DeltaR at 45° and 70°; at 0° (s ≡ p) the forward
is identically 0 (R) and the fold P is 0; wrap crossing (target near ±π)
gives the short-way residual; degenerate point gives zero contribution;
a retardance-only synthesis on a quarter-wave-retarder benchmark
(e.g. Ret_R = π/2 at 45°, 500-600 nm) converges.

**Risks.** At normal incidence retardance is identically zero and the
target is unreachable; warn (not refuse) when every angle of a retardance
demand is 0°. Near |r_s| → 0 (Brewster-adjacent for p, or a strong
anti-reflection for s) the phase gradient blows up as 1/|a|²; the floor
plus a per-point tolerance (W4) is the user-side mitigation.

---

## 4. Risks

| Risk | Where | Mitigation |
|---|---|---|
| Fingerprint and needle digests move (W1) | every needle bench and golden | one dedicated re-baseline commit; before/after table of final merits and film counts per bench problem; refusal to re-baseline anything else in that commit |
| Behaviour change for existing users at 0° (P magnitudes halve for Ru/s-only sets) | stagnation thresholds, convergence on summed P (`cycle.rs` gain shift, `stagnation.rs`) | check every P-scaled threshold; most compare relative gains; document in CHANGELOG "Changed (breaking)" |
| Gradient form is exact only at the op point | Ru, log, colour | already the colour contract; refold per cycle is default; initial sim in `from_spec` removes the static case |
| Back incidence from an absorbing substrate | W2 | refuse with a named error (k threshold) rather than return non-physical T_b |
| Field builds double for s+p sets vs today's shared buckets | W1 | unchanged count when both branches are used (today already builds both); s-only sets get faster |
| Fold dict layout change breaks Python callers of `needle_gradient` | `navette.smatrix.needle`, parity tests | hard change per project policy; update all in-repo callers and `.pyi` in the same commit |
| Weighted integral changes meaning of `integral` with weights | W4 | all-ones default is identical; documented formula |
| Aliased labels hide a user's intent to sweep two "different" materials | W7 | the run report names every alias; near-equality off by default |
| Incoherent cascade numerics (near-zero Re y, regularisation in `star_real`) | W6 | reuse the crate-shared regularisation; oracle tests vs the engine |
| Analytic deposits for back channels disagree with FD | W2.4 | FD fallback stays; FD-vs-analytic test per channel before enabling |

---

## 5. Further and future gaps (not in W1-W7)

1. **Ellipsometric and polarization targets in synthesis**: Psi,
   Rs/Rp ratio, extinction ratio, DOP (retardance and Delta moved into
   W8a). Need `CurveId`s, complex rows (W2.3 gives r_s, r_p), and
   gradient-form folds (W5). Natural next step after W5.
1b. **Differential phase on reflection (PDrs, PDrp) — deferred
   (decision 2026-11-09).** `reference_phase` already takes `passes = 2`;
   the open question is the reference plane of the round trip
   (implementation_plan_pd.md:721). Revisit once a user case fixes the
   convention; the fold is the existing phi pair plus the gain shift.
2. **Dispersion targets**: GD/GDD/TOD on reflection or transmission
   (chirped mirrors). The engine computes them (`dispmaps`); synthesis
   needs curve ids, a merit form and the needle derivative of GD (phase
   derivative in λ). Larger research item.
3. **Angle- and bandwidth-averaged targets**: cone of incidence, finite
   spectral bandwidth. Integration weights over (θ, λ) on the solver grid;
   W4's point weights are the building block.
4. **Colour on back incidence** (documented v1 limit) and colour with a
   spectral emphasis profile (W4b).
5. **Weight profiles per environment** (F2 multi-environment): a profile
   that differs per environment for the same demand.
6. **Merit vs fold parity as a standing property test**: a randomized
   target-set generator checking P = dM/dδ on random stacks (would have
   caught 1.1).
7. **General fit engine interplay** (`general_fit_plan.md`): its
   `MeritSpec` adapter (§3.5.3) needs W2's curve coverage; its parameter
   layer will make group `n_factor` a fit parameter, which W7.4's labelling
   must survive (index changes every evaluation, label must not).
8. **Error locations for notebook users**: W3's Rust messages name the
   target by section and index (`spectral[3]`); the Python wrapper only
   forwards them. No validation logic in Python (Rust-first rule).

---

## 6. Follow-ups after this plan

- Update `docs/spectralweave-target-kinds.md` (fold reference: per-branch
  buckets, gradient form, weights vs tolerances).
- Update `docs/needle-*` docs and the synthesis sections of the user docs
  with the supported-target matrix (generated from W3's capability sets so
  it cannot drift).
- CHANGELOG: W1 under "Fixed" + "Changed (breaking)" (fold dict, P
  magnitudes, re-baselined pins); W2/W3 under "Fixed"; W4-W7 under
  "Added".
- Memory/release ritual: the two fingerprints change once in step 3;
  record the new values.

## 7. Open questions

1. W1 first fold: compute an initial simulation in `from_spec` (proposed),
   or keep a static conservative fold?
2. W3: should `merit()` keep the 1e6 penalty at all after W3, or become a
   hard error everywhere?
3. W4: ship named set-level profiles in the first version (proposed), or
   per-target `weights` / `weight_profile` only?
4. W6: is a coated back side (two coated faces) needed, or is "uncoated
   back surface" enough for the first version?
5. W7.3: do a group's thickness factor/summand apply to a needle's
   optimized thickness (proposed: no, only in error draws), or to every
   member including needles?
6. W7.6: is fusing adjacent films of different labels but identical nk
   ever wanted (e.g. a user splitting one layer into two named halves)?
   Proposed: never fuse across labels.
