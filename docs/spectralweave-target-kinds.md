# SpectralWeave target kinds and the needle fold

Target semantics live in two Rust mirrors that must stay in lockstep:

- `rust/navette/src/spectralweave/targetweaver.rs` — ingestion
  (`TargetKind`, `TargetEntry`, `register_metadata`) and
- `rust/navette-py/src/spectralweave_target.rs` — `calculate_merit`.
- `rust/navette/src/smatrix/synthesis/merit.rs` — `ConstraintKind`,
  `MeritTarget`, `MeritSpec`: a verbatim lift of the merit kernel into
  residual space for the thickness optimizer and the needle pass.
- `rust/navette/src/smatrix/synthesis/needle_pass.rs` —
  `build_needle_targets`: folds a `MeritSpec` into flat
  `(targets_r, weights_r)` quadratics for the analytic needle operator.

Python surface: `src/navette/spectralweave/target.py`
(`SpectralTarget`, `AngularTarget`, `TargetCollection`).

## Normalization strategy (verified)

Each target curve is self-normalized at ingestion so heterogeneous
quantities (R ~ 1, weak T ~ 0.01, phases in radians) contribute comparably:
values are scaled to O(1), the sim is scaled identically, and the residual
is divided by the **raw** tolerance. Consequences:

| mode | norm_factor | scaled values | tolerance means |
|---|---|---|---|
| `linear` | $1/\|mean\|$ | raw·nf | **fraction of the curve mean** |
| `log` | $1/mean(\|log_{10}\|)$ | log₁₀·nf | fraction of the mean log magnitude |
| `phase` | $1$ | raw (radians) | absolute radians |
| `complex` | $1$ | raw (no normalization) | raw units |
| `auto` | log iff all-positive and max/min ≥ 100×, else linear | — | — |

So `tolerances=[0.05]` on an R ~ 0.5 curve is an effective absolute window
of 0.025 (5% of the mean) — tolerances are **relative to the curve level**,
not absolute, except in phase/complex modes. The same triple
(normalized, nf, floored tolerances) feeds `calculate_merit`, `MeritSpec`
and the fold, so all three agree by construction.

Resolution scope and guards:

- Normalization is **per curve**: one `(mode, nf)` for a whole spectral
  curve, and — since the angular fix — one shared pair for a whole
  angular curve (previously each angle normalized itself, weighting equal
  absolute deviations differently per point).
- **Zero-mean guard** (linear): mean-normalization would explode on
  zero-mean/cancelling data (old nf hit 1e12), so when
  $|mean| \le 10^{-9}\cdot spread$ the curve scales by the half-range
  instead (exact O(1) for symmetric data); all-constant data keeps raw
  scale (nf = 1). Non-degenerate curves are bit-for-bit legacy.
- **~= 1 guard** (log): all-$|log| \approx 0$ data (values ~= 1) falls
  back to raw log scale for the same reason.
- Tolerance floor (`tolerance_floor`, default 1e-12) clamps near-zero
  tolerances so the merit can never divide by zero.
- `complex` is currently **raw-scale linear** (real f64 data, nf = 1) —
  an opt-out of normalization, not complex-number support.

## Weights, count normalization, integral targets

Merit is a sum over points — a 200-point target outweighs a 10-point one
20:1 at equal residuals (measured), and a single angular point drowns next
to dense spectral curves. Two per-target knobs control this (all engines —
weaver, spec/LM, fold/needle — apply them identically):

- `weight` (default 1): multiplies the frame's merit sum (residuals scale
  by √weight, so LM Jacobians stay consistent). Relative importance across
  targets; finite and ≥ 0 (rejected otherwise, at every trust boundary).
- `normalize_count` (default off): divides the frame's sum by the
  target-level point count (spectral grid size, angular angle count —
  resolved at ingestion, since angular targets span many single-point
  entries). Turns the sum into a mean: equal say regardless of sampling.
  `tol·√N` per frame is exactly equivalent (pinned by test).
- `integral` (default off): constrains the MEAN of the scaled diffs —
  single residual `R = mean(d)/mean(tol)` with kinds applied once to the
  mean (integral-`a` = lower bound on the average). Rejects
  `normalize_count` (the mean already is one — the combo would
double-dilute). Regular and integral targets mix freely in one run.

Frame contribution = `weight × (Σr² / count)` pointwise,
`weight × kind(R)` integral. Missing-data penalties are unaffected by
weights (drop a target to silence those — weight 0 only mutes present
data). The CenterBand `+1` accounting becomes
$M_{true} = M_{folded} + \Sigma\, weight/count$ over violated `c` points.

The integral fold matches the mean-form merit's UNIFORM gradient exactly
at the operating point (`w_i = W/N²`, `t_i = s_i − N·G` per point) —
values differ by dropped constants (same loss class as the `+1`/overlap
terms), gradients superpose exactly across overlapping integral frames,
and the PD gain-shift formula needs no changes (it lands on
$−2W(\bar m−T)\overline{k_z}$ automatically).

## Kinds (per point, merit space)

Scaled residual `d = sim_scaled − target_scaled`, floored tolerance `tol`,
scaled band half-width `bw`. `d` is phase-wrapped to $[-\pi,\pi]$ in phase
mode; in log mode both `d` and `bw` live in log-normalized space.

| kind | code | merit contribution |
|---|---|---|
| Exact | `e` | $(d/tol)^2$ |
| Above (lower bound) | `a` | $d<0 \Rightarrow (d/tol)^2$, else $0$ |
| Below (upper bound) | `b` | $d>0 \Rightarrow (d/tol)^2$, else $0$ |
| Range (hard box) | `r` | $\|d\|\le bw \Rightarrow 0$, else $((\|d\|-bw)/tol)^2$ |
| CenterBand (soft box) | `c` | $\|d\|\le bw \Rightarrow (d/bw)^2$, else $((\|d\|-bw)/tol)^2 + 1$ |

`r` is exactly the combination of paired `a`/`b` targets at centre∓band.
`c` keeps an `e`-style centre with proportionally reduced weight inside:
$(d/bw)^2 = (d/tol)^2\cdot(tol/bw)^2$. A band of $1\%$ with tolerance
$0.1\%$ pulls 100× softer inside than `e`, reaches $1.0$ at the band edge
from both sides (continuous by construction), then grows with the outer
$0.1\%$ scale outside.

## The `band` parameter

Per-point half-width in **raw units** (same units as `values`); a Python
scalar broadcasts. Scaled at ingestion by the same `norm_factor` as the
targets — exact for linear/phase/complex, per-point exact on the upward
side for log (symmetric approximation otherwise; avoid huge relative
bands on log targets). Negative bands are rejected in Python and floored
in Rust. Omitting `band`:

- `r` falls back to `tol` as the half-width (dead-band of ±tolerance);
- `c` degrades gracefully to `e`.

`MeritSpec` accepts an empty `band` array meaning all-zero (unused).

## The needle fold

The analytic needle operator understands one merit shape per spectral
point — a homogeneous quadratic $f_k = w_k\cdot(R_k - t_k)^2$ — fed via
`targets_r`/`weights_r`. `build_needle_targets` converts each `MeritSpec`
entry at the current operating point (active-set linearization,
recomputed every iteration):

| kind | condition | folded $(t_k, w_k)$ |
|---|---|---|
| `e` | always | centre, $nf^2/tol^2$ |
| `a` / `b` | violated (at current sim) | centre, $nf^2/tol^2$; satisfied → skipped |
| `r` | violated | nearest band edge, $nf^2/tol^2$; in-band → skipped |
| `c` | inside | centre, **reduced** $nf^2/bw^2$ |
| `c` | outside | nearest band edge, $nf^2/tol^2$ |
| `r` / `c` | no sim yet (first iteration) | centre, $nf^2/tol^2$ (conservative) |
| any | multi-environment (K > 1) | folded **per environment**, against that environment's own simulation; K folds, one per roster entry |

Folding is per quantity: front R/T/A demands → `(targets_r, weights_r)` /
`(targets_t, weights_t)` / `(targets_a, weights_a)`; back demands → the
`rb`/`tb`/`ab` siblings; absorption (`As`/`Ap`/`Au`, `ABs`/`ABp`/`ABu`)
derives $A = 1 − R − T$ from the companion pairs — the kind table above
applies identically in every bucket. Missing companions/rows fold
conservatively (exact at centre).

Phase demands fold to one `(targets, weights)` pair per S-matrix channel
(`phi[0..=3]` → one `P_PHI` call each). Intensity/absorption demands
require linear normalization; phase demands accept linear or phase
(wrapped residuals mirror the evaluator); anything else is an `Err`.
Phase demands must carry raw values with `norm_factor == 1` (the phase
arm scales nothing — the converter unscales the resolved triple).

Spectral-label mapping for the converter (`TargetCollection` → `MeritSpec`):
`R/T/A` × s/p/u → front demands, `RB/TB/AB` × s/p/u → back demands;
`phase=True` targets become phase demands on the mapped curve's element
(R → r_front, T → t_fwd, RB → r_back, TB → t_back). Anything else raises
`ValueError`. Angular targets expand to one single-point demand per angle.

## Environments (multi-environment designs)

One coating, several surroundings. The same films sit under a bare
surface, under a laminate, behind a cover glass; each of those is an
**environment**, and the demands that describe it are tagged with its
name:

```python
tc.add(SpectralTarget(wl, R_bare,  tol, 0.0, "s", "R", environment="bare"))
tc.add(SpectralTarget(wl, R_lam,   tol, 0.0, "s", "R", environment="laminated"))
spec = build_merit_spec(tc, environments=["bare", "laminated"])
```

**The tag is a name, resolved once.** `environment=` takes a name from
the roster passed to `build_merit_spec`, and the resolution happens there
— the one place a demand and the roster are both in scope. An unknown
name refuses, quoting the demand and the known environments, because a
typo has no shape error of its own: it would silently resolve to some
other environment and the run would report a merit for a coating nobody
described.

**After resolution the roster is gone, so ORDER is the binding — pass
the collection, not a pre-built spec.** Resolution turns a name into an
index and the compiled `MeritSpec` keeps only `n_envs`, not the names.
The run door can therefore check that the spec and the design agree on
the *count* of environments, and nothing more. A roster written in a
different order in the two places passes that check and silently scores
every demand against the wrong surroundings:

```python
spec = build_merit_spec(tc, environments=["laminated", "bare"])   # order A
run_needle(design=..., environments=[bare_env, lam_env], targets=spec)
#                                    ^ order B — accepted, wrong answer
```

A *typo* does refuse, and refuses early, at `build_merit_spec`, because
the tag is not in the roster it was handed. A *permutation* has no such
shape error: both names exist, just at swapped indices.

The safe spelling is to hand `run_needle` the `TargetCollection` itself
and let it build the spec — it passes the request's own roster, so the
two cannot disagree. Reuse a pre-built spec only across runs that pass the same
`environments=` list in the same order. Making the run door compare
names rather than counts (a roster on `MeritSpec`) is open work, not a
property of the current build.

**Absent means the first environment.** Every target set written before
environments existed is therefore still meaningful and its JSON is
unchanged — the key is emitted only when the tag is set. On a
single-environment run the roster is still a roster and its one entry is
called `"default"`, so a tagged demand on a flat run gets a refusal that
can say what the name *could* have been.

**The merit is one number over all environments.** Residuals are ordered
environment-major: every demand of environment 0, then every demand of
environment 1. There is no per-environment merit to trade off — a joint
design is one optimization whose residual vector happens to span several
coatings, and weights are the only instrument for saying one surface
matters more than another.

**Parameter identity is the film NAME, not the position.** A design
segment shared by K environments is one set of thicknesses, and what ties
a row in one assembly to the same parameter in another is the film's name.
Two films of the same physical material are therefore two names carrying
identical tables — a repeated name is one parameter spelled twice and is
refused. This is also why a `material_code` in a program document's
`design:` section is both the material and the identity: a `LayerRow` has
no second field for a name.

**The needle sees all K at once.** Each environment scans its own
assembly with its own fold (the coverage row above); a site that lands in
the shared design is routed by design parameter into a shared bucket and
summed. The candidate the sweep inserts is the one that helps the joint
merit, not the one that helps environment 0 — and the insertion edits the
shared design and every environment's template at the same time, at each
environment's own row.

**Cost is ×K, in solves and in folds.** One evaluation assembles K stacks
and solves K of them; the needle sweep folds K times per cycle. The
scaling is linear in the number of environments and it is *not* amortized
— there is no shared factorization between environments, because a
different surrounding is a different stack, not a perturbation of one.
Two environments cost about twice one; six cost about six times one.

The cost is worse when a surrounding is *graded*, because a gradient span
expands to many rows and those rows are solved in every environment that
carries them. A 64-sublayer cover on three environments is 192 extra rows
per evaluation. Where the surroundings are thick and passive, an S-matrix
embedding (solve the surroundings once, reuse the result) would remove
most of that — it is a known follow-up, not part of v1.

**Surroundings are fixed, by two refusals.** A `layers` segment's
rows carry `optimize` and `needle` forced false, and an explicit `true`
refuses at compile. A `per_film_flags` override is keyed by material
code and applied *after* the row, so it could undo that forced false;
since 0.7.6 the assembler refuses it there instead, naming the
surrounding row it reached (`{env}.fixed[{seg}][{i}]`) rather than only
the material code the caller wrote. The refusal is worth its own
message because the failure it replaces was quiet: under K > 1 the LM
moved environment 0's copy of the surrounding against
environment-0-only residuals while every other environment kept the
compiled thickness — "one design, several surroundings" broken with no
error (review PB, M1; measured at 500.0 -> 534.0298 nm, with the shared
design film driven to its clamp floor).

`film_flags`, the global map, was never part of this and is untouched:
it is applied *before* the row, so the forced false wins on its own.
`film_flags={"optimize": True}` remains the ordinary way to say
"optimize the design" on a segmented run. If you meant to free a layer
that a surrounding currently holds, move it into a design segment —
then every environment shares it, which is what a design variable is.

**What a joint run gives up.** While K > 1 the thin-layer floor is a hard
bound on the optimizer rather than a post-hoc sweep, and the sweep clamps
up instead of removing. The reason is structural: eliminating a sub-floor
span would leave the compile carrying a parameter the design no longer
has. The visible consequence is that a rejected needle seed parks at the
floor instead of disappearing, so `thin_layer_policy="remove"` is refused
at the door on a multi-environment run rather than silently reinterpreted.

## Differential phase (`PDts`/`PDtp`)

`PDts`/`PDtp` demand the coating-induced transmitted phase: the design's
`arg(t)` minus the equivalent incidence-medium layer,

$$\Delta\varphi(\lambda) = \arg t(\lambda) - passes\cdot\frac{2\pi\,n_{inc}\,D\cos\theta_{inc}}{\lambda},$$

with `passes = 1` (single traversal; `passes = 2` covers a reflection
round trip if reflection labels are ever added), `D` the total coating
thickness and `n_inc` the real incidence index. GD/GDD (and TOD/FOD)
over $\Delta\varphi$ carry the reference's dispersion (PD3, decision:
option 1) — with $n_{inc}(\omega)$ dispersive, $\mathrm{GD}_{ref} =
(n_{inc} + \omega\,dn_{inc}/d\omega)\,D\cos\theta_{inc}/c$ is $\lambda$-dependent
and the reference's GDD contribution $2\,dn/d\omega + \omega\,d^2n/d\omega^2$
is non-zero, so finite differences of the corrected $\Delta\varphi$ (or
GD/GDD of the reference-only stack subtracted from the absolute keys)
are the right recipe. Before 0.6.45 a frozen scalar index contributed a
constant GD and exactly zero GDD, which is what made the old "finite
differences kill the reference anyway" remark true. The differential
form matters for absolute-phase targets and for correct needle-gain
bookkeeping; the absolute GD/GDD keys (`GD_T_s`, …) stay reference-free
— pinned by the PD3 parity twin.

Evaluation points (all in solver convention — see the sign note below):

- Ingestion forces phase normalization (raw radians, `nf = 1`); the
  converter maps `PDts` → (`Ts`, phase, passes 1), `PDtp` → (`Tp`, …).
  `phase=False` or a polarization mismatch raises `ValueError`.
- `SimCurves` carries `total_d`/`n_front_re`/`n_back_re` (defaults 0/1/1
  zero the reference); the thickness-optimizer evaluator fills them from
  the stack (the incidence/exit indices as per-λ columns — PD1; a
  dispersive medium now carries its real index at every λ, where the
  pre-0.6.45 scalar froze the centre λ and erred by
  `2π·D·cosθ·[n(λ) − n(λ_c)]/λ`) — but only when the spec asks
  (`uses_phase()` gates complex-row assembly, `uses_differential()` the
  metadata; intensity-only LM loops pay zero extra allocations).
  `total_d = 0` reproduces absolute phase bit-for-bit.
- The merit Phase arm subtracts the reference before wrapping; the fold
  passes PD demands to the `phi` buckets unchanged (same channel as the
  absolute element) and additionally accumulates the exact `dM/dD`
  correction `phi_gain_shift[ch] = Σ −2·kz·w·(s−rt)` (`kz =
  passes·reference_wavenumber`). Subtract it from the assembled `P_PHI`
  (`needle_gradient(…, gain_shift_phi=…)`); it is uniform in z, so the
  needle site (`argmax`) never moves — only predicted-gain bookkeeping.
- Manual-sim recipe (zero core involvement, doubles as test oracle):
  `apply_reference_rotation` multiplies complex rows by `e^{−i·ref}`;
  an absolute-phase demand on rotated rows is exactly a differential
  demand on raw rows.

Sign convention: the reference is `+kD`, matching this crate's
forward-propagation phase (an all-matched slab simulated by the solver
has `arg(tf) = +kD` — pinned by test). That is the conjugate of
Macleod/`e^{+iωt}` textbooks; the crate is self-consistent (absolute
phase demands, `P_PHI`, GD/GDD share it), so only textbook-imported
target numbers need conjugating — never solver-produced ones.

### The dropped `+1` level

Outside the band, `calculate_merit` contributes
$((|d|-bw)/tol)^2 + 1$ while the fold carries only
$w\cdot(R-t_{edge})^2 = ((|d|-bw)/tol)^2$. The $+1$ is a pure number:
independent of $R$, thicknesses, and everything the needle can perturb.
Hence needle gradients, the $P(z)$ profile, site ranking and the best-site
pick are **exact**; only the scalar merit value reads lower by exactly
$N_{outside}$ (the count of currently-violated `c` points):

$$M_{true} = M_{folded} + N_{outside}$$

Carrying the constant would need a third per-point array through the
`needle_gradient` FFI and every call site, for a term that cannot change
any decision — so the fold drops it by design. If a future consumer needs
true values from the folded path (e.g. a trust-region acceptance test),
return the outside-count alongside `(targets, weights)`; the fold already
knows it.

### Overlap constant (same class of loss)

When several demands fold onto one solver point, completing the square
$\sum w_i(s-t_i)^2 = W(s-\bar t)^2 + [\sum w_i t_i^2 - W\bar t^2]$ keeps
$(\bar t, W)$ and drops the bracket. Gradients stay exact; folded merit
*values* under-read by that constant. Exact value identity holds only for
non-overlapping demands (at most one demand per solver point per bucket).
Like the $+1$, this never affects needle gradients, site ranking, or
insertion decisions — only human-readable merit comparisons.

### Multiblock (incoherent) needle

`Pmb` differentiates the cascade totals through the adjoint weights
($g[p][e] = \partial v[e]/\partial param$): R via $v[0]$, T via $v[2]$,
A via $1 − v[0] − v[2]$ with negated summed weights. Request with
`P_MB`/`P_MB_T`/`P_MB_A` (reusing the R/T/A target arrays); phase has no
multiblock path (phases only live inside coherent blocks).

### Kinks

Value-continuity holds at every band edge (both sides read $1.0$ for `c`,
$0.0$ for `r`), but the *slope* has a kink there. A needle step that
crosses an edge is evaluated with the pre-step linearization and
self-corrects on the next re-fold; persistent oscillation around an edge
means the step size is too large — damp it.

## Color demands (`ColorTarget`, v1)

One demand = front R/T intensity curve × angle × illuminant × observer ×
quantity × distance × weight. Defaults D65 + 1931 2°, all overridable;
names (`"D65"`, `"1931_2deg"`) resolve to the embedded defaults, anything
else must be an explicit table (no registry to drift).

- **Quantities/distances:** Lab|XyY (P1), LCh|Oklab|Y (P2),
  DomWl|sRGB|Luv|XYZ|Din99|White|Yellow (P3) ×
  DeltaE2000|DeltaE76|Channels, with the compat matrix enforced at
  compile (DeltaE2000|DeltaE76 live in Lab|LCh only — everything else
  takes Channels; Y|Yellow are scalar Channels; DomWl|White are pair
  Channels). Lab white is always the demand illuminant's own white
  (computed at compile from native tables, never adapted).
- **Kinds/transforms:** `Exact` + `linear` only. No `integral` /
  `count_norm` / `phase` / `band` — a color demand already is integral
  (one residual per demand, not per wavelength).
- **Channels tolerances are unit and fixed (v1):** there is no
  per-channel tol — every Channels residual divides by 1.0 in native
  channel units, and `weight` scales the whole demand only. Effective
  weighting this implies: Lab/LCh/sRGB channels are comparable-by-luck
  (L ~ 0-100 vs a/b, hue degrees, rgb units); XyY/XYZ/Y compare
  commensurate 0-1 units; White W (0-100) vs Tw comparable-ish. The
  exception is **DomWl: wavelength (nm) vs purity (0-1) equally
  weighted, so 1 nm ≡ 1.0 purity — purity is ~100x mute in practice
  and the demand overwhelmingly optimizes wavelength.** A purity-first
  DomWl demand is not expressible in v1 (per-channel tol is parked as
  a non-breaking additive follow-up); use `weight` for demand-level
  balance only.
- **Wavelength window:** optional `wavelength_range: [lo, hi]` nm per
  demand restricts the *sample* integral only — the white stays the
  full-illuminant white and the DomWl locus stays full-CMF, so windowed
  numbers stay comparable to full-range numbers (`None` = full overlap).
  `lo < hi` + finiteness refused natively; a window missing the grids
  errors naming the range. The needle fold inherits it (covered points
  deposit, the rest keep zero buckets) with zero fold changes.
- **Merit path:** missing curve fails the whole key group (standard
  missing-penalty path, parity with pointwise); empty table/sim overlap
  errors instead of skipping (a demand that sees nothing is a spec bug).
  Coverage note: the integral runs over the OVERLAP SUBSET of sim grid
  and demand tables — partial overlap narrows the integral silently
  (deliberate: narrow-window twins use this on purpose). A typo'd table
  range therefore yields a finite, plausible, wrong merit; when a demand
  looks dead, check table-vs-sim coverage first.
  **This is the opposite of the pointwise rule, and the asymmetry is the
  porter trap.** A pointwise `SpectralTarget` whose grid does not line up
  with the sim grid is *interpolated* onto it (aligned fast path, else a
  two-pointer linear interpolation), so every target point still scores; a
  colour demand whose tables do not line up is *narrowed*, so the points
  outside the overlap silently score nothing. Same mismatch, two different
  answers. Both are right for their own maths — an interpolated illuminant
  or CMF table would quietly change the colorimetry, which is why the
  colour path refuses to guess — but code ported from the pointwise side on
  the assumption that "the grids get sorted out" gets a narrower integral
  and no warning. Empty overlap is the one case that errors.
- **Needle fold (Option B):** each demand deposits its analytic gradient
  g(point) = dF/dcurve (weight, residual, and the U-curve half included)
  into the `grad_r`/`grad_t` buckets — not a (target, weight) pair. `Ru`
  / `Tu` split half per polarization branch (Ru is the Rs/Rp mean); s/p
  demands ride the shared bucket under the same convention as pointwise
  targets. Zero points are skipped in the hot pass, so color-free sets
  are bitwise unaffected.

### Extended quantities (P2: LCh | Oklab | Y)

- **LCh:** references convert LCh→Lab once per eval; DeltaE is measured
  in Lab. Channels compare (L, C, h) with the hue difference wrapped to
  [−180, 180] *before* scaling (179 vs −179 is 2 deg, not 358).
  Near-neutral colors (C≈0) have ill-defined hue — the kernel gradient
  stays finite (FD class) but directionally unstable; weight hue
  tolerance generously near the achromatic axis.
- **Oklab:** D65-defined. Non-D65 demands Bradford-adapt XYZ→D65 first
  (same `adapt` the bindings use, unclipped so the FD gradient stays
  smooth). DeltaE on Oklab is refused — equal-tol Channels *is*
  unweighted Euclidean. Note: the in-tree D65 white constant and the
  Oklab matrix disagree ~1e-4 in b (pre-existing); systematic, negligible.
- **Y:** scalar reference, single residual off `tol[0]` (= 1.0 unit —
  see the unit-tol note above; the AR classic: dark residual, hue free).

### Dominant wavelength + purity (P3: `DomWl`)

- **Reference** is a `[wavelength_nm, purity]` pair; residual is
  Channels-style off `tol[0]` (nm) / `tol[1]` (purity) — both unit (see
  the unit-tol note above: purity is effectively ~100x mute). DeltaE
  refused.
- **Geometry:** forward ray (own-white → sample) vs the demand-locus
  polyline (monochromatic chromaticities of the demand CMF): a hit past
  the sample is spectral (`purity = 1/t`). A miss means a purple-direction
  ray; the backward extension hits the locus (complementary branch,
  purity **negative** — explicit in the value, never a hidden mode).
- **Achromatic samples** (|d| < 1e-9 in xy) return (0, 0): hue carries no
  information at white, purity does the work. Any rule kinks here
  (lambda is undefined at white) — same family as band-edge kinks; pair
  hue demands with a purity/Y floor, or use Lab near neutrals.
- Mixed-branch refs (spectral op vs complementary ref) stay continuous
  in purity but kink in lambda at white — pick refs in the right family.

### Coordinate quantities (P3: `sRGB` | `Luv` | `XYZ`)

- All three take Channels only (DeltaE76/2000 are Lab-space; Luv's own
  dE_uv is out of scope — Channels on (L,u,v) covers the need).
- **`sRGB`:** D65-adapted, gamma-encoded, UNCLIPPED. In-gamut values are
  exact standard sRGB; out-of-gamut extends linearly (the transfer's
  linear branch below 0.0031308 keeps it smooth and monotone). Display
  clipping is presentation — demands compare unclipped, so gradients
  never vanish outside the gamut.
- **`Luv`:** under the demand illuminant's own white (like Lab).
- **`XYZ`:** raw tristimulus, linear sensor-space matching.

### Index quantities (P3: `Din99` | `White` | `Yellow`)

- **`Din99`:** graphics DIN99 coordinates (ke = kch = 1) of Lab under the
  demand white; Channels compare. DIN99 is already euclidean — DeltaE
  refused (equal-tol Channels *is* unweighted Euclidean).
- **`White`:** CIE whiteness `[W, Tw]` pair ref (W on 0-100; Y here is
  0-1, hence x100). The formula needs only XYZ + demand white —
  illuminant-agnostic by construction. Perfect diffuser hits exactly
  (100, 0). No CIE validity-box enforcement (documented scope).
- **`Yellow`:** ASTM E313 scalar `100*(Cx*X - Cz*Z)/Y` off `tol[0]`
  (= 1.0 unit).
  Coefficients ride the demand (`yi_cx`/`yi_cz`, defaulting to the E313
  D65/10-deg table 1.3013/1.1498) — other geometries pass their own
  table values explicitly, never silently reuse D65/10-deg.
