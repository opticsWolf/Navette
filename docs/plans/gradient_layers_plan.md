# Gradient-index & inhomogeneous layers — implementation plan

## D0. Vocabulary (locked for this plan)

- **Inhomogeneous (exists, extended by this plan):** ONE material whose index
  drifts through the layer. `Layer { inhomogen: bool, inh_delta: f64 }` —
  expansion scales the base nk by a linear factor `1−δ … 1+δ` over
  `sub_layer_count()` slices (`thickness^0.4`-scaled). Phenomenological
  deposition drift. The scaling physics is FROZEN (bit-identical legacy
  path); what this plan adds is the application mode for δ:
  - **(a) Fixed (default, current behavior):** δ comes from `inh_delta`
    (+ group summand, §D2) — independent of thickness.
  - **(b) RateCapped:** δ grows with thickness,
    `δ_layer = min(rate · thickness / ref_thickness, cap)`
    (e.g. 0.05 per 100 nm capped at 0.30: a 50 nm film drifts ±0.025,
    a 600 nm film drifts ±0.30, saturated). Same ratio/delta
    canonicalization as gradients (§D0, gradients): one slope key.
- **Gradient-index (new):** a TWO-material mixture whose volume fraction
  `f(z)` varies through the layer, evaluated per sublayer with the native
  EMA kernels (`materials/ema.rs`: Bruggeman, MaxwellGarnett, Looyenga,
  Lichtenecker — all in-tree, all bound). Two sub-types:
  - **(a) FixedSpan:** `f` runs `f_start → f_end` linearly in physical `z`,
    independent of thickness. The spiritual successor of `inh_delta` for
    mixtures: 50 nm or 500 nm, endpoints are 0.1 → 0.8.
  - **(b) RateCapped:** `f(z) = f_start + rate · (z / ref_thickness)`,
    clamped to `[f_min, f_max]`. The slope is thickness-relative (e.g.
    0.3 per 100 nm); a 400 nm film saturates into a flat pure-material
    tail where the cap binds. `ref_thickness` user-configurable
    (default 100 nm), caps default `[0.0, 1.0]`.
- **Ratio vs delta (user wording, canonicalized):** both reduce to one
  slope + caps. `rate` IS the per-reference-thickness change of the mix
  fraction; there is no second degree of freedom. If a config supplies
  two slope-like keys for one layer → refuse (`"gradient: 'rate' and
  'delta' are the same slope — give one."`). No silent precedence.

## D1. Current-state inventory (what exists, with addresses)

- `structure/layer.rs`: `inhomogen`, `inh_delta` (default 0.1),
  `sub_layer_count()` (`thickness^0.4` rule), serde round-trip, state keys.
- `structure/expansion.rs` (~l.233): graded branch — factor profile,
  error application, inversion reversal, roughness only on first sublayer.
  Interface slices already use `looyenga_mix` of two nk spectra.
- `smatrix/synthesis/structure.rs::from_design`: non-background graded
  films HOMOGENIZE (base index, one row) with a warning; background-pinned
  graded films expand WITH profile (optimize/needle forced false).
  Merge is nk-keyed (distinct nk never collapses); thin-removal and
  needle/LM are flag-guarded.
- `materials/ema.rs` + `materials/__init__.py::evaluate`: homogeneous EMA
  composites via nested `MaterialSpec(host/inclusion, f)` — scalar `f`
  only, no spatial variation. **This is the kernel pool the gradient
  reuses; no new mixing math.**
- `smatrix/synthesis/design_config.rs` (l.84): `inhomogen` row flag —
  the gradient block docks next to it.
- Solver sees only uniform slabs (`SolverArrays` rows). **Zero solver /
  needle-operator change in this plan** — all work is expansion +
  pipeline guards. Gradient-span rows carry `optimize=False,
  needle=False` (background rule), so needle/LM skip them automatically.

## D2. Data model (Rust-first, single source)

```rust
pub enum EmaModel { Bruggeman, MaxwellGarnett, Looyenga, Lichtenecker }
pub enum GradientMode {
  /// (a) f_start → f_end over the layer, any thickness.
  FixedSpan { f_start: f64, f_end: f64 },
  /// (b) f(z) = f_start + rate·(z/ref_thickness), clamped to [f_min, f_max].
  RateCapped { f_start: f64, rate: f64, ref_thickness: f64, f_min: f64, f_max: f64 },
}
pub struct GradientSpec {
  pub material_a: String,   // f = 0 endpoint; ALWAYS the host (§D9.4)
  pub material_b: String,   // f = 1 endpoint (f = fraction of B in A)
  pub ema: EmaModel,        // default Bruggeman (§D9.1, locked)
  pub mode: GradientMode,
  pub shape: ProfileShape,  // Linear only v1; enum-ready (§D9.2, locked)
  pub sublayers: Option<u32>, // override; default = resolution rule §D3
}
/// Intra-span profile shape. V1 ships Linear only; SCurve/Exponential
/// (or similar) arrive as new variants — no schema break, `mode` and
/// `shape` are already enums. Each new variant gets its own twin batch.
pub enum ProfileShape { Linear }
// Layer gains: pub gradient: Option<GradientSpec>  (default None)
// Layer gains: pub inh_mode: InhMode                  (default Fixed)
/// Application mode for the legacy scaling drift (physics frozen).
pub enum InhMode {
  /// Current behavior, bit-identical: δ = inh_delta.
  Fixed,
  /// δ_layer = min(rate·thickness/ref_thickness, cap).
  RateCapped { rate: f64, ref_thickness: f64, cap: f64 },
}
```

Inhomogeneous validation + combination order (native):
- `rate` finite (sign free — negative inverts the drift direction);
  `ref_thickness > 0`; `cap` in `[0, 1]`.
- Combination (both modes, unchanged legacy arithmetic):
  `δ_nominal = (δ_layer + inh_delta_summand) * 0.5`, where `δ_layer`
  is `inh_delta` in Fixed mode, the rate formula in RateCapped mode.
  In RateCapped mode only, `δ_nominal` is clamped to `[0, cap]` — the
  cap binds the *nominal* value; the stochastic error draw applies
  after, unclamped (clamping draws would bias Monte-Carlo statistics;
  noise is allowed to exceed the nominal cap, documented).
- `sub_layer_count()` feeds δ through the same path: in RateCapped
  mode the count uses `δ_layer(thickness)` (pure function of thickness
  — no circularity, evaluated once).
- Fixed mode is byte-for-byte the legacy path (existing differential
  pins keep passing unmodified — the proof).

Validation (native, at parse + at expand):
- `f_*` in `[0,1]`; `ref_thickness > 0`; `f_min ≤ f_max`, both in `[0,1]`;
  `rate` finite (sign free — negative slopes allowed).
- `material_a != material_b` (refuse — a gradient of a material with
  itself is a spec bug; single-material drift is `inhomogen`).
- `gradient` + `inhomogen` on one layer → refuse (two profile engines;
  message names both and points at the alternative).
- `inh_delta_summand` / inh-errors are IGNORED with a warning when a
  gradient is set (error engine is factor-based; mixtures need no δ).
  Thickness/roughness errors apply normally (they scale sublayer
  geometry, untouched).
- Both materials must resolve in the SAME provider at expansion
  (`DictProvider` entries / `SpecProvider` names / weaver keys);
  missing → `Err` naming the layer + the absent material (never fall
  back to base index — a half-mixture is worse than a refusal).

Profile shape is LINEAR in physical z for v1 (documented; S-curve/
exponential shapes reserved — one enum variant later, no schema break
since `mode` is already an enum).

## D3. Expansion (the only physics code)

New branch in `expansion.rs` beside the `inhomogen` branch (mutually
exclusive by §D2 validation):

1. Resolve `nk_a(λ)`, `nk_b(λ)` from the provider (same grid assertion
   as every other lookup — length-only is silent-wrong, §grid rule).
2. Sublayer count: default rule
   `clamp(ceil(thickness / max_step), 3, 64)` with
   `max_step = min(ref_thickness_or_thickness/8, λ_min/10 in optical
   terms…)` — exact formula fixed at implementation, pinned by a
   differential test over randomized thicknesses (same technique as the
   `sub_layer_count` pin: integer boundaries must agree bit-exactly).
   `sublayers: Some(n)` overrides (clamped to `[2, 256]`, warning
   outside).
3. Per sublayer `i` (midpoint `z_i`): `f_i` from the mode formula
   (RateCapped clamped), `nk_i(λ) = EMA(nk_a, nk_b, f_i)` via the
   layer's `EmaModel`, `d_i = thickness / n` (last sublayer absorbs
   float remainder so Σd ≡ thickness exactly).
4. Inversion (`inv`) reverses row ORDER (same as legacy — physical flip),
   which equivalently swaps the profile direction; tested both ways.
5. Roughness/type on first sublayer only (legacy convention, keeps
   interface physics identical).

Type-(a) thickness-independence is structural: endpoints depend only on
`f_start/f_end`, never on thickness (only the COUNT changes). Type-(b)
saturation is structural: the clamp binds → flat tail rows of pure
`material_b` (or `_a` for negative rates). Both asserted in §D6.

Inhomogeneous modes share the branch: `δ_layer` is resolved first
(`inh_delta` vs rate formula), everything downstream — factor profile
`1−δ … 1+δ`, error draw, inversion, first-sublayer roughness — is the
frozen legacy code path, parameterized by one float. Fixed mode must
produce bitwise-identical rows to pre-plan builds (differential suite
is the gate, not review).

Cost note: expansion runs once per design AND once per error-draw
(thickness errors re-expand). Gradient EMA is `n_sublayers × n_λ`
complex ops per draw — trivial vs the solve. No memo needed; if the
bench (§D6) disagrees, memoize on `(nk_a, nk_b, f)` per draw.

## D4. Pipeline guards (span provenance — the one structural change)

`from_design` currently drops span identity after building rows. Gradient
spans need three protections no flag alone provides:

1. **Thin-removal exemption:** `remove_thin_layers` deletes thin rows;
   gradient sublayers ARE thin by construction. Carry `span_id` per row
   (`DesignStack` gains `spans: Vec<Span>` alongside rows — `Span`
   already exists in `expansion.rs`) and exempt non-singleton-span rows
   from thin-removal. (Also hardens legacy graded spans — same bug.)
2. **Needle-refusal:** inserting a needle seed INSIDE a gradient span
   would orphan the profile. `insert_layer`/`needle` hosts only rows
   whose span is a singleton uniform film; gradient spans refuse with
   the span's material + range named. (Interface slices already force
   `optimize/needle=False` — same rule, one place.)
3. **Homogenize path:** non-background gradient films CANNOT use the
   base index (there is none). Homogenize to `EMA(f_mid)` —
   `f_mid = (f_start+f_end)/2` for (a), `clamp(f(start…end) mean)` for
   (b) — single row, warning names the mixture + `f_mid`. Background
   rule unchanged: full profile, pinned, silent.

Merge needs NO change (nk-keyed — sublayer nk differ → preserved;
verified by test, not by reasoning).

## D5. Config + Python surface (thin, last)

- `design_config.rs::LayerRow` gains `gradient: Option<GradientJson>`
  (mirror struct + defaults; serde). `builders.py` passes through;
  `Layer` dataclass gains `gradient: GradientSpec | None` (frozen
  holder, validated natively — pydantic stays deleted).
- **Schema version:** new optional fields change state shape →
  `SCHEMA_VERSION 1 → 2` with default-`None` migration: `gradient`
  defaults to `None`, `inh_mode` to `Fixed`, `shape` to `Linear` when
  the keys are absent — so every v1 file parses unchanged and
  reproduces bit-identical results with zero user action (no file
  rewriting, no converter). The fingerprint test is extended to cover
  the new fields. No silent additive change — the fingerprint
  discipline exists for exactly this.
- Python `MaterialSpec` EMA nesting is REUSED for endpoint evaluation
  documentation (users already know host/inclusion/`f`); no new
  material classes.

## D6. Verification (twins per phase, §8-style)

- **G-expansion:** per-row nk vs direct `_native.ema_*` oracle on the
  same `(nk_a, nk_b, f_i)` — BITWISE (same kernel, same inputs).
  Profile test: extracted `f_i` (invert Looyenga? no — recompute
  expected `f` from the mode formula and compare EMA outputs) matches
  formula to 1e-15; caps produce exactly-flat tails (bitwise equal
  pure-material rows).
- **G-thickness:** type (a) at two thicknesses → endpoint rows bitwise
  equal, counts differ. Type (b): 200 nm film, start 0.1, rate 0.3/100
  nm → end 0.7; same spec at 400 nm → saturates to 1.0 from z=300 nm
  (flat tail = pure-`b` rows, bitwise).
- **G-inh-modes:** Fixed rows bitwise vs legacy oracle (existing pins +
  randomized differential); RateCapped δ formula vs hand-computed
  `min(rate·t/t_ref, cap)` over a thickness sweep (saturation knee exact);
  double-thickness → double-δ below the cap; error draws statistically
  unclamped (mean preserved across the clamp boundary, documented).
- **G-pipeline:** `simulate()` with background gradient vs manually
  built explicit sublayer stack → bit-equal spectra (proves the guards
  preserve physics end to end).
- **G-guards:** thin-removal keeps spans; needle-insert into span
  refuses naming material+range; homogenize warning names mixture +
  `f_mid`; both-materials-missing refuses naming the absent one;
  `gradient`+`inhomogen` refuses naming both.
- **G-differential:** sublayer counts pinned over randomized
  thicknesses (integer-boundary agreement, legacy technique).
- **G-bench:** expansion cost per draw with/without gradient
  (memo decision data).

## D7. Progression (after color P1 — no parallel schema churn)

- **G1 (0.4.28):** §D2 model + §D3 FixedSpan expansion + G-expansion/
  G-thickness twins, PLUS the inhomogeneous modes (§D2–D3, §D6
  G-inh-modes). Both are expansion-only changes sharing one twin
  harness — one patch, no pipeline change (gradient films homogenize
  via §D4.3 formula already in this patch — safe default).
- **G2 (0.4.29):** RateCapped mode + caps + saturation twins.
- **G3 (0.4.30):** §D4 span provenance (spans on DesignStack,
  thin-exemption, needle-refusal) + G-guards + G-pipeline twins.
- **G4 (0.4.31):** §D5 config/schema-v2/Python surface + full
  validation green.
- Each: implement → gates → twin batch → commit/push dev → ff main.

## D8. Risks

| Risk | Mitigation |
|---|---|
| Sublayer-count float boundaries diverge | differential pin test (legacy precedent) |
| EMA per-draw cost surprises | G-bench before/after; memo fallback ready |
| Users conflate `inh_delta` with gradient | refusal message cross-references; docs §D0 |
| Legacy `inhomogen` deprecation pressure | explicitly NOT deprecated — single-material drift is real (oxidation gradients, nitrides); coexistence rule §D2 |
| Schema bump churn vs color P2 (0.4.27) | G-series starts after 0.4.26 lands; P2 is quantity-only, no schema change |

## D9. Decisions (all locked 2026-09-06)

0. Inhomogeneous stays simple scaling; modes (a) Fixed default /
   (b) RateCapped — cap binds the nominal δ pre-error-draw (§D2).
   `inhomogen` + `gradient` on one layer → refused (two profile
   engines; message cross-references the alternative).
1. Default EMA model: **Bruggeman** (right for co-sputtered mixtures;
   symmetric — no host question until MG-family models are chosen).
2. **Linear-only profiles v1**, enum-ready: `shape: ProfileShape`
   ships with the single variant `Linear`; curved profiles (S-curve,
   exponential…) arrive as new variants — no schema break, each with
   its own twin batch.
3. `SCHEMA_VERSION → 2` with default-`None` migration (see §D5).
4. **Asymmetric EMA (MaxwellGarnett, MoriTanaka): `material_a` is
   ALWAYS the host, `f` is the volume fraction of B in A.**
   Effect: `f = 0.1` means 10% B inclusions in 90% A host — swapping
   A↔B is NOT the same physics (MG is host/inclusion-asymmetric by
   construction) and silently swapping would corrupt fits. This remark
   MUST be repeated in the `GradientSpec` docstring at implementation
   (G1) so the choice is visible at the call site, not just here.
   **Deferred follow-up:** revisit whether to allow per-layer host
   selection (e.g. an explicit `host: A|B` key or auto-switching near
   f=0.5). Options deferred until a real asymmetric-mixture use case
   lands — do NOT design it now.
   (Host-effectiveness comparison + full docstring checklist: needle
   plan §§N7–N8 — the G1 implementing patch answers to both.)
