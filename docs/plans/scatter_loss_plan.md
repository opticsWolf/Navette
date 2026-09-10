# Plan — Debye-Waller / TIS scatter-loss option

**Status:** proposed (not yet scheduled)
**Motivation:** every existing roughness model (`RoughnessType` 0–5) is a
**specular-only** treatment. NEVOT_CROCE (type 5) conserves specular energy by
construction (`R+T=1`, graded-interface picture); the graded types (1–4) damp
both beams but their `R+T<1` deficit is a form-factor artifact, not a derived
scatter loss. **None of them model diffuse scatter as a tracked quantity.** A
topographically rough interface scatters energy out of the specular beam into
the diffuse hemisphere (Total Integrated Scatter). This plan adds a physically-
grounded Debye-Waller / TIS option and — critically — a dedicated scatter
channel so the lost energy is reported rather than silently mislabeled.

References: Beckmann & Spizzichino, *The Scattering of EM Waves from Rough
Surfaces* (1963); Bennett & Mattsson, *Introduction to Surface Roughness and
Scattering* (OSA, 1989); Stover, *Optical Scattering* (SPIE, 3rd ed.); Névot &
Croce, Rev. Phys. Appl. 15, 761 (1980); de Boer, Phys. Rev. B 49, 5817 (1994).

---

## 0. The accounting problem this must solve first

The engine computes absorptance as a **residual**:

```
solver.rs:482   A_s = 1 - rs - ts
solver.rs:483   A_p = 1 - rp - tp
solver.rs:484   A_avg = 1 - 0.5*(rs+rp) - 0.5*(ts+tp)
```

So *any* energy missing from the specular beams lands in the `A_*` channel. If a
Debye-Waller model damps both `r` and `t` (as it must, to represent scatter
loss), that deficit will be reported as **absorption** — which is wrong and
actively misleading for a lossless-but-rough interface. Therefore a scatter
model is **not** allowed to ship without a dedicated scatter channel that
separates diffuse loss from true bulk absorption:

```
1 = R_spec + T_spec + A_bulk + S_diffuse
```

Two ways to split the residual `(1 - R - T)` into `A_bulk` and `S_diffuse`:

- **(A) Model-derived TIS (recommended).** The Debye-Waller factors give the
  per-interface specular attenuation in closed form, so the scattered fraction
  is computable directly from `sigma`, `kz`, and the smooth-interface `R0/T0`.
  Sum per-interface TIS into `S_diffuse`; then `A_bulk = (1 - R - T) - S_diffuse`.
  This is exact within the model and needs no field integration.
- **(B) Field-integration absorption.** Compute true bulk absorption from the
  internal field (Poynting divergence) and get scatter as the remainder. Much
  larger lift (the engine does not expose per-layer absorbed power today) and
  unnecessary if (A) is accepted. Deferred.

**Decision: (A).** The whole point of a factor model is that the loss is
analytic — expose it.

---

## 1. Physics — what the Debye-Waller factors are

For a single interface with RMS roughness `sigma`, wavevector normal components
`kz1 = (2*pi/lambda)*n1*cos(theta1)` and `kz2` in the two media:

**Specular reflection (amplitude):**
```
r~ = r0 * exp(-2 * kz1^2 * sigma^2)              (Debye-Waller, incident medium)
```
so `R_spec = R0 * exp(-4 * kz1^2 * sigma^2)`. Note this is the **single-surface
Debye-Waller** form (uses `kz1^2`), which differs from Névot-Croce's
`exp(-2*kz1*kz2*sigma^2)` (uses the cross term `kz1*kz2`). NC → 0 loss as
`kz → 0` near the critical angle; DW does not. Which is "correct" depends on the
regime (DW for optical scatter from a single surface; NC near total external
reflection in XRR). **Offer both reflection conventions behind a flag** rather
than hard-coding one (see §3, `RScatterConvention`).

**Specular transmission (amplitude):**
```
t~ = t0 * exp(-((kz1 - kz2) * sigma)^2 / 2)       (transmission Debye-Waller)
```
This is the **negative-exponent sibling** of the NC transmission factor: NC
*enhances* transmission (`exp(+...)`) to conserve specular energy; DW *damps* it
(`exp(-...)`) because the scattered energy genuinely leaves the specular beam.
This one line is the physical crux distinguishing the two models.

**Total Integrated Scatter (the reported diffuse loss):**
```
S_R = R0 * (1 - exp(-4 * kz1^2 * sigma^2))        # scattered out of spec. reflection
S_T = T0 * (1 - exp(-((kz1-kz2)*sigma)^2))        # scattered out of spec. transmission
```
For reflection at small `kz*sigma`, `S_R/R0 ≈ (2 kz1 sigma)^2 = (4*pi*sigma*cos(theta)/lambda)^2`
— the textbook TIS formula, a good sanity anchor for the validation.

**Energy bookkeeping (must hold to machine precision for lossless media):**
```
R0 + T0 = 1     (smooth Fresnel)
R_spec + T_spec + S_R + S_T = 1   for a single lossless interface
```
Unlike NC, this is exactly unitary by construction (the scattered fractions are
defined as the complement of the specular factors), so there is no
high-contrast overshoot to apologize for.

---

## 2. Multi-interface subtlety (must be designed, not discovered)

The engine builds the stack with Redheffer star products over per-interface
modified Fresnel coefficients (`coherent_block.rs`). For NC/graded types the
per-interface amplitude factor folds cleanly into the recursion. For scatter:

- **Specular amplitudes** are still just damped Fresnel coefficients — they fold
  into the existing recursion exactly like NC (one extra `exp` per interface, on
  both `r` and `t`). No structural change to the solver hot loop.
- **Scatter accounting is the new part.** `S_R`/`S_T` are per-interface and must
  be accumulated across the stack accounting for the specular throughput that
  reaches each interface and the fraction that returns. First cut: report the
  **total** diffuse loss as the residual split `S_total = (1 - R - T) - A_bulk`
  where `A_bulk` is computed from the *bulk* extinction only (k>0 layers), using
  the same per-point machinery that already produces `tau = exp(-2*Im(beta))`
  for incoherent spacers. This sidesteps per-interface scatter re-interception
  (a second-order effect) for v1 and keeps the invariant `R+T+A_bulk+S=1` exact.
- **Coherent vs incoherent, s vs p:** the damping is per-polarization already
  (kz is pol-independent but enters through `cos(theta)`); reuse the existing
  dual-solver path. Incoherent spacers: the scatter factors apply at the
  interface, the `tau` intensity decay is unchanged.

Document the v1 scope explicitly: **single-scatter (no re-interception of
diffusely scattered light), specular beam fully modeled, diffuse light reported
as an integrated fraction, not resolved into angle (no BRDF).**

---

## 3. API surface

### Rust engine
1. `RoughnessType::DEBYE_WALLER = 6` in `structure/types` (Python enum) and the
   Rust match arms in `coherent_block.rs` (both the single-pol and dual sites,
   mirroring the type-5 edit from R1.1).
2. New per-interface branch computing `(r*fr, r*fr, t*ft, t*ft)` with the DW
   amplitude factors, plus emitting the per-interface `(S_R, S_T)` into a new
   accumulation buffer threaded through the block solvers.
3. `enum RScatterConvention { DebyeWaller, NevotCroceCorrelated }` config on the
   solver (default `DebyeWaller`) selecting `kz1^2` vs `kz1*kz2` in the
   reflection factor — the one genuine physics choice.

### Observable channels (the accounting fix)
4. New `Request` bits: `SCATTER_R`, `SCATTER_T`, `SCATTER_AVG` (and update the
   `expected_keys` map + the Python `Request` IntFlag — guarded by the R2.5
   request-bit sync test once that lands, so the new bits can't desync).
5. **Redefine the absorption channel semantics:** keep `A_* = 1 - R - T` as the
   *total non-specular* residual for backward compatibility, but add
   `A_BULK_*` = `A_* - S_*` (true absorption) and document that, with a
   scatter model active, `A_*` includes scatter while `A_BULK_*` does not.
   Without a scatter model active, `S_* = 0` and `A_BULK_* == A_*` (no change
   for every existing user — a hard requirement).

### Python
6. `ScatterMatrix.scatter()` convenience view (like `.absorption()`), and a
   `scatter_convention=` kwarg on `__init__`. Docstrings cross-link the
   `RoughnessType` caveat already added in the docs pass.

---

## 4. Validation

- **Closed-form single interface (the primary anchor).** Mirror the R1.1
  literature test: assert engine `R_spec`, `T_spec`, `S_R`, `S_T` match the
  closed forms in §1 to machine precision across contrast × sigma × angle ×
  pol, and assert `R+T+S = 1` exactly for lossless media (no overshoot, unlike
  NC).
- **TIS textbook limit.** At small `kz*sigma`, assert `S_R/R0` → `(4*pi*sigma*
  cos(theta)/lambda)^2` to first order.
- **Absorption separation.** Lossless rough stack: `A_BULK ≈ 0`, `S > 0`,
  `A_total = S`. Absorbing rough stack: `A_BULK > 0`, `S > 0`, `A_total =
  A_BULK + S`. This is the test that proves scatter is no longer mislabeled as
  absorption.
- **Backward compatibility.** With no DW interface present, every existing
  golden/parity/differential test is bit-identical (S channel all-zero,
  `A_BULK == A`). This is the gate: the change is purely additive for types 0–5.
- **Cross-model sanity.** At matched sigma/contrast, DW specular R < NC specular
  R (NC has no scatter loss), and `NC(R+T) ≈ 1` while `DW(R+T) = 1 - S < 1` — a
  regression-style comparison table committed alongside.
- **Parity mirror.** Add the DW factors to `validation/parity/smatrix/refs/
  loom_matrix.py` in lockstep (as R1.1 did for NC) so Python↔Rust parity stays
  meaningful.

New test files: `validation/smoke/test_scatter_loss.py` (closed-form + TIS
limit + absorption separation), and a cargo unit test for the per-interface
factor and the `R+T+S=1` invariant.

---

## 5. Risk, effort, sequencing

- **Risk M.** (a) The absorption-channel semantics change is the sharp edge —
  mitigated by the hard backward-compat gate (`S=0`, `A_BULK==A` when no DW
  interface is present) and the request-bit sync test. (b) Multi-interface
  scatter re-interception is deliberately out of scope for v1 (documented as
  single-scatter); revisit only if a use case needs it. (c) The `kz1^2` vs
  `kz1*kz2` reflection convention is a real physics choice — expose both, pick
  DW as default, document the XRR-vs-optical distinction.
- **Effort M–L.** The specular amplitude factors are S (one match arm per site).
  The scatter accounting channel + the `A_BULK` split + request-bit plumbing +
  Python surface + tests are the bulk (M–L).
- **Sequencing.** After R2.5 (request-bit sync test) so the new `SCATTER_*` /
  `A_BULK_*` bits are guarded from day one. Independent of the physics-fix
  phase; a natural Phase-4 feature item (call it R4.7). Ships in a minor version
  with a changelog entry documenting the `A_*` semantics note.

---

## 6. Explicitly out of scope (v1)

- Angle-resolved diffuse scatter / BRDF (only integrated TIS is reported).
- Re-interception of diffusely scattered light by other interfaces (single-
  scatter assumption).
- Correlated roughness between interfaces (each interface independent).
- Field-integration bulk absorption (the model-derived TIS split is used).
