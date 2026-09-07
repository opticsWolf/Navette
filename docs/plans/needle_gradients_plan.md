# Needle + gradient/inhomogeneous layers — implementation plan

## N0. Verdict on the hypothesis (up front)

"Only thickness-dependent is stable" is **correct for insertion, with one
refinement — and thickness-optimization is stable for all four modes.**
The mechanism is *locality*: needle theory is a local sensitivity
analysis (merit change ≈ P(z*)·s for an infinitesimal seed at z*). The
finite step is only justified if the edit stays local. Consequences:

| profile | bisection | thickness-LM | verdict |
|---|---|---|---|
| Gradient RateCapped (b) | EXACT (absolute-z profile, §N2) | smooth (one kink at cap knee) | needle-stable, natural |
| Gradient FixedSpan (a) | EXACT only via endpoint subdivision (§N2); naive re-normalization is maximally non-local | smooth | needle-stable WITH the subdivision rule + defrag |
| Inhomogeneous either mode | NO exact bisection (symmetric ±δ can't represent an asymmetric sub-range) | smooth (same cap-knee note for (b)) | boundary-only insertion v1 |
| Uniform (today) | exact (identical halves) | smooth | unchanged |

So the refined rule: **insertion must be profile-preserving; thickness is
always a safe free parameter; inhomogeneous interiors are never split.**

## N1. Current-state inventory (addresses)

- `synthesis/needle_pass.rs`: `ScanSite { film_idx, depth_into_layer_nm,
  z_nm }`, `build_scan_sites` (interior multiples of `scan_step_nm` in
  every `needle=true` film — Python-parity grid), `NeedleTargets` fold,
  `run_needle_pass` (analytic sweep → most-negative-P site).
- `synthesis/cycle.rs`: one pass = repeated insertions (`max_needles`,
  `needle_seed_thickness_nm` = 5 nm default); convergence test technique:
  `MF_now − MF(1 nm test needle) ≈ −P(z*)·1 nm`.
- `synthesis/driver.rs`: `ArraySeed` (host → fresh seed carrier);
  seeds are catalog uniform slabs.
- Prerequisite: gradient plan §D4 (span provenance — `spans` on
  `DesignStack`). **N1 starts after G3 lands.** Without span identity,
  insertion can't distinguish "split a uniform row" from "subdivide a
  profile" — the whole plan keys on that distinction.

## N2. Profile-preserving bisection (the core specification)

Splitting a span at fraction α (position z*) MUST leave the realized
physical profile pointwise identical; only the bookkeeping changes:

- **RateCapped gradient:** left keeps `(f_start, rate, ref, caps)`;
  right gets `f_start' = f(z*)` with z measured from the split, same
  rate/ref/caps. Exact by construction (absolute-z profile).
- **FixedSpan gradient:** left `FixedSpan(f0, fα)`, right
  `FixedSpan(fα, f1)`, `fα = f0 + α(f1−f0)`. Exact (endpoint
  subdivision). NAIVE split (two re-normalized `f0→f1` halves) is
  FORBIDDEN — it fabricates a `f1|f0` discontinuity at z* and the P(z*)
  prediction becomes worthless (this is the instability the hypothesis
  senses; specified here so it can't be implemented by accident).
- **Inhomogeneous (both modes):** no split — interiors host no scan
  sites (§N3). Rationale: symmetric `±δ` cannot represent the
  asymmetric sub-ranges. (Future: an asymmetric-δ form would enable
  exact bisection — same deferred bucket as curved profiles.)
- **Uniform rows:** unchanged (identical halves; defrag re-merges).

Repeated subdivision fragments spans → post-pass **defrag sweep**:
adjacent gradient spans re-join iff subdivision-compatible (same
materials/EMA/mode, continuous endpoints within 1e-12, same rate/ref/
caps). Defrag is tested round-trip: subdivide-then-defrag ≡ original
spec field-for-field.

## N3. Scan sites inside spans (implementation in `needle_pass.rs`)

1. `build_scan_sites` becomes span-aware (takes spans + rows):
   - Span boundaries: always sites (as today — every interface).
   - Gradient-span interiors: sites at SUBLAYER boundaries of the
     discretized span (finite, consistent with the model the P-operator
     actually sees — not the blind `k·step` grid, which could land
     mid-slab and split uniformity pointlessly).
   - Inhomogeneous-span interiors: NO sites (boundary-only).
   - Uniform needle films: unchanged `k·step` grid.
2. Sensitivity at a span-interior site uses the LOCAL row's nk as the
   ambient for the infinitesimal insertion (the `p_coherent_*` machinery
   already operates per-row — no operator change, just correct indexing;
   same pattern as the Option-B color bucket: new call pattern, no new
   kernel).
3. Insertion at an interior site bisects per §N2 (subdivision), then
   inserts the catalog-material seed as a new uniform span between the
   halves. nk-keyed merge can't fuse seed↔span (nk differs — verified,
   not assumed); defrag must not fuse either (compatibility gate §N2).

## N4. Thickness-LM over spans (no new optimizer code)

- Span thickness is ONE free parameter; sublayer thicknesses scale
  proportionally (Σd ≡ thickness preserved exactly, remainder rule §D3
  of the gradient plan). Profile parameters (`f_start/f_end`, `rate`,
  `δ`) are NOT LM-free in N1 (design choices, not variables).
- Smoothness: FixedSpan/inhomogeneous-Fixed fully smooth (geometry
  scales, nk fixed). RateCapped predicts a kink where the cap
  binds/unbinds mid-optimization — piecewise-smooth, LM-tolerant (same
  class as box-constraint kinks everywhere else). Monitor via the N1
  convergence smoke; no special handling v1.
- **Follow-up (N2, separate patch): profile-LM** — `f_start/f_end` or
  `rate` as free parameters with FD derivatives (2 params, cheap; the
  color plan's FD-twin discipline applies). **Far follow-up (N3):
  gradient seeds** (inserting a profiled seed = rugate synthesis) —
  needs seed-profile policy + profile-LM first. Neither is designed
  here beyond naming the order.

## N5. Verification (twins + the stability acceptance test)

- **N-bisect identity:** bisect-then-`simulate()` ≡ pre-bisect spectra
  BITWISE (the profile-preservation proof; all three bisectable types).
- **N-midpoint continuity:** `f` continuous across a FixedSpan split to
  1e-15; RateCapped right-half `f_start'` equals hand-computed `f(z*)`.
- **N-site parity:** brute-force 1 nm test seed at a span-interior site
  reproduces `−P(z*)·1 nm` (the existing cycle.rs technique, extended
  into spans — proves the local-ambient indexing §N3.2).
- **N-absence:** inhomogeneous interiors emit zero sites; needle-insert
  into an inhomogeneous span refuses naming material+range (same
  message family as the G3 span guards).
- **N-defrag round-trip:** subdivide-then-defrag ≡ original spec.
- **N-stability (acceptance):** needle cycle on a design with one (a)-
  and one (b)-gradient span converges — MF non-increasing across
  passes, no site oscillation (same site re-picked with opposite sign
  = failure). This is the test that would catch naive re-normalization.

## N6. Progression (after G3 — span provenance is load-bearing)

- **N1 (0.4.31+, first patch after G3/G4):** §N2 bisection + §N3 sites +
  defrag + N-bisect/continuity/parity/absence/defrag twins. Boundary
  insertion into gradient spans works from day one of N1 (it's just
  span-aware bookkeeping); interior sites unlock in the same patch.
- **N2:** profile-LM (FD over profile params + twins). **N3:**
  gradient seeds (design when N2 lands).
- Color interplay: none — folded targets see rows; gradient rows are
  rows. No joint testing beyond the standard validation green.

## N7. Host effectiveness: gradients vs inhomogeneous (locked 2026-09-06)

Gradient spans are the better needle hosts — three reasons:
1. **More candidate sites:** N−1 interior sites + boundaries vs
   boundaries only. Needle is search over positions; finer search wins.
2. **More optimizable freedom:** mixture endpoints today, profile
   parameters (`f_start`, `rate`) as LM variables in N2; inhomogeneous
   carries one frozen δ.
3. **Moves mean something physical:** gradients are realizable
   deposition programs (EMA of two materials); hard optimization around
   a phenomenological scaling risks exploiting the fudge.
Caveats (also locked): gradient spans cost more rows per solve +
 fragmentation pressure (bench pairing quantifies); and the two model
 different physics — single-material drift (oxidation through a nitride)
 cannot be a mixture without a fake second material, so inhomogeneous
 keeps its domain. **Modeling guidance: where both could describe the
 same transition, prefer the gradient; reserve inhomogeneous for what it
 uniquely describes.**

## N8. Docstring contract (enforced at implementation, not review)

The following remarks MUST appear in the listed docstrings when N1/G1
 land (plan text alone is insufficient — the call site is where the
 choice bites). Review checklist for the implementing patch:
- `GradientSpec`: A-host rule + effect (`f = 0.1` ≡ 10% B in 90% A;
  A↔B swap is different physics) + deferred host-selection follow-up
  (gradient plan §D9.4 — repeated here so neither patch drops it).
- `GradientSpec::shape`: Linear-only v1, curved variants later.
- `InhMode`: Fixed = bit-identical legacy; RateCapped formula + cap
  binds nominal pre-error-draw; interiors never bisected (boundary-only).
- `build_scan_sites` / insertion path: site rules per span kind (§N3)
  + naive FixedSpan re-normalization forbidden (§N2) + §N7 modeling
  guidance (prefer gradients where both fit).
- Insertion refusal messages name material + range (G3 family).

## N9. Risks

| Risk | Mitigation |
|---|---|
| Naive FixedSpan split implemented by accident | §N2 forbids it normatively + N-bisect twin fails bitwise on any violation |
| Span fragmentation explodes row count | defrag sweep + `sublayers` cap (256); bench row counts in N1 |
| Cap-knee kink stalls LM | convergence smoke watches for it; kink is mild (flat-tail onset) |
| Site grid vs sublayer boundaries mismatch | sites ARE sublayer boundaries by construction — no second grid |
