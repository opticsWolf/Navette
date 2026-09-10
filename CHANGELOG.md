# Changelog

All notable changes to Navette are recorded here. Work items reference
`docs/remediation_plan.md` (Rx.y) and `docs/code_review.md` (§).

## [0.5.2] — energy_conservation() accepts 2-D input (R1.2, §15)

### Fixed

- **`ScatterMatrix.energy_conservation()` no longer raises `TypeError`.** The
  wrapper feeds 2-D `[n_angles, n_wavs]` arrays to the native
  `solver_energy_conservation`, which previously accepted only 1-D input, so
  the advertised API failed on every call. The binding now accepts 2-D, checks
  shape agreement, and reshapes the elementwise residual back to 2-D; the
  single-angle path still squeezes to 1-D.

### Added

- `validation/smoke/test_energy_conservation.py`: lossless (~0 residual),
  absorbing (residual == absorptance, monotone in thickness), 2-D shape, and
  1-D squeeze coverage.

### Changed

- `test_package_version` validates semver format instead of pinning a literal,
  so it survives per-item version bumps (release workflow still enforces sync).

## [0.5.1] — Névot-Croce transmission factor (R1.1, §3.2)

### Fixed

- **Energy conservation in the Névot-Croce roughness path (roughness type 5).**
  The type-5 branch previously applied the reflection Debye-Waller factor
  `f = exp(-2·kz1·kz2·σ²)` to the transmission amplitudes as well, silently
  destroying transmitted energy at rough interfaces. Transmission now takes the
  canonical wavevector-difference factor `exp(+((kz1-kz2)·σ)²/2)` (Névot & Croce
  1980; de Boer, Phys. Rev. B 44, 498 (1991); Stearns, J. Appl. Phys. 65, 491
  (1989)), so `R + T → 1` for lossless interfaces.

  Measured single interface, n = 1.0 → 1.5001, λ = 550 nm, normal incidence
  (no absorption — energy must be conserved):

  | σ (nm) | old R+T (buggy) | new R+T (fixed) |
  |--------|-----------------|-----------------|
  | 0      | 1.000000        | 1.000000        |
  | 3      | 0.992977        | 1.000001        |
  | 6      | 0.972202        | 1.000016        |
  | 10     | 0.924678        | 1.000125        |

  **This changes numerical results for every stack using roughness type 5**
  (XRR is the main consumer). Fixed at both call sites in
  `rust/navette/src/smatrix/coherent_block.rs`; the in-repo loom parity mirror
  was patched in lockstep.

### Added

- `validation/smoke/test_rtype5_energy.py`: energy-conservation sweep
  (rtype × σ × contrast × angle), a regression guard against re-damping
  transmission, and an independent closed-form Névot-Croce literature anchor
  checking engine R and T individually to machine precision.
