# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Shared dtypes, enums and solver-array containers for layer stacks.

All solver-facing arrays use fixed dtypes (``FLOAT_TYPE``/``COMPLEX_TYPE``/
``INT_TYPE``) so the Python stack model matches the native engine layouts
bit-for-bit. The ``*Mask`` enums index the per-layer flag/error vectors.
"""
from dataclasses import dataclass
from enum import IntEnum
from typing import NamedTuple
import numpy as np

# Standardised numeric types (match the native engine dtypes)
FLOAT_TYPE = np.float64
COMPLEX_TYPE = np.complex128
INT_TYPE = np.int32

class ErrorType(IntEnum):
    """Statistical law used when drawing fabrication errors."""
    GAUSSIAN = 0
    UNIFORM = 1
    COMBINED = 2

class RoughnessType(IntEnum):
    """Per-interface roughness form factor (solver contract, [nm] sigma).

    Canonical definition shared by the structure model and the smatrix
    engine (`navette.smatrix.RoughnessType` re-exports this): NONE is an
    ideal interface; LINEAR/STEP/EXPONENTIAL/GAUSSIAN are analytic
    graded-index profiles; NEVOT_CROCE is the Nevot-Croce X-ray factor.
    Stored on :class:`Layer.rough_type` and passed to the engine as int.

    Scatter-loss caveat (important)
    -------------------------------
    **All of these models are specular-beam (coherent) only. None of them
    account for diffuse scatter as a tracked quantity.** A topographically
    rough interface scatters energy out of the specular beam into the diffuse
    hemisphere (Total Integrated Scatter, TIS ~= (4*pi*sigma*cos(theta)/lambda)^2);
    that energy is not represented in R, T, A, or any other channel here.

    Two physically distinct pictures share this enum:

    * **Graded transition (LINEAR/STEP/EXPONENTIAL/GAUSSIAN, and NEVOT_CROCE).**
      Roughness is treated as a laterally-uniform index grading, so energy is
      only redistributed between the specular beams -- no diffuse scatter.
      NEVOT_CROCE is constructed to conserve specular energy: reflection is
      damped by exp(-2*kz1*kz2*sigma^2) and transmission is *enhanced* by
      exp(+((kz1-kz2)*sigma)^2/2) so that R_spec + T_spec = 1 to first order in
      sigma^2.

      It is only *perturbatively* unitary, and the cancellation is first-order
      only: it holds for |kz*sigma| << 1 and degrades outside that regime.
      Typical optical coatings (sigma <~ 5 nm, visible) overshoot by ~1e-4;
      but the transmission factor grows without bound, so at high contrast and
      large sigma the model produces flatly unphysical output. Measured, single
      lossless interface, normal incidence, 550 nm:

        n 1 -> 2.35, sigma = 20 nm:  R+T = 1.0206,  A = 1-R-T = -2.1e-2
        n 1 -> 4.28, sigma = 20 nm:  T   = 1.0768,  A = 1-R-T = -2.3e-1

      The controlling quantity is the index *contrast*, not sigma alone: the
      exponent is ((kz1-kz2)*sigma)^2/2, i.e. (2*pi*dn*sigma/lambda)^2/2 at
      normal incidence. That is why a model borrowed from X-ray reflectometry
      (dn ~ 1e-5, where the factor is 1.000000...) misbehaves at optical
      contrast. Usable budget -- sigma injecting at most a fraction eps of
      spurious energy at one interface:

        sigma_max = sqrt(ln(1+eps)) * lambda / (2*pi*dn)

      For a 1% budget (eps = 0.01, so sigma_max ~= 0.0159*lambda/dn):

        dn  = 1.35   lambda = 400 nm ->  4.7 nm  |  550 nm ->  6.5 nm
        dn  = 0.6    lambda = 400 nm -> 10.6 nm  |  550 nm -> 14.6 nm
        dn  = 0.3    lambda = 400 nm -> 21.2 nm  |  550 nm -> 29.1 nm

      Note the consequences: **T can exceed 1 and the residual absorptance
      A = 1 - R - T can go negative.** Neither is clamped. Note the direction
      too: real roughness scatters energy *out* of the specular beam, so the
      physically correct result is R+T slightly *below* 1. R+T > 1 has no
      physical mechanism behind it and is the unambiguous sign that the model
      has been pushed past its range. Check
      :meth:`~navette.smatrix.ScatterMatrix.energy_conservation` if you are
      near the validity edge -- but note it reports |1-R-T|, so an energy
      *gain* is indistinguishable from absorption by magnitude alone.
      The graded types (1-4) instead
      damp both coefficients, so R_spec + T_spec < 1, but the deficit is a
      form-factor artifact, not a derived scatter loss.

    * **True topographic roughness with diffuse loss** is NOT modeled by any
      type here. If you need specular attenuation from real scatter (R+T<1 with
      the deficit equal to TIS), a Debye-Waller / TIS option is required
      (planned; see docs/plans/scatter_loss_plan.md).
    """
    NONE = 0
    LINEAR = 1
    STEP = 2
    EXPONENTIAL = 3
    GAUSSIAN = 4
    NEVOT_CROCE = 5

class ErrorMask(IntEnum):
    """Slots of the per-layer error vector (thickness, n/k, roughness, ...)."""
    THICKNESS = 0
    N_REAL = 1
    N_IMAG = 2
    ROUGHNESS = 3
    INH_DELTA = 4
    INTERFACE = 5

# State-schema version (structure get_state dicts + config states).
# v1 = the current baseline. There is no past: untagged states are
# malformed, not legacy. Bump on any breaking change (removed/renamed
# keys, changed meaning of an existing key); purely additive keys are
# safe without a bump (readers ignore unknown keys) — the fingerprint
# test in validation/regression/structure/test_roundtrip.py enforces
# this decision on every key-set change. Readers refuse anything but
# the current version (no silent misreads).
SCHEMA_VERSION = 1


def check_schema_version(state, what: str) -> None:
    """Refuse states not written at the current schema version.

    Thin over the native gate (single canonical home)."""
    from navette._structure import check_schema_version as _check
    _check(dict(state), what)


WARNING_PREFIX = "warning: "


def is_warning(issue: str) -> bool:
    """True for advisory validation issues (reported, never solve-blocking)."""
    return isinstance(issue, str) and issue.startswith(WARNING_PREFIX)


class LayerMask(IntEnum):
    """Slots of the per-layer status mask produced by :meth:`Layer.mask`."""
    ACTIVE = 0
    COHERENT = 1
    INHOMOGEN = 2
    ROUGHNESS = 3

class OptMask(IntEnum):
    """Slots of the per-group optimization mask (`Group.optimization_mask`).

    Convention: 1 = the property may be optimized, 0 = fixed. The mask
    refines `Layer.optimize` per property for optimizers that support it;
    `get_optimization_parameters` already honors the THICKNESS slot.
    MATERIAL governs material-substitution moves (e.g. needle insertion).
    """
    THICKNESS = 0
    N = 1
    K = 2
    ROUGHNESS = 3
    INH_DELTA = 4
    INTERFACE = 5
    MATERIAL = 6

class LayerType(IntEnum):
    """Design role of a layer (vocabulary for the former open `layer_typ` int).

    Markers delimit stacks: a `STACK` block opens with `AMBIENT` and closes
    with `SUBSTRATE`; `FILM` rows are the thin-film sequence. `FILM = 1`
    keeps every pre-enum state file valid (1 was the only value in use).
    """
    AMBIENT = 0
    FILM = 1
    SUBSTRATE = 2

class BlockKind(IntEnum):
    """Composition role of an architect block (declared, never inferred).

    `STACK` spans half-space to half-space; `FILMS` is a thin-film-only
    sequence legal only between stacks (or films), taking its boundary
    media from its neighbors at expansion time.
    """
    STACK = 0
    FILMS = 1

@dataclass(frozen=True)
class InterpolationSettings:
    method: str = "linear"
    floater_hormann_d: int = 3
    robust: bool = False

class SolverArrays(NamedTuple):
    indices: np.ndarray          # complex128, shape (n_total, n_wavs)
    thicknesses: np.ndarray      # float64,    shape (n_total,)
    incoherent_flags: np.ndarray # bool,       shape (n_total,)
    rough_types: np.ndarray      # int32,      shape (n_total,)
    rough_vals: np.ndarray       # float64,    shape (n_total,)