# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""
Navette: Weaving the mathematics of light in thin film systems
Copyright (c) 2026 opticsWolf

Optimized Python Wrapper over `spectralweave` Rust extension for Targets.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Literal, Union

import json as _json

import numpy as np

# Import the ultra-fast Rust implementation
try:
  from navette._spectralweave import TargetWeaver
  from navette._spectralweave import OpticalWeaver
  from navette._spectralweave import calculate_merit as _calculate_merit
  from navette._smatrix import WeightProfile
except ImportError as exc:  # pragma: no cover - native not built yet
  raise ImportError(
    "Targets need the compiled `navette._spectralweave` extension. "
    "Build it with: maturin develop --release  # from the repo root"
  ) from exc

TargetType = Literal["a", "b", "e", "r", "c"]
#: Labels that are phase quantities by definition (raw values, phase=True).
_PHASE_ONLY = ("PDts", "PDtp", "RetR", "RetT", "DeltaR", "DeltaT",
               "GDR", "GDT", "GDDR", "GDDT", "TODR", "TODT", "FODR", "FODT")
NormalizationMode = Literal["auto", "linear", "log", "phase", "complex"]

# ---------------------------------------------------------------------------
# 1. Target data classes (user-facing inputs)
# ---------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class SpectralTarget:
    """One constraint curve: value vs wavelength at fixed angle.

    ``spectral`` + ``polarization`` name the curve. The vocabulary is::

        spectral:      R | T | A        reflectance / transmittance / absorptance
                       RB | TB | AB     the same, from the BACK side
                       PDts | PDtp      transmitted differential phase
                       RetR | RetT      retardance arg(p·conj(s)), r_f / t_f
                       DeltaR | DeltaT  ellipsometric Δ (RetR − π / RetT)
                       RatioR | RatioT  intensity ratio Rs/Rp, Ts/Tp
                       ExtR | ExtT      extinction ratio Rp/Rs, Tp/Ts
                       DiattR | DiattT  diattenuation (p − s)/(p + s)
                       PsiR | PsiT      ellipsometric Ψ = atan √(p/s)
                       GDR  | GDT       group delay dφ/dω (fs), r_f / t_f
                       GDDR | GDDT      group-delay dispersion (fs²)
                       TODR | TODT      third-order dispersion (fs³)
                       FODR | FODT      fourth-order dispersion (fs⁴)
                       S0R … S3R, DOPR  Stokes vector / DOP of r_f for an input
                       S0T … S3T, DOPT  the same, t_f
        polarization:  s | p | u        u = unpolarized (the s/p mean)

    giving ``Rs``…``Au``, ``RBs``…``ABu``, the two differential labels,
    the four retardance labels, the eight ratio labels, the eight
    dispersion labels and the ten Stokes labels.
    Any other pair raises at construction. ``A`` is derived as 1 − R − T on
    the shared grid, so a missing companion curve fails the whole key group.

    ``PDts``/``PDtp`` are the exception to the grid above: the label already
    encodes its polarization (``PDts`` is s, ``PDtp`` is p — a mismatch
    raises), it forces ``phase=True`` and raw-radian normalization, and it
    demands the coating-INDUCED phase, i.e. ``arg(t)`` minus the equivalent
    incidence-medium slab. **Transmission only.** The machinery carries a
    ``passes`` factor and ``reference_phase`` handles ``passes=2`` for a
    reflection round trip, but no reflection differential label exists, so
    there is no ``PDrs``/``PDrp`` to reach for. Absolute-phase targets on
    ``R``/``T`` with ``phase=True`` are unaffected and stay available.

    ``RetR``/``RetT``/``DeltaR``/``DeltaT`` are p-minus-s phase
    differences of one element, so they take ``polarization="u"`` (anything
    else raises) and, like ``PDts``, need ``phase=True`` and force
    raw-radian normalization. Their values are the engine's
    ``Retardance_R``/``Retardance_T`` and ``Delta_R``/``Delta_T``
    (``Delta_R = Retardance_R − π`` wrapped, ``Delta_T = Retardance_T``).
    At normal incidence s and p coincide and the retardance is identically
    0; a retardance demand at 0° only compiles with a warning.

    ``RatioR``/``RatioT``, ``ExtR``/``ExtT``, ``DiattR``/``DiattT`` and
    ``PsiR``/``PsiT`` are front-side functions of the s and p intensities:
    s/p, p/s, the engine's ``Diattenuation_R``/``Diattenuation_T`` and its
    ``Psi_R``/``Psi_T`` (radians, ``π/2`` below the engine's Rs/Ts floor).
    They take ``polarization="u"`` and ``phase=False`` (anything else
    raises), and each ratio's denominator carries the engine's 1e-20 guard.
    At normal incidence s and p are equal, the value is fixed (1, 1, 0 or
    π/4) and a demand at 0° only compiles with a warning.

    ``GDR``…``FODT`` are the engine's ``GD``/``GDD``/``TOD``/``FOD`` of
    the s or p amplitude (``polarization="s"`` or ``"p"``; ``"u"`` raises):
    the unwrapped ``arg`` of r_f / t_f differentiated 1–4 times in
    ω = 2πc/λ over the *simulated* wavelength grid, then sampled at the
    target points. Values are fs^n, need ``phase=True`` and are never
    normalized or wrapped. Like the engine's, they are coherent-stack
    quantities, and the simulated grid must be fine enough to resolve them
    (one wavelength gives zero).

    ``S0R``…``S3R``/``DOPR`` and ``S0T``…``S3T``/``DOPT`` are the Stokes
    vector and degree of polarization of the light r_f / t_f returns for
    the input state ``stokes_input``, a normalized Stokes vector
    ``(q1, q2, q3)`` with ``|q| <= 1`` (required for these labels, refused
    for any other). The axes are the engine's: ``(1, 0, 0)`` is p,
    ``(-1, 0, 0)`` s, and ``(0, 1, 0)`` the engine's own ``Es = Ep`` input,
    for which the values are half the engine's ``S0_R``…``S3_T`` (the input
    carries ``S0 = 1``) and the DOP is the engine's ``DOP_R``/``DOP_T``.
    Partially polarized input (``|q| < 1``) is what gives DOP its meaning:
    a coherent stack keeps a fully polarized input fully polarized. They
    take ``polarization="u"`` and ``phase=False`` and are normalized
    linearly (S1–S3 are signed). ``S2``/``S3`` read the retardance and are
    dropped across an incoherent film; ``S0``, ``S1`` and DOP read only
    the intensities and reach through one.
    See ``docs/spectralweave-target-kinds.md`` for the sign convention
    (the reference is ``+kD``, the conjugate of Macleod) and the needle
    ``gain_shift_phi`` bookkeeping.

    Kinds ``e``/``a``/``b`` are exact/above/below; ``r`` is a hard range box
    (zero merit inside ±``band``, quadratic exceedance outside) and ``c`` a
    soft center box (reduced ``(d/band)^2`` inside, exceedance + 1 outside).
    ``band`` is a per-point half-width in raw units (scalar broadcasts);
    omit it and ``r`` falls back to ±``tolerances``, ``c`` to ``e``.
    Tolerances are fractions of the curve level (linear/log), absolute
    radians (phase) — see ``docs/spectralweave-target-kinds.md``.
    ``weight`` multiplies the target's merit sum (default 1; relative
    importance across targets — the 200-pt vs 10-pt density question).
    ``weights`` (one per point) or ``weight_profile`` multiply ``weight``
    point by point: a residual is scaled by ``√(weight·w_i/count)``, and an
    integral target takes the weighted means ``Σw_i d_i/Σw_i`` over
    ``Σw_i tol_i/Σw_i``. Tolerances keep their meaning (scale, band);
    weights carry emphasis. A profile is evaluated on the target's axis:
    wavelength in nm here, angle in degrees on :class:`AngularTarget`.
    ``normalize_count=True`` additionally divides by the target-level point
    count (spectral grid size, angular angle count), turning the sum into a
    mean: equal say regardless of sampling density. Both flow verbatim into
    ``MeritSpec`` and the needle fold (missing-data penalties unaffected —
    drop a target to silence those, weight 0 only mutes present data).
    ``integral=True`` instead constrains the MEAN of the scaled diffs
    (single residual ``mean(d)/mean(tol)`` — kinds apply once to the mean,
    e.g. integral-``a`` is a lower bound on the average); it rejects
    ``normalize_count`` (the mean already is one). Regular and integral
    targets mix freely in one collection/run.
    Note for needle optimization: the folded needle evaluation drops the
    ``c`` +1 outside-level (constant offset — gradients exact, values read
    lower by the violated-point count); see
    ``docs/spectralweave-target-kinds.md``.
    """
    wavelengths:  np.ndarray
    values:       np.ndarray
    tolerances:   np.ndarray
    polarization: str
    spectral:     str
    kind:         TargetType = "e"
    normalization_mode: NormalizationMode = "auto"
    band:         Union[np.ndarray, float, None] = None
    phase:        bool = False
    weight:       float = 1.0
    normalize_count: bool = False
    integral:     bool = False
    environment:  Union[str, None] = None
    #: W4 point weights: one per point (``weights``) or a profile over
    #: the target's axis (``weight_profile``: a :class:`WeightProfile`, or
    #: the name of one in the collection's ``weight_profiles``). Give at
    #: most one.
    weights:      Union[np.ndarray, None] = None
    weight_profile: Union["WeightProfile", str, None] = None
    #: W8: the input polarization of a Stokes label (``S0R``…``DOPT``), as
    #: a normalized Stokes vector ``(q1, q2, q3)`` with ``|q| <= 1``.
    stokes_input: Union[tuple, None] = None
    #: The angle of incidence, in degrees or in radians: give exactly one
    #: (keyword-only). It is kept as given; the native layer converts.
    angle_deg:    Union[float, None] = field(default=None, kw_only=True)
    angle_rad:    Union[float, None] = field(default=None, kw_only=True)

    def _dump(self) -> dict:
        band = self.band
        if isinstance(band, (int, float)):
            band = float(band)
        elif band is not None:
            band = np.ascontiguousarray(np.asarray(band, dtype=np.float64)).ravel().tolist()
        return {
            **_dump_environment(self.environment),
            **_dump_weights(self.weights, self.weight_profile),
            **_dump_stokes_input(self.stokes_input),
            "wavelengths": np.ascontiguousarray(np.asarray(self.wavelengths, dtype=np.float64)).ravel().tolist(),
            "values": np.ascontiguousarray(np.asarray(self.values, dtype=np.float64)).ravel().tolist(),
            "tolerances": np.ascontiguousarray(np.asarray(self.tolerances, dtype=np.float64)).ravel().tolist(),
            **_dump_angle("angle", self.angle_deg, self.angle_rad),
            "polarization": str(self.polarization),
            "spectral": str(self.spectral), "kind": str(self.kind),
            "normalization_mode": str(self.normalization_mode),
            "band": band, "phase": bool(self.phase),
            "weight": float(self.weight) if isinstance(self.weight, (int, float)) else self.weight,
            "normalize_count": bool(self.normalize_count),
            "integral": bool(self.integral),
        }

    def __post_init__(self) -> None:
        # Shaping (broadcast/contiguity) stays; all CHECKS run natively,
        # once, via validate_targets (single validator, no duplication).
        from navette._structure import validate_targets as _validate
        import json as _json
        band = self.band
        if isinstance(band, (int, float)):
            band = np.ascontiguousarray(np.full(np.shape(self.values), float(band), dtype=np.float64))
        elif band is not None:
            band = np.ascontiguousarray(np.asarray(band, dtype=np.float64))
        # FFI Guard: Ensure contiguous float64 memory to prevent Rust segfaults.
        # This is an O(1) no-op if the array is already correctly formatted.
        object.__setattr__(self, 'wavelengths', np.ascontiguousarray(self.wavelengths, dtype=np.float64))
        object.__setattr__(self, 'values', np.ascontiguousarray(self.values, dtype=np.float64))
        object.__setattr__(self, 'tolerances', np.ascontiguousarray(self.tolerances, dtype=np.float64))
        object.__setattr__(self, 'band', band)
        _validate(_json.dumps({"spectral": [self._dump()], "angular": []}))


@dataclass(slots=True, frozen=True)
class AngularTarget:
    """One constraint curve: value vs angle at fixed wavelength.

    Same ``kind``/``band`` semantics as :class:`SpectralTarget`.
    """
    wavelength:   float
    values:       np.ndarray
    tolerances:   np.ndarray
    polarization: str
    spectral:     str
    kind:         TargetType = "e"
    normalization_mode: NormalizationMode = "auto"
    band:         Union[np.ndarray, float, None] = None
    phase:        bool = False
    weight:       float = 1.0
    normalize_count: bool = False
    integral:     bool = False
    environment:  Union[str, None] = None
    #: W4 point weights: one per point (``weights``) or a profile over
    #: the target's axis (``weight_profile``: a :class:`WeightProfile`, or
    #: the name of one in the collection's ``weight_profiles``). Give at
    #: most one.
    weights:      Union[np.ndarray, None] = None
    weight_profile: Union["WeightProfile", str, None] = None
    #: W8: the input polarization of a Stokes label, as on
    #: :class:`SpectralTarget`.
    stokes_input: Union[tuple, None] = None
    #: The angle grid, in degrees or in radians: give exactly one
    #: (keyword-only). It is kept as given; the native layer converts.
    angles_deg:   Union[np.ndarray, None] = field(default=None, kw_only=True)
    angles_rad:   Union[np.ndarray, None] = field(default=None, kw_only=True)

    def _dump(self) -> dict:
        band = self.band
        if isinstance(band, (int, float)):
            band = float(band)
        elif band is not None:
            band = np.ascontiguousarray(np.asarray(band, dtype=np.float64)).ravel().tolist()
        return {
            **_dump_environment(self.environment),
            **_dump_weights(self.weights, self.weight_profile),
            **_dump_stokes_input(self.stokes_input),
            "wavelength": float(self.wavelength),
            **_dump_angle("angles", self.angles_deg, self.angles_rad),
            "values": np.ascontiguousarray(np.asarray(self.values, dtype=np.float64)).ravel().tolist(),
            "tolerances": np.ascontiguousarray(np.asarray(self.tolerances, dtype=np.float64)).ravel().tolist(),
            "polarization": str(self.polarization),
            "spectral": str(self.spectral), "kind": str(self.kind),
            "normalization_mode": str(self.normalization_mode),
            "band": band, "phase": bool(self.phase),
            "weight": float(self.weight) if isinstance(self.weight, (int, float)) else self.weight,
            "normalize_count": bool(self.normalize_count),
            "integral": bool(self.integral),
        }

    def __post_init__(self) -> None:
        from navette._structure import validate_targets as _validate
        import json as _json
        band = self.band
        if isinstance(band, (int, float)):
            band = np.ascontiguousarray(np.full(np.shape(self.values), float(band), dtype=np.float64))
        elif band is not None:
            band = np.ascontiguousarray(np.asarray(band, dtype=np.float64))
        # FFI Guard
        for name in ('angles_deg', 'angles_rad'):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, np.ascontiguousarray(
                    np.asarray(getattr(self, name), dtype=np.float64)).ravel())
        object.__setattr__(self, 'values', np.ascontiguousarray(self.values, dtype=np.float64))
        object.__setattr__(self, 'tolerances', np.ascontiguousarray(self.tolerances, dtype=np.float64))
        object.__setattr__(self, 'band', band)
        _validate(_json.dumps({"spectral": [], "angular": [self._dump()]}))

ColorQuantity = Literal["Lab", "XyY", "LCh", "Oklab", "Y", "DomWl", "sRGB", "Luv", "XYZ",
                             "Din99", "White", "Yellow"]
ColorDistance = Literal["DeltaE2000", "DeltaE76", "Channels"]


@dataclass(slots=True, frozen=True)
class ColorTarget:
    """One colorimetric demand: front R/T intensity curve at fixed angle,
    integrated against an illuminant x observer pair.

    Defaults are D65 + 1931 2 deg; every field is overridable. ``reference``
    is a triple (e.g. Lab ``(62, 18, -34)``) — or a scalar for quantity
    ``"Y"`` (P2); mismatched shapes are refused natively. ``illuminant`` /
    ``observer`` are ``"D65"`` / ``"1931_2deg"`` or explicit table dicts
    (``{"wavelengths", "values"}`` / ``{"wavelengths", "xyz"}``);
    unknown names are refused natively (no registry to drift).
    Lab white point is always the demand illuminant's own white.
    Limits (all refused natively with named messages): front R/T curves
    only (no back, no absorption/phase curves), kind ``"Exact"`` only,
    ``"linear"`` only; quantities Lab|XyY|LCh|Oklab|Y|DomWl|sRGB|Luv|XYZ|
    Din99|White|Yellow — DeltaE2000|DeltaE76 in Lab|LCh only (everything
    else takes Channels; Y|Yellow are scalar Channels; DomWl|White are
    [2]-pair Channels: (nm, purity) | (whiteness, tint)). ``reference``
    is a triple, a scalar for ``"Y"``|``"Yellow"``, or a pair for
    ``"DomWl"``|``"White"``. Channels tolerances are unit and fixed
    (no per-channel tol in v1) — see
    ``docs/spectralweave-target-kinds.md`` for the effective weighting
    this implies per quantity (notably DomWl purity).
    Oklab references are D65-Oklab (non-D65 demands Bradford-adapt first;
    in-tree D65 white carries a pre-existing ~1e-4 b offset — systematic,
    far below demand relevance). LCh hue wraps before scaling; near-neutral
    hues are gradient-unstable by nature (weight hue tolerance generously
    near the achromatic axis). No ``integral`` /
    ``count_norm`` / ``phase`` / ``band`` — a color demand already is
    integral. ``weight`` scales the demand's merit sum (default 1).
    ``wavelength_range`` is an optional ``[lo, hi]`` nm window restricting
    the *sample* integral only — the white stays the full-illuminant white
    and the DomWl locus stays full-CMF, so windowed numbers stay comparable
    to full-range numbers (``None`` = full overlap).
    Needle note: ``Ru``/``Tu`` demands split half per polarization branch
    (Ru is the Rs/Rp mean); s/p demands ride their shared bucket, same
    convention as pointwise targets — see
    ``docs/plans/color_targets_optionB_plan.md`` D5.
    All checks run natively, once, via validate_targets (single validator,
    no duplication) — shaping here is normalization only.
    """
    curve:        str = "Ru"
    illuminant:   Union[str, dict] = "D65"
    observer:     Union[str, dict] = "1931_2deg"
    quantity:     str = "Lab"
    reference:    Union[tuple, list, float] = (62.0, 18.0, -34.0)
    distance:     str = "DeltaE2000"
    weight:       float = 1.0
    kind:         str = "Exact"
    transform:    str = "linear"
    wavelength_range: Union[tuple, list, None] = None
    yi_cx:        Union[float, None] = None
    yi_cz:        Union[float, None] = None
    environment:  Union[str, None] = None
    #: The angle of incidence, in degrees or in radians: give at most one
    #: (keyword-only); neither is normal incidence. It is kept as given;
    #: the native layer converts.
    angle_deg:    Union[float, None] = field(default=None, kw_only=True)
    angle_rad:    Union[float, None] = field(default=None, kw_only=True)

    def _dump(self) -> dict:
        return {
            **_dump_environment(self.environment),
            "curve": str(self.curve),
            **_dump_angle("angle", self.angle_deg, self.angle_rad, normal=True),
            "illuminant": _dump_table(self.illuminant, "illuminant"),
            "observer": _dump_table(self.observer, "observer"),
            "quantity": str(self.quantity),
            "reference": _dump_reference(self.reference),
            "distance": str(self.distance),
            "weight": float(self.weight) if isinstance(self.weight, (int, float)) else self.weight,
            "kind": str(self.kind), "transform": str(self.transform),
            "wavelength_range": _dump_range(self.wavelength_range),
            "yi_cx": self.yi_cx, "yi_cz": self.yi_cz,
        }

    def __post_init__(self) -> None:
        from navette._structure import validate_targets as _validate
        import json as _json
        _validate(_json.dumps({"spectral": [], "angular": [], "color": [self._dump()]}))


def _dump_angle(stem: str, deg, rad, normal: bool = False) -> dict:
    """The JSON angle key(s) as given: ``"<stem>_deg"`` and/or
    ``"<stem>_rad"``. The native validator refuses both or neither; with
    ``normal`` (an optional angle), neither is normal incidence."""
    out = {}
    for unit, value in (("deg", deg), ("rad", rad)):
        if value is not None:
            out[f"{stem}_{unit}"] = (float(value) if np.ndim(value) == 0 else
                                     np.asarray(value, dtype=np.float64).ravel().tolist())
    if normal and not out:
        out[f"{stem}_rad"] = 0.0
    return out


def _dump_environment(env):
    """F2.4: the `environment` tag, emitted only when it is set.

    Absent means "the first environment", which is what keeps every
    pre-F2.1 target set meaningful without an edit — and emitting the key
    as ``null`` would change the JSON of documents that never asked for
    environments, so it is left out entirely.

    The NAME is not checked here. ``__post_init__`` validates one target
    against no roster at all; the roster is a property of the run, and
    ``build_merit_spec`` is where both are in the same room. An unknown
    name refuses there, naming the demand and the known environments.
    """
    if env is None:
        return {}
    if not isinstance(env, str) or not env:
        raise ValueError("environment must be a non-empty name, or None.")
    return {"environment": env}


def _dump_stokes_input(q) -> dict:
    """W8: ``stokes_input``, emitted only when set (every other label's
    document is unchanged). Its checks (three finite numbers, ``|q| <= 1``,
    a Stokes label) run at compile."""
    if q is None:
        return {}
    return {"stokes_input": [float(x) for x in q]}


def _dump_weights(weights, profile) -> dict:
    """W4: ``weights`` / ``weight_profile``, emitted only when set (so a
    target without them dumps the JSON it always did). A
    :class:`WeightProfile` travels as its JSON document, a name as the
    string the compile resolves against ``weight_profiles``. No checks
    here: the native validator owns them."""
    out = {}
    if weights is not None:
        out["weights"] = np.ascontiguousarray(
            np.asarray(weights, dtype=np.float64)).ravel().tolist()
    if profile is not None:
        out["weight_profile"] = (profile if isinstance(profile, str)
                                 else _json.loads(profile.to_json()))
    return out


def _dump_table(spec: Union[str, dict], role: str):
    """Name passthrough (native resolves ``D65``/``1931_2deg`` from the
    embedded defaults) or explicit table dict with array normalization.
    No rule checks here — the native validator owns them."""
    if isinstance(spec, str):
        return spec
    if isinstance(spec, dict):
        out = {}
        for key, value in spec.items():
            if isinstance(value, np.ndarray):
                value = np.ascontiguousarray(np.asarray(value, dtype=np.float64)).ravel().tolist()
            out[key] = value
        return out
    return spec


def _dump_range(wr: Union[tuple, list, None]):
    """[lo, hi] -> [f, f], None -> None; anything else passes through
    for the native validator to refuse with a named message."""
    if wr is None:
        return None
    if isinstance(wr, (tuple, list)):
        return [float(v) for v in wr]
    return wr


def _dump_reference(ref: Union[tuple, list, float]):
    """Triple -> [f, f, f], scalar -> float; anything else passes through
    for the native validator to refuse with a named message."""
    if isinstance(ref, bool):
        return ref
    if isinstance(ref, (int, float)):
        return float(ref)
    if isinstance(ref, (tuple, list)):
        return [float(v) for v in ref]
    return ref


BaseTarget = Union[SpectralTarget, AngularTarget, ColorTarget]

# ---------------------------------------------------------------------------
# 2. TargetCollection (lightweight, standalone)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class TargetCollection:
    """
    User-facing container for mixed Spectral, Angular, and Color targets.
    """
    _spectral_targets: list[SpectralTarget] = field(default_factory=list)
    _angular_targets:  list[AngularTarget]  = field(default_factory=list)
    _color_targets:    list[ColorTarget]    = field(default_factory=list)
    #: W4: named profiles a target's ``weight_profile`` may name.
    weight_profiles:   dict[str, WeightProfile] = field(default_factory=dict)

    def add(self, target: BaseTarget) -> None:
        """Append a :class:`SpectralTarget`, :class:`AngularTarget`, or
        :class:`ColorTarget`."""
        if isinstance(target, SpectralTarget):
            self._spectral_targets.append(target)
        elif isinstance(target, AngularTarget):
            self._angular_targets.append(target)
        elif isinstance(target, ColorTarget):
            self._color_targets.append(target)
        else:
            raise TypeError(f"Unsupported target type: {type(target)}")

    def clear(self) -> None:
        """Remove all spectral, angular, and color targets."""
        self._spectral_targets.clear()
        self._angular_targets.clear()
        self._color_targets.clear()

    @property
    def spectral_targets(self) -> list[SpectralTarget]:
        """The ingested wavelength-domain targets."""
        return self._spectral_targets

    @property
    def angular_targets(self) -> list[AngularTarget]:
        """The ingested angle-domain targets."""
        return self._angular_targets

    @property
    def color_targets(self) -> list[ColorTarget]:
        """The ingested colorimetric demands."""
        return self._color_targets

    @property
    def count(self) -> int:
        """Total number of spectral plus angular plus color targets."""
        return len(self._spectral_targets) + len(self._angular_targets) + len(self._color_targets)

    def build_weaver(self, cache_size: int = 128, tolerance_floor: float = 1e-12) -> TargetWeaver:
        """
        Compiles the defined targets into a Rust-native TargetWeaver.
        """
        for t in (*self._spectral_targets, *self._angular_targets):
            if t.weights is not None or t.weight_profile is not None:
                raise ValueError(
                    "build_weaver: point weights (weights / weight_profile) are "
                    "read by the synthesis merit (build_merit_spec) only; "
                    "calculate_merit would ignore them. Drop them for this path.")
        weaver = TargetWeaver(cache_size=cache_size, tolerance_floor=tolerance_floor)

        # Iterating in Python and passing to Rust here is fine because it
        # only happens ONCE during initialization, not inside the hot loop.
        # PDts/PDtp and the retardance labels force phase normalization
        # (raw radians, nf == 1): they are phase quantities and any other
        # resolution would not unscale back to radians in the synthesis
        # converter.
        for t in self._spectral_targets:
            weaver.add_spectral_target(
                t.wavelengths, t.values, t.tolerances,
                t.polarization, t.spectral, t.kind,
                "phase" if t.spectral in _PHASE_ONLY else t.normalization_mode,
                t.band,
                weight=t.weight, normalize_count=t.normalize_count,
                integral=t.integral,
                angle_deg=t.angle_deg, angle_rad=t.angle_rad,
            )

        for t in self._angular_targets:
            weaver.add_angular_target(
                t.wavelength, t.values, t.tolerances,
                t.polarization, t.spectral, t.kind,
                "phase" if t.spectral in _PHASE_ONLY else t.normalization_mode,
                t.band,
                weight=t.weight, normalize_count=t.normalize_count,
                integral=t.integral,
                angles_deg=t.angles_deg, angles_rad=t.angles_rad,
            )

        return weaver


# ---------------------------------------------------------------------------
# 3. Merit Function
# ---------------------------------------------------------------------------
def calculate_merit(
    sim_weaver: OpticalWeaver,  # Will accept the unwrapped Rust object
    target_weaver: TargetWeaver,
    *,
    missing_penalty: float = 1e6,
) -> float:
    """
    Delegates the Merit Function computation completely to compiled Rust,
    bypassing the Python GIL entirely for true concurrency.
    """
    return _calculate_merit(sim_weaver, target_weaver, missing_penalty)