# -*- coding: utf-8 -*-
"""Construction-time validation of ``ScatterMatrix`` inputs (R3.1).

The engine is permissive on purpose -- it is also the optimizer's inner loop.
Before R3.1 that permissiveness reached the user unfiltered, and
``validation/review/garbage_in.py`` recorded what came back:

  * a negative thickness and a NaN thickness both produced *the same numbers
    as deleting the layer* -- a plausible answer to a question nobody asked;
  * 120 deg aliased onto 60 deg, because only ``sin(theta)`` reaches the
    engine;
  * a duplicated wavelength made every dispersion channel NaN (the GD/GDD
    kernels divide by the grid spacing);
  * a NaN index turned the whole output NaN with nothing naming the layer.

None raised. None warned. These tests hold the line on both sides of it: the
garbage must raise, with the offending index and value in the message, and
the legitimately odd inputs must still go through untouched.
"""

import numpy as np
import pytest

from navette.smatrix.smatrix import Request, ScatterMatrix

_N = np.array([1.0 + 0j, 2.35 + 0j, 1.46 + 0j, 1.52 + 0j])
_D = np.array([0.0, 120.0, 200.0, 0.0])
_WLS = np.linspace(450.0, 750.0, 31)
_ANG = [30.0]


def _sm(**kw):
    args = dict(layer_indices=_N, thicknesses=_D, wavelengths=_WLS, angles=_ANG)
    args.update(kw)
    layer_indices = args.pop("layer_indices")
    thicknesses = args.pop("thicknesses")
    return ScatterMatrix(layer_indices, thicknesses, **args)


def _with(arr, idx, value):
    out = np.array(arr, copy=True)
    out[idx] = value
    return out


# --------------------------------------------------------------------------
# Rejected inputs -- one row per silent-wrong-answer path
# --------------------------------------------------------------------------

_REJECT = [
    # (id, kwargs, fragments the message must contain)
    ("nan_index", dict(layer_indices=_with(_N, 1, np.nan)), ["layer_indices", "layer 1", "nan"]),
    ("inf_index", dict(layer_indices=_with(_N, 2, np.inf)), ["layer_indices", "layer 2", "inf"]),
    ("nan_imag_index", dict(layer_indices=_with(_N, 1, 2.0 + np.nan * 1j)), ["layer_indices", "layer 1"]),
    ("overflow_index", dict(layer_indices=_with(_N, 1, 1e308 + 0j)), ["layer_indices", "layer 1", "magnitude"]),
    ("nan_thickness", dict(thicknesses=_with(_D, 2, np.nan)), ["thicknesses", "layer 2", "nan"]),
    ("inf_thickness", dict(thicknesses=_with(_D, 2, np.inf)), ["thicknesses", "layer 2", "inf"]),
    ("negative_thickness", dict(thicknesses=_with(_D, 2, -50.0)), ["thicknesses", "layer 2", "-50.0"]),
    ("nan_wavelength", dict(wavelengths=_with(_WLS, 5, np.nan)), ["wavelengths", "index 5", "nan"]),
    ("zero_wavelength", dict(wavelengths=_with(_WLS, 0, 0.0)), ["wavelengths", "index 0", "0.0"]),
    ("negative_wavelength", dict(wavelengths=_with(_WLS, 3, -1.0)), ["wavelengths", "index 3", "-1.0"]),
    ("duplicate_wavelength", dict(wavelengths=np.array([500.0, 550.0, 550.0, 600.0])),
     ["wavelengths", "strictly increasing", "index 2", "550.0"]),
    ("descending_wavelengths", dict(wavelengths=_WLS[::-1].copy()),
     ["wavelengths", "strictly increasing", "index 1"]),
    ("nan_angle", dict(angles=[30.0, np.nan]), ["angles", "index 1", "nan"]),
    ("angle_above_90", dict(angles=[95.0]), ["angles", "index 0", "95.0"]),
    ("angle_retrograde", dict(angles=[120.0]), ["angles", "index 0", "120.0"]),
    ("angle_negative", dict(angles=[-30.0]), ["angles", "index 0", "-30.0"]),
    ("angle_radians_out_of_range", dict(angles=[2.0], angles_in_radians=True),
     ["angles", "index 0", "2.0", "rad"]),
]


@pytest.mark.parametrize("kwargs,fragments", [r[1:] for r in _REJECT],
                         ids=[r[0] for r in _REJECT])
def test_bad_input_raises_and_names_the_offender(kwargs, fragments):
    """Every rejected input must raise ``ValueError`` naming index and value.

    The index matters as much as the rejection: "a thickness is negative" sends
    the caller hunting through a 40-layer stack, "layer 2: -50.0" does not.
    """
    with pytest.raises(ValueError) as excinfo:
        _sm(**kwargs)
    message = str(excinfo.value)
    missing = [f for f in fragments if f not in message]
    assert not missing, f"message does not name {missing}: {message!r}"


# --------------------------------------------------------------------------
# Accepted inputs -- the checks must not overreach
# --------------------------------------------------------------------------

_ACCEPT = [
    ("grazing_90_deg", dict(angles=[90.0])),
    ("normal_0_deg", dict(angles=[0.0])),
    ("radians_pi_over_2", dict(angles=[np.pi / 2.0], angles_in_radians=True)),
    ("single_wavelength", dict(wavelengths=np.array([550.0]))),
    ("zero_thickness_everywhere", dict(thicknesses=np.zeros(4))),
    ("huge_thickness_1e9", dict(thicknesses=_with(_D, 2, 1e9))),
    ("strongly_absorbing_index", dict(layer_indices=_with(_N, 1, 1.2 + 7.5j))),
    ("index_below_one", dict(layer_indices=_with(_N, 1, 0.2 + 3.0j))),
    ("two_dim_index_array", dict(layer_indices=np.repeat(_N[:, None], _WLS.size, axis=1))),
]


@pytest.mark.parametrize("kwargs", [a[1] for a in _ACCEPT], ids=[a[0] for a in _ACCEPT])
def test_legitimate_input_still_constructs_and_solves(kwargs):
    """Odd but physical inputs must survive.

    Grazing incidence, a one-point grid, a metallic index with n < 1, a
    millimetre-thick layer: all legal, all previously working, and a
    validation layer that rejects any of them has broken the library to fix a
    bug. The solve is included because construction alone would not catch a
    check that mangles the arrays on the way through.
    """
    out = _sm(**kwargs).compute(Request.RS | Request.TS)
    assert set(out) == {"Rs", "Ts"}
    assert np.all(np.isfinite(np.asarray(out["Rs"], float)))


# --------------------------------------------------------------------------
# Positive controls: rejecting is not the same as repairing
# --------------------------------------------------------------------------

def test_ascending_grid_is_stored_verbatim():
    """Reject, never reorder.

    Silently sorting a descending grid would desynchronize it from the
    caller's own wavelength-indexed arrays -- their n(lambda) table would no
    longer line up with the grid the engine used, and nothing would say so.
    This asserts the grid comes out exactly as it went in.
    """
    st = _sm()
    assert np.array_equal(st.wavls, _WLS)
    assert st.wavls[0] < st.wavls[-1]


def test_a_rejected_construction_does_not_mutate_the_caller_arrays():
    """The validator reads; it must not write.

    A check that repaired its input in place would leave the caller holding a
    silently different array after catching the error.
    """
    wls = _WLS[::-1].copy()
    thick = _with(_D, 2, -50.0)
    idx = _with(_N, 1, np.nan)
    wls_before, thick_before = wls.copy(), thick.copy()
    idx_before = idx.copy()

    for kwargs in (dict(wavelengths=wls), dict(thicknesses=thick),
                   dict(layer_indices=idx)):
        with pytest.raises(ValueError):
            _sm(**kwargs)

    assert np.array_equal(wls, wls_before)
    assert np.array_equal(thick, thick_before)
    assert np.array_equal(idx[~np.isnan(idx)], idx_before[~np.isnan(idx_before)])
    assert np.isnan(idx[1])


def test_the_first_offender_is_the_one_reported():
    """With several bad entries, the message names the first, not an arbitrary one."""
    bad = _with(_with(_D, 1, -1.0), 3, -2.0)
    with pytest.raises(ValueError, match=r"layer 1: -1\.0"):
        _sm(thicknesses=bad)


def test_validation_runs_before_the_native_solver_is_built():
    """Bad input must not reach the engine at all.

    If the native Solver were constructed first, a rejected stack would still
    have paid for (and possibly warned from) the native setup, and a future
    native panic on garbage input would surface instead of the clear error.
    """
    import navette._smatrix as native

    calls = []
    real = native.Solver

    class _Spy:
        def __new__(cls, *a, **kw):
            calls.append(1)
            return real(*a, **kw)

    import navette.smatrix.smatrix as mod
    original = mod._NativeSolver
    mod._NativeSolver = _Spy
    try:
        with pytest.raises(ValueError):
            _sm(thicknesses=_with(_D, 2, -50.0))
        assert calls == [], "native Solver was constructed for a rejected stack"
        _sm()
        assert calls == [1], "native Solver was not constructed for a valid stack"
    finally:
        mod._NativeSolver = original
