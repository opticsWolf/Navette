# SPDX-License-Identifier: LGPL-3.0-or-later
"""ScatterMatrix.energy_conservation() coverage (R1.2, review §15).

The wrapper passes 2-D ``[n_angles, n_wavs]`` arrays (squeeze=False) to the
native ``solver_energy_conservation``, which previously accepted only 1-D input
and raised ``TypeError`` on every call. The binding now accepts 2-D and
reshapes the elementwise result back, so the advertised API works and the
single-angle path still squeezes to 1-D.
"""

import numpy as np

from navette.smatrix.smatrix import ScatterMatrix, Request

_WLS = np.linspace(450.0, 750.0, 7)


def test_lossless_stack_conserves():
    """3-layer lossless oxide -> conservation residual ~ 0 for both pols."""
    n = np.array([1.0 + 0j, 1.46 + 0j, 2.10 + 0j, 1.52 + 0j])
    d = np.array([0.0, 120.0, 80.0, 0.0])
    sm = ScatterMatrix(n, d, wavelengths=_WLS, angles=[0.0, 30.0, 60.0])
    cons = sm.energy_conservation()
    assert cons.shape == (3, len(_WLS))
    assert float(np.max(cons)) < 1e-12


def test_two_d_shape_preserved():
    """Multi-angle input keeps the 2-D [n_angles, n_wavs] shape."""
    n = np.array([1.0 + 0j, 1.46 + 0j, 1.52 + 0j])
    d = np.array([0.0, 120.0, 0.0])
    sm = ScatterMatrix(n, d, wavelengths=_WLS, angles=[0.0, 15.0])
    assert sm.energy_conservation().shape == (2, len(_WLS))


def test_single_angle_squeezes_to_1d():
    """A single angle returns a 1-D array (squeeze semantics preserved)."""
    n = np.array([1.0 + 0j, 1.46 + 0j, 1.52 + 0j])
    d = np.array([0.0, 120.0, 0.0])
    sm = ScatterMatrix(n, d, wavelengths=_WLS, angles=[0.0])
    cons = sm.energy_conservation()
    assert cons.shape == (len(_WLS),)


def test_absorbing_equals_absorptance():
    """For an absorbing stack the residual equals the absorptance A = 1-R-T,
    and A grows monotonically with the absorbing layer thickness within [0, 1]."""
    n = np.array([1.0 + 0j, 1.46 + 0.05j, 1.52 + 0j])
    prev = -1.0
    for thick in (0.0, 50.0, 200.0, 800.0):
        d = np.array([0.0, thick, 0.0])
        sm = ScatterMatrix(n, d, wavelengths=np.array([550.0]), angles=[0.0])
        cons = float(sm.energy_conservation()[0])
        out = sm.compute(Request.A_AVG | Request.A_S | Request.A_P)
        a_max = max(float(out["A_s"][0]), float(out["A_p"][0]))
        assert 0.0 <= cons <= 1.0
        # residual is max over pols of |1-R-T| == max pol absorptance
        assert abs(cons - a_max) < 1e-12
        assert cons >= prev - 1e-12
        prev = cons


# --------------------------------------------------------------------------
# The native binding's own contract (0.6.25). R1.2 widened it from 1-D to 2-D
# by changing the parameter type, which silently dropped the 1-D form -- and
# PyO3 reports the mismatch as "'ndarray' object is not an instance of
# 'ndarray'", which tells a caller nothing. It now takes both.
# --------------------------------------------------------------------------


def test_native_accepts_1d_and_returns_1d():
    from navette._smatrix import solver_energy_conservation as f

    rs = np.array([0.04, 0.10, 0.25])
    ts = np.array([0.96, 0.90, 0.75])
    out = f(rs, rs, ts, ts)
    assert out.shape == (3,)
    assert float(np.max(np.abs(out))) < 1e-15


def test_native_returns_the_callers_own_shape():
    from navette._smatrix import solver_energy_conservation as f

    for shape in ((3, 7), (5, 2), (1, 4), (6,)):
        a = np.full(shape, 0.25)
        b = np.full(shape, 0.75)
        assert f(a, a, b, b).shape == shape


def test_native_rejects_mismatched_shapes_and_ranks():
    import pytest
    from navette._smatrix import solver_energy_conservation as f

    a = np.zeros((2, 3))
    with pytest.raises(ValueError, match="shapes must match"):
        f(a, np.zeros((3, 2)), a, a)
    cube = np.zeros((2, 3, 4))
    with pytest.raises(ValueError, match="expected 1-D"):
        f(cube, cube, cube, cube)


def test_native_is_the_max_over_polarizations_not_a_row_per_pol():
    """The one thing the `.pyi` used to get wrong, pinned.

    The docstring claimed "`A = 1 - R - T` per polarization, as a `(2, n)`
    array". It is a single residual per grid point, maxed over the two
    polarizations, in the caller's own shape.
    """
    from navette._smatrix import solver_energy_conservation as f

    rs, ts = np.array([0.10]), np.array([0.80])   # |1-Rs-Ts| = 0.10
    rp, tp = np.array([0.05]), np.array([0.60])   # |1-Rp-Tp| = 0.35
    out = f(rs, rp, ts, tp)
    assert out.shape == (1,)
    assert abs(float(out[0]) - 0.35) < 1e-15
