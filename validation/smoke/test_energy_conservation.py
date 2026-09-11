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
