# SPDX-License-Identifier: LGPL-3.0-or-later
"""The roughness cap: each interface is solved with min(sigma, d / cap).

A film thinner than its roughness cannot carry two separate rough faces.
The cap makes a film that thins to zero vanish from the answer instead of
leaving a sharp-plus-rough interface pair behind; the stored roughness is
kept as given.
"""

import numpy as np
import pytest

from navette.smatrix.smatrix import ScatterMatrix, Request

WL = np.array([400.0, 550.0, 700.0])


def _rs(rows, **kw):
    # rows: (n, d, sigma) between air and glass; sigma on the row's top face
    n = np.array([1.0] + [r[0] for r in rows] + [1.52], dtype=complex)
    d = np.array([0.0] + [r[1] for r in rows] + [0.0])
    rv = np.array([0.0] + [r[2] for r in rows] + [0.0])
    rt = np.where(rv > 0, 5, 0).astype(np.int32)
    sm = ScatterMatrix(n, d, wavelengths=WL, angles_deg=[0.0],
                       roughness_types=rt, roughness_values=rv, **kw)
    return np.asarray(sm.compute(Request.RS)["Rs"], float).ravel(), sm


BASE = [(1.46, 100.0, 0.0), (2.35, 80.0, 10.0)]
SPLIT = [(1.46, 100.0, 0.0), (2.35, 40.0, 10.0), (1.8, 0.0, 0.0), (2.35, 40.0, 10.0)]


def test_a_vanishing_film_leaves_no_trace_with_the_cap():
    r0, sm = _rs(BASE)
    r1, _ = _rs(SPLIT)
    assert sm.roughness_cap == 2.0
    np.testing.assert_allclose(r1, r0, rtol=0, atol=1e-12)
    # The caller's roughness is kept as given.
    assert sm.roughness_values[2] == 10.0


def test_the_cap_can_be_turned_off():
    r0, _ = _rs(BASE, roughness_cap=0.0)
    r1, _ = _rs(SPLIT, roughness_cap=0.0)
    assert np.max(np.abs(r1 - r0)) > 1e-4


def test_thick_films_are_not_capped():
    r_on, _ = _rs(BASE)
    r_off, _ = _rs(BASE, roughness_cap=0.0)
    np.testing.assert_array_equal(r_on, r_off)


def test_a_non_finite_cap_is_refused():
    with pytest.raises(ValueError, match="roughness_cap"):
        _rs(BASE, roughness_cap=float("nan"))


def test_the_design_stack_carries_the_cap():
    from navette._smatrix import DesignStack, LayerSpec

    nk = lambda v: np.full(WL.shape, v, dtype=complex)  # noqa: E731
    films = [LayerSpec("H", nk(2.35), 3.0, rough_type=5, rough_val=4.0)]
    st = DesignStack(LayerSpec("air", nk(1.0), 0.0),
                     LayerSpec("glass", nk(1.52), 0.0), films)
    assert st.roughness_cap == 2.0
    st.roughness_cap = 0.0
    assert st.roughness_cap == 0.0
    with pytest.raises(ValueError, match="roughness_cap"):
        st.roughness_cap = float("inf")
