# SPDX-License-Identifier: LGPL-3.0-or-later
"""W7: a film's material is its label, not its index.

Two deposited materials can share one fitted index and still differ (another
group, other insertion rules). Films and pool entries carry the material's
label; an unlabelled film is matched by its index, as before.
"""

import json

import numpy as np
import pytest

from navette._smatrix import DesignStack, LayerSpec

WL = np.array([500.0, 600.0])


def _nk(v):
    return np.full(WL.shape, v, dtype=complex)


def _konstant(name, n):
    return {"name": name, "model": "Konstant", "params": {"n": n}}


def _film(code, d):
    return {
        "material_code": code, "thickness_nm": d, "coherent": True,
        "roughness_nm": 0.0, "rough_type": 0, "inhomogen": False,
        "inh_delta": 0.1, "interface": False, "interface_thickness_nm": 0.0,
        "optimize": True, "needle": True, "layer_type": 1, "gradient": None,
    }


def _design(codes, materials=None, groups=()):
    request = {
        "structure": {
            "label": "t",
            "layers": [_film(c, 50.0) for c in codes],
            "groups": list(groups),
        },
        "library": [_konstant("L", 1.46), _konstant("H", 2.1), _konstant("Hb", 2.1)],
        "materials": materials,
        "ambient_name": "air",
        "substrate_name": "sub",
    }
    return DesignStack.design_from_config(json.dumps(request), WL)


def test_a_layer_spec_carries_its_label():
    assert LayerSpec("H", _nk(2.1), 10.0).label is None
    assert LayerSpec("H#1", _nk(2.1), 10.0, label="H").label == "H"


def test_design_films_are_labelled_by_their_code():
    stack, _ = _design(["L", "H"])
    assert [f["label"] for f in stack.films()] == ["L", "H"]


def test_two_codes_of_one_index_are_two_pool_entries():
    _, pool = _design(["L", "H"], materials=["L", "H", "Hb", "H"])
    assert [m.label for m in pool] == ["L", "H", "Hb"]
    np.testing.assert_array_equal(pool[1].nk, pool[2].nk)


def test_a_given_pool_takes_the_group_factors():
    group = {"name": "H", "n_factor": 1.1}
    stack, default = _design(["L", "H"], groups=[group])
    _, given = _design(["L", "H"], materials=["L", "H"], groups=[group])
    for d, g in zip(default, given):
        assert d.label == g.label
        np.testing.assert_array_equal(d.nk, g.nk)
    np.testing.assert_allclose(given[1].nk.real, 2.1 * 1.1, rtol=0, atol=1e-12)


def test_an_unknown_pool_code_is_refused():
    with pytest.raises(ValueError, match='"Q"'):
        _design(["L"], materials=["Q"])
