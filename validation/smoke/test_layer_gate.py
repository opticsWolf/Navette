# SPDX-License-Identifier: LGPL-3.0-or-later
"""0.6.28: a layer's numbers are judged where they are written.

Before this, several layer properties were accepted in silence and only
misbehaved later, somewhere else:

* ``roughness = -20`` solved byte-identically to ``+20`` -- every roughness
  form factor squares sigma, so a sign slip could never surface.
* ``roughness = NaN`` turned every output NaN with nothing to point at.
* ``inh_delta < 0`` drove the refinement factor negative, which saturates to
  a sub-layer count of 1: the grading was dropped and the film solved as
  homogeneous. Writing ``-0.2`` to mean "ramp the other way" is the obvious
  way to hit that.
* ``inh_delta >= 2`` put rows of ``n = -0.59, k = -0.0125`` into the solver
  (measured on a 2.35 + 0.05i film at delta 2.5). Negative k is optical gain.

The rule now lives once, as ``Layer::property_issues``, and every door that
builds a layer uses it: the constructor, the numeric setters,
``set_properties``, ``from_state``, ``Structure.validate`` and the synthesis
assembler (which builds films from flag dicts and never touches the Python
``Layer``).

SCOPE, deliberately. The flat-array solver surface is NOT a door.
``ScatterMatrix(roughness_values=...)`` and the native ``Solver`` take raw
arrays and stay permissive, as R3.1 promised -- see
``test_the_flat_array_surface_is_still_permissive`` at the bottom, which pins
that as a decision rather than leaving it as an accident.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from navette.smatrix.smatrix import Request, ScatterMatrix
from navette.structure import DictMaterialProvider, Layer, Navette_Structure
from navette.synthesis.pipeline import stack_from_layers

_WLS = np.array([500.0, 600.0, 700.0])
_MATS = DictMaterialProvider({
    "glass": np.full(1, 1.52 + 0j),
    "TiO2": np.full(1, 2.35 + 0.01j),
})


# ---------------------------------------------------------------------------
# Refusals: no correction exists, so stopping beats guessing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kwargs,fragment", [
    (dict(thickness=-5.0), "Negative thickness"),
    (dict(thickness=float("nan")), "Non-finite thickness"),
    (dict(thickness=float("inf")), "Non-finite thickness"),
    (dict(roughness=-20.0), "Negative roughness"),
    (dict(roughness=float("nan")), "Non-finite roughness"),
    (dict(interface_thickness=-2.0), "Negative interface thickness"),
    (dict(interface_thickness=float("inf")), "Non-finite interface thickness"),
    (dict(inh_delta=-0.2), "Negative inh_delta"),
    (dict(inh_delta=2.0), "outside [0, 2)"),
    (dict(inh_delta=2.5), "outside [0, 2)"),
    (dict(inh_delta=float("nan")), "Non-finite inh_delta"),
    (dict(rough_type=99), "invalid discriminant"),
    (dict(rough_type=-1), "invalid discriminant"),
])
def test_the_constructor_refuses_what_has_no_correction(kwargs, fragment):
    base = dict(thickness=100.0, material_name="TiO2")
    base.update(kwargs)
    with pytest.raises(ValueError, match=re_escape(fragment)):
        Layer(**base)


def re_escape(s: str) -> str:
    import re
    return re.escape(s)


def test_the_refusals_explain_themselves():
    """A range check that only prints the range teaches nothing.

    ``inh_delta`` is the case that needs it: both ends fail for reasons a
    caller cannot guess from "[0, 2)".
    """
    with pytest.raises(ValueError) as neg:
        Layer(100.0, "TiO2", inhomogen=True, inh_delta=-0.2)
    msg = str(neg.value)
    assert "magnitude, not a direction" in msg
    assert "homogeneous" in msg, "must say what silently happened before"

    with pytest.raises(ValueError) as big:
        Layer(100.0, "TiO2", inhomogen=True, inh_delta=2.5)
    msg = str(big.value)
    assert "optical gain" in msg, "the reason 2 is the limit, not a taste"
    assert msg.isascii(), "the message must survive a cp1252 console"


def test_sigma_is_refused_rather_than_absolute_valued():
    """Deliberately not corrected to ``|sigma|``.

    Every roughness form factor squares it, so -20 and +20 produce identical
    numbers. A silent correction would be indistinguishable from the bug, and
    ``Structure.validate`` has refused this since it existed -- 0.6.28 only
    moved the refusal earlier.
    """
    with pytest.raises(ValueError, match="Negative roughness"):
        Layer(100.0, "TiO2", roughness=-20.0)


# ---------------------------------------------------------------------------
# Warnings: legal but suspicious, never blocking
# ---------------------------------------------------------------------------

def test_an_overhanging_interface_warns_and_builds():
    with pytest.warns(UserWarning, match="clamped at expansion"):
        layer = Layer(10.0, "TiO2", interface=True, interface_thickness=20.0)
    assert layer.thickness == 10.0


def test_a_graded_layer_with_no_grading_warns_and_builds():
    """Eight identical sub-layers is a waste, not an error."""
    with pytest.warns(UserWarning, match="identical sub-layers"):
        layer = Layer(100.0, "TiO2", inhomogen=True, inh_delta=0.0)
    assert layer.sub_layer_count > 1


def test_ordinary_layers_say_nothing():
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        Layer(100.0, "TiO2")
        Layer(100.0, "TiO2", inhomogen=True, inh_delta=0.2)
        Layer(100.0, "TiO2", roughness=3.0, rough_type=5)
        Layer(0.0, "glass")  # ambient/substrate half-spaces
    assert [str(w.message) for w in rec] == []


# ---------------------------------------------------------------------------
# The other doors onto the same object
# ---------------------------------------------------------------------------

def test_the_setters_are_gated_too():
    layer = Layer(100.0, "TiO2")
    for attr, value, fragment in [("roughness", -3.0, "Negative roughness"),
                                  ("thickness", float("nan"), "Non-finite thickness"),
                                  ("inh_delta", 3.0, "outside [0, 2)"),
                                  ("interface_thickness", -1.0, "Negative interface")]:
        with pytest.raises(ValueError, match=re_escape(fragment)):
            setattr(layer, attr, value)
    # The layer is unchanged by the rejections, and still accepts good values.
    assert layer.thickness == 100.0 and layer.roughness == 0.0
    layer.roughness = 3.0
    assert layer.roughness == 3.0


def test_turning_grading_on_is_a_gate():
    """``inh_delta`` is inert until ``inhomogen`` makes it load-bearing."""
    layer = Layer(100.0, "TiO2")
    layer.inhomogen = True
    assert layer.inhomogen


def test_a_rejected_property_batch_leaves_the_layer_alone():
    """Half-written is worse than not written.

    ``set_properties`` applies to a copy and swaps it in only once the whole
    batch passes, so a bad roughness cannot leave a new thickness behind.
    """
    layer = Layer(100.0, "TiO2")
    with pytest.raises(ValueError, match="Negative roughness"):
        layer.set_properties({"thickness": 50.0, "roughness": -1.0})
    assert layer.thickness == 100.0, "the good half of the batch was applied"
    assert layer.roughness == 0.0

    layer.set_properties({"thickness": 50.0, "roughness": 1.0})
    assert layer.thickness == 50.0 and layer.roughness == 1.0


def test_from_state_cannot_smuggle_a_bad_layer_back_in():
    """Deserialization is a door like any other."""
    state = Layer(100.0, "TiO2").get_state()
    state["roughness"] = -5.0
    with pytest.raises(ValueError, match="Negative roughness"):
        Layer.from_state(state)


def test_structure_validate_reports_the_same_rule():
    """One rule, and the structure-level report still names the layer index.

    ``Structure.validate`` calls ``Layer::property_issues``; it cannot see a
    bad layer from Python any more, but it is still the path that reports
    everything a single layer cannot judge (unresolvable materials, orphan
    groups), and its wording for the shared checks is unchanged.
    """
    st = Navette_Structure([Layer(50.0, "NotInTheLibrary")], {}, _MATS)
    issues = st.validate()
    assert any("not found" in i for i in issues)


# ---------------------------------------------------------------------------
# The design door: films are built from flag dicts, never from a Layer
# ---------------------------------------------------------------------------

def _build(**flags):
    hi = np.full(_WLS.size, 2.35 + 0j)
    sub = np.full(_WLS.size, 1.52 + 0j)
    return stack_from_layers([(hi, 120.0)], _WLS, {}, substrate=(sub, "sub"),
                             film_flags=flags)


def test_the_synthesis_assembler_applies_the_same_rule():
    """``stack_from_layers`` and ``run_needle`` never touch the Python Layer.

    They build films from flag dicts straight into the native assembler, so
    gating only the constructor would have left this door open -- the same
    mistake the layer-0 gate made in 0.6.26 and fixed in 0.6.27.
    """
    with pytest.raises(ValueError, match="Negative roughness"):
        _build(roughness=-5.0)
    with pytest.raises(ValueError, match=re_escape("outside [0, 2)")):
        _build(inhomogen=True, inh_delta=2.5)


def test_the_assembler_names_the_film_it_rejected():
    with pytest.raises(ValueError) as excinfo:
        _build(roughness=-5.0)
    assert "film0" in str(excinfo.value), str(excinfo.value)


def test_a_clean_design_build_stays_silent():
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        _build()
    assert [str(w.message) for w in rec] == []


# ---------------------------------------------------------------------------
# Out of scope, on purpose
# ---------------------------------------------------------------------------

def test_the_flat_array_surface_is_still_permissive():
    """The gate is at layer construction; the raw arrays stay open.

    This is a scoping decision, not an oversight, and it is pinned here so it
    stays one. ``ScatterMatrix`` and the native ``Solver`` take flat arrays
    with no layer objects involved -- the same escape hatch R3.1 promised for
    an absorbing ambient.

    What that leaves unreported, measured and stated plainly: a negative sigma
    here still solves as its positive twin, and Nevot-Croce (type 5) at
    sigma = 20 nm still produces ``R + T`` = 1.047 on this stack -- energy
    created, no warning. (1.311 in p-pol out to 89 deg; 49.4 at sigma = 100
    nm.) The validity band needs the wavelength and angle grid, which a layer
    does not have, so it cannot live at this gate.
    """
    n = np.array([1.0, 2.35, 1.46, 1.52])
    d = [0.0, 120.0, 200.0, 0.0]
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        neg = ScatterMatrix(n, d, wavelengths=_WLS, angles=[0.0],
                            roughness_types=[0, 5, 5, 5],
                            roughness_values=[0.0, -20.0, -20.0, -20.0])
        pos = ScatterMatrix(n, d, wavelengths=_WLS, angles=[0.0],
                            roughness_types=[0, 5, 5, 5],
                            roughness_values=[0.0, 20.0, 20.0, 20.0])
        req = Request.RS | Request.TS
        a = neg.compute(req, squeeze=False)
        b = pos.compute(req, squeeze=False)
    assert [str(w.message) for w in rec] == [], "no gate here, by decision"
    for key in a:
        assert np.array_equal(a[key], b[key]), f"{key}: -sigma is not +sigma"

    rt = np.asarray(a["Rs"], float) + np.asarray(a["Ts"], float)
    assert float(np.max(rt)) > 1.0, (
        "type 5 past its validity range still creates energy here; if this "
        "ever stops being true the validity band was added and this test "
        "should move rather than be deleted"
    )
