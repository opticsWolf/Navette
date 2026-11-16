# SPDX-License-Identifier: LGPL-3.0-or-later
"""A4: per-material thickness bounds on the keyword and document doors.

The config route (`pipeline_from_config`) has carried a material's
`thickness_bounds` since 0.11.0. These are the other two ways into a
needle run: `run_needle(...)` itself (both arms) and a program document
handed through `design_from_program`. Each must clamp a bounded
material's films into its interval, or refuse - never drop it.
"""
import json

import numpy as np
import pytest

from navette.config import load_program
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis.pipeline import (
    PipelineConfig,
    _environments_request,
    design_from_program,
    run_needle,
)

WL = np.linspace(450.0, 650.0, 11)
H, L, G = 2.2 + 0j, 1.46 + 0j, 1.52 + 0j
LO, HI = 60.0, 400.0
#: The flat arm takes evaluated curves (or specs), not bare indices.
LA, HA = np.full(WL.shape, L), np.full(WL.shape, H)

#: Enough cycles for needles to go in, short enough for a regression run.
CFG = dict(max_macro_cycles=3, needles_per_cycle=2)

DOCUMENT = {
    "schema_version": 2, "kind": "program", "name": "bounded",
    "sections": {
        "materials": [
            {"name": "L", "code": "L", "model": "Konstant",
             "params": {"n": 1.46, "k": 0.0}},
            {"name": "H", "code": "H", "model": "Konstant",
             "params": {"n": 2.2, "k": 0.0},
             "thickness_bounds": {"min_nm": LO, "max_nm": HI}},
        ],
        "design": {"coat": {"layers": [
            {"material_code": "L", "thickness_nm": 300.0},
            {"material_code": "H", "thickness_nm": 90.0},
        ]}},
        "environments": [{"name": "bare", "stack": [{"design": "coat"}]}],
    },
}


def _targets(env=None):
    tc = TargetCollection()
    kw = {"environment": env} if env else {}
    tc.add(SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, 1e-3),
                          "s", "R", angle_deg=0.0, **kw))
    return tc


def _write(tmp_path, doc):
    p = tmp_path / "prog.json"
    p.write_text(json.dumps(doc), encoding="utf-8")
    return str(p)


def _h_films(out):
    films = out["stack"].to_dict()["films"]
    return films, [f["thickness"] for f in films
                   if f["material"].split("#")[0] == "H"]


def _assert_bounded(out):
    films, hs = _h_films(out)
    assert len(films) > 2, films           # needles went in
    assert hs, films
    assert all(d == 0.0 or LO - 1e-9 <= d <= HI + 1e-9 for d in hs), hs


# --------------------------------------------------------------------------
# Program documents
# --------------------------------------------------------------------------

@pytest.mark.parametrize("context", [None, {"groups": {}}],
                         ids=["native-load", "context-load"])
def test_a_document_interval_reaches_the_request(tmp_path, context):
    """Both load paths keep the section's interval, and it lands on the
    library entry the core reads (`library[].thickness_bounds`)."""
    prog = load_program(_write(tmp_path, DOCUMENT), WL, context=context)
    assert prog.material_bounds == {"H": (LO, HI)}
    design, environments = design_from_program(prog)
    req = _environments_request(design, environments, None, WL,
                                None, None, None, None, None)
    lib = {e["code"]: e for e in req["library"]}
    assert lib["H"]["thickness_bounds"] == {"min_nm": LO, "max_nm": HI}
    assert "thickness_bounds" not in lib["L"]


def test_a_document_interval_holds_through_a_needle_run(tmp_path):
    prog = load_program(_write(tmp_path, DOCUMENT), WL)
    design, environments = design_from_program(prog)
    out = run_needle(targets=_targets("bare"), angles_deg=[0.0], wavelengths=WL,
                     design=design, environments=environments,
                     pipeline_config=PipelineConfig(**CFG), substrate=(G, "sub"))
    _assert_bounded(out)


# --------------------------------------------------------------------------
# The keyword, flat arm
# --------------------------------------------------------------------------

def _flat(**kw):
    return run_needle(layers=[(LA, 300.0), (HA, 90.0)], names=["L", "H"],
                      targets=_targets(), angles_deg=[0.0], wavelengths=WL,
                      pipeline_config=PipelineConfig(**CFG), **kw)


def test_the_keyword_bounds_a_flat_run():
    _assert_bounded(_flat(thickness_bounds={"H": (LO, HI)}))


def test_a_bounded_mapping_material_bounds_a_flat_run():
    h = {"model": "Konstant", "params": {"n": 2.2, "k": 0.0},
         "thickness_bounds": {"min_nm": LO, "max_nm": HI}}
    out = run_needle(layers=[(LA, 300.0), (h, 90.0)], names=["L", "H"],
                     targets=_targets(), angles_deg=[0.0], wavelengths=WL,
                     pipeline_config=PipelineConfig(**CFG))
    _assert_bounded(out)


def test_the_keyword_bounds_a_pool_label_on_a_flat_run():
    out = _flat(materials={"L": LA, "H": HA}, thickness_bounds={"H": (LO, HI)})
    _assert_bounded(out)


def test_an_unknown_label_refuses_on_the_flat_arm():
    with pytest.raises(ValueError, match="names no film and no needle material"):
        _flat(thickness_bounds={"X": (LO, HI)})


def test_two_intervals_for_one_label_refuse_on_the_flat_arm():
    h = {"model": "Konstant", "params": {"n": 2.2, "k": 0.0},
         "thickness_bounds": (LO, HI)}
    with pytest.raises(ValueError, match="two different intervals"):
        run_needle(layers=[(LA, 300.0), (h, 90.0)], names=["L", "H"],
                   targets=_targets(), angles_deg=[0.0], wavelengths=WL,
                   thickness_bounds={"H": (10.0, HI)})


def test_an_empty_interval_refuses():
    with pytest.raises(ValueError):
        _flat(thickness_bounds={"H": (HI, LO)})


# --------------------------------------------------------------------------
# The keyword, environments arm
# --------------------------------------------------------------------------

def _envs(**kw):
    return run_needle(targets=_targets("bare"), angles_deg=[0.0], wavelengths=WL,
                      design={"coat": [(L, 300.0, "L"), (H, 90.0, "H")]},
                      environments=[{"name": "bare", "stack": [{"design": "coat"}]}],
                      pipeline_config=PipelineConfig(**CFG), substrate=(G, "sub"),
                      **kw)


def test_the_keyword_bounds_an_environments_run():
    _assert_bounded(_envs(thickness_bounds={"H": (LO, HI)}))


def test_an_unknown_label_refuses_on_the_environments_arm():
    with pytest.raises(ValueError, match="names no film and no needle material"):
        _envs(thickness_bounds={"X": (LO, HI)})
