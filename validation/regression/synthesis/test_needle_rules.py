# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Needle rules: which pool materials go into which, and the edge layers.

`insert_rules`, `outer` and `inner` narrow a needle run on every door:
`run_needle` (flat and multi-environment), `NeedlePipeline`, and the
config builder. The mechanics (scan filter, cleanup and clamp guards)
are pinned in the core; these are the seams:

  * a rule is obeyed end to end, and no rules is the old run bit for bit;
  * the edges the caller named are the edges the run ends with; an edge
    film matches by label, the floor reports an edge it holds, and
    `seed_edges` puts a listed film on an edge that starts off its list;
  * every refusal is raised at the door, naming what was wrong.
"""
import numpy as np
import pytest

from navette._smatrix import NeedlePipeline
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec
from navette.synthesis.pipeline import PipelineConfig, run_needle, stack_from_layers

WL = np.linspace(450.0, 750.0, 31)
N = {"L": 1.46 + 0j, "M": 1.65 + 0j, "H": 2.35 + 0j}
MATERIALS = {k: np.full(WL.shape, v) for k, v in N.items()}
SEED = [(MATERIALS["L"], 300.0), (MATERIALS["M"], 200.0)]
NAMES = ["L", "M"]


def _targets():
    tc = TargetCollection()
    tc.add(SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, 0.01),
                          0.0, "s", "R"))
    return tc


def _cfg():
    return PipelineConfig(max_macro_cycles=3, needles_per_cycle=2)


def _run(**kw):
    return run_needle(layers=SEED, names=NAMES, targets=_targets(),
                      angles_deg=[0.0], wavelengths=WL, materials=MATERIALS,
                      pipeline_config=_cfg(), substrate=(1.52, "sub"), **kw)


def _inserted(res):
    return [nr["insertion"]["material"] for p in res["phases"]
            for nr in p["needle_results"] if nr["insertion"]]


def _labels(res):
    """Each film's material label, read from its index."""
    out = []
    for f in res["stack"].films():
        n = complex(np.asarray(f["nk"])[0])
        out.append(next(k for k, v in N.items() if v == n))
    return out


def test_no_rules_is_the_default_run_bit_for_bit():
    base = _run()
    assert _inserted(base), "the problem must insert something"
    for kw in (dict(insert_rules={}), dict(insert_rules=None, outer=None)):
        res = _run(**kw)
        assert res["final_mf"].hex() == base["final_mf"].hex()
        assert res["final_layer_count"] == base["final_layer_count"]


def test_a_closed_material_never_takes_a_needle():
    # M takes nothing and is never inserted: it stays one film, untouched
    # by any needle, while L and H trade needles.
    res = _run(insert_rules={"M": [], "L": ["H"], "H": ["L"]})
    ins = _inserted(res)
    assert ins and set(ins) <= {"L", "H"}, ins
    assert _labels(res).count("M") <= 1
    assert not any(f["material"].startswith("M#") for f in res["stack"].films())


def test_the_edges_hold_to_the_end():
    res = _run(outer="L", inner=["M", "H"])
    labels = _labels(res)
    assert labels[0] == "L" and labels[-1] in ("M", "H"), labels


def test_the_rules_reach_the_pipeline_class():
    stack, pool = stack_from_layers(SEED, WL, MATERIALS, names=NAMES)
    pipe = NeedlePipeline(stack, build_merit_spec(_targets()), np.array([0.0]),
                          WL, pool, _cfg(), insert_rules={"M": [], "L": ["H"]},
                          outer=["L"])
    res = pipe.run()
    assert set(_inserted(res)) <= {"H", "L"}
    assert _labels(res)[0] == "L"
    with pytest.raises(ValueError, match="outer"):
        NeedlePipeline(stack, build_merit_spec(_targets()), np.array([0.0]),
                       WL, pool, _cfg(), outer=["H"])


@pytest.mark.parametrize("kw, match", [
    (dict(outer="H"), r"outer \(ambient side\) edge film \"L\""),
    (dict(inner="L"), r"inner \(substrate side\)"),
    (dict(outer="Q"), r"\"Q\" is not a pool material"),
    (dict(insert_rules={"L": ["Q"]}), r"\"Q\" is not a pool material"),
    (dict(insert_rules={"L": ["L"]}), "lists itself"),
])
def test_bad_rules_refuse_at_the_door(kw, match):
    with pytest.raises(ValueError, match=match):
        _run(**kw)


def test_the_clamp_reports_an_edge_it_holds():
    # A pinned 1 nm L on top, under the 2 nm floor: removing it would put H
    # outside, so the floor holds it and says so.
    seed = [(MATERIALS["L"], 1.0), (MATERIALS["H"], 100.0), (MATERIALS["M"], 200.0)]
    res = run_needle(layers=seed, names=["L", "H", "M"], targets=_targets(),
                     angles_deg=[0.0], wavelengths=WL, materials=MATERIALS,
                     pipeline_config=_cfg(), substrate=(1.52, "sub"), outer="L",
                     per_film_flags={"L": {"optimize": False}})
    held = [h for p in res["phases"] if p.get("clamp_report")
            for h in p["clamp_report"]["edges_held"]]
    assert held and held[0] == "outer: L (1.0 nm)", held
    assert _labels(res)[0] == "L"


def test_an_edge_film_matches_by_label():
    # The film's index is off by a part in 1e12 (a rescaled spectrum): it
    # is still "L", so the edge list takes it.
    seed = [(MATERIALS["L"] * (1 + 1e-12), 300.0), (MATERIALS["M"], 200.0)]
    res = run_needle(layers=seed, names=NAMES, targets=_targets(),
                     angles_deg=[0.0], wavelengths=WL, materials=MATERIALS,
                     pipeline_config=_cfg(), substrate=(1.52, "sub"), outer="L")
    assert _labels_near(res)[0] == "L"


def _labels_near(res):
    out = []
    for f in res["stack"].films():
        n = complex(np.asarray(f["nk"])[0])
        out.append(min(N, key=lambda k: abs(N[k] - n)))
    return out


def test_seed_edges_put_a_listed_film_on_an_edge_off_its_list():
    with pytest.raises(ValueError, match="outer"):
        _run(outer="H")
    res = _run(outer="H", seed_edges=2.0)
    assert res["stack"].films()[0]["material"] == "H#1"
    assert _labels(res)[0] == "H"
    assert "M" in _labels(res)
    with pytest.raises(ValueError, match="seed_edges"):
        _run(outer="H", seed_edges=0.0)


def test_the_pipeline_class_seeds_edges_too():
    stack, pool = stack_from_layers(SEED, WL, MATERIALS, names=NAMES)
    pipe = NeedlePipeline(stack, build_merit_spec(_targets()), np.array([0.0]),
                          WL, pool, _cfg(), outer=["H"], seed_edges=2.0)
    assert _labels(pipe.run())[0] == "H"


def test_insert_rules_must_be_a_mapping():
    with pytest.raises(TypeError, match="insert_rules must be a mapping"):
        _run(insert_rules=[("L", "H")])


def test_the_environments_door_takes_the_rules():
    """The edges are the shared design's, whatever surrounds it."""
    G = 1.52 + 0j
    kw = dict(
        targets=_targets(), angles_deg=[0.0], wavelengths=WL,
        design={"coat": [(N["L"], 100.0, "L"), (N["H"], 80.0, "H")]},
        environments=[
            {"name": "bare", "stack": [{"design": "coat"}]},
            {"name": "laminated", "stack": [
                {"layers": [(G, 300.0, "G")]}, {"design": "coat"}]},
        ],
        pipeline_config=PipelineConfig(max_macro_cycles=1, needles_per_cycle=2,
                                       enable_cleanup=False, enable_inflate=False,
                                       thin_layer_policy="clamp_up_final"),
        substrate=(G, "sub"))
    # Only H takes needles (of L); L is closed. The joint run inserts.
    res = run_needle(outer="L", inner="H", insert_rules={"H": ["L"], "L": []}, **kw)
    ins = _inserted(res)
    assert ins and set(ins) == {"L"}, ins
    # Closing H instead leaves only H needles in L, and none of them
    # improves this merit: the rule changes the run, nothing lands.
    assert not _inserted(run_needle(insert_rules={"L": ["H"], "H": []}, **kw))
    with pytest.raises(ValueError, match="outer"):
        run_needle(outer="H", **kw)
    # Seeded, the H edge film goes on top in both environments.
    res = run_needle(outer="H", seed_edges=2.0, **kw)
    assert res["stack"].films()[0]["material"] == "H#1"


def test_edge_lists_belong_to_a_segment():
    """With two segments each one's lists are checked on its own films,
    in an order that differs per environment."""
    G = 1.52 + 0j
    cfg = PipelineConfig(max_macro_cycles=1, needles_per_cycle=1,
                         enable_cleanup=False, enable_inflate=False,
                         thin_layer_policy="clamp_up_final")

    def run(top_outer, **extra):
        design = {
            "top": {"layers": [(N["L"], 100.0, "L"), (N["H"], 80.0, "H")],
                    "outer": top_outer, "inner": "H"},
            "bottom": {"layers": [(N["M"], 60.0, "M")], "outer": ["M"]},
        }
        return run_needle(
            targets=_targets(), angles_deg=[0.0], wavelengths=WL,
            design=design,
            environments=[
                {"name": "a", "stack": [{"design": "top"}, {"design": "bottom"}]},
                {"name": "b", "stack": [{"layers": [(G, 300.0, "G")]},
                                        {"design": "bottom"}, {"design": "top"}]},
            ],
            pipeline_config=cfg, substrate=(G, "sub"), **extra)

    assert "stack" in run("L")
    with pytest.raises(ValueError, match=r'design segment "top".*outer'):
        run("H")
    with pytest.raises(ValueError, match="ambiguous"):
        run("L", outer="L")
    with pytest.raises(ValueError, match="'layers'"):
        run_needle(targets=_targets(), angles_deg=[0.0], wavelengths=WL,
                   design={"top": {"films": []}},
                   environments=[{"name": "a", "stack": [{"design": "top"}]}])
