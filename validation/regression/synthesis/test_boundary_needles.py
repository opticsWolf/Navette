# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Boundary needles: `NeedleCycleConfig(boundary_needles=True)`.

The scan also tries every span boundary, the outer edges included, and a
needle there is a new film between two layers. The mechanics (sites,
neighbour filter, span bookkeeping, P continuity) are pinned in the core;
these are the seams:

  * off is the old run bit for bit;
  * on, a film too thin to split is reached through its boundaries, and
    the report says which insertions were boundary needles;
  * the edge lists and the insertion rules hold on boundaries too;
  * the multi-environment door refuses the option rather than ignore it.
"""
import numpy as np
import pytest

from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis.pipeline import NeedleCycleConfig, PipelineConfig, run_needle

WL = np.linspace(450.0, 750.0, 31)
N = {"L": 1.46 + 0j, "M": 1.65 + 0j, "H": 2.35 + 0j}
MATERIALS = {k: np.full(WL.shape, v) for k, v in N.items()}


def _targets():
    tc = TargetCollection()
    tc.add(SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, 0.01),
                          0.0, "s", "R"))
    return tc


def _run(seed, names, boundary=None, **kw):
    if boundary is not None:
        kw["needle_config"] = NeedleCycleConfig(boundary_needles=boundary)
    return run_needle(layers=[(MATERIALS[n], d) for n, d in seed], names=names,
                      targets=_targets(), angles_deg=[0.0], wavelengths=WL,
                      materials=MATERIALS, substrate=(1.52, "sub"),
                      pipeline_config=PipelineConfig(max_macro_cycles=2,
                                                     needles_per_cycle=2),
                      **kw)


def _insertions(res):
    return [nr["insertion"] for p in res["phases"]
            for nr in p["needle_results"] if nr["insertion"]]


def _labels(res):
    out = []
    for f in res["stack"].films():
        n = complex(np.asarray(f["nk"])[0])
        out.append(next(k for k, v in N.items() if v == n))
    return out


AR = ([("L", 300.0), ("M", 200.0)], ["L", "M"])
THIN = ([("M", 3.0)], ["M"])
# Held at 3 nm: free, the optimizer would thin M (it only adds reflection
# on glass) and remove it, leaving no film to scan.
PIN = {"per_film_flags": {"M": {"optimize": False}}}


def test_off_is_the_old_run_bit_for_bit():
    base = _run(*AR)
    off = _run(*AR, boundary=False)
    assert off["final_mf"].hex() == base["final_mf"].hex()
    assert off["final_layer_count"] == base["final_layer_count"]
    assert all(i["boundary"] is False for i in _insertions(off))


def test_a_film_too_thin_to_split_is_reached_on_its_boundaries():
    assert not _insertions(_run(*THIN, boundary=False, **PIN))
    res = _run(*THIN, boundary=True, **PIN)
    ins = _insertions(res)
    assert ins, "no boundary needle landed on the 3 nm film"
    assert ins[0]["boundary"] and ins[0]["depth_into_layer_nm"] == 0.0
    assert res["final_mf"] < _run(*THIN, boundary=False, **PIN)["final_mf"]
    labels = _labels(res)
    assert "M" in labels
    assert all(a != b for a, b in zip(labels, labels[1:])), labels


def test_the_edge_lists_hold_with_boundary_needles():
    # The free run's first needle is L on the substrate side; listing only
    # M for the inner edge keeps the bottom M.
    free = _insertions(_run(*THIN, boundary=True, **PIN))
    assert (free[0]["material"], free[0]["film_idx"], free[0]["boundary"]) == ("L", 1, True)
    res = _run(*THIN, boundary=True, outer=["M", "H"], inner=["M"], **PIN)
    labels = _labels(res)
    assert labels[0] in ("M", "H") and labels[-1] == "M", labels


def test_the_insertion_rules_hold_on_boundaries():
    # M is closed, so nothing goes next to it: its two boundaries are the
    # only sites a 3 nm film has.
    assert not _insertions(_run(*THIN, boundary=True, insert_rules={"M": []}, **PIN))


def test_the_environments_door_refuses_the_option():
    G = 1.52 + 0j
    with pytest.raises(ValueError, match="boundary_needles"):
        run_needle(
            targets=_targets(), angles_deg=[0.0], wavelengths=WL,
            design={"coat": [(N["L"], 100.0, "L"), (N["H"], 80.0, "H")]},
            environments=[{"name": "bare", "stack": [{"design": "coat"}]},
                          {"name": "laminated", "stack": [
                              {"layers": [(G, 300.0, "G")]}, {"design": "coat"}]}],
            pipeline_config=PipelineConfig(max_macro_cycles=1, enable_cleanup=False,
                                           enable_inflate=False,
                                           thin_layer_policy="clamp_up_final"),
            needle_config=NeedleCycleConfig(boundary_needles=True),
            substrate=(G, "sub"))
