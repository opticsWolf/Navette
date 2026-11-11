# SPDX-License-Identifier: LGPL-3.0-or-later
"""W6: synthesis through an incoherent plate (coatings on both faces)."""
import numpy as np
import pytest

from navette.spectralweave.target import TargetCollection, SpectralTarget
from navette.synthesis.pipeline import run_needle

WL = np.linspace(450.0, 750.0, 31)
H = np.full(WL.shape, 2.35 + 0j)
L = np.full(WL.shape, 1.46 + 0j)
G = np.full(WL.shape, 1.52 + 0j)
PLATE = {"coherent": False, "optimize": False, "needle": False}


def _targets(phase=False):
    tc = TargetCollection()
    tc.add(SpectralTarget(WL, np.zeros(31), np.full(31, 0.01), "s", "R",
                          kind="e", angle_deg=0.0))
    if phase:
        tc.add(SpectralTarget(WL, np.zeros(31), np.full(31, 0.1), "s", "R",
                              kind="e", angle_deg=0.0, phase=True))
    return tc


def _run(targets):
    from navette._smatrix import PipelineConfig
    return run_needle(
        [(H, 100.0), (L, 80.0), (G, 1.0e6), (L, 90.0)], targets, WL,
        {"H": H, "L": L},
        pipeline_config=PipelineConfig(max_macro_cycles=2, needles_per_cycle=1),
        names=["H", "L", "plate", "Lb"], angles_deg=[0.0],
        substrate=(1.0, "air"), per_film_flags={"plate": PLATE})


def test_needle_run_coats_both_faces_of_a_fixed_plate():
    res = _run(_targets())
    assert res["phases"], res["termination"]
    films = res["stack"].films()
    plate = [f for f in films if not f["coherent"]]
    assert len(plate) == 1 and plate[0]["thickness"] == 1.0e6
    # The budget and the report count coating only.
    coating = sum(f["thickness"] for f in films if f["coherent"])
    assert res["final_total_thickness_nm"] == pytest.approx(coating)
    # Both faces stay coated, a needle went in, and the merit only falls.
    i = films.index(plate[0])
    assert 0 < i < len(films) - 1
    assert any(p["needle_results"] for p in res["phases"])
    assert res["final_mf"] <= res["phases"][0]["mf_after_needle"]


def test_phase_demand_through_a_plate_is_refused():
    with pytest.raises(Exception, match="phase"):
        _run(_targets(phase=True))
