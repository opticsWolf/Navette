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


def _joint(phase_in=()):
    # One coat, two surroundings in air: alone, and on both faces of a
    # 1 mm plate (the second face a copy of the first).
    from navette._smatrix import PipelineConfig
    plate = {"material": G, "thickness_nm": 1.0e6, "name": "plate",
             "coherent": False}
    envs = [
        {"name": "alone", "stack": [{"design": "coat"}]},
        {"name": "plate", "stack": [{"design": "coat"}, {"layers": [plate]},
                                    {"design": "coat"}]},
    ]
    tc = TargetCollection()
    for env in ("alone", "plate"):
        tc.add(SpectralTarget(WL, np.zeros(31), np.full(31, 0.01), "s", "R",
                              kind="e", angle_deg=0.0, environment=env))
    for env in phase_in:
        tc.add(SpectralTarget(WL, np.zeros(31), np.full(31, 0.1), "s", "T",
                              kind="e", angle_deg=0.0, phase=True,
                              environment=env))
    return run_needle(
        targets=tc, angles_deg=[0.0], wavelengths=WL, materials={"H": H, "L": L},
        design={"coat": [(H, 100.0, "H"), (L, 80.0, "L")]}, environments=envs,
        pipeline_config=PipelineConfig(max_macro_cycles=2, needles_per_cycle=1,
                                       enable_cleanup=False, enable_inflate=False,
                                       thin_layer_policy="clamp_up_final"),
        substrate=(1.0, "air"))


def test_phase_demand_is_dropped_only_in_the_plate_environment():
    # D4: the plate environment's phase demand goes, with a warning; the
    # coherent environment's stays and the run completes.
    with pytest.warns(UserWarning, match="phase demand.*'plate'"):
        res = _joint(phase_in=("alone", "plate"))
    assert res["phases"], res["termination"]


def test_phase_demand_only_behind_a_plate_is_refused():
    with pytest.raises(Exception, match="every phase demand"):
        _joint(phase_in=("plate",))


def test_joint_run_with_a_plate_in_one_environment():
    res = _joint()
    assert any(p["needle_results"] for p in res["phases"])
    assert res["final_mf"] < res["phases"][0]["mf_after_needle"] * 1.0001


def test_an_environment_keeps_its_own_exit_medium():
    # The run's substrate is the default exit; an environment that ends in
    # its own layer_type 2 row keeps that one instead of getting two.
    from navette._smatrix import PipelineConfig

    def run(exit_n):
        envs = [
            {"name": "in_air", "stack": [{"design": "coat"}]},
            {"name": "on_glass", "stack": [
                {"design": "coat"},
                {"layers": [{"material": exit_n, "thickness_nm": 0.0,
                             "name": "exit", "layer_type": 2}]}]},
        ]
        tc = TargetCollection()
        for env in ("in_air", "on_glass"):
            tc.add(SpectralTarget(WL, np.zeros(31), np.full(31, 0.01), "s", "R",
                                  kind="e", angle_deg=0.0, environment=env))
        return run_needle(
            targets=tc, angles_deg=[0.0], wavelengths=WL,
            materials={"H": H, "L": L},
            design={"coat": [(H, 100.0, "H"), (L, 80.0, "L")]}, environments=envs,
            pipeline_config=PipelineConfig(max_macro_cycles=1, needles_per_cycle=0,
                                           enable_cleanup=False, enable_inflate=False,
                                           thin_layer_policy="clamp_up_final"),
            substrate=(1.0, "air"))

    # The environment's own exit is the one solved against: glass behind
    # the coat scores differently from air behind it.
    assert run(1.52)["final_mf"] != run(1.0)["final_mf"]


def test_phase_demand_through_a_plate_is_refused():
    with pytest.raises(Exception, match="phase"):
        _run(_targets(phase=True))
