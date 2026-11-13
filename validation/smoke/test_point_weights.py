# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Per-point weights in synthesis (W4): ``weights`` / ``weight_profile``.

A pointwise residual is scaled by ``√w_i``; an integral target takes the
weighted means ``Σw d/Σw`` over ``Σw tol/Σw``. The merit, sensitivity and
needle fold against hand formulas and finite differences are cargo tests;
this file is the end-to-end part through the Python surface: the three
ways to give weights, the refusals, the profile door, and a weighted
needle run re-scored by the engine.
"""

import numpy as np
import pytest

from navette.spectralweave import TargetCollection, WeightProfile
from navette.spectralweave.target import AngularTarget, SpectralTarget
from navette.synthesis import build_merit_spec, sim_curves_from_arrays

_WL = np.linspace(450.0, 650.0, 11)
_RS = np.linspace(0.10, 0.30, _WL.size)
_SIM = sim_curves_from_arrays(_WL, {"Rs": _RS.reshape(1, -1)}, angles_deg=[0.0])


def _residuals(*targets, profiles=None):
    tc = TargetCollection()
    for t in targets:
        tc.add(t)
    for name, p in (profiles or {}).items():
        tc.weight_profiles[name] = p
    return np.asarray(build_merit_spec(tc).residuals(_SIM))


def _rs(tol=0.01, **kw):
    return SpectralTarget(_WL, np.full(_WL.size, 0.2), np.broadcast_to(tol, _WL.shape).copy(),
                          "s", "R", angle_deg=0.0, **kw)


def test_a_point_weight_is_a_tolerance_rescaling():
    w = np.linspace(0.5, 3.0, _WL.size)
    np.testing.assert_allclose(_residuals(_rs(weights=w)),
                               _residuals(_rs(tol=0.01 / np.sqrt(w))), rtol=1e-13)


def test_ones_are_no_weights():
    plain = _residuals(_rs())
    np.testing.assert_array_equal(_residuals(_rs(weights=np.ones(_WL.size))), plain)
    np.testing.assert_array_equal(_residuals(_rs(weight_profile=WeightProfile.constant(1.0))), plain)


def test_explicit_inline_and_named_weights_agree():
    p = WeightProfile.product(WeightProfile.ramp(450.0, 650.0, 0.5, 2.0),
                              WeightProfile.gaussian(550.0, 80.0, 2.0, base=1.0))
    explicit = _residuals(_rs(weights=p.evaluate(_WL)))
    np.testing.assert_array_equal(_residuals(_rs(weight_profile=p)), explicit)
    np.testing.assert_array_equal(
        _residuals(_rs(weight_profile="emph"), profiles={"emph": p}), explicit)


def test_an_integral_takes_the_weighted_mean():
    # Against the unweighted pointwise residuals, which carry d_i/tol times
    # the target's normalization scale.
    w = np.linspace(1.0, 4.0, _WL.size)
    pointwise = _residuals(_rs())
    got = _residuals(_rs(weights=w, integral=True))
    assert got.shape == (1,)
    assert got[0] == pytest.approx(np.sum(w * pointwise) / np.sum(w), rel=1e-12)


def test_angular_weights_take_degrees():
    angles = np.array([0.0, 20.0, 40.0])
    sim = sim_curves_from_arrays([550.0], {"Rs": np.full((3, 1), 0.3)}, angles_deg=angles)
    def res(**kw):
        tc = TargetCollection()
        tc.add(AngularTarget(550.0, np.full(3, 0.2), np.full(3, 0.1), "s", "R",
                             angles_deg=angles, **kw))
        return np.sort(np.asarray(build_merit_spec(tc).residuals(sim)))

    p = WeightProfile.table([0.0, 40.0], [1.0, 4.0], interp="linear")
    np.testing.assert_allclose(res(weight_profile=p), res() * np.sqrt([1.0, 2.5, 4.0]),
                               rtol=1e-12)


def test_profiles_evaluate_and_round_trip():
    b = WeightProfile.bands([(500.0, 550.0, 3.0), (550.0, 600.0, 2.0)], outside=0.5)
    np.testing.assert_array_equal(b.evaluate(np.array([450.0, 520.0, 580.0, 620.0])),
                                  [0.5, 3.0, 2.0, 0.5])
    assert WeightProfile.from_json(b.to_json()) == b
    t = WeightProfile.table([400.0, 500.0, 600.0], [1.0, 3.0, 1.0])
    v = t.evaluate(np.linspace(350.0, 650.0, 61))
    assert v.min() >= 1.0 and v.max() <= 3.0   # pchip stays inside the data; clamp outside


@pytest.mark.parametrize("make, match", [
    (lambda: WeightProfile.ramp(600.0, 500.0, 1.0, 2.0), "lo"),
    (lambda: WeightProfile.bands([(500.0, 560.0, 1.0), (550.0, 600.0, 1.0)]), "overlap"),
    (lambda: WeightProfile.table([400.0, 500.0], [1.0, 2.0], interp="cubic"), "interp"),
    (lambda: WeightProfile.from_json('{"kind": "constant", "value": 1, "x": 2}'), "unknown field"),
])
def test_bad_profiles_refuse(make, match):
    with pytest.raises(ValueError, match=match):
        make()


@pytest.mark.parametrize("kw, profiles, match", [
    ({"weights": np.ones(_WL.size), "weight_profile": WeightProfile.constant(1.0)}, None, "not both"),
    ({"weights": np.ones(_WL.size - 1)}, None, "shape mismatch"),
    ({"weights": np.zeros(_WL.size)}, None, "all zero"),
    ({"weights": np.r_[-1.0, np.ones(_WL.size - 1)]}, None, ">= 0|finite"),
    ({"weight_profile": "missing"}, {"emph": WeightProfile.constant(2.0)}, "not in weight_profiles"),
])
def test_bad_weights_refuse(kw, profiles, match):
    with pytest.raises(ValueError, match=match):
        _residuals(_rs(**kw), profiles=profiles)


def test_weights_on_an_integral_angular_target_refuse():
    tc = TargetCollection()
    tc.add(AngularTarget(550.0, np.full(3, 0.2), np.full(3, 0.1), "s", "R", integral=True,
                         weights=np.array([1.0, 2.0, 3.0]), angles_deg=np.array([0.0, 20.0, 40.0])))
    with pytest.raises(ValueError, match="integral angular"):
        build_merit_spec(tc)


def test_the_weaver_path_refuses_weights():
    tc = TargetCollection()
    tc.add(_rs(weights=np.ones(_WL.size)))
    with pytest.raises(ValueError, match="build_merit_spec"):
        tc.build_weaver()


def test_a_weighted_needle_run_is_scored_with_its_weights():
    # The 500-600 nm quarter-wave retarder with the band centre at double
    # weight (merit 0.002, 14 films); the run's merit is the engine's
    # weighted sum. A 4x centre or a 1-2 ramp stalls at 7 films from this
    # seed (163, 18): the landscape moves with the weights as it does with
    # the seed, and the weighted fold is the merit slope (cargo oracle).
    from navette.smatrix.smatrix import Request, ScatterMatrix
    from navette.synthesis.pipeline import PipelineConfig, run_needle

    wl = np.linspace(500.0, 600.0, 21)
    h, l, glass = (np.full(wl.size, n + 0j) for n in (2.30, 1.46, 1.52))
    p = WeightProfile.bands([(530.0, 570.0, 2.0)], outside=1.0)
    w = p.evaluate(wl)
    tc = TargetCollection()
    tc.add(SpectralTarget(wl, np.full(wl.size, np.pi / 2), np.full(wl.size, 0.05),
                          "u", "RetR", kind="e", phase=True, angle_deg=45.0, weight_profile=p))
    cfg = PipelineConfig(max_macro_cycles=60, needles_per_cycle=2, max_film_layers=40,
                         max_total_thickness_nm=6000.0)
    out = run_needle(layers=[(h, 600.0), (l, 600.0)], names=["H0", "L0"], targets=tc,
                     angles_deg=[45.0], wavelengths=wl, materials={"H": h, "L": l},
                     pipeline_config=cfg, substrate=(glass, "glass"))
    films = out["stack"].to_dict()["films"]
    rows = [np.ones(wl.size, complex)] + [np.asarray(f["nk"]) for f in films] + [glass]
    d = [0.0] + [f["thickness"] for f in films] + [0.0]
    sm = ScatterMatrix(np.vstack(rows), d, wavelengths=wl, angles_deg=[45.0])
    ret = np.asarray(sm.compute(Request.RETARD_R)["Retardance_R"]).ravel()
    assert np.sum(w * ((ret - np.pi / 2) / 0.05) ** 2) == pytest.approx(out["final_mf"], rel=1e-6)
    assert out["final_mf"] < 0.05, out["final_mf"]
