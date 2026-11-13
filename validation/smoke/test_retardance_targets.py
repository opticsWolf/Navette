# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Retardance targets in synthesis (W8a): ``RetR``/``RetT``/``DeltaR``/``DeltaT``.

``Ret = arg(z_p · conj(z_s))`` of one element (r_f for R, t_f for T), the
engine's ``Retardance_R``/``Retardance_T``; ``DeltaR = RetR − π`` wrapped and
``DeltaT = RetT`` are the engine's ``Delta_R``/``Delta_T``. The per-point
comparison against the engine itself, the sensitivity and Jacobian against
finite differences, and the needle fold against the merit slope are cargo
tests; this file is the end-to-end part through the Python surface.
"""

import warnings

import numpy as np
import pytest

from navette._smatrix import LmConfig, SmatrixContext
from navette.materials import MaterialSpec
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec, sim_curves_from_arrays
from navette.synthesis.pipeline import stack_from_layers

_WL = np.linspace(480.0, 640.0, 9)


def _spec(label, wl, values, tol, angle_deg, phase=True, pol="u"):
    tc = TargetCollection()
    tc.add(SpectralTarget(wl, values, np.full(wl.size, tol), pol, label,
                          kind="e", phase=phase, angle_deg=angle_deg))
    return build_merit_spec(tc)


def _wrap(x):
    return x - 2.0 * np.pi * np.round(x / (2.0 * np.pi))


def test_retardance_labels_take_u_and_need_phase():
    with pytest.raises(ValueError, match='polarization="u"'):
        _spec("RetR", _WL, np.zeros(_WL.size), 1.0, 45.0, pol="p")
    with pytest.raises(ValueError, match="phase=True"):
        _spec("RetT", _WL, np.zeros(_WL.size), 1.0, 45.0, phase=False)


def test_the_forward_is_p_minus_s_and_delta_is_its_offset():
    # Hand rows: a zero target with tolerance one leaves the residual equal
    # to the demand's own value, wrapped.
    rng = np.random.default_rng(3)
    rows = {k: rng.uniform(0.2, 0.9, _WL.size) * np.exp(1j * rng.uniform(-np.pi, np.pi, _WL.size))
            for k in ("Rs", "Rp", "Ts", "Tp")}
    sim = sim_curves_from_arrays(_WL, {}, {k: v.reshape(1, -1) for k, v in rows.items()},
                                 angles_deg=[45.0])
    ret_r = np.angle(rows["Rp"] * np.conj(rows["Rs"]))
    ret_t = np.angle(rows["Tp"] * np.conj(rows["Ts"]))
    for label, want in [("RetR", ret_r), ("RetT", ret_t),
                        ("DeltaR", _wrap(ret_r - np.pi)), ("DeltaT", ret_t)]:
        got = np.asarray(_spec(label, _WL, np.zeros(_WL.size), 1.0, 45.0).residuals(sim))
        np.testing.assert_allclose(_wrap(got - want), 0.0, atol=1e-12, err_msg=label)


def test_a_demand_at_normal_incidence_only_warns():
    with pytest.warns(UserWarning, match="normal incidence"):
        _spec("RetR", _WL, np.zeros(_WL.size), 0.1, 0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _spec("RetR", _WL, np.zeros(_WL.size), 0.1, 45.0)


_INDICES = [2.30, 1.46, 2.30]
_REFERENCE = [118.0, 203.0, 64.0]
_START = [122.0, 198.0, 67.0]


def _layers(thicknesses):
    idx = [_INDICES[i % len(_INDICES)] if len(thicknesses) == 3
           else (2.30 if i % 2 == 0 else 1.46) for i in range(len(thicknesses))]
    return [(MaterialSpec("Konstant", dict(n=n)), d) for n, d in zip(idx, thicknesses)]


def _ctx(spec, wl, jacobian="analytic"):
    return SmatrixContext(spec, wl, 1e-9, 1000.0, LmConfig(jacobian=jacobian),
                          angles_deg=np.array([45.0]))


def _stack(ds, wl):
    st, _ = stack_from_layers(_layers(ds), wl, {}, names=[f"L{i}" for i in range(len(ds))])
    return st


@pytest.mark.parametrize("label", ["RetR", "DeltaT"])
def test_retardance_takes_the_analytic_path_to_the_reference(label):
    # Targets read off a reference stack put a merit-zero optimum a few nm
    # from the start; both Jacobian modes must find it.
    flat = _spec(label, _WL, np.zeros(_WL.size), 1.0, 45.0)
    ret = np.asarray(flat.residuals(_ctx(flat, _WL).simulate(_stack(_REFERENCE, _WL))))
    spec = _spec(label, _WL, ret, 0.05, 45.0)
    for jacobian in ("analytic", "fd"):
        st = _stack(_START, _WL)
        mf, report = _ctx(spec, _WL, jacobian).refine_thicknesses_report(st)
        assert (report["analytic_jacobians"] > 0) == (jacobian == "analytic")
        assert mf < 1e-10, (jacobian, mf)
        for got, want in zip([f["thickness"] for f in st.films()], _REFERENCE):
            assert got == pytest.approx(want, abs=1e-3)


def test_a_quarter_wave_retarder_in_reflection():
    # RetR = pi/2 at 45 deg, 545-555 nm, on H/L layers over glass: a
    # retardance-only synthesis from seeded random starts. The demand is
    # reachable (merit ~1e-26); plain LM finds it from some starts and not
    # others, so the test takes the best of a few.
    wl = np.linspace(545.0, 555.0, 3)
    spec = _spec("RetR", wl, np.full(wl.size, np.pi / 2), 0.05, 45.0)
    ctx = _ctx(spec, wl)
    rng = np.random.default_rng(7)
    best = (np.inf, None)
    for _ in range(12):
        st = _stack(list(rng.uniform(20.0, 250.0, 8)), wl)
        mf, _ = ctx.refine_thicknesses_report(st)
        if mf < best[0]:
            best = (mf, st)
    assert best[0] < 1e-12
    probe = _spec("DeltaR", wl, np.zeros(wl.size), 1.0, 45.0)
    delta = np.asarray(probe.residuals(_ctx(probe, wl).simulate(best[1])))
    np.testing.assert_allclose(delta, -np.pi / 2, atol=1e-6)


def test_a_broadband_quarter_wave_retarder_needs_the_needle():
    # The same demand on 500-600 nm: thickness refinement alone stalls at
    # merit 900-2700 from quarter-wave and random starts. Needle synthesis
    # reaches it, but which start does is chaotic: moving the H seed by
    # 1e-12 nm turns merit 0.019 into 5378, so the last-bit differences
    # between platform maths libraries pick a different path (Linux CI
    # stalled at 130 from the start that converged on Windows). Over the
    # 25 two-film starts below, 11 reach merit < 0.1 on Windows, spread
    # through the grid. The test takes the first start that converges and
    # re-scores its design with the engine's own Retardance_R.
    from navette.smatrix.smatrix import Request, ScatterMatrix
    from navette.synthesis.pipeline import PipelineConfig, run_needle

    wl = np.linspace(500.0, 600.0, 21)
    h, l, glass = (np.full(wl.size, n + 0j) for n in (2.30, 1.46, 1.52))
    tc = TargetCollection()
    tc.add(SpectralTarget(wl, np.full(wl.size, np.pi / 2), np.full(wl.size, 0.05),
                          "u", "RetR", kind="e", phase=True, angle_deg=45.0))
    cfg = PipelineConfig(max_macro_cycles=60, needles_per_cycle=2, max_film_layers=40,
                         max_total_thickness_nm=6000.0)
    tried = []
    for dh in (300.0, 310.0, 320.0, 280.0, 290.0):
        for dl in (300.0, 310.0, 320.0, 280.0, 290.0):
            out = run_needle(layers=[(h, dh), (l, dl)], names=["H0", "L0"], targets=tc,
                             angles_deg=[45.0], wavelengths=wl,
                             materials={"H": h, "L": l}, pipeline_config=cfg,
                             substrate=(glass, "glass"))
            tried.append(out["final_mf"])
            if out["final_mf"] < 0.1:
                break
        else:
            continue
        break
    assert out["final_mf"] < 0.1, tried
    films = out["stack"].to_dict()["films"]
    rows = [np.ones(wl.size, complex)] + [np.asarray(f["nk"]) for f in films] + [glass]
    d = [0.0] + [f["thickness"] for f in films] + [0.0]
    sm = ScatterMatrix(np.vstack(rows), d, wavelengths=wl, angles_deg=[45.0])
    ret = np.asarray(sm.compute(Request.RETARD_R)["Retardance_R"]).ravel()
    # Merit < 0.1 at tolerance 0.05 bounds every point by 0.05 * sqrt(0.1).
    assert np.max(np.abs(ret - np.pi / 2)) < 0.05 * np.sqrt(0.1)
    assert np.sum(((ret - np.pi / 2) / 0.05) ** 2) == pytest.approx(out["final_mf"], rel=1e-6)
