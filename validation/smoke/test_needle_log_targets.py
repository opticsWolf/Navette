# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Log-normalized targets in needle synthesis (W5).

A Log demand folds the tangent of ``log10(s)·nf − target`` at the
operating point, so its needle gradient is the merit's (the fold against
the merit slope, for every kind, the integral and the back side, is the
cargo oracle ``evaluator::needle_oracle``). Before W5 the needle refused
a Log demand mid-run, and ``"auto"`` resolves to Log for any target with
a 100x range. This is the end-to-end part through the Python surface.
"""

import numpy as np
import pytest

from navette.smatrix.smatrix import Request, ScatterMatrix
from navette.spectralweave import TargetCollection
from navette.spectralweave.target import SpectralTarget
from navette.synthesis import build_merit_spec, sim_curves_from_arrays
from navette.synthesis.pipeline import PipelineConfig, run_needle

_WL = np.linspace(450.0, 650.0, 21)
_H, _L, _GLASS = (np.full(_WL.size, n + 0j) for n in (2.30, 1.46, 1.52))


def _engine_merit(tc, layers):
    rows = [np.ones(_WL.size, complex)] + [nk for nk, _ in layers] + [_GLASS]
    d = [0.0] + [t for _, t in layers] + [0.0]
    res = ScatterMatrix(np.vstack(rows), d, wavelengths=_WL, angles_deg=[0.0]).compute(
        Request.RS | Request.RP)
    rs, rp = (np.asarray(res[k]).reshape(1, -1) for k in ("Rs", "Rp"))
    sim = sim_curves_from_arrays(_WL, {"Rs": rs, "Rp": rp, "Ru": 0.5 * (rs + rp)},
                                 angles_deg=[0.0])
    return build_merit_spec(tc).merit(sim, 1e6)


def test_a_log_ar_target_runs_the_needle():
    # R = 1e-4 over 450-650 nm, 12 % tolerance in log10 units. From
    # 100 nm H + 100 nm L the needle reaches 9 films (R ~ 10^-3.6 in band,
    # merit 177, a plateau of these two materials), not a refusal.
    tc = TargetCollection()
    tc.add(SpectralTarget(_WL, np.full(_WL.size, 1e-4), np.full(_WL.size, 0.05), "u", "R",
                          kind="e", normalization_mode="log", angle_deg=0.0))
    seed = [(_H, 100.0), (_L, 100.0)]
    cfg = PipelineConfig(max_macro_cycles=30, needles_per_cycle=2, max_film_layers=20,
                         max_total_thickness_nm=2000.0)
    out = run_needle(layers=seed, names=["H0", "L0"], targets=tc, angles_deg=[0.0],
                     wavelengths=_WL, materials={"H": _H, "L": _L}, pipeline_config=cfg,
                     substrate=(_GLASS, "glass"))
    films = out["stack"].to_dict()["films"]
    assert len(films) > len(seed)
    final = _engine_merit(tc, [(np.asarray(f["nk"]), f["thickness"]) for f in films])
    assert final == pytest.approx(out["final_mf"], rel=1e-6)
    assert out["final_mf"] < 0.1 * _engine_merit(tc, seed)
