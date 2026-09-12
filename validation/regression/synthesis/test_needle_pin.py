# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""The needle-run fingerprint pin (ground rule 5 / R1, F0.1).

Until this file existed there was no in-tree test that pinned a bit-exact
needle run: the coverage was spread over ``validation/review/*.py`` and a
scratchpad harness outside the repo. This pin is committed BEFORE the
first line of ``DesignStack`` changes, so every later "fingerprint
unmoved" claim names a test that lives in the tree (F0.1, ground rule 5).

The run is deterministic by construction: fixed design, fixed targets,
fixed seed materials, no error groups (so `from_design` expands
deterministically), a deterministic P-function scan and deterministic LM.
The double-run assertion proves that at pin time; the recorded digest is
what every later item in the plan must hold.

Recorded at 0.6.32 (``72a2d4d``, release build, before F0.1).
"""

import hashlib

import numpy as np

from navette.spectralweave.target import TargetCollection, SpectralTarget, AngularTarget
from navette.synthesis.pipeline import run_needle

WL = np.linspace(450.0, 750.0, 31)
ANGS = [0.0, 30.0]

# air | H 100 | L 80 | H 150 | glass  -- constant-index, roughness-free.
_N_H = 2.35 + 0.01j
_N_L = 1.46 + 0.0j
_N_SUB = 1.52 + 0.0j

LAYERS = [
    (np.full(WL.shape, _N_H), 100.0),
    (np.full(WL.shape, _N_L), 80.0),
    (np.full(WL.shape, _N_H), 150.0),
]
NAMES = ["film0", "film1", "film2"]
CONTRAST = {
    "film0": np.full(WL.shape, _N_L),
    "film1": np.full(WL.shape, _N_H),
    "film2": np.full(WL.shape, _N_L),
}


def _targets():
    # Satisfiable mix: spectral exact + integral band + angular exact.
    tc = TargetCollection()
    tc.add(SpectralTarget(WL, np.full(31, 0.15), np.full(31, 0.02),
                          0.0, "s", "R", kind="e", weight=2.0))
    tc.add(SpectralTarget(WL, np.full(31, 0.85), np.full(31, 0.10),
                          0.0, "s", "T", kind="a", integral=True,
                          weight=3.0))
    tc.add(AngularTarget(600.0, np.array([0.0]), np.array([0.0]),
                         np.array([0.05]), "p", "R", kind="e",
                         normalize_count=True))
    return tc


def _cfg():
    from navette._smatrix import PipelineConfig
    return PipelineConfig(max_macro_cycles=3, needles_per_cycle=2,
                          enable_cleanup=True, enable_inflate=True,
                          stagnation_window=100)


def _canon(obj):
    """Recursively render floats bit-exactly (hex) for hashing."""
    if isinstance(obj, bool) or isinstance(obj, int):
        return obj
    if isinstance(obj, float):
        return obj.hex()
    if isinstance(obj, complex):
        return [obj.real.hex(), obj.imag.hex()]
    if isinstance(obj, dict):
        return {str(k): _canon(obj[k]) for k in sorted(obj, key=str)}
    if isinstance(obj, (list, tuple)):
        return [_canon(v) for v in obj]
    if isinstance(obj, str):
        return obj
    if obj is None:
        return None
    return repr(obj)


def _digest(res):
    payload = {
        "termination": res["termination"],
        "final_mf": res["final_mf"],
        "final_layer_count": res["final_layer_count"],
        "final_total_thickness_nm": res["final_total_thickness_nm"],
        "stagnation_detail": res["stagnation_detail"],
        "phases": res["phases"],
        "stack": res["stack"].to_dict(),
    }
    blob = repr(_canon(payload)).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


# Recorded on the 0.6.32 release build BEFORE F0.1 touched DesignStack.
# Placeholder on first run; the value below is what F0.1 and every later
# item must reproduce bit-exactly.
RECORDED_DIGEST = "cf753a910c7bf1d2e6c5ea5b9f5c696601ab675abaf497db28bac5340a63044d"


def test_needle_run_is_deterministic_and_matches_the_recorded_digest():
    res1 = run_needle(LAYERS, _targets(), ANGS, WL, CONTRAST,
                      pipeline_config=_cfg(), names=NAMES)
    res2 = run_needle(LAYERS, _targets(), ANGS, WL, CONTRAST,
                      pipeline_config=_cfg(), names=NAMES)
    d1, d2 = _digest(res1), _digest(res2)
    # Determinism first: two runs of the same inputs, one digest. If this
    # fails the pin is measuring noise, and the recorded value is void.
    assert d1 == d2, "needle run is not deterministic - pin is void"
    assert d1 == RECORDED_DIGEST, (
        f"needle fingerprint moved: recorded {RECORDED_DIGEST}, got {d1}"
    )


def test_needle_run_shape():
    # Sanity on the pin's problem: the run must actually exercise the
    # machinery the fingerprint is supposed to cover (insertions happened,
    # cleanup ran, more than one macro cycle).
    res = run_needle(LAYERS, _targets(), ANGS, WL, CONTRAST,
                     pipeline_config=_cfg(), names=NAMES)
    assert len(res["phases"]) == 3
    assert any(p["needle_results"] for p in res["phases"]), (
        "the pin must cover needle insertions - this design inserts nothing"
    )
    films = res["stack"].films()
    assert len(films) == res["final_layer_count"]
