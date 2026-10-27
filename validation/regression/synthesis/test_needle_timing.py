# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""`NeedleCycleConfig(record_timing=True)` observes the run and never steers it.

The timers exist for `validation/benches/synthesis/bench_needle_insertion.py`.
Two promises are pinned here, through the public `run_needle` surface:

* a timed run is the untimed run, bit for bit: same final stack, same merit
  trajectory, same insertion sites;
* an untimed record carries no `timing` key at all, so a caller that never
  asked for timings sees exactly the dict it saw before the timers existed.
"""

import numpy as np

from navette._smatrix import NeedleCycleConfig, PipelineConfig
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis.pipeline import run_needle

WL = np.linspace(450.0, 650.0, 11)
STEPS = ("refold", "scan", "merit", "lm", "initial_lm")


def _run(record_timing):
    hi, lo = np.full(WL.shape, 2.35 + 0j), np.full(WL.shape, 1.46 + 0j)
    tc = TargetCollection()
    tc.add(SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, 0.005),
                          0.0, "s", "R"))
    cfg = PipelineConfig(max_macro_cycles=3, needles_per_cycle=2)
    nc = NeedleCycleConfig(record_timing=record_timing)
    return run_needle([(lo, 200.0), (hi, 100.0)], tc, [0.0], WL,
                      {0: hi, 1: lo}, pipeline_config=cfg, needle_config=nc,
                      substrate=(np.full(WL.shape, 1.52 + 0j), "sub"))


def _records(res):
    return [r for p in res["phases"] for r in p["needle_results"]]


def _bits(x):
    return None if x is None else np.float64(x).tobytes()


def test_timed_run_is_the_untimed_run():
    off, on = _run(False), _run(True)
    assert _bits(off["final_mf"]) == _bits(on["final_mf"])
    thick = lambda r: [_bits(f["thickness"]) for f in r["stack"].films()]
    assert thick(off) == thick(on)
    a, b = _records(off), _records(on)
    assert any(r["insertion"] is not None for r in a), "no insertion: degenerate"
    assert len(a) == len(b)
    for ra, rb in zip(a, b):
        assert {k: v for k, v in rb.items() if k != "timing"}.keys() == ra.keys()
        for key in ("merit_before", "merit_after", "best_p",
                    "predicted_improvement"):
            assert _bits(ra[key]) == _bits(rb[key]), key
        assert ra["insertion"] == rb["insertion"]


def test_timing_key_only_when_asked():
    for r in _records(_run(False)):
        assert "timing" not in r
    for p in _run(True)["phases"]:
        recs = p["needle_results"]
        for i, r in enumerate(recs):
            t = r["timing"]
            assert set(t) == set(STEPS)
            assert all(np.isfinite(v) and v >= 0.0 for v in t.values()), t
            if r["insertion"] is not None:
                assert t["scan"] > 0.0 and t["lm"] > 0.0, t
            # The pass's opening optimization rides on its first record only.
            assert (t["initial_lm"] > 0.0) == (i == 0), (i, t)
