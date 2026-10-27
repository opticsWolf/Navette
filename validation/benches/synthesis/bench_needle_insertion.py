#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Benchmark: where a needle insertion spends its time, against film count.

`bench_refold.py` times the per-cycle re-fold on a 3-wavelength, one-film
toy, and says nothing about the insertion scan — the sweep of the needle
gradient over every candidate depth, which is the part of synthesis that
grows with the design. This bench grows a broadband AR coating from a
4-film seed to ~26 films (61 wavelengths, 400-700 nm, normal incidence)
under `NeedleCycleConfig(record_timing=True)` and reads the engine's own
per-insertion timers:

    refold   simulate + fold against the live operating point
    scan     scan-site build + one analytic sweep per host material
    merit    the convergence test's evaluate_merit
    lm       the thickness re-optimization after the insertion

The run is deterministic, so it is repeated and each insertion keeps its
fastest repeat; the trajectories are asserted identical first. The timers
are opt-in and observe only (`test_needle_timing.py` pins that).

Run explicitly:  python validation/benches/synthesis/bench_needle_insertion.py
"""

# --- bench preamble (R2.1/R2.2) ------------------------------------------
# UTF-8 console, repo `src/` on sys.path, and a hard gate against timing a
# debug build. Must precede any `navette` import.
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from _bench_common import require_release, setup_bench  # noqa: E402

setup_bench()
require_release()
# -------------------------------------------------------------------------

import sys
import time

import numpy as np

from navette._smatrix import NeedleCycleConfig, PipelineConfig
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis.pipeline import run_needle

WL = np.linspace(400.0, 700.0, 61)
H, L, GLASS = 2.35, 1.46, 1.52
TOL = 0.002
STEPS = ("refold", "scan", "merit", "lm")
BUCKETS = ((1, 10), (11, 20), (21, 30), (31, 10**6))
REPEATS = 3

FAILURES = []


def check(name, ok, detail=""):
    print(f"  {name}: {'OK' if ok else 'FAIL'} {detail}")
    if not ok:
        FAILURES.append(name)


def run_once():
    """One timed synthesis: (wall seconds, result, insertion records)."""
    hi, lo = np.full(WL.shape, H + 0j), np.full(WL.shape, L + 0j)
    tc = TargetCollection()
    tc.add(SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, TOL),
                          0.0, "s", "R"))
    # A deliberately poor seed; each film may host the other material.
    seed = [(lo, 250.0), (hi, 150.0), (lo, 250.0), (hi, 150.0)]
    # The site's AR problem, with the stagnation stop off: left on, it ends
    # the run at ~22 films once the merit plateaus, and the bench wants
    # the insertions past that point too. 30 cycles stay the bound.
    cfg = PipelineConfig(max_macro_cycles=30, needles_per_cycle=2,
                         max_film_layers=32, stagnation_window=100)
    nc = NeedleCycleConfig(record_timing=True)
    t0 = time.perf_counter()
    res = run_needle(seed, tc, [0.0], WL, {0: hi, 1: lo, 2: hi, 3: lo},
                     pipeline_config=cfg, needle_config=nc,
                     substrate=(np.full(WL.shape, GLASS + 0j), "sub"))
    wall = time.perf_counter() - t0
    recs = [r for p in res["phases"] for r in p["needle_results"]]
    return wall, res, recs


def main():
    print(f"--- needle insertion: {WL.size} wavelengths, 4-film seed, "
          f"best of {REPEATS} ---")
    runs = [run_once() for _ in range(REPEATS)]
    wall = min(r[0] for r in runs)
    res, recs = runs[0][1], runs[0][2]

    # Same trajectory every repeat, or the per-insertion minimum compares
    # different insertions.
    shape = lambda rs: [(r["layer_count"], r["insertion"]) for r in rs]
    check("deterministic_trajectory",
          all(shape(r[2]) == shape(recs) for r in runs[1:]))

    # Per record, the fastest repeat of each step.
    best = [{s: min(run[2][i]["timing"][s] for run in runs) for s in STEPS}
            for i in range(len(recs))]
    initial_lm = sum(min(run[2][i]["timing"]["initial_lm"] for run in runs)
                     for i in range(len(recs)))
    inserted = [(r, b) for r, b in zip(recs, best) if r["insertion"] is not None]
    # `layer_count` is recorded after the insertion; the scan saw one less
    # film-split, i.e. two fewer films.
    films_scanned = lambda r: r["layer_count"] - 2

    print(f"  {'films':>5}  {'scan ms':>9}  {'lm ms':>9}  {'refold ms':>9}  "
          f"{'merit ms':>9}  {'scan %':>6}")
    for r, b in inserted:
        tot = sum(b.values())
        print(f"  {films_scanned(r):5d}  {b['scan'] * 1e3:9.3f}  "
              f"{b['lm'] * 1e3:9.3f}  {b['refold'] * 1e3:9.3f}  "
              f"{b['merit'] * 1e3:9.3f}  {b['scan'] / tot * 100:5.1f}%")

    print("--- by film count at the scan (mean per insertion) ---")
    print(f"  {'films':>7}  {'n':>3}  {'scan ms':>9}  {'lm ms':>9}  "
          f"{'scan/lm':>7}")
    for lo_f, hi_f in BUCKETS:
        sel = [b for r, b in inserted if lo_f <= films_scanned(r) <= hi_f]
        if not sel:
            continue
        scan = np.mean([b["scan"] for b in sel])
        lm = np.mean([b["lm"] for b in sel])
        label = f"{lo_f}-{hi_f}" if hi_f < 10**6 else f"{lo_f}+"
        print(f"  {label:>7}  {len(sel):3d}  {scan * 1e3:9.3f}  "
              f"{lm * 1e3:9.3f}  {scan / lm:7.2f}")

    totals = {s: sum(b[s] for b in best) for s in STEPS}
    timed = sum(totals.values()) + initial_lm
    print("--- whole run ---")
    print(f"  wall:            {wall * 1e3:9.1f} ms  "
          f"({len(inserted)} insertions, {res['final_layer_count']} films, "
          f"mf {res['final_mf']:.4g}, {res['termination']})")
    for s in STEPS:
        print(f"  {s:<15}  {totals[s] * 1e3:9.1f} ms  "
              f"{totals[s] / wall * 100:5.1f}% of wall")
    print(f"  {'initial lm':<15}  {initial_lm * 1e3:9.1f} ms  "
          f"{initial_lm / wall * 100:5.1f}% of wall")
    print(f"  untimed rest     {(wall - timed) * 1e3:9.1f} ms  "
          f"{(wall - timed) / wall * 100:5.1f}% of wall  "
          "(cleanup, inflate, final optimize, Python glue)")

    counts = [films_scanned(r) for r, _ in inserted]
    check("reaches_10_films", bool(counts) and max(counts) >= 10,
          f"max {max(counts, default=0)} films scanned")
    check("insertions_timed",
          all(b["scan"] > 0.0 and b["lm"] > 0.0 for _, b in inserted))
    # Each repeat's timed steps are a subset of its wall; the per-step
    # minima sum to less than the fastest wall a fortiori.
    check("timed_within_wall", timed <= wall,
          f"{timed * 1e3:.1f} ms of {wall * 1e3:.1f} ms")

    print("ALL OK" if not FAILURES else f"MISMATCH {FAILURES}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
