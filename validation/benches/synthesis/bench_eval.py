#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Benchmark: the design pipeline's inner eval, assemble -> simulate -> merit.

The multi-environment phase (F2.2) turns one solve per eval into K, behind
one branch taken outside every loop. The gate on that branch is that K=1
costs what the pre-segment path cost -- threshold <1%
(``docs/plans/multi_environment_plan.md`` §4.6 (b)). A gate needs a
baseline, and a baseline taken *after* the change is not a baseline, so
F2.1 records one.

It also needed a harness. Nothing under ``validation/benches`` timed this:
``bench_refold.py`` is the closest and measures re-fold against static fold
with an LM cost breakdown; ``smatrix/`` times the solver core and the
backside path; ``structure/`` times a grid assert. So this file is F2.1's
deliverable as much as the segment schema is (plan §4, F2.1, AMENDED A7).

Three phases are timed separately, because F2.2 moves exactly one of them:

``assemble``
    ``stack_from_layers`` -- run start and after each needle insertion, NOT
    per eval. F2.2 makes this K concatenations.
``simulate``
    one solve. F2.2 makes this K solves; §6 limitation 6 says the ×K cost
    is inherent, not overhead.
``merit``
    residuals + merit over the demand set. F2.2 routes each demand to its
    environment's curves; at K=1 it must be the same call sequence.

``eval`` does NOT equal ``assemble + simulate + merit_only``, and this is a
property of the measurement rather than of the pipeline. Assembling a stack
and immediately solving it, in one loop that frees it again, measures about
1.7x the phases measured separately on this machine -- allocate/solve/free
in a tight loop is a different memory pattern from solving a stack that was
built before the loop started. Both numbers are honest; neither is the
other's sum. F2.2 must compare each phase against the SAME phase, never a
phase against a sum.

Run explicitly:  python validation/benches/synthesis/bench_eval.py
Write the committed baseline:  ... bench_eval.py --write
"""

# --- bench preamble (R2.1/R2.2) ------------------------------------------
# UTF-8 console, repo `src/` on sys.path, and a hard gate against timing a
# debug build. Must precede any `navette` import.
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from _bench_common import bench_provenance, require_release, setup_bench  # noqa: E402

setup_bench()
require_release()
# -------------------------------------------------------------------------

import argparse
import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from navette._smatrix import SmatrixContext
from navette.materials import MaterialSpec
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec
from navette.synthesis.pipeline import stack_from_layers

RESULTS = Path(__file__).resolve().parent / "results" / "eval_baseline.json"

# One fixed problem, stated here and nowhere else. It is deliberately a
# plain multilayer on a plain grid: the baseline exists to measure the
# K=1 branch, not the physics, and a problem with graded spans or color
# demands would make the number move for reasons F2.2 does not control.
WL = np.linspace(400.0, 700.0, 61)
ANGLES = np.array([0.0, 30.0])
N_H, N_L = 2.35, 1.46
D_H, D_L = 95.0, 160.0
N_PAIRS = 6


def problem():
    """(layers, contrast, targets) for the fixed baseline design."""
    layers = []
    for i in range(N_PAIRS):
        layers.append((MaterialSpec("Konstant", dict(n=N_H)), D_H))
        layers.append((MaterialSpec("Konstant", dict(n=N_L)), D_L))
    names = [f"film{i}" for i in range(len(layers))]
    contrast = {name: MaterialSpec("Konstant", dict(n=N_H if i % 2 else N_L))
                for i, name in enumerate(names)}
    tc = TargetCollection()
    n_demands = 0
    for ang in ANGLES:
        tc.add(SpectralTarget(WL, np.zeros(WL.size), np.full(WL.size, 0.01),
                              float(ang), "s", "R", kind="e"))
        tc.add(SpectralTarget(WL, np.ones(WL.size), np.full(WL.size, 0.01),
                              float(ang), "p", "T", kind="e"))
        n_demands += 2
    return layers, names, contrast, tc, n_demands


def timed(fn, repeats, inner):
    """Best-of-`repeats` mean seconds per call.

    Best-of, not mean-of: on a desktop the distribution's tail is other
    processes, and the <1% gate has to compare the machine's floor, not
    its mood.

    Every phase runs the same `inner` count. An earlier draft ran the
    allocating phases (assemble, eval) at a tenth of it and reported an
    `eval` 1.8x the sum of its parts -- allocation churn that had not
    amortized, not a cost. A phase measured on a different number of
    iterations than the phase it will be compared against is not
    comparable to it.
    """
    best = float("inf")
    for _ in range(repeats):
        t0 = time.perf_counter()
        for _ in range(inner):
            fn()
        dt = (time.perf_counter() - t0) / inner
        best = min(best, dt)
    return best


def measure(repeats=5, inner=50):
    layers, names, contrast, tc, n_demands = problem()
    spec = build_merit_spec(tc)

    stack, cmap = stack_from_layers(layers, WL, contrast, names=names)
    ctx = SmatrixContext(spec, ANGLES, WL)
    sim = ctx.simulate(stack)
    merit = float(spec.merit(sim, 1e6))
    n_res = len(np.asarray(spec.residuals(sim)).ravel())

    t_assemble = timed(
        lambda: stack_from_layers(layers, WL, contrast, names=names),
        repeats, inner)
    t_simulate = timed(lambda: ctx.simulate(stack), repeats, inner)
    t_merit = timed(lambda: spec.merit(ctx.simulate(stack), 1e6), repeats, inner)

    def full_eval():
        st, _ = stack_from_layers(layers, WL, contrast, names=names)
        return spec.merit(ctx.simulate(st), 1e6)

    t_eval = timed(full_eval, repeats, inner)

    return {
        "problem": {
            "n_films": len(layers),
            "n_wavelengths": int(WL.size),
            "n_angles": int(ANGLES.size),
            "n_demands": n_demands,
            "n_residuals": n_res,
            "n_contrast_seeds": len(cmap),
        },
        "merit": merit,
        "seconds": {
            "assemble": t_assemble,
            "simulate": t_simulate,
            # merit-only is the difference: `t_merit` times simulate+merit,
            # because the demand set is evaluated against a fresh solve.
            "merit_only": max(0.0, t_merit - t_simulate),
            "simulate_merit": t_merit,
            # The plan's phrase, measured whole: assemble + simulate +
            # merit. Not the per-eval cost (assembly is run-start and
            # post-insertion only, 4.6) -- it is here because F2.2 moves
            # assembly as well as the solve, and a gate that watched only
            # the solve would not see a concatenation that got expensive.
            "eval": t_eval,
        },
        "timing": {"repeats": repeats, "inner": inner, "reduction": "best-of"},
    }


def git_hash():
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                             text=True, check=True,
                             cwd=str(Path(__file__).resolve().parents[3]))
        return out.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true",
                    help="write results/eval_baseline.json (F2.2's gate reads it)")
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--inner", type=int, default=50)
    args = ap.parse_args()

    res = measure(args.repeats, args.inner)
    res["provenance"] = bench_provenance()
    res["provenance"]["commit"] = git_hash()
    res["provenance"]["recorded_utc"] = datetime.now(timezone.utc).isoformat(
        timespec="seconds")
    res["note"] = (
        "Pre-segment eval baseline (F2.1). F2.2's gate: K=1 segmented within "
        "1% of these numbers on the same machine and build profile "
        "(multi_environment_plan.md 4.6 (b)). Cross-machine comparison is "
        "meaningless -- re-record rather than compare. NOTE: 'eval' is "
        "measured whole and is not the sum of the other phases (see the "
        "module docstring); compare each phase against the same phase."
    )

    p = res["problem"]
    s = res["seconds"]
    print(f"design: {p['n_films']} films, {p['n_wavelengths']} wl x "
          f"{p['n_angles']} angles, {p['n_demands']} demands "
          f"({p['n_residuals']} residuals)")
    print(f"merit:  {res['merit']:.12e}")
    for name in ("assemble", "simulate", "merit_only", "simulate_merit", "eval"):
        print(f"  {name:11s} {s[name] * 1e6:10.2f} us")

    if args.write:
        RESULTS.parent.mkdir(parents=True, exist_ok=True)
        RESULTS.write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {RESULTS}")
    else:
        print("(not written -- pass --write to record the baseline)")


if __name__ == "__main__":
    main()
