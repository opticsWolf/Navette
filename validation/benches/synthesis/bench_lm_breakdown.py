#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Benchmark: where a thickness optimization spends its time.

`bench_needle_insertion.py` shows the LM re-optimization is most of a
needle run (the scan is ~1%). This bench opens the LM up. Its wall time is

    iterations x (Jacobian + QR + damped solves)  +  evals x eval

so it can be cut two ways: cheaper iterations, or fewer of them. Each case
is measured both ways on the same seed.

Cost per iteration (the main table):

``t_eval``
    one `ctx.simulate` + `spec.residuals`, best-of mean. Timed from
    Python, so it carries a few microseconds of call overhead the LM's
    internal evals do not: `eval%` is an upper bound.
``lm``
    `optimize_thicknesses_report` with the analytic Jacobian at the
    default budget: wall, iterations, evals, Jacobian builds. `rest%` is
    the wall not spent in evals, i.e. the Jacobian builds plus the linear
    algebra; `rest/it` is that per iteration in units of `t_eval`.
``fd``
    the same with `jacobian="fd"` (central differences, rayon across
    columns): its wall over the analytic one is what the analytic chain
    saves, or costs.

Iterations needed (the budget table): the merit after 20, 50, 200 (the
default) and 1000 iterations, as a fraction of the merit at 200, with the
termination at 1000. A case still falling at 1000 is limited by how far
each step goes, not by what a step costs.

Cases vary one thing at a time from a 16-film AR coating on 61
wavelengths at normal incidence: film count, wavelength count, angle
count, a layer feature (roughness, inhomogeneity, interface slices), and
a differential-phase target, which the analytic chain cannot express (the
whole run is then differenced, whatever `jacobian` says).

Run explicitly:  python validation/benches/synthesis/bench_lm_breakdown.py
Fewer cases:     ... bench_lm_breakdown.py --quick
Write results:   ... bench_lm_breakdown.py --write
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

from navette._smatrix import LmConfig, SmatrixContext
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec
from navette.synthesis.pipeline import stack_from_layers

RESULTS = Path(__file__).resolve().parent / "results" / "lm_breakdown.json"

N_H, N_L, N_SUB = 2.35, 1.46, 1.52
TOL = 0.002
BUDGETS = (20, 50, 200, 1000)


def case(name, n_films=16, n_wl=61, angles=(0.0,), flags=None, phase=False):
    return dict(name=name, n_films=n_films, n_wl=n_wl, angles=list(angles),
                flags=flags or {}, phase=phase)


CASES = [
    case("base_16f"),
    case("films_4", n_films=4),
    case("films_8", n_films=8),
    case("films_32", n_films=32),
    case("wl_201", n_wl=201),
    case("wl_601", n_wl=601),
    case("angles_3", angles=(0.0, 30.0, 60.0)),
    case("rough_gauss", flags=dict(roughness=1.0, rough_type=4)),
    case("rough_nevot", flags=dict(roughness=1.0, rough_type=5)),
    case("inhomogen", flags=dict(inhomogen=True, inh_delta=0.05)),
    case("interface", flags=dict(interface=True, interface_thickness=2.0)),
    case("phase_PDts", phase=True),
]
QUICK = ("base_16f", "films_32", "wl_201", "rough_gauss", "phase_PDts")


def problem(c):
    """(layers, names, targets, wl): a perturbed quarter-wave AR seed.

    The thicknesses are off their quarter-wave values by a fixed pattern,
    so every case starts away from an optimum and the LM has work to do.
    """
    wl = np.linspace(400.0, 700.0, c["n_wl"])
    hi, lo = np.full(wl.shape, N_H + 0j), np.full(wl.shape, N_L + 0j)
    layers, names = [], []
    for i in range(c["n_films"]):
        h = i % 2 == 0
        qw = 550.0 / (4.0 * (N_H if h else N_L))
        layers.append((hi if h else lo, qw * (1.0 + 0.3 * np.sin(1.7 * i + 0.4))))
        names.append(f"{'H' if h else 'L'}{i}")
    tc = TargetCollection()
    for a in c["angles"]:
        for pol in ("s", "p"):
            tc.add(SpectralTarget(wl, np.zeros_like(wl), np.full_like(wl, TOL),
                                  pol, "R", angle_deg=a))
        if c["phase"]:
            tc.add(SpectralTarget(wl, np.zeros_like(wl), np.full_like(wl, 0.05),
                                  "s", "PDts", phase=True, angle_deg=a))
    return layers, names, tc, wl


def build(c):
    """(stack, targets, wl): a fresh native stack from the seed."""
    layers, names, tc, wl = problem(c)
    stack, _ = stack_from_layers(layers, wl, None, names=names,
                                 substrate=(N_SUB, "sub"),
                                 film_flags=c["flags"] or None)
    return stack, tc, wl


def context(c, spec, wl, **lm):
    return SmatrixContext(spec, wl, angles_deg=np.asarray(c["angles"], dtype=float),
                          lm=LmConfig(**lm))


def best_mean(fn, repeats, inner):
    best = float("inf")
    for _ in range(repeats):
        t0 = time.perf_counter()
        for _ in range(inner):
            fn()
        best = min(best, (time.perf_counter() - t0) / inner)
    return best


def lm_run(c, spec, wl, repeats, **lm):
    """Best-of wall seconds and the report of one LM from the seed."""
    best, rep = float("inf"), None
    for _ in range(repeats):
        stack = build(c)[0]
        ctx = context(c, spec, wl, **lm)
        t0 = time.perf_counter()
        mf, r = ctx.optimize_thicknesses_report(stack)
        dt = time.perf_counter() - t0
        if dt < best:
            best, rep = dt, dict(r, merit=mf)
    return best, rep


def measure(c, repeats, inner):
    stack, tc, wl = build(c)
    spec = build_merit_spec(tc)
    ctx = context(c, spec, wl)
    n_res = len(spec.residuals(ctx.simulate(stack)))
    t_eval = best_mean(lambda: spec.residuals(ctx.simulate(stack)), repeats, inner)

    wall, rep = lm_run(c, spec, wl, repeats)
    wall_fd, rep_fd = lm_run(c, spec, wl, repeats, jacobian="fd")
    budget = {}
    for it in BUDGETS:
        w, r = lm_run(c, spec, wl, 1, max_iterations=it)
        budget[it] = dict(wall=w, merit=r["merit"], iterations=r["iterations"],
                          termination=r["termination"])

    keys = ("iterations", "evals", "analytic_jacobians", "termination", "merit")
    rest = max(0.0, wall - rep["evals"] * t_eval)
    return {
        "case": c,
        "n_residuals": n_res,
        "t_eval": t_eval,
        "lm": dict(wall=wall, **{k: rep[k] for k in keys}),
        "fd": dict(wall=wall_fd, **{k: rep_fd[k] for k in keys}),
        "eval_share": rep["evals"] * t_eval / wall,
        "rest_share": rest / wall,
        "rest_per_iter_in_evals": rest / max(rep["iterations"], 1) / t_eval,
        "budget": budget,
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
    ap.add_argument("--quick", action="store_true", help=f"only {', '.join(QUICK)}")
    ap.add_argument("--write", action="store_true",
                    help="write results/lm_breakdown.json")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--inner", type=int, default=20)
    args = ap.parse_args()

    cases = [c for c in CASES if not args.quick or c["name"] in QUICK]
    rows = [measure(c, args.repeats, args.inner) for c in cases]

    print(f"--- cost per iteration: best of {args.repeats}, "
          f"t_eval over {args.inner} calls ---")
    print(f"  {'case':12s} {'res':>5} {'t_eval':>9} {'lm ms':>8} {'it':>4} "
          f"{'evals':>5} {'jac':>4} {'eval%':>6} {'rest%':>6} {'rest/it':>7} "
          f"{'fd ms':>8} {'fd/an':>6}")
    for r in rows:
        c, lm, fd = r["case"], r["lm"], r["fd"]
        print(f"  {c['name']:12s} {r['n_residuals']:5d} {r['t_eval'] * 1e6:7.1f}us "
              f"{lm['wall'] * 1e3:8.2f} {lm['iterations']:4d} {lm['evals']:5d} "
              f"{lm['analytic_jacobians']:4d} {r['eval_share'] * 100:5.1f}% "
              f"{r['rest_share'] * 100:5.1f}% {r['rest_per_iter_in_evals']:7.2f} "
              f"{fd['wall'] * 1e3:8.2f} {fd['wall'] / lm['wall']:6.2f}")
    print("  rest = Jacobian builds + linear algebra; rest/it in units of t_eval.")

    print("--- iterations needed: merit / merit at 200 iterations ---")
    print(f"  {'case':12s} " + " ".join(f"{f'@{b}':>8}" for b in BUDGETS)
          + f"  {'ms @1000':>9}  stop @1000")
    for r in rows:
        ref = r["budget"][200]["merit"]
        b1000 = r["budget"][1000]
        print(f"  {r['case']['name']:12s} "
              + " ".join(f"{r['budget'][b]['merit'] / ref:8.4f}" for b in BUDGETS)
              + f"  {b1000['wall'] * 1e3:9.1f}  {b1000['termination']}"
              f" ({b1000['iterations']} it)")

    res = {"cases": rows, "provenance": bench_provenance()}
    res["provenance"]["commit"] = git_hash()
    res["provenance"]["recorded_utc"] = datetime.now(timezone.utc).isoformat(
        timespec="seconds")
    if args.write:
        RESULTS.parent.mkdir(parents=True, exist_ok=True)
        RESULTS.write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {RESULTS}")
    else:
        print("(not written -- pass --write to record the results)")


if __name__ == "__main__":
    main()
