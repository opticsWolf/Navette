#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compare optional solvers on the review harness's thin-film cases (issue #2).

Each (case, Jacobian mode, backend) is warmed up, then timed ``--repetitions``
times on a fresh stack. The committed result keeps a five-number summary of
those wall times -- min, quartiles, median, max, in milliseconds -- rather
than every sample: the median is the headline, the interquartile range says
how far to trust it, and min/max catch a run that was disturbed. Pass
``--raw`` to keep the samples as well (for a local investigation, not for the
committed file).

Solver work (iterations, evaluations, termination) is deterministic and is
recorded once per row. Timing is informational and gates nothing.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _bench_common import bench_provenance, require_release, setup_bench

setup_bench()
require_release()
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "review"))
from lm_check import basin_cases, case_arguments, run_engine

import argparse
import json
import os
import platform
import statistics
import subprocess
import time
import tomllib

import numpy as np
import scipy

from navette._smatrix import available_optimizers


def summarize_ms(samples):
    """Five-number summary of wall times in seconds, reported in ms."""
    ms = sorted(1e3 * s for s in samples)
    q1, median, q3 = statistics.quantiles(ms, n=4, method="inclusive")
    return {"n": len(ms), "min": ms[0], "q1": q1, "median": median, "q3": q3, "max": ms[-1]}


def rounded(value, digits=4):
    """Round floats in nested data to `digits` significant figures for timing."""
    if isinstance(value, float):
        return float(f"{value:.{digits}g}")
    if isinstance(value, dict):
        return {k: rounded(v, digits) for k, v in value.items()}
    if isinstance(value, list):
        return [rounded(v, digits) for v in value]
    return value


def dump(result):
    """Provenance indented; one line per result row, so the file diffs by row."""
    rows = ",\n".join("    " + json.dumps(row, allow_nan=False) for row in result["results"])
    provenance = json.dumps(result["provenance"], indent=2, allow_nan=False).replace("\n", "\n  ")
    return f'{{\n  "provenance": {provenance},\n  "results": [\n{rows}\n  ]\n}}\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repetitions", type=int, default=21,
                        help="timed runs per row (default 21; odd keeps the median a sample)")
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--raw", action="store_true", help="also keep every timing sample")
    args = parser.parse_args()
    if args.repetitions < 2:
        parser.error("--repetitions must be at least 2 for quartiles")
    repo = Path(__file__).resolve().parents[3]
    lock = tomllib.loads((repo / "Cargo.lock").read_text())
    provenance = bench_provenance()
    cpu = platform.processor()
    cpuinfo = Path("/proc/cpuinfo")
    if not cpu and cpuinfo.exists():
        cpu = next((line.split(":", 1)[1].strip() for line in cpuinfo.read_text().splitlines()
                    if line.startswith("model name")), platform.machine())
    provenance.update({
        "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
        "dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=repo, text=True)),
        "numpy": np.__version__, "scipy": scipy.__version__,
        "cpu": cpu, "logical_cpus": os.cpu_count(),
        "threads": {name: os.environ.get(name, "runtime default") for name in
                    ("RAYON_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")},
        "solver_versions": {p["name"]: p["version"] for p in lock["package"]
                            if p["name"] in {"basin", "argmin", "levenberg-marquardt"}},
        "backends": available_optimizers(),
        "warmups": args.warmups, "repetitions": args.repetitions,
        "timing": "wall ms per fresh-stack run: n, min, q1, median, q3, max "
                  "(inclusive quartiles), 4 significant figures",
        "max_iterations": 200, "max_evals": 100_000,
        "ftol": 1e-12, "xtol": 1e-12, "gtol": 1e-10,
    })
    rows = []
    for case in basin_cases():
        for jacobian in ("analytic", "fd"):
            for backend in available_optimizers():
                row = {"case": case["label"], "backend": backend, "jacobian": jacobian,
                       "start_nm": case["x0"].tolist(), "bounds_nm": [0.0, case["clamp_max"]]}
                try:
                    for _ in range(args.warmups):
                        run_engine(**case_arguments(case), optimizer=backend, jacobian=jacobian)
                    samples = []
                    for _ in range(args.repetitions):
                        start = time.perf_counter()
                        x, cost, report = run_engine(**case_arguments(case), optimizer=backend,
                                                    jacobian=jacobian)
                        samples.append(time.perf_counter() - start)
                    ms = summarize_ms(samples)
                    row.update(ms=rounded(ms), cost=cost, thicknesses_nm=x.tolist(),
                               iterations=report["iterations"], evals=report["evals"],
                               analytic_jacobians=report["analytic_jacobians"],
                               termination=report["termination"])
                    if args.raw:
                        row["samples_ms"] = [1e3 * s for s in samples]
                    print(f"{case['label']:22} {backend:20} {jacobian:8} "
                          f"{ms['median']:9.3f} ms  [{ms['q1']:.3f}, {ms['q3']:.3f}]  "
                          f"cost={cost:.8g}  evals={report['evals']}  {report['termination']}")
                except ValueError as error:
                    row["error"] = str(error)
                    print(f"{case['label']:22} {backend:20} {jacobian:8} ERROR {error}")
                rows.append(row)
    result = {"provenance": provenance, "results": rows}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(dump(result))


if __name__ == "__main__":
    main()
