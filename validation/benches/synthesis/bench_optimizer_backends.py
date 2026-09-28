#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compare optional solvers on the review harness's thin-film cases (issue #2)."""

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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
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
        "backends": available_optimizers(), "warmups": 1, "repetitions": 5,
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
                    run_engine(**case_arguments(case), optimizer=backend, jacobian=jacobian)
                    samples = []
                    for _ in range(5):
                        start = time.perf_counter()
                        x, cost, report = run_engine(**case_arguments(case), optimizer=backend,
                                                    jacobian=jacobian)
                        samples.append(time.perf_counter() - start)
                    row.update(seconds=samples, median_seconds=statistics.median(samples),
                               cost=cost, thicknesses_nm=x.tolist(),
                               iterations=report["iterations"], evals=report["evals"],
                               analytic_jacobians=report["analytic_jacobians"],
                               termination=report["termination"])
                    print(f"{case['label']:22} {backend:20} {jacobian:8} "
                          f"{1e3 * row['median_seconds']:9.3f} ms  "
                          f"cost={cost:.8g}  evals={report['evals']}  {report['termination']}")
                except ValueError as error:
                    row["error"] = str(error)
                    print(f"{case['label']:22} {backend:20} {jacobian:8} ERROR {error}")
                rows.append(row)
    result = {"provenance": provenance, "results": rows}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
