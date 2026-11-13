# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Fit a deposited coating back from its measured spectra.

A six-layer TiO2/SiO2 coating on glass was designed, deposited, and
measured: R at 8 deg, T at normal incidence and Rp at 45 deg, 400-1000 nm.
The coater's groups carry their error laws (2 % thickness, 1 % index,
relative), so a "deposited" truth is drawn from them, and the spectra are
simulated on it with 0.2 % noise. The fit recovers that truth from the
design, which is what a reverse-engineering run does with real data.

Steps:

  1. The problem: the program's model, the sample it was measured on, the
     data, and the scalars that may move (every film thickness, and each
     group's index factor).
  2. A two-stage strategy: thicknesses first, then everything, with a
     progress callback that prints each stage.
  3. A tournament of strategies on the same data, ranked by chi2, and a
     small benchmark that ranks them by how well they recover truths.
  4. The problem saved as a program file and read back; a callback that
     raises stops a fit.

Run:  python examples/fit_coating.py [--plot]
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

from navette.fit import (
    Dataset,
    Environment,
    FitProblem,
    Parameter,
    Stage,
    Strategy,
    benchmark,
    draw_truth,
    fit,
    orderings,
    plot_fit,
    tournament,
)

# -- 1. The problem ----------------------------------------------------------


def erring(name: str) -> dict:
    """A group whose thickness and index scatter by 2 % and 1 % per run."""
    return {
        "name": name,
        "error_mask": [1, 1, 0, 0, 0, 0],
        "thickness_error_params": {"abs_std_dev": 0.0, "abs_variance": 0.0,
                                   "rel_std_dev": 0.02},
        "n_error_params": {"abs_std_dev": 0.0, "abs_variance": 0.0, "rel_std_dev": 0.01},
    }


DESIGN_NM = [("H", 21.0), ("L", 35.0), ("H", 110.0), ("L", 22.0), ("H", 64.0), ("L", 95.0)]

MODEL = {
    "materials": [
        {"name": "air", "model": "Konstant", "params": {"n": 1.0}},
        {"name": "glass", "model": "Konstant", "params": {"n": 1.52}},
        {"name": "H", "model": "Cauchy", "params": {"A": 2.2, "B": 0.02, "C": 0.001}},
        {"name": "L", "model": "Cauchy", "params": {"A": 1.45, "B": 0.0036, "C": 0.0}},
    ],
    "groups": [erring("H"), erring("L")],
    "structures": [{
        "label": "coat",
        "layers": [{"material_code": m, "thickness_nm": d} for m, d in DESIGN_NM],
        "groups": [],
    }],
}

WL = np.linspace(400.0, 1000.0, 121)


def placeholder(id_: str, channel: str, angle_deg: float) -> Dataset:
    """A measurement slot; draw_truth fills in the values."""
    return Dataset(id_, channel, WL, np.zeros_like(WL), sigma=0.002, angle_deg=angle_deg)


def problem() -> FitProblem:
    return FitProblem(
        model=MODEL,
        environments=[Environment("sample", ambient="air", substrate="glass",
                                  stack=[{"structure": "coat"}])],
        data=[placeholder("R@8", "R_avg", 8.0), placeholder("T@0", "T_avg", 0.0),
              placeholder("Rp@45", "Rp", 45.0)],
        parameters=[
            Parameter(select="kind=thickness"),
            Parameter("group[H].n_factor"),
            Parameter("group[L].n_factor"),
        ],
    )


def main() -> int:
    nominal = problem()
    truth = draw_truth(nominal, seed=3)
    measured = problem()
    measured.data = truth.data
    print("deposited truth (drawn from the error laws):")
    for addr, v in truth.values.items():
        print(f"  {addr:<26} {v:9.4f}")

    # -- 2. A staged fit with progress -----------------------------------
    strategy = Strategy("thicknesses, then all", [
        Stage("thicknesses", ["set:thicknesses"]),
        Stage("all", ["*"]),
    ])
    improved = [0]

    def progress(event: dict) -> None:
        if event["event"] == "stage_start":
            print(f"stage {event['name']!r}: {len(event['free'])} free")
        elif event["event"] == "improved":
            improved[0] += 1
        else:
            print(f"  ... chi2 {event['chi2']:.4g} after {improved[0]} improvements")

    report = fit(measured, strategy, progress=progress)
    print()
    print(report.summary())
    worst = max(abs(report.values[a] / v - 1.0) for a, v in truth.values.items())
    print(f"\nworst relative error against the truth: {worst:.2e}")
    assert report.chi2_reduced < 2.0, report.chi2_reduced
    assert worst < 0.02, worst
    assert improved[0] > 0

    # -- 3. Which strategy works best? -----------------------------------
    strategies = [
        Strategy("all at once", [Stage("all", ["*"])]),
        Strategy("groups only", [Stage("groups", ["set:groups"])]),
        orderings(["thicknesses", "groups"]),
    ]
    ranked, _ = tournament(measured, strategies)
    print("\ntournament on the measured data:")
    for e in ranked:
        print(f"  {e.rank}. {e.name:<32} chi2_red {e.chi2_reduced:9.3f}  {e.evals:5d} evals")
    assert ranked[-1].name == "groups only"

    bench = benchmark(nominal, strategies, seeds=[1, 2])
    print("\nbenchmark over two drawn truths:")
    print(bench.summary())
    assert bench.entries[-1].name == "groups only"

    # -- 4. Saved, read back, stopped ------------------------------------
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "fit_coating.json"
        measured.strategy = strategy
        measured.save(path)
        again = FitProblem.load(path)
    assert again.to_document() == measured.to_document()
    assert fit(again).chi2 == report.chi2
    print("\nsaved as a program file and read back: the same fit")

    def impatient(event: dict) -> None:
        if event["event"] == "improved":
            raise KeyboardInterrupt("enough")

    try:
        fit(measured, progress=impatient)
    except ValueError as e:
        assert "enough" in str(e)
        print(f"a callback that raises stops the fit: {e}")
    else:
        raise AssertionError("the fit did not stop")

    if "--plot" in sys.argv:
        import matplotlib.pyplot as plt
        plot_fit(report)
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
