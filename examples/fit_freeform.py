# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Fit the optical constants of an unknown film, and the instrument too.

A 420 nm film of a material no dispersion formula describes was measured
as R at 8 deg, T at normal incidence, and Rs, Rp at 60 deg, 400-1000 nm.
The spectrophotometer's wavelength scale was off by 1.2 nm on the T
channel. The fit recovers:

  1. the film's n and k as a ``Spline`` on seven knots (pchip), every knot
     value a fit parameter (``material[X].n[3]``), with a smoothness prior
     on n (second differences, one row per inner knot);
  2. the film thickness;
  3. the T channel's wavelength shift, an instrument parameter
     (``data[T].wavelength_shift``).

``BSplineEps`` (eps2 as a B-spline, eps1 by Kramers-Kronig) and
``PointWise`` (one value per wavelength) are the other free-form models.

Run:  python examples/fit_freeform.py [--plot]
"""

from __future__ import annotations

import sys

import numpy as np

from navette.fit import Dataset, Environment, FitProblem, Parameter, evaluate, fit, plot_fit
from navette.materials import MODELS, evaluate as nk_of

KNOTS = np.linspace(400.0, 1000.0, 7)
TRUE_N = 2.05 + 0.025 / (KNOTS / 1000.0) ** 2
TRUE_K = 0.002 + 0.03 * np.exp(-(KNOTS - 400.0) / 80.0)
TRUE_D = 420.0
TRUE_SHIFT = 1.2
WL = np.linspace(400.0, 1000.0, 121)
SIGMA = 0.001


def model(n, k, d: float) -> dict:
    film = {"wavelengths": list(KNOTS), "n": list(n), "k": list(k), "method": "pchip"}
    return {
        "materials": [
            {"name": "air", "model": "Konstant", "params": {"n": 1.0}},
            {"name": "glass", "model": "Cauchy", "params": {"A": 1.5, "B": 0.004, "C": 0.0}},
            {"name": "X", "model": "Spline", "params": film},
        ],
        "structures": [{
            "label": "film",
            "layers": [{"material_code": "X", "thickness_nm": d}],
            "groups": [],
        }],
    }


SAMPLE = [Environment("sample", ambient="air", substrate="glass", stack=[{"structure": "film"}])]


CURVES = [("R", "R_avg", 8.0), ("T", "T_avg", 0.0), ("Rs", "Rs", 60.0), ("Rp", "Rp", 60.0)]


def measured() -> list:
    """The truth simulated, with the T channel's shift and seeded noise."""
    slots = [Dataset(i, ch, WL, np.zeros_like(WL), sigma=SIGMA, angle_deg=a,
                     instrument={"wavelength_shift": TRUE_SHIFT} if i == "T" else None)
             for i, ch, a in CURVES]
    truth = evaluate(FitProblem(model(TRUE_N, TRUE_K, TRUE_D), data=slots, environments=SAMPLE))
    rng = np.random.default_rng(7)
    for slot, curve in zip(slots, truth.curves):
        slot.values = curve.simulated + rng.normal(0.0, SIGMA, WL.size)
        slot.instrument = None  # the fit does not know the shift
    return slots


def main() -> int:
    assert {"PointWise", "Spline", "BSplineEps"} <= set(MODELS)
    at_knots = nk_of({"model": "Spline", "params": {
        "wavelengths": list(KNOTS), "n": list(TRUE_N), "k": list(TRUE_K)}}, KNOTS)
    assert np.allclose(at_knots, TRUE_N + 1j * TRUE_K), "a spline passes its knots"

    knots = [f"material[X].{q}[{i}]" for q in ("n", "k") for i in range(len(KNOTS))]
    problem = FitProblem(
        model(np.full(7, 2.1), np.full(7, 0.005), 410.0),
        data=measured(),
        environments=SAMPLE,
        parameters=[Parameter(a) for a in knots]
        + [Parameter("layer[film:0].thickness"), Parameter("data[T].wavelength_shift")],
        smoothness=[("X", "n", 0.01)],
    )
    report = fit(problem)
    print(f"chi2_red {report.chi2_reduced:.3f}  ({report.termination}, "
          f"{report.iterations} iterations) {report.covariance_note or ''}")
    truths = dict(zip(knots, list(TRUE_N) + list(TRUE_K)))
    truths["layer[film:0].thickness"] = TRUE_D
    truths["data[T].wavelength_shift"] = TRUE_SHIFT
    worst = 0.0
    for addr, want in truths.items():
        p = report.parameter(addr)
        z = abs(p.value - want) / p.std_error
        worst = max(worst, z)
        print(f"  {addr:<28} {p.value:9.4f} +- {p.std_error:.4f}   truth {want:9.4f}")
    print(f"worst deviation {worst:.2f} standard errors")
    if "--plot" in sys.argv:
        import matplotlib.pyplot as plt
        plot_fit(report)
        plt.show()
    ok = worst < 3.0 and abs(report.chi2_reduced - 1.0) < 0.3
    print("recovered" if ok else "NOT recovered")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
