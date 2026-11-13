# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Optical constants from measurements: a dispersion fit and a film fit.

Two common jobs:

  1. **A dispersion model for a measured n, k table.** A film's optical
     constants were measured on a grid (here a Sellmeier glass plus noise).
     :func:`table_to_model` ranks candidate models by BIC, and
     :func:`fit_dispersion` fits the chosen one to the table with its
     uncertainties.
  2. **A film's thickness and index from a reflectance file.** The
     spectrophotometer writes R in percent at 8 deg; :func:`to_solver`
     reads the file in that convention (a preset), and the fit recovers the
     film's thickness and its Cauchy A and B together.

Run:  python examples/fit_optical_constants.py
"""

from __future__ import annotations

import sys

import numpy as np

from navette.fit import (
    Dataset,
    Environment,
    FitProblem,
    NkDataset,
    Parameter,
    conventions,
    evaluate,
    fit,
    fit_dispersion,
    from_solver,
    table_to_model,
    to_solver,
)

RNG = np.random.default_rng(7)

# -- 1. A model for a measured table -------------------------------------------

SELLMEIER = {"model": "Sellmeier", "params": {
    "B1": 0.6961663, "C1": 0.0046791, "B2": 0.4079426, "C2": 0.0135121,
    "B3": 0.8974794, "C3": 97.934}}


def true_n(wl: np.ndarray) -> np.ndarray:
    """n of fused silica (Malitson's Sellmeier)."""
    um2 = (wl / 1000.0) ** 2
    p = SELLMEIER["params"]
    return np.sqrt(1.0 + sum(p[f"B{i}"] * um2 / (um2 - p[f"C{i}"]) for i in (1, 2, 3)))


def dispersion() -> None:
    # The measured table: 250-2000 nm with 2e-4 noise.
    wl = np.linspace(250.0, 2000.0, 90)
    n = true_n(wl) + RNG.normal(0.0, 2e-4, wl.size)
    table = {"model": "Table", "params": {"n_data": [wl.tolist(), n.tolist()]}}
    ranked = table_to_model(table, [
        {"model": "Konstant", "params": {"n": 1.46}},
        {"model": "Cauchy", "params": {"A": 1.45, "B": 0.004, "C": 0.0}},
        {"model": "Sellmeier", "params": {"B1": 0.7, "C1": 0.005, "B2": 0.4,
                                          "C2": 0.0135121, "B3": 0.9, "C3": 100.0}},
    ])
    print("models for the table, best first (BIC):")
    for c in ranked:
        print(f"  {c['model']:<10} BIC {c['bic']:10.1f}")
    assert [c["model"] for c in ranked] == ["Sellmeier", "Cauchy", "Konstant"]

    # The chosen model against the table, with the noise as sigma.
    data = [NkDataset("table", "silica", wl, {"n": (n, 2e-4)})]
    start = {"model": "Sellmeier", "params": dict(ranked[0]["spec"]["params"])}
    spec, report = fit_dispersion(start, data)
    print("\nSellmeier fitted to the table:")
    print(report.summary())
    assert abs(report.chi2_reduced - 1.0) < 0.5, report.chi2_reduced
    # Six terms over this range trade off against each other (see the
    # error bars and report.correlation), so the coefficients are not the
    # generating ones; the curve is, to within the noise.
    c = report.curves[0]
    rms = float(np.sqrt(np.mean((c.simulated - true_n(c.wavelengths)) ** 2)))
    print(f"  rms against the noiseless index: {rms:.1e} (noise 2e-4)")
    assert rms < 2e-4, rms
    assert spec["model"] == "Sellmeier"


# -- 2. A film from a reflectance file ---------------------------------------

TRUTH = {"layer[film:0].thickness": 412.0, "material[film].A": 1.95,
         "material[film].B": 0.012}


def film_problem(data) -> FitProblem:
    return FitProblem(
        model={
            "materials": [
                {"name": "air", "model": "Konstant", "params": {"n": 1.0}},
                {"name": "glass", "model": "Konstant", "params": {"n": 1.52}},
                {"name": "film", "model": "Cauchy",
                 "params": {"A": 1.9, "B": 0.01, "C": 0.0}},
            ],
            "structures": [{"label": "film", "groups": [],
                            "layers": [{"material_code": "film", "thickness_nm": 400.0}]}],
        },
        environments=[Environment("sample", ambient="air", substrate="glass",
                                  stack=[{"structure": "film"}])],
        data=data,
        parameters=[
            Parameter("layer[film:0].thickness", lb=350.0, ub=450.0),
            Parameter("material[film].A", lb=1.7, ub=2.2),
            Parameter("material[film].B", lb=0.0, ub=0.05),
        ],
    )


def instrument_file() -> dict:
    """What the spectrophotometer wrote: R in percent at 8 deg."""
    wl = np.linspace(400.0, 900.0, 201)
    slot = Dataset("R", "R_avg", wl, np.zeros_like(wl), sigma=0.002, angle_deg=8.0)
    truth = film_problem([slot])
    truth.parameters = [Parameter(a, value=v, fixed=True) for a, v in TRUTH.items()]
    r = evaluate(truth).curves[0].simulated + RNG.normal(0.0, 0.002, wl.size)
    return {"id": "R", "quantity": {"spectrum": "R_avg"}, "angle": 8.0,
            "abscissa": wl, "columns": [100.0 * r]}


def film() -> None:
    raw = instrument_file()
    convention = {**conventions()["nm_deg_percent"], "sigma": {"constant": 0.2}}
    data = to_solver(raw, convention)
    print(f"\n{data[0].id}: {data[0].source}")
    assert np.allclose(data[0].values, raw["columns"][0] / 100.0)
    assert np.allclose(data[0].sigma, 0.002)

    report = fit(film_problem(data))
    print(report.summary())
    for addr, v in TRUTH.items():
        p = report.parameter(addr)
        assert abs(p.value - v) < 3.0 * p.std_error, (addr, p)

    back = from_solver("R", data, convention)
    assert np.allclose(back["columns"][0], raw["columns"][0])
    print("written back in the instrument's convention: the same file")


def main() -> int:
    dispersion()
    film()
    return 0


if __name__ == "__main__":
    sys.exit(main())
