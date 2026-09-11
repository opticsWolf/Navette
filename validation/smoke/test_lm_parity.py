# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""The bounded LM against ``scipy.optimize.least_squares`` (R4.4).

``docs/`` has claimed scipy parity for ``synthesis/thick_opt.rs`` since it
replaced ``least_squares(method="trf")``; review §18.2 found no evidence it had
ever been reproduced. ``validation/review/lm_check.py`` is the full harness
(three film counts, far-start behaviour, the cargo tests' pinned optima
recomputed). This file is the part CI runs.

Two facts about the comparison, both learned by running it:

* **Reflectance against thickness is oscillatory**, so two local solvers
  started far from an optimum legitimately land in different basins. Cost
  parity is only meaningful when both start inside one, which is what the
  fixture arranges -- it anchors on scipy's own answer and perturbs.
* **The engine removes sub-minimum films** after optimizing (Navette's
  contract, not the LM's), so ``clamp_min`` is set to nothing here; a film
  scipy leaves at 0.0 is one the engine would delete, changing the parameter
  count and leaving nothing to compare.
"""

import numpy as np
import pytest
from scipy.optimize import least_squares

from navette._smatrix import SmatrixContext
from navette.materials import MaterialSpec
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec
from navette.synthesis.pipeline import stack_from_layers

_WL = np.linspace(450.0, 750.0, 13)
_ANGLES = np.array([0.0])
_CLAMP_MAX = 1000.0
_NO_REMOVAL_MIN = 1e-9

_FILMS = [(2.35, 80.0), (1.46, 140.0), (2.10, 70.0)]


@pytest.fixture(scope="module")
def bench():
    tc = TargetCollection()
    tc.add(SpectralTarget(_WL, np.zeros(_WL.size), np.full(_WL.size, 0.01),
                          0.0, "s", "R", kind="e", weight=2.0))
    spec = build_merit_spec(tc)
    layers = [(MaterialSpec("Konstant", dict(n=n)), d) for n, d in _FILMS]
    names = [f"L{i}" for i in range(len(_FILMS))]

    def ctx():
        return SmatrixContext(spec, _ANGLES, _WL, _NO_REMOVAL_MIN, _CLAMP_MAX)

    def scipy_from(x0):
        c = ctx()
        stack, _ = stack_from_layers(layers, _WL, {}, names=names)

        def fun(x):
            for i, v in enumerate(x):
                stack.set_thickness(i, float(v))
            return np.asarray(spec.residuals(c.simulate(stack)), dtype=float)

        sp = least_squares(fun, np.asarray(x0, float), method="trf",
                           bounds=(np.zeros(len(x0)), np.full(len(x0), _CLAMP_MAX)),
                           xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=20000)
        return sp.x.copy(), float(np.sum(sp.fun ** 2))

    def engine_from(x0):
        c = ctx()
        stack, _ = stack_from_layers(
            [(m, float(d)) for (m, _), d in zip(layers, x0)], _WL, {}, names=names)
        cost = c.optimize_thicknesses(stack)
        return np.array([f["thickness"] for f in stack.films()]), float(cost)

    anchor, _ = scipy_from([d for _, d in _FILMS])
    start = np.clip(anchor + np.array([3.0, -3.0, 2.0]), 5.0, _CLAMP_MAX - 5.0)
    return {"spec": spec, "layers": layers, "names": names, "ctx": ctx,
            "scipy_from": scipy_from, "engine_from": engine_from,
            "start": start}


def test_the_same_basin_gives_the_same_cost(bench):
    _, e_cost = bench["engine_from"](bench["start"])
    _, s_cost = bench["scipy_from"](bench["start"])
    assert e_cost == pytest.approx(s_cost, rel=1e-6)


def test_the_same_basin_gives_the_same_thicknesses(bench):
    e_x, _ = bench["engine_from"](bench["start"])
    s_x, _ = bench["scipy_from"](bench["start"])
    assert e_x.size == s_x.size
    assert np.max(np.abs(e_x - s_x)) < 1e-3 * max(1.0, float(e_x.max()))


def test_the_engine_answer_is_a_fixed_point(bench):
    """A solver that keeps stepping at a stationary point is not converged."""
    x1, c1 = bench["engine_from"]([d for _, d in _FILMS])
    _, c2 = bench["engine_from"](x1)
    assert c2 >= c1 - 1e-6 * abs(c1)


def test_a_far_start_still_improves_the_merit(bench):
    """Basin-free: whichever minimum it reaches, it may not end up worse."""
    stack, _ = stack_from_layers(bench["layers"], _WL, {}, names=bench["names"])
    before = bench["ctx"]().evaluate_merit(stack)
    _, after = bench["engine_from"]([d for _, d in _FILMS])
    assert after < before


@pytest.mark.parametrize("name,fun,x0,bounds,want,tol", [
    # linear_least_squares_exact_recovery
    ("linear",
     lambda p: p[0] * np.arange(5.0) + p[1] - np.array([-1.0, 1.0, 3.0, 5.0, 7.0]),
     [0.0, 0.0], (-1e3, 1e3), [2.0, -1.0], 1e-8),
    # exponential_fit_nonlinear
    ("exponential",
     lambda p: (p[0] * np.exp(-p[1] * np.arange(0.0, 3.5, 0.5))
                - 2.5 * np.exp(-0.7 * np.arange(0.0, 3.5, 0.5))),
     [1.0, 0.2], ([0.1, 0.1], [10.0, 5.0]), [2.5, 0.7], 1e-6),
    # bound_active_run_still_converges_to_the_boundary_optimum
    ("boundary",
     lambda p: np.array([p[0] - 5.0, p[1] + 5.0, 0.25 * (p[0] - p[1])]),
     [0.5, 0.5], ([0.0, 0.0], [1.0, 1.0]), [1.0, 0.0], 1e-9),
    # a_clamped_step_reports_the_gain_ratio_of_the_step_it_took
    ("clamped_corner",
     lambda p: p[0] * np.arange(5.0) + p[1] - (100.0 * np.arange(5.0) - 50.0),
     [0.5, 0.5], ([0.0, 0.0], [1.0, 1.0]), [1.0, 1.0], 1e-9),
])
def test_the_cargo_tests_pinned_optima_are_what_scipy_finds(
        name, fun, x0, bounds, want, tol):
    """The constants in ``thick_opt.rs``'s unit tests were written by hand.

    Without this, a wrong pin would be invisible: the solver agrees with it by
    construction. There is no Python entry point for the engine LM on an
    arbitrary residual closure (the FD Jacobian is rayon-parallel, so a Python
    callback would need the GIL inside every worker), so this checks scipy
    against the pinned values; the tests above are what exercise the Rust
    solver.
    """
    sp = least_squares(fun, np.array(x0, float), method="trf", bounds=bounds,
                       xtol=1e-15, ftol=1e-15, gtol=1e-15)
    assert np.max(np.abs(sp.x - np.array(want))) < tol, f"{name}: {sp.x}"
