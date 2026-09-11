# -*- coding: utf-8 -*-
"""The analytic Jacobian of the thickness optimizer (R4.5).

``J[i,k] = sum_terms dr_i/d(curve value) * d(curve value)/dd_k``: the merit
half from ``MeritSpec::curve_sensitivity``, the deposit half from the same
solver sweep that produces the curves. The per-element cross-check against
central differences is a cargo test (it needs the two halves separately);
this file is the end-to-end part CI runs, through the Python surface a user
actually touches.

Two facts about the comparison, both load-bearing:

* **The optimum has to be interior and unique.** Reflectance against
  thickness is oscillatory, so two solvers that differ only in their Jacobian
  can legitimately end in different local minima when started far from one.
  The targets here are read off a reference stack, which puts a merit-zero
  optimum a few nm from the start, in one basin, with no bound active.
* **The engine removes sub-minimum films** after optimizing (Navette's
  contract, not the LM's). ``clamp_min`` is floored to nothing so neither run
  can delete a layer and change the parameter count out from under the
  comparison.
"""

import numpy as np
import pytest

from navette._smatrix import LmConfig, SmatrixContext
from navette.materials import MaterialSpec
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec
from navette.synthesis.pipeline import stack_from_layers

_WL = np.linspace(480.0, 640.0, 9)
_ANGLES = np.array([0.0])
_CLAMP_MAX = 1000.0
_NO_REMOVAL_MIN = 1e-9

_INDICES = [2.30, 1.46, 2.30]
_REFERENCE = [118.0, 203.0, 64.0]
_START = [122.0, 198.0, 67.0]
_NAMES = ["H", "L", "H2"]


def _layers(thicknesses):
    return [(MaterialSpec("Konstant", dict(n=n)), d)
            for n, d in zip(_INDICES, thicknesses)]


def _ctx(spec, jacobian):
    return SmatrixContext(spec, _ANGLES, _WL, _NO_REMOVAL_MIN, _CLAMP_MAX,
                          LmConfig(jacobian=jacobian))


@pytest.fixture(scope="module")
def interior():
    """An R target the reference stack hits exactly, so the optimum is known."""
    # `SimCurves` has no Python reader for a raw curve, but a residual does
    # the job exactly: an Exact demand for zero with tolerance one returns
    # R itself, point for point.
    flat = TargetCollection()
    flat.add(SpectralTarget(_WL, np.zeros(_WL.size), np.ones(_WL.size),
                            0.0, "s", "R", kind="e"))
    flat_spec = build_merit_spec(flat)
    probe = _ctx(flat_spec, "analytic")
    ref, _ = stack_from_layers(_layers(_REFERENCE), _WL, {}, names=_NAMES)
    rs = np.asarray(flat_spec.residuals(probe.simulate(ref)), dtype=float)

    tc = TargetCollection()
    tc.add(SpectralTarget(_WL, rs, np.full(_WL.size, 0.02), 0.0, "s", "R",
                          kind="e"))
    return build_merit_spec(tc)


def _run(spec, jacobian, start=_START):
    ctx = _ctx(spec, jacobian)
    stack, _ = stack_from_layers(_layers(start), _WL, {}, names=_NAMES)
    mf, report = ctx.optimize_thicknesses_report(stack)
    d = [f["thickness"] for f in stack.films()]
    return mf, report, d


def test_the_analytic_path_is_the_default():
    assert "jacobian=Analytic" in repr(LmConfig())


def test_the_report_says_which_jacobian_actually_ran(interior):
    # Asserted, not inferred from a timing: a silent fallback would otherwise
    # look exactly like a fast analytic run.
    _, an, _ = _run(interior, "analytic")
    _, fd, _ = _run(interior, "fd")
    assert an["analytic_jacobians"] > 0
    assert fd["analytic_jacobians"] == 0


def test_the_two_modes_reach_the_same_optimum(interior):
    mf_an, _, d_an = _run(interior, "analytic")
    mf_fd, _, d_fd = _run(interior, "fd")
    assert len(d_an) == len(d_fd) == len(_REFERENCE)
    assert mf_an == pytest.approx(mf_fd, rel=1e-6, abs=1e-12)
    for a, b in zip(d_an, d_fd):
        assert a == pytest.approx(b, abs=1e-4)


def test_the_optimum_is_the_reference_stack(interior):
    # The targets were read off it, so merit zero is reachable and both runs
    # must find it. Without this the agreement test could pass on two runs
    # that agreed on the wrong answer.
    mf, _, d = _run(interior, "analytic")
    assert mf < 1e-12
    for got, want in zip(d, _REFERENCE):
        assert got == pytest.approx(want, abs=1e-3)


def test_the_analytic_path_spends_far_fewer_residual_evaluations(interior):
    # 2n central differences per iteration is what it removes; with three
    # films that is six evaluations an iteration, and the ratio grows with
    # the stack.
    _, an, _ = _run(interior, "analytic")
    _, fd, _ = _run(interior, "fd")
    assert an["evals"] < fd["evals"]


def test_a_phase_demand_falls_back_to_differences():
    # Phase rows have no curve sensitivity, so the source declines and the
    # driver differences — the run still has to work, and has to say so.
    tc = TargetCollection()
    tc.add(SpectralTarget(_WL, np.zeros(_WL.size), np.full(_WL.size, 0.05),
                          0.0, "s", "PDts", kind="e", phase=True))
    spec = build_merit_spec(tc)
    mf, report, d = _run(spec, "analytic")
    assert report["analytic_jacobians"] == 0
    assert np.isfinite(mf)
    assert len(d) == len(_START)


def test_an_unknown_jacobian_mode_is_refused():
    with pytest.raises(ValueError, match="analytic"):
        LmConfig(jacobian="exact")
