# -*- coding: utf-8 -*-
"""Which least-squares solver runs, and how a build says what it has (R4.4c).

The alternative backends are cargo features, so what a wheel can run is a
property of how it was built. That makes two things testable everywhere, and
one only where the feature is on:

* **everywhere** — the build can be *asked* what it has
  (``available_optimizers()``), every name it reports actually constructs, and
  a name it does not have is refused with the command that would supply it.
  Refusing is the contract: a missing backend must never fall back to a
  different solver, because then a comparison against "the reference LM" would
  quietly be a comparison against ourselves;
* **only with ``opt-minpack-lm``** — that the reference LM lands on the same
  optimum as the built-in.

Bounds are deliberately not compared. The built-in's veto+clamp lets a
thickness finish *on* a bound; the unbounded backends run on an interior
reparametrization and finish strictly inside one. That is a documented
contract difference, so the comparison below uses a target whose optimum is
interior for both.
"""

import numpy as np
import pytest

from navette._smatrix import LmConfig, SmatrixContext, available_optimizers
from navette.materials import MaterialSpec
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec
from navette.synthesis.pipeline import stack_from_layers

_WL = np.linspace(480.0, 640.0, 9)
_ANGLES = np.array([0.0])
_CLAMP_MAX = 1000.0
# The engine deletes sub-minimum films after optimizing, which would change
# the parameter count out from under a backend comparison.
_NO_REMOVAL_MIN = 1e-9

_INDICES = [2.30, 1.46, 2.30]
_REFERENCE = [118.0, 203.0, 64.0]
_START = [122.0, 198.0, 67.0]
_NAMES = ["H", "L", "H2"]

_HAVE_MINPACK = "minpack_lm" in available_optimizers()


def _layers(thicknesses):
    return [(MaterialSpec("Konstant", dict(n=n)), d)
            for n, d in zip(_INDICES, thicknesses)]


def _ctx(spec, **lm):
    return SmatrixContext(spec, _ANGLES, _WL, _NO_REMOVAL_MIN, _CLAMP_MAX,
                          LmConfig(**lm))


@pytest.fixture(scope="module")
def interior():
    """An R target the reference stack hits exactly: optimum interior, merit 0."""
    flat = TargetCollection()
    flat.add(SpectralTarget(_WL, np.zeros(_WL.size), np.ones(_WL.size),
                            0.0, "s", "R", kind="e"))
    flat_spec = build_merit_spec(flat)
    ref, _ = stack_from_layers(_layers(_REFERENCE), _WL, {}, names=_NAMES)
    rs = np.asarray(flat_spec.residuals(_ctx(flat_spec).simulate(ref)),
                    dtype=float)

    tc = TargetCollection()
    tc.add(SpectralTarget(_WL, rs, np.full(_WL.size, 0.02), 0.0, "s", "R",
                          kind="e"))
    return build_merit_spec(tc)


def _run(spec, optimizer):
    ctx = _ctx(spec, optimizer=optimizer)
    stack, _ = stack_from_layers(_layers(_START), _WL, {}, names=_NAMES)
    mf, report = ctx.optimize_thicknesses_report(stack)
    return mf, report, [f["thickness"] for f in stack.films()]


def test_the_build_can_be_asked_what_it_has():
    names = available_optimizers()
    assert "builtin" in names
    for n in names:
        assert repr(LmConfig(optimizer=n))


def test_the_builtin_backend_is_the_default():
    assert 'optimizer="builtin"' in repr(LmConfig())


def test_an_unknown_backend_is_refused_and_says_what_exists():
    with pytest.raises(ValueError, match="minpack_lm"):
        LmConfig(optimizer="scipy")


@pytest.mark.skipif(_HAVE_MINPACK,
                    reason="this build has opt-minpack-lm compiled in")
def test_a_backend_this_build_lacks_is_refused_with_the_rebuild_command():
    # Not a fallback. A wheel that silently ran the built-in here would make
    # every "reference LM" comparison a comparison with itself.
    with pytest.raises(ValueError, match="opt-minpack-lm"):
        LmConfig(optimizer="minpack_lm")


def test_the_report_says_which_backend_ran(interior):
    _, report, _ = _run(interior, "builtin")
    assert report["backend"] == "builtin"


@pytest.mark.skipif(not _HAVE_MINPACK,
                    reason="needs the opt-minpack-lm cargo feature")
def test_the_reference_lm_finds_the_same_interior_optimum(interior):
    mf_a, rep_a, d_a = _run(interior, "builtin")
    mf_b, rep_b, d_b = _run(interior, "minpack_lm")
    assert rep_b["backend"] == "minpack_lm"
    assert len(d_a) == len(d_b) == len(_REFERENCE)
    assert mf_b == pytest.approx(mf_a, rel=1e-5, abs=1e-12)
    for got, want in zip(d_b, _REFERENCE):
        assert got == pytest.approx(want, abs=1e-2)
