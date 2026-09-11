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

``trf`` (R4.6) is hand-rolled and needs no feature, so it is in the first
group: always present, always selectable. It is the one alternative backend
with bounds of its own, so unlike the reparametrized ones it is compared on a
*bound-active* problem as well as an interior one.

Bounds are otherwise deliberately not compared. The built-in's veto+clamp lets
a thickness finish *on* a bound; the unbounded backends run on an interior
reparametrization and finish strictly inside one. That is a documented
contract difference, so the shared comparison below uses a target whose
optimum is interior for every backend.
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


def test_trf_needs_no_feature_to_be_selectable():
    # Hand-rolled, so a standard wheel has it. If this ever starts skipping,
    # something moved it behind a cargo feature and the docs are now wrong.
    assert "trf" in available_optimizers()
    assert 'optimizer="trf"' in repr(LmConfig(optimizer="trf"))


def test_trf_finds_the_same_interior_optimum(interior):
    mf_a, _, d_a = _run(interior, "builtin")
    mf_b, rep_b, d_b = _run(interior, "trf")
    assert rep_b["backend"] == "trf"
    assert len(d_a) == len(d_b) == len(_REFERENCE)
    assert mf_b == pytest.approx(mf_a, rel=1e-5, abs=1e-12)
    for got, want in zip(d_b, _REFERENCE):
        assert got == pytest.approx(want, abs=1e-2)


def test_trf_reaches_a_bound_active_optimum_without_landing_on_the_bound():
    """The case TRF exists for, and the one visible contract difference.

    A single low-index film on glass improves an AR demand monotonically up
    to its quarter wave (~103 nm here), so a 50 nm clamp puts the optimum on
    the clamp. Both bounded backends must find it; only the built-in is
    allowed to return the bound itself. TRF's iterates are strictly interior
    by construction -- that is what keeps the Coleman-Li scaling
    differentiable -- so it stops an ulp short, and a caller testing
    ``x == ub`` would be surprised. The distance is asserted, not hand-waved.
    """
    clamp = 50.0
    tc = TargetCollection()
    tc.add(SpectralTarget(_WL, np.zeros(_WL.size), np.full(_WL.size, 0.01),
                          0.0, "s", "R", kind="e"))
    spec = build_merit_spec(tc)

    def run(optimizer):
        ctx = SmatrixContext(spec, _ANGLES, _WL, _NO_REMOVAL_MIN, clamp,
                             LmConfig(optimizer=optimizer, max_iterations=2000))
        layers = [(MaterialSpec("Konstant", dict(n=1.38)), 10.0)]
        stack, _ = stack_from_layers(layers, _WL, {}, names=["L"])
        mf, _ = ctx.optimize_thicknesses_report(stack)
        return mf, [f["thickness"] for f in stack.films()]

    mf_b, d_b = run("builtin")
    mf_t, d_t = run("trf")
    assert len(d_b) == len(d_t) == 1
    assert mf_t == pytest.approx(mf_b, rel=1e-9, abs=1e-14)
    assert d_t[0] == pytest.approx(d_b[0], abs=1e-4)
    # The optimum really is the bound, or this test proves nothing.
    assert d_b[0] == pytest.approx(clamp, abs=1e-9)
    # ...and TRF got there from strictly inside.
    assert 0.0 < d_t[0] < clamp
    assert clamp - d_t[0] < 1e-9


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
