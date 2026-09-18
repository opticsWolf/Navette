# SPDX-License-Identifier: LGPL-3.0-or-later
"""0.7.13 (C7): the incoherent cascade against something that is not itself.

The regression twins of ``validation/review/incoherent_check.py``. Everything
that pinned the incoherent path before C7 pinned it against a copy of itself:
``validation/parity/smatrix/refs/loom_matrix.py`` is a port of the same block
sweep citing the same paper, ``intensity_path_matches_full_path_bitwise``
compares two routes through one algorithm, and ``test_physics_mirror.py``
translates the Rust tests through the Python API. All useful; none of them
evidence that the physics is right.

Three oracles here, none of which know how Navette works:

* **A closed form.** A lossless slab in air with both surfaces incoherent has
  ``R = 2*R1/(1 + R1)`` -- the geometric sum over internal reflections with
  the phase discarded. One line of algebra.
* **Thickness independence.** The only thing an incoherent layer's thickness
  feeds is a phase that has been averaged away, so R must be *bit*-identical
  across three decades of thickness.
* **The phase average.** The definition: a layer is incoherent when its
  round-trip phase is uniformly distributed, so the flagged answer IS the
  coherent answer averaged over one phase period. Sweeping the layer's
  thickness across exactly one period with the flag OFF assumes nothing about
  the block sweep at all.

Each check carries a control that fails it, because an oracle no stack can
violate measures nothing.

The averages use 64 samples, not thousands. The coherent answer is an
analytic periodic function of the round-trip phase, so an equispaced Riemann
sum over it converges geometrically: 8 samples land at 1.7e-07 and 16 already
at 5e-15. 64 is machine precision with room to spare, and costs ~6 ms.

Two traps are pinned as tests in their own right, because both were live
mistakes and both look like engine bugs when you hit them:
``test_the_absorbing_case_needs_kd_held_fixed`` and
``test_averaging_dop_is_not_averaging_the_stokes_vector``.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from navette.smatrix.smatrix import CoherenceMode, Request, ScatterMatrix

LAM = 550.0
# air / 95 nm H / 50 um slab / 110 nm L / air, the slab flagged
IDX = [1.0, 2.35, 1.52, 1.46, 1.0]
THICK = [0.0, 95.0, 50_000.0, 110.0, 0.0]
FLAGS = [0, 0, 1, 0, 0]
SLAB = 2
STOKES = Request.S0_R | Request.S1_R | Request.S2_R | Request.S3_R
SAMPLES = 64


def _solve(idx, thick, req, *, flags=None, mode=CoherenceMode.FRONT_BLOCK,
           lam=LAM, angle=0.0):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # C3/C5 flag chatter is not the subject
        sm = ScatterMatrix(
            np.asarray(idx, dtype=complex), list(thick),
            wavelengths=np.asarray([lam], float), angles=[angle],
            incoherent_flags=flags, coherence_mode=mode,
        )
    out = sm.compute(req, squeeze=False)
    return {k: float(np.asarray(v, float).ravel()[0]) for k, v in out.items()}


def _averaged(idx, thick, req, *, angle=0.0, lam=LAM, hold_kd=False,
              n_samples=SAMPLES):
    """The coherent answer averaged over one round-trip phase period."""
    n_slab = complex(idx[SLAB])
    d0 = thick[SLAB]
    cos_t = np.sqrt(1.0 - (np.sin(np.radians(angle)) / n_slab.real) ** 2 + 0j)
    period = float(lam / (2.0 * np.real(n_slab.real * cos_t)))

    acc = None
    for i in range(n_samples):
        d = d0 + period * i / n_samples
        th = list(thick)
        th[SLAB] = d
        ix = list(idx)
        if hold_kd and n_slab.imag != 0.0:
            ix[SLAB] = complex(n_slab.real, n_slab.imag * d0 / d)
        vals = _solve(ix, th, req, lam=lam, angle=angle)
        acc = vals if acc is None else {k: acc[k] + vals[k] for k in acc}
    return {k: v / n_samples for k, v in acc.items()}


# ---- A. the closed form ----------------------------------------------------
def _closed_form(n_s):
    r1 = ((1.0 - n_s) / (1.0 + n_s)) ** 2
    return 2.0 * r1 / (1.0 + r1), (1.0 - r1) / (1.0 + r1)


@pytest.mark.parametrize(
    "mode", [CoherenceMode.FRONT_BLOCK, CoherenceMode.COHERENCY_MATRIX])
def test_a_lossless_incoherent_slab_matches_the_geometric_sum(mode):
    want_r, want_t = _closed_form(1.5)
    got = _solve([1.0, 1.5, 1.0], [0.0, 1e5, 0.0], Request.RS | Request.TS,
                 flags=[0, 1, 0], mode=mode)
    assert got["Rs"] == pytest.approx(want_r, abs=1e-12)
    assert got["Ts"] == pytest.approx(want_t, abs=1e-12)


def test_a_coherent_slab_does_not_satisfy_the_closed_form():
    """The control. Without it the check above would pass on any R near 0.077."""
    want_r, _ = _closed_form(1.5)
    got = _solve([1.0, 1.5, 1.0], [0.0, 1e5, 0.0], Request.RS)
    assert abs(got["Rs"] - want_r) > 1e-6


# ---- B. thickness independence ---------------------------------------------
def test_thickness_independence_is_bit_exact():
    """Not a tolerance: there is no numerical reason for a single bit to move.

    The thickness of an incoherent layer feeds exactly one quantity, the phase,
    and that has been averaged away. 10 um to 3.7 mm.
    """
    vals = [_solve([1.0, 1.5, 1.0], [0.0, d, 0.0], Request.RS,
                   flags=[0, 1, 0])["Rs"]
            for d in (10e3, 37e3, 150e3, 800e3, 3.7e6)]
    assert len({v.hex() for v in vals}) == 1, vals


def test_a_coherent_slab_fringes_over_the_same_span():
    """The control: sub-wavelength steps must move a coherent answer."""
    vals = [_solve([1.0, 1.5, 1.0], [0.0, d, 0.0], Request.RS)["Rs"]
            for d in (10e3, 10e3 + 91.0, 10e3 + 183.0)]
    assert max(vals) - min(vals) > 1e-3


# ---- C. the phase average, on R and T --------------------------------------
@pytest.mark.parametrize("angle", [0.0, 45.0])
@pytest.mark.parametrize("key", ["Rs", "Tp"])
def test_the_flagged_answer_is_the_coherent_answer_averaged(angle, key):
    req = Request.RS | Request.TP
    flagged = _solve(IDX, THICK, req, flags=FLAGS, angle=angle)
    avg = _averaged(IDX, THICK, req, angle=angle)
    assert flagged[key] == pytest.approx(avg[key], abs=1e-12)


def test_the_phase_average_converges_geometrically():
    """Why 64 samples is not a fudge.

    The coherent answer is analytic and periodic in the round-trip phase, so
    an equispaced Riemann sum over it converges faster than any power of the
    sample count. If this ever degrades to algebraic convergence, something
    has put a discontinuity into the sweep and the averages above are no
    longer the identity they claim to be.
    """
    req = Request.RS
    flagged = _solve(IDX, THICK, req, flags=FLAGS, angle=45.0)["Rs"]
    coarse = abs(_averaged(IDX, THICK, req, angle=45.0, n_samples=8)["Rs"] - flagged)
    fine = abs(_averaged(IDX, THICK, req, angle=45.0, n_samples=16)["Rs"] - flagged)
    assert coarse > 1e-9, "8 samples should not already be exact"
    assert fine < 1e-13, f"doubling to 16 should reach machine precision, got {fine:.1e}"


def test_the_absorbing_case_needs_kd_held_fixed():
    """A trap, pinned because the naive version looks like an engine bug.

    Sweeping ``d`` to turn the phase also sweeps ``tau = exp(-2*Im(beta))``,
    so the naive average is over a stack whose absorption changes underneath
    it -- it disagrees at 3.4e-04 and is measuring itself. Holding ``k*d``
    invariant while the phase turns drops the residual to ~1e-07, the
    remainder being the slab's own Fresnel coefficients moving as ``k`` is
    rescaled. Both bounds are asserted: the fix has to work AND the naive
    version has to fail, or the next reader will delete the correction.
    """
    idx = [1.0, 2.35, complex(1.52, 0.002), 1.46, 1.0]
    req = Request.RS | Request.TP
    flagged = _solve(idx, THICK, req, flags=FLAGS)
    naive = _averaged(idx, THICK, req)
    held = _averaged(idx, THICK, req, hold_kd=True)
    assert abs(flagged["Tp"] - naive["Tp"]) > 1e-5, "the trap has stopped biting"
    for key in ("Rs", "Tp"):
        assert flagged[key] == pytest.approx(held[key], abs=5e-6)


# ---- D. the phase average, on the Stokes vector -- the C2 check -------------
@pytest.mark.parametrize("key", ["S0_R", "S1_R", "S2_R", "S3_R"])
def test_the_coherency_matrix_stokes_vector_is_the_phase_average(key):
    """Mode B against the definition, not against another implementation."""
    b = _solve(IDX, THICK, STOKES, flags=FLAGS,
               mode=CoherenceMode.COHERENCY_MATRIX, angle=45.0)
    avg = _averaged(IDX, THICK, STOKES, angle=45.0)
    assert b[key] == pytest.approx(avg[key], abs=1e-12)


def test_averaging_dop_is_not_averaging_the_stokes_vector():
    """The other trap, and the physics it hides.

    ``S0..S3`` are bilinear in the fields, so incoherent superposition
    averages *them*. ``DOP = sqrt(S1^2+S2^2+S3^2)/S0`` is a nonlinear function
    OF them, and averaging it instead is a different quantity: this stack is
    non-depolarizing, so every coherent sample has DOP = 1 exactly and their
    mean is 1, while the averaged Stokes vector has DOP = 0.994186. The gap is
    not an error -- partial depolarization is precisely what incoherent
    superposition produces, and reproducing it is the point.

    Pinned because the first draft of the C7 harness averaged DOP and reported
    a 5.8e-03 "disagreement" that was entirely its own.
    """
    mean_of_dop = _averaged(IDX, THICK, Request.DOP_R, angle=45.0)["DOP_R"]
    assert mean_of_dop == pytest.approx(1.0, abs=1e-12)

    avg = _averaged(IDX, THICK, STOKES, angle=45.0)
    dop_of_mean = float(
        np.hypot(np.hypot(avg["S1_R"], avg["S2_R"]), avg["S3_R"]) / avg["S0_R"])
    assert dop_of_mean == pytest.approx(0.994186094, abs=1e-8)

    b = _solve(IDX, THICK, Request.DOP_R, flags=FLAGS,
               mode=CoherenceMode.COHERENCY_MATRIX, angle=45.0)
    assert b["DOP_R"] == pytest.approx(dop_of_mean, abs=1e-12), (
        "the engine must depolarize by the amount the definition demands")


def test_what_front_block_refuses_does_not_satisfy_the_average():
    """C2, measured from the definition rather than from a disagreement.

    Mode A's cross channel comes from the front block while its intensities
    are totals. 0.7.10 made it refuse a cross observable on a flagged stack;
    reaching past that door through the raw engine shows what it is refusing.
    Its S2/S3 miss the phase average by 2.3e-02 and 3.9e-03 -- while its S0
    and S1 pass the very same average to 1e-16. That split is the finding.
    """
    import navette._smatrix as native

    with pytest.raises(ValueError, match="front_block"):
        _solve(IDX, THICK, Request.DOP_R, flags=FLAGS,
               mode=CoherenceMode.FRONT_BLOCK, angle=45.0)

    cache = []
    for n in IDX:
        cache += [complex(n).real, complex(n).imag]
    raw = native.core_engine(
        np.array([LAM]), np.array([np.sin(np.radians(45.0))]), len(IDX),
        np.array(cache), np.array(THICK, float),
        np.array(FLAGS, np.int32), np.zeros(len(IDX), np.int32),
        np.zeros(len(IDX)), int(CoherenceMode.FRONT_BLOCK), int(STOKES),
    )
    a = {k: float(np.asarray(v, float).ravel()[0]) for k, v in raw.items()}
    avg = _averaged(IDX, THICK, STOKES, angle=45.0)

    assert abs(a["S2_R"] - avg["S2_R"]) > 1e-3
    assert abs(a["S3_R"] - avg["S3_R"]) > 1e-3
    for key in ("S0_R", "S1_R"):
        assert a[key] == pytest.approx(avg[key], abs=1e-12), (
            f"mode A's {key} is a total and must still satisfy the average")
