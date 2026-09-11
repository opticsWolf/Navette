# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""The color gradient on the Python needle path (R4.2).

``build_needle_targets`` -> ``needle_gradient`` is the documented way to
assemble a needle cycle in Python. It carried every pointwise demand
correctly and dropped color demands on the floor.

A color demand (Lab/DE2000, White, Yellow, dominant wavelength) integrates
the whole spectrum into a *single* residual, so it has no per-point target
to fold into the ``r``/``t`` pairs. The native fold emits its analytic
``dF/dcurve`` per solver point instead -- ``NeedleTargets.grad_r/grad_t``,
deposited into the same accumulator as the pointwise terms by
``needle_pass.rs``. The Python dict omitted both arrays and
``needle_gradient`` had no argument that could have taken them, so the
folded ``r`` pair was all zeros and the returned gradient was zero.

Not a crash and not a wrong shape: a silent wrong-optimizer. ``run_design``
was never affected -- it folds and deposits inside Rust.

``validation/review/color_grad_python.py`` is the full harness (FD against
the merit, the chain rule, additivity); this file is the part CI runs.
"""

import json

import numpy as np
import pytest

from navette._smatrix import SimCurves, build_needle_targets, compile_merit_spec
from navette.smatrix.needle import NeedleRequest, needle_gradient
from navette.smatrix.smatrix import Request, ScatterMatrix

_N = np.array([1.0 + 0j, 2.35 + 0j, 1.46 + 0j, 2.10 + 0j, 1.52 + 0j])
_D = np.array([0.0, 120.0, 200.0, 80.0, 0.0])
_WLS = np.linspace(450.0, 750.0, 31)
_ANGLES = np.array([30.0])
_NW = _WLS.size

# A needle of the host's own index, inside the host: insertion == growth.
_Z = 220.0
_HOST_N = 1.46 + 0j


def _spec(curve="Rs"):
    return compile_merit_spec(json.dumps({
        "spectral": [], "angular": [],
        "color": [{
            "curve": curve, "angle": 30.0,
            "illuminant": "D65", "observer": "1931_2deg",
            "quantity": "Lab", "reference": [55.0, 12.0, -18.0],
            "distance": "DeltaE2000", "weight": 1.7,
        }],
    }))


@pytest.fixture(scope="module")
def bench():
    st = ScatterMatrix(_N, _D, wavelengths=_WLS, angles=_ANGLES)
    out = st.compute(Request.RS | Request.TS)
    rs = np.atleast_2d(np.asarray(out["Rs"], float))
    ts = np.atleast_2d(np.asarray(out["Ts"], float))

    def sim_of(rs_row=None):
        s = SimCurves(_ANGLES, _WLS, float(_D[1:-1].sum()), 1.0, 1.52)
        s.set_curve("Rs", np.ascontiguousarray(
            (rs if rs_row is None else rs_row).ravel()))
        s.set_curve("Ts", np.ascontiguousarray(ts.ravel()))
        return s

    spec = _spec()
    return {
        "stack": st, "rs": rs, "ts": ts, "sim_of": sim_of, "spec": spec,
        "fold": build_needle_targets(spec, _ANGLES, _WLS, sim_of()),
    }


def _grad(request=NeedleRequest.P, **kw):
    st = ScatterMatrix(_N, _D, wavelengths=_WLS, angles=_ANGLES)
    return needle_gradient(st, np.full(_NW, _HOST_N), [_Z], request, **kw)


# --------------------------------------------------------------------------
# 1. The fold surfaces the buckets, and they are what they claim to be
# --------------------------------------------------------------------------

def test_the_fold_emits_the_color_buckets(bench):
    fold = bench["fold"]
    assert "grads_r" in fold and "grads_t" in fold, sorted(fold)
    assert np.asarray(fold["grads_r"]).shape == (_NW,)
    assert np.any(np.asarray(fold["grads_r"]) != 0.0)
    # An Rs demand deposits nowhere else.
    assert np.all(np.asarray(fold["grads_t"]) == 0.0)


def test_the_pointwise_pair_stays_empty_for_a_color_demand(bench):
    """The bug's fingerprint: a color spec folds to all-zero r targets.

    Nothing is wrong with that -- there *is* no per-point target -- but it
    is why the old path's zero gradient looked like an honest answer.
    """
    r = bench["fold"]["r"]
    assert np.all(np.asarray(r["weights"], float) == 0.0)
    assert np.all(np.asarray(r["targets"], float) == 0.0)


def test_grads_r_equals_dF_dR(bench):
    """The fold's analytic dF/dcurve against a finite difference of the merit."""
    spec, rs, sim_of = bench["spec"], bench["rs"], bench["sim_of"]
    g = np.asarray(bench["fold"]["grads_r"], float)

    h = 1e-6
    fd = np.empty(_NW)
    for i in range(_NW):
        rp = rs.copy(); rp[0, i] += h
        rm = rs.copy(); rm[0, i] -= h
        fd[i] = (spec.merit(sim_of(rp), 1e6) - spec.merit(sim_of(rm), 1e6)) / (2 * h)

    assert np.max(np.abs(g - fd)) / np.max(np.abs(fd)) < 1e-6


# --------------------------------------------------------------------------
# 2. The regression
# --------------------------------------------------------------------------

def test_without_the_buckets_the_color_gradient_is_zero(bench):
    """The old behaviour, pinned so its return would be visible."""
    r = bench["fold"]["r"]
    out = _grad(targets_r=np.asarray(r["targets"], float),
                weights_r=np.asarray(r["weights"], float))
    assert np.all(out["P_s"] == 0.0)


def test_with_the_buckets_the_gradient_is_the_merit_derivative(bench):
    """End to end: dF/d(thickness) through the Python path vs merit FD.

    ``P`` is the half-gradient (see ``needle_gradient``'s docstring), so
    ``dF/dd = 2*sum_k P_k``; a needle of the host's own index placed inside
    the host is exactly growing that layer.
    """
    r, spec, sim_of = bench["fold"]["r"], bench["spec"], bench["sim_of"]
    out = _grad(targets_r=np.asarray(r["targets"], float),
                weights_r=np.asarray(r["weights"], float),
                grads_r=np.asarray(bench["fold"]["grads_r"], float))
    analytic = 2.0 * float(np.sum(out["P_s"]))

    def merit_at(th):
        st = ScatterMatrix(_N, np.asarray(th, float), wavelengths=_WLS,
                           angles=_ANGLES)
        rs = np.atleast_2d(np.asarray(st.compute(Request.RS)["Rs"], float))
        return spec.merit(sim_of(rs), 1e6)

    step = np.eye(_D.size)[2] * 0.01
    fd = (merit_at(_D + step) - merit_at(_D - step)) / 0.02
    assert analytic == pytest.approx(fd, rel=2e-3)
    assert abs(analytic) > 1.0, "the FD agreement above would be vacuous at zero"


def test_the_deposit_is_the_pointwise_kernel_with_g_in_place_of_the_residual(bench):
    """1e-12, not 1e-6: the two kernels differ only in one scalar factor.

    ``p_coherent_from_fields`` multiplies ``Re{conj(r) dr/ddz}`` by
    ``2*w*(R - t)``; the color kernel multiplies it by ``g``. Driving the
    pointwise channel with target 0 / weight 1 gives ``2*R*Re{...}``, so
    the deposit must equal ``g * P_ref / (2*R)`` point for point. Anything
    looser would not notice the deposit landing on the wrong depth row.
    """
    g = np.asarray(bench["fold"]["grads_r"], float)
    ref = _grad(targets_r=np.zeros(_NW), weights_r=np.ones(_NW))["P_s"].ravel()
    got = _grad(targets_r=np.zeros(_NW), weights_r=np.zeros(_NW),
                grads_r=g)["P_s"].ravel()
    want = 0.5 * g * ref / bench["rs"].ravel()
    assert np.max(np.abs(got - want)) / np.max(np.abs(want)) < 1e-12


def test_pointwise_and_color_superpose(bench):
    """One accumulator, as in ``needle_pass.rs`` -- not one overwriting the other."""
    g = np.asarray(bench["fold"]["grads_r"], float)
    pair = dict(targets_r=np.full(_NW, 0.05), weights_r=np.full(_NW, 3.0))
    only_pt = _grad(**pair)["P_s"]
    only_col = _grad(targets_r=np.zeros(_NW), weights_r=np.zeros(_NW),
                     grads_r=g)["P_s"]
    both = _grad(grads_r=g, **pair)["P_s"]
    assert np.max(np.abs(both - (only_pt + only_col))) \
        / np.max(np.abs(both)) < 1e-12


# --------------------------------------------------------------------------
# 3. Channels stay separate
# --------------------------------------------------------------------------

def test_a_transmission_color_demand_rides_the_T_channel(bench):
    fold_t = build_needle_targets(_spec("Ts"), _ANGLES, _WLS, bench["sim_of"]())
    gr = np.asarray(fold_t["grads_r"], float)
    gt = np.asarray(fold_t["grads_t"], float)
    assert np.all(gr == 0.0) and np.any(gt != 0.0)

    out = _grad(NeedleRequest.P | NeedleRequest.P_T,
                targets_r=np.zeros(_NW), weights_r=np.zeros(_NW),
                targets_t=np.zeros(_NW), weights_t=np.zeros(_NW),
                grads_r=gr, grads_t=gt)
    assert np.max(np.abs(out["P_T_s"])) > 0.0
    assert np.all(out["P_s"] == 0.0), "an R-channel deposit appeared from nowhere"


# --------------------------------------------------------------------------
# 4. A bucket with nowhere to land is an error
# --------------------------------------------------------------------------

def test_grads_r_without_P_requested_is_refused(bench):
    """Silently dropping it would be the same bug under a new name."""
    g = np.asarray(bench["fold"]["grads_r"], float)
    with pytest.raises(ValueError, match="NREQ_P was not requested"):
        _grad(NeedleRequest.P_T, grads_r=g)


def test_grads_t_without_P_T_requested_is_refused(bench):
    fold_t = build_needle_targets(_spec("Ts"), _ANGLES, _WLS, bench["sim_of"]())
    with pytest.raises(ValueError, match="NREQ_P_T was not requested"):
        _grad(NeedleRequest.P, grads_t=np.asarray(fold_t["grads_t"], float))


def test_an_all_zero_bucket_is_not_a_demand(bench):
    """A color-free fold hands both arrays through; that must stay legal."""
    out = _grad(NeedleRequest.P, grads_t=np.zeros(_NW))
    assert np.all(np.isfinite(out["P_s"]))


def test_a_non_finite_bucket_names_its_index(bench):
    g = np.asarray(bench["fold"]["grads_r"], float).copy()
    g[7] = np.nan
    with pytest.raises(ValueError, match=r"grads_r\[7\]"):
        _grad(NeedleRequest.P, grads_r=g)


def test_a_wrong_length_bucket_is_refused(bench):
    with pytest.raises(ValueError, match="grads_r"):
        _grad(NeedleRequest.P, grads_r=np.ones(_NW + 1))
