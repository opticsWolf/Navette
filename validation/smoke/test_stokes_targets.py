# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Stokes / DOP targets in synthesis (W8.6): ``S0R``…``S3R``, ``DOPR``,
``S0T``…``S3T``, ``DOPT`` for an input state ``stokes_input``.

The comparison against the engine's own ``S0_R``…``DOP_T`` (and against the
coherency matrix), the merit rows against finite differences and the needle
fold against the merit slope are cargo tests; this file is the end-to-end
part through the Python surface.
"""

import numpy as np
import pytest

from navette._smatrix import LmConfig, SmatrixContext
from navette.materials import MaterialSpec
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec, sim_curves_from_arrays
from navette.synthesis.pipeline import stack_from_layers

_WL = np.linspace(480.0, 640.0, 17)
_LABELS = ["S0R", "S1R", "S2R", "S3R", "DOPR", "S0T", "S1T", "S2T", "S3T", "DOPT"]
_Q = (0.3, 0.5, 0.2)


def _spec(label, wl, values, tol, angle_deg, q=_Q, phase=False, pol="u", spectral=None):
    tc = TargetCollection()
    tc.add(SpectralTarget(wl, values, np.full(wl.size, tol), pol, spectral or label,
                          kind="e", phase=phase, angle_deg=angle_deg, stokes_input=q))
    return build_merit_spec(tc)


@pytest.mark.parametrize("label", _LABELS)
def test_stokes_labels_refuse_what_they_cannot_mean(label):
    z = np.zeros(_WL.size)
    with pytest.raises(ValueError, match="needs stokes_input"):
        _spec(label, _WL, z, 1.0, 45.0, q=None)
    with pytest.raises(ValueError, match='"u"'):
        _spec(label, _WL, z, 1.0, 45.0, pol="p")
    with pytest.raises(ValueError, match="phase=True"):
        _spec(label, _WL, z, 1.0, 45.0, phase=True)
    with pytest.raises(ValueError):
        _spec(label, _WL, z, 1.0, 45.0, q=(1.0, 1.0, 0.0))  # |q| > 1


def test_other_labels_take_no_stokes_input():
    with pytest.raises(ValueError, match="takes no stokes_input"):
        _spec(None, _WL, np.zeros(_WL.size), 1.0, 45.0, spectral="R")


def _hand(label, xs, xp, zs, zp, q):
    a, b = (1.0 + q[0]) / 2.0, (1.0 - q[0]) / 2.0
    sigma = -1.0 if label.endswith("R") else 1.0
    m = np.sqrt(xs * xp)
    d = np.angle(zp * np.conj(zs))
    s0, s1 = a * xp + b * xs, a * xp - b * xs
    s2 = sigma * m * (q[1] * np.cos(d) - q[2] * np.sin(d))
    s3 = sigma * m * (q[1] * np.sin(d) + q[2] * np.cos(d))
    q_ = label[:-1]
    if q_ == "DOP":
        return np.minimum(1.0, np.sqrt(s1 ** 2 + s2 ** 2 + s3 ** 2) / (s0 + 1e-20))
    return {"S0": s0, "S1": s1, "S2": s2, "S3": s3}[q_]


@pytest.mark.parametrize("q", [_Q, (0.0, 1.0, 0.0), (-0.6, 0.0, 0.8)])
def test_the_forward_is_the_stokes_vector_of_the_rows(q):
    # A target equal to the hand formula on hand rows leaves zero residuals.
    rng = np.random.default_rng(11)
    rows, ints = {}, {}
    for k in ("Rs", "Rp", "Ts", "Tp"):
        rows[k] = rng.uniform(0.3, 0.9, _WL.size) * np.exp(1j * rng.uniform(-3, 3, _WL.size))
        ints[k] = rng.uniform(0.05, 0.9, _WL.size)
    sim = sim_curves_from_arrays(_WL, {k: v.reshape(1, -1) for k, v in ints.items()},
                                 {k: v.reshape(1, -1) for k, v in rows.items()},
                                 angles_deg=[45.0])
    for label in _LABELS:
        side = label[-1]
        want = _hand(label, ints[side + "s"], ints[side + "p"],
                     rows[side + "s"], rows[side + "p"], q)
        got = np.asarray(_spec(label, _WL, want, 1.0, 45.0, q=q).residuals(sim))
        np.testing.assert_allclose(got, 0.0, atol=1e-10, err_msg=f"{label} {q}")


_REFERENCE = [118.0, 203.0, 64.0]
_START = [121.0, 199.0, 66.0]


def _ctx(spec, wl, jacobian="analytic"):
    return SmatrixContext(spec, wl, 1e-9, 1000.0, LmConfig(jacobian=jacobian),
                          angles_deg=np.array([60.0]))


def _stack(ds, wl):
    layers = [(MaterialSpec("Konstant", dict(n=n)), d)
              for n, d in zip([2.30, 1.46, 2.30], ds)]
    st, _ = stack_from_layers(layers, wl, {}, names=[f"L{i}" for i in range(len(ds))])
    return st


@pytest.mark.parametrize("label", ["S1R", "S2R", "S3T", "DOPR"])
def test_a_stokes_target_takes_the_analytic_path_to_the_reference(label):
    # Targets read off a reference stack put a merit-zero optimum a few nm
    # from the start; both Jacobian modes must find it.
    flat = _spec(label, _WL, np.zeros(_WL.size), 1.0, 60.0)
    value = np.asarray(flat.residuals(_ctx(flat, _WL).simulate(_stack(_REFERENCE, _WL))))
    spec = _spec(label, _WL, value, 0.01, 60.0)
    for jacobian in ("analytic", "fd"):
        st = _stack(_START, _WL)
        mf, report = _ctx(spec, _WL, jacobian).refine_thicknesses_report(st)
        assert (report["analytic_jacobians"] > 0) == (jacobian == "analytic")
        assert mf < 1e-10, (jacobian, mf)
        for got, want in zip([f["thickness"] for f in st.films()], _REFERENCE):
            assert got == pytest.approx(want, abs=1e-3)
