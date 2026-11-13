# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Intensity-ratio targets in synthesis (W8.3, W8.4): ``RatioR/T``, ``ExtR/T``,
``DiattR/T`` and the ellipsometric ``PsiR/T``.

``Ratio = s/p``, ``Ext = p/s`` and ``Diatt = (p − s)/(p + s)`` of one front
intensity, each denominator guarded by the engine's 1e-20, and
``Psi = atan √(p/s)``, the engine's ``Psi_R``/``Psi_T``. The comparison
against the engine's own maps, the merit rows against finite differences
and the needle fold against the merit slope are cargo tests; this file is
the end-to-end part through the Python surface.
"""

import warnings

import numpy as np
import pytest

from navette._smatrix import LmConfig, SmatrixContext
from navette.materials import MaterialSpec
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec, sim_curves_from_arrays
from navette.synthesis.pipeline import stack_from_layers

_WL = np.linspace(480.0, 640.0, 9)
_LABELS = ("RatioR", "RatioT", "ExtR", "ExtT", "DiattR", "DiattT", "PsiR", "PsiT")


def _spec(label, wl, values, tol, angle_deg, phase=False, pol="u"):
    tc = TargetCollection()
    tc.add(SpectralTarget(wl, values, np.full(wl.size, tol), pol, label,
                          kind="e", phase=phase, angle_deg=angle_deg,
                          normalization_mode="linear"))
    return build_merit_spec(tc)


@pytest.mark.parametrize("label", _LABELS)
def test_ratio_labels_take_u_and_refuse_phase(label):
    with pytest.raises(ValueError, match='"u"'):
        _spec(label, _WL, np.ones(_WL.size), 1.0, 45.0, pol="p")
    with pytest.raises(ValueError, match="phase"):
        _spec(label, _WL, np.ones(_WL.size), 1.0, 45.0, phase=True)


def test_the_forward_is_the_ratio_of_the_intensity_rows():
    # A target equal to the ratio of hand rows leaves zero residuals.
    rng = np.random.default_rng(5)
    rows = {k: rng.uniform(0.2, 0.9, _WL.size) for k in ("Rs", "Rp", "Ts", "Tp")}
    sim = sim_curves_from_arrays(_WL, {k: v.reshape(1, -1) for k, v in rows.items()},
                                 angles_deg=[45.0])
    for label in _LABELS:
        s, p = rows[label[-1] + "s"], rows[label[-1] + "p"]
        want = {"Ratio": s / p, "Ext": p / s, "Diatt": (p - s) / (p + s),
                "Psi": np.arctan(np.sqrt(p / s))}[label[:-1]]
        got = np.asarray(_spec(label, _WL, want, 0.01, 45.0).residuals(sim))
        np.testing.assert_allclose(got, 0.0, atol=1e-9, err_msg=label)


def test_a_demand_at_normal_incidence_only_warns():
    with pytest.warns(UserWarning, match="normal incidence"):
        _spec("ExtR", _WL, np.ones(_WL.size), 0.1, 0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _spec("ExtR", _WL, np.ones(_WL.size), 0.1, 45.0)


_REFERENCE = [118.0, 203.0, 64.0]
_START = [122.0, 198.0, 67.0]
# The reflection labels that are monotone in Rp/Rs (DiattR, PsiR) stall at
# the same local minimum from _START, analytic and FD alike (merit 716 and
# 94); they start closer.
_START_R = [120.0, 201.0, 65.5]


def _ctx(spec, wl, jacobian="analytic"):
    return SmatrixContext(spec, wl, 1e-9, 1000.0, LmConfig(jacobian=jacobian),
                          angles_deg=np.array([45.0]))


def _stack(ds, wl):
    layers = [(MaterialSpec("Konstant", dict(n=n)), d)
              for n, d in zip([2.30, 1.46, 2.30], ds)]
    st, _ = stack_from_layers(layers, wl, {}, names=[f"L{i}" for i in range(len(ds))])
    return st


@pytest.mark.parametrize("label, start", [("ExtT", _START), ("DiattT", _START),
                                          ("PsiT", _START), ("DiattR", _START_R),
                                          ("PsiR", _START_R)])
def test_a_ratio_takes_the_analytic_path_to_the_reference(label, start):
    # Targets read off a reference stack put a merit-zero optimum a few nm
    # from the start; both Jacobian modes must find it.
    flat = _spec(label, _WL, np.zeros(_WL.size), 1.0, 45.0)
    value = np.asarray(flat.residuals(_ctx(flat, _WL).simulate(_stack(_REFERENCE, _WL))))
    spec = _spec(label, _WL, value, 0.01, 45.0)
    for jacobian in ("analytic", "fd"):
        st = _stack(start, _WL)
        mf, report = _ctx(spec, _WL, jacobian).optimize_thicknesses_report(st)
        assert (report["analytic_jacobians"] > 0) == (jacobian == "analytic")
        assert mf < 1e-10, (jacobian, mf)
        for got, want in zip([f["thickness"] for f in st.films()], _REFERENCE):
            assert got == pytest.approx(want, abs=1e-3)
