# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Phase-dispersion targets in synthesis (W8.5): ``GDR/T``, ``GDDR/T``,
``TODR/T``, ``FODR/T`` × s/p.

The n-th ω derivative (fs^n) of the unwrapped phase of r_f / t_f over the
simulated grid, the engine's ``GD``…``FOD``. The comparison against the
engine's own maps, the merit rows against finite differences and the needle
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
_LABELS = {"GDR": 1, "GDT": 1, "GDDR": 2, "GDDT": 2,
           "TODR": 3, "TODT": 3, "FODR": 4, "FODT": 4}
_C_NM_PER_FS = 299.792458


def _spec(label, wl, values, tol, angle_deg, phase=True, pol="s"):
    tc = TargetCollection()
    tc.add(SpectralTarget(wl, values, np.full(wl.size, tol), pol, label,
                          kind="e", phase=phase, angle_deg=angle_deg))
    return build_merit_spec(tc)


@pytest.mark.parametrize("label", list(_LABELS))
def test_dispersion_labels_take_s_or_p_and_need_phase(label):
    with pytest.raises(ValueError, match='"s" or "p"'):
        _spec(label, _WL, np.zeros(_WL.size), 1.0, 45.0, pol="u")
    with pytest.raises(ValueError, match="phase=True"):
        _spec(label, _WL, np.zeros(_WL.size), 1.0, 45.0, phase=False)


def test_the_forward_is_the_unwrapped_phase_differentiated_in_omega():
    # A target equal to numpy's gradient chain of a hand row leaves zero
    # residuals: unwrap, then np.gradient n times against omega.
    rng = np.random.default_rng(7)
    omega = 2.0 * np.pi * _C_NM_PER_FS / _WL
    rows = {}
    for k in ("Rs", "Rp", "Ts", "Tp"):
        phase = np.cumsum(rng.uniform(1.0, 2.5, _WL.size))  # wraps several times
        rows[k] = rng.uniform(0.3, 0.9, _WL.size) * np.exp(1j * phase)
    sim = sim_curves_from_arrays(_WL, {}, {k: v.reshape(1, -1) for k, v in rows.items()},
                                 angles_deg=[45.0])
    for label, order in _LABELS.items():
        for pol in "sp":
            g = np.unwrap(np.angle(rows[label[-1] + pol]))
            for _ in range(order):
                g = np.gradient(g, omega)
            got = np.asarray(_spec(label, _WL, g, 10.0 ** order, 45.0, pol=pol).residuals(sim))
            np.testing.assert_allclose(got, 0.0, atol=1e-8, err_msg=f"{label}/{pol}")


_REFERENCE = [118.0, 203.0, 64.0]
_START = [121.0, 199.0, 66.0]


def _ctx(spec, wl, jacobian="analytic"):
    return SmatrixContext(spec, wl, 1e-9, 1000.0, LmConfig(jacobian=jacobian),
                          angles_deg=np.array([45.0]))


def _stack(ds, wl):
    layers = [(MaterialSpec("Konstant", dict(n=n)), d)
              for n, d in zip([2.30, 1.46, 2.30], ds)]
    st, _ = stack_from_layers(layers, wl, {}, names=[f"L{i}" for i in range(len(ds))])
    return st


@pytest.mark.parametrize("label, pol", [("GDR", "p"), ("GDT", "p"), ("GDDR", "s"),
                                        ("GDDT", "p")])
def test_a_dispersion_target_takes_the_analytic_path_to_the_reference(label, pol):
    # Targets read off a reference stack put a merit-zero optimum a few nm
    # from the start; both Jacobian modes must find it. (GDT/s stalls at a
    # local minimum from _START, merit 2.9e-4, analytic and FD alike.)
    flat = _spec(label, _WL, np.zeros(_WL.size), 1.0, 45.0, pol=pol)
    value = np.asarray(flat.residuals(_ctx(flat, _WL).simulate(_stack(_REFERENCE, _WL))))
    spec = _spec(label, _WL, value, 0.01 * 10.0 ** _LABELS[label], 45.0, pol=pol)
    for jacobian in ("analytic", "fd"):
        st = _stack(_START, _WL)
        mf, report = _ctx(spec, _WL, jacobian).refine_thicknesses_report(st)
        assert (report["analytic_jacobians"] > 0) == (jacobian == "analytic")
        assert mf < 1e-10, (jacobian, mf)
        for got, want in zip([f["thickness"] for f in st.films()], _REFERENCE):
            assert got == pytest.approx(want, abs=1e-3)
