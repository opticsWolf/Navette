# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Synthesis refuses, up front, a demand its evaluator cannot score (W3).

The synthesis evaluator fills the front s/p intensities (Rs, Rp, Ts, Tp)
and the complex forward t (Ts/Tp phase, PDts/PDtp). A demand on any other
row (unpolarized, back incidence, front reflection phase) used to reach
the merit as a constant 1e6 penalty and the LM as a mid-run "missing
curve". Every door that builds a synthesis context now refuses it at
construction, naming each offending curve once with its angles. The
check lives in Rust (``MeritSpec::check_supported``); these tests pin the
Python doors.
"""
import numpy as np
import pytest

from navette import _smatrix as S
from navette.spectralweave.target import AngularTarget, SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec
from navette.synthesis.pipeline import PipelineConfig, run_needle

WL = np.array([500.0, 600.0, 700.0])
ANG = np.array([45.0])
V, T = np.full(3, 0.1), np.full(3, 0.01)
MATERIALS = {k: np.full(WL.shape, v) for k, v in {"L": 1.46 + 0j, "H": 2.35 + 0j}.items()}
SEED = [(MATERIALS["L"], 100.0), (MATERIALS["H"], 80.0)]


def _spec(*targets):
  tc = TargetCollection()
  for t in targets:
    tc.add(t)
  return build_merit_spec(tc)


def _ctx(*targets):
  return S.SmatrixContext(_spec(*targets), WL, angles_deg=ANG)


@pytest.mark.parametrize("make", [
  lambda: SpectralTarget(WL, V, T, "s", "R", angle_deg=45.0),
  lambda: SpectralTarget(WL, V, T, "p", "T", angle_deg=45.0),
  lambda: SpectralTarget(WL, V, T, "s", "A", angle_deg=45.0),
  lambda: SpectralTarget(WL, np.full(3, 0.5), np.full(3, 0.1), "p", "T", angle_deg=45.0,
                         phase=True),
  lambda: SpectralTarget(WL, np.full(3, 0.5), np.full(3, 0.1), "s", "PDts", angle_deg=45.0,
                         phase=True),
], ids=["Rs", "Tp", "As", "Tp-phase", "PDts"])
def test_supported_demands_construct(make):
  _ctx(make())


@pytest.mark.parametrize("make, line", [
  (lambda: SpectralTarget(WL, V, T, "u", "R", angle_deg=45.0), "Ru intensity at 45 deg"),
  (lambda: SpectralTarget(WL, V, T, "u", "A", angle_deg=45.0), "Au intensity at 45 deg"),
  (lambda: SpectralTarget(WL, V, T, "s", "RB", angle_deg=45.0), "RBs intensity at 45 deg"),
  (lambda: SpectralTarget(WL, np.full(3, 0.5), np.full(3, 0.1), "s", "R", angle_deg=45.0,
                          phase=True), "Rs phase at 45 deg"),
], ids=["Ru", "Au", "RBs", "Rs-phase"])
def test_unsupported_demands_are_refused_by_name(make, line):
  with pytest.raises(ValueError, match="synthesis cannot score 1 demand") as exc:
    _ctx(make())
  assert line in str(exc.value)


def test_an_angular_demand_is_one_line():
  angles = np.array([0.0, 20.0, 40.0, 60.0])
  tc = TargetCollection()
  tc.add(AngularTarget(600.0, np.full(4, 0.1), np.full(4, 0.01), "u", "T", angles_deg=angles))
  with pytest.raises(ValueError) as exc:
    S.SmatrixContext(build_merit_spec(tc), WL, angles_deg=angles)
  msg = str(exc.value)
  assert "cannot score 1 demand" in msg
  assert "Tu intensity at 0 deg, 20 deg, 40 deg, 60 deg" in msg


def test_one_bad_demand_refuses_the_whole_set():
  with pytest.raises(ValueError, match="Ru intensity"):
    _ctx(SpectralTarget(WL, V, T, "s", "R", angle_deg=45.0),
         SpectralTarget(WL, V, T, "u", "R", angle_deg=45.0))


def test_run_needle_refuses_before_solving():
  tc = TargetCollection()
  tc.add(SpectralTarget(WL, V, T, "u", "R", angle_deg=45.0))
  calls = []
  with pytest.raises(ValueError, match="Ru intensity at 45 deg"):
    run_needle(layers=SEED, names=["L", "H"], targets=tc, wavelengths=WL,
               materials=MATERIALS, substrate=(1.52, "sub"), angles_deg=ANG,
               pipeline_config=PipelineConfig(max_macro_cycles=1, needles_per_cycle=1),
               callback=lambda cycle, phase: calls.append(cycle))
  assert calls == []
