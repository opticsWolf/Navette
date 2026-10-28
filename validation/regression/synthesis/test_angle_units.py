# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""`angles_in_radians` on every Python door that takes an angle.

Degrees are the default everywhere; ``angles_in_radians=True`` takes the
same angles in radians. A solver grid in radians goes straight into
``sin`` (or ``cos`` for the differential-phase reference); an angle that
only labels something (a target, a merit key, a ``SimCurves`` row) is
stored in degrees, the merit-key unit. The seams pinned here:

  * radians and degrees describe the same problem: every door agrees with
    its degree twin to rounding, and bit for bit at normal incidence;
  * the radians reach the solver as radians, not read as degrees: a
    radian grid matches the engine's own radian door;
  * leaving the flag out is the degree run bit for bit.

The degree path itself is pinned by the forward golden baseline and the
needle fingerprint.
"""
import numpy as np
import pytest

from navette import _smatrix as S
from navette.smatrix.smatrix import Request, ScatterMatrix
from navette.spectralweave.target import (AngularTarget, ColorTarget,
                                          SpectralTarget, TargetCollection)
from navette.synthesis import (apply_reference_rotation, build_merit_spec,
                               sim_curves_from_arrays)
from navette.synthesis.pipeline import PipelineConfig, run_needle

WL = np.linspace(450.0, 750.0, 31)
DEG = np.array([0.0, 30.0, 60.0])
RAD = np.radians(DEG)
N = {"L": 1.46 + 0j, "H": 2.35 + 0j}
MATERIALS = {k: np.full(WL.shape, v) for k, v in N.items()}
SEED = [(MATERIALS["L"], 300.0), (MATERIALS["H"], 200.0)]
CURVES = ("Rs", "Rp", "Ts", "Tp")


def _stack():
  films = [S.LayerSpec(name, np.full(WL.shape, n), d)
           for name, n, d in (("H", N["H"], 120.0), ("L", N["L"], 180.0))]
  return S.DesignStack(S.LayerSpec("air", np.ones(WL.size, complex), 0.0),
                       S.LayerSpec("sub", np.full(WL.shape, 1.52 + 0j), 0.0), films)


def _readback_spec(angles_deg):
  """Neutral targets: residuals are the simulated curves themselves."""
  spec = S.MeritSpec()
  zeros, ones = np.zeros(WL.size), np.ones(WL.size)
  for a in angles_deg:
    for c in CURVES:
      spec.add_target(spec.add_key(float(a), c), WL, zeros, ones, "e", "linear", 1.0)
  return spec


def _curves(spec, angles, **kw):
  sim = S.SmatrixContext(spec, np.asarray(angles, float), WL, **kw).simulate(_stack())
  return np.asarray(spec.residuals(sim)).reshape(-1, len(CURVES), WL.size)


def test_a_radian_grid_reads_each_angle_as_the_engine_does():
  got = _curves(_readback_spec(DEG), RAD, angles_in_radians=True)
  sm = ScatterMatrix(np.array([1.0, N["H"], N["L"], 1.52]), [0.0, 120.0, 180.0, 0.0],
                     wavelengths=WL, angles=RAD, angles_in_radians=True)
  eng = sm.compute(Request.RS | Request.RP | Request.TS | Request.TP, squeeze=False)
  for ci, c in enumerate(CURVES):
    assert np.max(np.abs(got[:, ci] - eng[c])) < 1e-13, c


def test_a_radian_grid_matches_its_degree_twin():
  spec = _readback_spec(DEG)
  deg = _curves(spec, DEG)
  rad = _curves(spec, RAD, angles_in_radians=True)
  assert np.max(np.abs(rad - deg)) < 1e-14
  # Normal incidence is 0 in both units: nothing to round.
  assert rad[0].tobytes() == deg[0].tobytes()


def test_the_default_is_degrees():
  spec = _readback_spec(DEG)
  assert _curves(spec, DEG).tobytes() == _curves(spec, DEG, angles_in_radians=False).tobytes()


def test_merit_keys_and_sim_rows_take_radians():
  spec_deg, spec_rad = _readback_spec(DEG), S.MeritSpec()
  zeros, ones = np.zeros(WL.size), np.ones(WL.size)
  for a in RAD:
    for c in CURVES:
      spec_rad.add_target(spec_rad.add_key(float(a), c, angles_in_radians=True),
                          WL, zeros, ones, "e", "linear", 1.0)
  rows = {c: np.linspace(0.1, 0.9, DEG.size * WL.size) for c in CURVES}
  sim_deg = sim_curves_from_arrays(DEG, WL, rows)
  sim_rad = sim_curves_from_arrays(RAD, WL, rows, angles_in_radians=True)
  ref = np.asarray(spec_deg.residuals(sim_deg))
  for spec in (spec_deg, spec_rad):
    for sim in (sim_deg, sim_rad):
      assert np.asarray(spec.residuals(sim)).tobytes() == ref.tobytes()


def test_targets_take_radians_and_store_degrees():
  vals, tol = np.zeros_like(WL), np.full_like(WL, 0.01)
  st = SpectralTarget(WL, vals, tol, np.pi / 6, "s", "R", angles_in_radians=True)
  assert st.angle == pytest.approx(30.0, abs=1e-12) and st.angles_in_radians is False
  at = AngularTarget(550.0, RAD, np.zeros(3), np.full(3, 0.01), "s", "R",
                     angles_in_radians=True)
  assert np.allclose(at.angles, DEG, atol=1e-12, rtol=0) and at.angles_in_radians is False
  ct = ColorTarget(angle=np.pi / 4, angles_in_radians=True)
  assert ct.angle == pytest.approx(45.0, abs=1e-12) and ct.angles_in_radians is False
  # Left out, the flag changes nothing.
  assert SpectralTarget(WL, vals, tol, 30.0, "s", "R").angle == 30.0


def test_the_reference_rotation_takes_radians():
  amp = np.exp(1j * np.linspace(0.0, 3.0, WL.size))[None, :]
  n_inc = np.ones(WL.size)
  for d, r in zip(DEG, RAD):
    deg = apply_reference_rotation(amp, WL, d, n_inc, 400.0)
    rad = apply_reference_rotation(amp, WL, r, n_inc, 400.0, angles_in_radians=True)
    assert np.max(np.abs(rad - deg)) < 1e-13
  # cos(pi/2) is not 0 in floating point, but it is cos of the radians given.
  rot = np.asarray(S.reference_rotation(WL, np.pi / 2, n_inc, 400.0, 1.0,
                                        angles_in_radians=True))
  ref = np.exp(-1j * 2 * np.pi * 400.0 * np.cos(np.pi / 2) / WL)
  assert np.max(np.abs(rot - ref)) < 1e-15


def test_the_needle_fold_takes_radians():
  spec = build_merit_spec(_targets(DEG))
  deg = S.build_needle_targets(spec, DEG, WL)
  rad = S.build_needle_targets(spec, RAD, WL, angles_in_radians=True)
  for k in ("r", "t"):
    for part in ("targets", "weights"):
      assert np.asarray(rad[k][part]).tobytes() == np.asarray(deg[k][part]).tobytes(), k


def _targets(angles_deg):
  tc = TargetCollection()
  for a in angles_deg:
    tc.add(SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, 0.01), float(a), "s", "R"))
  return tc


def _run(angles, **kw):
  return run_needle(layers=SEED, names=["L", "H"], targets=_targets(DEG),
                    wavelengths=WL, materials=MATERIALS, substrate=(1.52, "sub"),
                    angles_deg=angles,
                    pipeline_config=PipelineConfig(max_macro_cycles=2, needles_per_cycle=1),
                    **kw)


def test_run_needle_takes_radians():
  base = _run(DEG)
  rad = _run(RAD, angles_in_radians=True)
  assert rad["final_layer_count"] == base["final_layer_count"]
  assert rad["final_mf"] == pytest.approx(base["final_mf"], rel=1e-9)
  # At normal incidence the two runs are one run.
  assert (_run([0.0], angles_in_radians=True)["final_mf"].hex()
          == _run([0.0])["final_mf"].hex())


def test_the_environments_door_takes_radians():
  kw = dict(targets=_targets(DEG), wavelengths=WL,
            design={"coat": [(N["L"], 300.0, "L"), (N["H"], 200.0, "H")]},
            environments=[{"name": "bare", "stack": [{"design": "coat"}]}],
            pipeline_config=PipelineConfig(max_macro_cycles=1, needles_per_cycle=1),
            substrate=(1.52 + 0j, "sub"))
  base = run_needle(angles_deg=DEG, **kw)
  rad = run_needle(angles_deg=RAD, angles_in_radians=True, **kw)
  assert rad["final_mf"] == pytest.approx(base["final_mf"], rel=1e-9)


def test_the_pipeline_class_takes_radians():
  from navette.synthesis.pipeline import stack_from_layers
  stack, pool = stack_from_layers(SEED, WL, MATERIALS, names=["L", "H"])
  spec = build_merit_spec(_targets(DEG))
  cfg = PipelineConfig(max_macro_cycles=1, needles_per_cycle=1)
  base = S.NeedlePipeline(stack, spec, DEG, WL, pool, cfg).run()
  rad = S.NeedlePipeline(stack, spec, RAD, WL, pool, cfg, angles_in_radians=True).run()
  assert rad["final_mf"] == pytest.approx(base["final_mf"], rel=1e-9)
