# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Every angle door names its unit: ``*_deg=`` or ``*_rad=``, exactly one.

There is no default unit. A grid is ``angles_deg=`` / ``angles_rad=``, a
single angle ``angle_deg=`` / ``angle_rad=``; both or neither is a
``ValueError`` naming the pair, and an angle outside [0, 90] deg /
[0, pi/2] rad is refused with the value and the unit it was read in. The
native layer converts, once; everything past the door is radians. The
seams pinned here:

  * the two units describe the same problem: every door agrees with its
    degree twin to rounding, and bit for bit at normal incidence;
  * a radian grid reaches the solver as radians (it matches the engine);
  * both, neither and out of range raise, on every door;
  * the old spellings (``angles=``, ``angle=``, ``angles_in_radians=``,
    the ``"angle"`` JSON key) are refused, never read as degrees;
  * target-weaver keys hold radians, and a simulation keyed otherwise
    warns instead of silently costing ``missing_penalty``.

The degree path itself is pinned by the forward golden baseline and the
needle fingerprint.
"""
import json
import warnings

import numpy as np
import pytest

from navette import _smatrix as S
from navette._spectralweave import OpticalWeaver, calculate_merit
from navette.smatrix.smatrix import Request, ScatterMatrix
from navette.spectralweave.target import (AngularTarget, ColorTarget,
                                          SpectralTarget, TargetCollection)
from navette.structure import solve_structure
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
IDX = np.array([1.0, N["H"], N["L"], 1.52])
THICK = [0.0, 120.0, 180.0, 0.0]


def _stack():
  films = [S.LayerSpec(name, np.full(WL.shape, n), d)
           for name, n, d in (("H", N["H"], 120.0), ("L", N["L"], 180.0))]
  return S.DesignStack(S.LayerSpec("air", np.ones(WL.size, complex), 0.0),
                       S.LayerSpec("sub", np.full(WL.shape, 1.52 + 0j), 0.0), films)


def _readback_spec(**angle):
  """Neutral targets: residuals are the simulated curves themselves."""
  (unit, grid), = angle.items()
  spec = S.MeritSpec()
  zeros, ones = np.zeros(WL.size), np.ones(WL.size)
  for a in grid:
    for c in CURVES:
      key = spec.add_key(c, **{unit.replace("angles_", "angle_"): float(a)})
      spec.add_target(key, WL, zeros, ones, "e", "linear", 1.0)
  return spec


def _curves(spec, **angles):
  sim = S.SmatrixContext(spec, WL, **angles).simulate(_stack())
  return np.asarray(spec.residuals(sim)).reshape(-1, len(CURVES), WL.size)


def _targets(**angle):
  (unit, grid), = angle.items()
  tc = TargetCollection()
  for a in grid:
    tc.add(SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, 0.01), "s", "R",
                          **{unit.replace("angles_", "angle_"): float(a)}))
  return tc


# --------------------------------------------------------------------------
# The two units agree
# --------------------------------------------------------------------------

def test_a_radian_grid_reads_each_angle_as_the_engine_does():
  got = _curves(_readback_spec(angles_deg=DEG), angles_rad=RAD)
  sm = ScatterMatrix(IDX, THICK, wavelengths=WL, angles_rad=RAD)
  eng = sm.compute(Request.RS | Request.RP | Request.TS | Request.TP, squeeze=False)
  for ci, c in enumerate(CURVES):
    assert np.max(np.abs(got[:, ci] - eng[c])) < 1e-13, c


def test_the_scatter_matrix_agrees_across_units():
  deg = ScatterMatrix(IDX, THICK, wavelengths=WL, angles_deg=DEG)
  rad = ScatterMatrix(IDX, THICK, wavelengths=WL, angles_rad=RAD)
  assert np.max(np.abs(deg.sin_theta - rad.sin_theta)) < 1e-15
  r_deg = deg.compute(Request.RS, squeeze=False)["Rs"]
  r_rad = rad.compute(Request.RS, squeeze=False)["Rs"]
  assert np.max(np.abs(r_deg - r_rad)) < 1e-14
  assert r_deg[0].tobytes() == r_rad[0].tobytes()


def _structure():
  from navette.structure import Layer, Navette_Structure
  from navette.structure.materials import DictMaterialProvider
  mats = {"glass": np.full(WL.size, 1.52 + 0j), "TiO2": np.full(WL.size, 2.35 + 0j)}
  return Navette_Structure([Layer(0.0, "glass"), Layer(50.0, "TiO2"), Layer(0.0, "glass")],
                           {}, DictMaterialProvider(mats, wavelength=WL))


def test_solve_structure_agrees_across_units():
  deg = solve_structure(_structure(), WL, angles_deg=DEG)
  rad = solve_structure(_structure(), WL, angles_rad=RAD)
  assert np.max(np.abs(np.asarray(deg["Rs"]) - np.asarray(rad["Rs"]))) < 1e-14


def test_a_radian_grid_matches_its_degree_twin():
  spec = _readback_spec(angles_deg=DEG)
  deg = _curves(spec, angles_deg=DEG)
  rad = _curves(spec, angles_rad=RAD)
  assert np.max(np.abs(rad - deg)) < 1e-14
  # Normal incidence is 0 in both units: nothing to round.
  assert rad[0].tobytes() == deg[0].tobytes()


def test_merit_keys_and_sim_rows_agree_across_units():
  spec_deg, spec_rad = _readback_spec(angles_deg=DEG), _readback_spec(angles_rad=RAD)
  rows = {c: np.linspace(0.1, 0.9, DEG.size * WL.size) for c in CURVES}
  sim_deg = sim_curves_from_arrays(WL, rows, angles_deg=DEG)
  sim_rad = sim_curves_from_arrays(WL, rows, angles_rad=RAD)
  ref = np.asarray(spec_deg.residuals(sim_deg))
  for spec in (spec_deg, spec_rad):
    for sim in (sim_deg, sim_rad):
      assert np.asarray(spec.residuals(sim)).tobytes() == ref.tobytes()


def test_targets_keep_the_unit_they_were_given():
  vals, tol = np.zeros_like(WL), np.full_like(WL, 0.01)
  st = SpectralTarget(WL, vals, tol, "s", "R", angle_rad=np.pi / 6)
  assert st.angle_rad == np.pi / 6 and st.angle_deg is None
  assert st._dump()["angle_rad"] == np.pi / 6 and "angle_deg" not in st._dump()
  at = AngularTarget(550.0, np.zeros(3), np.full(3, 0.01), "s", "R", angles_deg=DEG)
  assert at.angles_deg.tobytes() == DEG.tobytes() and at.angles_rad is None
  assert ColorTarget(angle_deg=45.0)._dump()["angle_deg"] == 45.0
  # ColorTarget's angle is optional: left out, it is normal incidence.
  assert ColorTarget()._dump()["angle_rad"] == 0.0


def test_targets_compile_alike_in_either_unit():
  rows = {c: np.linspace(0.1, 0.9, DEG.size * WL.size) for c in CURVES}
  sim = sim_curves_from_arrays(WL, rows, angles_deg=DEG)
  deg = build_merit_spec(_targets(angles_deg=DEG)).residuals(sim)
  rad = build_merit_spec(_targets(angles_rad=RAD)).residuals(sim)
  assert np.max(np.abs(np.asarray(deg) - np.asarray(rad))) == 0.0


def test_the_reference_rotation_agrees_across_units():
  amp = np.exp(1j * np.linspace(0.0, 3.0, WL.size))[None, :]
  n_inc = np.ones(WL.size)
  for d, r in zip(DEG, RAD):
    deg = apply_reference_rotation(amp, WL, n_inc, 400.0, angle_deg=d)
    rad = apply_reference_rotation(amp, WL, n_inc, 400.0, angle_rad=r)
    assert np.max(np.abs(rad - deg)) < 1e-13
  # cos(pi/2) is not 0 in floating point, but it is cos of the radians given.
  rot = np.asarray(S.reference_rotation(WL, n_inc, 400.0, 1.0, angle_rad=np.pi / 2))
  ref = np.exp(-1j * 2 * np.pi * 400.0 * np.cos(np.pi / 2) / WL)
  assert np.max(np.abs(rot - ref)) < 1e-15


def test_the_needle_fold_agrees_across_units():
  spec = build_merit_spec(_targets(angles_deg=DEG))
  deg = S.build_needle_targets(spec, WL, angles_deg=DEG)
  rad = S.build_needle_targets(spec, WL, angles_rad=RAD)
  for k in ("r", "t"):
    for part in ("targets", "weights"):
      assert np.asarray(rad[k][part]).tobytes() == np.asarray(deg[k][part]).tobytes(), k


def _run(**angles):
  return run_needle(layers=SEED, names=["L", "H"], targets=_targets(angles_deg=DEG),
                    wavelengths=WL, materials=MATERIALS, substrate=(1.52, "sub"),
                    pipeline_config=PipelineConfig(max_macro_cycles=2, needles_per_cycle=1),
                    **angles)


def test_run_needle_agrees_across_units():
  base = _run(angles_deg=DEG)
  rad = _run(angles_rad=RAD)
  assert rad["final_layer_count"] == base["final_layer_count"]
  assert rad["final_mf"] == pytest.approx(base["final_mf"], rel=1e-9)
  # At normal incidence the two runs are one run.
  assert _run(angles_rad=[0.0])["final_mf"].hex() == _run(angles_deg=[0.0])["final_mf"].hex()


def _env_kw():
  return dict(targets=_targets(angles_deg=DEG), wavelengths=WL,
              design={"coat": [(N["L"], 300.0, "L"), (N["H"], 200.0, "H")]},
              environments=[{"name": "bare", "stack": [{"design": "coat"}]}],
              pipeline_config=PipelineConfig(max_macro_cycles=1, needles_per_cycle=1),
              substrate=(1.52 + 0j, "sub"))


def test_the_environments_door_agrees_across_units():
  base = run_needle(angles_deg=DEG, **_env_kw())
  rad = run_needle(angles_rad=RAD, **_env_kw())
  assert rad["final_mf"] == pytest.approx(base["final_mf"], rel=1e-9)


def _pipeline(**angles):
  from navette.synthesis.pipeline import stack_from_layers
  stack, pool = stack_from_layers(SEED, WL, MATERIALS, names=["L", "H"])
  spec = build_merit_spec(_targets(angles_deg=DEG))
  cfg = PipelineConfig(max_macro_cycles=1, needles_per_cycle=1)
  return S.NeedlePipeline(stack, spec, WL, pool, cfg, **angles)


def test_the_pipeline_class_agrees_across_units():
  base = _pipeline(angles_deg=DEG).run()
  rad = _pipeline(angles_rad=RAD).run()
  assert rad["final_mf"] == pytest.approx(base["final_mf"], rel=1e-9)


# --------------------------------------------------------------------------
# Exactly one unit, in range, on every door
# --------------------------------------------------------------------------

def _spec():
  return _readback_spec(angles_deg=[0.0])


#: door name -> call taking the angle keywords (grid doors).
GRID_DOORS = {
  "ScatterMatrix": lambda **a: ScatterMatrix(IDX, THICK, wavelengths=WL, **a),
  "Solver": lambda **a: S.Solver(WL, np.ascontiguousarray(
    np.repeat(IDX[:, None], WL.size, axis=1)).ravel(), IDX.size, **a),
  "SimCurves": lambda **a: S.SimCurves(WL, **a),
  "sim_curves_from_arrays": lambda **a: sim_curves_from_arrays(WL, {}, **a),
  "build_needle_targets": lambda **a: S.build_needle_targets(_spec(), WL, **a),
  "SmatrixContext": lambda **a: S.SmatrixContext(_spec(), WL, **a),
  "NeedlePipeline": lambda **a: _pipeline(**a),
  "run_needle": lambda **a: _run(**a),
  "run_needle(environments)": lambda **a: run_needle(**_env_kw(), **a),
  "AngularTarget": lambda **a: AngularTarget(550.0, np.zeros(2), np.full(2, 0.01),
                                             "s", "R", **a),
  "add_angular_target": lambda **a: TargetCollection().build_weaver().add_angular_target(
    550.0, np.zeros(2), np.full(2, 0.01), "s", "R", "e", "auto", **a),
}

#: door name -> call taking the single-angle keywords.
ONE_DOORS = {
  "add_key": lambda **a: S.MeritSpec().add_key("Rs", **a),
  "reference_rotation": lambda **a: S.reference_rotation(WL, np.ones(WL.size), 400.0, 1.0, **a),
  "apply_reference_rotation": lambda **a: apply_reference_rotation(
    np.ones((1, WL.size), complex), WL, np.ones(WL.size), 400.0, **a),
  "SpectralTarget": lambda **a: SpectralTarget(WL, np.zeros_like(WL),
                                               np.full_like(WL, 0.01), "s", "R", **a),
  "ColorTarget": lambda **a: ColorTarget(**a),
  "add_spectral_target": lambda **a: TargetCollection().build_weaver().add_spectral_target(
    WL, np.zeros_like(WL), np.full_like(WL, 0.01), "s", "R", "e", "auto", **a),
}


@pytest.mark.parametrize("door", list(GRID_DOORS))
def test_a_grid_door_takes_exactly_one_unit_in_range(door):
  call = GRID_DOORS[door]
  # Two angles: the two-point grids the target doors need fit everywhere.
  grid = np.array([0.0, 30.0])
  with pytest.raises(ValueError, match="angles_deg.*angles_rad.*not both"):
    call(angles_deg=grid, angles_rad=np.radians(grid))
  with pytest.raises(ValueError, match="angles_deg.*angles_rad"):
    call()
  # Degrees handed over as radians: 30 rad is out of range.
  with pytest.raises(ValueError, match="30 rad"):
    call(angles_rad=grid)
  with pytest.raises(ValueError, match="91 deg"):
    call(angles_deg=np.array([0.0, 91.0]))
  # A NaN is refused too; a target refuses it already as JSON, which has
  # no NaN, so only the array doors can name it.
  with pytest.raises(ValueError, match=None if door == "AngularTarget" else "NaN"):
    call(angles_deg=np.array([0.0, np.nan]))


@pytest.mark.parametrize("door", list(ONE_DOORS))
def test_a_single_angle_door_takes_exactly_one_unit_in_range(door):
  call = ONE_DOORS[door]
  with pytest.raises(ValueError, match="angle_deg.*angle_rad.*not both"):
    call(angle_deg=30.0, angle_rad=np.pi / 6)
  if door != "ColorTarget":  # the one door whose angle is optional
    with pytest.raises(ValueError, match="angle_deg.*angle_rad"):
      call()
  with pytest.raises(ValueError, match="30 rad"):
    call(angle_rad=30.0)
  with pytest.raises(ValueError, match="-1 deg"):
    call(angle_deg=-1.0)
  # The edges are legal: normal and grazing incidence.
  call(angle_deg=90.0)
  call(angle_rad=np.pi / 2)


def test_solve_structure_takes_exactly_one_unit():
  with pytest.raises(ValueError, match="angles_deg.*angles_rad.*not both"):
    solve_structure(_structure(), WL, angles_deg=[0.0], angles_rad=[0.0])
  with pytest.raises(ValueError, match="angles_deg.*angles_rad"):
    solve_structure(_structure(), WL)
  with pytest.raises(ValueError, match="30 rad"):
    solve_structure(_structure(), WL, angles_rad=[0.0, 30.0])
  with pytest.raises(TypeError, match="angles_in_radians"):
    solve_structure(_structure(), WL, angles_deg=[0.0], angles_in_radians=True)


# --------------------------------------------------------------------------
# The old spellings are refused, never read as degrees
# --------------------------------------------------------------------------

@pytest.mark.parametrize("call", [
  lambda: ScatterMatrix(IDX, THICK, wavelengths=WL, angles=[0.0]),
  lambda: ScatterMatrix(IDX, THICK, wavelengths=WL, angles_deg=[0.0],
                        angles_in_radians=False),
  lambda: SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, 0.01), "s", "R",
                         angle=0.0),
  lambda: ColorTarget(angle=0.0),
  lambda: S.MeritSpec().add_key("Rs", angle=0.0),
  lambda: sim_curves_from_arrays(WL, {}, angles=[0.0]),
], ids=["ScatterMatrix(angles=)", "ScatterMatrix(angles_in_radians=)",
        "SpectralTarget(angle=)", "ColorTarget(angle=)", "add_key(angle=)",
        "sim_curves_from_arrays(angles=)"])
def test_the_old_keywords_are_refused(call):
  with pytest.raises(TypeError):
    call()


def test_run_needle_refuses_the_old_flag():
  with pytest.raises(TypeError, match="angles_in_radians"):
    _run(angles_deg=DEG, angles_in_radians=False)


@pytest.mark.parametrize("section,target", [
  ("spectral", {"wavelengths": [500.0], "values": [0.5], "tolerances": [0.1],
                "angle": 30.0, "polarization": "s", "spectral": "R"}),
  ("angular", {"wavelength": 500.0, "values": [0.5], "tolerances": [0.1],
               "angles": [30.0], "polarization": "s", "spectral": "R"}),
])
def test_the_old_json_key_is_refused(section, target):
  stem = "angle" if section == "spectral" else "angles"
  with pytest.raises(ValueError, match=f'"{stem}" is no longer read.*"{stem}_deg" or "{stem}_rad"'):
    S.compile_merit_spec(json.dumps({section: [target]}))


def test_json_keys_name_the_unit():
  base = {"wavelengths": list(WL), "values": [0.1] * WL.size, "tolerances": [0.01] * WL.size,
          "polarization": "s", "spectral": "R"}
  deg = S.compile_merit_spec(json.dumps({"spectral": [{**base, "angle_deg": 30.0}]}))
  rad = S.compile_merit_spec(json.dumps({"spectral": [{**base, "angle_rad": np.pi / 6}]}))
  rows = {"Rs": np.full(WL.size, 0.2)}
  sim = sim_curves_from_arrays(WL, rows, angles_deg=[30.0])
  assert np.asarray(deg.residuals(sim)).tobytes() == np.asarray(rad.residuals(sim)).tobytes()
  with pytest.raises(ValueError, match="not both"):
    S.compile_merit_spec(json.dumps({"spectral": [{**base, "angle_deg": 30.0,
                                                   "angle_rad": 0.5}]}))


# --------------------------------------------------------------------------
# Target-weaver keys hold radians
# --------------------------------------------------------------------------

def _weaver():
  tc = TargetCollection()
  tc.add(SpectralTarget(WL, np.full_like(WL, 0.5), np.full_like(WL, 0.01), "s", "R",
                        angle_deg=30.0))
  return tc.build_weaver()


def test_weaver_keys_hold_radians():
  (entry,) = _weaver().export_entries()
  assert entry["angle_rad"] == np.radians(30.0) and "angle" not in entry


def test_a_simulation_keyed_otherwise_warns():
  sim = OpticalWeaver()
  sim.set_data((np.radians(30.0), "s", "R"), np.full_like(WL, 0.5), WL)
  with warnings.catch_warnings():
    warnings.simplefilter("error")
    assert calculate_merit(sim, _weaver()) == pytest.approx(0.0, abs=1e-12)
  stale = OpticalWeaver()
  stale.set_data((30.0, "s", "R"), np.full_like(WL, 0.5), WL)
  with pytest.warns(UserWarning, match="no simulated curve.*radians"):
    assert calculate_merit(stale, _weaver(), 7.0) == 7.0
