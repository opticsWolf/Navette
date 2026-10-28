# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Case definitions and the evaluator for the forward-model golden baseline.

Shared by `test_forward_golden.py` (recompute and compare),
`test_forward_physics.py` (independent references) and
`gen_forward_golden.py` (the only writer). Everything here is
deterministic: analytic dispersion formulas, fixed grids, literal stacks.
The survey matrix these cases span is `docs/forward-golden.md`.

Two forward paths are pinned, because they are two different solvers:

* the engine -- `Structure.get_solver_inputs()` expanded to solver rows,
  then `ScatterMatrix.compute(mask)` for every request bit (0..50);
* the synthesis evaluator -- `SmatrixContext.simulate(DesignStack)`, the
  dual coherent-block solve every optimisation loop calls, read back
  through `MeritSpec.residuals` with neutral targets (target 0, tolerance
  1, linear or phase transform), one context per angle.
"""
from __future__ import annotations

import warnings
from typing import Any, Dict, List

import numpy as np

FORMAT_VERSION = 1

# Three wavelength ranges, solved separately, each on its own fixed grid.
RANGES = {"uv": (250.0, 400.0), "vis": (400.0, 800.0), "nir": (800.0, 2000.0)}
N_PER_RANGE = 8
# Normal incidence to near grazing.
ANGLES = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 89.0)

N_REQUEST_BITS = 51
ALL_REQUEST = (1 << N_REQUEST_BITS) - 1


def grid(range_name: str) -> np.ndarray:
  lo, hi = RANGES[range_name]
  return np.linspace(lo, hi, N_PER_RANGE)


# ---------------------------------------------------------------------------
# Materials: closed-form dispersion, lambda in nm, Im(n) >= 0 is loss.
# ---------------------------------------------------------------------------
def _drude_ag(wl: np.ndarray) -> np.ndarray:
  # eps = eps_inf - (lambda/lambda_p)^2 / (1 + i*lambda/lambda_gamma)
  x = wl / 20000.0
  eps = 3.7 - (wl / 136.0) ** 2 / (1.0 + 1j * x)
  return np.sqrt(eps)


MATERIALS = {
  "air": lambda wl: np.full(wl.shape, 1.0 + 0.0j),
  "glass": lambda wl: (1.5046 + 4200.0 / wl ** 2) + 0.0j,
  "SiO2": lambda wl: (1.450 + 3600.0 / wl ** 2) + 0.0j,
  "Ta2O5": lambda wl: (2.050 + 25000.0 / wl ** 2) + 0.0j,
  "TiO2": lambda wl: (2.220 + 40000.0 / wl ** 2)
  + 1j * (2e-4 + 0.02 * np.exp(-(wl - 250.0) / 40.0)),
  "Ag": _drude_ag,
  "Si": lambda wl: (3.42 + 1.2e5 / wl ** 2) + 1j * 0.02 * (400.0 / wl) ** 3,
}
LOSSLESS = frozenset({"air", "glass", "SiO2", "Ta2O5"})


def material_table(wl: np.ndarray) -> Dict[str, np.ndarray]:
  return {k: np.ascontiguousarray(f(wl), dtype=np.complex128) for k, f in MATERIALS.items()}


# ---------------------------------------------------------------------------
# Cases. Layer dicts are `navette.structure.Layer` keyword arguments plus
# `material`; the first row is the ambient, the last the substrate.
# ---------------------------------------------------------------------------
def _amb(m: str = "air") -> Dict[str, Any]:
  return {"material": m, "thickness": 0.0, "layer_type": 0}


def _sub(m: str = "glass", **kw: Any) -> Dict[str, Any]:
  return {"material": m, "thickness": 0.0, "layer_type": 2, **kw}


def _f(m: str, d: float, **kw: Any) -> Dict[str, Any]:
  return {"material": m, "thickness": d, **kw}


QW_H, QW_L = 64.47, 94.05  # quarter waves of Ta2O5 / SiO2 at 550 nm


def _qw_mirror(**kw_h: Any) -> List[Dict[str, Any]]:
  films = []
  for i in range(11):
    films.append(_f("Ta2O5", QW_H, **kw_h) if i % 2 == 0 else _f("SiO2", QW_L))
  return films


_EDGE = [41.2, 77.9, 58.3, 101.4, 63.0, 88.8, 49.7, 120.5, 66.1, 93.3,
         55.0, 84.2, 70.6, 97.1, 61.9, 110.0, 52.4, 79.5, 68.8, 104.7, 30.3]

_GRAD_BRUGGEMAN = {"material_a": "Ta2O5", "material_b": "SiO2", "ema": "Bruggeman",
                   "mode": {"FixedSpan": {"f_start": 0.0, "f_end": 1.0}}}
_GRAD_LOOYENGA = {"material_a": "TiO2", "material_b": "SiO2", "ema": "Looyenga",
                  "mode": {"RateCapped": {"f_start": 0.1, "rate": 0.3, "ref_thickness": 100.0,
                                          "f_min": 0.0, "f_max": 0.8}}}
_GRAD_MG = {"material_a": "Ta2O5", "material_b": "SiO2", "ema": "MaxwellGarnett",
            "mode": {"FixedSpan": {"f_start": 0.2, "f_end": 0.7}}, "sublayers": 8}
_RATE_INH = {"RateCapped": {"rate": 0.05, "ref_thickness": 100.0, "cap": 0.3}}

_GROUPS_SCALED = {
  "Ta2O5": {"thick_factor": 1.03, "thick_summand": -0.7, "n_factor": 0.99,
            "k_factor": 1.0, "inh_delta_summand": 0.01, "roughness_summand": 0.4,
            "interface_summand": 0.5},
  "SiO2": {"thick_factor": 0.98, "thick_summand": 0.3, "n_factor": 1.01},
}


def _combined() -> List[Dict[str, Any]]:
  return [
    _amb(),
    _f("TiO2", 85.0, inhomogen=True, inh_delta=0.05, roughness=1.0, rough_type=5),
    _f("SiO2", 110.0, interface=True, interface_thickness=2.0, roughness=1.5, rough_type=4),
    _f("Ta2O5", 140.0, gradient=_GRAD_BRUGGEMAN),
    _f("Ag", 8.0),
    _f("SiO2", 60.0, roughness=0.8, rough_type=1),
    _f("glass", 1.0e6, coherent=False),
    _sub("air"),
  ]


def _rough(rt: int) -> List[Dict[str, Any]]:
  return [_amb(), _f("Ta2O5", 120.0, roughness=2.0, rough_type=rt), _f("SiO2", 90.0), _sub()]


CASES: Dict[str, Dict[str, Any]] = {
  # -- bare substrates (Fresnel) --
  "bare_glass": {"layers": [_amb(), _sub()]},
  "bare_si": {"layers": [_amb(), _sub("Si")]},
  "bare_ag": {"layers": [_amb(), _sub("Ag")]},
  "bare_tir": {"layers": [_amb("glass"), _sub("air")]},
  # -- single films (Airy) --
  "film_sio2": {"layers": [_amb(), _f("SiO2", 250.0), _sub()]},
  "film_tio2_abs": {"layers": [_amb(), _f("TiO2", 80.0), _sub()]},
  "film_ag": {"layers": [_amb(), _f("Ag", 25.0), _sub()]},
  "film_ftir": {"layers": [_amb("glass"), _f("air", 300.0), _sub()]},
  # -- multilayers --
  "qw_mirror": {"layers": [_amb(), *_qw_mirror(), _sub()]},
  "ar_broadband": {"layers": [_amb(), _f("SiO2", 92.0), _f("Ta2O5", 120.0),
                              _f("SiO2", 16.0), _f("Ta2O5", 13.0), _sub()]},
  "mdm_ag": {"layers": [_amb(), _f("SiO2", 80.0), _f("Ag", 15.0), _f("SiO2", 80.0), _sub()]},
  "edge_filter": {"layers": [_amb(), *[_f("Ta2O5" if i % 2 == 0 else "SiO2", d)
                                       for i, d in enumerate(_EDGE)], _sub()]},
  # -- interface roughness, one model each, then mixed --
  "rough_linear": {"layers": _rough(1)},
  "rough_step": {"layers": _rough(2)},
  "rough_exponential": {"layers": _rough(3)},
  "rough_gaussian": {"layers": _rough(4)},
  "rough_nevot_croce": {"layers": _rough(5)},
  "rough_mixed_ml": {"layers": [_amb(), *[
    dict(l, roughness=0.5 + 0.25 * i, rough_type=1 + i % 5) for i, l in enumerate(_qw_mirror())],
    _sub(roughness=1.0, rough_type=4)]},
  "rough_cap_binds": {"layers": [_amb(), _f("Ta2O5", 3.0, roughness=4.0, rough_type=4),
                                 _f("SiO2", 100.0), _sub()], "roughness_cap": 2.0},
  "rough_cap_off": {"layers": [_amb(), _f("Ta2O5", 3.0, roughness=4.0, rough_type=4),
                               _f("SiO2", 100.0), _sub()], "roughness_cap": 0.0},
  # -- inhomogeneity / graded index --
  "inh_fixed": {"layers": [_amb(), _f("TiO2", 150.0, inhomogen=True, inh_delta=0.1), _sub()]},
  "inh_ratecapped": {"layers": [_amb(), _f("Ta2O5", 180.0, inhomogen=True, inh_mode=_RATE_INH),
                                _sub()]},
  "inh_multilayer": {"layers": [_amb(), *_qw_mirror(inhomogen=True, inh_delta=0.04), _sub()]},
  "gradient_bruggeman": {"layers": [_amb(), _f("Ta2O5", 150.0, gradient=_GRAD_BRUGGEMAN), _sub()]},
  "gradient_looyenga_rate": {"layers": [_amb(), _f("TiO2", 200.0, gradient=_GRAD_LOOYENGA),
                                        _sub()]},
  "gradient_maxwell_garnett": {"layers": [_amb(), _f("Ta2O5", 120.0, gradient=_GRAD_MG),
                                          _f("SiO2", 90.0), _sub()]},
  # -- interface slices --
  "interface_slices": {"layers": [_amb(), _f("Ta2O5", 120.0, interface=True, interface_thickness=5.0),
                                  _f("SiO2", 90.0, interface=True, interface_thickness=3.0),
                                  _f("Ta2O5", 120.0, interface=True, interface_thickness=5.0),
                                  _sub()]},
  # -- group scaling (thickness, n/k, and every summand made load-bearing) --
  "group_scaled": {"layers": [_amb(), *[
    dict(l, inhomogen=True, inh_delta=0.02, interface=True, interface_thickness=1.0,
         roughness=0.5, rough_type=4) if l["material"] == "Ta2O5" else l
    for l in _qw_mirror()], _sub()], "groups": _GROUPS_SCALED},
  # -- incoherent / thick layers and the three coherence modes --
  "incoh_substrate_fb": {"layers": [_amb(), _f("Ta2O5", 120.0), _f("SiO2", 90.0),
                                    _f("glass", 1.0e6, coherent=False), _sub("air")],
                         "coherence_mode": 0},
  "incoh_substrate_cm": {"layers": [_amb(), _f("Ta2O5", 120.0), _f("SiO2", 90.0),
                                    _f("glass", 1.0e6, coherent=False), _sub("air")],
                         "coherence_mode": 1},
  "incoh_interior_cm": {"layers": [_amb(), _f("TiO2", 80.0), _f("SiO2", 50000.0, coherent=False),
                                   _f("Ta2O5", 100.0), _sub()], "coherence_mode": 1},
  "fully_coherent_thick": {"layers": [_amb(), _f("Ta2O5", 120.0),
                                      _f("glass", 20000.0, coherent=False), _sub("air")],
                           "coherence_mode": 2},
  # -- everything at once --
  "combined_cm": {"layers": _combined(), "groups": {"SiO2": {"thick_factor": 1.02,
                                                              "roughness_summand": 0.2}},
                  "coherence_mode": 1},
  "combined_fb": {"layers": _combined(), "groups": {"SiO2": {"thick_factor": 1.02,
                                                              "roughness_summand": 0.2}},
                  "coherence_mode": 0},
}


def case_options(case: Dict[str, Any]) -> Dict[str, Any]:
  return {"coherence_mode": int(case.get("coherence_mode", 0)),
          "roughness_cap": float(case.get("roughness_cap", 2.0))}


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------
def build_structure(case: Dict[str, Any], wl: np.ndarray):
  from navette.structure import DictMaterialProvider, Group, Layer, Navette_Structure

  layers = []
  for spec in case["layers"]:
    kw = dict(spec)
    kw["material_name"] = kw.pop("material")
    layers.append(Layer(**kw))
  groups = {name: Group(name, **kw) for name, kw in case.get("groups", {}).items()}
  provider = DictMaterialProvider(material_table(wl), wavelength=wl)
  return Navette_Structure(layers, groups, provider)


def expand(case: Dict[str, Any], wl: np.ndarray) -> Dict[str, np.ndarray]:
  """Solver rows for `case` on `wl` (the structure layer's output)."""
  sa = build_structure(case, wl).get_solver_inputs()
  return {
    "indices": np.ascontiguousarray(sa.indices, dtype=np.complex128),
    "thicknesses": np.ascontiguousarray(sa.thicknesses, dtype=np.float64),
    "flags": np.ascontiguousarray(sa.incoherent_flags, dtype=np.int32),
    "rough_types": np.ascontiguousarray(sa.rough_types, dtype=np.int32),
    "rough_vals": np.ascontiguousarray(sa.rough_vals, dtype=np.float64),
  }


def request_mask(case: Dict[str, Any], rows: Dict[str, np.ndarray]) -> int:
  """Every request bit, less the cross channel where the door refuses it
  (FRONT_BLOCK over an interior incoherent layer, C2)."""
  from navette._smatrix import NEEDS_CROSS

  opts = case_options(case)
  if opts["coherence_mode"] == 0 and np.any(rows["flags"][1:-1] != 0):
    return ALL_REQUEST & ~int(NEEDS_CROSS)
  return ALL_REQUEST


def scatter_matrix(case: Dict[str, Any], rows: Dict[str, np.ndarray], wl: np.ndarray,
                   angles=ANGLES):
  from navette.smatrix.smatrix import ScatterMatrix

  opts = case_options(case)
  with warnings.catch_warnings():
    # The thin-flag advisory (C3) is expected on the 20 um coherent case.
    warnings.simplefilter("ignore", UserWarning)
    return ScatterMatrix(
      rows["indices"], rows["thicknesses"], wavelengths=wl,
      angles_deg=np.asarray(angles, dtype=np.float64),
      incoherent_flags=rows["flags"], roughness_types=rows["rough_types"],
      roughness_values=rows["rough_vals"], coherence_mode=opts["coherence_mode"],
      roughness_cap=opts["roughness_cap"])


def engine_outputs(case: Dict[str, Any], rows: Dict[str, np.ndarray],
                   wl: np.ndarray, angles=ANGLES) -> Dict[str, np.ndarray]:
  sm = scatter_matrix(case, rows, wl, angles)
  out = sm.compute(request_mask(case, rows), squeeze=False)
  return {k: np.asarray(v) for k, v in out.items()}


def synthesis_eligible(case: Dict[str, Any], rows: Dict[str, np.ndarray]) -> bool:
  """The synthesis evaluator solves one coherent block, ambient to substrate."""
  return not np.any(rows["flags"] != 0)


# The baseline pins what the code computes today, known defects included:
# its job is to catch drift. Each known defect is a strict xfail in
# `test_forward_physics.py`; fixing one is an intended change, so the fix
# regenerates the baseline (`gen_forward_golden.py --write`).


SIM_CURVES = ("Rs", "Rp", "Ts", "Tp")
SIM_PHASES = ("Ts", "Tp")


def design_stack(case: Dict[str, Any], rows: Dict[str, np.ndarray]):
  from navette import _smatrix as S

  n, d = rows["indices"], rows["thicknesses"]
  rt, rv = rows["rough_types"], rows["rough_vals"]
  last = n.shape[0] - 1

  def spec(i: int, thick: float):
    return S.LayerSpec(f"row{i}", np.ascontiguousarray(n[i]), thick,
                       rough_type=int(rt[i]), rough_val=float(rv[i]))

  films = [spec(i, float(d[i])) for i in range(1, last)]
  stack = S.DesignStack(spec(0, 0.0), spec(last, 0.0), films)
  stack.roughness_cap = case_options(case)["roughness_cap"]
  return stack


def synthesis_outputs(case: Dict[str, Any], rows: Dict[str, np.ndarray],
                      wl: np.ndarray, angles=ANGLES) -> Dict[str, np.ndarray]:
  """`SmatrixContext.simulate` read back per angle (one context each)."""
  from navette import _smatrix as S

  stack = design_stack(case, rows)
  zeros, ones = np.zeros(wl.size), np.ones(wl.size)
  names = ([f"sim_{c}" for c in SIM_CURVES] + [f"sim_phi_{c[0].lower()}{c[1]}" for c in SIM_PHASES]
           + [f"sim_PD{c[0].lower()}{c[1]}" for c in SIM_PHASES])
  out = {k: np.empty((len(angles), wl.size)) for k in names}
  for ai, a in enumerate(angles):
    spec = S.MeritSpec()
    for c in SIM_CURVES:
      spec.add_target(spec.add_key(c, angle_deg=float(a)), wl, zeros, ones, "e", "linear", 1.0)
    for c in SIM_PHASES:
      spec.add_target(spec.add_key(c, angle_deg=float(a)), wl, zeros, ones, "e", "phase", 1.0, phase=True)
    for c in SIM_PHASES:
      spec.add_target(spec.add_key(c, angle_deg=float(a)), wl, zeros, ones, "e", "phase", 1.0,
                      phase=True, differential_passes=1.0)
    sim = S.SmatrixContext(spec, wl, angles_deg=np.array([float(a)])).simulate(stack)
    res = np.asarray(spec.residuals(sim)).reshape(len(names), wl.size)
    for k, row in zip(names, res):
      out[k][ai] = row
  return out


def evaluate(name: str) -> Dict[str, np.ndarray]:
  """Every golden array of one case: `{range}.{in|out|sim}.{key}`."""
  case = CASES[name]
  arrays: Dict[str, np.ndarray] = {}
  for rng in RANGES:
    wl = grid(rng)
    rows = expand(case, wl)
    for k, v in rows.items():
      arrays[f"{rng}.in.{k}"] = v
    for k, v in engine_outputs(case, rows, wl).items():
      arrays[f"{rng}.out.{k}"] = v
    if synthesis_eligible(case, rows):
      for k, v in synthesis_outputs(case, rows, wl).items():
        arrays[f"{rng}.sim.{k}"] = v
  return arrays


# ---------------------------------------------------------------------------
# Storage: one stacked block per (case, range, path, dtype) instead of one
# zip member per key, which cost more in member headers than the data.
# The key order of each block lives in the manifest's `layout`.
# ---------------------------------------------------------------------------
def pack(name: str, flat: Dict[str, np.ndarray]):
  """`(npz members, layout)` for one case's `evaluate` output."""
  members: Dict[str, np.ndarray] = {}
  layout: Dict[str, List[str]] = {}
  for rng in RANGES:
    groups: Dict[str, List[str]] = {}
    for k in sorted(flat):
      r, path, q = k.split(".", 2)
      if r != rng:
        continue
      v = flat[k]
      if path == "in":
        members[f"{name}.{k}"] = v
        continue
      groups.setdefault(f"{path}.{v.dtype.str}", []).append(q)
    for block, keys in groups.items():
      if layout.setdefault(block, keys) != keys:
        raise ValueError(f"{name}: {block} keys differ between ranges")
      path = block.split(".")[0]
      members[f"{name}.{rng}.{block}"] = np.stack([flat[f"{rng}.{path}.{q}"] for q in keys])
  return members, layout


def unpack(npz, name: str, layout: Dict[str, List[str]]) -> Dict[str, np.ndarray]:
  """Inverse of `pack`: the flat `{range}.{path}.{key}` dict of one case."""
  flat: Dict[str, np.ndarray] = {}
  for rng in RANGES:
    for k in ("indices", "thicknesses", "flags", "rough_types", "rough_vals"):
      flat[f"{rng}.in.{k}"] = npz[f"{name}.{rng}.in.{k}"]
    for block, keys in layout.items():
      stack = npz[f"{name}.{rng}.{block}"]
      path = block.split(".")[0]
      for q, v in zip(keys, stack):
        flat[f"{rng}.{path}.{q}"] = v
  return flat


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------
# Angle-valued keys are compared as wrapped differences, so a value sitting
# on the +-pi branch cut cannot fail by 2*pi.
def is_angle_key(key: str) -> bool:
  q = key.rsplit(".", 1)[-1]
  return (q.startswith(("phi_", "Delta_", "Retardance_", "PD", "sim_phi_", "sim_PD")))


def is_dispersion_key(key: str) -> bool:
  return key.rsplit(".", 1)[-1].startswith(("GD_", "GDD_", "TOD_", "FOD_"))


# Ill-conditioned ratios of small quantities (a phase of a near-zero
# amplitude, a diattenuation of two near-zero reflectances) get the looser
# absolute tolerance; dispersion is relative to the channel's own scale.
_RATIO_PREFIXES = ("Psi_", "DOP_", "Diattenuation_")
TOL_ABS = 1e-12
TOL_RATIO = 1e-9
TOL_DISP_REL = 1e-8


def tolerance(key: str, golden: np.ndarray) -> float:
  q = key.rsplit(".", 1)[-1]
  if is_dispersion_key(key):
    finite = np.abs(golden[np.isfinite(golden)])
    return TOL_DISP_REL * max(1.0, float(finite.max()) if finite.size else 1.0)
  if is_angle_key(key) or q.startswith(_RATIO_PREFIXES):
    return TOL_RATIO
  return TOL_ABS


def deviation(key: str, got: np.ndarray, golden: np.ndarray) -> float:
  """Worst |got - golden| (wrapped for angles); inf on shape/NaN mismatch."""
  if got.shape != golden.shape or got.dtype.kind != golden.dtype.kind:
    return float("inf")
  if got.dtype.kind in "iub":
    return 0.0 if np.array_equal(got, golden) else float("inf")
  nan_g, nan_r = np.isnan(golden), np.isnan(got)
  if not np.array_equal(nan_g, nan_r):
    return float("inf")
  ok = ~nan_g
  if not ok.any():
    return 0.0
  diff = got[ok] - golden[ok]
  if is_angle_key(key):
    diff = np.angle(np.exp(1j * diff))
  return float(np.max(np.abs(diff)))
