# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Independent references for the forward-model golden baseline.

`test_forward_golden.py` proves the outputs did not move; this module proves
the values it pins are right where an independent answer exists, so the
baseline is not merely self-consistent. Dense grids, nothing stored:

* bare interfaces against the Fresnel equations (absorbing substrates and
  total internal reflection included), on all three wavelength ranges;
* single films against the analytic Airy summation (complex r and t);
* energy bookkeeping on every golden case: A = 1 - R - T, unpolarised
  averages, R + T = 1 for lossless stacks, R + T <= 1 for graded roughness;
* the ellipsometric and Stokes identities against the complex amplitudes;
* phase continuity on a dense grid and GD against d(phi)/d(omega);
* feature limits that must collapse to the plain stack;
* the two engine doors (`solve_structure`, `ScatterMatrix`) and the
  synthesis evaluator against each other.

Navette's conventions, as measured here: fields carry exp(+i*beta) with
beta = 2*pi*q*d/lambda, q = sqrt(n^2 - (n0*sin(theta))^2) on the branch
Im(q) >= 0; `rp_c` is the negative of the textbook r_p (so rs = rp at
normal incidence); `tp_c` is the textbook t_p times cos(theta_exit) /
cos(theta_0); Psi is in radians; Delta_R = arg(-rp/rs), Delta_T =
arg(tp/ts); S0 = Rs + Rp, S1 = Rp - Rs, S2_R = -2 Re(cross_R), S3_R =
-2 Im(cross_R) with cross_R = rp rs*; GD = +d(phi)/d(omega) in fs.

Two known defects are pinned as strict xfails (and by the baseline, which
pins today's numbers); fixing one flips its xfail and needs a regenerate.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _cases as C  # noqa: E402

from navette.smatrix.smatrix import CoherenceMode, Request, ScatterMatrix  # noqa: E402

DENSE = {rng: np.linspace(lo, hi, 121) for rng, (lo, hi) in C.RANGES.items()}
DENSE_ANGLES = np.linspace(0.0, 89.5, 40)
C_NM_PER_FS = 299.792458


def _q(n, n0, theta_deg):
  """Normal wavevector factor, decaying branch: q = n*cos(theta_n)."""
  kx = n0 * np.sin(np.radians(theta_deg))
  q = np.sqrt(np.asarray(n, dtype=complex) ** 2 - kx ** 2 + 0j)
  return np.where(q.imag < 0, -q, q)


def _fresnel_std(na, nb, qa, qb, pol):
  """Textbook (r, t): s; p in the total-field convention."""
  if pol == "s":
    return (qa - qb) / (qa + qb), 2 * qa / (qa + qb)
  den = nb ** 2 * qa + na ** 2 * qb
  return (nb ** 2 * qa - na ** 2 * qb) / den, 2 * na * nb * qa / den


def _stack_reference(ns, ds, wl, theta):
  """Airy summation over one coherent block -> Navette-convention dict.

  `ns` is a list of per-wavelength index arrays, ambient first; `ds` the
  thicknesses of the interior layers. Recursive from the substrate up.
  """
  n0 = ns[0]
  qs = [_q(n, n0, theta) for n in ns]
  out = {}
  for pol in "sp":
    r, t = _fresnel_std(ns[-2], ns[-1], qs[-2], qs[-1], pol)
    for j in range(len(ns) - 2, 0, -1):
      e1 = np.exp(1j * 2 * np.pi * qs[j] * ds[j - 1] / wl)
      r0, t0 = _fresnel_std(ns[j - 1], ns[j], qs[j - 1], qs[j], pol)
      den = 1 + r0 * r * e1 ** 2
      r, t = (r0 + r * e1 ** 2) / den, t0 * t * e1 / den
    cos0, cos_e = qs[0] / ns[0], qs[-1] / ns[-1]
    if pol == "s":
      out["rs_c"], out["ts_c"] = r, t
      out["Rs"] = np.abs(r) ** 2
      out["Ts"] = np.abs(t) ** 2 * qs[-1].real / qs[0].real
    else:
      out["rp_c"], out["tp_c"] = -r, t * cos_e / cos0
      out["Rp"] = np.abs(r) ** 2
      out["Tp"] = (np.abs(t) ** 2 * (np.conj(ns[-1]) * cos_e).real
                   / (np.conj(ns[0]) * cos0).real)
  return out


REF_KEYS = ("rs_c", "rp_c", "ts_c", "tp_c", "Rs", "Rp", "Ts", "Tp")
REF_REQUEST = (Request.RS | Request.RP | Request.TS | Request.TP | Request.RS_C
               | Request.RP_C | Request.TS_C | Request.TP_C)


def _engine(ns, ds, wl, angles=DENSE_ANGLES, **kw):
  idx = np.array(ns)
  d = [0.0, *ds, 0.0]
  with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)
    sm = ScatterMatrix(idx, d, wavelengths=wl, angles=angles, **kw)
  return sm


def _compare_reference(ns, ds, wl, tol):
  out = _engine(ns, ds, wl).compute(REF_REQUEST, squeeze=False)
  worst = {}
  for ai, a in enumerate(DENSE_ANGLES):
    ref = _stack_reference(ns, ds, wl, a)
    for k in REF_KEYS:
      worst[k] = max(worst.get(k, 0.0), float(np.max(np.abs(out[k][ai] - ref[k]))))
  bad = {k: v for k, v in worst.items() if v > tol}
  assert not bad, f"worst deviation from the analytic reference: {worst}"


M = C.MATERIALS


@pytest.mark.parametrize("rng", list(C.RANGES))
@pytest.mark.parametrize("pair", [("air", "glass"), ("air", "Si"), ("air", "Ag"),
                                  ("glass", "air")])
def test_bare_interface_matches_fresnel(pair, rng):
  wl = DENSE[rng]
  _compare_reference([M[pair[0]](wl), M[pair[1]](wl)], [], wl, 1e-13)


@pytest.mark.parametrize("rng", list(C.RANGES))
def test_total_internal_reflection_is_total(rng):
  wl = DENSE[rng]
  ng, na = M["glass"](wl), M["air"](wl)
  out = _engine([ng, na], [], wl).compute(Request.RS | Request.RP | Request.TS | Request.TP,
                                          squeeze=False)
  crit = np.degrees(np.arcsin(na.real / ng.real))
  beyond = DENSE_ANGLES[:, None] > crit[None, :] + 0.1
  for k in ("Rs", "Rp"):
    assert np.max(np.abs(out[k][beyond] - 1.0)) < 1e-14, k
  for k in ("Ts", "Tp"):
    assert np.max(np.abs(out[k][beyond])) < 1e-14, k


@pytest.mark.parametrize("rng", list(C.RANGES))
@pytest.mark.parametrize("film", [("air", "SiO2", 250.0, "glass"), ("air", "TiO2", 80.0, "glass"),
                                  ("air", "Ag", 25.0, "glass"), ("glass", "air", 300.0, "glass"),
                                  ("air", "Ta2O5", 120.0, "Si")])
def test_single_film_matches_airy(film, rng):
  wl = DENSE[rng]
  amb, mat, d, sub = film
  _compare_reference([M[amb](wl), M[mat](wl), M[sub](wl)], [d], wl, 1e-12)


@pytest.mark.parametrize("rng", list(C.RANGES))
def test_quarter_wave_mirror_matches_recursive_airy(rng):
  wl = DENSE[rng]
  layers = C.CASES["qw_mirror"]["layers"]
  ns = [M[l["material"]](wl) for l in layers]
  ds = [l["thickness"] for l in layers[1:-1]]
  _compare_reference(ns, ds, wl, 1e-11)


# ---------------------------------------------------------------------------
# Energy bookkeeping, on the golden cases themselves
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def evaluated():
  return {name: C.evaluate(name) for name in C.CASES}


def _per_range(flat, rng):
  return {k.split(".", 2)[2]: v for k, v in flat.items() if k.startswith(f"{rng}.out.")}


def _rows(flat, rng):
  return {k.split(".", 2)[2]: v for k, v in flat.items() if k.startswith(f"{rng}.in.")}


@pytest.mark.parametrize("name", list(C.CASES))
def test_energy_bookkeeping(name, evaluated):
  for rng in C.RANGES:
    o, rows = _per_range(evaluated[name], rng), _rows(evaluated[name], rng)
    for p in ("s", "p"):
      R, T, A = o[f"R{p}"], o[f"T{p}"], o[f"A_{p}"]
      assert np.max(np.abs(A - (1.0 - R - T))) < 1e-15, (rng, p)
    assert np.max(np.abs(o["R_avg"] - 0.5 * (o["Rs"] + o["Rp"]))) < 1e-15
    assert np.max(np.abs(o["T_avg"] - 0.5 * (o["Ts"] + o["Tp"]))) < 1e-15
    assert np.max(np.abs(o["A_avg"] - 0.5 * (o["A_s"] + o["A_p"]))) < 1e-15
    lossless = not np.any(rows["indices"].imag != 0.0)
    rough = rows["rough_types"][rows["rough_vals"] > 0.0]
    if lossless and rough.size == 0:
      for p in ("s", "p"):
        dev = np.max(np.abs(o[f"R{p}"] + o[f"T{p}"] - 1.0))
        assert dev < 1e-12, f"{rng} {p}: lossless stack, |R+T-1| = {dev:.2e}"
    if not np.any(rough == 5):
      # Graded roughness and absorption only ever remove specular energy;
      # Nevot-Croce may add some (documented in RoughnessType).
      for p in ("s", "p"):
        assert np.min(o[f"A_{p}"]) > -1e-12, (rng, p)
        assert np.max(o[f"R{p}"]) <= 1.0 + 1e-12 and np.min(o[f"T{p}"]) >= -1e-15


# ---------------------------------------------------------------------------
# Ellipsometry and Stokes against the complex amplitudes (single block)
# ---------------------------------------------------------------------------
def _wrapped(x):
  return np.abs(np.angle(np.exp(1j * x)))


@pytest.mark.parametrize("name", [n for n in C.CASES if C.CASES[n].get("coherence_mode", 0) != 2])
@pytest.mark.filterwarnings("ignore:invalid value encountered:RuntimeWarning")
def test_ellipsometric_and_stokes_identities(name, evaluated):
  for rng in C.RANGES:
    o, rows = _per_range(evaluated[name], rng), _rows(evaluated[name], rng)
    if np.any(rows["flags"] != 0):
      if "Psi_R" in o:  # amplitude ratio of totals: still defined
        ok = o["Rs"] > 1e-10
        assert np.max(np.abs(o["Psi_R"] - np.arctan(np.sqrt(o["Rp"] / o["Rs"])))[ok]) < 1e-12
      if "DOP_R" in o:  # partially coherent echoes depolarise, never over-polarise
        assert np.max(o["DOP_R"]) <= 1.0 + 1e-12 and np.max(o["DOP_T"]) <= 1.0 + 1e-12
      continue
    rs, rp, ts, tp = o["rs_c"], o["rp_c"], o["ts_c"], o["tp_c"]
    ok_r = o["Rs"] > 1e-10
    assert np.max(np.abs(o["Psi_R"] - np.arctan(np.sqrt(o["Rp"] / o["Rs"])))[ok_r]) < 1e-12
    assert np.max(_wrapped(o["Delta_R"] - np.angle(-rp / rs))[ok_r]) < 1e-9
    ok_t = o["Ts"] > 1e-10
    assert np.max(np.abs(o["Psi_T"] - np.arctan(np.sqrt(o["Tp"] / o["Ts"])))[ok_t]) < 1e-12
    assert np.max(_wrapped(o["Delta_T"] - np.angle(tp / ts))[ok_t]) < 1e-9
    assert np.max(np.abs(o["cross_R"] - rp * np.conj(rs))) < 1e-15
    assert np.max(np.abs(o["cross_T"] - tp * np.conj(ts))) < 1e-14
    assert np.max(np.abs(o["S0_R"] - (o["Rs"] + o["Rp"]))) < 1e-15
    assert np.max(np.abs(o["S1_R"] - (o["Rp"] - o["Rs"]))) < 1e-15
    assert np.max(np.abs(o["S2_R"] + 2 * o["cross_R"].real)) < 1e-15
    assert np.max(np.abs(o["S3_R"] + 2 * o["cross_R"].imag)) < 1e-15
    assert np.max(np.abs(o["S0_T"] - (o["Ts"] + o["Tp"]))) < 1e-15
    assert np.max(np.abs(o["S1_T"] - (o["Tp"] - o["Ts"]))) < 1e-15
    pure = o["S1_R"] ** 2 + o["S2_R"] ** 2 + o["S3_R"] ** 2
    assert np.max(np.abs(pure - o["S0_R"] ** 2)) < 1e-14, "coherent reflection is a pure state"
    assert np.max(np.abs(o["DOP_R"] - 1.0)[o["S0_R"] > 1e-8]) < 1e-9
    d_r = (o["Rp"] - o["Rs"]) / (o["Rp"] + o["Rs"])
    assert np.max(np.abs(o["Diattenuation_R"] - d_r)[o["S0_R"] > 1e-8]) < 1e-12
    for p in ("rs", "rp", "ts", "tp", "rbs", "rbp", "tbs", "tbp"):
      assert np.max(_wrapped(o[f"phi_{p}"] - np.angle(o[f"{p}_c"]))) < 1e-15, p
    # Below the Psi/Delta floor (T = 0 past TIR) Delta takes its degenerate
    # value and the retardance does not; compare where there is signal.
    assert np.max(_wrapped(o["Retardance_R"] - o["Delta_R"] + np.pi)[ok_r]) < 1e-12
    assert np.max(_wrapped(o["Retardance_T"] - o["Delta_T"])[ok_t]) < 1e-12


@pytest.mark.parametrize("name", ["film_sio2", "qw_mirror", "combined_cm"])
def test_differential_phase_subtracts_the_equivalent_medium(name, evaluated):
  case = C.CASES[name]
  for rng in C.RANGES:
    o, rows = _per_range(evaluated[name], rng), _rows(evaluated[name], rng)
    wl = C.grid(rng)
    D = float(np.sum(rows["thicknesses"][1:-1]))
    n0 = rows["indices"][0].real
    cos0 = np.cos(np.radians(np.asarray(C.ANGLES)))[:, None]
    rot = np.exp(-1j * 2 * np.pi * n0[None, :] * D * cos0 / wl[None, :])
    for p in ("s", "p"):
      ref = np.angle(o[f"t{p}_c"] * rot)
      assert np.max(_wrapped(o[f"PDt{p}"] - ref)) < 1e-9, (name, rng, p, case.get("coherence_mode"))


# ---------------------------------------------------------------------------
# Phase continuity and group delay
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("stack", ["film_sio2", "qw_mirror", "mdm_ag"])
def test_phase_is_continuous_and_gd_is_its_derivative(stack):
  wl = np.linspace(400.0, 800.0, 4001)
  layers = C.CASES[stack]["layers"]
  ns = [M[l["material"]](wl) for l in layers]
  ds = [l["thickness"] for l in layers[1:-1]]
  sm = _engine(ns, ds, wl, angles=np.array([0.0, 45.0]))
  out = sm.compute(Request.PHI_RS | Request.PHI_TP | Request.DISP_R_S | Request.DISP_T_P,
                   squeeze=False)
  omega = 2 * np.pi * C_NM_PER_FS / wl
  for phi, gd in (("phi_rs", "GD_R_s"), ("phi_tp", "GD_T_p")):
    for ai in range(2):
      up = np.unwrap(out[phi][ai])
      step = np.max(np.abs(np.diff(up)))
      assert step < 0.2, f"{stack} {phi}: jump {step:.3f} rad between neighbours"
      ref = np.gradient(up, omega)
      scale = np.max(np.abs(ref))
      inner = slice(5, -5)
      dev = np.max(np.abs(out[gd][ai][inner] - ref[inner])) / scale
      assert dev < 1e-4, f"{stack} {gd}: |GD - dphi/domega| / max = {dev:.2e}"


# ---------------------------------------------------------------------------
# Feature limits
# ---------------------------------------------------------------------------
def _solve_layers(layers, wl, groups=None, **opts):
  case = {"layers": layers, "groups": groups or {}, **opts}
  rows = C.expand(case, wl)
  return C.engine_outputs(case, rows, wl)


@pytest.mark.filterwarnings("ignore:.*inh_delta 0 expands:UserWarning")
@pytest.mark.parametrize("variant", [
  {"inhomogen": True, "inh_delta": 0.0},
  {"roughness": 0.0, "rough_type": 4},
  {"roughness": 0.0, "rough_type": 5},
])
def test_inert_feature_settings_collapse_to_the_plain_film(variant):
  wl = C.grid("vis")
  base = [{"material": "air", "thickness": 0.0, "layer_type": 0},
          {"material": "TiO2", "thickness": 150.0},
          {"material": "glass", "thickness": 0.0, "layer_type": 2}]
  feat = [dict(l) for l in base]
  feat[1].update(variant)
  a, b = _solve_layers(base, wl), _solve_layers(feat, wl)
  for k in ("Rs", "Rp", "Ts", "Tp", "rs_c", "tp_c"):
    assert np.max(np.abs(a[k] - b[k])) < 1e-12, (variant, k)


def test_identity_group_is_inert():
  wl = C.grid("vis")
  layers = C.CASES["qw_mirror"]["layers"]
  a = _solve_layers(layers, wl)
  b = _solve_layers(layers, wl, groups={"Ta2O5": {}, "SiO2": {}})
  for k in a:
    assert np.array_equal(a[k], b[k], equal_nan=True), k


# ---------------------------------------------------------------------------
# Door and path agreement
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", list(C.CASES))
def test_solve_structure_door_matches_scatter_matrix(name):
  from navette.structure import solve_structure

  case = C.CASES[name]
  opts = C.case_options(case)
  if opts["roughness_cap"] != 2.0:
    pytest.skip("solve_structure has no roughness_cap option")
  for rng in C.RANGES:
    wl = C.grid(rng)
    rows = C.expand(case, wl)
    mask = C.request_mask(case, rows)
    want = C.engine_outputs(case, rows, wl)
    with warnings.catch_warnings():
      warnings.simplefilter("ignore", UserWarning)
      got = solve_structure(C.build_structure(case, wl), wl, list(C.ANGLES), request=mask,
                            coherence_mode=opts["coherence_mode"])
    assert sorted(got) == sorted(want)
    for k in want:
      assert np.array_equal(np.asarray(got[k]), want[k], equal_nan=True), (rng, k)


@pytest.mark.parametrize("name", [n for n in C.CASES
                                  if not any(not l.get("coherent", True)
                                             for l in C.CASES[n]["layers"])])
def test_synthesis_evaluator_matches_the_engine(name, evaluated):
  """Two solvers, one answer: the dual block solve vs the engine (R, phases)."""
  for rng in C.RANGES:
    flat = evaluated[name]
    absorbing = np.any(flat[f"{rng}.in.indices"][-1].imag > 0.0)
    for c in C.SIM_CURVES:
      if absorbing and c.startswith("T"):
        continue  # known defect, see the xfail below
      dev = np.max(np.abs(flat[f"{rng}.sim.sim_{c}"] - flat[f"{rng}.out.{c}"]))
      assert dev < 1e-13, (rng, c, dev)
    for p in ("ts", "tp"):
      for sim, eng in ((f"sim_phi_{p}", f"phi_{p}"), (f"sim_PD{p}", f"PD{p}")):
        dev = np.max(_wrapped(flat[f"{rng}.sim.{sim}"] - flat[f"{rng}.out.{eng}"]))
        assert dev < 1e-12, (rng, sim, dev)


# ---------------------------------------------------------------------------
# Known defects (found while building this baseline). Strict: a fix flips
# them to XPASS-fail, which is the cue to drop the mark and regenerate.
# ---------------------------------------------------------------------------
def test_multi_angle_synthesis_reads_each_angle():
  from navette import _smatrix as S

  wl = C.grid("vis")
  case = C.CASES["qw_mirror"]
  rows = C.expand(case, wl)
  stack = C.design_stack(case, rows)
  angles = np.array([0.0, 30.0, 60.0])
  spec = S.MeritSpec()
  for a in angles:
    spec.add_target(spec.add_key(float(a), "Rs"), wl, np.zeros(wl.size), np.ones(wl.size),
                    "e", "linear", 1.0)
  got = np.asarray(spec.residuals(S.SmatrixContext(spec, angles, wl).simulate(stack)))
  want = C.engine_outputs(case, rows, wl, angles)["Rs"]
  assert np.max(np.abs(got.reshape(want.shape) - want)) < 1e-13


@pytest.mark.xfail(strict=True, reason=(
  "SmatrixContext.simulate's Ts/Tp into an absorbing substrate disagree with the "
  "Poynting transmittance the engine returns (T = 84 on bare Ag at 400 nm)."))
@pytest.mark.parametrize("sub", ["Ag", "Si"])
def test_synthesis_transmittance_into_absorbing_substrate(sub):
  wl = C.grid("vis")
  case = {"layers": [{"material": "air", "thickness": 0.0, "layer_type": 0},
                     {"material": sub, "thickness": 0.0, "layer_type": 2}]}
  rows = C.expand(case, wl)
  sim = C.synthesis_outputs(case, rows, wl)
  eng = C.engine_outputs(case, rows, wl)
  assert np.max(np.abs(sim["sim_Ts"] - eng["Ts"])) < 1e-13
  assert np.max(np.abs(sim["sim_Tp"] - eng["Tp"])) < 1e-13


@pytest.mark.parametrize("n_exit", [
  1.0,
  pytest.param(1.52, marks=pytest.mark.xfail(strict=True, reason=(
    "S2_T/S3_T are built from tp*ts* without the exit-medium factor that S0_T/S1_T "
    "carry, so DOP_T of a coherent stack is < 1 whenever n_exit != n_ambient "
    "(0.658 for a film on glass at normal incidence)."))),
])
def test_coherent_transmission_is_a_pure_state(n_exit):
  wl = DENSE["vis"]
  sm = _engine([np.full(wl.size, 1.0 + 0j), np.full(wl.size, 2.0 + 0j),
                np.full(wl.size, n_exit + 0j)], [120.0], wl, angles=np.array([0.0, 30.0, 60.0]))
  o = sm.compute(Request.DOP_T, squeeze=False)
  assert np.max(np.abs(o["DOP_T"] - 1.0)) < 1e-12
