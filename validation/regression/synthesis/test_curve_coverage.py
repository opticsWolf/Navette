# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Every target type the compiler accepts is scored by synthesis (W2, W3).

The synthesis evaluator used to fill only Rs, Rp, Ts, Tp and the complex
forward t. A demand on anything else (unpolarized, back incidence, front
or back reflection phase) cost a constant 1e6 in the merit and stopped
the LM with "missing curve". The evaluator now fills every row a spec
reads, and these tests pin it:

  * every label x polarization x phase flag constructs, scores a finite
    merit (not the penalty), runs the thickness LM and folds for the
    needle pass;
  * each new row equals the forward engine (``ScatterMatrix``): targets
    set to the engine's own values score a merit of ~0;
  * back incidence from an absorbing exit medium is refused, naming it;
  * a needle run on an unpolarized target completes.

Unpolarized phase and colour on back curves stay refused at compile;
``MeritSpec::check_supported`` (Rust) refuses anything an evaluator does
not produce before a run starts.
"""
import numpy as np
import pytest

from navette import _smatrix as S
from navette.smatrix.smatrix import Request, ScatterMatrix
from navette.spectralweave.target import AngularTarget, ColorTarget, SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec
from navette.synthesis.pipeline import PipelineConfig, run_needle

WL = np.array([500.0, 560.0, 620.0, 680.0])
DEG = 45.0
ANG = np.array([DEG])
N = {"air": 1.0 + 0j, "H": 2.35 + 0.01j, "L": 1.46 + 0j, "sub": 1.52 + 0j}
D = {"H": 80.0, "L": 100.0}


def _stack(sub=N["sub"]):
  films = [S.LayerSpec(k, np.full(WL.shape, N[k]), D[k]) for k in ("H", "L")]
  return S.DesignStack(S.LayerSpec("air", np.full(WL.shape, N["air"]), 0.0, optimize=False, needle=False),
                       S.LayerSpec("sub", np.full(WL.shape, sub), 0.0, optimize=False, needle=False),
                       films)


def _spec(*targets):
  tc = TargetCollection()
  for t in targets:
    tc.add(t)
  return build_merit_spec(tc)


def _ctx(*targets):
  return S.SmatrixContext(_spec(*targets), WL, angles_deg=ANG)


V, T = np.full(WL.size, 0.1), np.full(WL.size, 0.01)
CASES = {}
for sp in ("R", "T", "A", "RB", "TB", "AB"):
  for pol in ("s", "p", "u"):
    CASES[f"{sp}{pol}"] = (lambda sp=sp, pol=pol: SpectralTarget(WL, V, T, pol, sp, angle_deg=DEG))
for sp in ("R", "T", "RB", "TB"):
  for pol in ("s", "p"):
    CASES[f"{sp}{pol}-phase"] = (lambda sp=sp, pol=pol: SpectralTarget(
      WL, np.full(WL.size, 0.5), np.full(WL.size, 0.1), pol, sp, angle_deg=DEG, phase=True))
CASES["angular-Ru"] = lambda: AngularTarget(560.0, V[:1], T[:1], "u", "R", angles_deg=ANG)
CASES["colour-Ru"] = lambda: ColorTarget(curve="Ru", angle_deg=DEG, illuminant="D65",
                                         observer="1931_2deg", quantity="Lab",
                                         reference=[50.0, 0.0, 0.0], distance="DeltaE76")


@pytest.mark.parametrize("name", list(CASES))
def test_every_target_type_scores_optimizes_and_folds(name):
  spec = _spec(CASES[name]())
  ctx = S.SmatrixContext(spec, WL, angles_deg=ANG)
  stack = _stack()
  m = ctx.evaluate_merit(stack)
  assert np.isfinite(m) and m != 1e6, m
  m_opt = ctx.optimize_thicknesses(stack)
  assert np.isfinite(m_opt) and m_opt <= m * (1 + 1e-12), (m, m_opt)
  S.build_needle_targets(spec, WL, ctx.simulate(stack), angles_deg=ANG)


def _engine():
  n = np.array([np.full(WL.shape, N[k]) for k in ("air", "H", "L", "sub")])
  st = ScatterMatrix(n, np.array([0.0, D["H"], D["L"], 0.0]), wavelengths=WL, angles_deg=[DEG])
  req = (Request.R_AVG | Request.T_AVG | Request.PHI_RS | Request.PHI_RP
         | Request.PHI_RBS | Request.PHI_RBP | Request.PHI_TBS | Request.PHI_TBP
         | Request.RBS_C | Request.RBP_C | Request.TBS_C | Request.TBP_C | Request.TS | Request.TP)
  return st.compute(req)


def _back_flux_factor(pol):
  """Re(y_ambient)/Re(y_exit) for back incidence, as the solver's finalize."""
  sin0 = np.sin(np.radians(DEG))
  cos = lambda n: np.sqrt(1 - (sin0 / n) ** 2 + 0j)
  y = (lambda n: n * cos(n)) if pol == "s" else (lambda n: n / cos(n))
  return y(N["air"]).real / y(N["sub"]).real


def _oracle_rows():
  e = _engine()
  rows = {
    ("R", "u", False): e["R_avg"], ("T", "u", False): e["T_avg"],
    ("R", "s", True): e["phi_rs"], ("R", "p", True): e["phi_rp"],
    ("RB", "s", True): e["phi_rbs"], ("RB", "p", True): e["phi_rbp"],
    ("TB", "s", True): e["phi_tbs"], ("TB", "p", True): e["phi_tbp"],
  }
  for pol in ("s", "p"):
    rb, tb = np.abs(e[f"rb{pol}_c"]) ** 2, np.abs(e[f"tb{pol}_c"]) ** 2 * _back_flux_factor(pol)
    rows[("RB", pol, False)] = rb
    rows[("TB", pol, False)] = tb
    # Reciprocity: back transmission equals front transmission.
    assert np.max(np.abs(tb - e[f"T{pol}"])) < 1e-12
  rows[("RB", "u", False)] = 0.5 * (rows[("RB", "s", False)] + rows[("RB", "p", False)])
  rows[("TB", "u", False)] = 0.5 * (rows[("TB", "s", False)] + rows[("TB", "p", False)])
  return rows


@pytest.mark.parametrize("key", list(_oracle_rows()), ids=lambda k: f"{k[0]}{k[1]}{'-phase' if k[2] else ''}")
def test_new_rows_equal_the_forward_engine(key):
  sp, pol, phase = key
  values = np.asarray(_oracle_rows()[key], dtype=float)
  tol = np.full(WL.size, 1e-3)
  ctx = _ctx(SpectralTarget(WL, values, tol, pol, sp, angle_deg=DEG, phase=phase))
  # Residual (sim - engine)/1e-3 per point: merit < 1e-14 is agreement to ~1e-10.
  assert ctx.evaluate_merit(_stack()) < 1e-14


def test_back_demands_refuse_an_absorbing_exit_medium():
  ctx = _ctx(SpectralTarget(WL, V, T, "s", "TB", angle_deg=DEG))
  with pytest.raises(ValueError, match="exit medium absorbs at 500 nm"):
    ctx.evaluate_merit(_stack(sub=0.05 + 3.0j))
  # A front demand on the same metal substrate is fine.
  _ctx(SpectralTarget(WL, V, T, "s", "R", angle_deg=DEG)).evaluate_merit(_stack(sub=0.05 + 3.0j))


def test_a_run_drops_back_demands_on_an_absorbing_exit_once():
  # The run door checks the complete stack once, at setup: the back demand
  # goes with ONE warning, the front demand still runs, and no solve inside
  # the pipeline reaches the evaluator's per-call refusal.
  import warnings
  tc = TargetCollection()
  tc.add(SpectralTarget(WL, np.zeros(WL.size), np.full(WL.size, 0.01), "s", "R", angle_deg=DEG))
  tc.add(SpectralTarget(WL, V, T, "s", "TB", angle_deg=DEG))
  mats = {k: np.full(WL.shape, N[k]) for k in ("H", "L")}
  with warnings.catch_warnings(record=True) as seen:
    warnings.simplefilter("always")
    out = run_needle(layers=[(mats["L"], 100.0), (mats["H"], 80.0)], names=["L", "H"], targets=tc,
                     wavelengths=WL, materials=mats, substrate=(np.full(WL.shape, 0.05 + 3.0j), "metal"),
                     angles_deg=ANG,
                     pipeline_config=PipelineConfig(max_macro_cycles=2, needles_per_cycle=1))
  dropped = [str(w.message) for w in seen if "back-incidence" in str(w.message)]
  assert len(dropped) == 1 and "'metal'" in dropped[0], [str(w.message) for w in seen]
  assert np.isfinite(out["final_mf"])


def test_a_needle_run_on_an_unpolarized_target_completes():
  tc = TargetCollection()
  tc.add(SpectralTarget(WL, np.zeros(WL.size), np.full(WL.size, 0.01), "u", "R", angle_deg=DEG))
  mats = {k: np.full(WL.shape, N[k]) for k in ("H", "L")}
  out = run_needle(layers=[(mats["L"], 100.0), (mats["H"], 80.0)], names=["L", "H"], targets=tc,
                   wavelengths=WL, materials=mats, substrate=(1.52, "sub"), angles_deg=ANG,
                   pipeline_config=PipelineConfig(max_macro_cycles=2, needles_per_cycle=1))
  assert np.isfinite(out["final_mf"]) and out["final_mf"] < 1e6
