# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Gradient mixture films in the synthesis pipeline (F1.1, A5/R4).

Design-path twins: a background gradient film expands WITH the profile
(pinned, silent) and its rows are bit-equal to the same sublayers built
explicitly from the native EMA oracle; a non-background gradient film
homogenizes to ONE row that is bitwise the direct EMA call at f_mid,
announced by a warning naming the mixture and f_mid; the span
bookkeeping survives merge/optimizer passes.
"""
import numpy as np
import pytest

from navette.materials import MaterialSpec, evaluate
from navette.synthesis.pipeline import stack_from_layers

WL = np.linspace(400., 800., 8)
TIO2 = {"model": "Konstant", "params": {"n": 2.35, "k": 0.01}}
SIO2 = {"model": "Konstant", "params": {"n": 1.46}}


def films(stack):
  return stack.to_dict()["films"]


def _manual_sublayers(t_nm, f_start, f_end, host, incl, model="Bruggeman"):
  """The explicit twin rows: n sublayers, inclusive endpoint sampling,
  the last row absorbing the remainder - the expansion branch's exact
  conventions, fed by the native EMA oracle."""
  # The count is read off a first assembly (the count differential is
  # the Rust twin's job; here the VALUES are the contract).
  probe, _ = stack_from_layers(
    [(host, t_nm)], WL, {}, names=["probe"],
    per_film_flags={"probe": {"gradient": {"material_b": incl,
                                            "f_start": f_start,
                                            "f_end": f_end},
                                "optimize": False, "needle": False}})
  n = len(films(probe))
  rows = []
  step = t_nm / n
  for i in range(n):
    d_i = t_nm - step * (n - 1) if i + 1 == n else step
    if i == 0:
      f_i = f_start
    elif i + 1 == n:
      f_i = f_end
    else:
      f_i = f_start + (f_end - f_start) * i / (n - 1)
    nk = evaluate(MaterialSpec(model=model,
                               params={"host": host, "inclusion": incl,
                                       "fraction": f_i}), WL)
    rows.append((nk, d_i))
  return rows


def _merit_ctx():
  from navette._smatrix import MeritSpec, SmatrixContext
  spec = MeritSpec()
  k = spec.add_key(0.0, "Rs")
  spec.add_target(k, WL.copy(), np.zeros(len(WL)),
                  np.full(len(WL), 0.01), "e", "linear", 1.0)
  return SmatrixContext(spec, np.array([0.0]), WL)


def test_gradient_background_rows_match_manual_stack():
  """G-pipeline: a background gradient film's rows are bit-equal to the
  same profile built explicitly from the native EMA oracle - the guards
  preserve physics end to end."""
  import warnings
  grad = {"gradient": {"material_b": TIO2, "f_start": 0.0, "f_end": 1.0},
            "optimize": False, "needle": False}
  with warnings.catch_warnings():
    warnings.simplefilter("error")  # background must NOT warn
    st, _ = stack_from_layers(
      [(SIO2, 120.0)], WL, {}, names=["graded"],
      per_film_flags={"graded": grad})
  subs = [f for f in films(st) if f["material"] == "graded"]
  assert len(subs) > 1
  assert all(not f["optimize"] and not f["needle"] for f in subs)
  assert sum(f["thickness"] for f in subs) == pytest.approx(120.0)

  manual = _manual_sublayers(120.0, 0.0, 1.0, SIO2, TIO2)
  assert len(manual) == len(subs)
  for (nk, d), f in zip(manual, subs):
    assert np.array_equal(np.asarray(f["nk"]), np.asarray(nk)), "row nk bitwise"
    assert f["thickness"] == d
  # Inclusive endpoint sampling: first row IS f_start, last IS f_end.
  assert np.array_equal(
    np.asarray(subs[0]["nk"]),
    np.asarray(evaluate(MaterialSpec(
      model="Bruggeman",
      params={"host": SIO2, "inclusion": TIO2, "fraction": 0.0}), WL)))
  assert np.array_equal(
    np.asarray(subs[-1]["nk"]),
    np.asarray(evaluate(MaterialSpec(
      model="Bruggeman",
      params={"host": SIO2, "inclusion": TIO2, "fraction": 1.0}), WL)))

  # The explicit twin stack simulates to the same merit, bitwise.
  ctx = _merit_ctx()
  st_manual, _ = stack_from_layers(
    [(np.array(nk, dtype=np.complex128), d) for (nk, d) in manual],
    WL, {}, names=[f"m{i}" for i in range(len(manual))],
    film_flags={"optimize": False, "needle": False})
  assert ctx.evaluate_merit(st) == ctx.evaluate_merit(st_manual)


def test_gradient_homogenizes_to_the_ema_at_f_mid():
  """R4: the non-background gradient homogenizes to ONE row, bitwise the
  direct EMA at f_mid, with a warning naming the mixture and f_mid."""
  with pytest.warns(UserWarning, match="carries a gradient"):
    st, _ = stack_from_layers(
      [(TIO2, 100.0)], WL, {}, names=["H"],
      per_film_flags={"H": {"gradient": {"material_b": SIO2,
                                          "f_start": 0.2, "f_end": 0.9,
                                          "ema": "Looyenga"}}})
  fs = films(st)
  assert len(fs) == 1
  assert fs[0]["optimize"] and fs[0]["needle"]
  oracle = evaluate(MaterialSpec(model="Looyenga",
                                 params={"host": TIO2, "inclusion": SIO2,
                                         "fraction": 0.55}), WL)
  assert np.array_equal(np.asarray(fs[0]["nk"]), np.asarray(oracle))
  assert fs[0]["thickness"] == pytest.approx(100.0)


def test_gradient_span_survives_merge_and_optimizer():
  """A background gradient span: nk-keyed merge keeps it (sublayer nk
  differ), the optimizer skips its pinned rows, nothing moves."""
  import warnings
  grad = {"gradient": {"material_b": TIO2, "f_start": 0.0, "f_end": 1.0},
            "optimize": False, "needle": False}
  with warnings.catch_warnings():
    warnings.simplefilter("error")
    st, _ = stack_from_layers(
      [(SIO2, 120.0), (TIO2, 50.0)], WL, {}, names=["graded", "H"],
      per_film_flags={"graded": grad})
  def snapshot(stack):
    return [(f["material"], f["thickness"], f["optimize"], f["needle"],
             np.asarray(f["nk"]).tobytes()) for f in films(stack)]

  before = snapshot(st)
  graded_before = [s for s in before if s[0] == "graded"]
  assert st.merge_adjacent() == 0
  assert [s for s in snapshot(st) if s[0] == "graded"] == graded_before
  ctx = _merit_ctx()
  m0 = ctx.evaluate_merit(st)
  ctx.optimize_thicknesses(st)
  # The gradient span is pinned: its rows are byte-identical after the
  # pass (the optimizable H film is free to move - that is its job).
  assert [s for s in snapshot(st) if s[0] == "graded"] == graded_before


def test_gradient_python_door_refusals():
  """The design door's refusals, each naming the film."""
  with pytest.raises(ValueError, match="material_b"):
    stack_from_layers([(TIO2, 50.0)], WL, {}, names=["H"],
                      per_film_flags={"H": {"gradient": {"f_start": 0.,
                                                          "f_end": 1.}}})
  with pytest.raises(ValueError, match="must be a material"):
    stack_from_layers(
      [(TIO2, 50.0)], WL, {}, names=["H"],
      per_film_flags={"H": {"gradient": {"material_b": "SiO2",
                                          "f_start": 0., "f_end": 1.}}})
  with pytest.raises(ValueError, match="unknown mixing rule"):
    stack_from_layers(
      [(TIO2, 50.0)], WL, {}, names=["H"],
      per_film_flags={"H": {"gradient": {"material_b": SIO2,
                                          "ema": "Wiener",
                                          "f_start": 0., "f_end": 1.}}})
  with pytest.raises(ValueError, match="exactly one key"):
    stack_from_layers(
      [(TIO2, 50.0)], WL, {}, names=["H"],
      per_film_flags={"H": {"gradient": {"material_b": SIO2,
                                          "ema": {"Looyenga": None,
                                                   "Bruggeman": {}},
                                          "f_start": 0., "f_end": 1.}}})
  with pytest.raises(TypeError, match="must be a mapping"):
    stack_from_layers([(TIO2, 50.0)], WL, {}, names=["H"],
                      per_film_flags={"H": {"gradient": 3}})
  # The f-window refusals surface from the native gate.
  with pytest.raises(ValueError, match=r"\[0, 1\]"):
    stack_from_layers(
      [(TIO2, 50.0)], WL, {}, names=["H"],
      per_film_flags={"H": {"gradient": {"material_b": SIO2,
                                          "f_start": -0.2, "f_end": 1.}}})
  with pytest.raises(ValueError, match="f_start == f_end"):
    stack_from_layers(
      [(TIO2, 50.0)], WL, {}, names=["H"],
      per_film_flags={"H": {"gradient": {"material_b": SIO2,
                                          "f_start": 0.5, "f_end": 0.5}}})
  with pytest.raises(ValueError, match="identical"):
    stack_from_layers(
      [(TIO2, 50.0)], WL, {}, names=["TiO2"],
      per_film_flags={"TiO2": {"gradient": {"material_b": TIO2,
                                             "f_start": 0., "f_end": 1.}}})


def test_flat_films_stay_silent_with_gradient_key_present():
  """gradient=None (the _FILM_DEFAULTS default) leaves the path
  untouched: no warning, no extra rows."""
  import warnings
  with warnings.catch_warnings():
    warnings.simplefilter("error")
    st, _ = stack_from_layers([(TIO2, 50.0)], WL, {}, names=["H"],
                              film_flags={"gradient": None})
  assert len(films(st)) == 1
