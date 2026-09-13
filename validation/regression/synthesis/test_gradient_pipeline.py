# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Gradient mixture films in the synthesis pipeline (F1.1/F1.2, A5/R4).

Design-path twins: a background gradient film expands WITH the profile
(pinned, silent) and its rows are bit-equal to the same sublayers built
explicitly from the native EMA oracle; a non-background gradient film
homogenizes to ONE row that is bitwise the direct EMA call at f_mid,
announced by a warning naming the mixture and f_mid; the span
bookkeeping survives merge/optimizer passes; RateCapped saturates into
a bitwise pure-material tail (F1.2).
"""
import itertools

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


def test_gradient_optimize_keeps_the_profile_and_homogenize_path_still_warns():
  """F1.6 licence item 1/2 (the one non-additive change): a FixedSpan
  gradient with optimize=true no longer homogenizes - it keeps the
  profile as ONE scalable span, silently. The R4 homogenize (warning,
  EMA at f_mid) survives for the posture that always homogenized:
  optimize=false + needle=true (not background)."""
  import warnings
  # The F1.6 posture: full profile, no warning.
  with warnings.catch_warnings():
    warnings.simplefilter("error")  # must NOT warn
    st, _ = stack_from_layers(
      [(TIO2, 100.0)], WL, {}, names=["H"],
      per_film_flags={"H": {"gradient": {"material_b": SIO2,
                                          "f_start": 0.2, "f_end": 0.9,
                                          "ema": "Looyenga"}}})
  fs = films(st)
  assert len(fs) > 1, "the profile kept its sublayers"
  assert all(f["optimize"] for f in fs)

  # The stable homogenize posture still warns, and the row is still
  # bitwise the direct EMA at f_mid.
  with pytest.warns(UserWarning, match="carries a gradient"):
    st2, _ = stack_from_layers(
      [(TIO2, 100.0)], WL, {}, names=["H"],
      per_film_flags={"H": {"gradient": {"material_b": SIO2,
                                          "f_start": 0.2, "f_end": 0.9,
                                          "ema": "Looyenga"},
                             "optimize": False, "needle": True}})
  fs2 = films(st2)
  assert len(fs2) == 1
  assert not fs2[0]["optimize"] and fs2[0]["needle"]
  oracle = evaluate(MaterialSpec(model="Looyenga",
                                 params={"host": TIO2, "inclusion": SIO2,
                                         "fraction": 0.55}), WL)
  assert np.array_equal(np.asarray(fs2[0]["nk"]), np.asarray(oracle))
  assert fs2[0]["thickness"] == pytest.approx(100.0)


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

  graded_before = [s for s in snapshot(st) if s[0] == "graded"]
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
                                          "f_end": 1.}}})
  with pytest.raises(ValueError, match="exactly one key"):
    stack_from_layers(
      [(TIO2, 50.0)], WL, {}, names=["H"],
      per_film_flags={"H": {"gradient": {"material_b": SIO2,
                                          "ema": {"Looyenga": None,
                                                   "Bruggeman": {}},
                                          "f_end": 1.}}})
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
                                             "f_end": 1.}}})


def test_flat_films_stay_silent_with_gradient_key_present():
  """gradient=None (the _FILM_DEFAULTS default) leaves the path
  untouched: no warning, no extra rows."""
  import warnings
  with warnings.catch_warnings():
    warnings.simplefilter("error")
    st, _ = stack_from_layers([(TIO2, 50.0)], WL, {}, names=["H"],
                              film_flags={"gradient": None})
  assert len(films(st)) == 1


# ---------------------------------------------------------------------------
# F1.2 - RateCapped: slope with caps, saturation, the N2 merge twin.
# ---------------------------------------------------------------------------


def _oracle(host, incl, f, model="Bruggeman"):
  return evaluate(MaterialSpec(
    model=model, params={"host": host, "inclusion": incl, "fraction": f}), WL)


def test_rate_capped_saturates_to_bitwise_pure_tail():
  """F1.2 (D6 G-thickness, type b): a rate-driven film saturates where
  the cap binds - the tail rows are BITWISE pure material_b (EMA at
  exactly f_max), the head rows follow the formula, and a 200 nm film
  of the same spec never reaches the cap."""
  import warnings
  rate, f_start = 0.25, 0.1
  grad = {"gradient": {"material_b": TIO2, "rate": rate, "f_start": f_start},
            "optimize": False, "needle": False}
  with warnings.catch_warnings():
    warnings.simplefilter("error")
    st, _ = stack_from_layers(
      [(SIO2, 400.0)], WL, {}, names=["g"], per_film_flags={"g": grad})
  rows = [f for f in films(st) if f["material"] == "g"]
  n = len(rows)
  assert sum(f["thickness"] for f in rows) == pytest.approx(400.0)
  tail_seen = False
  for i, f in enumerate(rows):
    z = 400.0 * i / (n - 1)
    raw = f_start + rate * (z / 100.0)
    f_i = min(max(raw, 0.0), 1.0)
    assert np.array_equal(np.asarray(f["nk"]), np.asarray(_oracle(SIO2, TIO2, f_i))), \
      f"row {i} (f={f_i})"
    tail_seen = tail_seen or raw >= 1.0
  assert tail_seen, "the 400 nm film must saturate"
  # The cap rows are bitwise pure material_b - the saturation proof.
  pure_b = _oracle(SIO2, TIO2, 1.0)
  assert any(np.array_equal(np.asarray(f["nk"]), np.asarray(pure_b))
             for f in rows)
  # 200 nm of the same spec (rate 0.3): the far face sits at
  # 0.1 + 0.3*2 = 0.7 - below the cap, no saturation.
  grad200 = {"gradient": {"material_b": TIO2, "rate": 0.3, "f_start": 0.1},
               "optimize": False, "needle": False}
  st2, _ = stack_from_layers(
    [(SIO2, 200.0)], WL, {}, names=["g"], per_film_flags={"g": grad200})
  rows2 = [f for f in films(st2) if f["material"] == "g"]
  f_last = min(0.1 + 0.3 * (200.0 / 100.0), 1.0)
  assert np.array_equal(np.asarray(rows2[-1]["nk"]),
                        np.asarray(_oracle(SIO2, TIO2, f_last)))
  assert not any(np.array_equal(np.asarray(f["nk"]), np.asarray(pure_b))
                 for f in rows2), "the 200 nm film must NOT saturate"


def test_rate_capped_negative_rate_saturates_low():
  """A negative rate saturates at the f_min cap at the far end; with
  the default caps the tail is the f = 0 mixture (the EMA kernels at
  f = 0 return the host permittivity, so the tail rows are bitwise
  the oracle's f = 0 row)."""
  import warnings
  grad = {"gradient": {"material_b": SIO2, "rate": -0.25, "f_start": 0.9},
            "optimize": False, "needle": False}
  with warnings.catch_warnings():
    warnings.simplefilter("error")
    st, _ = stack_from_layers(
      [(TIO2, 400.0)], WL, {}, names=["g"], per_film_flags={"g": grad})
  rows = [f for f in films(st) if f["material"] == "g"]
  tail = _oracle(TIO2, SIO2, 0.0)
  assert any(np.array_equal(np.asarray(f["nk"]), np.asarray(tail))
             for f in rows), "the low tail must be the f = 0 mixture"
  # And the head row is f_start (inside the caps here).
  assert np.array_equal(np.asarray(rows[0]["nk"]),
                        np.asarray(_oracle(TIO2, SIO2, 0.9)))


def test_rate_capped_merge_twin_n2():
  """N2: the saturated film through the nk-keyed merge - the tail rows
  collapse, the span stays one contiguous run of the carrier material,
  and the simulated merit is unchanged to float precision.

  NOTE (CORRECTIONS to the plan's wording): the plan says the spectra
  are BITWISE equal before and after the merge; they are not - the
  merge folds two rows into one, and exp(i d1) exp(i d2) differs from
  exp(i (d1 + d2)) at the last ulp. The merge is a no-op to ~1e-15
  relative, and that is what this twin pins.
  """
  import warnings
  grad = {"gradient": {"material_b": TIO2, "rate": 0.25, "f_start": 0.1},
            "optimize": False, "needle": False}
  with warnings.catch_warnings():
    warnings.simplefilter("error")
    st, _ = stack_from_layers(
      [(SIO2, 400.0)], WL, {}, names=["g"], per_film_flags={"g": grad})
  ctx = _merit_ctx()
  m0 = ctx.evaluate_merit(st)
  graded0 = [f for f in films(st) if f["material"] == "g"]
  merged = st.merge_adjacent()
  assert merged > 0, "the saturated tail must merge"
  graded1 = [f for f in films(st) if f["material"] == "g"]
  assert len(graded1) < len(graded0), "the tail rows collapsed"
  # The span stays one contiguous run of the carrier material.
  seq = [f["material"] for f in films(st)]
  runs = [(k, len(list(g))) for k, g in itertools.groupby(seq)]
  graded_runs = [r for r in runs if r[0] == "g"]
  assert len(graded_runs) == 1, f"span split by the merge: {runs}"
  # The optics are untouched (to float precision).
  m1 = ctx.evaluate_merit(st)
  assert m1 == pytest.approx(m0, rel=1e-12, abs=1e-12)


def test_rate_capped_door_refusals():
  """The A9.2 refusal (two slope spellings) and the RateCapped native
  refusals surface at the door."""
  both = {"material_b": SIO2, "rate": 0.3, "f_end": 0.9}
  with pytest.raises(ValueError, match="same slope"):
    stack_from_layers([(TIO2, 50.0)], WL, {}, names=["H"],
                      per_film_flags={"H": {"gradient": both}})
  neither = {"material_b": SIO2}
  with pytest.raises(ValueError, match="requires 'f_end'"):
    stack_from_layers([(TIO2, 50.0)], WL, {}, names=["H"],
                      per_film_flags={"H": {"gradient": neither}})
  # Native RateCapped validation: NaN rate, non-positive ref, empty cap
  # window, zero rate (degenerate), f_* outside [0, 1].
  for bad, match in [
    ({"material_b": SIO2, "rate": float("nan")}, "finite"),
    ({"material_b": SIO2, "rate": 0.3, "ref_thickness": 0.0}, "ref_thickness"),
    ({"material_b": SIO2, "rate": 0.3, "f_min": 0.8, "f_max": 0.2}, "f_min"),
    ({"material_b": SIO2, "rate": 0.0}, "degenerate"),
    ({"material_b": SIO2, "rate": 0.3, "f_max": 1.5}, r"\[0, 1\]"),
  ]:
    with pytest.raises(ValueError, match=match):
      stack_from_layers(
        [(TIO2, 50.0)], WL, {}, names=["H"],
        per_film_flags={"H": {"gradient": dict(bad)}})


def test_rate_capped_optimize_keeps_the_profile_f1_7():
  """F1.7 completes F1.6's licence: optimize=true on a rate-driven
  gradient film also keeps the full profile (the mode reads absolute
  depth, so the span carries a recipe and is refreshed at the
  construction points; the LM sees ONE parameter). The rows are the
  hand-computed rate ramp at the authored thickness, same oracle as
  the saturation twin."""
  import warnings
  rate, f_start = 0.25, 0.1
  grad = {"gradient": {"material_b": TIO2, "rate": rate, "f_start": f_start}}
  with warnings.catch_warnings():
    warnings.simplefilter("error")  # must NOT warn
    st, _ = stack_from_layers(
      [(SIO2, 200.0)], WL, {}, names=["g"], per_film_flags={"g": grad})
  rows = [f for f in films(st) if f["material"] == "g"]
  assert len(rows) > 1, "the rate profile kept its sublayers"
  assert all(f["optimize"] for f in rows)
  n = len(rows)
  assert sum(f["thickness"] for f in rows) == pytest.approx(200.0)
  for i, f in enumerate(rows):
    z = 200.0 * i / (n - 1)
    f_i = min(max(f_start + rate * (z / 100.0), 0.0), 1.0)
    assert np.array_equal(np.asarray(f["nk"]), np.asarray(_oracle(SIO2, TIO2, f_i))), \
      f"row {i} (f={f_i})"

  # The stable posture still homogenizes loudly (EMA at f_mid).
  with pytest.warns(UserWarning, match="carries a gradient"):
    st2, _ = stack_from_layers(
      [(SIO2, 200.0)], WL, {}, names=["g"],
      per_film_flags={"g": {"gradient": {"material_b": TIO2, "rate": rate,
                                          "f_start": f_start},
                             "optimize": False, "needle": True}})
  fs2 = [f for f in films(st2) if f["material"] == "g"]
  assert len(fs2) == 1
  # R4: EMA at f_mid, the mode's own mid-depth value (z = t/2 = 100).
  f_mid = min(max(f_start + rate * (200.0 / 2.0 / 100.0), 0.0), 1.0)
  assert np.array_equal(np.asarray(fs2[0]["nk"]), np.asarray(_oracle(SIO2, TIO2, f_mid)))
  assert fs2[0]["thickness"] == pytest.approx(200.0)
