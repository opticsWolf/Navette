# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Graded films in the synthesis pipeline: homogenize-with-warning by
default, pinned background on request (Rust `from_design` twins)."""
import numpy as np
import pytest

from navette.synthesis.pipeline import stack_from_layers

WL = np.linspace(400., 800., 8)
TIO2 = {"model": "Konstant", "params": {"n": 2.35, "k": 0.01}}
SIO2 = {"model": "Konstant", "params": {"n": 1.46}}


def films(stack):
  return stack.to_dict()["films"]


def test_graded_optimize_keeps_the_profile_and_homogenize_path_still_warns():
  """F1.6 licence item 1/2 (the one non-additive change): optimize=true
  on a profiled film used to mean "homogenize me"; it now keeps the
  profile as ONE scalable span, silently. The homogenize-with-warning
  path survives for the posture that always homogenized: optimize=false
  + needle=true (not background - background requires needle=false)."""
  # The F1.6 posture: full profile, no warning, one multi-row span.
  import warnings
  with warnings.catch_warnings():
    warnings.simplefilter("error")  # must NOT warn
    st, _ = stack_from_layers(
      [(TIO2, 50.0)], WL, {}, names=["TiO2"],
      per_film_flags={"TiO2": {"inhomogen": True, "inh_delta": 0.2}})
  fs = films(st)
  assert len(fs) > 1  # the profile kept its sublayers
  assert all(f["optimize"] for f in fs)  # every bulk row free

  # The stable homogenize posture still warns.
  with pytest.warns(UserWarning, match="homogeneous"):
    st2, _ = stack_from_layers(
      [(TIO2, 50.0)], WL, {}, names=["TiO2"],
      per_film_flags={"TiO2": {"inhomogen": True, "inh_delta": 0.2,
                                "optimize": False, "needle": True}})
  fs2 = films(st2)
  assert len(fs2) == 1  # single base-index row, profile dropped loudly
  assert fs2[0]["thickness"] == pytest.approx(50.0)
  assert not fs2[0]["optimize"] and fs2[0]["needle"]


def test_background_pins_profile_silently():
  import warnings
  with warnings.catch_warnings():
    warnings.simplefilter("error")  # background must NOT warn
    st, cmap = stack_from_layers(
      [(SIO2, 100.0), (TIO2, 50.0)], WL, {"TiO2": SIO2},
      names=["sub", "TiO2"],
      per_film_flags={"sub": {"inhomogen": True, "inh_delta": 0.2,
                                "optimize": False, "needle": False}})
  fs = films(st)
  subs = [f for f in fs if f["material"] == "sub"]
  assert len(subs) > 1  # profile expanded, not flattened
  assert sum(f["thickness"] for f in subs) == pytest.approx(100.0)
  assert all(not f["optimize"] and not f["needle"] for f in subs)
  design = [f for f in fs if f["material"] == "TiO2"]
  assert len(design) == 1 and design[0]["optimize"] and design[0]["needle"]
  # Pinned span survives merge (nk-keyed, not name-keyed).
  assert st.merge_adjacent() == 0
  assert len(films(st)) == len(fs)


def test_flat_stacks_stay_silent():
  import warnings
  with warnings.catch_warnings():
    warnings.simplefilter("error")
    st, _ = stack_from_layers([(TIO2, 50.0)], WL, {}, names=["TiO2"])
  assert len(films(st)) == 1

