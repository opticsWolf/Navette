# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Regression: state round-trips and validation (BUG-1, BUG-2 accept).

BUG-1: `from_state(get_state())` preserves every field — material names
included — for Layer, Group, Navette_Structure and Navette_Architect
(shared-block aliasing and material re-attachment included).
BUG-2: `validate()` never crashes; missing materials are reported.
"""

import numpy as np
import pytest

from navette.structure import (
  DictMaterialProvider,
  Group,
  Layer,
  Navette_Architect,
  Navette_Structure,
)

WL = np.array([1000.0])
MATS = DictMaterialProvider({
  "glass": np.full(1, 1.52 + 0j),
  "TiO2": np.full(1, 2.35 + 0.01j),
})


def _full_layer() -> Layer:
  return Layer(
    thickness=50.0, material_name="TiO2", coherent=False,
    roughness=2.0, rough_type=5, inhomogen=True, inh_delta=0.2,
    interface=True, interface_thickness=3.0, optimize=False,
    needle=False, layer_type=1,
  )


def test_layer_roundtrip_all_fields():
  original = _full_layer()
  back = Layer.from_state(original.get_state())
  assert back.get_state() == original.get_state()
  assert back.material == "TiO2"  # BUG-1: was '' before STRUCT-2


def test_layer_roundtrip_minimal():
  original = Layer(thickness=10.0, material_name="glass")
  back = Layer.from_state(original.get_state())
  assert (back.material, back.thickness) == ("glass", 10.0)


def test_v1_state_fixture_loads_and_expands_bit_identically():
  """F1.4/R5: the committed v1 fixture is the acceptance oracle.

  The fixture was dumped from a live build at 0.6.40 (schema_version 1)
  and committed BEFORE the range gate existed, so the "v1 stays
  readable" half of the contract is pinned by a file that predates it -
  not by a test written alongside the change. It must load and expand
  bit-identically to the live equivalent: every field that drives the
  expansion rides the state, so dict equality plus bitwise solver
  inputs is the full proof. Note the loaded state re-tags itself at the
  current SCHEMA_VERSION on write (the tag is checked, not stored), so
  this test is stable across the bump.
  """
  import json
  from pathlib import Path

  fixture = (Path(__file__).resolve().parents[2] / "fixtures" / "state"
             / "v1_architect.json")
  state = json.loads(fixture.read_text(encoding="utf-8"))
  assert state["schema_version"] == 1

  loaded = Navette_Architect.from_state(state, materials=MATS)

  layers = [
    Layer(thickness=100.0, material_name="TiO2", inhomogen=True,
          inh_delta=0.2, optimize=True, needle=False),
    Layer(thickness=50.0, material_name="glass"),
    Layer(thickness=30.0, material_name="TiO2", inhomogen=True,
          inh_mode={"RateCapped": {"rate": 0.05, "ref_thickness": 100.0,
                                    "cap": 0.3}}),
  ]
  live = Navette_Architect(materials=MATS)
  live.add_structure(Navette_Structure(
    layers, {"TiO2": Group("TiO2", n_factor=1.1)}, MATS))

  assert loaded.get_state() == live.get_state()
  a, b = loaded.get_solver_inputs(), live.get_solver_inputs()
  assert np.array_equal(np.asarray(a.thicknesses), np.asarray(b.thicknesses))
  assert np.array_equal(np.asarray(a.indices), np.asarray(b.indices))
  # And the graded expansion is real (not a degenerate 3-row stack).
  assert len(np.asarray(a.thicknesses)) > len(layers)


def test_newer_build_state_refused():
  """F1.4/R5: a state naming a NEWER build is refused.

  Written against the un-bumped point gate (which refuses version 3 as
  stale); the F1.4 range gate keeps the refusal and sharpens the reason
  to 'newer build' - the assertion here is refusal-first by design, so
  the oracle exists before the change that must satisfy it.
  """
  arch = Navette_Architect(materials=MATS)
  arch.add_structure(Navette_Structure([Layer(10.0, "glass")], {}, MATS))
  future = arch.get_state()
  future["schema_version"] = 999
  with pytest.raises(ValueError):
    Navette_Architect.from_state(future, materials=MATS)


def test_stale_schema_versions_refused():
  """F1.4: the gate is a readable RANGE, and the refusals name direction.

  The pre-F1.4 test asserted SCHEMA_VERSION - 1 is refused; after the
  bump that is v1, which the range MUST accept (the committed v1
  fixture is the oracle). Inverted deliberately: v1 accepted, v0
  refused as stale, a newer-build version refused with the newer-build
  reason, an untagged state still refused as malformed.
  """
  from navette.structure.types import MIN_READABLE_SCHEMA_VERSION, SCHEMA_VERSION
  layer_state = _full_layer().get_state()
  assert layer_state["schema_version"] == SCHEMA_VERSION
  # The readable past: a v1-tagged state (this test's pre-F1.4
  # 'SCHEMA_VERSION - 1' case) loads.
  older = dict(layer_state)
  older["schema_version"] = MIN_READABLE_SCHEMA_VERSION
  assert Layer.from_state(older).get_state() == layer_state
  untagged = dict(layer_state)
  del untagged["schema_version"]  # no past before v1: untagged is malformed
  with pytest.raises(ValueError):
    Layer.from_state(untagged)
  stale = dict(layer_state)
  stale["schema_version"] = MIN_READABLE_SCHEMA_VERSION - 1
  with pytest.raises(ValueError, match="stale"):
    Layer.from_state(stale)
  arch = Navette_Architect(materials=MATS)
  arch.add_structure(Navette_Structure([Layer(10.0, "glass")], {}, None))
  future = arch.get_state()
  future["schema_version"] = SCHEMA_VERSION + 999
  with pytest.raises(ValueError, match="newer build"):
    Navette_Architect.from_state(future, materials=MATS)


# Fingerprint: the exact serialized key set per entity, at the version it
# was recorded at. If this fails, the key set changed — classify FIRST:
# any key-set change                            -> bump SCHEMA_VERSION,
#   re-record the fingerprint below (F1.4 rewrote this policy: the
#   newer-writer hazard means even additive keys bump the version - the
#   readable PAST widens via MIN_READABLE_SCHEMA_VERSION instead);
# a key older writers leave implicit with a reconstruct-on-read default
#   (F1.3's inh_mode, F1.4's gradient)           -> the key list here
#   grows a conditional entry, and the plain key set stays byte-
#   identical to every earlier release.
FINGERPRINT = {
  "version": 2,
  "Layer": ["coherent", "inh_delta", "inhomogen", "interface",
              "interface_thickness", "layer_type", "material_name",
              "needle", "optimize", "rough_type", "roughness",
              "schema_version", "thickness"],
  # F1.4: `gradient` rides only when Some, so a no-gradient stack's v2
  # state is byte-identical to its v1 state apart from the tag. The
  # nested key set gets its own entry - a nested object would otherwise
  # have weaker protection than every top-level key.
  "GradientSpec": ["ema", "material_a", "material_b", "mode", "shape",
                    "sublayers"],
  "Group": ["group_name", "inh_delta_summand", "inh_delta_error_params",
              "inh_delta_error_type", "interface_error_params",
              "interface_error_type", "interface_summand", "k_error_params",
              "k_error_type", "k_factor", "n_error_params", "n_error_type",
              "n_factor", "roughness_error_params", "roughness_error_type",
              "roughness_summand", "schema_version", "thick_factor",
              "thick_summand", "thickness_error_params",
              "thickness_error_type", "error_mask", "optimization_mask"],
  "Navette_Structure": ["groups", "layers", "schema_version"],
  "Navette_Architect": ["blocks", "schema_version", "structures"],
}


def test_state_fingerprint():
  from navette.structure.types import SCHEMA_VERSION
  assert FINGERPRINT["version"] == SCHEMA_VERSION, \
    "fingerprint recorded at a different version — re-classify (see comment)"
  assert sorted(_full_layer().get_state()) == sorted(FINGERPRINT["Layer"])
  assert sorted(Group("g").get_state()) == sorted(FINGERPRINT["Group"])
  # F1.3, additive key: a RateCapped layer's state gains `inh_mode`
  # (only when the mode is not the Fixed default - a Fixed layer's key
  # set above is unchanged, byte-identical to every pre-F1.3 build).
  capped = Layer(50.0, "TiO2", inh_mode={"RateCapped": {"rate": 0.05,
                                                       "ref_thickness": 100.0,
                                                       "cap": 0.3}})
  assert sorted(capped.get_state()) ==     sorted(FINGERPRINT["Layer"] + ["inh_mode"])
  assert capped.get_state()["inh_mode"]["RateCapped"]["rate"] == 0.05
  back = Layer.from_state(capped.get_state())
  assert back.inh_mode == {"RateCapped": {"rate": 0.05,
                                          "ref_thickness": 100.0,
                                          "cap": 0.3}}
  # F1.4, additive key: a gradient layer's state gains `gradient`
  # (only when Some - the plain list above stays byte-identical to
  # every pre-F1.4 build). The layer is built through `from_state`
  # with a hand-written v2 dict: the constructor surface for gradient
  # is F1.5, but the STATE surface exists now and gets pinned now.
  graded_state = _full_layer().get_state()
  graded_state["gradient"] = {
    "material_a": "TiO2", "material_b": "glass",
    "ema": {"Bruggeman": {"max_iter": 100, "tol": 1e-9}},
    "mode": {"FixedSpan": {"f_start": 0.0, "f_end": 1.0}},
    "shape": "Linear", "sublayers": None,
  }
  graded = Layer.from_state(graded_state)
  gs = graded.get_state()
  assert sorted(gs) == sorted(FINGERPRINT["Layer"] + ["gradient"])
  assert sorted(gs["gradient"]) == sorted(FINGERPRINT["GradientSpec"])
  assert Layer.from_state(gs).get_state() == gs  # v2 -> v2 preserves it
  st = Navette_Structure([Layer(10.0, "glass")], {"glass": Group("glass")}, MATS)
  assert sorted(st.get_state()) == sorted(FINGERPRINT["Navette_Structure"])
  arch = Navette_Architect(materials=MATS)
  arch.add_structure(st)
  assert sorted(arch.get_state()) == sorted(FINGERPRINT["Navette_Architect"])


def test_group_roundtrip():
  original = Group("TiO2", thick_factor=1.1, n_factor=0.9, k_factor=1.2)
  original.error_mask = [1, 0, 1, 0, 0, 0]
  original.optimization_mask = [0, 1, 1, 1, 1, 1, 1]
  back = Group.from_state(original.get_state())
  assert back.get_state() == original.get_state()
  # Masks are independent copies, not aliased dicts/lists.
  back.error_mask[0] = 9
  assert original.error_mask[0] == 1


def test_structure_roundtrip():
  original = Navette_Structure(
    [_full_layer(), Layer(100.0, "glass")],
    {"TiO2": Group("TiO2", n_factor=1.1)},
    MATS,
  )
  back = Navette_Structure.from_state(original.get_state(), materials=MATS)
  assert [l.get_state() for l in back] == [l.get_state() for l in original]
  assert set(back.group_dict) == {"TiO2"}
  assert back.group_dict["TiO2"].n_factor == 1.1
  assert back.validate() == []


def test_architect_roundtrip_shared_aliasing():
  shared = Navette_Structure([Layer(50.0, "TiO2")], {}, MATS)
  original = Navette_Architect(materials=MATS)
  original.add_structure(shared, label="a")
  original.add_structure(shared, inverted=True, label="b")
  back = Navette_Architect.from_state(original.get_state(), materials=MATS)
  assert back.block_count == 2
  # Shared-block aliasing survives the trip (one state, two refs).
  assert back.blocks[0].structure is back.blocks[1].structure
  assert back.blocks[1].inverted is True
  assert back.validate() == []


def test_validate_missing_material_no_crash():
  arch = Navette_Architect(materials=MATS)
  arch.add_structure(Navette_Structure([Layer(50.0, "Nope")], {}, None))
  issues = arch.validate()  # BUG-2: AttributeError before STRUCT-6
  assert issues == ["Layer 0: Material 'Nope' not found in material provider."]


def test_validate_without_materials_skips_coverage():
  arch = Navette_Architect()
  arch.add_structure(Navette_Structure([Layer(50.0, "Nope")], {}, None))
  assert arch.validate() == []


def test_validate_catches_solver_blockers():
  """0.6.28: a negative thickness is refused at construction instead.

  The architect's solve gate is still the backstop and still raises; it is
  demonstrated here with an unresolvable material, which the layer gate has
  no way to judge.
  """
  import pytest

  with pytest.raises(ValueError, match="Negative thickness"):
    Layer(-5.0, "TiO2")

  bad = Navette_Structure([Layer(5.0, "NotInTheLibrary")], {}, MATS)
  assert any("not found" in i for i in bad.validate())
  arch = Navette_Architect(materials=MATS)
  arch.add_structure(bad)
  try:
    arch.get_solver_inputs()
  except ValueError:
    pass
  else:
    raise AssertionError("solve gate did not raise on invalid structure")
