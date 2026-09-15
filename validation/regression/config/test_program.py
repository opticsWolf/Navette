# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Phase 1 program documents: envelope gate, partial + full restore,
prefix namespaces, legacy flat files."""
import pathlib

import numpy as np
import pytest

from navette.config import load_document, load_program
from navette.config.program import PROGRAM_SCHEMA_VERSION

WL = np.array([400., 600., 800.])
EXAMPLE = "src/navette/config/example_program.yaml"


def test_version_gate(tmp_path):
  import json
  from navette._structure import gate_document
  def gate(doc):
    p = tmp_path / "d.json"
    p.write_text(json.dumps(doc), encoding="utf-8")
    return load_document(str(p))
  with pytest.raises(ValueError, match="unsupported"):
    gate({"kind": "materials", "schema_version": 99, "materials": []})
  with pytest.raises(ValueError, match="unsupported"):
    gate({"kind": "materials", "materials": []})  # missing
  with pytest.raises(ValueError, match="unknown"):
    gate({"kind": "pipeline", "schema_version": 1})
  kind, name, payload = gate({"schema_version": 1, "kind": "groups", "groups": []})
  assert (kind, name, payload) == ("groups", None, [])


def test_full_restore():
  prog = load_program(EXAMPLE, WL)
  assert prog.name == "AR demo"
  assert prog.materials.contains("L") and prog.materials.contains("H")
  assert set(prog.groups) == {"H"}
  assert set(prog.structures) == {"ar"}
  st = prog.structures["ar"]
  assert [layer.material for layer in st.layer_list] == ["L", "H"]
  assert prog.architect is not None
  assert len(prog.architect) == 1
  assert prog.architect.blocks[0].label == "main"
  # Graded layer survived the trip with its flags.
  graded = st.layer_list[1]
  assert graded.inhomogen and not graded.optimize and not graded.needle


def test_prefix_namespace():
  a = load_program(EXAMPLE, WL)
  b = load_program(EXAMPLE, WL, prefix="run2_")
  assert set(b.structures) == {"run2_ar"}
  assert b.materials.contains("run2_L")
  assert [layer.material for layer in
          b.structures["run2_ar"].layer_list] == ["run2_L", "run2_H"]
  assert set(b.groups) == {"run2_H"}
  # Unprefixed load is unaffected.
  assert set(a.structures) == {"ar"}


def test_partial_restore_and_context():
  kind, name, payload = load_document(EXAMPLE)
  assert kind == "program" and name == "AR demo"
  import navette.config.program as P
  mats = P.load_materials(payload["materials"], WL)
  assert mats.contains("H")
  groups = P.load_groups(payload["groups"])
  assert set(groups) == {"H"}
  # Structures resolve against a context provider (file section absent).
  single = {"label": "solo", "layers": [{"material_code": "H", "thickness_nm": 10.0}]}
  st = P.load_structure(single, mats, groups)
  assert [layer.material for layer in st.layer_list] == ["H"]
  with pytest.raises(KeyError, match="needs materials"):
    P.load_structure(single, None, groups)


def test_missing_ref_names_section():
  import navette.config.program as P
  mats = P.load_materials([
    {"name": "H", "code": "H", "model": "Konstant", "params": {"n": 2.0, "k": 0.0}}], WL)
  bad = {"label": "x", "layers": [{"material_code": "NOPE", "thickness_nm": 5.0}]}
  with pytest.raises(KeyError, match="NOPE"):
    P.load_structure(bad, mats, {})


def test_legacy_flat_files(tmp_path):
  import yaml
  flat_mat = {"materials": [
    {"name": "H", "code": "H", "model": "Konstant", "params": {"n": 2.0, "k": 0.0}}]}
  p = tmp_path / "m.yaml"
  p.write_text(yaml.safe_dump(flat_mat))
  kind, _, payload = load_document(str(p))
  assert kind == "materials"
  import navette.config.program as P
  assert P.load_materials(payload, WL).contains("H")
  flat_stack = {"layers": [{"material_code": "H", "thickness_nm": 5.0}]}
  q = tmp_path / "s.yaml"
  q.write_text(yaml.safe_dump(flat_stack))
  kind, _, payload = load_document(str(q))
  assert kind == "structure"
  prog = load_program(str(q), WL,
                      context={"materials": P.load_materials(flat_mat["materials"], WL)})
  assert set(prog.structures) == {"stack"}


def _assert_shared_blocks(prog, label):
  named = prog.structures[label]
  block_shell = prog.architect.blocks[0].structure
  assert block_shell._inner.core_id() == named._inner.core_id()
  # Edits propagate through the shared handle.
  before = prog.architect.get_global_layer_count()
  named._inner.append_layer(named._inner.layer_list[0])
  try:
    assert prog.architect.get_global_layer_count() == before + 1
  finally:
    named._inner.remove_layer(len(named._inner.layer_list) - 1)


def test_restored_blocks_alias_structures():
  # Whole-document native path.
  _assert_shared_blocks(load_program(EXAMPLE, WL), "ar")


def test_context_path_blocks_alias_structures():
  # Section-wise context path keeps the same invariant.
  from navette.config.program import load_materials
  kind, name, payload = load_document(EXAMPLE)
  mats = load_materials(payload["materials"] if kind == "program"
                        else payload, WL)
  prog = load_program(EXAMPLE, WL, context={"materials": mats})
  _assert_shared_blocks(prog, "ar")


# -- F2.4: the envelope became a range, and gained two sections --------------

def _stamp(tmp_path, version, extra=None):
  """`example_program.yaml`, re-stamped, optionally with extra sections."""
  import json
  import yaml
  doc = yaml.safe_load(pathlib.Path(EXAMPLE).read_text(encoding="utf-8"))
  doc["schema_version"] = version
  if extra:
    doc["sections"].update(extra)
  p = tmp_path / f"v{version}.json"
  p.write_text(json.dumps(doc), encoding="utf-8")
  return str(p)


def _fingerprint(prog):
  """Everything a restored program asserts about ITSELF, as bits.

  Thicknesses and evaluated indices, not a summary: a resample that moved
  one sample, or a flag that defaulted the other way, is invisible to
  every comparison coarser than this one.
  """
  import struct
  out = [prog.name or "", sorted(prog.groups), sorted(prog.structures)]
  for label in sorted(prog.structures):
    st = prog.structures[label]
    for layer in st.layer_list:
      out.append((layer.material,
                  struct.pack("<d", float(layer.thickness)).hex(),
                  bool(layer.inhomogen), bool(layer.optimize),
                  bool(layer.needle)))
  for code in sorted(getattr(prog.materials, "_dict", {})):
    nk = prog.materials.get_nk(code)
    out.append((code, b"".join(struct.pack("<dd", z.real, z.imag)
                               for z in nk).hex()))
  return out


def test_a_v1_program_assembles_bit_identically_under_the_v2_build(tmp_path):
  """The oldest readable envelope is READ, not merely accepted.

  A range gate that let a v1 document through but assembled it differently
  would be worse than one that refused it: the refusal is visible and the
  drift is not. The fixture is stamped v1 and v2 and compared bit for bit
  -- the bump added sections, and adding a section must not change what a
  document that has none restores to.
  """
  a = load_program(_stamp(tmp_path, 1), WL)
  b = load_program(_stamp(tmp_path, 2), WL)
  assert _fingerprint(a) == _fingerprint(b)
  # And the real v1 fixture on disk, unedited, agrees with both.
  assert _fingerprint(load_program(EXAMPLE, WL)) == _fingerprint(a)


def test_the_envelope_gate_is_a_range_with_two_sided_messages(tmp_path):
  """Too old and too new are different problems, so they read differently.

  Mirrors the state gate (F1.4). A single "unsupported" for both ends
  leaves the user unable to tell whether to upgrade navette or regenerate
  the file, which are the only two remedies and are not interchangeable.
  """
  from navette.config.program import MIN_READABLE_PROGRAM_SCHEMA_VERSION
  assert (MIN_READABLE_PROGRAM_SCHEMA_VERSION, PROGRAM_SCHEMA_VERSION) == (1, 2)
  with pytest.raises(ValueError, match="newer build"):
    load_program(_stamp(tmp_path, PROGRAM_SCHEMA_VERSION + 1), WL)
  with pytest.raises(ValueError, match=r"reads 1..=2"):
    load_program(_stamp(tmp_path, MIN_READABLE_PROGRAM_SCHEMA_VERSION - 1), WL)


def test_an_unknown_section_refuses_by_name(tmp_path):
  """Silent-drop is the bug the bump exists to work around; close it too.

  Every section lookup is a bare `get`, so an unrecognised name used to
  vanish. That is exactly how a v2 program loses its environments on a v1
  build -- and it would have happened again for the next section added.
  """
  with pytest.raises(ValueError, match="enviroments"):
    load_program(_stamp(tmp_path, 2, {"enviroments": []}), WL)


def test_a_design_section_needs_an_environment_to_reach_it(tmp_path):
  """An unreferenced design segment is unreachable, not merely unused."""
  with pytest.raises(ValueError, match="environments"):
    load_program(_stamp(tmp_path, 2, {
      "design": {"coat": {"layers": [
        {"material_code": "L", "thickness_nm": 100.0}]}}}), WL)


def test_the_two_new_sections_round_trip(tmp_path):
  """Carried verbatim: the loader checks shape, the compiler checks meaning.

  Building them here would mean building them twice -- once for the file
  path and once for `run_needle(design=...)` -- with two chances to
  disagree about what a segment is.
  """
  path = _stamp(tmp_path, 2, {
    "design": {"coat": {"layers": [
      {"material_code": "L", "thickness_nm": 100.0},
      {"material_code": "H", "thickness_nm": 80.0}]}},
    "environments": [{"name": "bare", "stack": [{"design": "coat"}]}]})
  prog = load_program(path, WL)
  assert list(prog.design) == ["coat"]
  assert [r["material_code"] for r in prog.design["coat"]["layers"]] == ["L", "H"]
  assert [e["name"] for e in prog.environments] == ["bare"]
  assert prog.environments[0]["stack"][0]["design"] == "coat"
