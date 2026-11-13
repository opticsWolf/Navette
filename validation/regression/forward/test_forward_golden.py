# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Forward-model golden baseline: every spectral and phase output, recomputed.

Each case in `_cases.CASES` is expanded through the structure layer, solved
for every request bit on three wavelength ranges (UV 250-400, VIS 400-800,
NIR 800-2000 nm) at seven angles (0-89 deg), and -- for single-block stacks
-- re-solved through the synthesis evaluator. Every array is compared with
`golden/forward_golden.npz`. The survey of what is covered, and why the
tolerances are what they are, is `docs/forward-golden.md`.

A failure lists the worst deviation per quantity for that case. Knobs:

  NAVETTE_GOLDEN_EXACT=1   demand bit-identity instead of the tolerances
  NAVETTE_GOLDEN_REPORT=1  print the worst deviation per quantity over all
                           cases at the end of the module

The baseline is rewritten only by `gen_forward_golden.py --write`, never by
this test.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _cases as C  # noqa: E402

NPZ = HERE / "golden" / "forward_golden.npz"
MANIFEST = HERE / "golden" / "manifest.json"
EXACT = os.environ.get("NAVETTE_GOLDEN_EXACT") == "1"
REPORT = os.environ.get("NAVETTE_GOLDEN_REPORT") == "1"

# quantity -> (worst deviation, case.range) over the module, for the report.
_WORST: dict = {}


@pytest.fixture(scope="module")
def manifest():
  return json.loads(MANIFEST.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def golden():
  with np.load(NPZ) as npz:
    yield {name: npz[name] for name in npz.files}


def test_manifest_matches_the_case_definitions(manifest):
  """The stored baseline describes exactly the cases this code defines."""
  assert manifest["format_version"] == C.FORMAT_VERSION
  assert manifest["angles_deg"] == list(C.ANGLES)
  for rng in C.RANGES:
    assert manifest["grids"][rng] == C.grid(rng).tolist(), rng
  assert sorted(manifest["cases"]) == sorted(C.CASES), \
    "case set changed; regenerate with gen_forward_golden.py --write"
  for name, case in C.CASES.items():
    stored = manifest["cases"][name]
    assert stored["definition"] == json.loads(json.dumps(case)), name
    assert stored["options"] == C.case_options(case), name


def test_baseline_file_is_the_one_the_manifest_names(manifest):
  sha = hashlib.sha256(NPZ.read_bytes()).hexdigest()
  assert sha == manifest["npz_sha256"], "forward_golden.npz edited outside the generator"


def _quantity(key: str) -> str:
  _, path, q = key.split(".", 2)
  return f"{path}.{q}"


@pytest.mark.parametrize("name", list(C.CASES))
def test_case_reproduces_the_baseline(name, manifest, golden):
  gold = C.unpack(golden, name, manifest["cases"][name]["layout"])
  got = C.evaluate(name)
  assert sorted(got) == sorted(gold), (
    f"{name}: keys differ -- new {sorted(set(got) - set(gold))}, "
    f"missing {sorted(set(gold) - set(got))}")
  for rng in C.RANGES:
    mask = C.request_mask(C.CASES[name], {"flags": got[f"{rng}.in.flags"]})
    assert mask == manifest["cases"][name]["request_mask"][rng], (name, rng)

  per_q: dict = {}
  failures = []
  for key in sorted(gold):
    dev = C.deviation(key, got[key], gold[key])
    tol = 0.0 if EXACT else C.tolerance(key, gold[key])
    q = _quantity(key)
    if dev > per_q.get(q, (-1.0, ""))[0]:
      per_q[q] = (dev, key.split(".", 1)[0], tol)
    if dev > tol:
      failures.append(key)
    if dev > _WORST.get(q, (-1.0, ""))[0]:
      _WORST[q] = (dev, f"{name}.{key.split('.', 1)[0]}")
  if failures:
    rows = sorted(per_q.items(), key=lambda kv: -kv[1][0])
    lines = [f"  {q:24s} {d:10.3e} (tol {t:.0e}, {rng}){'  <-- FAIL' if d > t else ''}"
             for q, (d, rng, t) in rows if d > 0.0]
    pytest.fail(f"{name}: {len(failures)} array(s) out of tolerance; worst per quantity:\n"
                + "\n".join(lines))


@pytest.mark.skipif(not REPORT, reason="set NAVETTE_GOLDEN_REPORT=1 for the summary")
def test_report_worst_deviation_per_quantity(capsys):
  with capsys.disabled():
    print("\nforward golden: worst deviation per quantity over all cases")
    for q, (d, where) in sorted(_WORST.items(), key=lambda kv: -kv[1][0]):
      print(f"  {q:24s} {d:10.3e}  {where}")
