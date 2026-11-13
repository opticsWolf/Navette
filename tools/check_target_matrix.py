# SPDX-License-Identifier: LGPL-3.0-or-later
"""CI sync guard: the supported-target matrix in the target-kinds doc.

The table under *What synthesis can target* in
``docs/spectralweave-target-kinds.md`` used to be written by hand from
``SmatrixContext::provides`` and the converter's label tables. Here every
cell is measured instead: each row's labels are compiled
(``build_merit_spec``) and scored (``SmatrixContext.evaluate_merit``) on

* a coherent stack (column "coherent stack"),
* the same coating around an incoherent plate (column "stack with a
  plate"; a refusal there is footnote 2, phase across a plate), and
* the coherent stack on an absorbing exit medium (a refusal there earns
  the row footnote 1).

A row whose labels split (some accepted, some not) is reported as drift,
not rendered: a row means one answer. Labels outside every row are probed
too, so a label that starts compiling (a new W8 target) fails the check
until it gets a row.

The table lives between ``<!-- target-matrix:begin -->`` and
``<!-- target-matrix:end -->``. ``--write`` regenerates it; without it,
exit 1 on any difference. Needs the built extension.
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DOC = ROOT / "docs" / "spectralweave-target-kinds.md"
BEGIN, END = "<!-- target-matrix:begin -->", "<!-- target-matrix:end -->"

#: (demand, label cell, [(spectral, pol, phase)]) - the doc's rows.
ROWS = [
    ("intensity, front", "`R`, `T` × s/p/u",
     [(sp, pol, False) for sp in ("R", "T") for pol in "spu"]),
    ("intensity, back (light from the exit side)", "`RB`, `TB` × s/p/u",
     [(sp, pol, False) for sp in ("RB", "TB") for pol in "spu"]),
    ("absorptance, front", "`A` × s/p/u",
     [("A", pol, False) for pol in "spu"]),
    ("absorptance, back", "`AB` × s/p/u",
     [("AB", pol, False) for pol in "spu"]),
    ("phase of one element, front", "`R`, `T` × s/p, `phase=True`",
     [(sp, pol, True) for sp in ("R", "T") for pol in "sp"]),
    ("phase of one element, back", "`RB`, `TB` × s/p, `phase=True`",
     [(sp, pol, True) for sp in ("RB", "TB") for pol in "sp"]),
    ("retardance", "`RetR`, `RetT`, `DeltaR`, `DeltaT`, u, `phase=True`",
     [(sp, "u", True) for sp in ("RetR", "RetT", "DeltaR", "DeltaT")]),
    ("intensity ratio", "`RatioR`, `RatioT`, `ExtR`, `ExtT`, `DiattR`, `DiattT`, u",
     [(sp, "u", False) for sp in ("RatioR", "RatioT", "ExtR", "ExtT", "DiattR", "DiattT")]),
    ("differential phase", "`PDts`, `PDtp`, `phase=True`",
     [("PDts", "s", True), ("PDtp", "p", True)]),
]
#: Colour: `ColorTarget` curves.
COLOR_ROW = ("colour (`ColorTarget`)", "front `R` / `T`",
             ["Rs", "Rp", "Ru", "Ts", "Tp", "Tu"])

#: Labels the engine knows that no row lists: each must stay refused.
UNLISTED = [("A", "s", True), ("AB", "u", True), ("R", "u", True),
            ("PDrs", "s", True), ("PDrp", "p", True),
            ("RatioR", "u", True), ("DiattT", "s", False)]

WL = np.linspace(400.0, 700.0, 16)
ANG = np.array([30.0])


def _stacks():
  from navette import _smatrix as S

  def layer(name, n, d, **kw):
    return S.LayerSpec(name, np.full(WL.shape, n, dtype=complex), d, **kw)

  def half(name, n):
    return layer(name, n, 0.0, optimize=False, needle=False)

  coat = [layer("H", 2.35 + 0.0j, 80.0), layer("L", 1.46 + 0.0j, 100.0)]
  plate = layer("G", 1.52 + 0.0j, 1.0e6, coherent=False, optimize=False, needle=False)
  air, glass = half("air", 1.0), half("glass", 1.52)
  return {
      "coherent": S.DesignStack(air, glass, coat),
      "plate": S.DesignStack(air, half("air", 1.0), coat + [plate, layer("L", 1.46, 90.0)]),
      "absorbing exit": S.DesignStack(air, half("metal", 0.2 + 3.0j), coat),
  }


def _scores(make, stack) -> bool:
  """True when the target `make()` builds compiles and scores on `stack`."""
  from navette import _smatrix as S
  from navette.spectralweave.target import TargetCollection
  from navette.synthesis import build_merit_spec
  tc = TargetCollection()
  try:
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      tc.add(make())
      spec = build_merit_spec(tc)
      ctx = S.SmatrixContext(spec, WL, angles_deg=ANG)
      return bool(np.isfinite(ctx.evaluate_merit(stack)))
  except (ValueError, TypeError):
    return False


def _spectral(sp, pol, phase):
  from navette.spectralweave.target import SpectralTarget
  return lambda: SpectralTarget(WL, np.zeros(WL.size), np.full(WL.size, 0.1), pol, sp,
                        kind="e", angle_deg=float(ANG[0]), phase=phase)


def _color(curve):
  from navette.spectralweave.target import ColorTarget
  return lambda: ColorTarget(curve=curve, angle_deg=float(ANG[0]))


def _cell(name, results, problems, yes="yes", no="no"):
  if all(results):
    return yes
  if not any(results):
    return no
  problems.append(f"{name}: its labels split ({sum(results)}/{len(results)} accepted)")
  return "split"


def render(problems: list) -> str:
  stacks = _stacks()
  lines = ["| demand | labels (`spectral`, pol.) | coherent stack | stack with a plate |",
           "|---|---|---|---|"]
  rows = [(d, lab, [_spectral(*c) for c in cs]) for d, lab, cs in ROWS]
  rows.append((COLOR_ROW[0], COLOR_ROW[1], [_color(c) for c in COLOR_ROW[2]]))
  for demand, labels, targets in rows:
    coh = [_scores(t, stacks["coherent"]) for t in targets]
    pla = [_scores(t, stacks["plate"]) for t in targets]
    exit_ok = [_scores(t, stacks["absorbing exit"]) for t in targets]
    note = "" if all(exit_ok) else "¹"
    if any(exit_ok) and not all(exit_ok):
      problems.append(f"{demand}: the absorbing exit splits its labels")
    # A cell that is already "no" needs no exit footnote.
    c = _cell(demand, coh, problems, yes="yes" + note)
    p = _cell(demand, pla, problems, yes="yes" + note, no="no²")
    lines.append(f"| {demand} | {labels} | {c} | {p} |")
  for c in UNLISTED:
    if _scores(_spectral(*c), stacks["coherent"]):
      problems.append(f"{c} compiles and scores but no row lists it")
  return "\n".join(lines)


def main() -> int:
  ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
  ap.add_argument("--write", action="store_true", help="regenerate the table in the doc")
  args = ap.parse_args()
  sys.stdout.reconfigure(encoding="utf-8")  # the table holds × and ¹
  text = DOC.read_text(encoding="utf-8")
  if text.count(BEGIN) != 1 or text.count(END) != 1:
    print(f"check_target_matrix: {DOC.name} needs exactly one {BEGIN} ... {END} block")
    return 1
  head, rest = text.split(BEGIN)
  old, tail = rest.split(END)
  problems: list = []
  table = render(problems)
  if problems:
    print("check_target_matrix: the measured matrix is inconsistent:")
    for p in problems:
      print("  -", p)
    return 1
  new = f"\n{table}\n"
  if args.write:
    DOC.write_text(head + BEGIN + new + END + tail, encoding="utf-8", newline="")
    print(f"check_target_matrix: wrote {DOC.relative_to(ROOT)}")
    return 0
  if old.strip() != table.strip():
    print(f"check_target_matrix: {DOC.relative_to(ROOT)} is out of step with the code.")
    print("Measured table:\n" + table)
    print("Run: python tools/check_target_matrix.py --write")
    return 1
  print("check_target_matrix: the supported-target matrix matches the code")
  return 0


if __name__ == "__main__":
  sys.exit(main())
