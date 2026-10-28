#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Regenerate the forward-model golden baseline -- on purpose only.

    python validation/regression/forward/gen_forward_golden.py          # dry run
    python validation/regression/forward/gen_forward_golden.py --write  # rewrite

The dry run recomputes every case and prints the worst deviation from the
stored baseline per case; nothing is written. `--write` rewrites
`golden/forward_golden.npz` and `golden/manifest.json` (case definitions,
grids, request masks, block layouts, Navette version, git commit, file
hash). Rewrite only after an intended change to the forward model, and say
why in the commit; a performance change must pass against the baseline
unchanged. Needs a `--release` extension (`maturin develop --release`).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _cases as C  # noqa: E402

GOLDEN_DIR = HERE / "golden"
NPZ = GOLDEN_DIR / "forward_golden.npz"
MANIFEST = GOLDEN_DIR / "manifest.json"


def _git(*args: str) -> str:
  try:
    return subprocess.run(["git", *args], cwd=HERE, capture_output=True, text=True,
                          check=True).stdout.strip()
  except (OSError, subprocess.CalledProcessError):
    return "unknown"


def compute_all() -> dict:
  """`{case: flat evaluate() dict}`."""
  return {name: C.evaluate(name) for name in C.CASES}


def build_manifest(flats: dict, layouts: dict, sha256: str) -> dict:
  import navette

  cases = {}
  for name, case in C.CASES.items():
    masks = {}
    for rng in C.RANGES:
      rows = {k: flats[name][f"{rng}.in.{k}"]
              for k in ("indices", "thicknesses", "flags", "rough_types", "rough_vals")}
      masks[rng] = C.request_mask(case, rows)
    cases[name] = {
      "definition": case,
      "options": C.case_options(case),
      "request_mask": masks,
      "layout": layouts[name],
    }
  status = _git("status", "--porcelain", "--untracked-files=no")
  return {
    "format_version": C.FORMAT_VERSION,
    "navette_version": navette.__version__,
    "git_commit": _git("rev-parse", "HEAD"),
    "git_dirty": bool(status) if status != "unknown" else None,
    "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
    "numpy_version": np.__version__,
    "python_version": platform.python_version(),
    "platform": platform.platform(),
    "npz_sha256": sha256,
    "grids": {rng: C.grid(rng).tolist() for rng in C.RANGES},
    "angles_deg": list(C.ANGLES),
    "tolerances": {"abs": C.TOL_ABS, "ratio_and_angle": C.TOL_RATIO,
                   "dispersion_rel": C.TOL_DISP_REL},
    "cases": cases,
  }


def report(flats: dict) -> int:
  if not (NPZ.exists() and MANIFEST.exists()):
    print(f"no baseline in {GOLDEN_DIR}; run with --write to create it")
    return 1
  manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
  bad = 0
  with np.load(NPZ) as npz:
    for name, flat in flats.items():
      if name not in manifest["cases"]:
        print(f"{name:26s} NEW (not in baseline)")
        bad += 1
        continue
      gold = C.unpack(npz, name, manifest["cases"][name]["layout"])
      worst, worst_key, n_exact = 0.0, "", 0
      for k in sorted(set(flat) | set(gold)):
        if k not in flat or k not in gold:
          print(f"  {name}.{k}: {'MISSING' if k in gold else 'NEW'}")
          bad += 1
          continue
        dev = C.deviation(k, flat[k], gold[k])
        n_exact += dev == 0.0
        bad += dev > C.tolerance(k, gold[k])
        if dev > worst:
          worst, worst_key = dev, k
      print(f"{name:26s} worst {worst:9.2e} ({worst_key or '-'}) bit-exact {n_exact}/{len(gold)}")
    for name in sorted(set(manifest["cases"]) - set(flats)):
      print(f"{name:26s} MISSING (in baseline, not defined)")
      bad += 1
  print("OK: baseline reproduced" if bad == 0 else f"{bad} out-of-tolerance or unmatched arrays")
  return 0 if bad == 0 else 1


def main() -> int:
  ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
  ap.add_argument("--write", action="store_true", help="rewrite the baseline files")
  args = ap.parse_args()
  flats = compute_all()
  if not args.write:
    return report(flats)
  members, layouts = {}, {}
  for name, flat in flats.items():
    m, layouts[name] = C.pack(name, flat)
    members.update(m)
  GOLDEN_DIR.mkdir(exist_ok=True)
  np.savez_compressed(NPZ, **members)
  sha = hashlib.sha256(NPZ.read_bytes()).hexdigest()
  manifest = build_manifest(flats, layouts, sha)
  MANIFEST.write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8", newline="\n")
  size = NPZ.stat().st_size / 1024
  print(f"wrote {NPZ.name} ({len(members)} members, {size:.0f} KiB) and {MANIFEST.name}")
  print(f"navette {manifest['navette_version']} @ {manifest['git_commit'][:10]}"
        f"{' (dirty)' if manifest['git_dirty'] else ''}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
