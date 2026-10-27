# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""A dual-band anti-reflection coating: a VIS band and a NIR band.

The coating must stay non-reflecting over the visible, 420-680 nm, and
more tightly over a NIR band, 1044-1084 nm. It is built from three
oxides on fused silica:

  SiO2   low index, Malitson's fused-silica dispersion
  Al2O3  medium index
  HfO2   high index

The Al2O3 and HfO2 dispersions are representative Cauchy fits for
evaporated films; put your own coater's fits here.

The demand is two targets in one merit: R = 0 over the VIS band with a
0.5 % tolerance, and R = 0 over the NIR band with 0.1 %, so a NIR point
weighs five times as much as a VIS point.

Two stages:

  1. **Synthesis.** Needle synthesis from 1100 nm of SiO2 over 250 nm of
     Al2O3, with all three oxides in the needle pool: every film accepts
     a needle of either other oxide, the needles themselves included,
     and each insertion is the best (position, material) pair of them
     all. The thick start gives the needles room: from a three-layer
     quarter-half-quarter seed the run stalls in a thin local minimum.
  2. **The deposition floor.** The same design re-optimized with every
     film held at 10 nm or more, and no new needles. Synthesis cannot run
     at that floor: a needle is born a few nanometres thick and grows,
     and with a 10 nm floor from the start no insertion survives its
     re-optimization. So the search runs at the default 2 nm floor, and
     the floor is applied afterwards; films that fall below it are
     removed and the rest re-optimized.

Run:  python examples/dual_band_ar.py
"""

from collections import Counter

import numpy as np

from navette.materials import MaterialSpec, evaluate
from navette.smatrix.smatrix import Request, ScatterMatrix
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis.pipeline import PipelineConfig, run_needle

# --- the problem ----------------------------------------------------------
VIS = np.arange(420.0, 680.1, 5.0)
NIR_CENTER = 1064.0
NIR = np.arange(NIR_CENTER - 20.0, NIR_CENTER + 20.1, 5.0)
WL = np.concatenate([VIS, NIR])
TOL_VIS, TOL_NIR = 0.005, 0.001

SIO2 = MaterialSpec("Sellmeier", dict(B1=0.6961663, C1=0.0684043 ** 2,
                                      B2=0.4079426, C2=0.1162414 ** 2,
                                      B3=0.8974794, C3=9.896161 ** 2))
AL2O3 = MaterialSpec("Cauchy", dict(A=1.615, B=0.0083, C=0.0))
HFO2 = MaterialSpec("Cauchy", dict(A=1.95, B=0.016, C=0.0004))
MATERIALS = {"SiO2": SIO2, "Al2O3": AL2O3, "HfO2": HFO2}
SUBSTRATE = SIO2                # fused silica

SEED_NM = {"SiO2": 1100.0, "Al2O3": 250.0}
MAX_FILMS = 30
FLOOR_NM = 10.0


def nk(material, wl=None):
    return evaluate(material, WL if wl is None else wl)


def demands():
    tc = TargetCollection()
    tc.add(SpectralTarget(VIS, np.zeros_like(VIS), np.full_like(VIS, TOL_VIS),
                          0.0, "s", "R"))
    tc.add(SpectralTarget(NIR, np.zeros_like(NIR), np.full_like(NIR, TOL_NIR),
                          0.0, "s", "R"))
    return tc


def synthesize():
    """Stage 1: needle synthesis at the default 2 nm floor."""
    seed = [(nk(SIO2), SEED_NM["SiO2"]), (nk(AL2O3), SEED_NM["Al2O3"])]
    cfg = PipelineConfig(max_macro_cycles=60, needles_per_cycle=2,
                         max_film_layers=MAX_FILMS,
                         max_total_thickness_nm=3000.0)
    return run_needle(layers=seed, names=["SiO2", "Al2O3"], targets=demands(),
                      angles_deg=[0.0], wavelengths=WL,
                      materials={m: nk(spec) for m, spec in MATERIALS.items()},
                      pipeline_config=cfg, substrate=(nk(SUBSTRATE), "sub"))


def apply_floor(result):
    """Stage 2: re-optimize with every film at FLOOR_NM or more.

    One macro cycle with no needles (an empty pool) is an optimization
    plus cleanup, and the cleanup's thin-layer sweep now runs at the
    floor. A film removed between two of one material leaves them
    touching, and the engine merges them: layers merge by index. The
    cycle repeats until the floor removes nothing more.
    """
    stack = films(result)
    while True:
        layers = [(nk(MATERIALS[name]), t) for name, _, t in stack]
        cfg = PipelineConfig(max_macro_cycles=1, needles_per_cycle=0,
                             clamp_min_nm=FLOOR_NM, max_film_layers=MAX_FILMS,
                             max_total_thickness_nm=3000.0)
        out = run_needle(layers=layers, targets=demands(), angles_deg=[0.0],
                         wavelengths=WL, materials={}, pipeline_config=cfg,
                         substrate=(nk(SUBSTRATE), "sub"))
        stack = films(out)
        if len(stack) == len(layers):
            return out


def films(result):
    """[(material name, nk on WL, thickness)], ambient side first.

    Needles are named "<material>#<n>" and the floor stage's films
    "film<i>", so each film is identified by which material's spectrum
    it carries; `check` then proves the reading right.
    """
    out = []
    for f in result["stack"].to_dict()["films"]:
        n = np.asarray(f["nk"])
        name = min(MATERIALS, key=lambda m: np.max(np.abs(nk(MATERIALS[m]) - n)))
        out.append((name, n, float(f["thickness"])))
    return out


def reflectance(stack_films, wl):
    """Rs at normal incidence of `stack_films` on the substrate, over `wl`."""
    rows = ([np.ones_like(wl)]
            + [nk(MATERIALS[name], wl) for name, _, _ in stack_films]
            + [nk(SUBSTRATE, wl)])
    d = [0.0] + [t for _, _, t in stack_films] + [0.0]
    sm = ScatterMatrix(np.vstack(rows), d, wavelengths=wl, angles=[0.0])
    return np.asarray(sm.compute(Request.RS)["Rs"]).ravel()


def figures(stack_films):
    """What the spec sheet would say, on the design grid."""
    R = reflectance(stack_films, WL)
    vis, nir = R[:len(VIS)], R[len(VIS):]
    d = [t for _, _, t in stack_films]
    return dict(vis_mean=float(vis.mean()), vis_max=float(vis.max()),
                nir_center=float(reflectance(stack_films, np.array([NIR_CENTER]))[0]),
                nir_max=float(nir.max()), films=len(d), total_nm=float(sum(d)),
                thinnest_nm=float(min(d)))


def check(result, stack_films):
    """Re-score the design by material name: it must reproduce the merit.

    If a film were matched to the wrong material, or the order read the
    wrong way round, this would not agree with the run.
    """
    R = reflectance(stack_films, WL)
    tol = np.concatenate([np.full(len(VIS), TOL_VIS), np.full(len(NIR), TOL_NIR)])
    merit = float(np.sum((R / tol) ** 2))
    assert abs(merit - result["final_mf"]) <= 1e-6 * max(1.0, result["final_mf"]), \
        (merit, result["final_mf"])


def main() -> None:
    s1 = synthesize()
    s2 = apply_floor(s1)
    f1, f2 = films(s1), films(s2)
    check(s1, f1)
    check(s2, f2)
    a, b = figures(f1), figures(f2)

    print(f"VIS {VIS[0]:.0f}-{VIS[-1]:.0f} nm (tol {TOL_VIS:.1%}), "
          f"NIR {NIR[0]:.0f}-{NIR[-1]:.0f} nm (tol {TOL_NIR:.1%}), fused silica")
    print(f"synthesis: {s1['termination']} after {len(s1['phases'])} cycles")
    print()
    print(f"  {'stage':<18}{'merit':>8}{'films':>7}{'nm':>7}{'thinnest':>10}"
          f"{'VIS mean':>10}{'VIS max':>9}{'R@1064':>9}")
    for label, res, m in (("needle, 2 nm", s1, a), (f"floor {FLOOR_NM:.0f} nm", s2, b)):
        print(f"  {label:<18}{res['final_mf']:>8.2f}{m['films']:>7}"
              f"{m['total_nm']:>7.0f}{m['thinnest_nm']:>10.1f}"
              f"{m['vis_mean']:>10.3%}{m['vis_max']:>9.3%}{m['nir_center']:>9.3%}")
    print()
    print("  the design, ambient side first:")
    for i, (name, n, t) in enumerate(f2, 1):
        print(f"  {i:>4}  {name:<6}{t:>8.1f} nm")
    counts = {m: sum(1 for name, _, _ in f2 if name == m) for m in MATERIALS}
    inserted = Counter(nr["insertion"]["material"] for p in s1["phases"]
                       for nr in p["needle_results"] if nr["insertion"])
    print()
    print("THREE MATERIALS. " + ", ".join(f"{c} {m}" for m, c in counts.items())
          + " films. The run inserted")
    print(", ".join(f"{n} {m}" for m, n in inserted.most_common())
          + " needles: each insertion is the best")
    print("(position, material) pair over every film, needles included, and every")
    print("oxide of another index than the film it lands in.")
    print()
    print(f"THE FLOOR. Lifting the minimum film from {a['thinnest_nm']:.1f} to "
          f"{FLOOR_NM:.0f} nm moved the merit")
    print(f"from {s1['final_mf']:.2f} to {s2['final_mf']:.2f}; that is the "
          "price of a design a coater can")
    print("actually lay down.")


if __name__ == "__main__":
    main()
