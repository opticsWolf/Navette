# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""One AR coating, two surroundings, optimized jointly.

A coating is designed for the surface it will be measured on, and then it
gets laminated. The bare-surface design is not the laminated design: an
adhesive at n ~ 1.52 against the top film removes most of the index step
the AR was built around, so a stack that reaches 0.2 % bare can sit at 2 %
under glass. The usual fix is to design for whichever case matters more
and accept the other one.

A joint design does not remove that trade -- nothing can, the two surfaces
want different films -- it makes the trade explicit and lets the optimizer
price it. The films are defined ONCE, as a named design segment; each
environment is that segment plus whatever sits around it; and one
optimization moves one set of thicknesses against the demands of every
environment at the same time. The films cannot drift apart, because there
is only one of each.

What this script shows, in order:

  1. the bare-only design, scored in both environments -- the trade, measured;
  2. the joint design, scored in both environments;
  3. the same joint design expressed as a program DOCUMENT, to show that the
     two spellings are one request (identical merit bits).

Note on the numbers: this is a short run (three macro cycles, 21
wavelengths, a deliberately poor seed) against a 0.2 % tolerance, so
neither design finishes anywhere near that tolerance. What is worth
reading is the SHAPE of the two rows, not their absolute size -- a real
design would run longer on a denser grid.

Run:  python examples/multi_environment_ar.py
"""

import json
import struct
import tempfile
from pathlib import Path

import numpy as np

from navette.config import load_program
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec
from navette.synthesis.pipeline import (
    PipelineConfig,
    SmatrixContext,
    design_from_program,
    run_needle,
    stack_from_layers,
)

# --- the problem ----------------------------------------------------------
WL = np.linspace(450.0, 650.0, 21)
ANGLES = [0.0]

H = 2.35 + 0j      # TiO2-ish
L = 1.46 + 0j      # SiO2-ish
GLASS = 1.52 + 0j  # substrate, and the laminate above
TOL = 0.002        # 0.2 % reflectance, the demand's own scale

#: A four-layer QW-ish starting stack. Deliberately unoptimized -- the
#: point is what the optimizer does with it under one surrounding versus
#: two, not how good the seed was.
SEED = [(L, 90.0, "a"), (H, 60.0, "b"), (L, 120.0, "c"), (H, 40.0, "d")]

#: Needle seeds: each film may host the other material. Keyed by film NAME,
#: which on this surface is the design parameter's identity.
CONTRAST = {"a": H, "b": L, "c": H, "d": L}

CFG = PipelineConfig(max_macro_cycles=3, needles_per_cycle=1,
                     enable_cleanup=False, enable_inflate=False,
                     thin_layer_policy="clamp_up_final")


def demands(*names):
    """Zero reflectance across the band, once per named environment."""
    tc = TargetCollection()
    for name in names:
        tc.add(SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, TOL),
                              0.0, "s", "R", environment=name))
    return tc


def environments():
    """Bare surface, and the same films under a 300 nm laminate.

    Both entries reference the SAME design segment. That is the whole
    mechanism: `"coat"` is defined once, and an environment is an ordered
    list of things, one of which is a reference to it.
    """
    return [
        {"name": "bare", "stack": [{"design": "coat"}]},
        {"name": "laminated", "stack": [
            {"layers": [(GLASS, 300.0, "lam")]},   # fixed, not optimized
            {"design": "coat"},
        ]},
    ]


def final_films(result):
    """``[(nk, d_nm, name)]`` of the finished design, in order."""
    return [(np.asarray(f["nk"]), float(f["thickness"]), f["material"])
            for f in result["stack"].to_dict()["films"]]


def score_in(films, env_name):
    """One finished design, scored ALONE against one surrounding.

    The two runs' own merits are not comparable -- the joint run carries
    twice the residuals -- so the trade has to be measured the way a
    reviewer would measure it: build the stack this environment actually
    presents, score it against that environment's demands only, and
    report the two numbers side by side.
    """
    layers = list(films)
    if env_name == "laminated":
        # Broadcast: the films carry evaluated nk arrays, and `stack_from_layers`
        # wants every film on the run grid (it does not broadcast constants).
        layers = [(np.full(WL.shape, GLASS), 300.0, "lam")] + layers
    # A needle insertion splits one film into three, and the two outer
    # pieces keep the host's material name -- which is exactly right for
    # a DESIGN (they are one parameter's material) and wrong for a flat
    # nk table, where a name is a key. Number them for this scoring pass
    # only; the design itself is untouched.
    names = [f"{n}{i}" for i, (_, _, n) in enumerate(layers)]
    tc = TargetCollection()
    tc.add(SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, TOL),
                          0.0, "s", "R"))
    spec = build_merit_spec(tc)
    stack, _ = stack_from_layers([(nk, d) for nk, d, _ in layers], WL, {},
                                 names=names, substrate=(GLASS, "sub"))
    ctx = SmatrixContext(spec, np.asarray(ANGLES, dtype=np.float64), WL)
    return float(spec.merit(ctx.simulate(stack), 1e6))


def bits(x):
    return struct.pack("<d", float(x)).hex()


def main() -> None:
    print(__doc__.split("Run:")[0].strip().splitlines()[0])
    print()

    # --- 1. the trade, measured ------------------------------------------
    # Optimize against the BARE demands only, with the laminated
    # environment still present so the run reports both. The laminated
    # demands are simply absent from the target set, so nothing pulls on
    # that surface -- this is the "design for one case" baseline.
    bare_only = run_needle(
        targets=demands("bare"), angles_deg=ANGLES, wavelengths=WL,
        contrast=CONTRAST, design={"coat": SEED},
        environments=environments(),
        pipeline_config=CFG, substrate=(GLASS, "sub"))

    # --- 2. the joint design ---------------------------------------------
    joint = run_needle(
        targets=demands("bare", "laminated"), angles_deg=ANGLES,
        wavelengths=WL, contrast=CONTRAST, design={"coat": SEED},
        environments=environments(),
        pipeline_config=CFG, substrate=(GLASS, "sub"))

    # --- the trade, scored the same way for both designs -----------------
    print("merit of each finished design, scored alone in each surrounding")
    print("(lower is better; same demands, same tolerance, K = 1 both times)")
    print()
    print(f"  {'design':<22}{'bare':>14}{'laminated':>14}   films")
    for label, res in (("optimized bare only", bare_only),
                       ("optimized jointly", joint)):
        films = final_films(res)
        print(f"  {label:<22}{score_in(films, 'bare'):>14.2f}"
              f"{score_in(films, 'laminated'):>14.2f}"
              f"   {len(films)}")
    print()
    print("The bare-only design wins bare by a wide margin and is an order")
    print("of magnitude worse laminated. The joint design is worse bare and")
    print("much better laminated: its WORST case is the number that moved.")
    print("Which of the two is right is a decision, not a computation --")
    print("what a joint run removes is the option of making it by accident.")
    print()
    for label, res in (("bare only", bare_only), ("joint", joint)):
        d = ", ".join(f"{n}={t:.1f}" for _, t, n in final_films(res))
        print(f"  {label:<12} {d}")
    print()

    # --- 3. the same design, as a document -------------------------------
    # The keyword call above and this document are the same request. The
    # surface does not build a stack and hand it over: it shapes the
    # `DesignRequest` a program file already is, and one compiler consumes
    # both. So the merits agree bit for bit, not to four decimals.
    doc = {
        "schema_version": 2, "kind": "program", "name": "joint AR",
        "sections": {
            "materials": [
                {"name": "a", "code": "a", "model": "Konstant",
                 "params": {"n": 1.46, "k": 0.0}},
                {"name": "b", "code": "b", "model": "Konstant",
                 "params": {"n": 2.35, "k": 0.0}},
                {"name": "c", "code": "c", "model": "Konstant",
                 "params": {"n": 1.46, "k": 0.0}},
                {"name": "d", "code": "d", "model": "Konstant",
                 "params": {"n": 2.35, "k": 0.0}},
                {"name": "lam", "code": "lam", "model": "Konstant",
                 "params": {"n": 1.52, "k": 0.0}},
            ],
            # One segment, referenced by both environments below.
            "design": {"coat": {"layers": [
                {"material_code": "a", "thickness_nm": 90.0},
                {"material_code": "b", "thickness_nm": 60.0},
                {"material_code": "c", "thickness_nm": 120.0},
                {"material_code": "d", "thickness_nm": 40.0},
            ]}},
            "environments": [
                {"name": "bare", "stack": [{"design": "coat"}]},
                {"name": "laminated", "stack": [
                    {"layers": [{"material_code": "lam",
                                 "thickness_nm": 300.0}]},
                    {"design": "coat"},
                ]},
            ],
        },
    }
    # Note the four material codes a/b/c/d: two of them are the same
    # physical material. On this surface a film NAME is the design
    # parameter's identity across environments, so two independent SiO2
    # layers are two codes carrying identical tables -- not one code used
    # twice, which would tie them together.

    path = Path(tempfile.mkdtemp()) / "joint_ar.json"
    path.write_text(json.dumps(doc, indent=2), encoding="utf-8")
    prog = load_program(str(path), WL)
    design, envs = design_from_program(prog)

    from_doc = run_needle(
        targets=demands("bare", "laminated"), angles_deg=ANGLES,
        wavelengths=WL, contrast=CONTRAST, design=design, environments=envs,
        pipeline_config=CFG, substrate=(GLASS, "sub"))

    print("the same joint design, two spellings:")
    print(f"  keywords  {bits(joint['final_mf'])}  {joint['final_mf']!r}")
    print(f"  document  {bits(from_doc['final_mf'])}  {from_doc['final_mf']!r}")
    same = bits(joint["final_mf"]) == bits(from_doc["final_mf"])
    print(f"  -> {'identical bits' if same else 'DIFFERENT BITS'}")
    print()
    print(f"document written to {path}")


if __name__ == "__main__":
    main()
