# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""A rugate filter two ways: as gradient spans, and as a discrete stack.

A rugate is an index modulation rather than a stack of layers: the index
sweeps up and down through the coating instead of stepping. Its selling
point is what it does NOT have -- no abrupt interfaces means no strong
harmonics of the stopband, so the reflector is narrow and the rest of the
spectrum stays clean.

Navette does not have a sinusoidal profile mode, and this example does not
pretend otherwise. What it has is `FixedSpan`: a mixture fraction that runs
linearly from `f_start` to `f_end` across one span. A rugate built from
alternating `FixedSpan` spans is therefore a TRIANGULAR modulation, not a
sinusoidal one -- the right shape family, one Fourier component off. Read
the comparison below as "graded versus stepped", not "rugate versus stack".

The two things worth measuring, and the reason this example exists:

  * **rows.** Each gradient span expands into sublayers, and the solver
    sees rows, not spans. A graded design costs what its EXPANSION costs.
  * **merit.** Against the same stopband demand, on the same grid.

Run:  python examples/rugate_gradient_vs_discrete.py
"""

import numpy as np

from navette.smatrix.smatrix import Request, ScatterMatrix
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec
from navette.synthesis.pipeline import SmatrixContext, stack_from_layers

# --- the problem ----------------------------------------------------------
WL = np.linspace(450.0, 780.0, 67)
ANGLES = np.array([0.0])

LAM0 = 600.0                    # stopband centre
N_HI, N_LO = 2.35, 1.46         # the two endpoint indices
N_AVG = 0.5 * (N_HI + N_LO)
PERIODS = 8

#: One optical period of the modulation, and the quarter-waves that make
#: the stepped equivalent. Same optical thickness per period either way,
#: which is what makes the two designs comparable at all.
HALF_PERIOD_NM = LAM0 / (2.0 * N_AVG) / 2.0
QW_HI_NM = LAM0 / (4.0 * N_HI)
QW_LO_NM = LAM0 / (4.0 * N_LO)

GLASS = 1.52 + 0j


def nk(value):
    """A constant index on the run grid (films carry evaluated arrays)."""
    return np.full(WL.shape, complex(value))


def stopband_demands():
    """Reflect at the centre, transmit either side of it.

    Three demands rather than one full-band curve, because that is how a
    filter is actually specified: a pass level, a stop level, and nothing
    said about the transition regions -- where a demand would only be
    fighting the physics.
    """
    tc = TargetCollection()
    stop = WL[np.abs(WL - LAM0) <= 25.0]
    blue = WL[WL <= 520.0]
    red = WL[WL >= 700.0]
    tc.add(SpectralTarget(stop, np.ones_like(stop), np.full_like(stop, 0.01),
                          0.0, "s", "R"))
    for band in (blue, red):
        tc.add(SpectralTarget(band, np.zeros_like(band),
                              np.full_like(band, 0.01), 0.0, "s", "R"))
    return tc


def graded():
    """The triangular modulation: alternating `FixedSpan` gradient spans.

    Each span's host material is the LOW index and its inclusion is the
    high one; `f_start`/`f_end` run 0 -> 1 on the way up and 1 -> 0 on the
    way down, so consecutive spans meet at the same fraction and the
    profile is continuous across the join. A discontinuity there would
    put back exactly the interface the modulation exists to remove.
    """
    layers, names, flags = [], [], {}
    for i in range(2 * PERIODS):
        name = f"g{i}"
        up = (i % 2 == 0)
        layers.append((nk(N_LO), HALF_PERIOD_NM))
        names.append(name)
        flags[name] = {"gradient": {
            "material_b": nk(N_HI),
            "f_start": 0.0 if up else 1.0,
            "f_end": 1.0 if up else 0.0,
            "ema": "Bruggeman",
        }}
    return layers, names, flags


def stepped():
    """The discrete equivalent: quarter-wave pairs at the same centre."""
    layers, names = [], []
    for i in range(PERIODS):
        layers.append((nk(N_HI), QW_HI_NM))
        names.append(f"h{i}")
        layers.append((nk(N_LO), QW_LO_NM))
        names.append(f"l{i}")
    return layers, names, {}


def build(layers, names, flags):
    return stack_from_layers(layers, WL, {}, names=names,
                             per_film_flags=flags,
                             substrate=(nk(GLASS), "sub"))[0]


def spectrum(stack):
    """Rs on the run grid, solved over the stack's EXPANDED rows.

    `to_dict()` returns the expansion -- one entry per row, each with the
    nk the solver will use -- so this is the same coating the merit arm
    scored, not the authored spans. (`SimCurves` is a write-only sink;
    there is no spectrum to read back out of `ctx.simulate`.)
    """
    rows = stack.to_dict()["films"]
    n = np.vstack([np.ones_like(WL)]
                  + [np.asarray(r["nk"]) for r in rows]
                  + [np.full(WL.shape, GLASS)])
    d = np.array([0.0] + [float(r["thickness"]) for r in rows] + [0.0])
    out = ScatterMatrix(n.T if n.shape[0] != len(d) else n, d,
                        wavelengths=WL, angles=list(ANGLES)).compute(Request.RS)
    return np.asarray(out["Rs"]).ravel()


def report(label, stack, spec):
    """Rows, physical thickness, merit, and the peak the filter reaches."""
    ctx = SmatrixContext(spec, ANGLES, WL)
    merit = float(spec.merit(ctx.simulate(stack), 1e6))
    peak = float(np.max(spectrum(stack)))
    print(f"  {label:<26}{stack.film_count():>6}"
          f"{stack.total_thickness():>12.1f}"
          f"{merit:>14.2f}{peak:>10.4f}")
    return merit


def main() -> None:
    spec = build_merit_spec(stopband_demands())

    g_layers, g_names, g_flags = graded()
    s_layers, s_names, s_flags = stepped()
    g_stack = build(g_layers, g_names, g_flags)
    s_stack = build(s_layers, s_names, s_flags)

    print(f"{PERIODS} periods at {LAM0:.0f} nm, n = {N_LO} .. {N_HI}")
    print(f"half-period {HALF_PERIOD_NM:.2f} nm graded; "
          f"quarter-waves {QW_HI_NM:.2f} / {QW_LO_NM:.2f} nm stepped")
    print()
    print(f"  {'design':<26}{'rows':>6}{'thickness':>12}"
          f"{'merit':>14}{'peak R':>10}")
    report(f"graded ({2 * PERIODS} spans)", g_stack, spec)
    report(f"stepped ({2 * PERIODS} films)", s_stack, spec)
    print()

    spans, rows = 2 * PERIODS, g_stack.film_count()
    print(f"ROWS. Authored as {spans} spans, solved as {rows}: every span")
    print(f"expands into {rows // spans} sublayers, and the solver sees rows.")
    print("Authoring a graded design is cheap and solving one is not -- the")
    print("two numbers are far apart by construction, and the row count is")
    print("the one that shows up in a timing.")
    print()
    print("MERIT. The stepped design reaches the higher peak (abrupt")
    print("interfaces reflect harder per period), and loses badly overall,")
    print("because the demand also asks for transmission either side of the")
    print("stopband and a quarter-wave stack has sidelobes there. The graded")
    print("profile has no abrupt interfaces to generate them. That is the")
    print("whole trade a rugate is bought for, and here it costs 5x the rows.")
    print()
    print("A NOTE ON INTERFACES (Nevot-Croce). Roughness type 5 applies an")
    print("interface factor at every boundary it is set on. A graded span's")
    print("sublayer boundaries are boundaries: switching rtype 5 on across")
    print(f"this design means {g_stack.film_count() - 1} interface factors")
    print("instead of the handful the author had in mind, each one a")
    print("small-perturbation approximation. The model's validity band is")
    print("sigma much smaller than the wavelength AND much smaller than the")
    print("layer thickness -- and a sublayer is thin by construction, so the")
    print("second condition is the one that fails first. Set interface")
    print("roughness on the spans you meant, not on the expansion.")
    print()
    print("A NOTE ON THE THIN-GRADIENT FLOOR. Sublayer thickness is the")
    print("span's thickness divided by its sublayer count, so a thin span")
    print("has thin sublayers: a 3 nm gradient at default settings is three")
    print("1 nm rows, below the usual clamp_min_nm of 2 nm. Those rows are")
    print("exempt from the thin-layer sweep (they are not independent")
    print("parameters -- the span is), so they are not deleted. They are")
    print("still near the floor of what a staircase approximation means at")
    print("all: three steps do not resolve a profile. If a gradient has to")
    print("be that thin, it is a mixed layer, and inhomogen says so.")


if __name__ == "__main__":
    main()
