# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""A rugate notch filter, built as gradient spans, against a stepped stack.

A rugate is an index modulation rather than a stack of layers: the index
sweeps up and down through the coating instead of stepping. Its selling
point is what it does NOT have. No abrupt interfaces means no higher-order
stopbands, so a notch at 1550 nm leaves the visible alone; and with the
modulation faded in and out (apodization) and the ends graded onto the
substrate and the top material (index matching), the passband has no
sidelobe ripple either.

The coating here is one a sputter coater could make: SiO2 and Ta2O5
co-deposited, the mixture fraction varied through the film, both with
dispersion, on BK7. The mixture index is Bruggeman's, the same rule the
gradient spans use.

How the profile is built. Navette's gradient spans are linear in mixture
fraction (`FixedSpan`), so a smooth profile is sampled at
`SPANS_PER_PERIOD` points per period and joined by linear spans. Twelve
per period is converged: doubling it moves no reflectance below by more
than 0.01, nor the OD by more than 0.1. The profile itself is the
textbook one: ln n is sinusoidal in OPTICAL thickness, because the
reflection of a slowly varying medium is driven by d(ln n)/dx. The
obvious alternative, n sinusoidal in physical depth, is not
harmonic-free: at this contrast it reflects over half the light at
lambda0/2.
The fraction for each sampled index is solved at the notch wavelength;
everywhere else the mixture disperses as its two materials do.

Three designs, the same index range (n = 1.55 .. 1.95) and period count:

  * **stepped.** Quarter-wave pairs of the two extreme mixtures.
  * **sinusoid.** The bare modulation: abrupt ends, constant amplitude.
  * **rugate.** The same modulation, apodized and index-matched.

Run:  python examples/rugate_gradient_vs_discrete.py
"""

import time

import numpy as np

from navette.materials import MaterialSpec, evaluate
from navette.smatrix.smatrix import Request, ScatterMatrix
from navette.synthesis.pipeline import stack_from_layers

# --- the problem ----------------------------------------------------------
WL = np.linspace(450.0, 1700.0, 501)

LAM0 = 1550.0                   # notch centre
N_HI, N_LO = 1.95, 1.55         # the index range the modulation spans
PERIODS = 40
SPANS_PER_PERIOD = 12
APODIZE = 0.35                  # share of the modulated core faded at each end
MATCH_PERIODS = 1.0             # optical length of each index-matching ramp

#: Fused silica (Malitson) and a representative sputtered Ta2O5; put your
#: own coater's fits here. BK7 is the substrate.
SIO2 = MaterialSpec("Sellmeier", dict(B1=0.6961663, C1=0.0684043 ** 2,
                                      B2=0.4079426, C2=0.1162414 ** 2,
                                      B3=0.8974794, C3=9.896161 ** 2))
TA2O5 = MaterialSpec("Cauchy", dict(A=2.05, B=0.016, C=0.0003))
BK7 = MaterialSpec("Sellmeier", dict(B1=1.03961212, C1=0.00600069867,
                                     B2=0.231792344, C2=0.0200179144,
                                     B3=1.01046945, C3=103.560653))

#: The optical period: lambda0 / 2, whatever the index.
PERIOD_OPT = LAM0 / 2.0


def n_at(spec, lam):
    return float(evaluate(spec, np.array([lam]))[0].real)


N_SIO2, N_TA2O5, N_SUB = (n_at(m, LAM0) for m in (SIO2, TA2O5, BK7))


def fraction(n):
    """Ta2O5 volume fraction in SiO2 whose Bruggeman index at LAM0 is `n`.

    Two-phase Bruggeman is linear in f once the target permittivity is
    fixed, so this is closed-form, not a root search.
    """
    eps = n * n
    a = (N_SIO2 ** 2 - eps) / (N_SIO2 ** 2 + 2 * eps)
    b = (N_TA2O5 ** 2 - eps) / (N_TA2O5 ** 2 + 2 * eps)
    return float(np.clip(a / (a - b), 0.0, 1.0))


def mixture(n):
    """A homogeneous co-deposited film of index `n` at LAM0."""
    return MaterialSpec("Bruggeman", dict(host=SIO2, inclusion=TA2O5,
                                          fraction=fraction(n)))


def smoothstep(u):
    """Quintic 0 -> 1 with zero slope and curvature at both ends."""
    u = np.clip(u, 0.0, 1.0)
    return u ** 3 * (10.0 - 15.0 * u + 6.0 * u ** 2)


def index_profile(shaped):
    """n(x) over optical depth x from the ambient side, and the total x.

    ln n = ln n0 + A w(x) sin(2 pi x / period), n0 the geometric mean of
    the range so the extremes are exactly N_LO and N_HI. `shaped` adds
    the envelope w (quintic fades over APODIZE of the core at each end)
    and the quintic ramps in ln n from the top material and onto the
    substrate; without it w = 1 and the modulation starts and stops dead.
    """
    ln0, amp = 0.5 * np.log(N_HI * N_LO), 0.5 * np.log(N_HI / N_LO)
    core = PERIODS * PERIOD_OPT
    ramp = MATCH_PERIODS * PERIOD_OPT if shaped else 0.0

    def n(x):
        if x < ramp:
            return np.exp(np.log(N_SIO2) + (ln0 - np.log(N_SIO2)) * smoothstep(x / ramp))
        if x > ramp + core:
            u = (x - ramp - core) / ramp
            return np.exp(ln0 + (np.log(N_SUB) - ln0) * smoothstep(u))
        u = (x - ramp) / core
        w = smoothstep(u / APODIZE) * smoothstep((1 - u) / APODIZE) if shaped else 1.0
        return np.exp(ln0 + amp * w * np.sin(2 * np.pi * (x - ramp) / PERIOD_OPT))

    return n, core + 2 * ramp


def graded(shaped):
    """The profile as linear `FixedSpan` spans between sampled points.

    Each span is a SiO2 host with Ta2O5 inclusions, its fraction running
    from one sample to the next, so consecutive spans meet at the same
    fraction and the profile is continuous across every join. A
    discontinuity there would put back exactly the interface the
    modulation exists to remove. The spans are pinned background
    (`optimize=False, needle=False`), which is what makes them expand
    with the profile instead of homogenizing.

    A span's physical thickness is its optical length over the index,
    integrated by the trapezoid rule across the span.
    """
    n_of, total = index_profile(shaped)
    k = int(round(total / PERIOD_OPT * SPANS_PER_PERIOD))
    x = np.linspace(0.0, total, k + 1)
    n = np.array([n_of(v) for v in x])
    d = np.diff(x) * 0.5 * (1.0 / n[1:] + 1.0 / n[:-1])
    f = [fraction(v) for v in n]
    host, inclusion = evaluate(SIO2, WL), evaluate(TA2O5, WL)
    layers, names, flags = [], [], {}
    for i in range(k):
        name = f"g{i}"
        layers.append((host, d[i]))
        names.append(name)
        flags[name] = {"optimize": False, "needle": False, "gradient": {
            "material_b": inclusion,
            "f_start": f[i],
            "f_end": f[i + 1],
            "ema": "Bruggeman",
        }}
    return layers, names, flags


def stepped():
    """Quarter-wave pairs of the two extreme mixtures, same optical period."""
    hi, lo = mixture(N_HI), mixture(N_LO)
    layers, names = [], []
    for i in range(PERIODS):
        layers.append((hi, LAM0 / (4.0 * N_HI)))
        names.append(f"h{i}")
        layers.append((lo, LAM0 / (4.0 * N_LO)))
        names.append(f"l{i}")
    return layers, names, {}


def build(layers, names, flags):
    return stack_from_layers(layers, WL, {}, names=names,
                             per_film_flags=flags,
                             substrate=(evaluate(BK7, WL), "sub"))[0]


def solver(stack):
    """A solver over the stack's EXPANDED rows.

    `to_dict()` returns the expansion -- one entry per row, each with the
    nk the solver will use -- so this is the coating the engine builds,
    not the authored spans.
    """
    rows = stack.to_dict()["films"]
    n = np.vstack([np.ones_like(WL)]
                  + [np.asarray(r["nk"]) for r in rows]
                  + [evaluate(BK7, WL)])
    d = np.array([0.0] + [float(r["thickness"]) for r in rows] + [0.0])
    return ScatterMatrix(n, d, wavelengths=WL, angles=[0.0])


def spectrum(stack):
    """Rs at normal incidence on the run grid."""
    return np.asarray(solver(stack).compute(Request.RS)["Rs"]).ravel()


def figures(R):
    """What a notch filter is specified by, read off one spectrum.

    od      optical density at the notch: -log10 of the worst T within
            +/-15 nm of LAM0 (lossless, so T = 1 - R)
    fwhm    width of the R > 0.5 band around LAM0 [nm]
    worst   the highest R anywhere in the passband, 450..1400 nm
    third   the highest R near LAM0 / 3, where a stepped stack's
            third-order band sits (dispersion moves it to ~525 nm)
    """
    stop = np.abs(WL - LAM0) <= 15.0
    od = -np.log10(max(1.0 - float(R[stop].min()), 1e-12))
    i = int(np.argmin(np.abs(WL - LAM0)))
    lo = hi = i
    while lo > 0 and R[lo - 1] > 0.5:
        lo -= 1
    while hi < len(WL) - 1 and R[hi + 1] > 0.5:
        hi += 1
    band = WL <= 1400.0
    third = (WL >= 490.0) & (WL <= 560.0)
    return dict(od=od, fwhm=float(WL[hi] - WL[lo]), worst=float(R[band].max()),
                third=float(R[third].max()))


def solve_ms(stack, repeats=3):
    """Best-of-`repeats` wall time of one solve over the grid [ms]."""
    sm = solver(stack)
    best = float("inf")
    for _ in range(repeats):
        t0 = time.perf_counter()
        sm.compute(Request.RS)
        best = min(best, time.perf_counter() - t0)
    return 1e3 * best


def designs():
    """The three designs, built: (label, stack, spans authored)."""
    out = []
    for label, shaped in (("sinusoid", False), ("rugate", True)):
        g = graded(shaped)
        out.append((label, build(*g), len(g[1])))
    s = stepped()
    out.insert(0, ("stepped", build(*s), len(s[1])))
    return out


def main() -> None:
    print(f"notch at {LAM0:.0f} nm, {PERIODS} periods, n = {N_LO} .. {N_HI} "
          f"(SiO2 {N_SIO2:.3f} : Ta2O5 {N_TA2O5:.3f} at {LAM0:.0f} nm)")
    print(f"fractions {fraction(N_LO):.3f} .. {fraction(N_HI):.3f} Ta2O5; "
          f"{SPANS_PER_PERIOD} linear spans per period")
    print()
    print(f"  {'design':<10}{'spans':>7}{'rows':>7}{'um':>7}{'solve ms':>10}"
          f"{'OD':>7}{'FWHM':>7}{'worst R':>9}{'3rd R':>8}")
    rows = {}
    for label, stack, spans in designs():
        m = figures(spectrum(stack))
        rows[label] = (stack, spans, m)
        print(f"  {label:<10}{spans:>7}{stack.film_count():>7}"
              f"{stack.total_thickness() / 1e3:>7.1f}{solve_ms(stack):>10.1f}"
              f"{m['od']:>7.2f}{m['fwhm']:>7.0f}{m['worst']:>9.3f}"
              f"{m['third']:>8.3f}")
    print()

    st, sn, rg = rows["stepped"][2], rows["sinusoid"][2], rows["rugate"][2]
    g_stack, g_spans, _ = rows["rugate"]
    print("SHAPE. The stepped stack reflects at LAM0 / 3 as well as at LAM0:")
    print(f"R = {st['third']:.3f} near 525 nm, a green mirror nobody asked for. A")
    print("square wave has odd harmonics and a sinusoid has none, so neither")
    print("graded design has a band there. The rugate sits at the bare-surface")
    print(f"level ({rg['third']:.3f}); the sinusoid's {sn['third']:.3f} is the ripple of its")
    print("abrupt ends, the same ripple it shows across the whole passband.")
    print("The stepped stack's fundamental is stronger and wider -- its first")
    print("Fourier component is 4/pi of the sinusoid's at the same index range.")
    print()
    print("SIDELOBES. The bare sinusoid still rings: its abrupt ends are")
    print(f"interfaces, and the passband peaks at R = {sn['worst']:.3f}. Apodization")
    print("and the matching ramps take that to "
          f"{rg['worst']:.3f}, at the price of rejection")
    print(f"(OD {sn['od']:.1f} -> {rg['od']:.1f}): the faded periods reflect less.")
    print("What is left in the passband is the outer surface: no index in")
    print("the mix is lower than SiO2's, so the top cannot be matched to air.")
    print()
    print(f"ROWS. Authored as {g_spans} spans, solved as "
          f"{g_stack.film_count()}: each span expands")
    print("into sublayers at assembly, and the solver sees rows. Authoring a")
    print("graded design is cheap and solving one is not; the row count is")
    print("the number that shows up in the solve time above.")
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
