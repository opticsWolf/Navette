# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compute every figure on the project site with Navette itself.

Each figure is a real engine run, not an illustration: the site's charts are
drawn from the JSON this script writes, so a change in the physics shows up
as a change on the page the next time the script runs.

Run:  python site/scripts/build_data.py      (writes site/data/figures.json)
"""

import json
import math
from pathlib import Path

import numpy as np

import navette
from navette.materials import MaterialSpec, evaluate
from navette.smatrix.smatrix import CoherenceMode, Request, ScatterMatrix
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis.pipeline import PipelineConfig, run_needle

OUT = Path(__file__).resolve().parents[1] / "data" / "figures.json"

H, L = 2.35, 1.46                       # TiO2-ish / SiO2-ish, as in examples/
BK7 = MaterialSpec("Sellmeier", dict(B1=1.03961212, C1=0.00600069867,
                                     B2=0.231792344, C2=0.0200179144,
                                     B3=1.01046945, C3=103.560653))
GOLD = MaterialSpec("Drude", dict(omega_p=9.0, gamma=0.07, epsilon_inf=9.5))
# Amorphous silicon, Tauc-Lorentz (Jellison & Modine 1996 parameters).
A_SI = MaterialSpec("TaucLorentz", dict(Eg=1.20, osc=[(122.0, 3.45, 2.54)],
                                        epsilon_inf=1.15))
# Crystalline silicon, rounded from Green (2008); k from alpha = 4*pi*k/lambda.
C_SI_TABLE = np.array([
    # nm     n      k
    [300, 5.00, 4.18], [350, 5.48, 2.94], [370, 6.70, 1.30], [400, 5.57, 0.387],
    [450, 4.67, 0.135], [500, 4.30, 0.073], [550, 4.08, 0.041],
    [600, 3.94, 0.025], [650, 3.85, 0.016], [700, 3.78, 0.011],
    [800, 3.68, 5.0e-3], [900, 3.61, 2.2e-3], [1000, 3.57, 5.1e-4],
    [1100, 3.54, 3.1e-5], [1200, 3.52, 2.0e-7],
])


def r(x, nd=6):
    """Round for JSON: the page draws, it does not re-derive."""
    a = np.asarray(x, dtype=float)
    out = np.where(np.isfinite(a), np.round(a, nd), np.nan)
    return [None if math.isnan(v) else float(v) for v in out.ravel()]


def solve(n_rows, d, wl, angles, **kw):
    """``n_rows``: one entry per row, each a scalar or an nk array on ``wl``."""
    wl = np.atleast_1d(np.asarray(wl, dtype=float))
    n = np.vstack([np.broadcast_to(np.asarray(v, dtype=complex), wl.shape)
                   for v in n_rows])
    return ScatterMatrix(n, list(d), wavelengths=list(wl),
                         angles=list(np.atleast_1d(angles)), **kw)


def qw_mirror(pairs, lam0):
    rows = [H, L] * pairs
    ds = [lam0 / 4 / H, lam0 / 4 / L] * pairs
    return rows, ds


# --- 1. a quarter-wave mirror over wavelength and angle ---------------------
def bragg():
    wl = np.linspace(380.0, 780.0, 201)
    ang = np.linspace(0.0, 85.0, 35)
    rows, ds = qw_mirror(8, 560.0)
    sub = evaluate(BK7, wl)
    sm = solve([1.0] + rows + [sub], [0.0] + ds + [0.0], wl, ang)
    out = sm.compute(Request.RS | Request.RP, squeeze=False)
    rs, rp = np.asarray(out["Rs"]), np.asarray(out["Rp"])
    return dict(wl=r(wl, 2), angles=r(ang, 2), Rs=[r(x, 4) for x in rs],
                Rp=[r(x, 4) for x in rp], pairs=8, lam0=560.0)


# --- 2. stability: a silicon wafer as a finite layer -------------------------
def tmm_naive(n_rows, ds, wl):
    """The textbook characteristic matrix (Abeles), normal incidence."""
    M = np.eye(2, dtype=complex)
    for n, d in zip(n_rows[1:-1], ds[1:-1]):
        delta = 2 * np.pi * n * d / wl
        M = M @ np.array([[np.cos(delta), -1j * np.sin(delta) / n],
                          [-1j * n * np.sin(delta), np.cos(delta)]])
    n0, ns = n_rows[0], n_rows[-1]
    B = M[0, 0] + M[0, 1] * ns
    C = M[1, 0] + M[1, 1] * ns
    return abs((n0 * B - C) / (n0 * B + C)) ** 2


def stability():
    wl = np.linspace(300.0, 1200.0, 451)
    t = C_SI_TABLE
    si = np.interp(wl, t[:, 0], t[:, 1]) + 1j * np.exp(
        np.interp(wl, t[:, 0], np.log(t[:, 2])))
    d_wafer = 525_000.0  # nm
    with np.errstate(all="ignore"):
        naive = np.array([tmm_naive([1.0, s, 1.0], [0, d_wafer, 0], w)
                          for s, w in zip(si, wl)])
    sm = solve([1.0, si, 1.0], [0.0, d_wafer, 0.0], wl, [0.0],
               incoherent_flags=[0, 1, 0])
    rt = sm.compute(Request.RS | Request.TS)
    fresnel = np.abs((1 - si) / (1 + si)) ** 2
    first_ok = float(wl[np.argmax(np.isfinite(naive))])
    return dict(wl=r(wl, 2), naive=r(naive, 5), R=r(rt["Rs"], 5),
                T=r(rt["Ts"], 5), fresnel=r(fresnel, 5),
                nan_fraction=float(np.mean(~np.isfinite(naive))),
                first_finite_nm=first_ok, d_um=d_wafer / 1000)


# --- 3. partial coherence: a 1 mm slide under a TiO2 film --------------------
def coherence():
    wl = np.linspace(600.0, 601.0, 1201)
    glass = 1.515
    rows = [1.0, H, glass, 1.0]
    ds = [0.0, 80.0, 1_000_000.0, 0.0]
    coh = solve(rows, ds, wl, [0.0]).compute(Request.RS)["Rs"]
    inc = solve(rows, ds, wl, [0.0], incoherent_flags=[0, 0, 1, 0]).compute(
        Request.RS)["Rs"]
    return dict(wl=r(wl, 4), coherent=r(coh, 5), incoherent=r(inc, 5),
                period_nm=float(600.0 ** 2 / (2 * glass * 1e6)))


# --- 4. graded versus stepped, from examples/rugate_gradient_vs_discrete.py --
def rugate():
    import importlib.util
    path = Path(__file__).resolve().parents[2] / "examples" / \
        "rugate_gradient_vs_discrete.py"
    spec = importlib.util.spec_from_file_location("rugate_example", path)
    ex = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ex)
    wl_dense = np.linspace(450.0, 780.0, 331)
    ex.WL = wl_dense  # the example's builders read the module grid
    g = ex.build(*ex.graded())
    s = ex.build(*ex.stepped())
    return dict(wl=r(wl_dense, 2), graded=r(ex.spectrum(g), 5),
                stepped=r(ex.spectrum(s), 5), rows_graded=g.film_count(),
                rows_stepped=s.film_count(), lam0=ex.LAM0)


# --- 5. needle synthesis of a broadband AR coating ---------------------------
def needle():
    band = (400.0, 700.0)
    wl = np.linspace(*band, 61)
    glass = 1.52 + 0j
    tol = 0.002
    tc = TargetCollection()
    tc.add(SpectralTarget(wl, np.zeros_like(wl), np.full_like(wl, tol),
                          0.0, "s", "R"))
    hi, lo = np.full(wl.shape, H + 0j), np.full(wl.shape, L + 0j)
    # A deliberately poor four-layer seed; each film may host the other one.
    seed = [(lo, 250.0), (hi, 150.0), (lo, 250.0), (hi, 150.0)]
    cfg = PipelineConfig(max_macro_cycles=30, needles_per_cycle=2)
    res = run_needle(layers=seed, targets=tc, angles_deg=[0.0],
                     wavelengths=wl, contrast={0: hi, 1: lo, 2: hi, 3: lo},
                     pipeline_config=cfg,
                     substrate=(np.full(wl.shape, glass), "sub"))
    films = res["stack"].to_dict()["films"]
    # Films are listed from the ambient side, as stack_from_layers takes them.
    ns = [complex(np.asarray(f["nk"]).ravel()[0]).real for f in films]
    ds = [float(f["thickness"]) for f in films]
    rows = [1.0] + [complex(n) for n in ns] + [glass]
    # Re-score the finished design on the run grid: if the film order were
    # read the wrong way round this would not reproduce the run's merit.
    R_run = np.asarray(solve(rows, [0.0] + ds + [0.0], wl, [0.0]).compute(
        Request.RS)["Rs"])
    check = float(np.sum((R_run / tol) ** 2))
    assert abs(check - res["final_mf"]) <= 1e-6 * max(1.0, res["final_mf"]), \
        (check, res["final_mf"])
    wl_dense = np.linspace(380.0, 780.0, 201)
    R = solve(rows, [0.0] + ds + [0.0], wl_dense, [0.0]).compute(
        Request.RS)["Rs"]
    R_bare = solve([1.0, glass], [0.0, 0.0], wl_dense, [0.0]).compute(
        Request.RS)["Rs"]
    R_seed = solve([1.0, L, H, L, H, glass],
                   [0.0, 250.0, 150.0, 250.0, 150.0, 0.0],
                   wl_dense, [0.0]).compute(Request.RS)["Rs"]
    phases = res["phases"]
    first = phases[0]["needle_results"][0]["merit_before"]
    history = [dict(cycle=0, merit=first, layers=len(seed),
                    thickness=800.0)]
    history += [dict(cycle=int(p["macro_cycle"]), merit=float(p["mf_end"]),
                     layers=int(p["layer_count"]),
                     thickness=round(float(p["total_thickness_nm"]), 2))
                for p in phases]
    return dict(wl=r(wl_dense, 2), R=r(R, 6), R_bare=r(R_bare, 6),
                R_seed=r(R_seed, 6),
                layers=[dict(n=round(n, 3), d=round(d, 2)) for n, d in zip(ns, ds)],
                band=list(band), tol=tol, termination=str(res["termination"]),
                final_mf=float(res["final_mf"]), history=history,
                rms_in_band=float(np.sqrt(np.mean(R_run ** 2))),
                max_in_band=float(np.max(R_run)))


# --- 6. surface plasmon resonance (Kretschmann) ------------------------------
def spr():
    wl = 633.0
    ang = np.linspace(38.0, 52.0, 281)
    n_prism = complex(evaluate(BK7, np.array([wl]))[0])
    gold = complex(evaluate(GOLD, np.array([wl]))[0])
    curves = {}
    for d in (30, 40, 50, 60):
        out = solve([n_prism, gold, 1.0], [0.0, float(d), 0.0], wl, ang).compute(
            Request.RP)
        curves[str(d)] = r(out["Rp"], 5)
    return dict(angles=r(ang, 3), curves=curves, n_prism=n_prism.real,
                gold=[gold.real, gold.imag], wl=wl)


# --- 7. ellipsometric trajectory of a growing oxide --------------------------
def ellipsometry():
    wl, aoi = 633.0, 70.0
    si = 3.875 + 0.018j
    thick = np.arange(0.0, 601.0, 2.0)
    psi, delta = [], []
    for d in thick:
        e = solve([1.0, 1.457, si], [0.0, float(d), 0.0], wl, [aoi]).ellipsometry()
        psi.append(float(np.ravel(e["Psi_R"])[0]))
        delta.append(float(np.ravel(e["Delta_R"])[0]))
    # Delta comes back in (-pi, pi]; the classic trajectory plot uses 0..360.
    psi, delta = np.degrees(psi), np.mod(np.degrees(delta), 360.0)
    return dict(d=r(thick, 1), psi=r(psi, 4), delta=r(delta, 4), aoi=aoi, wl=wl)


# --- 8. group delay and GDD of the mirror ------------------------------------
def group_delay():
    wl = np.linspace(470.0, 690.0, 441)
    rows, ds = qw_mirror(8, 560.0)
    sm = solve([1.0] + rows + [1.52], [0.0] + ds + [0.0], wl, [0.0])
    out = sm.compute(Request.DISP_R_S | Request.RS)
    return dict(wl=r(wl, 2), GD=r(out["GD_R_s"], 4), GDD=r(out["GDD_R_s"], 3),
                R=r(out["Rs"], 5), lam0=560.0)


# --- 9. dispersion models ----------------------------------------------------
def dispersion():
    wl = np.linspace(300.0, 1000.0, 141)
    mats = {"BK7 (Sellmeier)": BK7, "Gold (Drude)": GOLD,
            "a-Si (Tauc-Lorentz)": A_SI}
    return dict(wl=r(wl, 2), materials={
        name: dict(n=r(nk.real, 4), k=r(nk.imag, 4))
        for name, spec in mats.items()
        for nk in [evaluate(spec, wl)]})


def main():
    figs = {}
    for name, fn in [("bragg", bragg), ("stability", stability),
                     ("coherence", coherence), ("rugate", rugate),
                     ("needle", needle), ("spr", spr),
                     ("ellipsometry", ellipsometry), ("group_delay", group_delay),
                     ("dispersion", dispersion)]:
        figs[name] = fn()
        print(f"  {name:<14} ok")
    figs["meta"] = dict(navette=navette.__version__)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(figs, separators=(",", ":"), allow_nan=False),
                   encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size / 1024:.0f} KiB)")


if __name__ == "__main__":
    main()
