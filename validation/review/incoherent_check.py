#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""The incoherent cascade, checked against something that is not itself (C7).

Everything that pinned the incoherent path before this pinned it against a
copy of itself. `validation/parity/smatrix/refs/loom_matrix.py` is a port of
the same block sweep, citing the same paper. `core_engine.rs`'s
`intensity_path_matches_full_path_bitwise` compares the full walk against
the intensity-only walk -- two routes through one algorithm.
`test_physics_mirror.py` translates the Rust tests through the Python API:
same stacks, same conventions, same hand-written values. All useful, none
of them evidence that the physics is right.

Four checks that do not depend on Navette's implementation:

  A. **Closed form.** A lossless slab in air, both surfaces incoherent, has
     `R = 2*R1/(1 + R1)` and `T = (1 - R1)/(1 + R1)` -- the geometric sum
     over multiple internal reflections with phase discarded. One line of
     algebra, no transfer matrix anywhere.

  B. **Thickness independence.** The same slab must return bit-identical R
     across thicknesses spanning three decades. A coherent slab cannot do
     this (it fringes); an incoherent one must, because the only thing its
     thickness feeds is a phase that has been averaged away. Asserted on
     `to_bits()`, not a tolerance, because there is no numerical reason for
     a single bit to differ.

  C. **The phase average, on R and T.** This is the definition, not an
     agreement. A layer is incoherent when its round-trip phase is
     uniformly distributed, so the incoherent answer IS the coherent answer
     averaged over one phase period. Sweep the flagged layer's thickness
     across exactly one period with the flag OFF, average, and compare to
     the flagged answer. Nothing about the block sweep is assumed.

  D. **The phase average, on the Stokes vector.** The same identity applied
     to S0..S3 rather than R and T -- the check that would have caught C2.
     Averaging must be done on the Stokes COMPONENTS, which are bilinear in
     the fields; DOP and Delta are nonlinear functions of them, and
     averaging those instead measures something else (see the note in
     `part_d`). Mode B reproduces the averaged vector to 1e-16, including
     the partial depolarization it implies: every coherent sample here has
     DOP = 1 exactly, while their averaged vector has DOP = 0.994186, which
     is the physics -- incoherent superposition depolarizes. Mode A's S0
     and S1 pass the same average to 1e-16 while its S2 and S3 miss by
     2.3e-2 and 3.9e-3. That contrast IS C2, measured against the
     definition rather than against another implementation.

Absorbing layers need care in C/D, and round 1 recorded why: sweeping `d`
to turn the phase also sweeps `tau = exp(-2*Im(beta))`, so the naive test
disagrees at 6e-4 and is measuring itself. Holding `k*d` invariant while
the phase turns drops the residual to ~4e-7, the remainder being the
slab's own Fresnel coefficients moving as `k` is rescaled. The lossless
cases carry the exact claim; the absorbing one is run with the correction
and a looser bound, and is labelled as such.

Exit code 0 if every check passes.
"""
from __future__ import annotations

import sys
import warnings

import numpy as np

from navette.smatrix.smatrix import CoherenceMode, Request, ScatterMatrix

FAILURES: list[str] = []


def check(name, ok, detail=""):
    print(f"  {'OK  ' if ok else 'FAIL'}  {name}{'  ' + detail if detail else ''}")
    if not ok:
        FAILURES.append(name)


def _solve(idx, thick, req, *, flags=None, mode=CoherenceMode.FRONT_BLOCK,
           wls=None, angles=(0.0,)):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # C3/C5 thin-flag chatter is not the point here
        sm = ScatterMatrix(
            np.asarray(idx, dtype=complex), list(thick),
            wavelengths=np.asarray(wls if wls is not None else [550.0], float),
            angles=list(angles), incoherent_flags=flags, coherence_mode=mode,
        )
        return sm.compute(req, squeeze=False)


# ---------------------------------------------------------------------------
# A. closed form
# ---------------------------------------------------------------------------
def part_a():
    print("A. lossless slab in air against the closed-form geometric sum")
    n_s = 1.5
    r1 = ((1.0 - n_s) / (1.0 + n_s)) ** 2
    want_r = 2.0 * r1 / (1.0 + r1)
    want_t = (1.0 - r1) / (1.0 + r1)
    print(f"     R1 = {r1:.12f} -> R = {want_r:.12f}, T = {want_t:.12f}")

    for mode in (CoherenceMode.FRONT_BLOCK, CoherenceMode.COHERENCY_MATRIX):
        out = _solve([1.0, n_s, 1.0], [0.0, 1e5, 0.0],
                     Request.RS | Request.TS, flags=[0, 1, 0], mode=mode)
        got_r = float(np.asarray(out["Rs"]).ravel()[0])
        got_t = float(np.asarray(out["Ts"]).ravel()[0])
        check(f"{mode.name} Rs", abs(got_r - want_r) < 1e-12,
              f"{got_r:.12f} vs {want_r:.12f}")
        check(f"{mode.name} Ts", abs(got_t - want_t) < 1e-12,
              f"{got_t:.12f} vs {want_t:.12f}")

    # The coherent stack must NOT match -- otherwise the test is blind.
    out = _solve([1.0, n_s, 1.0], [0.0, 1e5, 0.0], Request.RS)
    coh = float(np.asarray(out["Rs"]).ravel()[0])
    check("a coherent slab does not satisfy it (control)",
          abs(coh - want_r) > 1e-6, f"coherent Rs = {coh:.12f}")


# ---------------------------------------------------------------------------
# B. thickness independence
# ---------------------------------------------------------------------------
def part_b():
    print("\nB. a lossless incoherent slab is bit-identical across thicknesses")
    bits = []
    for d in (10e3, 37e3, 150e3, 800e3, 3.7e6):  # 10 um -> 3.7 mm
        out = _solve([1.0, 1.5, 1.0], [0.0, d, 0.0],
                     Request.RS | Request.TS, flags=[0, 1, 0])
        bits.append((d, float(np.asarray(out["Rs"]).ravel()[0])))
    ref = bits[0][1]
    same = all(v.hex() == ref.hex() for _, v in bits)
    check("Rs bit-identical over 10 um -> 3.7 mm", same,
          f"{ref!r} " + ("" if same else str(bits)))

    # Control: with the flag cleared the same sweep must fringe.
    coh = []
    for d in (10e3, 10e3 + 91.0, 10e3 + 183.0):
        out = _solve([1.0, 1.5, 1.0], [0.0, d, 0.0], Request.RS)
        coh.append(float(np.asarray(out["Rs"]).ravel()[0]))
    check("a coherent slab fringes over the same span (control)",
          max(coh) - min(coh) > 1e-3, f"spread = {max(coh) - min(coh):.6f}")


# ---------------------------------------------------------------------------
# C / D. the phase average
# ---------------------------------------------------------------------------
def _period(n_slab, d0, lam, sin_th0, n_amb=1.0):
    """Thickness change that advances the round-trip phase by exactly 2*pi."""
    cos_t = np.sqrt(1.0 - (n_amb * sin_th0 / n_slab) ** 2 + 0j)
    return float(lam / (2.0 * np.real(n_slab * cos_t)))


def _averaged(idx, thick, req, slab, lam, angle_deg, n_samples=2048,
              mode=CoherenceMode.FRONT_BLOCK, hold_kd=False):
    """Coherent answer averaged over one round-trip phase period."""
    n_slab = complex(idx[slab])
    d0 = thick[slab]
    per = _period(n_slab.real, d0, lam, np.sin(np.radians(angle_deg)))
    acc = None
    for i in range(n_samples):
        d = d0 + per * i / n_samples
        th = list(thick)
        th[slab] = d
        ix = list(idx)
        if hold_kd and n_slab.imag != 0.0:
            # Sweeping d also sweeps tau; hold k*d so only the phase turns.
            ix[slab] = complex(n_slab.real, n_slab.imag * d0 / d)
        out = _solve(ix, th, req, mode=mode, wls=[lam], angles=(angle_deg,))
        vals = {k: np.asarray(v, float).ravel()[0] for k, v in out.items()}
        acc = vals if acc is None else {k: acc[k] + vals[k] for k in acc}
    return {k: v / n_samples for k, v in acc.items()}


def part_c():
    print("\nC. the flagged answer IS the coherent answer averaged over one period")
    # air / 95 nm H / 50 um slab / 110 nm L / air, slab flagged
    idx = [1.0, 2.35, 1.52, 1.46, 1.0]
    thick = [0.0, 95.0, 50_000.0, 110.0, 0.0]
    flags = [0, 0, 1, 0, 0]
    req = Request.RS | Request.TP
    for angle in (0.0, 45.0):
        inc = _solve(idx, thick, req, flags=flags, wls=[550.0], angles=(angle,))
        inc = {k: float(np.asarray(v, float).ravel()[0]) for k, v in inc.items()}
        avg = _averaged(idx, thick, req, 2, 550.0, angle)
        for key in ("Rs", "Tp"):
            d = abs(inc[key] - avg[key])
            check(f"{angle:g} deg {key}", d < 5e-9,
                  f"incoherent = {inc[key]:.9f}  averaged = {avg[key]:.9f}  "
                  f"diff = {d:.1e}")

    print("   absorbing slab, k*d held invariant (see module docstring)")
    idx_a = [1.0, 2.35, complex(1.52, 0.002), 1.46, 1.0]
    inc = _solve(idx_a, thick, req, flags=flags, wls=[550.0], angles=(0.0,))
    inc = {k: float(np.asarray(v, float).ravel()[0]) for k, v in inc.items()}
    avg = _averaged(idx_a, thick, req, 2, 550.0, 0.0, hold_kd=True)
    for key in ("Rs", "Tp"):
        d = abs(inc[key] - avg[key])
        check(f"absorbing {key}", d < 5e-6,
              f"incoherent = {inc[key]:.9f}  averaged = {avg[key]:.9f}  "
              f"diff = {d:.1e}")


def part_d():
    print("\nD. the same average on the Stokes vector -- the C2 check")
    idx = [1.0, 2.35, 1.52, 1.46, 1.0]
    thick = [0.0, 95.0, 50_000.0, 110.0, 0.0]
    flags = [0, 0, 1, 0, 0]
    # The STOKES COMPONENTS, not DOP. S0..S3 are bilinear in the fields, so
    # incoherent superposition averages them; DOP = sqrt(S1^2+S2^2+S3^2)/S0
    # and Delta = arg(cross) are nonlinear functions OF them, and averaging
    # those instead is a different quantity. The first draft of this check
    # averaged DOP and reported a 5.8e-3 "disagreement" that was entirely its
    # own: for this stack every coherent sample has DOP = 1 exactly (a
    # non-depolarizing stack), so their mean is 1, while the averaged Stokes
    # vector legitimately has DOP < 1 -- partial depolarization is precisely
    # what incoherent superposition produces. Pinned here because the mistake
    # is easy and the wrong version looks like an engine bug.
    req = Request.S0_R | Request.S1_R | Request.S2_R | Request.S3_R
    angle = 45.0

    avg = _averaged(idx, thick, req, 2, 550.0, angle)
    print("     averaged coherent Stokes: " + "  ".join(
        f"{k} = {avg[k]:.9f}" for k in ("S0_R", "S1_R", "S2_R", "S3_R")))

    b = _solve(idx, thick, req, flags=flags, mode=CoherenceMode.COHERENCY_MATRIX,
               wls=[550.0], angles=(angle,))
    b = {k: float(np.asarray(v, float).ravel()[0]) for k, v in b.items()}
    for key in ("S0_R", "S1_R", "S2_R", "S3_R"):
        d = abs(b[key] - avg[key])
        check(f"COHERENCY_MATRIX {key} matches the average", d < 5e-9,
              f"{b[key]:.9f} vs {avg[key]:.9f}  diff = {d:.1e}")

    dop_avg = np.hypot(np.hypot(avg["S1_R"], avg["S2_R"]), avg["S3_R"]) / avg["S0_R"]
    check("the averaged Stokes vector is partially depolarized",
          dop_avg < 1.0 - 1e-6,
          f"DOP of the averaged vector = {dop_avg:.9f} (each coherent sample "
          f"has DOP = 1)")

    # Mode A refuses to answer at all (C2), which is the fix for exactly the
    # failure this part would otherwise measure. Reaching past the door shows
    # what it is refusing.
    cross_req = Request.DOP_R
    try:
        _solve(idx, thick, cross_req, flags=flags,
               mode=CoherenceMode.FRONT_BLOCK, wls=[550.0], angles=(angle,))
        check("FRONT_BLOCK refuses a cross observable on a flagged stack",
              False, "it answered")
    except ValueError as exc:
        check("FRONT_BLOCK refuses a cross observable on a flagged stack",
              "front_block" in str(exc), str(exc)[:60])

    import navette._smatrix as native
    lam, n_layers = 550.0, len(idx)
    cache = []
    for n in idx:
        cache += [complex(n).real, complex(n).imag]
    raw = native.core_engine(
        np.array([lam]), np.array([np.sin(np.radians(angle))]), n_layers,
        np.array(cache), np.array(thick, float),
        np.array(flags, np.int32), np.zeros(n_layers, np.int32),
        np.zeros(n_layers), int(CoherenceMode.FRONT_BLOCK), int(req),
    )
    a = {k: float(np.asarray(v, float).ravel()[0]) for k, v in raw.items()}
    gaps = {k: abs(a[k] - avg[k]) for k in ("S2_R", "S3_R")}
    check("and what it refuses does NOT satisfy the average (C2, measured)",
          max(gaps.values()) > 1e-3,
          "  ".join(f"{k}: mode A = {a[k]:.9f} vs averaged {avg[k]:.9f} "
                    f"(gap {gaps[k]:.2e})" for k in gaps))
    for key in ("S0_R", "S1_R"):
        d = abs(a[key] - avg[key])
        check(f"  while its {key} passes the same average", d < 5e-9,
              f"{a[key]:.9f} vs {avg[key]:.9f}  diff = {d:.1e}")


def main():
    part_a()
    part_b()
    part_c()
    part_d()
    print("\nALL OK" if not FAILURES else f"\nFAILURES: {FAILURES}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
