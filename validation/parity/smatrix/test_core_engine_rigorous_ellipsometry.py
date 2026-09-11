#!/usr/bin/env python3
"""Comparison test for core_engine_rigorous_ellipsometry — numba vs rust (R2.4a).

The legacy kernel returned a fixed 13-tuple and took a ``debug_flag`` to
decide whether the last element was filled. Its successor is the
request-driven ``navette._smatrix.core_engine``: twelve of those thirteen are
channels in a mask, and the thirteenth — the energy-conservation diagnostic —
became a separate helper (``solver_energy_conservation``, R1.2). The whole
tuple is still compared here; nothing is dropped because it moved.

**The mode mapping is the one thing this port could get wrong**, so it is
stated rather than assumed. The legacy engine's own docstring describes its
partial-coherence treatment as: reflection S₂/S₃ from the *first coherent
block's* field amplitudes, transmission S₂/S₃ from a Mueller cross-term
product accumulated across blocks. In `core_engine.rs` that is exactly the
``track_cross_channel == false`` path — ``cross_r = rp₀·conj(rs₀)`` and
``cross_t = cross_t_acc`` — which is ``FRONT_BLOCK`` (mode A).
``COHERENCY_MATRIX`` (mode B) is a *different, later* treatment that cascades
the complex p-s channel through the incoherent echoes, and it is not what
this kernel did. Half the cases below flag an incoherent layer so the
distinction is actually exercised: with every layer coherent, all three modes
collapse to the same answer and the mapping would be untested.

Delta is compared on the circle, and the last case is built to sit on top of
±π (transparent, low contrast, near normal incidence) so that "on the circle"
is not a hypothetical. Stated plainly: **the two engines pick the same side of
the cut**, there and everywhere else, so a plain difference would pass this
file today. The circle comparison is a guard, not a finding -- it is here so
that if the two ever do disagree about the branch, the report is the 1e-16
disagreement it really is rather than a spurious 2π.
"""

import sys, os, time
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
UNIT_NAME = "core_engine_rigorous_ellipsometry"

NUM_WAVS = 50
NUM_ANGLES = 10
N_LAYERS = 6

# Dual-shaped: runs as a script and is collected by pytest (R2.4).
sys.path.insert(0, PROJECT_ROOT)
from _parity import console_utf8, report, require_reference  # noqa: E402

console_utf8()

numba_imported = False
try:
    sys.path.insert(0, os.path.join(SCRIPT_DIR, "refs"))
    from loom_matrix import core_engine_rigorous_ellipsometry as numba_func
    numba_imported = True
except ImportError:
    numba_func = None

require_reference(numba_imported, "numba reference loom_matrix",
                  "pip install numba; see validation/parity/smatrix/refs/")

# Unified layout: kernels live in navette._smatrix (aggregated extension).
import navette._smatrix as rust_mod  # noqa: E402
from navette.smatrix.smatrix import CoherenceMode, Request  # noqa: E402

rust_func = rust_mod.core_engine
energy_conservation = rust_mod.solver_energy_conservation

REQUESTED = Request.PHOTOMETRY | Request.ELLIPSOMETRY

# The legacy return tuple, in order. `None` is the conservation diagnostic:
# not a channel any more, reconstructed below from the four intensities.
LEGACY_TUPLE = [
    "Psi_R", "Delta_R", "DOP_R", "Rs", "Rp", "R_avg",
    "Psi_T", "Delta_T", "DOP_T", "Ts", "Tp", "T_avg",
    None,
]
ANGLES = {"Delta_R", "Delta_T"}


def angular_diff(a, b):
    """|a − b| on the circle, so ±π does not read as a 2π disagreement."""
    return np.abs(np.angle(np.exp(1j * (np.asarray(a) - np.asarray(b)))))


def generate_test_data(rng, incoherent_layer=None, rough=False,
                       angles=(5.0, 70.0), lossless=False):
    """One stack, in both engines' input layouts (see the photometry port)."""
    wavls = np.linspace(400.0, 800.0, NUM_WAVS).astype(np.float64)
    # Away from exactly 0°, where Psi/Delta are degenerate by construction
    # (r_p and r_s coincide) and the comparison would say nothing.
    angles_deg = np.linspace(angles[0], angles[1], NUM_ANGLES).astype(np.float64)
    sin_theta_arr = np.sin(np.radians(angles_deg)).astype(np.float64)

    if lossless:
        # Transparent, low contrast, near normal incidence: Delta_R sits on
        # top of ±π and individual grid points land on either side of the
        # branch cut. This is the case that makes the angular comparison
        # load-bearing -- without it, comparing Delta with a plain difference
        # would pass just as well.
        n_real_base = [1.0, 1.46, 1.38, 1.46, 1.38, 1.52]
        n_imag_base = [0.0] * N_LAYERS
    else:
        n_real_base = [1.0, 1.5 + rng.random() * 0.3, 2.0 + rng.random() * 0.5,
                       1.8 + rng.random() * 0.4, 2.5 + rng.random() * 0.3, 1.52]
        n_imag_base = [0.0, 0.001 + rng.random() * 0.01, 0.01 + rng.random() * 0.05,
                       0.001 + rng.random() * 0.02, 0.02 + rng.random() * 0.08, 0.0]

    n_complex = np.empty((NUM_WAVS, N_LAYERS), dtype=np.complex128)
    for wv in range(NUM_WAVS):
        wav_factor = (wavls[wv] - 400.0) / 400.0
        for li in range(N_LAYERS):
            jitter = 0.0 if lossless else wav_factor * rng.uniform(-0.1, 0.1)
            n_complex[wv, li] = complex(n_real_base[li] + jitter,
                                        max(0.0, n_imag_base[li]))
    n_flat = np.empty(NUM_WAVS * N_LAYERS * 2, dtype=np.float64)
    n_flat[0::2] = n_complex.real.ravel()
    n_flat[1::2] = n_complex.imag.ravel()

    thicknesses = np.array([50.0, 100.0, 200.0, 150.0, 80.0, 500.0],
                           dtype=np.float64)
    incoherent_flags = np.zeros(N_LAYERS, dtype=np.int32)
    if incoherent_layer is not None:
        incoherent_flags[incoherent_layer] = 1
    rough_types = np.zeros(N_LAYERS, dtype=np.int32)
    rough_vals = np.zeros(N_LAYERS, dtype=np.float64)
    if rough:
        rough_types[1:N_LAYERS - 1] = [1, 2, 5, 4][:N_LAYERS - 2]
        rough_vals[1:N_LAYERS - 1] = rng.uniform(1.0, 4.0, N_LAYERS - 2)

    return dict(wavls=wavls, sin_theta_arr=sin_theta_arr,
                n_complex=n_complex, n_flat=n_flat, thicknesses=thicknesses,
                incoherent_flags=incoherent_flags, rough_types=rough_types,
                rough_vals=rough_vals)


def call_numba(case, debug_flag=1):
    return numba_func(
        case["wavls"], case["sin_theta_arr"], N_LAYERS, case["n_complex"],
        case["thicknesses"], case["incoherent_flags"], case["rough_types"],
        case["rough_vals"], np.int32(debug_flag),
    )


def call_rust(case):
    return rust_func(
        case["wavls"], case["sin_theta_arr"], N_LAYERS, case["n_flat"],
        case["thicknesses"], case["incoherent_flags"], case["rough_types"],
        case["rough_vals"], int(CoherenceMode.FRONT_BLOCK), int(REQUESTED),
    )


rng = np.random.default_rng(seed=1234)
CASES = [
    ("all coherent",        generate_test_data(rng)),
    ("all coherent, rough", generate_test_data(rng, rough=True)),
    ("incoherent layer 4",  generate_test_data(rng, incoherent_layer=4)),
    ("incoherent + rough",  generate_test_data(rng, incoherent_layer=2, rough=True)),
    ("Delta on the branch cut",
     generate_test_data(rng, angles=(1.0, 15.0), lossless=True)),
]

print(f"=== CORRECTNESS TEST: {UNIT_NAME} ===")
all_pass = True
for label, case in CASES:
    numba_out = call_numba(case)
    out = call_rust(case)
    # The 13th element: a separate helper since R1.2, not a channel.
    got_all = dict(out)
    got_all[None] = energy_conservation(
        out["Rs"], out["Rp"], out["Ts"], out["Tp"])

    for key, ref in zip(LEGACY_TUPLE, numba_out):
        ref = np.asarray(ref)
        got = np.asarray(got_all[key])
        name = key if key is not None else "conservation_err"
        if got.shape != ref.shape:
            print(f"  {label} {name}: FAIL | shape {got.shape} vs {ref.shape}")
            all_pass = False
            continue
        if key in ANGLES:
            diff = float(np.max(angular_diff(ref, got)))
            ok = diff < 1e-9
        else:
            diff = float(np.max(np.abs(ref - got)))
            ok = np.allclose(ref, got, rtol=1e-9, atol=1e-12)
        print(f"  {label} {name}: {'PASS' if ok else 'FAIL'} | diff_max={diff:.2e}")
        all_pass &= ok

# The mode mapping, asserted rather than assumed: on a stack with an
# incoherent layer, FRONT_BLOCK must reproduce the legacy kernel (it does,
# above) and COHERENCY_MATRIX must NOT -- otherwise "mode A is the legacy
# treatment" is an untested claim that happens to hold because the two modes
# never differ here.
inc_case = CASES[2][1]
mode_b = rust_func(
    inc_case["wavls"], inc_case["sin_theta_arr"], N_LAYERS, inc_case["n_flat"],
    inc_case["thicknesses"], inc_case["incoherent_flags"],
    inc_case["rough_types"], inc_case["rough_vals"],
    int(CoherenceMode.COHERENCY_MATRIX), int(REQUESTED),
)
ref_delta = np.asarray(call_numba(inc_case)[LEGACY_TUPLE.index("Delta_R")])
b_delta = np.asarray(mode_b["Delta_R"])
differs = float(np.max(angular_diff(ref_delta, b_delta))) > 1e-9
# The branch-cut case only earns its name if it actually straddles the cut.
cut_delta = np.asarray(call_rust(CASES[4][1])["Delta_R"])
straddles = cut_delta.max() > 3.0 and cut_delta.min() < -3.0
print(f"  branch-cut case: {'PASS' if straddles else 'FAIL'} | Delta_R spans "
      f"[{cut_delta.min():.4f}, {cut_delta.max():.4f}]")
all_pass &= straddles

print(f"  mode mapping: {'PASS' if differs else 'FAIL'} | "
      f"COHERENCY_MATRIX gives a different Delta_R, so FRONT_BLOCK is a real choice")
all_pass &= differs

def cooldown(seconds=1.0):
    """Let the other engine's thread pool go quiet before timing this one.

    numba's default threading layer leaves its workers spin-waiting after a
    call returns, and a spinning worker holds a core while the *other* engine
    is being timed -- so timing the two back to back charges whichever runs
    second. Measured at 20 000 points: the same Rust call is 1.52 ms on its
    own, 2.00 ms immediately after a block of numba calls, 1.58 ms again after
    a one-second pause. The speed line below was reporting that artefact as a
    property of the engine (R5.3).
    """
    time.sleep(seconds)

print(f"\n=== SPEED BENCHMARK: {UNIT_NAME} ===")
bench_case = CASES[-1][1]
call_numba(bench_case)      # warm the jit
call_rust(bench_case)
NUM_RUNS = 50
numba_times, rust_times = [], []
cooldown()
for _ in range(NUM_RUNS):
    t0 = time.perf_counter()
    call_numba(bench_case)
    numba_times.append(time.perf_counter() - t0)
cooldown()
for _ in range(NUM_RUNS):
    t0 = time.perf_counter()
    call_rust(bench_case)
    rust_times.append(time.perf_counter() - t0)
numba_mean_ms = float(np.mean(numba_times)) * 1000
rust_mean_ms = float(np.mean(rust_times)) * 1000
speedup = numba_mean_ms / rust_mean_ms if rust_mean_ms > 0 else float("inf")
print(f"  Numba avg: {numba_mean_ms:.3f} ms")
print(f"  Rust  avg: {rust_mean_ms:.3f} ms")
print(f"  Speedup:   {speedup:.2f}x")
print(f"\nBENCH_RESULT {UNIT_NAME} numba_time_ms={numba_mean_ms:.3f} "
      f"rust_time_ms={rust_mean_ms:.3f} speedup={speedup:.2f}")

report(UNIT_NAME, all_pass)


def test_parity():
    """pytest entry point: the module body above ran the comparisons."""
    assert all_pass, f"{UNIT_NAME}: parity against the numba reference FAILED"


if __name__ == "__main__":
    sys.exit(0 if all_pass else 1)
