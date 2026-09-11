#!/usr/bin/env python3
"""Comparison test for core_engine_photometry_only — numba vs rust (R2.4a).

The legacy kernel took ``calc_s, calc_p`` and always returned four arrays,
zero-filling the polarization it was told to skip. Its successor is the
request-driven ``navette._smatrix.core_engine``, where the same choice is two
bits in a mask and a skipped channel is *absent from the result* rather than
present and zero. Porting therefore means three things, and each is a claim
the comparisons below have to earn:

1. **The mask means what ``calc_s``/``calc_p`` meant.** ``Request.RS|TS``
   against ``calc_s=1, calc_p=0``, and so on for the three combinations.
2. **The zeros were never data.** A channel the legacy filled with zeros is a
   key the new engine does not emit, so the test asserts absence rather than
   comparing against zeros — which would have passed for the wrong reason.
3. **Asking for less does not change what you get.** `Rs` from the s-only
   request must be bit-identical to `Rs` from the both-polarizations request;
   the fast path is a fast path, not a different physics.

``coherence_mode`` is ``FRONT_BLOCK`` (0), the legacy engine's own
partial-coherence treatment: incoherent boundaries cascade as intensities and
only the front block keeps its phase. Stated precisely, because it is easy to
overclaim here: for *pure intensities* ``FRONT_BLOCK`` and
``COHERENCY_MATRIX`` are the same computation — the coherency channel is only
tracked when something asks for it, and photometry never does. The mode
distinction this file can make is against ``FULLY_COHERENT``, which ignores
the flags altogether, and it is asserted at the bottom: identical to the other
two when every layer is coherent, different as soon as one is not. (The A-vs-B
distinction is real for ellipsometry and is asserted in that port.)
"""

import sys, os, time
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
UNIT_NAME = "core_engine_photometry_only"

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
    from loom_matrix import core_engine_photometry_only as numba_func
    numba_imported = True
except ImportError:
    numba_func = None

require_reference(numba_imported, "numba reference loom_matrix",
                  "pip install numba; see validation/parity/smatrix/refs/")

# Unified layout: kernels live in navette._smatrix (aggregated extension).
import navette._smatrix as rust_mod  # noqa: E402
from navette.smatrix.smatrix import CoherenceMode, Request  # noqa: E402

rust_func = rust_mod.core_engine

# The legacy positional `calc_s, calc_p` as request bits. `R_avg`/`T_avg` are
# deliberately not asked for: the legacy kernel does not compute them, so
# requesting them would compare a channel with no reference.
POL_CASES = [
    ("s+p", 1, 1, Request.RS | Request.RP | Request.TS | Request.TP),
    ("s",   1, 0, Request.RS | Request.TS),
    ("p",   0, 1, Request.RP | Request.TP),
]

# Which dict key carries each element of the legacy return tuple.
LEGACY_TUPLE = ["Rs", "Rp", "Ts", "Tp"]
# ...and which polarization each of them belongs to, so a case that switched
# one off knows which keys must be missing.
KEY_POL = {"Rs": "s", "Ts": "s", "Rp": "p", "Tp": "p"}


def generate_test_data(rng, incoherent_layer=None, rough=False):
    """One stack, in both engines' input layouts.

    The numba kernel wants a complex (n_wavs, n_layers) index cache; the
    binding wants the same numbers interleaved re/im and flattened wav-major.
    They are built from one source here so a layout bug cannot masquerade as
    a physics difference.
    """
    wavls = np.linspace(400.0, 800.0, NUM_WAVS).astype(np.float64)
    angles_deg = np.linspace(-30.0, 30.0, NUM_ANGLES).astype(np.float64)
    sin_theta_arr = np.sin(np.radians(angles_deg)).astype(np.float64)

    n_real_base = [1.0, 1.5 + rng.random() * 0.3, 2.0 + rng.random() * 0.5,
                   1.8 + rng.random() * 0.4, 2.5 + rng.random() * 0.3, 1.0]
    n_imag_base = [0.0, 0.001 + rng.random() * 0.01, 0.01 + rng.random() * 0.05,
                   0.001 + rng.random() * 0.02, 0.02 + rng.random() * 0.08, 0.0]

    n_complex = np.empty((NUM_WAVS, N_LAYERS), dtype=np.complex128)
    for wv in range(NUM_WAVS):
        wav_factor = (wavls[wv] - 400.0) / 400.0
        for li in range(N_LAYERS):
            n_complex[wv, li] = complex(
                n_real_base[li] + wav_factor * rng.uniform(-0.1, 0.1),
                max(0.0, n_imag_base[li]),
            )
    n_flat = np.empty(NUM_WAVS * N_LAYERS * 2, dtype=np.float64)
    n_flat[0::2] = n_complex.real.ravel()
    n_flat[1::2] = n_complex.imag.ravel()

    thicknesses = np.array([50.0, 100.0, 200.0, 150.0, 80.0, 500.0],
                           dtype=np.float64)
    # Length n_layers: the binding requires one flag per layer, and the numba
    # kernel only ever reads the first n_layers-1 of them.
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


def call_numba(case, calc_s, calc_p):
    return numba_func(
        case["wavls"], case["sin_theta_arr"], N_LAYERS, case["n_complex"],
        case["thicknesses"], case["incoherent_flags"], case["rough_types"],
        case["rough_vals"], np.int32(calc_s), np.int32(calc_p),
    )


def call_rust(case, requested, mode=CoherenceMode.FRONT_BLOCK):
    return rust_func(
        case["wavls"], case["sin_theta_arr"], N_LAYERS, case["n_flat"],
        case["thicknesses"], case["incoherent_flags"], case["rough_types"],
        case["rough_vals"], int(mode), int(requested),
    )


rng = np.random.default_rng(seed=42)
CASES = [
    ("all coherent",        generate_test_data(rng)),
    ("all coherent, rough", generate_test_data(rng, rough=True)),
    ("incoherent layer 3",  generate_test_data(rng, incoherent_layer=3)),
    ("incoherent + rough",  generate_test_data(rng, incoherent_layer=2, rough=True)),
]

print(f"=== CORRECTNESS TEST: {UNIT_NAME} ===")
all_pass = True
for label, case in CASES:
    for pol_label, calc_s, calc_p, requested in POL_CASES:
        numba_out = call_numba(case, calc_s, calc_p)
        rust_out = call_rust(case, requested)
        tag = f"{label} [{pol_label}]"

        for key, ref in zip(LEGACY_TUPLE, numba_out):
            wanted = (KEY_POL[key] == "s" and calc_s) or (KEY_POL[key] == "p" and calc_p)
            if not wanted:
                # The legacy kernel zero-filled this; the mask omits it. Never
                # compare against those zeros -- that passes for the wrong
                # reason and would keep passing if the mask were ignored.
                absent = key not in rust_out
                print(f"  {tag} {key}: {'PASS' if absent else 'FAIL'} | "
                      f"not requested, {'absent' if absent else 'PRESENT'}")
                all_pass &= absent
                continue
            got = np.asarray(rust_out[key])
            ref = np.asarray(ref)
            ok = got.shape == ref.shape and np.allclose(ref, got, rtol=1e-9, atol=1e-12)
            diff = float(np.max(np.abs(ref - got))) if got.shape == ref.shape else float("nan")
            print(f"  {tag} {key}: {'PASS' if ok else 'FAIL'} | diff_max={diff:.2e}")
            all_pass &= ok

    # Asking for one polarization must not change the other's numbers.
    both = call_rust(case, POL_CASES[0][3])
    s_only = call_rust(case, POL_CASES[1][3])
    p_only = call_rust(case, POL_CASES[2][3])
    for key, sub in (("Rs", s_only), ("Ts", s_only), ("Rp", p_only), ("Tp", p_only)):
        same = np.array_equal(np.asarray(both[key]), np.asarray(sub[key]))
        print(f"  {label} {key}: {'PASS' if same else 'FAIL'} | "
              f"single-pol request is bit-identical to the full one")
        all_pass &= same

# What the coherence mode does to photometry, asserted rather than assumed.
for label, case in CASES:
    flagged = bool(case["incoherent_flags"].any())
    req = POL_CASES[0][3]
    a = call_rust(case, req, CoherenceMode.FRONT_BLOCK)
    b = call_rust(case, req, CoherenceMode.COHERENCY_MATRIX)
    c = call_rust(case, req, CoherenceMode.FULLY_COHERENT)
    # No cross channel is requested, so A and B are the same computation.
    same_ab = all(np.array_equal(np.asarray(a[k]), np.asarray(b[k]))
                  for k in LEGACY_TUPLE)
    print(f"  {label} modes A/B: {'PASS' if same_ab else 'FAIL'} | "
          f"identical without a cross-channel request")
    all_pass &= same_ab
    # C ignores the flags: identical when there are none, different otherwise.
    same_c = all(np.array_equal(np.asarray(a[k]), np.asarray(c[k]))
                 for k in LEGACY_TUPLE)
    ok = same_c != flagged
    print(f"  {label} mode C: {'PASS' if ok else 'FAIL'} | fully-coherent "
          f"{'differs' if not same_c else 'matches'}, flags "
          f"{'set' if flagged else 'clear'}")
    all_pass &= ok

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
requested = POL_CASES[0][3]
call_numba(bench_case, 1, 1)          # warm the jit
call_rust(bench_case, requested)
NUM_RUNS = 50
numba_times, rust_times = [], []
cooldown()
for _ in range(NUM_RUNS):
    t0 = time.perf_counter()
    call_numba(bench_case, 1, 1)
    numba_times.append(time.perf_counter() - t0)
cooldown()
for _ in range(NUM_RUNS):
    t0 = time.perf_counter()
    call_rust(bench_case, requested)
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
