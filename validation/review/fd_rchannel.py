# SPDX-License-Identifier: LGPL-3.0-or-later
import numpy as np
from navette.smatrix.smatrix import ScatterMatrix, Request
from navette.smatrix.needle import NeedleRequest, needle_gradient

N = np.array([1.0+0j, 2.35+0j, 1.46+0j, 2.10+0j, 1.52+0j])
D = np.array([0.0, 120.0, 200.0, 80.0, 0.0])
WLS = np.linspace(450.0, 750.0, 31)
st0 = ScatterMatrix(N, D, wavelengths=WLS, angles=[0.0, 30.0])
NN = np.full(WLS.size, 1.8+0.05j, dtype=np.complex128)

def stack_ins(delta):
    idx = np.array([1.0+0j, 2.35+0j, 1.46+0j, 1.8+0.05j, 1.46+0j, 2.10+0j, 1.52+0j])
    thick = np.array([0.0, 120.0, 100.0, delta, 100.0, 80.0, 0.0])
    return ScatterMatrix(idx, thick, wavelengths=WLS, angles=[0.0, 30.0])

# one-sided 2nd-order stencil; converge h until the O(h^2) plateau is
# reached (h=2nm sits above the stencil truncation on this oscillatory
# merit — see docs/code_review.md section 19.1)
analytic_ref = needle_gradient(st0, NN, [220.0], NeedleRequest.P, pol="s")
analytic = 2.0 * np.asarray(analytic_ref["P_s"]).ravel().sum()
h = 2.0
prev = None
rates = []
rels = []
for _ in range(8):
    R0 = np.asarray(st0.compute(Request.RS)["Rs"], float)
    Rp = np.asarray(stack_ins(+h).compute(Request.RS)["Rs"], float)
    R2 = np.asarray(stack_ins(+2 * h).compute(Request.RS)["Rs"], float)
    fdF = np.sum((-3 * R0 ** 2 + 4 * Rp ** 2 - R2 ** 2) / (2 * h))
    rel = abs(analytic - fdF) / max(abs(fdF), abs(analytic))
    print(f"h={h:7.4f}nm: fd={fdF:.10g} rel={rel:.3e}"
          + (f" rate={prev / rel:.1f}" if prev else ""))
    rates.append(prev / rel if prev else 0.0)
    rels.append(rel)
    prev = rel
    h /= 2.0
# pass = clean O(h^2) (rate ~4 on two consecutive halvings) AND the
# converged value agrees with the analytic slope
ok = (len(rates) >= 3 and rates[1] >= 3.5 and rates[2] >= 3.5
      and rels[2] < 1e-3)
print("R-channel slope with distinct needle material:",
      "OK" if ok else "FAIL",
      f"(rate h2->h1={rates[1]:.1f}, h1->h0.5={rates[2]:.1f}, rel={rels[2]:.2e})")
