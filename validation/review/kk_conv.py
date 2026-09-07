import numpy as np
from kk_validate import eps1_pair, eps2_lor
# convergence of MY quadrature at the TL comparison point E=4.13 (Lorentz proxy)
def eps2_tl(E, Eg, oscs):
    E = np.asarray(E, float)
    out = np.zeros_like(E)
    for (A, E0, C) in oscs:
        m = E > Eg
        out[m] += (A * E0 * C * (E[m] - Eg) ** 2
                   / (((E[m] ** 2 - E0 ** 2) ** 2 + C ** 2 * E[m] ** 2) * E[m]))
    return out
args = (3.0, [(25.0, 4.0, 1.2), (8.0, 9.0, 2.5)])
prev = None
for h in (2e-3, 1e-3, 5e-4, 2.5e-4):
    v = eps1_pair(eps2_tl, args, 4.13, h=h)
    print(f"h={h:.1e}: mine={v:.6f}" + (f"  step={abs(v-prev):.2e}" if prev else ""))
    prev = v
# Lorentz analytic error at same h values
print("Lorentz analytic check at E=4.13:")
for h in (2e-3, 5e-4, 1.25e-4):
    m = eps1_pair(eps2_lor, (4.0, 1.2, 25.0), 4.13, h=h)
    ex = 1.0 + 25.0 * (16.0 - 4.13**2) / ((4.13**2 - 16.0)**2 + 1.44 * 4.13**2)
    print(f"  h={h:.1e}: mine={m:.6f} exact={ex:.6f} rel={abs(m-ex)/abs(ex):.2e}")
