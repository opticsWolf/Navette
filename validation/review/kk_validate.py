import numpy as np

# analytic Lorentz oscillator: eps2 = f*G*E/((E^2-E0^2)^2 + G^2 E^2)
# KK: eps1 = 1 + f*(E0^2 - E^2)/((E^2-E0^2)^2 + G^2 E^2)   [standard result]
def eps2_lor(E, E0, G, f):
    return f * G * E / ((E * E - E0 * E0) ** 2 + G * G * E * E)
def eps1_lor_exact(E, E0, G, f):
    return 1.0 + f * (E0 * E0 - E * E) / ((E * E - E0 * E0) ** 2 + G * G * E * E)

def eps2_of(fn, x, *a):
    return fn(x, *a)

def eps1_pair(fn, args, E, eps_inf=1.0, E_hi=80.0, h=2e-3):
    # + branch: xi in (E, E_hi); - branch: xi in (0, E). Symmetric pairs
    # (same u) cancel the PV singular part exactly; beyond umin = min(E, E_hi-E)
    # the - branch ends, so the remaining + branch integrates alone.
    umin = min(E, E_hi - E)
    npair = int(umin / h)
    up = (np.arange(npair) + 0.5) * h
    xip, xim = E + up, E - up
    gp, gm = fn(xip, *args), fn(xim, *args)
    dp, dm = xip * xip - E * E, xim * xim - E * E
    pair = gp * xip / dp + gm * xim / dm
    total = np.trapezoid(pair, up)
    rest = int((E_hi - E - umin) / h)
    if rest > 2:
        ur = umin + (np.arange(rest) + 0.5) * h
        xr = E + ur
        gr = fn(xr, *args)
        dr = xr * xr - E * E
        total += np.trapezoid(gr * xr / dr, ur)
    return eps_inf + (2.0 / np.pi) * total


def main():
    print("validate my KK against analytic Lorentz:")
    for E in (2.0, 3.5, 4.13, 6.0):
        args = (4.0, 1.2, 25.0)
        mine = eps1_pair(eps2_lor, args, E)
        exact = eps1_lor_exact(E, *args)
        print(f"  E={E}: mine={mine:.6f} exact={exact:.6f} rel={abs(mine-exact)/abs(exact):.2e}")
    
    # now map the TL discrepancy around E0=4.0
    from navette.materials import MaterialSpec, evaluate
    HC = 1239.841984
    def eps2_tl(E, Eg, oscs):
        E = np.asarray(E, float)
        out = np.zeros_like(E)
        for (A, E0, C) in oscs:
            m = E > Eg
            out[m] += (A * E0 * C * (E[m] - Eg) ** 2
                       / (((E[m] ** 2 - E0 ** 2) ** 2 + C ** 2 * E[m] ** 2) * E[m]))
        return out
    spec = MaterialSpec("TaucLorentz", {"Eg": 3.0,
        "osc": [(25.0, 4.0, 1.2), (8.0, 9.0, 2.5)], "epsilon_inf": 1.0})
    Es = np.array([3.2, 3.5, 3.8, 4.0, 4.05, 4.13, 4.3, 5.0, 6.0, 10.0])
    wls = HC / Es
    nk = evaluate(spec, wls)
    e1 = np.real(nk) ** 2 - np.imag(nk) ** 2
    e2 = 2.0 * np.real(nk) * np.imag(nk)
    e2_mine = eps2_tl(Es, 3.0, [(25.0, 4.0, 1.2), (8.0, 9.0, 2.5)])
    print("\nTL engine vs my pair-KK (same eps2):")
    print("  eps2 agree:", np.allclose(e2, e2_mine, rtol=1e-7, atol=1e-10))
    for i, E in enumerate(Es):
        mine = eps1_pair(eps2_tl, (3.0, [(25.0, 4.0, 1.2), (8.0, 9.0, 2.5)]), E)
        print(f"  E={E:5.2f}: engine={e1[i]:.6f} mine={mine:.6f} rel={abs(e1[i]-mine)/abs(mine):.2e}")


if __name__ == "__main__":
    main()
