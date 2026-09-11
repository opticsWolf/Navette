# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""spectralweave: weave fragments into one curve, then unweave a target back.

A simulation rarely produces one spectrum on one grid. It produces pieces --
a dense scan across a resonance, a coarse sweep either side of it, a handful
of points somebody measured -- each on its own wavelength grid. `spectralweave`
holds those pieces as *frames*, hands you the concatenated curve to work with,
and slices an edited curve back into the frames it came from.

Run:  python examples/spectralweave_example.py
"""

import numpy as np

from navette.spectralweave import OpticalFragment, SimulationWeaver


def main() -> None:
    weaver = SimulationWeaver(cache_size=64)

    # --- weave: three fragments of the same curve, on three grids ----------
    # They must not overlap: each wavelength belongs to exactly one frame.
    grids = [
        np.linspace(400.0, 499.0, 100),   # coarse, blue side
        np.linspace(500.0, 600.0, 401),   # dense, across the feature
        np.linspace(601.0, 800.0, 200),   # coarse, red side
    ]
    for wl in grids:
        weaver.add_fragment(OpticalFragment(
            wavelengths=wl,
            values=0.5 + 0.4 * np.sin(wl / 40.0),
            base_wavelength=550.0,      # \
            data_type="R",              #  > the three fields that form the key
            polarization="s",           # /
        ))

    print(f"frames held: {weaver.frame_count}")

    # --- the continuous curve ---------------------------------------------
    wl, values = weaver.get_continuous_curve(550.0, "R", "s")
    print(f"weaved curve: {wl.size} points, {wl[0]:.1f}-{wl[-1]:.1f} nm, "
          f"R in [{values.min():.3f}, {values.max():.3f}]")

    # --- unweave: push an edited target back into the same frames ----------
    # `wl` is the weaved grid, so every frame is covered exactly. Distribution
    # is by exact wavelength match -- a fresh linspace over the same range
    # would share almost no points and be rejected (see the guard below).
    target = np.clip(values * 1.05, 0.0, 1.0)

    template = OpticalFragment(
        wavelengths=wl, values=target,
        base_wavelength=550.0, data_type="R", polarization="s",
    )
    written = weaver.unweave(template, wl, target)
    print(f"unweave wrote {written} fragment(s) back")

    wl_after, values_after = weaver.get_continuous_curve(550.0, "R", "s")
    assert np.array_equal(wl_after, wl), "the grid must survive a round trip"
    assert np.allclose(values_after, target), "values must survive a round trip"
    print("round trip: grid and values match")

    # --- what happens when the grid does not line up ----------------------
    stranger = np.linspace(400.0, 800.0, 500)
    try:
        weaver.unweave(template, stranger, np.cos(stranger / 100.0))
    except ValueError as exc:
        print(f"rejected, as it should be: {exc}")

    # --- batch form: many keys, one grid, one FFI crossing ------------------
    keys = {
        OpticalFragment(wavelengths=wl, values=target, base_wavelength=550.0,
                        data_type="R", polarization=pol): target * scale
        for pol, scale in (("s", 1.0), ("p", 0.9))
    }
    print(f"unweave_batch wrote {weaver.unweave_batch(wl, keys)} fragment(s)")


if __name__ == "__main__":
    main()
