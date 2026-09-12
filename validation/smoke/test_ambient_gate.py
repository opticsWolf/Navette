# SPDX-License-Identifier: LGPL-3.0-or-later
"""R3.4 (0.6.27): the layer-0 gate, checked at every door that opens on it.

0.6.26 taught ``ScatterMatrix`` to drop absorption in the incident medium and
say so. That covered one of five ways into a solve. The synthesis and design
surface never builds a ``ScatterMatrix`` at all -- it hands its own ``ambient``
straight to the native assembly -- so an absorbing ambient went in silently and
the optimizer then spent thousands of merit evaluations fitting a coating
against physics that does not hold.

The rule now lives once, in the engine
(``optics_core::sanitize_incident_index``), and is applied at each Python-facing
door:

* ``ScatterMatrix``                       -- Python, ``_sanitize_incident_medium``
* ``stack_from_layers``                   -- native ``assemble_design``
* ``run_needle``                          -- native ``run_design``
* ``navette._smatrix.DesignStack(...)``   -- native ``with_films``
* ``DesignStack.from_design(...)``        -- native ``from_design``

Three of those reach ``DesignStack::from_design``, which is the single
production constructor and runs once per assembly. That placement is the point:
the ambient is private with no mutator, so needle insertion, merging, clamping
and thickness steps cannot put the absorption back, and nothing has to be
re-checked inside the optimizer.

What is deliberately NOT gated: the native ``Solver`` and ``core_engine``. R3.1
promised that escape hatch and ``test_the_native_solver_is_still_permissive``
pins it.
"""
from __future__ import annotations

import warnings

import numpy as np

from navette._smatrix import DesignStack, LayerSpec
from navette.smatrix.smatrix import Request, ScatterMatrix
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis.pipeline import run_needle, stack_from_layers

_WLS = np.array([900.0, 1000.0, 1100.0])
_ANGS = np.array([0.0])
_NW = _WLS.size

_AIR = np.full(_NW, 1.0 + 0.0j)
_SUB = np.full(_NW, 1.52 + 0.0j)
_HI = np.full(_NW, 2.35 + 0.0j)
_LO = np.full(_NW, 1.46 + 0.0j)
_ABSORBING_AIR = np.full(_NW, 1.0 + 0.05j)

# The sentence both implementations must carry. Not hard-coded here: it is
# lifted from whichever warning the Python door produces and then required of
# the Rust one, so the test compares the two implementations against each
# other rather than against a third copy that could drift from both.
_SHARED_FROM = "Its absorption has been dropped"


def _only_ambient_warning(record) -> str:
    """The one warning about layer 0, or a failure naming what was said."""
    msgs = [str(w.message) for w in record]
    hits = [m for m in msgs if "incident medium" in m and "absorbing" in m]
    assert len(hits) == 1, f"expected exactly one layer-0 warning, got {msgs!r}"
    return hits[0]


def _catch(fn):
    """Run ``fn`` capturing warnings; return ``(result, [messages])``."""
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        out = fn()
    return out, rec


# ---------------------------------------------------------------------------
# Door 1: the wrapper (0.6.26, re-checked here as the parity reference)
# ---------------------------------------------------------------------------

def _scatter_matrix_warning() -> str:
    idx = np.stack([_ABSORBING_AIR, _HI, _SUB])
    _, rec = _catch(lambda: ScatterMatrix(idx, [0.0, 120.0, 0.0],
                                          wavelengths=_WLS, angles=_ANGS))
    return _only_ambient_warning(rec)


def test_scatter_matrix_still_warns():
    msg = _scatter_matrix_warning()
    assert _SHARED_FROM in msg
    assert msg.isascii(), "the message must survive a cp1252 console"


# ---------------------------------------------------------------------------
# Door 2: stack_from_layers -> assemble_design
# ---------------------------------------------------------------------------

def test_stack_from_layers_drops_the_absorption_and_says_so():
    """The door the 0.6.26 fix missed entirely.

    Before 0.6.27 this returned a stack carrying ``nk = 1+0.05j`` on the
    ambient, with nothing said.
    """
    def build():
        return stack_from_layers([(_HI, 120.0), (_LO, 200.0)], _WLS, {},
                                 ambient=(_ABSORBING_AIR, "air"),
                                 substrate=(_SUB, "sub"))

    (stack, _cmap), rec = _catch(build)
    msg = _only_ambient_warning(rec)
    assert "3 of 3 wavelengths" in msg, msg

    nk = np.asarray(stack.to_dict()["ambient"]["nk"])
    assert np.all(nk.imag == 0.0), f"ambient still absorbs: {nk[:2]}"
    assert np.allclose(nk.real, 1.0), "the real part must be kept, not replaced"


def test_a_transparent_ambient_stays_silent():
    """The common path pays nothing and says nothing."""
    def build():
        return stack_from_layers([(_HI, 120.0)], _WLS, {},
                                 ambient=(_AIR, "air"), substrate=(_SUB, "sub"))

    _, rec = _catch(build)
    assert [str(w.message) for w in rec] == []


def test_only_layer_zero_is_touched():
    """An absorbing substrate or film is ordinary physics and must survive.

    A gate that quietly flattened those would be far worse than no gate.
    """
    abs_sub = np.full(_NW, 1.52 + 0.1j)
    abs_film = np.full(_NW, 2.35 + 0.3j)

    def build():
        return stack_from_layers([(abs_film, 120.0)], _WLS, {},
                                 ambient=(_ABSORBING_AIR, "air"),
                                 substrate=(abs_sub, "sub"))

    (stack, _cmap), rec = _catch(build)
    _only_ambient_warning(rec)
    d = stack.to_dict()
    assert np.allclose(np.asarray(d["substrate"]["nk"]).imag, 0.1)
    assert np.allclose(np.asarray(d["films"][0]["nk"]).imag, 0.3)


# ---------------------------------------------------------------------------
# Door 3: run_needle -> run_design (the path that swallowed its warnings)
# ---------------------------------------------------------------------------

def test_a_full_design_run_is_not_the_silent_one():
    """``run_design`` discarded its assembly warnings until 0.6.27.

    That made the most expensive path the quietest one: a needle design would
    correct the ambient and never mention it, having already spent thousands
    of merit evaluations. The same line was swallowing the graded-film
    homogenization warning, which is why this is a fix and not a new feature.
    """
    from navette._smatrix import NeedleCycleConfig, PipelineConfig

    tc = TargetCollection()
    tc.add(SpectralTarget(_WLS, np.zeros(_NW), np.full(_NW, 0.01), 0.0,
                          "s", "R", kind="e", weight=1.0))
    cfg = PipelineConfig(max_macro_cycles=1, needles_per_cycle=1,
                         enable_cleanup=False, enable_inflate=False,
                         stagnation_window=100)
    nc = NeedleCycleConfig(scan_step_nm=20.0, refold_per_cycle=False)

    def run():
        return run_needle([(_LO, 120.0)], tc, _ANGS, _WLS,
                          {"film0": _HI}, pipeline_config=cfg, needle_config=nc,
                          names=["L"], ambient=(_ABSORBING_AIR, "air"),
                          substrate=(_SUB, "sub"))

    res, rec = _catch(run)
    _only_ambient_warning(rec)
    nk = np.asarray(res["stack"].to_dict()["ambient"]["nk"])
    assert np.all(nk.imag == 0.0)


# ---------------------------------------------------------------------------
# Doors 4 and 5: the native DesignStack constructors
# ---------------------------------------------------------------------------

def test_the_plain_design_stack_constructor_is_gated():
    """``DesignStack(ambient, substrate, films)`` reaches ``with_films``.

    It does not pass through ``from_design``, so it needs the rule applied at
    the boundary. Missing this one is the same mistake as 0.6.26's: gating the
    door you were looking at rather than all of them.
    """
    def build():
        return DesignStack(
            LayerSpec("air", _ABSORBING_AIR, 0.0, optimize=False, needle=False),
            LayerSpec("sub", _SUB, 0.0, optimize=False, needle=False),
            [LayerSpec("H", _HI, 120.0)],
        )

    stack, rec = _catch(build)
    _only_ambient_warning(rec)
    assert np.all(np.asarray(stack.to_dict()["ambient"]["nk"]).imag == 0.0)


def test_the_expanded_design_stack_constructor_is_gated():
    from navette._structure import Layer

    def build():
        return DesignStack.from_design(
            LayerSpec("air", _ABSORBING_AIR, 0.0, optimize=False, needle=False),
            LayerSpec("sub", _SUB, 0.0, optimize=False, needle=False),
            [Layer(thickness=120.0, material_name="H")],
            {"H": _HI},
            {},
            _WLS,
        )

    stack, rec = _catch(build)
    _only_ambient_warning(rec)
    assert np.all(np.asarray(stack.to_dict()["ambient"]["nk"]).imag == 0.0)


# ---------------------------------------------------------------------------
# The two implementations must not drift apart
# ---------------------------------------------------------------------------

def test_both_doors_explain_it_the_same_way():
    """One rule, two implementations (Python wrapper, Rust engine).

    They are kept separate on purpose -- the Python warning's ``stacklevel``
    points at the caller's own constructor, which a warning raised from Rust
    cannot do -- so the explanation is the thing that has to stay identical.
    Whichever one is edited, this fails until the other follows.
    """
    py_msg = _scatter_matrix_warning()
    shared = py_msg[py_msg.index(_SHARED_FROM):]
    assert len(shared) > 300, f"parity slice looks truncated: {shared!r}"

    def build():
        return stack_from_layers([(_HI, 120.0)], _WLS, {},
                                 ambient=(_ABSORBING_AIR, "air"),
                                 substrate=(_SUB, "sub"))

    _, rec = _catch(build)
    rust_msg = _only_ambient_warning(rec)
    assert shared in rust_msg, (
        "the two doors have drifted apart\n"
        f"python: {shared!r}\n"
        f"rust  : {rust_msg!r}"
    )
    assert rust_msg.isascii(), "the message must survive a cp1252 console"

    # Both must still name what was dropped, not merely that something was.
    for msg in (py_msg, rust_msg):
        for fragment in ("incident medium", "layer 0", "wavelengths",
                         "largest |Im(n)|", "transparent ambient", "substrate"):
            assert fragment in msg, f"{fragment!r} missing from {msg!r}"


# ---------------------------------------------------------------------------
# The gate is at assembly, not in the optimizer
# ---------------------------------------------------------------------------

def test_the_correction_survives_every_stack_mutation():
    """Once at assembly is enough, which is why it is not on the hot path.

    The ambient is private with no mutator, so the operations the pipeline
    runs thousands of times per design cannot reintroduce the absorption.
    Checking per merit evaluation would re-establish something that cannot
    have changed.
    """
    def build():
        return stack_from_layers([(_LO, 200.0), (_HI, 100.0)], _WLS, {},
                                 ambient=(_ABSORBING_AIR, "air"),
                                 substrate=(_SUB, "sub"))

    (stack, _cmap), rec = _catch(build)
    _only_ambient_warning(rec)

    stack.insert_needle_seed(0, 100.0, LayerSpec("H", _HI, 5.0))
    stack.merge_adjacent()
    stack.clamp_all(1.0, 1000.0)
    stack.set_thickness(0, 33.0)
    nk = np.asarray(stack.to_dict()["ambient"]["nk"])
    assert np.all(nk.imag == 0.0), f"absorption came back: {nk[:2]}"


def test_the_gated_stack_matches_its_transparent_twin():
    """The correction is exact: same numbers as building it transparent.

    Bit-for-bit, not ``allclose`` -- "we dropped k" has one right answer, and
    any drift means something other than k was touched.
    """
    def build(ambient):
        return stack_from_layers([(_HI, 120.0), (_LO, 200.0)], _WLS, {},
                                 ambient=(ambient, "air"), substrate=(_SUB, "sub"))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        corrected, _ = build(_ABSORBING_AIR)
    twin, _ = build(_AIR)

    a = np.asarray(corrected.to_dict()["ambient"]["nk"])
    b = np.asarray(twin.to_dict()["ambient"]["nk"])
    assert np.array_equal(a, b), f"{a[:2]} != {b[:2]}"


def test_the_native_solver_is_still_permissive():
    """R3.1's escape hatch stays open.

    The gate is a property of the design and ``ScatterMatrix`` surfaces, not
    of the engine. Someone who knows what an absorbing ambient does to their
    amplitudes can still reach the raw ones.
    """
    from navette._smatrix import Solver

    idx = np.concatenate([_ABSORBING_AIR, _HI, _SUB])
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        s = Solver(_WLS, _ANGS, idx, 3, thicknesses=np.array([0.0, 120.0, 0.0]))
        out = s.solve(int(Request.RS))
    assert [str(w.message) for w in rec] == []
    assert np.all(np.isfinite(np.asarray(out["Rs"], dtype=float)))
