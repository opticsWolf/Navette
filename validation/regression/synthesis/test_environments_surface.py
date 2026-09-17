# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""F2.4: the Python multi-environment surface and the program document.

`run_needle(design=..., environments=[...])` and a program file's
`design:` / `environments:` sections describe the same thing, and the
whole point of F2.4's shape is that they describe it in the same words:
the keyword call builds the `DesignRequest` a document already is, and
one compiler consumes both. Nothing here tests the physics -- F2.1-F2.3
own that -- these are the seams:

  * the two spellings land on the SAME BITS, not merely on the same
    answer to four decimals;
  * every refusal is a `ValueError`, raised where the mistake was made,
    naming the thing that was wrong and readable on an ASCII terminal;
  * the signature change (`layers` became optional) refuses both and
    neither rather than quietly preferring one.
"""
import json
import struct

import numpy as np
import pytest

from navette.config import load_program
from navette.spectralweave.target import SpectralTarget, TargetCollection
from navette.synthesis import build_merit_spec
from navette.synthesis.pipeline import (
    PipelineConfig,
    design_from_program,
    design_from_program as _bridge,
    run_needle,
)

WL = np.linspace(450.0, 650.0, 11)
H, L, G = 2.1 + 0j, 1.45 + 0j, 1.52 + 0j

#: One macro cycle, one needle, no cleanup and no inflate: the shortest
#: run that still exercises a joint insertion, which is the part of the
#: pipeline a design/environment split actually changes.
CFG = dict(max_macro_cycles=1, needles_per_cycle=1, enable_cleanup=False,
           enable_inflate=False, thin_layer_policy="clamp_up_final")

DOCUMENT = {
    "schema_version": 2, "kind": "program", "name": "two envs",
    "sections": {
        "materials": [
            {"name": "L", "code": "L", "model": "Konstant",
             "params": {"n": 1.45, "k": 0.0}},
            {"name": "H", "code": "H", "model": "Konstant",
             "params": {"n": 2.1, "k": 0.0}},
            {"name": "G", "code": "G", "model": "Konstant",
             "params": {"n": 1.52, "k": 0.0}},
        ],
        "design": {"coat": {"layers": [
            {"material_code": "L", "thickness_nm": 100.0},
            {"material_code": "H", "thickness_nm": 80.0},
        ]}},
        "environments": [
            {"name": "bare", "stack": [{"design": "coat"}]},
            {"name": "laminated", "stack": [
                {"layers": [{"material_code": "G", "thickness_nm": 300.0}]},
                {"design": "coat"},
            ]},
        ],
    },
}


def _targets():
    """The same two demands, one per environment, every time."""
    tc = TargetCollection()
    for env in ("bare", "laminated"):
        tc.add(SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, 0.01),
                              0.0, "s", "R", environment=env))
    return tc


def _hex(x):
    """Little-endian IEEE-754 bits, the project's fingerprint spelling."""
    return struct.pack("<d", float(x)).hex()


def _write(tmp_path, doc):
    p = tmp_path / "prog.json"
    p.write_text(json.dumps(doc), encoding="utf-8")
    return str(p)


# --------------------------------------------------------------------------
# 1. Surface == JSON, hex-compared
# --------------------------------------------------------------------------

def test_a_document_and_a_keyword_call_are_the_same_run(tmp_path):
    """Same design, two spellings, identical bits.

    This is the gate the whole transport design exists to pass. The
    keyword surface does not build a stack and hand it over -- it shapes
    the `DesignRequest` a program file already is, so the file path and
    the keyword path differ only in where the material tables came from.
    If they ever differ in anything else, the merit moves, and a float
    comparison with a tolerance would hide exactly the drift that matters
    (a resample, a reordered library, a flag defaulted differently).
    """
    prog = load_program(_write(tmp_path, DOCUMENT), WL)
    design, environments = design_from_program(prog)

    from_file = run_needle(
        targets=_targets(), angles_deg=[0.0], wavelengths=WL,
        contrast={"L": H, "H": L}, design=design, environments=environments,
        pipeline_config=PipelineConfig(**CFG), substrate=(G, "sub"))

    from_kwargs = run_needle(
        targets=_targets(), angles_deg=[0.0], wavelengths=WL,
        contrast={"L": H, "H": L},
        design={"coat": [(L, 100.0, "L"), (H, 80.0, "H")]},
        environments=[
            {"name": "bare", "stack": [{"design": "coat"}]},
            {"name": "laminated", "stack": [
                {"layers": [(G, 300.0, "G")]}, {"design": "coat"}]},
        ],
        pipeline_config=PipelineConfig(**CFG), substrate=(G, "sub"))

    assert _hex(from_file["final_mf"]) == _hex(from_kwargs["final_mf"]), (
        f"document {_hex(from_file['final_mf'])} "
        f"!= keywords {_hex(from_kwargs['final_mf'])}")
    assert from_file["final_layer_count"] == from_kwargs["final_layer_count"]
    assert from_file["termination"] == from_kwargs["termination"]


def test_the_document_sections_survive_the_load(tmp_path):
    """`design:` / `environments:` come back, in order, unedited.

    A section that parses and then evaporates is the failure the F2.4
    schema bump exists to prevent; asserting the bump without asserting
    the payload would leave that hole open one layer up.
    """
    prog = load_program(_write(tmp_path, DOCUMENT), WL)
    assert list(prog.design) == ["coat"]
    assert [e["name"] for e in prog.environments] == ["bare", "laminated"]
    assert [r["material_code"] for r in prog.design["coat"]["layers"]] == ["L", "H"]
    # Order is evaluation order and is load-bearing: the roster index a
    # demand's `environment=` tag resolves to is a position in this list.
    assert prog.environments[1]["stack"][0]["layers"][0]["material_code"] == "G"


def test_a_program_without_environments_restores_unchanged(tmp_path):
    """The flat case is every pre-F2.4 document, and it stays flat."""
    doc = json.loads(json.dumps(DOCUMENT))
    del doc["sections"]["design"]
    del doc["sections"]["environments"]
    prog = load_program(_write(tmp_path, doc), WL)
    assert prog.design == {} and prog.environments == []


# --------------------------------------------------------------------------
# 2. The signature change (N8)
# --------------------------------------------------------------------------

@pytest.mark.parametrize("kwargs, wanted", [
    (dict(layers=[(L, 100.0)], design={"coat": [(L, 100.0)]}), "both were given"),
    (dict(), "neither was given"),
])
def test_run_needle_refuses_both_and_neither(kwargs, wanted):
    """`layers` became optional; exactly one of the two is still required.

    Defaulting to `layers` when both are present would run the flat path
    and ignore the environments in silence -- an answer for a different
    coating. Defaulting to an empty stack when neither is present would
    optimize nothing and report success.
    """
    with pytest.raises(ValueError, match=wanted) as exc:
        run_needle(targets=None, wavelengths=WL, angles_deg=[0.0], **kwargs)
    msg = str(exc.value)
    assert "run_needle" in msg and "layers" in msg and "design" in msg
    assert msg.isascii(), msg


# --------------------------------------------------------------------------
# 3. Refusals: at construction, ValueError, named, ASCII
# --------------------------------------------------------------------------

def _refuses(fn, *needles):
    with pytest.raises(ValueError) as exc:
        fn()
    msg = str(exc.value)
    assert msg.isascii(), f"non-ASCII refusal: {msg!r}"
    for needle in needles:
        assert needle in msg, f"{needle!r} missing from {msg!r}"
    return msg


def test_an_unknown_environment_tag_refuses_naming_the_roster():
    """Scoring against the wrong surroundings looks like a physics result.

    A typo in an `environment=` tag has no shape error and no exception
    of its own -- it would simply resolve to some other environment, or
    to the first one, and the run would report a merit for a coating the
    user never described. So the roster and the demand meet at
    `build_merit_spec`, and a name that is not on the roster stops there.
    """
    tc = TargetCollection()
    tc.add(SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, 0.01),
                          0.0, "s", "R", environment="laminted"))
    _refuses(lambda: build_merit_spec(tc, environments=["bare", "laminated"]),
             "laminted", "bare", "laminated")


def test_a_tagged_demand_on_a_flat_run_still_names_what_it_knows():
    """One environment is still a roster, and it is called "default".

    The alternative -- no roster at all on the flat path -- would make
    this refusal say that the name is unknown and that there is nothing
    it could have been, which tells the user nothing about the fix.
    """
    tc = TargetCollection()
    tc.add(SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, 0.01),
                          0.0, "s", "R", environment="bare"))
    _refuses(lambda: build_merit_spec(tc), "bare", "default")


def test_an_empty_environment_tag_refuses_at_construction():
    """`environment=""` is a mistake, not a way of saying "the first one".

    Absent means the first environment and is how every pre-F2.1 target
    set stays meaningful. An empty string would resolve to nothing and
    has no second reading worth guessing at.
    """
    # At CONSTRUCTION, not at dump: `__post_init__` validates the target
    # by serializing it, so the tag is checked at the moment the mistake
    # is made rather than at the moment the run starts.
    _refuses(lambda: SpectralTarget(WL, np.zeros_like(WL),
                                    np.full_like(WL, 0.01), 0.0, "s", "R",
                                    environment=""),
             "environment", "non-empty")


def test_an_untagged_demand_emits_no_environment_key():
    """Absent stays absent: a pre-F2.1 target set's JSON does not move."""
    t = SpectralTarget(WL, np.zeros_like(WL), np.full_like(WL, 0.01),
                       0.0, "s", "R")
    assert "environment" not in t._dump()


def test_a_design_row_naming_an_unknown_material_refuses(tmp_path):
    """A `material_code` resolves against the program's own library.

    Same contract a structure's layers have had since Phase 1, and the
    same failure: an unresolvable code is a typo, and the fix is visible
    only if the refusal quotes it.
    """
    doc = json.loads(json.dumps(DOCUMENT))
    doc["sections"]["design"]["coat"]["layers"][0]["material_code"] = "Ll"
    prog = load_program(_write(tmp_path, doc), WL)
    _refuses(lambda: _bridge(prog), "Ll", "coat")


def test_a_design_section_without_environments_refuses(tmp_path):
    """A segment no environment references is unreachable, not optional.

    The loader could carry it and let the compiler find nothing to do
    with it, but by then the filename is gone from the error.
    """
    doc = json.loads(json.dumps(DOCUMENT))
    del doc["sections"]["environments"]
    _refuses(lambda: load_program(_write(tmp_path, doc), WL),
             "design", "environments")


def _surrounded(**extra):
    """The file's request with a surrounding in EVERY environment.

    Environment 0 carries one too, because the returned `stack` is
    environment 0's: a surrounding only environment 1 has cannot be read
    back off the result, and M1's control has to read the flag off the
    very row it is about.
    """
    return dict(
        targets=_targets(), angles_deg=[0.0], wavelengths=WL,
        contrast={"L": H, "H": L},
        design={"coat": [(L, 100.0, "L"), (H, 80.0, "H")]},
        environments=[
            {"name": "bare", "stack": [
                {"layers": [(G, 300.0, "G")]}, {"design": "coat"}]},
            {"name": "laminated", "stack": [
                {"layers": [(G, 600.0, "G2")]}, {"design": "coat"}]},
        ],
        pipeline_config=PipelineConfig(**CFG), substrate=(G, "sub"), **extra)


@pytest.mark.parametrize("flag", ["optimize", "needle"])
def test_a_per_film_flag_cannot_free_a_fixed_surrounding(flag):
    """M1 (review PB, 0.7.6): an override cannot undo a forced flag.

    Flag application runs global map -> row -> per-film override, and the
    override is keyed by MATERIAL CODE, which a surrounding row carries
    just as a design row does. So the override used to land after the row
    schema had forced `optimize`/`needle` false and undo it, in silence:
    measured on 0.7.5, a 500 nm surrounding moved to 534.0298 nm as a
    free variable that only environment 0's residuals pulled on, while
    every other environment kept the compiled thickness and the shared
    design film was driven to its clamp floor. A surrounding is not a
    design variable, so this refuses instead -- at compile, before any
    evaluation, the same refusal `to_row` gives for an explicit `true`.
    """
    msg = _refuses(
        lambda: run_needle(**_surrounded(per_film_flags={"G": {flag: True}})),
        "per_film_flags", flag, "FIXED surrounding row", "design segment")
    # The row it names is the auto-generated surrounding, not only the
    # material code the caller wrote: the caller addressed a material and
    # has to be told which row that reached.
    assert "bare.fixed[0][0]" in msg, msg


def test_the_global_flag_map_still_loses_to_a_fixed_row():
    """M1's control: `film_flags` was never the hole, and stays open.

    The global map is applied BEFORE the row, so the row schema's forced
    false wins on its own and no refusal is needed -- which matters,
    because `film_flags={'optimize': True}` is the ordinary way to say
    "optimize the design" and refusing it would break every segmented run
    that asks for that. Measured, not assumed: the surrounding keeps its
    compiled 300.0 nm and its false, and the design films are free.
    """
    out = run_needle(**_surrounded(film_flags={"optimize": True}))

    films = out["stack"].to_dict()["films"]
    surrounding = [f for f in films if abs(f["thickness"] - 300.0) < 1e-12]
    assert len(surrounding) == 1, [f["thickness"] for f in films]
    assert surrounding[0]["optimize"] is False, surrounding[0]
    assert any(f["optimize"] for f in films), "the map landed on nothing"


def test_duplicate_film_names_refuse_because_a_name_is_an_identity():
    """Two rows, one name: one parameter spelled twice.

    On this surface a film name IS the cross-environment parameter
    identity -- that is how a design segment shared by K environments
    stays one set of thicknesses. Two films of the same physical material
    are therefore two codes carrying identical tables, and letting a name
    repeat would silently tie two layers together.
    """
    _refuses(lambda: run_needle(
        targets=_targets(), angles_deg=[0.0], wavelengths=WL,
        design={"coat": [(L, 100.0, "same"), (H, 80.0, "same")]},
        environments=[{"name": "bare", "stack": [{"design": "coat"}]}],
        pipeline_config=PipelineConfig(**CFG)),
        "same", "unique")
