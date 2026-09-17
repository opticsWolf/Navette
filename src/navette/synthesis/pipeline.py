# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""navette.synthesis.pipeline — needle synthesis driver.

Thin wrapper over the native pipeline (``navette._smatrix``): target
definition stays on the spectralweave surface — pass a
:class:`~navette.spectralweave.target.TargetCollection` (or a ready-made
``MeritSpec``) and a layer list; this module evaluates materials onto the
simulation grid, builds the native ``DesignStack`` + ``NeedlePipeline``,
and runs the macro-loop. All target options (kind/band/phase/weight/
count/integral, PD labels) flow through ``build_merit_spec`` untouched.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

from navette._smatrix import (
    DesignStack as _NativeDesignStack,
    LayerSpec as _NativeLayerSpec,
    LmConfig as _NativeLmConfig,
    NeedleCycleConfig as _NativeNeedleCycleConfig,
    NeedlePipeline as _NativePipeline,
    PipelineConfig as _NativePipelineConfig,
    SmatrixContext as _NativeContext,
)
from navette.materials import MaterialSpec, evaluate
from navette.spectralweave.target import TargetCollection
from navette.synthesis import build_merit_spec as _build_merit_spec

__all__ = [
    "LayerSpec",
    "DesignStack",
    "SmatrixContext",
    "LmConfig",
    "PipelineConfig",
    "NeedleCycleConfig",
    "NeedlePipeline",
    "layer_from_material",
    "stack_from_layers",
    "run_needle",
    "design_from_program",
]

# Re-export native classes under friendly names.
LayerSpec = _NativeLayerSpec
DesignStack = _NativeDesignStack
SmatrixContext = _NativeContext
LmConfig = _NativeLmConfig
PipelineConfig = _NativePipelineConfig
NeedleCycleConfig = _NativeNeedleCycleConfig
NeedlePipeline = _NativePipeline


def _eval_nk(material: Any, wavelengths: np.ndarray) -> np.ndarray:
    """MaterialSpec (or mapping) → complex nk on ``wavelengths``."""
    if isinstance(material, MaterialSpec):
        return np.ascontiguousarray(evaluate(material, wavelengths))
    if isinstance(material, Mapping):
        return np.ascontiguousarray(
            evaluate(MaterialSpec(model=material["model"],
                                  params=dict(material.get("params", {}))),
                     wavelengths)
        )
    arr = np.ascontiguousarray(np.asarray(material, dtype=np.complex128))
    if arr.shape != wavelengths.shape:
        raise ValueError(
            f"nk array shape {arr.shape} != wavelengths shape {wavelengths.shape}."
        )
    return arr


def layer_from_material(material: Any, thickness: float, wavelengths,
                        name: str = "", **flags) -> Any:
    """One native ``LayerSpec`` from a spec/mapping/nk array + thickness.

    ``flags``: coherent, rough_type, rough_val, optimize, needle.
    """
    wl = np.ascontiguousarray(np.asarray(wavelengths, dtype=np.float64))
    return LayerSpec(str(name), _eval_nk(material, wl), float(thickness), **flags)


_LAYER_KEYS = ("roughness", "rough_type", "interface", "interface_thickness")
_FILM_DEFAULTS = dict(coherent=True, optimize=True, needle=True,
                      roughness=0.0, rough_type=0, inhomogen=False,
                      inh_delta=0.1, interface=False,
                      interface_thickness=0.0,
                      gradient=None)

_EMA_RULE_NAMES = ("Bruggeman", "MaxwellGarnett", "Looyenga",
                   "Lichtenecker", "MoriTanaka", "PowerLaw")


def _norm_gradient(g, film_name, wl, film_nk):
    """Shape one film's ``gradient`` flag into the native dict (F1.1/A5).

    The design path has no materials library - every nk is film-supplied
    - so the inclusion spectrum is evaluated HERE and rides the film
    dict (``nk_b``), registered under a provider key of its own
    (``<film>~b``, overridable via ``b_name``). ``material_a`` defaults
    to the film's own name: its registered nk IS the host spectrum.
    ``ema`` is a rule name (defaults per variant) or a one-key
    ``{name: params}`` map. The profile mode needs exactly one slope
    spelling (A9.2): ``f_end`` (FixedSpan - the F1.1 mode) or ``rate``
    (RateCapped - F1.2: ``f(z) = f_start + rate*z/ref_thickness``
    clamped to ``[f_min, f_max]``; supplying both is refused).
    """
    if not isinstance(g, Mapping):
        raise TypeError(f"film {film_name!r}: gradient must be a mapping.")
    g = dict(g)
    if "material_b" not in g:
        raise ValueError(f"film {film_name!r}: gradient requires 'material_b'.")
    mb = g["material_b"]
    if isinstance(mb, str):
        raise ValueError(
            f"film {film_name!r}: gradient material_b must be a material "
            "(spec/mapping/nk array) - the design path has no materials "
            "library to resolve a bare name; the host is the film itself "
            "(material_a defaults to the film's own name)."
        )
    ema = g.get("ema", "Bruggeman")
    if isinstance(ema, str):
        if ema not in _EMA_RULE_NAMES:
            raise ValueError(
                f"film {film_name!r}: unknown mixing rule {ema!r} "
                f"(one of {', '.join(_EMA_RULE_NAMES)})."
            )
    elif isinstance(ema, Mapping):
        if len(ema) != 1:
            raise ValueError(
                f"film {film_name!r}: ema map must have exactly one key "
                "(the rule name)."
            )
        ((name, params),) = ema.items()
        if name not in _EMA_RULE_NAMES:
            raise ValueError(
                f"film {film_name!r}: unknown mixing rule {name!r} "
                f"(one of {', '.join(_EMA_RULE_NAMES)})."
            )
        if params is not None and not isinstance(params, Mapping):
            raise ValueError(f"film {film_name!r}: ema params must be a mapping.")
    else:
        raise TypeError(f"film {film_name!r}: ema must be a name or a one-key mapping.")
    nk_b = _eval_nk(mb, wl)
    if np.array_equal(nk_b, film_nk):
        raise ValueError(
            f"film {film_name!r}: gradient materials are identical "
            "(the film's own nk IS material_b's) - a gradient of a material "
            "with itself is a spec bug; single-material drift is inhomogen."
        )
    # F1.2: the profile mode - exactly one slope spelling (A9.2: two
    # slope-like keys are refused, no silent precedence).
    has_rate = "rate" in g
    has_f_end = "f_end" in g
    if has_rate and has_f_end:
        raise ValueError(
            f"film {film_name!r}: gradient: 'rate' and 'f_end' are the same "
            "slope - give one."
        )
    if has_rate:
        mode = {"RateCapped": {
            "f_start": g.get("f_start", 0.0),
            "rate": g["rate"],
            "ref_thickness": g.get("ref_thickness", 100.0),
            "f_min": g.get("f_min", 0.0),
            "f_max": g.get("f_max", 1.0),
        }}
        for key in ("rate", "ref_thickness", "f_min", "f_max", "f_start"):
            if key in mode["RateCapped"] and not isinstance(mode["RateCapped"][key], (int, float)):
                raise TypeError(f"film {film_name!r}: gradient {key} must be a number.")
    elif has_f_end:
        mode = {"FixedSpan": {"f_start": g.get("f_start", 0.0), "f_end": g["f_end"]}}
    else:
        raise ValueError(
            f"film {film_name!r}: gradient requires 'f_end' (FixedSpan) or "
            "'rate' (RateCapped)."
        )
    g["ema"] = ema
    g["mode"] = mode
    g["material_a"] = g.get("material_a", film_name)
    g["material_b"] = g.get("b_name", f"{film_name}~b")
    g["nk_b"] = [(float(z.real), float(z.imag)) for z in nk_b]
    return g


def _film_dicts(layers, names, wl, film_flags=None, per_film_flags=None):
    """Shape ``[(material, d)]`` films into native ``ArrayFilm`` dicts.

    Merge order per film: defaults < global ``film_flags`` < per-film
    override (``rough_val`` aliases ``roughness``); materials evaluate
    to nk on ``wl``. Shared by ``stack_from_layers`` and ``run_needle``.

    F1.1: the ``gradient`` flag carries a mixture gradient. Its
    ``material_b`` is an evaluable material (evaluated HERE - the design
    path has no materials library, so the inclusion spectrum rides the
    film dict as ``nk_b``, registered under a ``<film>~b`` provider
    key); ``material_a`` defaults to the film's own name, whose
    registered nk IS the host spectrum. Background-pinned gradient
    films (``optimize=False, needle=False``) expand WITH the profile
    (inclusive endpoint sampling: the first sublayer IS ``f_start``,
    the last IS ``f_end``); any other gradient film homogenizes to the
    EMA at ``f_mid`` with a warning (R4 - a gradient film has no base
    index to flip to).
    """
    flags = dict(film_flags or {})
    if "rough_val" in flags:
        flags["roughness"] = flags.pop("rough_val")
    unknown = set(flags) - set(_FILM_DEFAULTS)
    if unknown:
        raise TypeError(f"unknown film flags {sorted(unknown)}.")
    overrides = per_film_flags or {}
    out = []
    for (mat, d), nm in zip(layers, names):
        fd = dict(_FILM_DEFAULTS)
        fd.update(flags)
        local = dict(overrides.get(nm, {}))
        if "rough_val" in local:
            local["roughness"] = local.pop("rough_val")
        unknown = set(local) - set(_FILM_DEFAULTS)
        if unknown:
            raise TypeError(f"film {nm!r}: unknown flags {sorted(unknown)}.")
        fd.update(local)
        row_nk = _eval_nk(mat, wl)
        row = {
            "name": str(nm), "nk": row_nk, "d_nm": float(d),
            "coherent": bool(fd["coherent"]),
            "roughness": float(fd["roughness"]),
            "rough_type": int(fd["rough_type"]),
            "inhomogen": bool(fd["inhomogen"]),
            "inh_delta": float(fd["inh_delta"]),
            "interface": bool(fd["interface"]),
            "interface_thickness": float(fd["interface_thickness"]),
            "optimize": bool(fd["optimize"]),
            "needle": bool(fd["needle"]),
        }
        if fd["gradient"] is not None:
            row["gradient"] = _norm_gradient(fd["gradient"], str(nm), wl, row_nk)
        out.append(row)
    return out


def _half_nk(value, wl):
    """Half-space nk: constant indices broadcast (mirrors old ``fixed``)."""
    if isinstance(value, (int, float, complex)):
        return np.full(wl.shape, complex(value), dtype=np.complex128)
    return _eval_nk(value, wl)


def _group_map(groups):
    """Bound groups pass through; param dicts construct natively."""
    from navette._structure import Group as _RsGroup
    return {str(name): (g if isinstance(g, _RsGroup) else _RsGroup(str(name), **dict(g)))
            for name, g in (groups or {}).items()}


def stack_from_layers(layers: Sequence[Tuple[Any, float]],
                      wavelengths, contrast: Mapping[Any, Any],
                      ambient: Tuple[Any, str] = (1.0, "air"),
                      substrate: Tuple[Any, str] = (1.52, "sub"),
                      names: Optional[Sequence[str]] = None,
                      film_flags: Optional[Dict[str, Any]] = None,
                      groups: Optional[Mapping[str, Any]] = None,
                      per_film_flags: Optional[Mapping[str, Dict[str, Any]]] = None,
                      ) -> Tuple[Any, Dict[str, Any]]:
    """Native ``(DesignStack, contrast_map)`` from ``[(material, d_nm)]`` films.

    Materials are ``MaterialSpec`` / mappings / nk arrays / constants
    (plain names are not resolvable here — pass specs). ``names`` gives
    the film material names (default ``film0…``); ``contrast`` maps host
    name *or film index* → seed material. ``ambient``/``substrate`` are
    ``(material, name)``; always fixed, non-hosts. ``film_flags`` are
    per-film ``LayerSpec`` defaults (``optimize``/``needle``/…, and the
    F1.1 ``gradient`` mixture profile - see ``_film_dicts``).
    ``groups`` maps material name → bound ``Group`` (or param dict):
    thickness/nk scaling, roughness and interface policy expand here —
    the silent-drop limitation is gone. Graded films take one of two
    paths, implied by their flags (no separate declaration): a graded
    film with ``optimize=False`` *and* ``needle=False`` expands WITH
    the profile as pinned background (silent; excluded from needle
    candidacy and thickness optimization; merge/cleanup preserve the
    span). Any other graded film homogenizes with a warning (base
    index — pin it to keep the profile). ``per_film_flags`` overrides
    ``film_flags`` per film name (e.g. grade only the substrate).
    Span-aware graded optimization is future work (D2).
    """
    # Thin over native assemble_design: shaping here, expansion in Rust.
    from navette._smatrix import assemble_design as _assemble
    wl = np.ascontiguousarray(np.asarray(wavelengths, dtype=np.float64))
    names = list(names) if names is not None else [f"film{i}" for i in range(len(layers))]
    if len(names) != len(layers):
        raise ValueError("names length must match layers length.")
    films = _film_dicts(layers, names, wl, film_flags, per_film_flags)
    seeds = [(str(k), f"{k}_seed", _eval_nk(v, wl)) for k, v in contrast.items()]
    return _assemble(
        _half_nk(ambient[0], wl), str(ambient[1] if len(ambient) > 1 else "air"),
        _half_nk(substrate[0], wl), str(substrate[1] if len(substrate) > 1 else "sub"),
        films, _group_map(groups), seeds, wl)



# --- F2.4: the multi-environment request -----------------------------------

_ROW_KEYS = ("coherent", "roughness_nm", "rough_type", "inhomogen", "inh_delta",
             "interface", "interface_thickness_nm", "layer_type", "gradient",
             "optimize", "needle")


def _tab(code: str, material: Any, wl: np.ndarray) -> Dict[str, Any]:
    """One library entry: the material, evaluated onto the run's own grid.

    The design path has no materials library of its own — every nk is
    caller-supplied — so a material arrives here as a spec, a mapping, an
    nk array or a constant, and leaves as a `TableMaterial` on exactly the
    wavelengths the run uses. No resampling happens later because there is
    nothing left to resample onto.
    """
    # `_half_nk`, not `_eval_nk`: a constant index is a legal material on
    # this surface (and on the flat one), and broadcasting it here is what
    # makes `(1.52, "sub")` mean the same thing in both.
    nk = _half_nk(material, wl)
    wls = [float(x) for x in wl]
    return {
        "name": code, "code": code, "model": "TableMaterial", "params": {},
        "n_data": {"wavelengths": wls, "values": [float(z.real) for z in nk]},
        "k_data": {"wavelengths": wls, "values": [float(z.imag) for z in nk]},
    }


def _entry_parts(entry, default_name: str):
    """One layer entry -> (material, thickness, name, row flags).

    Three spellings, in order of how much they say: ``(material, d_nm)``,
    ``(material, d_nm, name)``, and a mapping for everything else
    (roughness, interface, ``layer_type`` for an ambient or substrate
    row, a gradient spec).
    """
    if isinstance(entry, Mapping):
        e = dict(entry)
        try:
            material = e.pop("material")
            thickness = e.pop("thickness_nm")
        except KeyError as exc:
            raise ValueError(
                f"layer entry {default_name!r}: a mapping needs 'material' and "
                f"'thickness_nm' (missing {exc.args[0]!r})."
            ) from None
        name = str(e.pop("name", default_name))
        unknown = set(e) - set(_ROW_KEYS)
        if unknown:
            raise ValueError(
                f"layer entry {name!r}: unknown row keys {sorted(unknown)} "
                f"(one of {', '.join(_ROW_KEYS)})."
            )
        return material, thickness, name, e
    seq = tuple(entry)
    if len(seq) == 2:
        return seq[0], seq[1], default_name, {}
    if len(seq) == 3:
        return seq[0], seq[1], str(seq[2]), {}
    raise ValueError(
        f"layer entry {default_name!r}: expected (material, d_nm), "
        f"(material, d_nm, name) or a mapping - got {len(seq)} items."
    )


def _program_material(materials, code: str, where: str):
    """One `material_code` -> something `_tab` can evaluate.

    The SPEC is preferred over the provider's evaluated curve, because the
    provider was built on the grid the program was LOADED with and the run
    has a grid of its own. Handing the spec through means the material is
    evaluated once, on the run's grid, exactly as an inline material is —
    which is what lets a document and a keyword call land on the same bits
    instead of on the same bits plus a resample.
    """
    if materials is None:
        raise ValueError(
            f"{where}: material_code {code!r} cannot be resolved - the "
            "program has no materials section and none was supplied."
        )
    shelf = getattr(materials, "_dict", None)
    if shelf is not None and code in shelf:
        return shelf[code]
    if getattr(materials, "contains", lambda _c: False)(code):
        # A provider without a spec shelf (a weaver, say). Its curve is on
        # its own grid; `_eval_nk` will refuse a length mismatch, which is
        # the right failure - silently interpolating here would be the
        # resample the docstring above exists to avoid.
        return materials.get_nk(code)
    raise ValueError(
        f"{where}: unknown material_code {code!r} - a design row resolves "
        "against the program's materials section, like a structure's layers do."
    )


def _program_rows(rows, materials, where: str):
    """`[LayerRow]` from a document -> `[mapping]` for `_entry_parts`.

    `material_code` is the film NAME as well as the material: on this
    surface a film name is the parameter's identity across environments,
    and a document has no second field to spell it with. Two rows of the
    same code in one run are therefore one parameter named twice, which
    `_environments_request` refuses by name.
    """
    out = []
    for i, row in enumerate(rows or ()):
        if not isinstance(row, Mapping) or "material_code" not in row:
            raise ValueError(
                f"{where}: layer {i} needs a 'material_code' "
                "(a design document's rows are LayerRows)."
            )
        r = dict(row)
        code = str(r.pop("material_code"))
        thickness = r.pop("thickness_nm", 0.0)
        # Keep only the flags this surface carries; a document may spell
        # defaults the request shape leaves implicit.
        flags = {k: v for k, v in r.items() if k in _ROW_KEYS}
        unknown = set(r) - set(_ROW_KEYS)
        if unknown:
            raise ValueError(
                f"{where}: layer {code!r} has unknown row keys "
                f"{sorted(unknown)} (one of {', '.join(_ROW_KEYS)})."
            )
        out.append({"material": _program_material(materials, code, where),
                    "thickness_nm": float(thickness), "name": code, **flags})
    return out


def design_from_program(program, materials=None):
    """``(design, environments)`` kwargs from a loaded program document.

    ``run_needle(design=..., environments=...)`` takes materials, because
    the design path has no library of its own; a program document takes
    material CODES, because it has one. This is the bridge, and it is the
    only one — both halves then go through `_environments_request`, so a
    run described in a file and the same run described in Python are the
    same request.

    ``materials`` overrides the program's own provider (the ``context=``
    case, where the library was supplied rather than parsed).
    """
    mats = materials if materials is not None else getattr(program, "materials", None)
    design_doc = getattr(program, "design", None) or {}
    env_doc = getattr(program, "environments", None) or []
    if not design_doc:
        raise ValueError(
            "this program has no 'design' section - nothing to run as a "
            "multi-environment design."
        )
    design = {
        str(seg): _program_rows(body.get("layers"), mats,
                                f"design segment {seg!r}")
        for seg, body in design_doc.items()
    }
    environments = []
    for env in env_doc:
        name = str(env.get("name", ""))
        stack = []
        for j, part in enumerate(env.get("stack") or ()):
            # The document spells the union with both keys present and one
            # of them null (serde `Option`); the run door spells it with
            # exactly one key. Translate, rather than teaching the door a
            # second spelling.
            ref = part.get("design")
            rows = part.get("layers")
            if ref is not None:
                stack.append({"design": str(ref)})
            elif rows is not None:
                stack.append({"layers": _program_rows(
                    rows, mats, f"environment {name!r} segment {j}")})
            else:
                raise ValueError(
                    f"environment {name!r}: stack entry {j} gives neither "
                    "'design' nor 'layers'."
                )
        environments.append({"name": name, "stack": stack})
    return design, environments


def _environments_request(design, environments, contrast, wl,
                          ambient, substrate, film_flags, per_film_flags,
                          groups) -> Dict[str, Any]:
    """The `DesignRequest` a multi-environment run is, as plain JSON data.

    This is the whole of F2.4's surface: `run_needle(design=...,
    environments=[...])` shapes this document and hands it to the core,
    which compiles it exactly as it compiles the same document loaded from
    a program file. One schema, one compiler, one set of refusals — the
    Python surface cannot drift from the file surface because there is
    nothing here for it to drift with.

    Film NAMES are the cross-environment parameter identity, and they are
    the library codes too: two films of the same physical material are two
    parameters and get two codes carrying identical tables. Unnamed films
    are named after their segment and position (``coat0``, ``coat1``), so
    a design that never mentions a name still has stable ones.
    """
    library: Dict[str, Dict[str, Any]] = {}

    def register(code: str, material) -> str:
        if code in library:
            raise ValueError(
                f"duplicate film name {code!r} - a film name is a design "
                "parameter's identity across environments, so it must be "
                "unique within the run."
            )
        library[code] = _tab(code, material, wl)
        return code

    if not isinstance(design, Mapping) or not design:
        raise ValueError(
            "design must be a non-empty mapping of segment id -> layer list."
        )
    design_doc: Dict[str, Any] = {}
    for seg, layers in design.items():
        rows = []
        for i, entry in enumerate(layers):
            material, d, name, flags = _entry_parts(entry, f"{seg}{i}")
            register(name, material)
            rows.append({"material_code": name, "thickness_nm": float(d), **flags})
        design_doc[str(seg)] = {"layers": rows}

    env_doc = []
    seen_env = set()
    for env in environments:
        if not isinstance(env, Mapping):
            raise ValueError(
                "each environment must be a mapping with 'name' and 'stack'."
            )
        name = str(env.get("name", ""))
        if not name:
            raise ValueError("each environment needs a non-empty 'name'.")
        if name in seen_env:
            raise ValueError(f"duplicate environment name {name!r}.")
        seen_env.add(name)
        unknown = set(env) - {"name", "stack"}
        if unknown:
            raise ValueError(
                f"environment {name!r}: unknown keys {sorted(unknown)} "
                "(only 'name' and 'stack')."
            )
        stack = []
        for j, part in enumerate(env.get("stack") or ()):
            if not isinstance(part, Mapping):
                raise ValueError(
                    f"environment {name!r}: stack entry {j} must be "
                    "{'design': id} or {'layers': [...]}."
                )
            has_design, has_layers = "design" in part, "layers" in part
            if has_design == has_layers:
                raise ValueError(
                    f"environment {name!r}: stack entry {j} must give exactly "
                    "one of 'design' or 'layers'."
                )
            if has_design:
                stack.append({"design": str(part["design"])})
                continue
            rows = []
            for i, entry in enumerate(part["layers"]):
                material, d, fname, flags = _entry_parts(entry, f"{name}~{j}~{i}")
                register(fname, material)
                rows.append({"material_code": fname,
                             "thickness_nm": float(d), **flags})
            stack.append({"layers": rows})
        env_doc.append({"name": name, "stack": stack})

    # Ambient and substrate: half-space rows, one per environment, unless
    # the caller placed their own (`layer_type` 0/2 inside a fixed
    # segment). Absent entirely, the core's n = 1.0 / 1.52 apply.
    for slot, value, ltype, default_name in (
        ("ambient", ambient, 0, "air"), ("substrate", substrate, 2, "sub"),
    ):
        if value is None:
            continue
        material = value[0] if isinstance(value, (tuple, list)) else value
        code = str(value[1]) if isinstance(value, (tuple, list)) and len(value) > 1 \
            else default_name
        library.setdefault(code, _tab(code, material, wl))
        row = {"material_code": code, "thickness_nm": 0.0, "layer_type": ltype}
        for env in env_doc:
            if slot == "ambient":
                env["stack"].insert(0, {"layers": [row]})
            else:
                env["stack"].append({"layers": [row]})

    # Contrast: host film name -> seed material. The seed becomes a
    # library entry of its own, the way it does on the flat path.
    cmap = {}
    for host, seed in (contrast or {}).items():
        host = str(host)
        if host not in library:
            raise ValueError(
                f"contrast: {host!r} is not a film in this run "
                f"(films: {', '.join(sorted(library))})."
            )
        code = f"{host}_seed"
        library.setdefault(code, _tab(code, seed, wl))
        cmap[host] = code

    return {
        "structure": {"label": "design", "layers": [], "groups": list(groups or ())},
        "library": list(library.values()),
        "contrast": cmap,
        "film_flags": dict(film_flags or {}),
        "per_film_flags": {str(k): dict(v) for k, v in (per_film_flags or {}).items()},
        "ambient_name": "air",
        "substrate_name": "sub",
        "design": design_doc,
        "environments": env_doc,
    }


def run_needle(layers: Optional[Sequence[Tuple[Any, float]]] = None,
               targets: Any = None,
               angles_deg=None, wavelengths=None,
               contrast: Optional[Mapping[Any, Any]] = None,
               pipeline_config=None, needle_config=None, lm_config=None,
               callback: Optional[Callable[[int, Dict], None]] = None,
               *,
               design: Optional[Mapping[str, Any]] = None,
               environments: Optional[Sequence[Mapping[str, Any]]] = None,
               **stack_kwargs) -> Dict[str, Any]:
    """Design a coating with the needle pipeline, end to end.

    Parameters
    ----------
    layers : [(material, d_nm)] films (ambient/substrate via ``stack_kwargs``).
        Materials are ``MaterialSpec`` / mappings / nk arrays (constant
        complex also accepted).
    targets : TargetCollection or MeritSpec (native).
    angles_deg / wavelengths : solver grid (degrees / nm).
    contrast : {host name, film index, or "film{i}" name: seed material}.
        Hosts without an entry are never split (empty needle history when
        nothing matches — check names when a run inserts nothing).
    pipeline_config / needle_config / lm_config : native configs (or None).
    callback : ``(macro_cycle, phase_dict)``; raising aborts (USER_ABORT).
    stack_kwargs : ambient, substrate, per-film flags (see
        :func:`stack_from_layers`).

    design : {segment id: [(material, d_nm)]} — named design segments,
        each defined ONCE and shared by every environment that references
        it. Mutually exclusive with ``layers``: one of the two, never both
        and never neither. Entries also accept ``(material, d_nm, name)``
        or a mapping (``material``, ``thickness_nm``, ``name``, plus row
        flags including ``layer_type`` for a half-space row).
    environments : [{"name": ..., "stack": [{"design": id} | {"layers": [...]}]}]
        Each environment's ordered segment list. The shared design films
        are the free variables; everything in a ``layers`` segment is
        fixed surroundings, and the row schema forces ``optimize`` and
        ``needle`` false there (an explicit ``true`` refuses). **Known
        gap (review PB, M1):** a ``per_film_flags`` entry is applied
        after the row and is keyed by material code, so it can re-enable
        either flag on a surrounding — and under K > 1 the LM then moves
        environment 0's surrounding against environment-0-only residuals
        while every other environment keeps the compiled thickness. Do
        not name a surrounding's material in ``per_film_flags`` until the
        assembler refuses it. Requires ``design``; the roster must match
        the one given to ``build_merit_spec(environments=...)`` **in the
        same order** — the core compares counts, not names (M7), so a
        permuted roster is accepted and scores against the wrong
        surroundings. Passing the ``TargetCollection`` instead of a
        pre-built spec avoids that entirely.

    With ``design``/``environments`` the contrast map is keyed by FILM
    NAME (the design parameter's identity), not by index.

    Returns the native result dict (``termination``, ``final_mf``,
    ``phases``, final ``stack``). With environments the returned stack is
    environment 0's assembly, which IS the shared design object.
    """
    if (layers is None) == (design is None):
        raise ValueError(
            "run_needle: give exactly one of `layers` (a flat film list) or "
            "`design` (named segments shared across `environments`) - "
            + ("both were given." if layers is not None else "neither was given.")
        )
    for name, value in (("targets", targets), ("angles_deg", angles_deg),
                        ("wavelengths", wavelengths)):
        if value is None:
            raise ValueError(f"run_needle: `{name}` is required.")
    if design is not None:
        return _run_needle_environments(
            design, environments, targets, angles_deg, wavelengths, contrast,
            pipeline_config, needle_config, lm_config, callback, stack_kwargs)
    if contrast is None:
        raise ValueError("run_needle: `contrast` is required.")
    if environments is not None:
        raise ValueError(
            "run_needle: `environments` needs `design` - a flat film list has "
            "no named segments for an environment to reference."
        )
    # Thin over native run_design: evaluate + key/flag shaping here,
    # assembly + macro-loop in Rust. Contrast-key normalization stays
    # (presentation over the film order).
    from navette._smatrix import run_design as _run_design
    wl = np.ascontiguousarray(np.asarray(wavelengths, dtype=np.float64))
    angs = np.ascontiguousarray(np.asarray(angles_deg, dtype=np.float64))
    names = stack_kwargs.pop("names", None)
    if names is None:
        names = [f"film{i}" for i in range(len(layers))]
    names = list(names)
    if len(names) != len(layers):
        raise ValueError("names length must match layers length.")
    def _host_key(k):
        if isinstance(k, bool):
            return str(k)
        if isinstance(k, int) and 0 <= k < len(names):
            return names[k]
        s = str(k)
        if s in names:
            return s
        if s.startswith("film") and s[4:].isdigit() and int(s[4:]) < len(names):
            return names[int(s[4:])]
        return s
    ambient = stack_kwargs.pop("ambient", (1.0, "air"))
    substrate = stack_kwargs.pop("substrate", (1.52, "sub"))
    film_flags = stack_kwargs.pop("film_flags", None)
    groups = stack_kwargs.pop("groups", None)
    per_film_flags = stack_kwargs.pop("per_film_flags", None)
    if stack_kwargs:
        raise TypeError(f"run_needle: unknown stack options {sorted(stack_kwargs)}.")
    films = _film_dicts(layers, names, wl, film_flags, per_film_flags)
    seeds = [(_host_key(k), f"{_host_key(k)}_seed", _eval_nk(v, wl))
             for k, v in contrast.items()]
    spec = (_build_merit_spec(targets) if isinstance(targets, TargetCollection)
            else targets)
    amb_nk = _half_nk(ambient[0], wl)
    sub_nk = _half_nk(substrate[0], wl)
    gmap = _group_map(groups)
    return _run_design(
        amb_nk, str(ambient[1] if len(ambient) > 1 else "air"),
        sub_nk, str(substrate[1] if len(substrate) > 1 else "sub"),
        films, gmap, seeds, wl, angs, spec,
        pipeline_config=pipeline_config, needle_config=needle_config,
        lm=lm_config, callback=callback)


def _run_needle_environments(design, environments, targets, angles_deg,
                             wavelengths, contrast, pipeline_config,
                             needle_config, lm_config, callback,
                             stack_kwargs) -> Dict[str, Any]:
    """The `design=`/`environments=` arm of :func:`run_needle` (F2.4).

    Shapes the request document and hands it to the core in one call. The
    flat arm above is untouched by this: it still goes through
    ``run_design`` with evaluated arrays, which is what keeps old calls on
    the bitwise path rather than on a re-implementation of it.
    """
    import json as _json
    from navette._smatrix import run_design_environments as _run_envs
    wl = np.ascontiguousarray(np.asarray(wavelengths, dtype=np.float64))
    angs = np.ascontiguousarray(np.asarray(angles_deg, dtype=np.float64))
    if not environments:
        raise ValueError(
            "run_needle: `design` needs `environments` - a design segment is "
            "only reachable through an environment that references it."
        )
    ambient = stack_kwargs.pop("ambient", None)
    substrate = stack_kwargs.pop("substrate", None)
    film_flags = stack_kwargs.pop("film_flags", None)
    groups = stack_kwargs.pop("groups", None)
    per_film_flags = stack_kwargs.pop("per_film_flags", None)
    if stack_kwargs:
        raise TypeError(f"run_needle: unknown stack options {sorted(stack_kwargs)}.")
    request = _environments_request(
        design, environments, contrast, wl, ambient, substrate,
        film_flags, per_film_flags, groups)
    roster = [e["name"] for e in request["environments"]]
    spec = (_build_merit_spec(targets, environments=roster)
            if isinstance(targets, TargetCollection) else targets)
    return _run_envs(
        _json.dumps(request), wl, angs, spec,
        pipeline_config=pipeline_config, needle_config=needle_config,
        lm=lm_config, callback=callback)
