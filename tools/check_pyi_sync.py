# SPDX-License-Identifier: LGPL-3.0-or-later
"""CI sync guard: the `.pyi` stubs vs the PyO3 registration lists.

Two passes, because they catch different rot:

1. **Surface** — parse each `#[pymodule]` body in `rust/navette-py/src/` for
   its `wrap_pyfunction!` / `add_class::<>` / `add(...)` entries and diff the
   names against the stub's top-level `def`s, `class`es and annotated
   assignments. Needs no build; this is the §24.1 audit, automated.
2. **Signatures** — if the extension imports, compare each function's and
   method's parameter *names* (from PyO3's `__text_signature__`) against the
   stub's. A renamed keyword argument is invisible to pass 1 and breaks every
   caller that used it.

Stdlib only. Exit 1 on any drift, loudly.
"""
from __future__ import annotations

import ast
import inspect
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BINDINGS = ROOT / "rust" / "navette-py" / "src"
STUBS = ROOT / "src" / "navette"

#: Submodules with no stub yet. Adding a module without a stub should fail,
#: so this list is explicit rather than inferred from which files exist.
UNSTUBBED = {"_color", "_interpolate", "_structure"}

# ---------------------------------------------------------------------------
# Pass 1: the registration list
# ---------------------------------------------------------------------------

_PYMODULE = re.compile(
    r"#\[pymodule\]\s*(?:pub\s+)?fn\s+(_\w+)\s*(?:<[^>]*>)?\s*\([^)]*\)\s*"
    r"->\s*PyResult<\(\)>\s*\{",
)
_FUNC = re.compile(r"wrap_pyfunction!\s*\(\s*(?:[\w:]*::)?(\w+)\s*,")
_CLASS = re.compile(r"add_class::<\s*(?:[\w:]*::)?(\w+)\s*>")
_CONST = re.compile(r'\bm\.add\s*\(\s*"(\w+)"')
#: `#[pyclass(name = "Foo")]` — the Rust type is `PyFoo`, Python sees `Foo`.
_PYCLASS_NAME = re.compile(r'#\[pyclass\([^)]*name\s*=\s*"(\w+)"')
#: The `struct`/`enum` the attribute lands on, possibly past derives and docs.
_NEXT_ITEM = re.compile(r"\b(?:struct|enum)\s+(\w+)")


def _brace_body(text: str, open_idx: int) -> str:
    """Source between the brace at `open_idx` and its match."""
    depth, i = 0, open_idx
    while i < len(text):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[open_idx + 1 : i]
        i += 1
    raise ValueError("unbalanced braces in a #[pymodule] body")


def rust_surface() -> tuple[dict[str, set[str]], dict[str, str]]:
    """`{module: {exported names}}` plus `{RustType: PythonName}`."""
    class_names: dict[str, str] = {}
    for src in sorted(BINDINGS.glob("*.rs")):
        text = src.read_text(encoding="utf-8")
        for m in _PYCLASS_NAME.finditer(text):
            item = _NEXT_ITEM.search(text, m.end())
            if item:
                class_names[item.group(1)] = m.group(1)

    modules: dict[str, set[str]] = {}
    for src in sorted(BINDINGS.glob("*.rs")):
        text = src.read_text(encoding="utf-8")
        for m in _PYMODULE.finditer(text):
            body = _brace_body(text, m.end() - 1)
            names = set(_FUNC.findall(body)) | set(_CONST.findall(body))
            for rust_ty in _CLASS.findall(body):
                names.add(class_names.get(rust_ty, rust_ty))
            modules[m.group(1)] = names
    return modules, class_names


# ---------------------------------------------------------------------------
# Stub parsing
# ---------------------------------------------------------------------------


def stub_surface(path: Path) -> tuple[set[str], dict[str, list[str]]]:
    """`({top-level names}, {qualified name: [parameter names]})`."""
    tree = ast.parse(path.read_text(encoding="utf-8"), str(path))
    names: set[str] = set()
    params: dict[str, list[str]] = {}

    def sig(fn: ast.FunctionDef, drop_self: bool) -> list[str]:
        a = fn.args
        got = [p.arg for p in (*a.posonlyargs, *a.args)]
        if drop_self and got and got[0] in ("self", "cls"):
            got = got[1:]
        return got + [p.arg for p in a.kwonlyargs]

    for node in tree.body:
        if isinstance(node, ast.FunctionDef):
            names.add(node.name)
            params[node.name] = sig(node, False)
        elif isinstance(node, ast.ClassDef):
            names.add(node.name)
            for sub in node.body:
                if not isinstance(sub, ast.FunctionDef):
                    continue
                decos = {
                    d.id if isinstance(d, ast.Name) else getattr(d, "attr", "")
                    for d in sub.decorator_list
                }
                if "property" in decos or "setter" in decos:
                    continue  # a getset descriptor has no comparable signature
                params[f"{node.name}.{sub.name}"] = sig(sub, "staticmethod" not in decos)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
        elif isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name):
                    names.add(t.id)
    return names, params


# ---------------------------------------------------------------------------
# Pass 2: live signatures
# ---------------------------------------------------------------------------


def live_params(obj) -> list[str] | None:
    """Parameter names from PyO3's `__text_signature__`, or None."""
    ts = getattr(obj, "__text_signature__", None)
    if not ts:
        return None
    try:
        sig = inspect.Signature.from_callable(obj)
    except (ValueError, TypeError):
        return None
    return [p for p in sig.parameters if p not in ("self", "cls")]


def live_surface(mod) -> dict[str, list[str]]:
    """`{qualified name: [parameter names]}` for one imported submodule."""
    out: dict[str, list[str]] = {}
    for name in dir(mod):
        if name.startswith("__"):
            continue
        obj = getattr(mod, name)
        if isinstance(obj, type):
            got = live_params(obj)
            if got is not None:
                out[f"{name}.__init__"] = got
            for attr in dir(obj):
                if attr.startswith("__"):
                    continue
                got = live_params(getattr(obj, attr))
                if got is not None:
                    out[f"{name}.{attr}"] = got
        elif callable(obj):
            got = live_params(obj)
            if got is not None:
                out[name] = got
    return out


# ---------------------------------------------------------------------------

def main() -> int:
    modules, _ = rust_surface()
    if not modules:
        print("FAILED: no #[pymodule] found under rust/navette-py/src/")
        return 1

    drift = 0
    checked = 0
    for mod_name in sorted(modules):
        if mod_name == "_navette":
            continue  # the aggregator: it registers the submodules, not names
        stub = STUBS / f"{mod_name}.pyi"
        if not stub.exists():
            if mod_name in UNSTUBBED:
                print(f"note  {mod_name}: no stub (known gap)")
                continue
            print(f"DRIFT {mod_name}: registered but has no {stub.name}")
            drift += 1
            continue
        if mod_name in UNSTUBBED:
            print(f"DRIFT {mod_name}: has a stub now — drop it from UNSTUBBED")
            drift += 1

        checked += 1
        want = modules[mod_name]
        got, stub_params = stub_surface(stub)
        # Aliases and type variables the stub defines for its own use are not
        # part of the surface; only names the module exports must match.
        for missing in sorted(want - got):
            print(f"DRIFT {mod_name}.{missing}: registered in Rust, absent from {stub.name}")
            drift += 1

        # --- pass 2 ---
        try:
            mod = getattr(__import__("navette._navette", fromlist=[mod_name]), mod_name)
        except Exception as exc:  # pragma: no cover - unbuilt tree
            print(f"note  {mod_name}: extension not importable ({exc.__class__.__name__});"
                  " signature pass skipped")
            continue
        live = live_surface(mod)
        for qual, live_names in sorted(live.items()):
            if qual not in stub_params:
                continue  # covered by the surface pass, or a property
            if stub_params[qual] != live_names:
                print(f"DRIFT {mod_name}.{qual}:")
                print(f"        native {live_names}")
                print(f"        stub   {stub_params[qual]}")
                drift += 1

    if drift:
        print(f"pyi sync FAILED ({drift} finding(s))")
        return 1
    print(f"pyi sync OK ({checked} stubbed module(s), {len(UNSTUBBED)} known gap(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
