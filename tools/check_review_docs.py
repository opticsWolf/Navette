# SPDX-License-Identifier: LGPL-3.0-or-later
"""CI guard: `validation/review/` vs the table in its README.

The battery runs every harness in that directory, but the README is the only
map of what each one proves and which review or plan item asked for it. Two
of them (`lm_check.py`, `color_grad_python.py`) sat in the battery undocumented
for several versions, because the table was written for `code_review.md`
sections 19-23 and nobody widened it when the R4 plan items added harnesses of
their own. Nothing failed; the map was just quietly incomplete.

Three passes:

1. **Coverage** — every `*.py` in `validation/review/` has a table row, and
   every row names a file that exists. A harness worth running is worth a
   sentence saying what it proves.
2. **Provenance** — each row's middle cell names the review section or plan
   item that asked for the check. A harness with no origin is one nobody can
   decide to retire.
3. **The README's own claim** — "The `loom` reference is never imported."
   These harnesses exist precisely because the parity suite compares the
   engine against a port of itself; one `import` away and they would be doing
   the same thing. Checked on the AST rather than by grep, so the word may
   still appear in prose (it does, in two docstrings explaining exactly this).

Stdlib only. Exit 1 on any finding, loudly.
"""
from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REVIEW = ROOT / "validation" / "review"
README = REVIEW / "README.md"

_HEADER = "| Script | Review section | What it verifies |"
#: A row's first cell: `| `name.py` | ...`
_ROW = re.compile(r"^\|\s*`([^`]+)`\s*\|([^|]*)\|")
#: Anything that reads as a review section or plan item: 19, 19.2, R4.6, C7.
_ORIGIN = re.compile(r"(§\s*\d|\bR\d+\.\d+|\bC\d+\b)")


def table_rows(text: str) -> dict[str, str]:
    """Map script name -> the row's 'Review section' cell."""
    lines = text.splitlines()
    try:
        start = next(i for i, ln in enumerate(lines) if ln.strip() == _HEADER)
    except StopIteration:
        print(f"MISSING the harness table header in {README.name}")
        print(f"        expected a line reading exactly: {_HEADER}")
        sys.exit(1)

    rows: dict[str, str] = {}
    for ln in lines[start + 2:]:          # skip the header and its |---| rule
        if not ln.startswith("|"):
            break                          # the table ends at the first non-row
        m = _ROW.match(ln)
        if m:
            rows[m.group(1)] = m.group(2).strip()
    return rows


def imports_loom(path: Path) -> bool:
    """Does this harness actually import the parity reference?"""
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except SyntaxError as exc:            # a broken harness is its own failure
        print(f"UNPARSEABLE {path.name}: {exc}")
        return False
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            if any("loom" in a.name for a in node.names):
                return True
        elif isinstance(node, ast.ImportFrom):
            if "loom" in (node.module or "") or any(
                    "loom" in a.name for a in node.names):
                return True
    return False


def main() -> int:
    rows = table_rows(README.read_text(encoding="utf-8"))
    scripts = {p.name: p for p in sorted(REVIEW.glob("*.py"))
               if not p.name.startswith("_")}

    findings = 0

    for name in sorted(set(scripts) - set(rows)):
        print(f"UNDOCUMENTED validation/review/{name}: the battery runs it, "
              f"{README.name} does not list it")
        findings += 1

    for name in sorted(set(rows) - set(scripts)):
        print(f"STALE {README.name} lists `{name}`, which is not in "
              f"validation/review/")
        findings += 1

    for name in sorted(set(rows) & set(scripts)):
        if not _ORIGIN.search(rows[name]):
            print(f"NO ORIGIN `{name}`: its 'Review section' cell "
                  f"({rows[name]!r}) names no section or plan item")
            findings += 1

    for name, path in sorted(scripts.items()):
        if imports_loom(path):
            print(f"LOOM {name} imports the parity reference; these harnesses "
                  f"exist to be independent of it")
            findings += 1

    if findings:
        print(f"review docs FAILED ({findings} finding(s))")
        return 1
    print(f"review docs OK ({len(scripts)} harness(es), all documented "
          f"and loom-free)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
