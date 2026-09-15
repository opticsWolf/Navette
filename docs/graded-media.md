# Graded media: spans, expansion, and the two things that bite

A mixture gradient is a layer **property**, not a hand-built stack of
sublayers. One span names its two endpoint materials, a mixing rule, and a
profile mode; the expansion into sublayers happens inside, at assembly.

```python
per_film_flags = {"grad": {"gradient": {
    "material_b": nk_high,     # the inclusion (material_a is the film itself)
    "f_start": 0.0, "f_end": 1.0,   # FixedSpan: linear across the span
    "ema": "Bruggeman",
}}}
```

Two profile modes, and exactly one slope spelling is accepted (giving both
is refused rather than silently preferred):

| mode | keys | meaning |
|---|---|---|
| `FixedSpan` | `f_start`, `f_end` | composition runs linearly across the span |
| `RateCapped` | `f_start`, `rate`, `ref_thickness`, `f_min`, `f_max` | `f(z) = f_start + rate·z/ref_thickness`, clamped — a thick span saturates into a pure-material tail |

Single-material drift (`inhomogen` + `inh_delta`) is the same idea with one
material and has the same two modes.

## Spans are authored; rows are solved

This is the number that matters and the one that is easy to miss. A span
expands into sublayers at assembly, and everything downstream — the
solver, the needle operator, a timing — sees **rows**.

`examples/rugate_gradient_vs_discrete.py` measures it: 16 authored spans
reach the solver as 80 rows (5 sublayers each), against 16 rows for the
stepped equivalent. Same optical thickness per period, same demands,
5× the rows. The graded design wins that comparison on merit — no abrupt
interfaces means no sidelobes where the demand asks for transmission,
which is the whole reason a rugate is bought — but it wins it at 5× the
solve cost, and the cost is linear in rows.

Two consequences worth stating plainly:

- **In a multi-environment run the row cost is ×K.** Each environment
  assembles and solves its own stack, so a 64-sublayer cover carried by
  three environments is 192 extra rows *per evaluation*. Where the
  surroundings are thick and passive, an S-matrix embedding (solve them
  once, reuse) would remove most of that; it is a known follow-up, not
  part of v1. See the environment section in
  [spectralweave-target-kinds.md](spectralweave-target-kinds.md).
- **A span's thickness can be an optimizer parameter.** The span is the
  parameter, not its sublayers — which is also why the sublayers are
  exempt from the thin-layer sweep (below).

## Névot-Croce on a graded span: read this before setting rtype 5

Interface roughness type 5 applies a Névot-Croce factor at **every
boundary it is set on**. A graded span's sublayer boundaries are
boundaries. Switching rtype 5 on across a 16-span rugate therefore means
79 interface factors, not the handful the author had in mind.

That is not merely expensive, it is outside the model. Névot-Croce is a
small-perturbation result, valid for

> σ ≪ λ **and** σ ≪ the layer thickness

and a sublayer is thin *by construction* — the span's thickness divided by
its sublayer count — so the second condition is the one that fails first.
A 78 nm span at 5 sublayers gives 16 nm rows; a σ of 3 nm is already a
fifth of the row.

The rule that follows: **set interface roughness on the spans you meant,
not on the expansion.** Roughness belongs on the real boundaries of the
coating — substrate, cover, the top surface — and a graded span's internal
boundaries are an artifact of the discretization, not features of the
film. Putting a Névot-Croce factor on them models a physical roughness
that is not there and does it with an approximation that does not hold.

## The thin-gradient floor

Sublayer thickness is the span's thickness divided by its sublayer count,
so a thin span has thin sublayers. At default settings a **3 nm gradient
is three 1 nm rows**, below the usual `clamp_min_nm` of 2 nm.

Those rows are **not deleted**. The thin-layer sweep filters on
`optimize`, and a sublayer is not an independent parameter — the span is —
so the expansion is exempt. (Before F0.1 it was not, and `clamp_all`
deleted pinned sublayers; that is fixed.)

What the exemption does *not* do is make a 3 nm gradient meaningful. Three
steps do not resolve a profile: the staircase approximation needs enough
rows that the index change per row is small compared to what a wave at λ
can distinguish, and at three rows it is not an approximation of a
gradient so much as a three-layer stack with unusual indices. If a
gradient has to be that thin, it is a **mixed layer** — author it as one,
with `inhomogen` if it drifts at all, and say so.

Rule of thumb: below roughly 6 nm at default sublayer counts you are at
the floor, and the honest description of the film is a mixture rather than
a gradient.

## See also

- `examples/rugate_gradient_vs_discrete.py` — the row-count and merit
  comparison above, runnable.
- `examples/multi_environment_ar.py` — one design, several surroundings,
  where the ×K row cost applies.
- [spectralweave-target-kinds.md](spectralweave-target-kinds.md) — the
  demand side, including the environment tag and the ×K cost note.
