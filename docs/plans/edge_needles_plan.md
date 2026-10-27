# Edge needles — plan (not built)

Status: proposed as a later step after the needle rules (`insert_rules`,
`outer`, `inner`). Nothing here is implemented.

## E0. What it adds

Today every scan site is strictly inside a film: `build_scan_sites`
(`synthesis/needle_pass.rs`) emits `k * scan_step_nm` for `k >= 1` and
stops before the film's far boundary. A needle therefore always lands
between two halves of a host, and the edge films of a design (the one
facing the ambient, the one on the substrate) can never change by
insertion. The needle rules rely on that: `outer` / `inner` are only
guarded where removal could expose a new edge (cleanup, the clamp).

Edge needles add two scan sites per run:

- the **top of the outer film**, the boundary between the ambient and
  the first film;
- the **bottom of the inner film**, the boundary between the last film
  and the substrate.

A needle chosen at one of these sites is not a split. It creates a new
film of the needle material, outside the current edge film, which then
becomes the new edge.

## E1. Why the method supports it

The needle function P(z) is the merit change per unit thickness of an
infinitesimal layer at z. At an outer boundary it is as well defined as
inside a film: the fields on both sides are known from the same solve,
and the perturbation is a thin slab between two half-spaces' worth of
known fields. No new physics, only two more evaluation points in the
existing sweep. The value at the boundary is the limit of the interior
value from the film side when the needle sits on the film's side, which
gives a direct test (E5).

## E2. Scope and switch

- One option, off by default: `edge_needles: bool` on the needle cycle
  config (names to settle with the rest of the config). Off means the
  scan grid, every candidate, and the run are bit-identical to today.
- Which edges: both, when on. A side whose edge film is pinned
  (`needle = false` on that film and the run's policy says the edge is
  fixed) gets no site. Whether `needle = false` on the edge film alone
  should close the edge site is an open question (E6).
- Sites go at the end of the site list (outer first, inner second), so
  the interior sites keep their indices and tie-breaking among interior
  sites is unchanged.

## E3. What a site may insert

- The candidate materials at an edge site are the pool, filtered by:
  - `outer` (for the top site) or `inner` (for the bottom site), when
    set: the new film becomes the edge, so it must be on that list;
  - the index rule: a material equal in index to the current edge film
    is skipped (it would merge into it and is a thickness change, which
    the optimizer already covers).
- `insert_rules` do not apply as written: there is no host on the far
  side. Two choices:
  1. edge sites ignore `insert_rules` and only `outer` / `inner` filter
     (simplest, and the edge lists already say what may go there);
  2. a reserved key per side (`"^"` ambient, `"$"` substrate) in
     `insert_rules` names the materials allowed there.
  Proposal: 1. Revisit if a real design needs edge rules that differ
  from the edge lists.
- After an edge insertion the old edge film is interior. The edge check
  (`EdgeMaterials::check`) runs again on the new stack; by construction
  it passes.

## E4. Mechanics to touch

- `needle_pass.rs`: `ScanSite` needs a kind (interior, top, bottom), or
  a sentinel `film_idx`; `build_scan_sites` emits the two edge sites
  when the option is on; the sweep evaluates P there.
- `structure.rs`: an insertion entry point that prepends or appends a
  film (span bookkeeping: a new singleton-bulk span at the front or the
  back; every later span index shifts by one when prepending).
- `cycle.rs`: `MaterialPool::sites_for` learns the edge filter; the
  candidate record carries the site kind so the report says "outer
  edge" rather than a film index.
- Multi-environment (`environments.rs`): the edges are the shared
  design's. Each environment's assembly sees the design's top and bottom
  at different depths; the edge site is routed to the design segment by
  parameter, like any interior design site. An edge that touches an
  environment's own surrounding layer (not the ambient) is still the
  design's edge, so the same lists apply.
- Reports: the insertion record names the site kind; `needle_results`
  gains nothing else.

## E5. Tests

- Off by default: the needle pin digest and every existing synthesis
  test unchanged (bit identity, not tolerance).
- P at the edge: equals the limit of interior P at the boundary from the
  film side (finite-difference check, a few points approaching it).
- Merit check: inserting a thin edge film of thickness s changes the
  merit by about P(edge) * s, as for interior needles.
- Edge lists: with `outer = [A]` and an edge site that would win with
  B, the run inserts A or nothing at the top, never B.
- Index skip: a material of the edge film's index is never a candidate
  at that edge.
- Multi-environment: over two identical environments, an edge needle
  lands where a single environment's run lands it.
- Span bookkeeping: prepend and append keep `spans` consistent
  (`DesignStack` invariants), and graded edge films are not split
  (they are not hosts here; the new film goes outside them).

## E6. Open questions

- Should `needle = false` on an edge film close the edge site next to
  it? The film flag says "do not split me"; adding a neighbour outside
  it is a different act. Proposal: no, the site belongs to the boundary,
  and only the option and the edge lists decide.
- Does an edge needle on a design with no `outer` / `inner` set need any
  extra guard? Without lists any pool material may become the edge,
  which matches how removal already behaves.
- Incoherent or thick edge films: the site is still well defined, but
  check the scan's coherent assumptions before enabling there.
