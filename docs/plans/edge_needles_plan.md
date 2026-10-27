# Boundary needles — built for flat runs; joint runs to do

Status: `NeedleCycleConfig(boundary_needles=True)` is built for flat
(single-environment) runs. A multi-environment run refuses the option.
§B5 is the remaining step; §B7 records the agreed follow-up work and
its order.

This replaces the earlier "edge needles" plan (a site at the top of the
outer film and at the bottom of the inner film only): every span
boundary is a site, and the two outer edges are two of them.

## B1. What it adds

Interior sites are `k * scan_step_nm` for `k = 1 .. floor(d / step) - 1`
inside each host. The range leaves out both film boundaries, so a film
thinner than `2 * scan_step_nm` has no site at all, a fresh 5 nm needle
has one, and the outer edges can never change by insertion.

With the option on, `build_boundary_sites` (`synthesis/needle_pass.rs`)
adds one site per span boundary, the ambient and substrate sides
included, wherever at least one neighbouring span is a needle host
(singleton-bulk, `needle` set). A boundary touching a zero-thickness span
is skipped. A needle chosen there is a new film between the two
neighbours; nothing is split
(`DesignStack::insert_film_at_boundary`). Boundary sites follow the
interior ones in the site list, so an exact tie goes to the interior site.

## B2. Why one site per boundary, and no new kernel

The needle kernel adds a slab without thinning the host, and the
tangential fields are continuous across an interface, so P has one value
on a boundary: the bottom of film i and the top of film i+1 are the same
site. The scan's own depth locator puts an exact boundary depth in the
row below at depth 0 (the substrate edge: in the last row, at its
thickness), so boundary sites run through the existing kernels.
`build_boundary_sites` sums the rows in the locator's order, so the depth
lands exactly. `p_is_continuous_across_a_boundary` pins the claim:
approached from either film, P reaches the boundary value.

Roughness sits on each row's top interface (`sigma = rv[i + 1]`). On a
rough boundary the locator's side is below the rough interface, while
the inserted film takes no roughness and the interface stays on the row
below it. So the prediction there is first order in a slightly different
model; the re-optimization sees the real stack.

## B3. What a boundary site may insert

`MaterialPool::boundary_accepts` (`synthesis/cycle.rs`):

- **Index:** the needle differs from both neighbours, the adjacent row's
  index and a layer's own bulk-row index. The ambient and the substrate
  count. A needle of a neighbour's index is that neighbour grown thicker,
  which the optimizer already covers (and the layers would merge).
- **Insertion rules:** every neighbour that is a needle host must accept
  the needle. A closed host (`[]`) closes its boundaries too. Non-host
  neighbours (`needle = false`, graded spans, ambient, substrate) only
  exclude by index.
- **Edge lists:** on the ambient side the needle becomes the edge layer,
  so it must be on `outer` when that is set; on the substrate side, on
  `inner`. `insert_film_at_boundary` checks the same at the door.

## B4. Cost

The sweep builds the stack fields once per (material, wavelength, angle,
polarization), at a cost linear in the row count. That part is
unchanged. Each site then costs a small fixed amount per spectral point.
A design with L films gains at most L + 1 sites against roughly
total thickness / step interior ones: about +3 % of the scan for a
14-film, 1134 nm design at the 2 nm default, more for a stack of many
thin films. The scan is a few percent of an insertion, the
re-optimization the rest. To be measured with `bench_needle_insertion`.

## B5. Remaining: joint (multi-environment) runs

The joint sweep routes each site to a shared design parameter (slot) and
a step ordinal, and `CompiledEnvironments::insert_seed` splits one slot in
every environment. A boundary needle needs:

- a boundary locus key: "before slot s" inside a segment, plus the
  segment's outer edges, agreed on by every environment;
- `CompiledEnvironments::insert_film_at_boundary`: insert a new slot in
  every environment's template (routing shifts by one, not two), with the
  design stack's own insertion kept in step (the alignment check catches
  any drift);
- on a segment edge, the neighbour outside the design differs per
  environment and is no design parameter, so it brings no insertion
  rules. The needle must differ from the inside neighbour and be on the
  segment's edge list. The outside neighbour excludes it only when the
  needle has its index in every environment (it would do nothing
  anywhere); an insertion that matters in at least one environment is
  allowed. §B7 lists the cases.

Tests to add then: over two identical environments a boundary needle
lands where the single-environment run lands it; the alignment check
holds after a boundary insertion; the refusal test turns into a run.

## B6. Tests in place

- Core: `boundary_sites_sit_on_every_host_boundary`,
  `p_is_continuous_across_a_boundary`,
  `a_boundary_needle_differs_from_both_neighbours`,
  `a_boundary_needle_needs_both_hosts_rules`,
  `a_boundary_needle_on_an_edge_is_on_its_list`,
  `boundary_needles_reach_films_too_thin_to_split`,
  `a_boundary_film_splits_nothing`, `a_boundary_film_obeys_the_edge_lists`.
- Python: `validation/regression/synthesis/test_boundary_needles.py`
  (off is bit-identical, a pinned 3 nm film reached on its boundaries,
  edge lists, rules, the environments refusal).
- Off by default: every existing test and the needle pin are unchanged.

## B7. Agreed next steps (2026-10-27)

Decisions from the review of edge insertion, `outer`/`inner` and
multi-environment runs. Order of work at the end.

### Roughness

- **Interior split stays.** Both halves of a split host keep the host's
  roughness. That is the intended behaviour, not a defect.
- **Boundary insertion keeps the interface with the film below.** The
  existing rough interface keeps its sigma and moves down with the lower
  film; the needle's top interface takes the needle material's own
  roughness (0 for pool materials). Same pattern as the interior split.
- **Roughness cap (missing today).** No code caps sigma against
  thickness: `Layer::property_issues` refuses only negative and
  non-finite values, and the solver reads `rv[i + 1]` as stored. Store
  sigma as given; at solve entry build one effective array,
  `sigma_eff = min(sigma, d / c)` over the two films that touch each
  interface (half-spaces do not cap, `c` configurable, default 2), and
  have every solve path read it (solver, needle kernels, coherent
  blocks). Designs without films thin against their roughness do not
  move.
- **Rough boundary sites.** Névot-Croce depends on both media at an
  interface, so `A | N | B` never reduces to the rough `A | B` as the
  needle thins: the first-order prediction misses a jump about the size
  of a 1 nm needle's own effect (scratch probe, sigma 3-10 nm). Score
  these sites, and only these, by a merit evaluation at a start
  thickness `d0 = max(min thickness, c * sigma)` and insert at `d0`.

### Flat runs

- **F1** Report edge films the clamp holds at the floor (`edges_held` in
  `ClampReport`); optionally let a held edge go when the film behind it
  is also on the list.
- **F2** Skip boundary sites that can accept nothing (an edge list that
  holds only the current edge material, which the index rule excludes).
- **F3** Match edge-list materials by label (the `#n` suffix stripped),
  nk equality as the fallback: bit-exact nk fails under scaled or
  re-evaluated spectra.
- **F4** Optional `seed_edges`: a starting design off its edge lists
  gets a floor-thickness film of the first listed material at the
  boundary (`insert_film_at_boundary`), instead of an error.
- **F5** Document for later: graded and interface-sliced edge films can
  never be on an edge list (`edge_allows` wants a singleton bulk).

### Environments

- **J1 Edge lists per segment.** `outer` faces the incoming light,
  `inner` faces away; the names stay. The joint check today takes the
  first and last design slots of the whole env-0 assembly, which mixes
  segments once there is more than one. Lists become per segment,
  checked on each segment's own top and bottom.
- **Segments reuse the Architect.** A segment reference gains `inverted`
  and `repeat`, with the Architect block's meaning. Each environment's
  entries come from the Architect's block composition
  (`Architect::iter_entries`) and its expansion, so inverted copies move
  roughness with the interface (`owner_of`), reverse grading and carve
  interface slices the way blocks already do. `CompiledEnvironments`
  keeps the routing, `insert_seed`, the thickness copy and the alignment
  check; the Architect lives at compile time only (its handles are
  `Rc<RefCell>`).
  - `DesignStack::from_design` asserts `inv = false` because an inverted
    slice reaches into the previous span. Each inverted copy expands as
    its own unit, and the per-span refresh respects that.
  - Edge lists stay with the segment: inverting a copy inverts its
    `outer`/`inner` with it.
  - Repeated use relaxes the exact-once rule to at-least-once; routing
    becomes one slot to a list of spans per environment.
- **Needles land in every environment.** One insertion edits the shared
  design and appears at every place the segment is used; the gain is
  summed over all of them. In an inverted copy depth `z` becomes
  `d - z`, the split lands in reverse order, and "before slot s" becomes
  "after slot s".
- **Segment-edge cases** (N needle, I inside design film, O_k the
  medium outside in environment k):

  | Case | Rule |
  |---|---|
  | O is the same half-space in every environment | exclude N of O's index |
  | O is a fixed film | allow (N thickens it there only) |
  | O differs per environment | allow |
  | O has N's index in every environment | exclude |
  | two segments touch in some environment | score both owners; the other segment's film counts as O |
  | O a half-space here, a film there | per environment, as above |
  | repeated or inverted copies | over every placement |
  | O incoherent | read the fields at the bottom of I, not in O |

  Always: N differs from I, I's insertion rules apply, the segment's
  edge list applies. Fixed films are never hosts.
- **No fusing across segments.** Touching films of the same material in
  neighbouring segments stay two layers unless the user collapses the
  segments. Joint runs refuse cleanup today, so nothing merges; if joint
  cleanup is allowed later, merging stops at segment and slot
  boundaries.
- If joint removal is ever allowed, `removal_keeps_edges` checks segment
  edges, not assembly edges.

### Order

1. Roughness cap. Built: `optics_core::effective_roughness`, applied at
   the `Solver` constructor and in `DesignStack::solver_arrays`;
   `roughness_cap` on `ScatterMatrix`, `Solver` and `DesignStack`.
2. Start-thickness scoring for rough boundary sites.
3. F1-F4, with F5 in the docs.
4. J1: per-segment edge lists.
5. Segment `inverted` and `repeat` through the Architect.
6. Joint boundary needles (§B5 with the cases above).
