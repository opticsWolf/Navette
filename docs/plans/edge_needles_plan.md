# Boundary needles — built for flat runs; joint runs to do

Status: `NeedleCycleConfig(boundary_needles=True)` is built for flat
(single-environment) runs. A multi-environment run refuses the option.
§B5 is the remaining step.

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
  environment and is no design parameter, so only the edge lists and the
  inside neighbour apply; no index exclusion against the surrounding.

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
