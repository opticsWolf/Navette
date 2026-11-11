# Synthesis: open items, inconsistencies and decisions

Status 2026-11-11, branch `needle-material-pool`, after W6 (first cut) and
the integral phase fold. Sources: `synthesis_targets_plan.md` (W1-W8, §4-§7),
the W6 work and a read of the code paths it touched. Each item says where
it lives, what is wrong or missing, and what would close it.

Severity: **bug** = a wrong answer today; **gap** = a refusal or a missing
feature, said out loud; **debt** = code, docs or pins out of step.

---

## 1. Bugs and inconsistencies found

### 1.1 `p_function_multiblock` misaligns its output under a host mask (bug)

`needle_operator.rs`, `locate_hosts_multiblock` + `p_function_multiblock`.
The locator `continue`s past a masked z without pushing an entry, so
`locs` is shorter than `z_grid`. `p_multiblock_point` returns one value
per `loc`, and `p_function_multiblock` adds `contrib[zi]` into `p_out[zi]`:
every z after the first masked one gets the P of a later site, and the
tail stays 0.

- Failure: `z_grid = [z_a (masked), z_b, z_c]` returns `[P(z_b), P(z_c), 0]`.
- Reach: the public Rust function only. It is not bound to Python, and
  synthesis does not call it (W6 uses its own site-to-host map in
  `needle_pass.rs`). Its tests pass no mask.
- Fix: carry the z index in `locs` (or push a sentinel) and scatter by it,
  and add a masked case to the multiblock FD test.

### 1.2 Linux needle-pin digest predates W1 (debt)

`validation/regression/synthesis/test_needle_pin.py` holds a Linux digest
from before the W1 fold; Windows was re-recorded at W1. W4, W5, W7, the
integral fold and W6 did not move the Windows digest. The Linux value
must be re-read from CI (`ubuntu-24.04`, release) and recorded with
provenance, or the Linux job fails on the first CI run of this branch.

### 1.3 Two notions of "total thickness" (debt, W6)

`DesignStack::total_thickness_nm` sums every film;
`coating_thickness_nm` (new) leaves out incoherent films. The thickness
budget, the pipeline's `total_thickness_nm` and `final_total_thickness_nm`
and the inflate report use the coating sum (equal bit for bit on coherent
stacks). The Python `DesignStack.total_thickness()` (and the stack's
`repr`) still call `total_thickness_nm`, so they count a plate. Decide whether it should be renamed, or the
report keys say "coating".

### 1.4 The PD reference thickness counts an incoherent plate (latent)

`simulate_inner` takes `D = Σ thicknesses[1..nl-1]` for differential
phase. On a stack with a plate this would include it. It cannot be
reached today: every phase demand is refused on an incoherent stack.
It matters if phase demands inside the first block are ever allowed
(§3, W6 follow-up 3).

### 1.5 Plan text out of step (debt)

- `synthesis_targets_plan.md` §1.6 still says construction refuses
  `coherent = false`; W6 lifted that for fixed films. Updated in this
  change, together with the W6 status note.
- §1.3b's row "any channel through an incoherent layer: blocked by W6"
  is now intensity-only done, phase refused.
- §6 follow-ups (target-kinds doc, needle docs with a generated
  supported-target matrix) have not been done for any workstream.

---

## 2. Deferred and not-done items, by workstream

### W3, capability check
- **Item 2 not done, deliberately.** `merit()` keeps the 1e6
  missing-curve penalty for runtime failures (open decision D2).

### W5, transforms in the needle
- `complex` transform: refused at the first fold (no fold form exists).
- Ratio, extinction and contrast targets: not started; would reuse the
  tangent-pair construction.

### W6, incoherent synthesis (first cut, this change)
Done: fixed incoherent films anywhere in the film list; forward through
the engine's block cascade (front rows bit-identical to `solve_point`);
needle P through the cascade adjoint for every intensity and colour
bucket (FD oracle, 0° and 30°); LM by finite differences; budget,
clamp, merge and inflate leave the plate alone. Not done:

1. **Analytic Jacobian through the cascade.** LM differences on any
   stack with a plate (n_par extra solves per iteration). The adjoint
   weights exist (`CascadeTrack.g`); deposits per block times them would
   close it.
2. **Multi-environment runs** with an incoherent film are refused
   (`simulate_all`). The joint door's gap places, backside segments and
   per-environment expansion have not been checked against a spacer.
3. **Phase demands** are refused on any stack with a spacer, including
   phases that live entirely in the first block (front `r`) or the last
   one. The engine reports first-block amplitudes; whether synthesis
   should accept "phase of the front coating alone" is decision D4.
4. **An uncoated face cannot grow a coating.** The plate is never a host,
   and a boundary site with no coherent film above it gets P = 0, so a
   face with no film on it stays bare. A face that should be coated
   needs at least one seed film. Documented here; not yet in user docs.
5. **No explicit exit medium.** The plan proposed
   `DesignStack.exit: Option<LayerSpec>` with the substrate as a finite
   incoherent layer. Implemented instead as "the plate is a film, the
   exit medium is the substrate". Same physics, no new field; the
   Python doors already pass `coherent` per film. Record it as the
   adopted model or schedule the field (decision D5).
6. **Edge materials** (`outer`/`inner`) next to a plate: the edge rule
   names the first/last film, which may now be the plate. Not tested.
7. **Roughness at a plate face**: carried to the solver unchanged, as the
   engine does; no synthesis test covers a rough plate interface.
8. **Back-incidence from an absorbing plate**: `check_back_exit` looks at
   the exit medium only. A strongly absorbing plate makes back demands
   ill-conditioned (tiny T_b) but not undefined; no guard.

### W7, material identity
- `pool_share_tol` (near-equal indices sharing a sweep): not done.
- Shared-sweep and merge **reports** (which entries shared a sweep,
  which films fused): not done.
- The joint (F2) scan sweeps each entry separately (its boundary loci
  differ per entry); only the flat scan shares sweeps.

### W8, new target types
Done: W8a retardance (RetR/RetT, DeltaR/DeltaT). Open, in plan order:
1. Rs/Rp ratio, extinction ratio, diattenuation (tangent pair or grad
   bucket).
2. Psi_R/Psi_T: needs a complex gradient operator
   (`p_coherent_grad_c_from_fields`).
3. DOP / Stokes for a given input state: needs an input-polarization
   field on the target.
4. GD…FOD: port the engine's dispersion forward into `simulate_inner`,
   add curve ids, fold onto the existing Q-row sensitivities; LM by FD
   first.
5. PDrs/PDrp: deferred (reference plane of the round trip, D3).

### Further gaps (plan §5)
- Angle-cone and bandwidth-averaged targets (W4 point weights are the
  building block).
- Colour on back incidence (compile refusal); colour with a spectral
  emphasis profile.
- Weight profiles per environment.
- A randomized standing property test P = dM/dδ (the oracle covers fixed
  cases only).
- General-fit interplay: group `n_factor` as a fit parameter must keep
  W7 labels stable.

---

## 3. Decisions required

| # | Question | Options | Default if nobody decides |
|---|---|---|---|
| D1 | First fold of a run: static (today, `from_spec` folds with no simulation) or an initial simulation in `from_spec`? Only matters with `refold_per_cycle = false`; nonlinear forms (unpolarized pairs, Log, phase) fold nothing or a normal-incidence approximation on the static fold. | static / initial sim | static (refold is the default) |
| D2 | `merit()`'s 1e6 missing-curve penalty after W3: keep for runtime failures, or a hard error everywhere? | keep / error | keep |
| D3 | PDrs/PDrp round-trip reference plane. | wait for a user case / pick the coating front | wait |
| D4 | W6: allow phase demands that live in one coherent block (front reflection of the first block, back reflection of the last), or keep refusing every phase on a stack with a spacer? | allow per block / refuse | refuse |
| D5 | W6 model: keep "plate as an incoherent film, exit medium = substrate", or add `DesignStack.exit` as the plan proposed? | keep / add field | keep |
| D6 | W6 follow-ups: analytic Jacobian through the cascade, and multi-environment support. Which first, and is either needed before a release? | Jacobian / environments / neither | neither |
| D7 | `total_thickness_nm` vs `coating_thickness_nm` (1.3): rename or keep both? | rename the Python method / keep both | keep both |
| D8 | Fix 1.1 (`p_function_multiblock` with a host mask) now, or refuse `host_mask` until it is fixed? | fix / refuse | fix (small) |
| D9 | W7 `pool_share_tol` and the sweep/merge reports: still wanted? | do / drop | drop until asked |

Plan open question 4 (coated back side) is answered by W6: both faces
can be coated, provided each starts with at least one film.

---

## 4. Housekeeping before merge

- Re-read the Linux needle-pin digest from CI (1.2).
- CHANGELOG entries for W7, the integral fold and W6 are in place;
  release ritual (7 version sites, fingerprints) not yet run.
- User docs: the supported-target matrix, the "seed a film on each face"
  rule for plates (W6 item 4), and the refused combinations (phase
  through a plate, multi-environment with a plate).
