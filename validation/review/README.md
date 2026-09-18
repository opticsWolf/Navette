# Independent verification harnesses

Reproduction scripts for the reviews and the remediation plan —
`docs/code_review.md` (sections 18–23), `docs/remediation_plan.md`
(R4.2, R4.4, R4.6) and `docs/physics_review_coherence.md` (C7). Every
check here is deliberately **independent of the parity reference** — the
oracles are finite differences through the real solver, hand-rolled
Python reimplementations, analytic closed forms, and first-principles
quadrature. The `loom` reference is never imported.

| Script | Review section | What it verifies |
|---|---|---|
| `fd_step1.py` | §19 | merit-fold evaluator (kinds/transforms/weights) vs independent reimplementation; end-to-end `dF/dθ` (fold → needle → solver FD); dispersion ladder (`dφ` vs physical FD, `dgd…dfod` vs independent `np.gradient`); P(z) constancy for host-material needles |
| `fd_rchannel.py` | §19.2 | the pure-reflectance needle channel with a *distinct* needle material — never FD-covered by the parity suite; pass = clean O(h²) convergence (rate ≈ 4) to the analytic slope |
| `color_merit_check.py` | §20 | color-demand merit vs my own XYZ→Lab→ΔE2000 chain (incl. the native-grid white rule); end-to-end color gradient through the solver via a hand-assembled chain rule; CIE whiteness / E313 yellowness / DomWl+purity / XyY / Oklab / Luv |
| `garbage_in.py` | §21 | NaN/inf/duplicate/descending inputs, sinθ>1 aliasing, negative & NaN thickness clamping, needle z range, UniInterpolator knot validation |
| `weaver_race.py` | §22 | OpticalWeaver LRU plan eviction (bit-exact reassembly under thrash) and an 8-thread mixed-op race (no exceptions, no torn writes) |
| `tauc_check.py` | §23 | Tauc–Lorentz ε₂ closed form + ε₁ via independent pair-sampled PV Kramers–Kronig |
| `kk_validate.py` | §23 | the PV-KK quadrature itself, validated against an analytic Lorentz oscillator ε₁ |
| `kk_conv.py` | §23 | h-refinement study attributing the ~1 % near-resonance residual to the FFT-KK grid, not the quadrature |
| `incoherent_check.py` | coherence review C7 | the incoherent cascade against three oracles that are not Navette: the lossless-slab closed form `R = 2R₁/(1+R₁)`, bit-exact thickness independence over three decades, and the phase average itself (the flagged answer IS the coherent answer averaged over one round-trip period) — on R/T and on the Stokes vector, where mode A's refused cross channel is measured missing it |
| `lm_check.py` | R4.4 / §18.2, R4.6 | the scipy parity the plan docs had been *claiming* for the bounded LM in `synthesis/thick_opt.rs` without ever reproducing it: the thin-film case against `least_squares(method="trf")`, the eight hand-written cargo-test optima recomputed by scipy so a wrong pin cannot hide behind a solver that agrees with it, termination reasons, and — part D, the tightest comparison in the plan — `optimizer="trf"` against scipy TRF for the same *optimum*, not merely the same merit (1e-9 on cost, 1e-4 nm on thicknesses) |
| `color_grad_python.py` | R4.2 | the color gradient on the *Python* needle path (`build_needle_targets` → `needle_gradient`), which returned zero for every color demand with no error because the Python dict dropped the fold's `grad_r`/`grad_t`: the fold against FD of the merit, the deposit against a needle call it can be re-derived from (1e-12, same kernel), the whole chain against a thickness FD, and that a color bucket aimed at an uncomputed channel is an error rather than a silent drop |

Run from the repo root with the project venv:

```
.venv/Scripts/python.exe validation/review/<script>.py   # Windows
.venv/bin/python validation/review/<script>.py           # POSIX
```

These are review harnesses, not the regression suite — they are not
collected by `pytest validation` (no `test_*` names). The known
solver behaviors they depend on (negative-thickness clamping, the
native-grid color white) are documented findings, not assumptions to
preserve blindly; re-read §19–§23 before trusting them after API changes.
