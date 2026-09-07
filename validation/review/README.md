# Independent verification harnesses

Reproduction scripts for `docs/code_review.md` (sections 19–23). Every
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
