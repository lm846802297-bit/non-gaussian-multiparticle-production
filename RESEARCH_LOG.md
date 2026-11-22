
## 2026-11-22 — All 16 CTP contour components

- Implemented a branch-state overlap construction: each subset of field insertions is evolved on the forward branch, and the CTP correlator is the overlap of backward- and forward-branch subset states.
- Computed all 16 connected CTP components for four spacetime configurations (64 complex components total) at Kcut=1, Nmax=12 and 16.
- Added a spatial-separation scan r=1,2,3 for a fixed unequal-time pattern and a second later-time configuration.
- Verified the largest-time identity to a maximum absolute residual of 4.62e-17 and complement conjugation at stored floating-point precision.
- The maximum relative Nmax=12 to 16 change across all 64 components is 4.60e-12.
- Cross-checked the all-forward components against the earlier independent time-ordered implementation; maximum relative difference is below 4.31e-13.
- Six automated tests pass. Current CTP all-branch results remain at Kcut=1 and cover selected configurations, not a full four-coordinate grid or continuum limit.

## 2026-11-22 — Composite six-field source and SD residual diagnostic

- Implemented the finite-mode normal-ordered composite operator `:Phi^3:_m` and independently checked `:Phi^3:=Phi^3-3 Delta_m(0) Phi` on interior Fock columns.
- Corrected an important conceptual point: the composite cumulant `K_{3;3}=<:phi_1^3: phi_2 phi_3 phi_4>_c` is nonzero in a free Gaussian vacuum because the three constituents of `:phi^3:` can contract with the three external fields; it equals `6 F12 F13 F14`. The elementary connected six-point function still vanishes in a Gaussian state.
- Evaluated the all-forward connected four-point Schwinger-Dyson residual with finite differences. At Kcut=1 the relative residual is about 3.1-3.6% across the tested Nmax=4,6,8 and h=0.1 to 0.0125 scans.
- A Kcut=2,Nmax=4 exploratory run gives a larger residual (~18.3%) but is not Fock-converged.
- Eight combined unit tests pass. The residual is classified as preliminary and not as a successful continuum SD identity check; the projected Fock boundary terms remain to be derived.
