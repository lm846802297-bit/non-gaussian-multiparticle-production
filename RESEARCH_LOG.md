
## 2026-11-22 — All 16 CTP contour components

- Implemented a branch-state overlap construction: each subset of field insertions is evolved on the forward branch, and the CTP correlator is the overlap of backward- and forward-branch subset states.
- Computed all 16 connected CTP components for four spacetime configurations (64 complex components total) at Kcut=1, Nmax=12 and 16.
- Added a spatial-separation scan r=1,2,3 for a fixed unequal-time pattern and a second later-time configuration.
- Verified the largest-time identity to a maximum absolute residual of 4.62e-17 and complement conjugation at stored floating-point precision.
- The maximum relative Nmax=12 to 16 change across all 64 components is 4.60e-12.
- Cross-checked the all-forward components against the earlier independent time-ordered implementation; maximum relative difference is below 4.31e-13.
- Six automated tests pass. Current CTP all-branch results remain at Kcut=1 and cover selected configurations, not a full four-coordinate grid or continuum limit.
