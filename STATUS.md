# Research Status

Last updated: 2026-11-22

## Current phase

Phase III: all 16 CTP contour components of the connected four-point function for selected spacetime configurations in the finite-volume 1+1D normal-ordered time-dependent `lambda(t) phi^4` model.

## Completed

- Phase I: lowest-order Keldysh cut derivation and finite-volume perturbative benchmark.
- Phase II: exact two- and connected four-point Schwinger-Dyson hierarchy; finite-volume all-orders-in-projected-Hamiltonian dynamics; selected equal-time and unequal-time all-forward correlators.
- Phase III: direct branch-state overlap construction for all 16 CTP branch assignments; four spacetime configurations evaluated; three-point spatial-separation scan at fixed unequal times; Nmax=12 to 16 convergence; CTP largest-time and conjugation checks.

## Phase III quantitative validation

- Four configurations × 16 branch assignments = 64 connected complex components.
- Maximum relative Nmax=12 to Nmax=16 change: `4.60e-12`.
- Maximum largest-time identity residual: `4.62e-17` absolute.
- Complement-conjugation residual: zero at stored floating-point precision.
- Independent all-forward comparison with the previous time-ordered implementation: relative difference below `4.31e-13`.
- Six automated tests pass.
- Nmax=16 basis dimension: 969 for Kcut=1.

## Remaining work

- Extend from selected configurations to a structured four-coordinate spacetime grid and additional time orderings.
- Increase momentum-mode cutoff and assess Kcut convergence; current CTP all-branch results use Kcut=1.
- Numerically evaluate the exact four-point Schwinger-Dyson residual, including the six-field composite source.
- Derive and compare controlled six-point closures (2PI/Bethe-Salpeter versus 4PI) against the finite-Hilbert reference.
- Complete continuum and infinite-volume extrapolations before making continuum-QFT claims.

## Scientific status

All reported Phase-III numbers are finite-volume, finite-mode, finite-Fock-space results. The 16 contour components are complete for the four stored spacetime configurations, not for the entire spacetime domain. No final novelty claim is made.
