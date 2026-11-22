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

## Phase IV preliminary: six-field composite source and SD residual

- Implemented the regulated normal-ordered composite operator `:Phi^3:_m` and checked its matrix identity against `Phi^3 - 3 Delta_m(0) Phi` on interior Fock columns (<2e-14).
- Verified in the free Gaussian vacuum that the composite cumulant is `K_{3;3}=6 F12 F13 F14`, even though the elementary connected four-point function vanishes.
- Evaluated the all-forward connected four-point SD residual for Kcut=1, Nmax=4,6,8 and finite-difference steps h=0.1 to 0.0125. The relative residual remains approximately 3.1-3.6%.
- A Kcut=2, Nmax=4 exploratory run gives a much larger residual (~18.3%); this low Fock cutoff is not converged.
- Eight combined automated tests pass.

## Phase IV decision

The composite source implementation is validated, but the finite-projection SD residual is not zero. The projected Fock-space operator algebra must be derived explicitly, including truncation boundary terms, before using this residual as a test of the continuum Schwinger-Dyson equation. No continuum closure or exact-equation verification is claimed.
