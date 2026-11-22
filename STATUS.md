# Research Status

Last updated: 2026-11-22

## Current phase

Phase IV: direct six-field composite source and numerical validation of the connected four-point Schwinger-Dyson equation in the finite-volume 1+1D normal-ordered time-dependent lambda(t) phi^4 model.

## Completed

- Phase I: lowest-order Keldysh cut derivation and finite-volume perturbative benchmark.
- Phase II: exact two- and connected four-point Schwinger-Dyson hierarchy; finite-volume all-orders-in-projected-Hamiltonian dynamics; selected equal-time and unequal-time correlators.
- Phase III: all 16 CTP branch assignments for four selected spacetime configurations; Fock-cutoff convergence and CTP identities.
- Phase IV: direct evaluation of the projected normal-ordered six-field source M33, the three M31*F2 subtraction terms, and the connected four-point Schwinger-Dyson residual for one noncoincident all-forward configuration.

## Phase III quantitative validation

- Four configurations x 16 branch assignments = 64 connected complex components.
- Maximum relative Nmax=12 to Nmax=16 change: 4.60e-12.
- Maximum largest-time identity residual: 4.62e-17 absolute.
- Complement-conjugation residual: zero at stored floating-point precision.
- Independent all-forward comparison with the previous time-ordered implementation: relative difference below 4.31e-13.
- Six Phase-III automated tests pass.

## Phase IV quantitative validation

- Model parameters: m=1, L=20, Kcut=1, lambda0=0.1, T=1; event times are noncoincident and event 1 is latest.
- The source S_4,c is the fourth cumulant with A_1=:phi_1^3:_m treated as one composite argument and phi_2, phi_3, phi_4 as the other arguments. It is not the elementary connected six-field cumulant.
- The projected normal-ordered cubic operator matches both the Hamiltonian commutator and Fourier-projected phi^3 - 3 Delta_m(0) phi.
- At Nmax=8, h=0.00625, the relative Schwinger-Dyson residual is 2.43384e-5.
- Halving the time stencil twice gives 6.12712e-6; the residual follows the expected second-order stencil scaling.
- Richardson fit of the complex residual versus h^2 gives a zero-step intercept magnitude 8.86e-13.
- Nmax=8 to 10 relative changes: F4c 2.96e-9, M33 6.66e-10, connected source 6.58e-10.
- Spatial finite-difference derivatives agree with the spectral derivative within 2.06e-5 to 3.29e-4 over a=0.05 to 0.2.
- Six Phase-IV automated tests pass; free-vacuum absolute SD residual is 2.30e-15.

## Remaining work

- Repeat the six-field source and SD-residual test for additional spacetime configurations and contour assignments.
- Increase momentum-mode cutoff and assess Kcut convergence; current Phase-III/IV calculations use Kcut=1.
- Extend from selected configurations to a structured four-coordinate spacetime grid.
- Derive and compare controlled six-point closures (2PI/Bethe-Salpeter versus 4PI) against the finite-Hilbert reference.
- Complete continuum and infinite-volume extrapolations before making continuum-QFT claims.

## Scientific status

The Phase-IV result is a successful numerical check of one noncoincident all-forward connected four-point equation in a finite-volume, finite-mode Hamiltonian model. It is not a continuum or infinite-volume result. The exact hierarchy has not yet been closed or solved over the full spacetime domain. No final novelty claim is made.
