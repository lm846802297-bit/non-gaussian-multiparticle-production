# Research Status

Last updated: 2026-11-22

## Current phase

Phase IV: multi-configuration validation of the connected four-point Schwinger-Dyson equation in the finite-volume 1+1D normal-ordered time-dependent lambda(t) phi^4 model.

## Completed

- Phase I: lowest-order Keldysh cut derivation and finite-volume perturbative benchmark.
- Phase II: exact two- and connected four-point Schwinger-Dyson hierarchy; finite-volume all-orders-in-projected-Hamiltonian dynamics; selected equal-time and unequal-time correlators.
- Phase III: all 16 CTP branch assignments for four selected spacetime configurations; Fock-cutoff convergence and CTP identities.
- Phase IV: projected normal-ordered six-field composite source and connected four-point SD residual for one reference configuration.
- Phase IV extension: repeated the all-forward SD residual at four noncoincident spacetime configurations, checked spatial translation invariance, and refined time stencil and Fock cutoff.

## Phase IV multi-configuration results

Parameters: m=1, L=20, Kcut=1, lambda0=0.1, T=1, Nmax=8, h=0.003125, a=0.1; DOP853 with rtol=2e-9 and atol=2e-11.

- A reference: relative SD residual 6.12712e-6.
- B spatially translated by 2.3: 6.13134e-6.
- C changed time separations: 6.10741e-6.
- D earlier pulse region: 5.37600e-6.
- Translation invariance: relative differences in F4c, M33, and connected source below 2.1e-13.
- Halving h from 0.003125 to 0.0015625 reduces the reference residual to 1.89285e-6.
- Increasing Nmax from 8 to 10 changes F4c, M33, and source by less than 3e-9 relative.
- Five multi-configuration tests pass.

## Remaining work

- Repeat the six-field source and SD-residual test for additional CTP contour assignments.
- Increase momentum-mode cutoff and assess Kcut convergence; current Phase-IV calculations use Kcut=1.
- Extend from selected configurations to a structured four-coordinate spacetime grid.
- Derive and compare controlled six-point closures (2PI/Bethe-Salpeter versus 4PI) against the finite-Hilbert reference.
- Complete continuum and infinite-volume extrapolations before making continuum-QFT claims.

## Scientific status

The all-forward connected four-point SD equation has passed numerical consistency checks at four noncoincident configurations in a finite-volume, finite-mode Hamiltonian model. This does not yet establish the full CTP hierarchy or continuum limit. No final novelty claim is made.