# Research Status

Last updated: 2026-11-22

## Current phase

Phase VI-A: extended momentum-mode cutoff and Fock-space convergence study for the finite-volume 1+1D normal-ordered time-dependent `lambda(t) phi^4` model.

## Completed

- Phase I: lowest-order Keldysh cut derivation and finite-volume perturbative benchmark.
- Phase II: exact two- and connected four-point Schwinger-Dyson hierarchy; finite-volume all-orders-in-projected-Hamiltonian dynamics; selected equal-time and unequal-time correlators.
- Phase III: all 16 CTP branch assignments for four selected spacetime configurations; Fock-cutoff convergence and CTP identities.
- Phase IV/V: projected normal-ordered six-field composite source and all-forward connected four-point SD residual; four noncoincident spacetime configurations; spatial translation check.
- Phase VI: all 16 CTP source/residual components at K=1 for two configurations; preliminary K=1,2,3 all-forward scan.
- Phase VI-A: matched Nmax=5 all-forward K=1,2,3,4 scan; K=3 Nmax=6 refinement; K=4 Nmax=4,5 comparison; optimized prefix-state memoization implementation cross-checked against the original engine.

## Phase VI-A results

- Matched K=1..4 scan is non-monotonic; K=4 values are below K=3 for both |F4c| and |S4c|.
- K=3 Nmax=5 to 6 changes |F4c| by 0.1269%, and |S4c| by 0.00256%.
- K=3 Nmax=6 SD residual improves from 2.65531e-5 at h=0.00625 to 7.15760e-6 at h=0.003125.
- K=4 Nmax=4 to 5 changes |F4c| by 0.4333%, and |S4c| by 0.1379%.
- K=4 Nmax=4 SD residual is about 2.7e-3; Nmax=5 reduces it to about 1.3e-5.
- K=4 Nmax=5 residual does not decrease when h is halved (1.30105e-5 to 1.40361e-5), so its SD stencil convergence is not established.
- K=4 Nmax=6 is not computed.
- Six Phase VI-A automated tests pass.

## Remaining work

- Improve the K=4 solver/moment evaluation to obtain Nmax=6 and a stable time-stencil residual.
- Extend the K scan to K=5 only after K=4 Fock and stencil convergence are controlled.
- Establish a controlled momentum-cutoff extrapolation before continuum claims.
- Compare controlled 2PI/Bethe-Salpeter and 4PI closures against the finite-Hilbert reference.

## Scientific status

The finite-mode results show significant and non-monotonic K dependence. No momentum-cutoff convergence, continuum limit, or infinite-volume extrapolation has been established. K=4 Nmax=5 SD residual is not time-stencil converged. No final novelty claim is made.