# Research Status

Last updated: 2026-11-22

## Current phase

Phase VI: all-16-branch connected four-point Schwinger-Dyson validation and preliminary momentum-cutoff scan in the finite-volume 1+1D normal-ordered time-dependent `lambda(t) phi^4` model.

## Completed

- Phase I: lowest-order Keldysh cut derivation and finite-volume perturbative benchmark.
- Phase II: exact two- and connected four-point Schwinger-Dyson hierarchy; finite-volume all-orders-in-projected-Hamiltonian dynamics; selected equal-time and unequal-time correlators.
- Phase III: all 16 CTP branch assignments for four selected spacetime configurations; Fock-cutoff convergence and CTP identities.
- Phase IV: direct projected normal-ordered six-field composite source and connected four-point SD residual for one noncoincident all-forward configuration.
- Phase V: repeated the six-field-source SD check for multiple noncoincident all-forward configurations, including spatial translation and changed insertion-time patterns.
- Phase VI: computed the six-field composite source and connected four-point SD residual for all 16 CTP assignments at two noncoincident spacetime configurations; compared Nmax=8 and 10; performed an all-forward K=1,2,3 momentum-cutoff scan at selected Fock cutoffs.

## Phase VI quantitative validation

- Configuration A: all 16 branch components; maximum relative SD residual `6.12712e-6` at Nmax=8 and `6.08935e-6` at Nmax=10.
- Configuration C: all 16 branch components; maximum relative SD residual `6.29661e-6` at Nmax=8.
- Maximum largest-time identity residual of F4c: `3.78e-17` absolute.
- Complement-branch conjugation residual: zero at stored floating-point precision for F4c and the SD residual.
- Nmax=8 to 10 maximum relative changes across all branches: F4c `1.56e-8`, M33 `1.67e-9`, composite source `1.66e-9`.
- Free-vacuum F4c is below `2e-17` for all branches; maximum absolute SD residual is below `2.7e-12`.
- All-forward Kcut scan at common Nmax=5 shows substantial cutoff dependence: relative to K=1, |F4c| changes by about 147% at K=2 and 196% at K=3. Momentum-mode convergence is not established.
- Combined Phase-IV/V/VI automated suite: 20/20 passed.

## Remaining work

- Extend all-16-branch SD checks to additional event configurations and larger K.
- Establish matched Fock-space convergence at K=3 and above; K=3, Nmax=6 is not yet available.
- Establish a controlled momentum-mode cutoff study and then finite-volume/infinite-volume extrapolations.
- Compare controlled six-point closures (2PI/Bethe-Salpeter versus 4PI) against the finite-Hilbert reference.
- Derive and test the full contour contact terms before evaluating coincident-time configurations.

## Scientific status

The all-16-branch source/residual test is complete at K=1 for two selected noncoincident configurations. The all-forward K=1,2,3 scan is a preliminary regulator study and shows that the observables are not yet cutoff-converged. These are finite-volume, finite-mode results, not continuum or infinite-volume results. No final novelty claim is made.

## Phase VII update (2026-11-23)

- K=1,...,5 common-Nmax=5 scan is nonmonotonic; Lambda/m reaches 1.571 at K=5.
- K=3 all-16 CTP maximum relative SD residuals: 1.28665e-3 (Nmax=5), 7.04053e-6 (Nmax=6), 5.50396e-6 (Nmax=7).
- K=3 Nmax=6->7 maximum relative changes across all branches: 2.101e-5 (F4c), 2.668e-6 (M33), 2.688e-6 (source).
- K=4 Nmax=5->6 changes: 1.819e-3 (F4c), 3.622e-5 (source).
- K=5 Nmax=4->5 changes: 5.193e-3 (F4c), 9.177e-4 (source).
- Phase VII tests: 13/13 PASS.
- No continuum or infinite-volume limit is established. K=6,Nmax=4 is exploratory.
