# Research Status

Last updated: 2026-11-22

## Current phase

Phase II: exact CTP hierarchy plus finite-volume real-time Hamiltonian dynamics for the normal-ordered 1+1D time-dependent `lambda(t) phi^4` model.

## Completed in this phase

- Derived the exact CTP two-point Schwinger-Dyson equation in terms of the connected four-point function.
- Derived the exact connected four-point equation and identified its coupling to the six-field composite moment / connected six-point sector.
- Implemented a finite-volume normal-ordered Hamiltonian in the zero-total-momentum, even-Z2 Fock sector.
- Evolved the free Gaussian vacuum through a `sech^2` interaction pulse with a DOP853 real-time solver.
- Computed `N1`, `N2`, all Gaussian two-point contractions, `kappa4`, and `R2` for generic, opposite-momentum, and coincident-mode channels.
- Computed equal-time spatial connected four-field correlator slices `F4c(t;0,0,0,r)`, time-resolved mode-projected cumulants, and selected unequal-time all-forward-branch CTP four-point components.
- Scanned coupling strength, pulse duration, Fock cutoff, momentum-mode cutoff, and solver maximum step.
- Passed eight automated tests, including Hermiticity, vacuum vertex normalization for repeated/distinct occupations, free-vacuum connected four-point vanishing, interaction-generated connected four-point, weak-coupling matching, and coincident-mode factorial moments.

## Main quantitative findings

For the generic `(n_k,n_q)=(1,2)` channel at `m=1`, `L=20`, `T=1`, `K=2`, `Nmax=10`, the exact finite-Hilbert result relative to the leading `O(lambda0^2)` spectrum decreases from 0.98954 at `lambda0/m^2=0.02` to 0.77127 at `lambda0/m^2=0.5`.

The Gaussian reconstruction is far below the full two-particle spectrum in all three measured channels. At `lambda0/m^2=0.2`, its fraction of `N2` is approximately `4.14e-8` (generic), `3.69e-3` (opposite momenta), and `5.43e-8` (coincident mode).

Fock-space convergence is strong at `K=2`: for `lambda0/m^2=0.5`, the generic-channel `N2` changes by about `1.1e-7` relative between `Nmax=10` and 12.

## Remaining work

- Momentum-mode cutoff convergence is not yet sufficient for a continuum claim; the `K=6` to `K=7` change is about 0.65% and the sequence is not monotonic.
- The numerical CTP four-point result covers equal-time spatial slices and selected unequal-time all-forward-branch configurations; a full spacetime grid and all 16 contour-index components remain to be computed.
- The exact CTP hierarchy couples to six-point information; a controlled 4PI / Bethe-Salpeter closure and its truncation error remain to be developed.
- Higher-order perturbative coefficients have not yet been isolated analytically; the finite-Hilbert real-time evolution includes them implicitly within its truncation.
- No infinite-volume or continuum extrapolation has been completed.

## Scientific status

The results establish a finite-volume all-orders-in-the-projected-Hamiltonian benchmark and demonstrate large non-Gaussian contributions for the stated pulse and parameters. They do not establish a continuum QFT prediction. No final novelty claim is made.