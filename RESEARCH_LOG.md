# Research Status

Last updated: 2026-11-22

## Current phase

Phase I: lowest-order Keldysh cut verification and finite-volume perturbative benchmark for the time-dependent, normal-ordered 1+1D `lambda(t) phi^4` model.

## Completed in the research branch

- Derived the leading `O(lambda0^2)` cross-branch contribution to the distinct-mode two-particle inclusive spectrum.
- Derived the normalized finite-volume vacuum-to-four-particle amplitude, including repeated-mode Bose factors.
- Derived the four-cut-Wightman representation and its residual symmetry factor `1/2`.
- Independently evaluated the same leading spectrum by an unordered normalized-Fock-state sum.
- Verified cut-Wightman and Fock-state sums for generic, opposite-momentum, and coincident-mode channels; the coincident channel uses `n_k (n_k - 1)`.
- Selected `lambda(t) = lambda0 sech^2(t/T)` and derived its analytic Fourier transform.
- Checked the Fourier transform against numerical Simpson integration at four frequencies.
- Evaluated finite-volume momentum sums and demonstrated cutoff convergence for the chosen parameters.
- Added three automated tests; all passed.

## Quantitative benchmark

For `m=1`, `L=20`, `T=1`, `lambda0/m^2=0.1`, and measured mode indices `(1,2)`, `(1,-1)`, `(1,1)`, the cutoff-24 leading spectra are recorded in `results/keldysh-cut-001/result.json` and explained in `theory/05-keldysh-cut-propagators-and-pulse-benchmark.md`.

## Not completed

- Higher-order corrections and the complete next-to-leading-order cumulant.
- Full nonperturbative finite-volume Hamiltonian evolution.
- Infinite-volume and continuum extrapolations.
- Independent comparison with an exact Hamiltonian solver.
- Complete citation-chain novelty audit.

## Scope of validation

The completed checks validate the leading perturbative cut/Fock equivalence, the Bose combinatorics, the analytic pulse transform, and momentum-sum convergence for this benchmark. They do not validate nonperturbative dynamics or establish a final novelty claim.
