# Research Status

Last updated: 2026-11-22

## Current phase

Phase I: lowest-order Keldysh/Dyson derivation for the time-dependent 1+1D normal-ordered lambda(t) phi^4 model.

## Completed in this branch

- Re-derived the in-in expansion for the normal-ordered out-mode pair number operator.
- Identified the first nonzero contribution as the cross-branch term with one interaction vertex on each Keldysh branch, at O(lambda_0^2).
- Derived the normalized finite-volume first-order vacuum-to-four-particle amplitude, including repeated-mode Bose factors.
- Derived the finite-volume occupation-state sum for the leading inclusive N2(k,q).
- Established the leading order separation N2_G = O(lambda_0^4) versus N2 = O(lambda_0^2) for distinct modes under the stated normal-ordering, free-vacuum, and homogeneous-pulse assumptions.
- Distinguished the O(lambda_0) anomalous four-annihilator coherence from the O(lambda_0^2) number-type connected cumulant.
- Independently checked the quartic-vertex occupation-number combinatorics for five representative repeated/distinct-mode patterns; all five tests passed.

## Not completed

- Full independent derivation of all contour-index and Wightman-projection factors.
- Explicit evaluation for the selected smooth pulse.
- Higher-order corrections and complete next-to-leading-order cumulant.
- Numerical finite-volume Hamiltonian implementation.
- Cutoff, volume, and time-step convergence.
- Independent comparison with a full quantum solver.
- Complete citation-chain novelty audit.

## Current decision

The leading finite-volume expression is explicit. The Bose factors in the four-particle vertex matrix element passed an independent combinatorial check. Proceed next to an independent Keldysh/Dyson comparison and an explicit finite-pulse momentum-sum evaluation. No full physical spectrum has yet been computed.
