# Research Log

## 2026-11-22 — Project bootstrap

- Established the working title and Phase-I model: finite-duration homogeneous interaction pulse in 1+1D lambda phi^4.
- Fixed the initial target observables: N1(k), N2(k,q), Gaussian reconstruction N2_G, connected cumulant kappa4, and relative diagnostic R2.
- Recorded prior work showing that general spectrum mappings, Gaussian HBT reconstruction, non-Gaussian initial conditions, four-point kinetic corrections, and MPS real-time phi^4 studies are not new by themselves.
- Provisional research gap: momentum-resolved connected outgoing two-particle spectrum generated dynamically from a Gaussian in-state under a finite interaction pulse. Priority remains subject to a deeper citation-chain audit.
- No physical numerical computation has been run and no physical result is claimed.

## 2026-11-22 — Lowest-order Keldysh/Dyson derivation

- Derived the in-in expansion for N2(k,q); the first nonzero contribution is the cross-branch term with one quartic vertex on each Keldysh branch.
- Derived the normalized finite-volume 0 -> 4 amplitude, including the factorials for repeated momentum occupations.
- Derived the leading finite-volume occupation-state sum for inclusive N2(k,q).
- Clarified that the anomalous four-annihilator coherence can be O(lambda), while the number-type connected out-mode cumulant kappa4(k,q) starts at O(lambda^2) for the normal-ordered free-vacuum pulse.
- A standalone combinatorial test checked five occupation patterns (4; 3+1; 2+2; 2+1+1; 1+1+1+1). All passed at absolute tolerance 1e-14.
- The test is a combinatorial validation only; no physical spectrum or time evolution was computed.
