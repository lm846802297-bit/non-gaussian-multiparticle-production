# Research Log

## 2026-11-22 — Project bootstrap

- Established the working title and Phase-I model: finite-duration homogeneous interaction pulse in 1+1D `lambda phi^4`.
- Fixed the initial target observables: `N1(k)`, `N2(k,q)`, Gaussian reconstruction `N2_G`, connected cumulant `kappa4`, and relative diagnostic `R2`.
- Recorded prior work showing that general spectrum mappings, Gaussian HBT reconstruction, non-Gaussian initial conditions, four-point kinetic corrections, and MPS real-time `phi^4` studies are not new by themselves.
- Provisional research gap: momentum-resolved connected outgoing two-particle spectrum generated dynamically from a Gaussian in-state under a finite interaction pulse. Priority remains subject to a deeper citation-chain audit.
- No physical numerical computation was run in the bootstrap phase.

## 2026-11-22 — Lowest-order Keldysh/Dyson derivation

- Derived the in-in expansion for `N2(k,q)`; the first nonzero contribution is the cross-branch term with one quartic vertex on each Keldysh branch.
- Derived the normalized finite-volume `0 -> 4` amplitude, including repeated-mode Bose factors.
- Derived the leading finite-volume occupation-state sum for inclusive `N2(k,q)`.
- Clarified that anomalous four-annihilator coherence can be `O(lambda)`, while the number-type connected out-mode cumulant `kappa4(k,q)` starts at `O(lambda^2)` for the normal-ordered free-vacuum pulse.
- Independently checked the Bose combinatoric factor for repeated and distinct occupations.

## 2026-11-22 — Exact CTP hierarchy and finite-volume real-time dynamics

- Derived the exact CTP two-point Schwinger-Dyson equation in terms of `F4,c`.
- Derived the exact connected four-point equation, including the normal-ordering subtraction and the coupling to the six-field composite moment / connected six-point sector.
- Implemented the projected normal-ordered Hamiltonian `H(t)=H0+lambda(t)V` in the total-momentum-zero, even-Z2 Fock sector.
- Evolved the Gaussian free vacuum with `lambda(t)=lambda0 sech(t/T)^2` and computed `N1`, `N2`, `N2_G`, `kappa4`, and `R2` for generic, opposite, and coincident channels.
- Computed equal-time spatial `F4,c(t;0,0,0,r)`, time-resolved reference-mode cumulants, and selected unequal-time all-forward-branch CTP four-point components.
- Coupling scan: for generic `(1,2)`, `N2/N2_LO` falls from 0.98954 at `lambda0/m^2=0.02` to 0.77127 at 0.5.
- Fock-space convergence at `K=2` is demonstrated through `Nmax=12`; momentum cutoff remains not fully converged through `K=7`.
- Nine automated tests pass. Selected unequal-time components converge between `Nmax=12` and 16 below `8e-13` relative. No continuum result or final novelty claim is made.