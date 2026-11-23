# Research Log

## 2026-11-22 — Project bootstrap

- Established the working title and Phase-I model: finite-duration homogeneous interaction pulse in 1+1D `lambda phi^4`.
- Fixed the initial target observables: `N1(k)`, `N2(k,q)`, Gaussian reconstruction `N2_G`, connected cumulant `kappa4`, and relative diagnostic `R2`.
- Recorded prior work showing that general spectrum mappings, Gaussian HBT reconstruction, non-Gaussian initial conditions, four-point kinetic corrections, and MPS real-time `phi^4` studies are not new by themselves.
- Provisional research gap: momentum-resolved connected outgoing two-particle spectrum generated dynamically from a Gaussian in-state under a finite interaction pulse. Priority remains subject to a deeper citation-chain audit.
- No numerical computation has been run and no physics result is claimed.

## 2026-11-22 — Four-cut Keldysh verification and pulse benchmark

- Derived the cross-branch cut expression with four Wightman lines: two measured lines at `k,q` and two unobserved lines summed over `r,s`.
- Counted the contractions explicitly: `(4*3)^2 * 2! / (4!)^2 = 1/2`.
- Rewrote the cut expression as a finite-volume momentum sum and independently compared it with an unordered normalized-Fock-state enumeration.
- Selected `lambda(t)=lambda0 sech^2(t/T)` and used `lambda_tilde(Omega)=lambda0*pi*T^2*Omega/sinh(pi*T*Omega/2)`.
- Evaluated generic `(k,q)=(1,2)`, opposite `(1,-1)`, and coincident `(1,1)` channels for `m=1`, `L=20`, `T=1`, `lambda0/m^2=0.1`.
- Cut and Fock sums agree at relative precision better than `5e-16` for the reported cutoffs. Momentum-sum changes from cutoff 20 to 24 are at or below `2.94e-15` for the three channels.
- Analytic pulse Fourier transform agrees with 200,000-interval Simpson integration to relative errors between `2.44e-14` and `5.42e-12` at the tested frequencies.
- Three automated tests passed. The test suite verifies cut/Fock equivalence, the coincident-mode factorial moment, and the pulse Fourier transform.
- During exploratory validation, an initial Fock enumerator used `n_k^2` in the coincident channel. This was corrected to the operator-required `n_k(n_k-1)`; the corrected independent enumeration agrees with the cut expression. The earlier discrepancy was a validation-code error, not a change to the analytic cut formula.
- No nonperturbative real-time simulation or continuum extrapolation has been performed.
