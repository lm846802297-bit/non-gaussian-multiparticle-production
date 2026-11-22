
## 2026-11-22 — Phase IV: six-field source and Schwinger-Dyson residual

- Implemented the projected normal-ordered cubic composite operator P_K :phi_K^3:_m, matching the finite-mode Hamiltonian rather than using an unprojected pointwise cube.
- Clarified that S_4,c is the fourth cumulant with the composite A_1=:phi_1^3:_m treated as one argument; it is not the elementary connected six-field cumulant. In a Gaussian vacuum M33 can be nonzero (3! F12 F13 F14) while M31 vanishes.
- Computed M33=<T P_K:phi_1^3:_m phi_2 phi_3 phi_4> and the three M31(1;j) F2 terms for a noncoincident all-forward configuration.
- Evaluated the time derivative by centered finite differences and the spatial derivative spectrally; independently compared spatial finite differences at a=0.05, 0.1, 0.2.
- At Nmax=8, h=0.00625, the relative SD residual is 2.43384e-5; at h=0.003125 it is 6.12712e-6, consistent with O(h^2). A fit versus h^2 gives a zero-step residual intercept magnitude 8.86e-13.
- Nmax=8 to 10 changes F4c, M33, and the connected source by 2.96e-9, 6.66e-10, and 6.58e-10 respectively.
- Independently checked the projected cubic operator against the Hamiltonian commutator and against Fourier-projected phi^3 - 3 Delta_m(0) phi.
- Six Phase-IV tests pass. The result covers one noncoincident all-forward configuration at Kcut=1; no continuum extrapolation is claimed.
