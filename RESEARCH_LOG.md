
## 2026-11-22 — Phase IV multi-configuration SD validation

- Repeated the direct projected six-field-source calculation at four noncoincident all-forward spacetime configurations, holding m=1, L=20, Kcut=1, lambda0=0.1, T=1, Nmax=8 fixed.
- Relative SD residuals: A reference 6.12712e-6; B common spatial translation 6.13134e-6; C changed times/separations 6.10741e-6; D earlier pulse region 5.37600e-6.
- A common spatial translation by 2.3 leaves F4c, M33, and the connected composite source invariant to relative differences below 2.1e-13.
- Halving the time stencil from 0.003125 to 0.0015625 reduces the reference residual to 1.89285e-6.
- Increasing Nmax from 8 to 10 changes the reference F4c, M33, and connected source by less than 3e-9 relative.
- Five dedicated tests pass. The test suite checks residuals across configurations, translation invariance, time-step refinement, Fock-cutoff refinement, and cancellation of the two SD equation terms.
- The result remains Kcut=1 and all-forward only; no continuum or all-contour source validation is claimed.