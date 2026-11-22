
## 2026-11-22 — Phase V: multi-configuration connected four-point SD validation

- Evaluated the projected normal-ordered six-field composite source and connected four-point SD residual at three noncoincident all-forward spacetime configurations.
- Changed both spatial positions and the ordering of the three non-differentiated insertions while keeping the differentiated event latest.
- For each configuration, scanned centered time steps h=0.0125, 0.00625, 0.003125. Relative residuals decrease by factors close to four under h -> h/2, consistent with O(h^2) time-stencil error.
- Spatial finite differences at a=0.1 agree with the retained-mode spectral derivative at relative discrepancy about 8.22e-5.
- Repeated configuration B at Nmax=10; F4c, M33, and composite source change by less than 8e-9 relative to Nmax=8.
- Re-ran all three configurations through the committed single-configuration driver; all reported F4c, M33, source, and residual values match the aggregate results exactly.
- Combined Phase-IV and Phase-V test suite passes 12/12.
- Results remain Kcut=1 finite-mode, finite-volume and all-forward only. Other contour assignments, higher K, and continuum extrapolation remain open.

## 2026-11-22 — Phase VI: all-CTP six-field source and momentum-cutoff scan

- Generalized the branch-state construction to evaluate M33 with the composite insertion A1=:phi_1^3:_m assigned to the same contour branch as event 1, and evaluated all three M31*F2 connected subtractions branch by branch.
- Computed all 16 CTP connected four-point SD residuals for two noncoincident event configurations (A and C), with event 1 latest in physical time.
- At K=1, Nmax=8, h=0.003125, the maximum relative residual is 6.12712e-6 for A and 6.29661e-6 for C. At Nmax=10, configuration A has maximum residual 6.08935e-6.
- Verified largest-time F4c identity (maximum absolute residual 3.78e-17), complement-branch conjugation (zero at stored precision), and agreement of the ++++ component with the independent all-forward engine.
- Across all 16 branches, Nmax=8 to 10 changes are below 1.56e-8 for F4c and below 1.67e-9 for M33 and the composite source.
- Free-vacuum F4c is below 2e-17; absolute SD residuals are below 2.7e-12. Relative residuals are not used in the zero-coupling limit because both equation terms vanish.
- Performed an all-forward Kcut scan at Nmax=5 for K=1,2,3. The magnitudes of F4c and the composite source change substantially with K; no momentum-cutoff convergence is claimed. Fock-cutoff refinements are recorded separately for K=2 and K=3.
- Combined Phase-IV/V/VI automated suite passes 20/20.