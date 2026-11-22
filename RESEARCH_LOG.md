
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

## 2026-11-23 — Phase VII

- Restored exact Phase-VI sources and verified Git blob identities before optimization.
- Optimized sparse Hamiltonian/cubic-source construction; matrix-element comparisons pass below 2e-14.
- Batched CTP states reproduce the shared-prefix branch tree below 1.7e-15 relative norm error in the K=1,Nmax=4 check.
- Reproduced Phase-VI K=1,Nmax=5 all-forward results to 2.21e-10 (F4c), 6.96e-12 (M33), and 6.99e-12 (source).
- Extended common-Nmax=5 scan to K=5; cutoff dependence is nonmonotonic.
- Extended K=3 all-16 CTP through Nmax=7; maximum relative SD residuals are 1.28665e-3, 7.04053e-6, and 5.50396e-6 at Nmax=5,6,7.
- K=3 Nmax=6->7 maximum relative changes are 2.101e-5 (F4c), 2.668e-6 (M33), 2.688e-6 (source).
- K=4 Nmax=5->6 changes are 0.1819% (F4c), 0.00362% (source); K=5 Nmax=4->5 changes are 0.5193% (F4c), 0.0918% (source).
- K=4,Nmax=6 time-stencil refinement lowers SD residual from 2.70945e-5 to 9.16434e-6 with stable undifferentiated observables.
- K=3,Nmax=7 CTP ++++ agrees with independent all-forward path at 8.00e-10 (F4c), 8.96e-12 (M33), 1.40e-11 (source).
- Phase-VII tests pass 13/13. No continuum/infinite-volume claim.
