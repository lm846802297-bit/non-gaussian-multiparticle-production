
## 2026-11-22 — Phase VI-A: extended K and Fock-cutoff scan

- Reran K=1,2,3 at common Nmax=5, h=0.00625, a=0.1, DOP853 rtol=2e-7, atol=2e-9, max_step=0.04; added K=4 at the same settings.
- Matched scan is non-monotonic: K=4 |F4c| and |S4c| are below K=3; K dependence remains large relative to K=1.
- Added K=3, Nmax=6 at h=0.00625 and h=0.003125. Nmax 5 to 6 changes |F4c| by 0.1269% and |S4c| by 0.00256%; time-stencil residual improves from 2.65531e-5 to 7.15760e-6.
- Added K=4, Nmax=4 and 5. Nmax 4 to 5 changes |F4c| by 0.4333% and |S4c| by 0.1379%. Nmax=4 has an inadequate 2.7e-3 SD residual; Nmax=5 gives 1.30105e-5 at h=0.00625.
- K=4, Nmax=5 residual does not improve under h halving (1.30105e-5 to 1.40361e-5); record as unresolved. K=4 Nmax=6 was not computed.
- Added prefix-state memoization in `sd_residual_cached.py`. At K=1, Nmax=5, cached engine reproduces original engine outputs exactly at stored double precision.
- Six Phase VI-A data-integrity and numerical-status tests pass. No continuum extrapolation is claimed.