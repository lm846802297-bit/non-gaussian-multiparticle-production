# Phase-II figures

All figures are generated from `results/phase2-dynamics-001/result.json` by `plot_phase2_results.py`.

- `phase2_n2_over_leading.svg/png`: exact finite-Fock `N2` divided by the leading `O(lambda0^2)` result versus peak coupling.
- `phase2_R2_vs_coupling.svg/png`: non-Gaussian connected correction relative to the Gaussian reconstruction.
- `phase2_kappa4_time.svg/png`: time-resolved reference-mode connected number cumulant during and after the pulse.
- `phase2_spatial_G4c.svg/png`: equal-time spatial connected four-field correlator at selected times.

SVG is the vector master; PNG is a raster preview. Re-running the plotting script regenerates both formats.