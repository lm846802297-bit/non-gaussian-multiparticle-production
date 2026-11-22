# Phase III CTP figures

- `phase3-ctp-components-vs-r.svg`: vector master showing the real and imaginary parts of selected connected CTP branch components as the spatial extent parameter varies from 1 to 3.
- `plot_ctp_components_vs_r.py`: deterministic plotting source. It reads the Nmax=16 JSON outputs in `results/phase3-ctp-all-001/`.

The figure is a finite-volume, finite-mode result at Kcut=1, not a continuum prediction. A PNG preview is retained in the local run archive but is not part of the Git commit; the SVG is the archival figure.
