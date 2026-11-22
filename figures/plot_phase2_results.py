#!/usr/bin/env python3
"""Generate figures from the immutable phase2-dynamics-001 JSON result."""
import json
from pathlib import Path
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = json.loads((ROOT / "results/phase2-dynamics-001/result.json").read_text())
OUT = ROOT / "figures"

coupling = DATA["coupling_scan"]
lam = [r["lambda0_over_m2"] for r in coupling]
ratio = [r["channels"]["generic"]["N2_over_LO"] for r in coupling]
R2 = [r["channels"]["generic"]["R2"] for r in coupling]

fig, ax = plt.subplots()
ax.plot(lam, ratio, marker="o")
ax.set_xlabel(r"$\lambda_0/m^2$")
ax.set_ylabel(r"$N_2/N_2^{(2)}$")
ax.set_title("Exact finite-Fock result relative to leading perturbation")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "phase2_n2_over_leading.png", dpi=180)
fig.savefig(OUT / "phase2_n2_over_leading.svg")
plt.close(fig)

fig, ax = plt.subplots()
ax.semilogy(lam, R2, marker="o")
ax.set_xlabel(r"$\lambda_0/m^2$")
ax.set_ylabel(r"$\mathcal{R}_2$")
ax.set_title("Non-Gaussian correction relative to Gaussian reconstruction")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "phase2_R2_vs_coupling.png", dpi=180)
fig.savefig(OUT / "phase2_R2_vs_coupling.svg")
plt.close(fig)

history = DATA["time_resolved_correlations"]
times = [r["time_over_T"] for r in history]
kappa = [r["mode_cumulants"]["generic"]["kappa4"] for r in history]
fig, ax = plt.subplots()
ax.plot(times, kappa, marker="o", markersize=3)
ax.set_xlabel(r"$t/T$")
ax.set_ylabel(r"$\kappa_4(k=1,q=2;t)$")
ax.set_title("Time-resolved reference-mode connected number cumulant")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "phase2_kappa4_time.png", dpi=180)
fig.savefig(OUT / "phase2_kappa4_time.svg")
plt.close(fig)

# Spatial equal-time connected field correlator at selected times.
fig, ax = plt.subplots()
selected = {0.0: None, 4.0: None, 8.0: None, 12.0: None}
for row in history:
    t = row["time_over_T"]
    if t in selected and "equal_time_G4c_x1112" in row:
        selected[t] = row["equal_time_G4c_x1112"]
for t, records in selected.items():
    if records is None:
        continue
    xs = [r["positions"][3] for r in records]
    ys = [r["G4_connected"][0] for r in records]
    ax.plot(xs, ys, marker="o", label=fr"$t/T={t:g}$")
ax.set_xlabel(r"$r$ (length units)")
ax.set_ylabel(r"$F_{4,c}(t;0,0,0,r)$")
ax.set_title("Equal-time spatial connected four-field correlator")
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "phase2_spatial_G4c.png", dpi=180)
fig.savefig(OUT / "phase2_spatial_G4c.svg")
plt.close(fig)