#!/usr/bin/env python3
"""Run multi-configuration connected four-point Schwinger-Dyson checks."""
from __future__ import annotations
import json
from pathlib import Path
from sd_residual import evaluate

CONFIGS = {
    "A_reference": [(1.0, 0.0), (0.5, 1.0), (-0.5, 2.0), (-1.0, 3.0)],
    "B_spatial_translation": [(1.0, 2.3), (0.5, 3.3), (-0.5, 4.3), (-1.0, 5.3)],
    "C_changed_separations_times": [(1.2, 0.0), (0.4, 0.7), (-0.3, 1.7), (-1.1, 3.2)],
    "D_early_pulse_side": [(0.4, 0.0), (-0.1, 0.6), (-0.8, 1.4), (-1.5, 2.1)],
}


def main():
    outdir = Path("results/phase4-multiconfig-001")
    outdir.mkdir(parents=True, exist_ok=True)
    summary = []
    for name, events in CONFIGS.items():
        result = evaluate(
            kcut=1, nmax=8, mass=1.0, L=20.0, lam=0.1, T=1.0,
            h=0.003125, a=0.1, events=events,
            rtol=2e-9, atol=2e-11, max_step=0.02,
        )
        result["configuration_name"] = name
        (outdir / f"{name}.json").write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n"
        )
        summary.append({
            "configuration_name": name,
            "events": events,
            "F4c": result["F4c"],
            "M33_projected": result["M33_projected"],
            "connected_source": result["connected_source"],
            "box_plus_mass_F4c": result["box_plus_mass_F4c"],
            "source_term": result["source_term"],
            "SD_residual": result["SD_residual"],
            "SD_residual_abs": result["SD_residual_abs"],
            "SD_relative_residual": result["SD_relative_residual"],
        })
    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()