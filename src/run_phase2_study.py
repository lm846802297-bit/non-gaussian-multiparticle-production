#!/usr/bin/env python3
"""Run coupling, truncation, pulse-duration, and equal-time G4c studies."""
from __future__ import annotations
import hashlib
import json
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from real_time_fock_dynamics import (
    evolve, mode_observables, leading_n2, equal_time_connected_four_field
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase2-dynamics-001" / "result.json"
MASS = 1.0
LENGTH = 20.0
CHANNELS = [("generic", 1, 2), ("opposite", 1, -1), ("coincident", 1, 1)]


def obs_record(result, lam0, duration, kcut, channel):
    name, nk, nq = channel
    o = mode_observables(result["modes"], result["basis"], result["psi"], nk, nq)
    lo = leading_n2(kcut, MASS, LENGTH, lam0, duration, nk, nq)
    o["N2_LO"] = lo
    o["N2_over_LO"] = o["N2"] / lo if lo else None
    o["channel"] = name
    return o


def run_single(kcut, nmax, lam0, duration, channels=CHANNELS, save_times=0,
               rtol=5e-10, atol=5e-12, max_step_scale=1.0):
    res = evolve(kcut, nmax, MASS, LENGTH, lam0, duration,
                 rtol=rtol, atol=atol, save_times=save_times,
                 max_step_scale=max_step_scale)
    return res, {name: obs_record(res, lam0, duration, kcut, ch)
                 for name, *rest in channels
                 for ch in [(name, *rest)]}


def main():
    started = time.time()
    output = {
        "study": "phase2-real-time-fock-dynamics",
        "model": "1+1D real scalar, normal-ordered lambda(t) phi^4, free Gaussian vacuum, sech^2 pulse",
        "units": "hbar=c=1; m=1; L=20",
        "channels": {name: {"mode_indices": [k, q],
                             "momenta_over_m": [2*np.pi*k/LENGTH, 2*np.pi*q/LENGTH]}
                     for name, k, q in CHANNELS},
        "coupling_scan": [],
        "fock_cutoff_convergence": [],
        "momentum_cutoff_convergence": [],
        "pulse_duration_scan": [],
        "time_resolved_correlations": [],
        "solver_step_convergence": [],
    }

    # 1. Coupling scan: all channels from the same exact finite-dimensional run.
    for lam0 in [0.02, 0.05, 0.1, 0.2, 0.3, 0.5]:
        res, obs = run_single(2, 10, lam0, 1.0)
        output["coupling_scan"].append({
            "lambda0_over_m2": lam0,
            "basis_dimension": res["dimension"],
            "state_norm": res["norm"],
            "hamiltonian_hermiticity_error": res["hermiticity_error"],
            "solver_nfev": res["nfev"],
            "channels": obs,
        })

    # 2. Fock-space truncation convergence at the strongest scanned coupling.
    for nmax in [4, 6, 8, 10, 12]:
        res, obs = run_single(2, nmax, 0.5, 1.0)
        output["fock_cutoff_convergence"].append({
            "Nmax": nmax, "basis_dimension": res["dimension"],
            "state_norm": res["norm"], "channels": obs,
        })

    # 3. Momentum-mode cutoff convergence. Nmax=6 allows Kcut up to 6 cheaply.
    for kcut in [2, 3, 4, 5, 6]:
        res, obs = run_single(kcut, 6, 0.2, 1.0, channels=[CHANNELS[0]])
        output["momentum_cutoff_convergence"].append({
            "Kcut": kcut, "basis_dimension": res["dimension"],
            "state_norm": res["norm"], "generic_channel": obs["generic"],
        })

    # 4. Pulse-duration scan at fixed peak coupling.
    for duration in [0.5, 0.75, 1.0, 1.25, 1.5]:
        res, obs = run_single(2, 10, 0.2, duration)
        output["pulse_duration_scan"].append({
            "T_times_m": duration, "basis_dimension": res["dimension"],
            "state_norm": res["norm"], "channels": obs,
        })

    # 5. Time-resolved number cumulants and equal-time spatial connected G4.
    res = evolve(2, 8, MASS, LENGTH, 0.2, 1.0, save_times=25,
                 rtol=5e-10, atol=5e-12)
    selected_indices = [12, 16, 20, 24]  # t/T = 0, 4, 8, 12
    positions = [0.0, 2.0, 5.0, 10.0]
    for it, t in enumerate(res["times"]):
        state = res["states"][:, it]
        mode_data = {
            name: mode_observables(res["modes"], res["basis"], state, k, q)
            for name, k, q in CHANNELS
        }
        row = {"time_over_T": float(t), "mode_cumulants": mode_data}
        if it in selected_indices:
            row["equal_time_G4c_x1112"] = []
            for x in positions:
                row["equal_time_G4c_x1112"].append(
                    equal_time_connected_four_field(
                        2, 8, MASS, LENGTH, state, (0.0, 0.0, 0.0, x)
                    )
                )
        output["time_resolved_correlations"].append(row)

    # 6. Solver time-step convergence; tolerances are held fixed.
    for scale in [1.0, 0.5, 0.25]:
        res, obs = run_single(2, 10, 0.5, 1.0, channels=[CHANNELS[0]],
                              max_step_scale=scale)
        output["solver_step_convergence"].append({
            "max_step_scale": scale, "max_step": res["max_step"],
            "nfev": res["nfev"], "elapsed_seconds": res["elapsed"],
            "state_norm": res["norm"], "generic_channel": obs["generic"],
        })

    output["runtime"] = {
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "cpu_count": os.cpu_count(),
    }
    output["elapsed_total_seconds"] = time.time() - started
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(OUT),
        "elapsed_total_seconds": output["elapsed_total_seconds"],
        "coupling_scan": [{
            "lambda": x["lambda0_over_m2"],
            "N2_generic": x["channels"]["generic"]["N2"],
            "R2_generic": x["channels"]["generic"]["R2"],
            "N2_over_LO": x["channels"]["generic"]["N2_over_LO"],
        } for x in output["coupling_scan"]],
        "nmax_convergence": [{"Nmax":x["Nmax"],"N2":x["channels"]["generic"]["N2"]}
                             for x in output["fock_cutoff_convergence"]],
        "kcut_convergence": [{"Kcut":x["Kcut"],"N2":x["generic_channel"]["N2"]}
                             for x in output["momentum_cutoff_convergence"]],
        "duration_scan": [{"T":x["T_times_m"],"N2":x["channels"]["generic"]["N2"],
                           "N2_over_LO":x["channels"]["generic"]["N2_over_LO"]}
                          for x in output["pulse_duration_scan"]],
    }, indent=2))


if __name__ == "__main__":
    main()