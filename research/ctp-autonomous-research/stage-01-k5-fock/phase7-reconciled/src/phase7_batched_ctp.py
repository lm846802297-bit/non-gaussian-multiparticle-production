#!/usr/bin/env python3
"""Batched exact finite-Fock CTP propagation for Phase VII.

The 16 plus-branch subset kets are propagated as 16 columns of one state
matrix. This is algebraically equivalent to the shared-prefix tree in the
Phase-VI implementation, but reduces repeated ODE integrations.
"""
from __future__ import annotations
import argparse, json, math, os, platform, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import sd_residual
from sd_residual import Correlator, pulse
from phase7_fast_ops import make_hamiltonian_fast, projected_normal_ordered_cubic_fast
sd_residual.make_hamiltonian = make_hamiltonian_fast
sd_residual.projected_normal_ordered_cubic = projected_normal_ordered_cubic_fast
from ctp_sd_all_branches import ctp_f4c, ctp_source, derivative_dxx_f4c


def evolve_batch(Y, t0, t1, calc):
    if abs(t1-t0) < 1e-14:
        return Y.copy()
    D, ncol = Y.shape
    E, V = calc.E, calc.V
    def rhs(t, y):
        M = y.reshape(D, ncol)
        dM = -1j * (E[:, None] * M + pulse(t, calc.lam, calc.T) * (V @ M))
        return dM.reshape(-1)
    sol = solve_ivp(rhs, (t0, t1), Y.reshape(-1), method="DOP853",
                   rtol=calc.rtol, atol=calc.atol, max_step=calc.max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:, -1].reshape(D, ncol)


def build_branch_states_batched(calc, events, kinds, masks=None):
    n = len(events)
    D = len(calc.basis)
    if masks is None:
        masks = list(range(1 << n))
    else:
        masks = list(masks)
    Y = np.repeat(calc.psi0[:, None], len(masks), axis=1)
    order = sorted(range(n), key=lambda i: events[i][0])
    now = calc.t0
    for i in order:
        t, x = events[i]
        Y = evolve_batch(Y, now, t, calc)
        kind = kinds[i]
        if kind == "phi": op = calc.phi(t, x)
        elif kind == "cubic": op = calc.cubic(t, x)
        elif kind == "dxx": op = calc.dxxphi(t, x)
        else: raise ValueError(kind)
        cols = [j for j, mask in enumerate(masks) if (mask >> i) & 1]
        if cols:
            Y[:, cols] = op @ Y[:, cols]
        now = t
    Y = evolve_batch(Y, now, calc.tf, calc)
    return {mask: Y[:, j].copy() for j, mask in enumerate(masks)}


def evaluate_all(kcut=1, nmax=8, mass=1.0, length=20.0, lam=0.1, duration=1.0,
                 h=0.003125, rtol=2e-9, atol=2e-11, max_step=0.02,
                 events=None):
    if events is None:
        events = [(1.0, 0.0), (0.5, 1.0), (-0.5, 2.0), (-1.0, 3.0)]
    if len(events) != 4 or not (events[0][0] > max(t for t, _ in events[1:])):
        raise ValueError("event 0 must be the latest physical-time event")
    if len({float(t) for t, _ in events}) != 4:
        raise ValueError("all event times must be distinct")
    tic = time.perf_counter()
    calc = Correlator(kcut, nmax, mass, length, lam, duration, rtol, atol, max_step)
    e0 = [(float(t), float(x)) for t, x in events]
    ep = [(e0[0][0] + h, e0[0][1]), *e0[1:]]
    em = [(e0[0][0] - h, e0[0][1]), *e0[1:]]
    kinds_phi = ["phi"] * 4
    kinds_cubic = ["cubic", "phi", "phi", "phi"]
    kinds_dxx = ["dxx", "phi", "phi", "phi"]
    states_phi = build_branch_states_batched(calc, e0, kinds_phi)
    states_plus = build_branch_states_batched(calc, ep, kinds_phi)
    states_minus = build_branch_states_batched(calc, em, kinds_phi)
    # Event 0 is the differentiated/composite insertion. For masks not
    # containing event 0, the dxx/cubic ket is exactly the ordinary phi ket.
    # Build only the eight nontrivial event-0-plus columns and reuse states_phi
    # for their complements. This is an exact algebraic reduction.
    event0_masks = [mask for mask in range(16) if mask & 1]
    states_dxx = dict(states_phi)
    states_dxx.update(build_branch_states_batched(calc, e0, kinds_dxx, event0_masks))
    states_cubic = dict(states_phi)
    states_cubic.update(build_branch_states_batched(calc, e0, kinds_cubic, event0_masks))
    results = {}
    for mask in range(16):
        f0, _ = ctp_f4c(states_phi, mask)
        ftp, _ = ctp_f4c(states_plus, mask)
        ftm, _ = ctp_f4c(states_minus, mask)
        dtt = (ftp - 2.0*f0 + ftm)/(h*h)
        dxx = derivative_dxx_f4c(states_dxx, states_phi, mask)
        box_mass = dtt - dxx + mass*mass*f0
        source, m33, terms = ctp_source(states_phi, states_cubic, mask)
        source_term = pulse(e0[0][0], lam, duration)*source/6.0
        residual = box_mass + source_term
        denom = abs(box_mass) + abs(source_term) + 1e-30
        label = "".join("+" if (mask >> i) & 1 else "-" for i in range(4))
        results[label] = {
            "F4c": [float(f0.real), float(f0.imag)],
            "M33": [float(m33.real), float(m33.imag)],
            "M31_F2_terms": [[float(z.real), float(z.imag)] for z in terms],
            "connected_composite_source": [float(source.real), float(source.imag)],
            "box_plus_mass_F4c": [float(box_mass.real), float(box_mass.imag)],
            "source_term": [float(source_term.real), float(source_term.imag)],
            "SD_residual": [float(residual.real), float(residual.imag)],
            "SD_residual_abs": float(abs(residual)),
            "SD_relative_residual": float(abs(residual)/denom),
        }
    largest_errors, conjugation_errors, residual_conjugation = [], [], []
    for mask in range(8):
        label0 = "".join("+" if (mask >> i) & 1 else "-" for i in range(4))
        label1 = "".join("+" if ((mask | 1) >> i) & 1 else "-" for i in range(4))
        z0 = complex(*results[label0]["F4c"]); z1 = complex(*results[label1]["F4c"])
        r0 = complex(*results[label0]["SD_residual"]); r1 = complex(*results[label1]["SD_residual"])
        largest_errors.append(abs(z0-z1))
        residual_conjugation.append(abs(r0-r1))
    for mask in range(16):
        label = "".join("+" if (mask >> i) & 1 else "-" for i in range(4))
        comp = 15 ^ mask
        labelc = "".join("+" if (comp >> i) & 1 else "-" for i in range(4))
        z = complex(*results[label]["F4c"]); zc = complex(*results[labelc]["F4c"])
        rr = complex(*results[label]["SD_residual"]); rrc = complex(*results[labelc]["SD_residual"])
        conjugation_errors.append(abs(zc-z.conjugate()))
        residual_conjugation.append(abs(rrc-rr.conjugate()))
    return {
        "model": "1+1D normal-ordered real scalar lambda(t) phi^4; finite volume and momentum projection",
        "parameters": {"Kcut": kcut, "Nmax": nmax, "m": mass, "L": length,
                       "lambda0": lam, "T": duration, "h": h,
                       "rtol": rtol, "atol": atol, "max_step": max_step,
                       "events": e0, "differentiated_event": 0,
                       "contour_convention": "+ branch chronological ket; - branch anti-chronological bra"},
        "basis_dimension": len(calc.basis),
        "interaction_hermiticity_error": calc.herm,
        "vacuum_norm_out": float(np.vdot(calc.vac_out, calc.vac_out).real),
        "normal_order_delta": calc.delta,
        "components": results,
        "checks": {
            "number_of_components": len(results),
            "max_largest_time_F4c_residual": float(max(largest_errors)),
            "max_F4c_complement_conjugation_error": float(max(conjugation_errors)),
            "max_residual_largest_time_branch_difference": float(max(residual_conjugation[:8])),
            "max_residual_complement_conjugation_error": float(max(residual_conjugation[8:])),
            "max_SD_relative_residual": float(max(v["SD_relative_residual"] for v in results.values())),
            "min_SD_relative_residual": float(min(v["SD_relative_residual"] for v in results.values())),
        },
        "runtime": {"python": sys.version, "numpy": np.__version__,
                    "scipy": __import__("scipy").__version__,
                    "platform": platform.platform(), "machine": platform.machine(),
                    "processor": platform.processor(), "cpu_count": os.cpu_count(),
                    "wall_seconds": time.perf_counter()-tic},
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kcut", type=int, default=1)
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--mass", type=float, default=1.0)
    ap.add_argument("--length", type=float, default=20.0)
    ap.add_argument("--lambda0", type=float, default=0.1)
    ap.add_argument("--duration", type=float, default=1.0)
    ap.add_argument("--h", type=float, default=0.003125)
    ap.add_argument("--rtol", type=float, default=2e-9)
    ap.add_argument("--atol", type=float, default=2e-11)
    ap.add_argument("--max-step", type=float, default=0.02)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = evaluate_all(kcut=args.kcut, nmax=args.nmax, mass=args.mass,
                       length=args.length, lam=args.lambda0, duration=args.duration,
                       h=args.h, rtol=args.rtol, atol=args.atol, max_step=args.max_step)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True)+"\n")
    print(json.dumps(out, indent=2, sort_keys=True))

if __name__ == "__main__": main()
