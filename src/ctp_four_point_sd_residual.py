#!/usr/bin/env python3
"""Finite-volume numerical check of the connected four-point CTP SD equation.

For the all-forward component, evaluate the connected four-field correlator,
the composite connected source < :phi^3(1): phi(2) phi(3) phi(4) >_c,
and a central finite-difference approximation to (Box_1 + m^2) F4c.
This is a regulated finite-Hilbert-space residual check, not a continuum proof.
"""
from __future__ import annotations
import argparse
import json
import math
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

from ctp_all_components import build_hamiltonian, field_matrix, branch_state, pulse

AKEY = "COMPOSITE_A1"


def composite_phi3_matrix(modes, basis, index, mass, length, x):
    """Matrix of :phi(x)^3:_m, normal ordered in the regulated free basis."""
    nm = len(modes)
    omega = [math.sqrt((2.0 * math.pi * p / length) ** 2 + mass ** 2) for p in modes]
    u = [1.0 / math.sqrt(2.0 * w * length) for w in omega]
    monomials = []
    for ncre in range(4):
        for ct in __import__('itertools').combinations_with_replacement(range(nm), ncre):
            from collections import Counter
            cc = Counter(ct)
            c_mom = sum(modes[i] * n for i, n in cc.items())
            c_factor = math.prod(u[i] ** n / math.factorial(n) for i, n in cc.items())
            for at in __import__('itertools').combinations_with_replacement(range(nm), 3 - ncre):
                ac = Counter(at)
                a_mom = sum(modes[i] * n for i, n in ac.items())
                a_factor = math.prod(u[i] ** n / math.factorial(n) for i, n in ac.items())
                phase = np.exp(1j * (a_mom - c_mom) * 2.0 * math.pi * x / length)
                coeff = 6.0 * c_factor * a_factor * phase
                cv = tuple(cc.get(i, 0) for i in range(nm))
                av = tuple(ac.get(i, 0) for i in range(nm))
                monomials.append((cv, av, coeff))
    rows, cols, vals = [], [], []
    for j, state in enumerate(basis):
        for cv, av, coeff in monomials:
            new = list(state)
            amp = coeff
            valid = True
            for i, nann in enumerate(av):
                ni = new[i]
                if ni < nann:
                    valid = False
                    break
                amp *= math.sqrt(math.factorial(ni) / math.factorial(ni - nann))
                new[i] -= nann
            if not valid:
                continue
            for i, ncre in enumerate(cv):
                ni = new[i]
                amp *= math.sqrt(math.factorial(ni + ncre) / math.factorial(ni))
                new[i] += ncre
            row = index.get(tuple(new))
            if row is not None and amp:
                rows.append(row); cols.append(j); vals.append(amp)
    mat = coo_matrix((vals, (rows, cols)), shape=(len(basis), len(basis)), dtype=np.complex128).tocsr()
    mat.sum_duplicates()
    anti = mat - mat.getH()
    herm = float(np.max(np.abs(anti.data))) if anti.nnz else 0.0
    return mat, herm


def subset_state(mask, events, psi0, t0, tf, operators, energies, V, lam, duration,
                 rtol=2e-9, atol=2e-11):
    return branch_state(mask, events, psi0, t0, tf, operators, energies, V, lam,
                        duration, rtol, atol)[0]


def connected_all_plus(events, psi0, t0, tf, operators, energies, V, lam, duration):
    # All + branch: <chi_empty | chi_subset> for the required subsets.
    needed = {0, 15, 3, 12, 5, 10, 9, 6}
    states = {mask: subset_state(mask, events, psi0, t0, tf, operators, energies, V, lam, duration)
              for mask in needed}
    f4 = np.vdot(states[0], states[15])
    pairs = []
    for a, b, c, d in [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)]:
        mask1 = (1 << a) | (1 << b)
        mask2 = (1 << c) | (1 << d)
        pairs.append(np.vdot(states[0], states[mask1]) * np.vdot(states[0], states[mask2]))
    return f4 - sum(pairs)


def composite_connected_all_plus(events, psi0, t0, tf, operators, energies, V, lam, duration):
    # Event 0 is the composite A=:phi^3(x1):; events 1..3 are elementary fields.
    needed = {0, 15, 3, 5, 9, 12, 10, 6}
    states = {mask: subset_state(mask, events, psi0, t0, tf, operators, energies, V, lam, duration)
              for mask in needed}
    moment = np.vdot(states[0], states[15])
    pair_products = []
    for j, klmask in [(1, (1 << 2) | (1 << 3)),
                      (2, (1 << 1) | (1 << 3)),
                      (3, (1 << 1) | (1 << 2))]:
        ajmask = 1 | (1 << j)
        pair_products.append(np.vdot(states[0], states[ajmask]) * np.vdot(states[0], states[klmask]))
    return moment - sum(pair_products)


def evaluate_at(t1, x1, events_template, psi0, t0, tf, modes, basis, index,
                energies, V, mass, length, lam, duration):
    events = [(t1, AKEY), (events_template[1][0], events_template[1][1]),
              (events_template[2][0], events_template[2][1]),
              (events_template[3][0], events_template[3][1])]
    operators = {x: field_matrix(modes, basis, index, mass, length, x)
                 for x in sorted(set(float(ev[1]) for ev in events[1:]))}
    operators[AKEY], comp_herm = composite_phi3_matrix(modes, basis, index, mass, length, x1)
    # Replace the first event key with A for the composite source; for F4 use phi(x1).
    f_events = [(t1, x1), *events_template[1:]]
    f_ops = {x: field_matrix(modes, basis, index, mass, length, x)
             for x in sorted(set(float(ev[1]) for ev in f_events))}
    f4c = connected_all_plus(f_events, psi0, t0, tf, f_ops, energies, V, lam, duration)
    k_events = [(t1, AKEY), *events_template[1:]]
    k3c = composite_connected_all_plus(k_events, psi0, t0, tf, operators, energies, V, lam, duration)
    return f4c, k3c, comp_herm


def calculate(kcut=1, nmax=8, mass=1.0, length=20.0, lam=0.2, duration=1.0,
              h_values=(0.1, 0.05, 0.025)):
    t0, tf = -12.0 * duration, 12.0 * duration
    modes, basis, index, V, energies, omega, h_herm = build_hamiltonian(kcut, nmax, mass, length)
    psi0 = np.zeros(len(basis), dtype=np.complex128)
    psi0[index[tuple(0 for _ in modes)]] = 1.0
    template = [(1.0, 0.0), (0.5, 1.0), (-0.5, 2.0), (-1.0, 3.0)]
    t_center, x_center = template[0]
    f_center, k_center, comp_herm = evaluate_at(t_center, x_center, template, psi0, t0, tf,
                                                 modes, basis, index, energies, V,
                                                 mass, length, lam, duration)
    records = []
    for h in h_values:
        fp, _, hp = evaluate_at(t_center + h, x_center, template, psi0, t0, tf,
                                modes, basis, index, energies, V, mass, length, lam, duration)
        fm, _, hm = evaluate_at(t_center - h, x_center, template, psi0, t0, tf,
                                modes, basis, index, energies, V, mass, length, lam, duration)
        xp, _, hxp = evaluate_at(t_center, x_center + h, template, psi0, t0, tf,
                                 modes, basis, index, energies, V, mass, length, lam, duration)
        xm, _, hxm = evaluate_at(t_center, x_center - h, template, psi0, t0, tf,
                                 modes, basis, index, energies, V, mass, length, lam, duration)
        dtt = (fp - 2.0 * f_center + fm) / (h * h)
        dxx = (xp - 2.0 * f_center + xm) / (h * h)
        box = dtt - dxx
        source = pulse(t_center, lam, duration) * k_center / 6.0
        residual = box + mass * mass * f_center + source
        scale = max(abs(box), abs(mass * mass * f_center), abs(source), 1e-300)
        records.append({
            "h": h,
            "F4c_center": [float(f_center.real), float(f_center.imag)],
            "K3c_composite_source": [float(k_center.real), float(k_center.imag)],
            "dtt_F4c": [float(dtt.real), float(dtt.imag)],
            "dxx_F4c": [float(dxx.real), float(dxx.imag)],
            "box_plus_mass_F4c": [float((box + mass*mass*f_center).real), float((box + mass*mass*f_center).imag)],
            "lambda_over_6_K3c": [float(source.real), float(source.imag)],
            "SD_residual": [float(residual.real), float(residual.imag)],
            "relative_residual_to_largest_term": float(abs(residual) / scale),
            "composite_hermiticity_error_max": max(comp_herm, hp, hm, hxp, hxm),
        })
    return {
        "run_id": "phase4-ctp-sd-residual-001",
        "model": "1+1D normal-ordered real scalar lambda(t) phi^4, finite full Fock basis",
        "parameters": {"Kcut": kcut, "Nmax": nmax, "m": mass, "L": length,
                       "lambda0": lam, "T": duration, "t0": t0, "tf": tf,
                       "events": template},
        "interaction_hermiticity_error": h_herm,
        "results": records,
        "runtime": {"python": sys.version, "numpy": np.__version__,
                    "scipy": __import__('scipy').__version__, "platform": platform.platform(),
                    "machine": platform.machine(), "cpu_count": os.cpu_count()},
        "limitations": [
            "Only the all-forward CTP component is used in this residual check.",
            "Derivatives are central finite differences in t1 and x1.",
            "Finite mode and Fock cutoffs; no continuum extrapolation.",
            "The result is a numerical residual diagnostic, not a proof of the continuum Schwinger-Dyson identity."
        ]
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nmax', type=int, default=8)
    ap.add_argument('--kcut', type=int, default=1)
    ap.add_argument('--lambda0', type=float, default=0.2)
    ap.add_argument('--output', type=Path, default=Path('results/phase4-ctp-sd-residual-001/result.json'))
    args = ap.parse_args()
    out = calculate(kcut=args.kcut, nmax=args.nmax, lam=args.lambda0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    print(json.dumps(out, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
