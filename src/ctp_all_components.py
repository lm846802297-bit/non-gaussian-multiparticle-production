#!/usr/bin/env python3
"""All 16 CTP branch components of a connected four-field correlator.

Finite-volume 1+1D normal-ordered lambda(t) phi^4 theory.  For a fixed set
of four spacetime events, each branch subset is represented by a forward
Schrodinger state with field insertions in chronological order.  The CTP
component for plus-mask S is <chi[complement S] | chi[S]>.

This computes selected event configurations, not a continuum spacetime grid.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
import os
import platform
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix


def full_basis(kcut: int, nmax: int):
    modes = tuple(range(-kcut, kcut + 1))
    basis = []
    def rec(i, remaining, occ):
        if i == len(modes):
            basis.append(tuple(occ))
            return
        for n in range(remaining + 1):
            occ.append(n)
            rec(i + 1, remaining - n, occ)
            occ.pop()
    rec(0, nmax, [])
    basis.sort(key=lambda s: (sum(s), s))
    return modes, basis, {s: i for i, s in enumerate(basis)}


def build_hamiltonian(kcut: int, nmax: int, mass: float, length: float):
    modes, basis, index = full_basis(kcut, nmax)
    nm = len(modes)
    omega = [math.sqrt((2 * math.pi * p / length) ** 2 + mass ** 2) for p in modes]
    u = [1.0 / math.sqrt(2.0 * w * length) for w in omega]
    monomials = []
    for ncre in range(5):
        for ct in itertools.combinations_with_replacement(range(nm), ncre):
            cc = Counter(ct)
            cp = sum(modes[i] * n for i, n in cc.items())
            cf = math.prod(u[i] ** n / math.factorial(n) for i, n in cc.items())
            for at in itertools.combinations_with_replacement(range(nm), 4 - ncre):
                ac = Counter(at)
                ap = sum(modes[i] * n for i, n in ac.items())
                if ap != cp:
                    continue
                af = math.prod(u[i] ** n / math.factorial(n) for i, n in ac.items())
                cv = tuple(cc.get(i, 0) for i in range(nm))
                av = tuple(ac.get(i, 0) for i in range(nm))
                monomials.append((cv, av, length * cf * af))
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
                rows.append(row)
                cols.append(j)
                vals.append(amp)
    V = coo_matrix((vals, (rows, cols)), shape=(len(basis), len(basis)), dtype=float).tocsr()
    V.sum_duplicates()
    anti = V - V.T
    herm = float(np.max(np.abs(anti.data))) if anti.nnz else 0.0
    energies = np.array([sum(n * w for n, w in zip(s, omega)) for s in basis])
    return modes, basis, index, V, energies, omega, herm


def pulse(t: float, lam: float, duration: float) -> float:
    z = t / duration
    if abs(z) > 350:
        return 0.0
    c = math.cosh(z)
    return lam / (c * c)


def propagate(v, t0, t1, energies, V, lam, duration, rtol=2e-9, atol=2e-11):
    if abs(t1 - t0) < 1e-14:
        return v.copy(), 0
    def rhs(t, y):
        return -1j * (energies * y + pulse(t, lam, duration) * (V @ y))
    max_step = min(duration / 24.0, 0.02)
    sol = solve_ivp(rhs, (t0, t1), v, method="DOP853", rtol=rtol,
                    atol=atol, max_step=max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:, -1], sol.nfev


def field_matrix(modes, basis, index, mass, length, x):
    rows, cols, vals = [], [], []
    for j, state in enumerate(basis):
        for i, nmode in enumerate(modes):
            p = 2.0 * math.pi * nmode / length
            w = math.sqrt(p * p + mass * mass)
            u = 1.0 / math.sqrt(2.0 * w * length)
            if state[i]:
                ns = list(state)
                ns[i] -= 1
                row = index.get(tuple(ns))
                if row is not None:
                    rows.append(row); cols.append(j)
                    vals.append(u * np.exp(1j * p * x) * math.sqrt(state[i]))
            ns = list(state)
            ns[i] += 1
            row = index.get(tuple(ns))
            if row is not None:
                rows.append(row); cols.append(j)
                vals.append(u * np.exp(-1j * p * x) * math.sqrt(state[i] + 1))
    return coo_matrix((vals, (rows, cols)), shape=(len(basis), len(basis)),
                      dtype=np.complex128).tocsr()


def branch_state(mask, events, psi0, t0, tf, operators, energies, V, lam, duration,
                 rtol=2e-9, atol=2e-11):
    """Ket generated by forward evolution with the events in mask inserted."""
    selected = [i for i in range(4) if (mask >> i) & 1]
    selected.sort(key=lambda i: events[i][0])
    state = psi0.copy()
    now = t0
    nfev = 0
    for i in selected:
        ti, xi = events[i]
        state, nf = propagate(state, now, ti, energies, V, lam, duration, rtol, atol)
        nfev += nf
        state = operators[xi] @ state
        now = ti
    state, nf = propagate(state, now, tf, energies, V, lam, duration, rtol, atol)
    nfev += nf
    return state, nfev


def ctp_components(events, kcut=1, nmax=12, mass=1.0, length=20.0,
                   lam=0.2, duration=1.0, rtol=2e-9, atol=2e-11):
    if len(events) != 4:
        raise ValueError("exactly four events are required")
    t0, tf = -12.0 * duration, 12.0 * duration
    if any(not (t0 <= t <= tf) for t, _ in events):
        raise ValueError("all event times must lie inside the evolution interval")
    modes, basis, index, V, energies, omega, herm = build_hamiltonian(kcut, nmax, mass, length)
    psi0 = np.zeros(len(basis), dtype=np.complex128)
    psi0[index[tuple(0 for _ in modes)]] = 1.0
    positions = sorted(set(float(x) for _, x in events))
    operators = {x: field_matrix(modes, basis, index, mass, length, x) for x in positions}
    states = {}
    total_nfev = 0
    start = time.perf_counter()
    for mask in range(16):
        states[mask], nf = branch_state(mask, events, psi0, t0, tf, operators,
                                        energies, V, lam, duration, rtol, atol)
        total_nfev += nf
    state_elapsed = time.perf_counter() - start

    def pair_value(i, j, plus_mask):
        # plus_mask is restricted to i,j; the complement pair is on the minus branch.
        pair_mask = (1 << i) | (1 << j)
        plus = plus_mask & pair_mask
        minus = pair_mask ^ plus
        return np.vdot(states[minus], states[plus])

    components = {}
    for mask in range(16):
        plus_state = states[mask]
        minus_state = states[15 ^ mask]
        f4 = np.vdot(minus_state, plus_state)
        pairs = []
        for i, j, k, l in [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)]:
            gij = pair_value(i, j, mask)
            gkl = pair_value(k, l, mask)
            pairs.append(gij * gkl)
        connected = f4 - sum(pairs)
        labels = "".join("+" if (mask >> i) & 1 else "-" for i in range(4))
        components[labels] = {
            "F4": [float(f4.real), float(f4.imag)],
            "pair_products": [[float(z.real), float(z.imag)] for z in pairs],
            "F4_connected": [float(connected.real), float(connected.imag)],
        }

    conjugation_errors = []
    for mask in range(16):
        label = "".join("+" if (mask >> i) & 1 else "-" for i in range(4))
        complement = 15 ^ mask
        label_c = "".join("+" if (complement >> i) & 1 else "-" for i in range(4))
        z = complex(*components[label]["F4_connected"])
        zc = complex(*components[label_c]["F4_connected"])
        conjugation_errors.append(abs(zc - z.conjugate()))

    out = {
        "model": "1+1D normal-ordered real scalar lambda(t) phi^4, finite full Fock basis",
        "parameters": {"Kcut": kcut, "Nmax": nmax, "m": mass, "L": length,
                       "lambda0": lam, "T": duration, "t0": t0, "tf": tf},
        "events": [[float(t), float(x)] for t, x in events],
        "basis_dimension": len(basis),
        "interaction_hermiticity_error": herm,
        "branch_state_propagation_nfev": total_nfev,
        "branch_state_wall_seconds": state_elapsed,
        "components": components,
        "checks": {
            "number_of_components": len(components),
            "max_complement_conjugation_error": float(max(conjugation_errors)),
            "all_plus_equals_conjugate_all_minus": conjugation_errors[0],
        },
        "runtime": {"python": sys.version, "numpy": np.__version__,
                    "scipy": scipy.__version__, "platform": platform.platform(),
                    "machine": platform.machine(), "cpu_count": os.cpu_count()},
    }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=12)
    ap.add_argument("--kcut", type=int, default=1)
    ap.add_argument("--lambda0", type=float, default=0.2)
    ap.add_argument("--duration", type=float, default=1.0)
    ap.add_argument("--output", type=Path, default=Path("results/phase3-ctp-all-001/result.json"))
    args = ap.parse_args()
    events = [(1.0, 0.0), (0.5, 1.0), (-0.5, 2.0), (-1.0, 3.0)]
    out = ctp_components(events, kcut=args.kcut, nmax=args.nmax,
                         lam=args.lambda0, duration=args.duration)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()
