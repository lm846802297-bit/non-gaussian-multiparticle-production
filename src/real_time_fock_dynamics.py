#!/usr/bin/env python3
"""Exact finite-dimensional real-time dynamics for 1+1D normal-ordered phi^4.

The model is a real scalar field in a periodic box with a finite momentum-mode
cutoff. Evolution is restricted to total momentum zero and even Z2 parity.
The interaction is normal ordered relative to the free mass-m vacuum. This is
an exact Hamiltonian evolution inside the stated finite-dimensional subspace;
it is not a continuum-limit result.
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
from scipy.sparse import coo_matrix, csr_matrix


def omega(n: int, mass: float, length: float) -> float:
    p = 2.0 * math.pi * n / length
    return math.sqrt(p * p + mass * mass)


def zero_momentum_basis(kcut: int, nmax: int):
    modes = tuple(range(-kcut, kcut + 1))
    mcount = len(modes)
    basis = []

    def rec(i, remaining, occ, total_p):
        if i == mcount:
            if total_p == 0 and sum(occ) % 2 == 0:
                basis.append(tuple(occ))
            return
        p = modes[i]
        for ni in range(remaining + 1):
            occ.append(ni)
            rec(i + 1, remaining - ni, occ, total_p + p * ni)
            occ.pop()

    rec(0, nmax, [], 0)
    basis.sort(key=lambda x: (sum(x), x))
    return modes, basis, {state: i for i, state in enumerate(basis)}


def build_interaction(kcut: int, nmax: int, mass: float, length: float):
    """Build the projected V=(1/4!) int dx :phi^4: matrix.

    The monomial coefficient for a^dagger_c a_a is obtained by summing all
    ordered field factors and the 1/4! vertex factor analytically. For creator
    and annihilator occupation multisets c,a, it is
       L prod_p u_p^(c_p+a_p)/(c_p! a_p!), u_p=(2 omega_p L)^(-1/2),
    provided sum_p p c_p = sum_p p a_p.
    """
    modes, basis, index = zero_momentum_basis(kcut, nmax)
    mcount = len(modes)
    ws = [omega(n, mass, length) for n in modes]
    u = [1.0 / math.sqrt(2.0 * w * length) for w in ws]
    monomials = []
    for ncre in range(5):
        for c_tuple in itertools.combinations_with_replacement(range(mcount), ncre):
            c_count = Counter(c_tuple)
            c_mom = sum(modes[i] * n for i, n in c_count.items())
            c_factor = math.prod(u[i] ** n / math.factorial(n) for i, n in c_count.items())
            for ann_tuple in itertools.combinations_with_replacement(range(mcount), 4 - ncre):
                a_count = Counter(ann_tuple)
                a_mom = sum(modes[i] * n for i, n in a_count.items())
                if a_mom != c_mom:
                    continue
                a_factor = math.prod(u[i] ** n / math.factorial(n) for i, n in a_count.items())
                coeff = length * c_factor * a_factor
                c_vec = tuple(c_count.get(i, 0) for i in range(mcount))
                a_vec = tuple(a_count.get(i, 0) for i in range(mcount))
                monomials.append((c_vec, a_vec, coeff))

    rows, cols, vals = [], [], []
    for j, state in enumerate(basis):
        for c_vec, a_vec, coeff in monomials:
            new = list(state)
            amp = coeff
            valid = True
            for i, nann in enumerate(a_vec):
                ni = new[i]
                if ni < nann:
                    valid = False
                    break
                amp *= math.sqrt(math.factorial(ni) / math.factorial(ni - nann))
                new[i] -= nann
            if not valid:
                continue
            for i, ncre in enumerate(c_vec):
                ni = new[i]
                amp *= math.sqrt(math.factorial(ni + ncre) / math.factorial(ni))
                new[i] += ncre
            irow = index.get(tuple(new))
            if irow is not None and amp:
                rows.append(irow)
                cols.append(j)
                vals.append(amp)
    V = coo_matrix((vals, (rows, cols)), shape=(len(basis), len(basis)), dtype=float).tocsr()
    V.sum_duplicates()
    anti = V - V.T
    herm_error = float(np.max(np.abs(anti.data))) if anti.nnz else 0.0
    energies = np.array([sum(n * w for n, w in zip(state, ws)) for state in basis])
    return modes, basis, index, V, energies, ws, herm_error, len(monomials)


def pulse(t: float, lam0: float, duration: float) -> float:
    z = t / duration
    if abs(z) > 350.0:
        return 0.0
    c = math.cosh(z)
    return lam0 / (c * c)


def evolve(kcut: int, nmax: int, mass: float, length: float, lam0: float,
           duration: float, rtol: float = 2e-10, atol: float = 2e-12,
           save_times: int = 0, max_step_scale: float = 1.0):
    modes, basis, index, V, energies, ws, herm_error, nmon = build_interaction(
        kcut, nmax, mass, length
    )
    psi0 = np.zeros(len(basis), dtype=np.complex128)
    psi0[index[tuple(0 for _ in modes)]] = 1.0
    tmax = 12.0 * duration
    max_step = max_step_scale * min(duration / 24.0, 0.025 / max(ws))

    def rhs(t, psi):
        return -1j * (energies * psi + pulse(t, lam0, duration) * (V @ psi))

    t_eval = np.linspace(-tmax, tmax, save_times) if save_times else None
    tic = time.time()
    sol = solve_ivp(rhs, (-tmax, tmax), psi0, method="DOP853", rtol=rtol,
                    atol=atol, max_step=max_step, t_eval=t_eval)
    elapsed = time.time() - tic
    if not sol.success:
        raise RuntimeError(sol.message)
    states = sol.y
    psi = states[:, -1]
    norm = float(np.vdot(psi, psi).real)
    return {
        "modes": modes, "basis": basis, "index": index, "V": V,
        "energies": energies, "ws": ws, "psi": psi,
        "times": sol.t, "states": states, "norm": norm,
        "hermiticity_error": herm_error, "dimension": len(basis),
        "monomials": nmon, "nfev": sol.nfev, "elapsed": elapsed,
        "max_step": max_step, "tmax": tmax,
    }


def mode_observables(modes, basis, psi, nk: int, nq: int):
    i_k, i_q = modes.index(nk), modes.index(nq)
    probs = np.abs(psi) ** 2
    n1k = n1q = n2 = 0.0
    for state, prob in zip(basis, probs):
        ck, cq = state[i_k], state[i_q]
        n1k += prob * ck
        n1q += prob * cq
        n2 += prob * (ck * (ck - 1) if nk == nq else ck * cq)

    # Translation invariance in the total-P=0 sector makes C_kq diagonal.
    Ckq = complex(n1k if nk == nq else 0.0)
    Akq = 0.0j
    if nk + nq == 0:
        index = {s: i for i, s in enumerate(basis)}
        for j, state in enumerate(basis):
            occ = list(state)
            if nk == nq:
                c = occ[i_k]
                if c < 2:
                    continue
                amp = math.sqrt(c * (c - 1))
                occ[i_k] -= 2
            else:
                if occ[i_k] == 0 or occ[i_q] == 0:
                    continue
                amp = math.sqrt(occ[i_k] * occ[i_q])
                occ[i_k] -= 1
                occ[i_q] -= 1
            ibra = index.get(tuple(occ))
            if ibra is not None:
                Akq += np.conjugate(psi[ibra]) * amp * psi[j]

    if nk == nq:
        n2g = 2.0 * n1k * n1k + abs(Akq) ** 2
    else:
        n2g = n1k * n1q + abs(Ckq) ** 2 + abs(Akq) ** 2
    kappa = n2 - n2g
    return {
        "N1_k": float(n1k), "N1_q": float(n1q), "N2": float(n2),
        "N2_G": float(n2g), "kappa4": float(kappa),
        "R2": float(kappa / n2g) if n2g > 0 else None,
        "C_kq": [Ckq.real, Ckq.imag], "A_kq": [Akq.real, Akq.imag],
    }


def lambda_tilde(energy: float, lam0: float, duration: float) -> float:
    x = math.pi * duration * energy / 2.0
    if abs(x) < 1e-7:
        return 2.0 * lam0 * duration * (1.0 - x*x/6.0 + 7*x**4/360.0)
    return 2.0 * lam0 * duration * x / math.sinh(x)


def leading_n2(kcut: int, mass: float, length: float, lam0: float,
               duration: float, nk: int, nq: int) -> float:
    wk, wq = omega(nk, mass, length), omega(nq, mass, length)
    total = 0.0
    for nr in range(-kcut, kcut + 1):
        ns = -(nk + nq + nr)
        if abs(ns) > kcut:
            continue
        wr, ws = omega(nr, mass, length), omega(ns, mass, length)
        E = wk + wq + wr + ws
        den = (2*wk)*(2*wq)*(2*wr)*(2*ws)
        total += lambda_tilde(E, lam0, duration)**2 / den
    return total / (2.0 * length * length)


def full_basis(kcut: int, nmax: int):
    """All Fock states, all momentum sectors, total N <= nmax."""
    modes = tuple(range(-kcut, kcut + 1))
    mcount = len(modes)
    basis = []
    def rec(i, remaining, occ):
        if i == mcount:
            basis.append(tuple(occ))
            return
        for ni in range(remaining + 1):
            occ.append(ni)
            rec(i + 1, remaining - ni, occ)
            occ.pop()
    rec(0, nmax, [])
    basis.sort(key=lambda x: (sum(x), x))
    return modes, basis, {s: i for i, s in enumerate(basis)}


def field_matrix(modes, basis, index, mass, length, x):
    rows, cols, vals = [], [], []
    for j, state in enumerate(basis):
        for i, nmode in enumerate(modes):
            w = omega(nmode, mass, length)
            p = 2.0 * math.pi * nmode / length
            u = 1.0 / math.sqrt(2.0 * w * length)
            if state[i] > 0:
                new = list(state); new[i] -= 1
                row = index.get(tuple(new))
                if row is not None:
                    rows.append(row); cols.append(j)
                    vals.append(u * np.exp(1j * p * x) * math.sqrt(state[i]))
            new = list(state); new[i] += 1
            row = index.get(tuple(new))
            if row is not None:
                rows.append(row); cols.append(j)
                vals.append(u * np.exp(-1j * p * x) * math.sqrt(state[i] + 1))
    return coo_matrix((vals, (rows, cols)), shape=(len(basis), len(basis)), dtype=np.complex128).tocsr()


def equal_time_connected_four_field(kcut, nmax_state, mass, length, psi_state,
                                     positions):
    """Compute equal-time C4,c at four points using a larger observable basis.

    The observable basis is extended by four particles so four successive
    field insertions are not cut off at the state-space boundary.
    """
    modes, basis, index = full_basis(kcut, nmax_state + 4)
    state_modes, state_basis, _ = zero_momentum_basis(kcut, nmax_state)
    assert modes == state_modes
    embed = np.zeros(len(basis), dtype=np.complex128)
    full_index = index
    for occ, amp in zip(state_basis, psi_state):
        embed[full_index[occ]] = amp
    phis = {x: field_matrix(modes, basis, index, mass, length, x) for x in set(positions)}

    def moment(xs):
        v = embed
        for x in reversed(xs):
            v = phis[x] @ v
        return np.vdot(embed, v)

    c4 = moment(positions)
    x1, x2, x3, x4 = positions
    c12, c34 = moment((x1, x2)), moment((x3, x4))
    c13, c24 = moment((x1, x3)), moment((x2, x4))
    c14, c23 = moment((x1, x4)), moment((x2, x3))
    c4c = c4 - c12*c34 - c13*c24 - c14*c23
    return {
        "positions": list(positions),
        "G4": [float(c4.real), float(c4.imag)],
        "G4_connected": [float(c4c.real), float(c4c.imag)],
        "pair_products": [[float((c12*c34).real), float((c12*c34).imag)],
                          [float((c13*c24).real), float((c13*c24).imag)],
                          [float((c14*c23).real), float((c14*c23).imag)]],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kcut", type=int, default=2)
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--mass", type=float, default=1.0)
    ap.add_argument("--length", type=float, default=20.0)
    ap.add_argument("--lambda0", type=float, default=0.1)
    ap.add_argument("--duration", type=float, default=1.0)
    ap.add_argument("--nk", type=int, default=1)
    ap.add_argument("--nq", type=int, default=2)
    ap.add_argument("--save-times", type=int, default=0)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    result = evolve(args.kcut, args.nmax, args.mass, args.length, args.lambda0,
                    args.duration, save_times=args.save_times)
    obs = mode_observables(result["modes"], result["basis"], result["psi"], args.nk, args.nq)
    lo = leading_n2(args.kcut, args.mass, args.length, args.lambda0,
                    args.duration, args.nk, args.nq)
    output = {
        "parameters": vars(args) | {"output": str(args.output) if args.output else None},
        "basis_dimension": result["dimension"], "interaction_monomials": result["monomials"],
        "interaction_hermiticity_error": result["hermiticity_error"],
        "state_norm": result["norm"], "solver_nfev": result["nfev"],
        "elapsed_seconds": result["elapsed"], "max_step": result["max_step"],
        "observables": obs, "N2_leading_perturbative": lo,
        "relative_N2_to_LO": obs["N2"] / lo if lo else None,
        "runtime": {"python": sys.version, "numpy": np.__version__,
                    "scipy": scipy.__version__, "platform": platform.platform(),
                    "machine": platform.machine(), "processor": platform.processor(),
                    "cpu_count": os.cpu_count()},
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()