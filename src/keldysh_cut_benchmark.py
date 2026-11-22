#!/usr/bin/env python3
"""Finite-volume leading-order cut-Keldysh benchmark for a sech^2 pulse.

Model: 1+1D normal-ordered lambda(t) phi^4, free Gaussian in-vacuum.
The script compares the four-cut-Wightman mode sum with an independent
normalized-Fock-occupation enumeration for generic, opposite, and coincident
measured momentum modes. It also checks the analytic pulse Fourier transform
against numerical Simpson integration.

This is a perturbative benchmark, not a real-time interacting simulation.
"""
from __future__ import annotations

import json
import math
import platform
import sys
from collections import defaultdict
from itertools import combinations_with_replacement
from pathlib import Path

M = 1.0
L = 20.0
T = 1.0
LAMBDA0 = 0.1  # lambda_0 / m^2 = 0.1


def omega(n: int) -> float:
    momentum = 2.0 * math.pi * n / L
    return math.sqrt(momentum * momentum + M * M)


def lambda_tilde(energy: float) -> float:
    """Fourier transform of lambda(t)=lambda0*sech(t/T)^2."""
    x = math.pi * T * energy / 2.0
    if abs(x) < 1.0e-7:
        ratio = 1.0 - x*x/6.0 + 7.0*x**4/360.0
    else:
        ratio = x / math.sinh(x)
    return 2.0 * LAMBDA0 * T * ratio


def n2_cut_sum(nk: int, nq: int, cutoff: int) -> float:
    """Four-cut-line Keldysh sum; valid for nk != nq and nk == nq."""
    wk, wq = omega(nk), omega(nq)
    total = 0.0
    for nr in range(-cutoff, cutoff + 1):
        ns = -(nk + nq + nr)
        if abs(ns) > cutoff:
            continue
        wr, ws = omega(nr), omega(ns)
        energy = wk + wq + wr + ws
        lt = lambda_tilde(energy)
        denominator = (2*wk) * (2*wq) * (2*wr) * (2*ws)
        total += lt*lt / denominator
    return total / (2.0 * L * L)


def n2_fock_sum(nk: int, nq: int, cutoff: int) -> float:
    """Independent sum over unordered normalized four-particle Fock states."""
    total = 0.0
    modes = range(-cutoff, cutoff + 1)
    for state in combinations_with_replacement(modes, 4):
        if sum(state) != 0:
            continue
        occupations: dict[int, int] = defaultdict(int)
        for mode in state:
            occupations[mode] += 1
        ck, cq = occupations[nk], occupations[nq]
        # a_k^dagger a_q^dagger a_q a_k gives n_k*n_q for k != q,
        # and n_k*(n_k-1) for k == q.
        weight = ck * (ck - 1) if nk == nq else ck * cq
        if weight == 0:
            continue
        bose_factor = math.prod(math.factorial(n) for n in occupations.values())
        energy_denominator = math.prod((2.0 * omega(p))**n for p, n in occupations.items())
        total_energy = sum(n * omega(p) for p, n in occupations.items())
        lt = lambda_tilde(total_energy)
        probability = lt*lt / (L*L*bose_factor*energy_denominator)
        total += weight * probability
    return total


def n1_fock_sum(nk: int, cutoff: int) -> float:
    """Leading single-particle spectrum from the same four-particle states."""
    total = 0.0
    modes = range(-cutoff, cutoff + 1)
    for state in combinations_with_replacement(modes, 4):
        if sum(state) != 0:
            continue
        occupations: dict[int, int] = defaultdict(int)
        for mode in state:
            occupations[mode] += 1
        count = occupations[nk]
        if count == 0:
            continue
        bose_factor = math.prod(math.factorial(n) for n in occupations.values())
        energy_denominator = math.prod((2.0 * omega(p))**n for p, n in occupations.items())
        total_energy = sum(n * omega(p) for p, n in occupations.items())
        lt = lambda_tilde(total_energy)
        probability = lt*lt / (L*L*bose_factor*energy_denominator)
        total += count * probability
    return total


def sech2(x: float) -> float:
    z = x / T
    if abs(z) > 350.0:
        return 0.0
    c = math.cosh(z)
    return 1.0 / (c*c)


def simpson_fourier(energy: float, extent_in_T: float = 20.0, steps: int = 200_000) -> complex:
    if steps % 2:
        steps += 1
    extent = extent_in_T * T
    h = 2.0 * extent / steps
    acc = 0.0j
    for i in range(steps + 1):
        time = -extent + i*h
        value = LAMBDA0 * sech2(time) * complex(math.cos(energy*time), math.sin(energy*time))
        coefficient = 1 if i in (0, steps) else (4 if i % 2 else 2)
        acc += coefficient * value
    return acc * h / 3.0


def main() -> None:
    output: dict[str, object] = {
        "parameters": {
            "m_over_m": M,
            "L_times_m": L*M,
            "T_times_m": T*M,
            "lambda0_over_m2": LAMBDA0/(M*M),
            "delta_k_over_m": 2.0*math.pi/(L*M),
            "pulse": "lambda(t)=lambda0*sech(t/T)^2",
            "cutoffs": [4, 6, 8, 10, 12, 16, 20, 24, 32, 40],
        },
        "cases": {},
        "fourier_checks": [],
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "visible_cpu_count": __import__("os").cpu_count(),
        },
    }

    cases = [("generic", 1, 2), ("opposite", 1, -1), ("coincident", 1, 1)]
    for name, nk, nq in cases:
        records = []
        previous = None
        for cutoff in output["parameters"]["cutoffs"]:
            cut = n2_cut_sum(nk, nq, cutoff)
            fock = n2_fock_sum(nk, nq, cutoff)
            relative_difference = abs(cut-fock)/max(abs(cut), 1e-300)
            step_change = None if previous is None else abs(cut-previous)/abs(cut)
            assert cut >= 0.0 and fock >= 0.0
            assert relative_difference < 5.0e-14, (name, cutoff, relative_difference)
            records.append({
                "cutoff": cutoff,
                "N2_cut": cut,
                "N2_fock": fock,
                "relative_cut_fock_difference": relative_difference,
                "relative_change_from_previous_cutoff": step_change,
            })
            previous = cut
        output["cases"][name] = {
            "mode_indices": [nk, nq],
            "momenta_over_m": [2*math.pi*nk/(L*M), 2*math.pi*nq/(L*M)],
            "N1_k_cutoff40": n1_fock_sum(nk, 40),
            "cutoff_scan": records,
        }

    for energy in [2.0, 4.0, 6.0, 8.0]:
        analytic = lambda_tilde(energy)
        numerical = simpson_fourier(energy)
        relative_error = abs(abs(numerical)-analytic)/analytic
        assert relative_error < 1.0e-10, (energy, relative_error)
        output["fourier_checks"].append({
            "Omega_over_m": energy/M,
            "analytic_lambda_tilde_over_m": analytic/M,
            "numeric_abs_lambda_tilde_over_m": abs(numerical)/M,
            "relative_error": relative_error,
        })

    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
