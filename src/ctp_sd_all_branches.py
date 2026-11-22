#!/usr/bin/env python3
"""All 16 CTP components of the connected four-point SD equation.

Finite-volume 1+1D normal-ordered lambda(t) phi^4 Hamiltonian.
The differentiated event is event 0 and must be the latest physical-time event.
The six-field source uses the same projected normal-ordered cubic operator as
src/sd_residual.py. No contact terms occur for the noncoincident events used.
"""
from __future__ import annotations
import argparse, json, math, platform, sys, time
from pathlib import Path
import numpy as np
from sd_residual import Correlator, pulse, evolve


def build_branch_states(calc, events, kinds):
    """Return 16 CTP branch kets via a shared-prefix chronological tree.

    For a subset S, the ket contains operators assigned to the + branch in
    chronological order. The complementary ket is used as the bra, giving
    anti-time ordering on the - branch through the overlap <chi[bar S]|chi[S]>.
    """
    n = len(events)
    order = sorted(range(n), key=lambda i: events[i][0])
    rank = {idx: r for r, idx in enumerate(order)}
    after = {0: (calc.psi0.copy(), calc.t0)}
    states = {}
    # Build each subset after its latest-in-time selected event.
    masks = sorted(range(1, 1 << n), key=lambda m: (m.bit_count(), max(rank[i] for i in range(n) if (m >> i) & 1)))
    for mask in masks:
        selected = [i for i in range(n) if (mask >> i) & 1]
        last = max(selected, key=lambda i: rank[i])
        parent = mask ^ (1 << last)
        state, tprev = after[parent]
        t, x = events[last]
        state = evolve(state, tprev, t, calc.E, calc.V, calc.lam, calc.T,
                       calc.rtol, calc.atol, calc.max_step)
        kind = kinds[last]
        if kind == "phi": op = calc.phi(t, x)
        elif kind == "cubic": op = calc.cubic(t, x)
        elif kind == "dxx": op = calc.dxxphi(t, x)
        else: raise ValueError(kind)
        state = op @ state
        after[mask] = (state, t)
    for mask in range(1 << n):
        state, tprev = after[mask]
        states[mask] = evolve(state, tprev, calc.tf, calc.E, calc.V, calc.lam,
                              calc.T, calc.rtol, calc.atol, calc.max_step)
    return states


def overlap(states, mask, subset_mask):
    plus = mask & subset_mask
    minus = subset_mask ^ plus
    return np.vdot(states[minus], states[plus])


def ctp_pair(states, mask, i, j):
    pair = (1 << i) | (1 << j)
    return overlap(states, mask, pair)


def ctp_f4c(states_phi, mask):
    allmask = 0b1111
    f4 = np.vdot(states_phi[allmask ^ mask], states_phi[mask])
    pairs = (
        ((0, 1), (2, 3)),
        ((0, 2), (1, 3)),
        ((0, 3), (1, 2)),
    )
    sub = 0j
    for (i, j), (k, l) in pairs:
        sub += ctp_pair(states_phi, mask, i, j) * ctp_pair(states_phi, mask, k, l)
    return f4 - sub, f4


def ctp_source(states_phi, states_cubic, mask):
    allmask = 0b1111
    m33 = np.vdot(states_cubic[allmask ^ mask], states_cubic[mask])
    terms = [
        overlap(states_cubic, mask, 0b0011) * ctp_pair(states_phi, mask, 2, 3),
        overlap(states_cubic, mask, 0b0101) * ctp_pair(states_phi, mask, 1, 3),
        overlap(states_cubic, mask, 0b1001) * ctp_pair(states_phi, mask, 1, 2),
    ]
    return m33 - sum(terms), m33, terms


def derivative_dxx_f4c(states_dxx, states_phi, mask):
    allmask = 0b1111
    raw = np.vdot(states_dxx[allmask ^ mask], states_dxx[mask])
    pairs = (
        ((0, 1), (2, 3)),
        ((0, 2), (1, 3)),
        ((0, 3), (1, 2)),
    )
    sub = 0j
    for (i, j), (k, l) in pairs:
        if i == 0:
            d2 = ctp_pair(states_dxx, mask, i, j)
        elif j == 0:
            d2 = ctp_pair(states_dxx, mask, j, i)
        else:
            raise AssertionError("differentiated event must be in one pair")
        sub += d2 * ctp_pair(states_phi, mask, k, l)
    return raw - sub


def evaluate_all(kcut=1, nmax=8, mass=1.0, length=20.0, lam=0.1, duration=1.0,
                 h=0.003125, rtol=2e-9, atol=2e-11, max_step=0.02,
                 events=None):
    if events is None:
        events = [(1.0, 0.0), (0.5, 1.0), (-0.5, 2.0), (-1.0, 3.0)]
    if len(events) != 4 or not (events[0][0] > max(t for t, _ in events[1:])):
        raise ValueError("event 0 must be the latest of four noncoincident events")
    if len({float(t) for t, _ in events}) != 4:
        raise ValueError("all event times must be distinct; contact terms are excluded")
    tic = time.perf_counter()
    calc = Correlator(kcut, nmax, mass, length, lam, duration, rtol, atol, max_step)
    e0 = [(float(t), float(x)) for t, x in events]
    ep = [(e0[0][0] + h, e0[0][1]), *e0[1:]]
    em = [(e0[0][0] - h, e0[0][1]), *e0[1:]]
    kinds_phi = ["phi"] * 4
    kinds_cubic = ["cubic", "phi", "phi", "phi"]
    kinds_dxx = ["dxx", "phi", "phi", "phi"]
    states_phi = build_branch_states(calc, e0, kinds_phi)
    states_plus = build_branch_states(calc, ep, kinds_phi)
    states_minus = build_branch_states(calc, em, kinds_phi)
    states_dxx = build_branch_states(calc, e0, kinds_dxx)
    states_cubic = build_branch_states(calc, e0, kinds_cubic)
    results = {}
    for mask in range(16):
        f0, f4 = ctp_f4c(states_phi, mask)
        ftp, _ = ctp_f4c(states_plus, mask)
        ftm, _ = ctp_f4c(states_minus, mask)
        dtt = (ftp - 2.0 * f0 + ftm) / (h*h)
        dxx = derivative_dxx_f4c(states_dxx, states_phi, mask)
        box_mass = dtt - dxx + mass*mass*f0
        source, m33, terms = ctp_source(states_phi, states_cubic, mask)
        lam1 = pulse(e0[0][0], lam, duration)
        source_term = lam1 * source / 6.0
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
    # CTP identities for the connected correlator and the residual.
    largest_errors, conjugation_errors, residual_conjugation = [], [], []
    for mask in range(8):
        label0 = "".join("+" if (mask >> i) & 1 else "-" for i in range(4))
        label1 = "".join("+" if ((mask | 1) >> i) & 1 else "-" for i in range(4))
        z0 = complex(*results[label0]["F4c"]); z1 = complex(*results[label1]["F4c"])
        r0 = complex(*results[label0]["SD_residual"]); r1 = complex(*results[label1]["SD_residual"])
        largest_errors.append(abs(z0-z1));
        # Event 0 is the latest physical-time insertion; the SD equation should
        # also be independent of its contour branch away from contact surfaces.
        # Residual comparison is recorded separately from F4c's largest-time identity.
        residual_conjugation.append(abs(r0-r1))
    for mask in range(16):
        label = "".join("+" if (mask >> i) & 1 else "-" for i in range(4))
        comp = 15 ^ mask
        labelc = "".join("+" if (comp >> i) & 1 else "-" for i in range(4))
        z = complex(*results[label]["F4c"]); zc = complex(*results[labelc]["F4c"])
        rr = complex(*results[label]["SD_residual"]); rrc = complex(*results[labelc]["SD_residual"])
        conjugation_errors.append(abs(zc-z.conjugate()))
        # Residuals transform by complex conjugation under complementing all branches.
        residual_conjugation.append(abs(rrc-rr.conjugate()))
    runtime = {
        "python": sys.version, "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "platform": platform.platform(), "machine": platform.machine(),
        "processor": platform.processor(), "cpu_count": __import__("os").cpu_count(),
        "wall_seconds": time.perf_counter()-tic,
    }
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
        "runtime": runtime,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kcut", type=int, default=1)
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--lambda0", type=float, default=0.1)
    ap.add_argument("--duration", type=float, default=1.0)
    ap.add_argument("--h", type=float, default=0.003125)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = evaluate_all(kcut=args.kcut, nmax=args.nmax, lam=args.lambda0,
                       duration=args.duration, h=args.h)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True)+"\n")
    print(json.dumps(out, indent=2, sort_keys=True))

if __name__ == "__main__": main()