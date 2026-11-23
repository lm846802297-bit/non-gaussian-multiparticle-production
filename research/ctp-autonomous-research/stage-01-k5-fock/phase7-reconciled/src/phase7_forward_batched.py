#!/usr/bin/env python3
"""Batched all-forward correlator evaluator for Phase VII cutoff scans.

Computes only the moments needed for F4c, its differentiated connected
subtractions, M33, and the three M31*F2 source subtractions. Each requested
T-ordered event subset is a column; all columns are propagated together.
"""
from __future__ import annotations
import argparse, json, os, platform, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import sd_residual
from sd_residual import Correlator, pulse
from phase7_fast_ops import make_hamiltonian_fast, projected_normal_ordered_cubic_fast
sd_residual.make_hamiltonian = make_hamiltonian_fast
sd_residual.projected_normal_ordered_cubic = projected_normal_ordered_cubic_fast


def evolve_batch(Y, t0, t1, calc):
    if abs(t1-t0) < 1e-14:
        return Y.copy()
    D, ncol = Y.shape
    E, V = calc.E, calc.V
    def rhs(t, y):
        M = y.reshape(D, ncol)
        return (-1j * (E[:, None] * M + pulse(t, calc.lam, calc.T) * (V @ M))).reshape(-1)
    sol = solve_ivp(rhs, (t0, t1), Y.reshape(-1), method="DOP853",
                   rtol=calc.rtol, atol=calc.atol, max_step=calc.max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:, -1].reshape(D, ncol)


def forward_moments(calc, events, kinds, masks):
    masks = list(masks)
    D = len(calc.basis)
    Y = np.repeat(calc.psi0[:, None], len(masks), axis=1)
    order = sorted(range(len(events)), key=lambda i: events[i][0])
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
    return {mask: np.vdot(calc.vac_out, Y[:, j]) for j, mask in enumerate(masks)}


def evaluate_forward(kcut=1, nmax=8, mass=1.0, length=20.0, lam=0.1,
                     duration=1.0, h=0.00625, a=0.1, rtol=2e-9,
                     atol=2e-11, max_step=0.1, events=None):
    if events is None:
        events = [(1.0, 0.0), (0.5, 1.0), (-0.5, 2.0), (-1.0, 3.0)]
    tic = time.perf_counter()
    calc = Correlator(kcut, nmax, mass, length, lam, duration, rtol, atol, max_step)
    e0 = [tuple(map(float, e)) for e in events]
    ep = [(e0[0][0] + h, e0[0][1]), *e0[1:]]
    em = [(e0[0][0] - h, e0[0][1]), *e0[1:]]
    pair_masks = (0b0011, 0b0101, 0b1001, 0b0110, 0b1010, 0b1100)
    dxx_masks = (0b1111, 0b0011, 0b0101, 0b1001)
    cubic_masks = (0b1111, 0b0011, 0b0101, 0b1001)
    phi_masks = (0b1111, *pair_masks)
    phi0 = forward_moments(calc, e0, ["phi"]*4, phi_masks)
    phip = forward_moments(calc, ep, ["phi"]*4, phi_masks)
    phim = forward_moments(calc, em, ["phi"]*4, phi_masks)
    dxx = forward_moments(calc, e0, ["dxx", "phi", "phi", "phi"], dxx_masks)
    cubic = forward_moments(calc, e0, ["cubic", "phi", "phi", "phi"], cubic_masks)
    pair_complements = ((0b1100, 0b0011), (0b1010, 0b0101), (0b0110, 0b1001))
    def connected_f4(mom):
        return mom[0b1111] - sum(mom[a]*mom[b] for a,b in pair_complements)
    F0 = connected_f4(phi0)
    Fp = connected_f4(phip)
    Fm = connected_f4(phim)
    dtt = (Fp - 2*F0 + Fm)/(h*h)
    dxx_raw = dxx[0b1111]
    dxx_sub = (dxx[0b0011]*phi0[0b1100]
               + dxx[0b0101]*phi0[0b1010]
               + dxx[0b1001]*phi0[0b0110])
    dxx_connected = dxx_raw - dxx_sub
    src_raw = cubic[0b1111]
    src_sub = (cubic[0b0011]*phi0[0b1100]
               + cubic[0b0101]*phi0[0b1010]
               + cubic[0b1001]*phi0[0b0110])
    source = src_raw - src_sub
    lam1 = pulse(e0[0][0], lam, duration)
    box_mass = dtt - dxx_connected + mass*mass*F0
    source_term = lam1*source/6.0
    residual = box_mass + source_term
    denom = abs(box_mass) + abs(source_term) + 1e-30
    return {
        "model": "1+1D normal-ordered real scalar lambda(t) phi^4; finite volume and momentum projection",
        "parameters": {"Kcut": kcut, "Nmax": nmax, "m": mass, "L": length,
                       "lambda0": lam, "T": duration, "h": h, "a": a,
                       "rtol": rtol, "atol": atol, "max_step": max_step, "events": e0},
        "basis_dimension": len(calc.basis),
        "interaction_hermiticity_error": calc.herm,
        "vacuum_norm_out": float(np.vdot(calc.vac_out, calc.vac_out).real),
        "normal_order_delta": calc.delta,
        "F4c": [float(F0.real), float(F0.imag)],
        "dtt_F4c": [float(dtt.real), float(dtt.imag)],
        "dxx_F4c_spectral": [float(dxx_connected.real), float(dxx_connected.imag)],
        "M33_projected": [float(src_raw.real), float(src_raw.imag)],
        "M31_F2_terms": [[float(z.real), float(z.imag)] for z in (
            cubic[0b0011]*phi0[0b1100], cubic[0b0101]*phi0[0b1010], cubic[0b1001]*phi0[0b0110])],
        "connected_source": [float(source.real), float(source.imag)],
        "lambda_at_t1": lam1,
        "source_term": [float(source_term.real), float(source_term.imag)],
        "box_plus_mass_F4c": [float(box_mass.real), float(box_mass.imag)],
        "SD_residual": [float(residual.real), float(residual.imag)],
        "SD_residual_abs": float(abs(residual)),
        "SD_relative_residual": float(abs(residual)/denom),
        "runtime": {"python": sys.version, "numpy": np.__version__,
                    "scipy": __import__("scipy").__version__,
                    "platform": platform.platform(), "machine": platform.machine(),
                    "processor": platform.processor(), "cpu_count": os.cpu_count(),
                    "wall_seconds": time.perf_counter()-tic},
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--kcut',type=int,default=1); ap.add_argument('--nmax',type=int,default=8)
    ap.add_argument('--mass',type=float,default=1.0); ap.add_argument('--length',type=float,default=20.0)
    ap.add_argument('--lambda0',type=float,default=0.1); ap.add_argument('--duration',type=float,default=1.0)
    ap.add_argument('--h',type=float,default=0.00625); ap.add_argument('--max-step',type=float,default=0.1)
    ap.add_argument('--rtol',type=float,default=2e-9); ap.add_argument('--atol',type=float,default=2e-11)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    d=evaluate_forward(kcut=args.kcut,nmax=args.nmax,mass=args.mass,length=args.length,
                       lam=args.lambda0,duration=args.duration,h=args.h,
                       rtol=args.rtol,atol=args.atol,max_step=args.max_step)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps(d,indent=2,sort_keys=True))
if __name__=='__main__': main()
