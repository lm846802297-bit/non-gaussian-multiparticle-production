# Phase VI-A: Momentum-Cutoff and Fock-Space Convergence Extension

Status: K=1..4 all-forward scan completed at matched Nmax=5; K=3 Fock refinement completed through Nmax=6; K=4 Fock refinement completed through Nmax=5 but its SD time-stencil residual is not converged. No continuum result is claimed.

## 1. Objective

Extend the finite-volume connected four-point Schwinger-Dyson benchmark to higher momentum cutoffs while separating momentum-mode cutoff effects from Fock-space truncation effects. The calculation uses the same normal-ordered finite-mode Hamiltonian and projected cubic source as the Phase-IV/V/VI work.

## 2. Matched K scan

The common scan uses

\[
m=1,\quad L=20,\quad \lambda_0=0.1,\quad T=1,
\]

\[
N_{\max}=5,\quad h=0.00625,\quad a=0.1,
\]

with DOP853, `rtol=2e-7`, `atol=2e-9`, and `max_step=0.04`. K=1,2,3 were rerun with these matched solver settings, and K=4 was added.

| K | Fock basis dimension | \(|F_{4,c}|\) | \(|M_{33}|\) | \(|\mathcal S_{4,c}|\) | relative SD residual |
|---:|---:|---:|---:|---:|---:|
| 1 | 56 | 1.51968034e-5 | 1.31475828e-3 | 1.31527954e-3 | 2.04180e-5 |
| 2 | 252 | 3.74294121e-5 | 3.32905683e-3 | 3.33007072e-3 | 1.77581e-5 |
| 3 | 792 | 4.49069636e-5 | 3.80564137e-3 | 3.80412374e-3 | 1.29011e-5 |
| 4 | 2002 | 3.92867110e-5 | 2.88608715e-3 | 2.88171695e-3 | 1.30116e-5 |

The K dependence is non-monotonic. K=4 is 12.5% below K=3 in \(|F_{4,c}|\), and 24.2% below K=3 in \(|\mathcal S_{4,c}|\). Relative to K=1, K=4 is still 158.5% higher in \(|F_{4,c}|\) and 119.1% higher in source magnitude. These values do not demonstrate cutoff convergence.

## 3. Fock-cutoff refinement at K=3

At K=3, Nmax=5 to 6 gives:

- \(|F_{4,c}|\): relative change 1.26874e-3 (0.1269%).
- \(|M_{33}|\): relative change 2.56386e-5.
- \(|\mathcal S_{4,c}|\): relative change 2.56403e-5.

The Nmax=6 SD residual decreases from 2.65531e-5 at h=0.00625 to 7.15760e-6 at h=0.003125, consistent with second-order time-stencil improvement.

## 4. Fock-cutoff refinement at K=4

At K=4, Nmax=4 to 5 gives:

- \(|F_{4,c}|\): relative change 4.33309e-3 (0.4333%).
- \(|M_{33}|\): relative change 1.38033e-3.
- \(|\mathcal S_{4,c}|\): relative change 1.37907e-3.

The Nmax=4 relative SD residual is approximately 2.7e-3, so Nmax=4 is inadequate for the equation-residual test. Nmax=5 reduces it to 1.30105e-5 at h=0.00625. However, halving h to 0.003125 gives 1.40361e-5 rather than a smaller value. This indicates that the K=4, Nmax=5 residual is not yet stencil-converged. K=4, Nmax=6 has not been computed.

## 5. Prefix-state memoization

A new implementation, `src/sd_residual_cached.py`, memoizes intermediate Schrödinger states for shared insertion prefixes. This avoids repeating identical early-time propagations across the 31 moment evaluations used by the SD residual. It does not change the Hamiltonian, regulator, operators, or solver tolerances.

At K=1, Nmax=5, the cached implementation reproduces the original engine's stored `F4c`, `M33`, connected source, kinetic term, and residual exactly at the stored double-precision values. This cross-check is included in the automated tests.

## 6. Validation and provenance

Six Phase VI-A tests pass. They cover cached/original agreement, K=3 Fock stability, K=3 time-stencil refinement, non-monotonic K dependence, the K=4 Nmax=4 residual failure indicator, and norm/Hermiticity checks.

The raw JSON results, source code, test output, environment description, and SHA-256 manifest are archived in `results/phase6-kcut-extension-001/`, `src/`, `tests/`, and `provenance/`.

## 7. Current conclusion and next work

The K=1..4 scan remains non-monotonic and unconverged. K=3 has a controlled Nmax=5 to 6 check for the reported source and correlator. K=4 has a preliminary Nmax=4 to 5 comparison, but its SD residual does not pass the h-refinement check, and Nmax=6 is missing.

The next calculation should improve the K=4 propagation strategy sufficiently to obtain Nmax=6 and repeat the time-stencil scan. Only after this is stable should K=5 and continuum extrapolation be attempted. No continuum or infinite-volume claim is made.