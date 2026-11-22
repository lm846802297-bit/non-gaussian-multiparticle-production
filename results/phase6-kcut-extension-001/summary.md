# Phase VI-A: Extended momentum-cutoff and Fock-cutoff scan

## Scope

This run extends the all-forward projected six-field-source Schwinger-Dyson benchmark from K=1,2,3 to K=4, and adds the previously missing K=3, Nmax=6 result. It uses the finite-volume 1+1D normal-ordered real scalar `lambda(t) phi^4` model with `lambda(t)=lambda0 sech(t/T)^2`.

## Matched K scan

For K=1,2,3,4 at common Nmax=5, the raw scan uses `m=1`, `L=20`, `lambda0=0.1`, `T=1`, `h=0.00625`, `a=0.1`, DOP853, `rtol=2e-7`, `atol=2e-9`, and `max_step=0.04`. The K=1,2,3 rows were rerun with this matched solver configuration; K=4 was computed at the same settings.

| K | Fock basis dimension | |F4c| | |M33| | |S4c| | relative SD residual |
|---:|---:|---:|---:|---:|
| 1 | 56 | 1.51968034e-5 | 1.31475828e-3 | 1.31527954e-3 | 2.04180e-5 |
| 2 | 252 | 3.74294121e-5 | 3.32905683e-3 | 3.33007072e-3 | 1.77581e-5 |
| 3 | 792 | 4.49069636e-5 | 3.80564137e-3 | 3.80412374e-3 | 1.29011e-5 |
| 4 | 2002 | 3.92867110e-5 | 2.88608715e-3 | 2.88171695e-3 | 1.30105e-5 |

The sequence is non-monotonic: K=4 is below K=3 for both `|F4c|` and `|S4c|`. The K=1 to K=4 changes remain large. This is not evidence of UV convergence.

## Fock-cutoff checks

At K=3, increasing Nmax from 5 to 6 changes:

- `|F4c|` by 0.1269%.
- `|M33|` by 0.00256%.
- `|S4c|` by 0.00256%.

At K=4, increasing Nmax from 4 to 5 changes:

- `|F4c|` by 0.4333%.
- `|M33|` by 0.1380%.
- `|S4c|` by 0.1379%.

The K=4, Nmax=4 relative SD residual is 2.70e-3, while K=4, Nmax=5 is 1.30e-5. This demonstrates that Nmax=4 is inadequate for the equation-residual test at K=4.

## Time-stencil check

At K=3, Nmax=6, halving h from 0.00625 to 0.003125 reduces the relative residual from 2.65531e-5 to 7.15760e-6, consistent with the expected second-order finite-difference trend.

At K=4, Nmax=5, the corresponding residual changes from 1.30105e-5 to 1.40361e-5 rather than decreasing. The raw correlators and source are stable to the stored precision, but the residual is not yet stencil-converged. This is retained as an unresolved validation issue; no SD-convergence claim is made for K=4.

## Software optimization check

A prefix-state memoization version of the correlator engine reuses intermediate Schrödinger states for shared insertion histories. At K=1, Nmax=5, its computed `F4c`, `M33`, connected source, kinetic term, and residual match the original implementation exactly at the stored double-precision values. This optimization changes the execution path only, not the model, regulator, or tolerances.

## Validation

Six Phase VI-A tests pass. They check exact cached/original agreement, K=3 Fock-cutoff stability, K=3 time-stencil refinement, non-monotonic K dependence, the K=4 Nmax=4 failure indicator, norm conservation, and Hamiltonian Hermiticity.

## Scientific conclusion

The K scan is still not converged. K=3 has a controlled Nmax=5 to 6 check for the reported correlators and source. K=4 has a useful Nmax=4 to 5 comparison, but its SD residual fails the h-refinement test and Nmax=6 has not been computed. The next required calculation is a more efficient or higher-resource K=4, Nmax=6 run, followed by a stable time-stencil residual check. No continuum or infinite-volume result is claimed.