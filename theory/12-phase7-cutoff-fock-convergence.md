# Phase VII — Matched Fock/momentum cutoff study

## Result

Phase VII extends the all-forward cutoff scan through K=5 and extends the all-16 CTP six-field-source Schwinger–Dyson validation at K=3 through Nmax=7. The calculation remains finite-volume and finite-mode.

### Common-Nmax=5 momentum scan

Fixed parameters: m=1, L=20, lambda0=0.1, T=1, Nmax=5, h=0.00625. The largest retained momentum is Lambda=2 pi K/L.

| K | Fock dimension | Lambda/m | |F4,c| | |S4,c| |
|---:|---:|---:|---:|---:|
| 1 | 56 | 0.314 | 1.520e-5 | 1.315e-3 |
| 2 | 252 | 0.628 | 3.743e-5 | 3.330e-3 |
| 3 | 792 | 0.942 | 4.491e-5 | 3.804e-3 |
| 4 | 2,002 | 1.257 | 3.929e-5 | 2.882e-3 |
| 5 | 4,368 | 1.571 | 4.156e-5 | 3.677e-3 |

The sequence is nonmonotonic. Adjacent source changes are +153.18%, +14.24%, -24.25%, and +27.59%. At K=5 the physical cutoff is only 1.571 m, so the scan does not yet reach a UV-asymptotic regime.

### All-16 CTP validation at K=3

| Nmax | Fock dimension | Maximum relative SD residual |
|---:|---:|---:|
| 5 | 792 | 1.28665e-3 |
| 6 | 1,716 | 7.04053e-6 |
| 7 | 3,432 | 5.50396e-6 |

At Nmax=7 the largest-time F4,c identity error is 2.025e-11; complementary-branch conjugation is exact at stored precision. From Nmax=6 to 7, the maximum relative changes across all 16 branches are 2.101e-5 for F4,c, 2.668e-6 for M33, and 2.688e-6 for the connected composite source.

The Nmax=5 mixed-branch residual of order 1e-3 falls by over two orders of magnitude at Nmax=6. This identifies the earlier mixed-branch failure as a Fock-boundary truncation effect, not a failure of the CTP branch construction.

### Fock and time-stencil checks

- K=4, Nmax=5 to 6: F4,c changes by 0.1819%; the connected source changes by 0.00362%.
- K=5, Nmax=4 to 5: F4,c changes by 0.5193%; the connected source changes by 0.0918%.
- K=4, Nmax=6: halving h from 0.00625 to 0.003125 lowers the all-forward SD residual from 2.70945e-5 to 9.16434e-6, while F4,c and the source are stable to below 1e-12 relative.
- At K=3, Nmax=7, the CTP ++++ component agrees with the independent all-forward engine to 8.00e-10 in F4,c, 8.96e-12 in M33, and 1.40e-11 in the source.

## Numerical implementation and independent checks

The optimized sparse Hamiltonian and cubic-source builders preserve the original Phase-VI monomial sums and projections. Their matrix elements agree with the original implementation below 2e-14 in the tested cases. The batched CTP state vectors agree with the original shared-prefix tree at relative norm error below 1.7e-15 for K=1, Nmax=4. The K=1, Nmax=5 all-forward output reproduces the archived Phase-VI result at relative differences 2.21e-10 (F4,c), 6.96e-12 (M33), and 6.99e-12 (source).

The Phase-VII automated suite passes 13/13 tests. Source, input, environment, result hashes, and the local native-output bundle are recorded in the provenance manifest.

## Scientific boundary and next work

The fixed free-mass Wick-ordering prescription is retained at every K; no K-dependent retuning of m or lambda0 is introduced. The increasing normal-order contraction diagnostic is recorded but is not applied as an extra counterterm.

No continuum or infinite-volume limit is established. K=6, Nmax=4 is exploratory and has an SD residual of 1.874e-3. K=5, Nmax=6, K=6 with higher Fock cutoffs, and K=4, Nmax=6 all-16 CTP remain open. The next sequence is to finish matched Fock convergence, extend the cutoff to larger Lambda/m, increase L separately, and only then compare 2PI/Bethe–Salpeter and 4PI closures against the finite-Hilbert reference.
