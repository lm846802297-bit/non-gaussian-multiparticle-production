# CTP closure benchmarks

This package contains **two distinct tests**, not one high-momentum 2PI/4PI result:

1. **K=5, configuration A:** independent strict-order Born and one-loop (all three bubble channels) for all 16 connected CTP four-point components, against the pre-existing finite-Fock reference `references/k5-n7.json` (K=5, Nmax=7, t_init=-12, λ0=.1). No nPI resummation is applied at K=5.
2. **K=0:** an actual self-consistent three-loop 2PI calculation with consistent bilocal Bethe–Salpeter derivative, and a separate self-consistent four-loop 4PI calculation, both compared to independent oscillator matrix evolution at matching t_init=-4, m=1, L=20, T=1 and external times (1,.5,-.5,-1). A single momentum mode makes space positions immaterial. λ0=.1,.5,2 map to oscillator couplings g0=λ0/L=.005,.025,.1.

No continuum, infinite-volume, finite-time-tail, or inclusive out-particle closure result is claimed here. In particular the K=0 closure comparison cannot predict the quality of either closure for K=5 or for an out-number cumulant.

## Equations and conventions

`F4c` is the elementary fourth contour cumulant with no extra factor of i. The initial state is the free vacuum; (p_n=2πn/L) and (|n|≤K). Both vertices have signed contour measure (∫_C=∫_+-∫_-). Because the interaction is free-vacuum normal ordered, no self-contraction appears in the strict (O(λ_0)) and (O(λ_0^2)) diagrams.

- Born: (-i ∫_C dz,λ(z) ∏_{j=1}^4 G_0(j,z)).
- One loop: (-\frac12 \sum_{(ij|kl)} ∫_C dz\,dw\, λ(z)λ(w) G_0(i,z)G_0(j,z)G_0(z,w)^2 G_0(k,w)G_0(l,w)), where the sum is over the three pairings, and (z,w) include space integration.
- The attached fixed-order code does spatial sums by conserved pair momentum and integrates time with Gauss–Legendre, splitting diagonal panels into Duffy triangles so the time-order cusp is never crossed within a panel.
- At K=0 use (q=\sqrt{L}\,\phi), (g=λ/L), and (G=\langle T_Cqq\rangle). Three-loop 2PI solves (G=G_0+G_0K[G]G) with (K_{\rm loc}=-(i/2)g(G(z,z)-1/(2m))δ_C) and (K_{\rm nl}=-(1/6)g(z)g(w)G(z,w)^3) in this Dyson-insertion convention. A normalized response to a symmetric bilocal source obeys (X=X_0+G K'[G;X]G), with (X_0=G_{\cdot k}G_{l\cdot}+G_{\cdot l}G_{k\cdot}); (F_{4,c}=(X-X_0)/L^2) in the external field normalization.
- Four-loop 4PI independently solves a symmetric proper vertex (V=V_0+\frac12\sum_{ij|kl} V (G⊗G)V) and its Dyson equation in a weighted symmetric contour pair basis. (V_0) carries the local (-ig) and the branch/time weight. The nonlocal self-energy has one bare vertex and one self-consistent vertex, with three joining propagators and symmetry factor 1/6. (F_{4,c}) is the stationary vertex dressed by the four external propagators, divided by (L^2). The raw solver self-consistently iterates *both* (G) and (V).
- Plain trapezoid values and the extrapolation are retained in raw JSON. The postprocessor subtracts the same-grid Born and one-loop terms from each closure and adds the converged quadrature values for those terms, leaving only the (O(λ_0^3))-and-higher closure part to extrapolate in (h^2). This is a quadrature correction, not a change to the closure equations; checks compare the two finest corrected meshes. Relative error means ( |F_{4,c}^{\rm test}-F_{4,c}^{\rm Fock}|/|F_{4,c}^{\rm Fock}| ) per branch.

Primary derivation references: [Berges, hep-ph/0409233](https://arxiv.org/abs/hep-ph/0409233); [Carrington et al., arXiv:1310.4352, Eqs. 10–13 and 14–19](https://arxiv.org/abs/1310.4352); [Carrington et al., arXiv:1310.3295](https://arxiv.org/abs/1310.3295). Those equilibrium 4PI calculations are not themselves nonequilibrium pulse results; the contour continuation is checked here against independently computed weak-coupling coefficients.

## Reproduce

Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 were used. Work in the directory containing this README. Source digests and the K=5 input SHA-256 are recorded in the corresponding native outputs.

```sh
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
python -m unittest -v test_born_bubble test_2pi_k0 test_4pi_k0 test_discrete_weak_order
python born_bubble.py --orders 40 48 56 --fock references/k5-n7.json k5-bubble-vs-fock.json
python closure_pilot.py k0-closure-pilot.json
python closure_pilot.py --h .125 --couplings 2 -- k0-closure-pilot-l2-finer.json
python summarize_closure.py k0-closure-pilot.json k0-closure-pilot-l2-finer.json k0-closure-summary.json
```

The Fock reference in `closure_pilot.py` is built and evaluated as an independent explicit finite oscillator matrix with DOP853. It archives Nmax=12, 16 and a tighter integrator comparison. The K=5 reference was generated from the finite-Fock Phase VII GitHub baseline `eab1f5b3198d03eafc3f162299c9038956beaa0b`, not by either of the two closure codes.

## Next gates before a higher-K closure claim

- Preserve the common finite momentum regulator, free normal ordering and pulse, including two independently validated insertion configurations; Fock and CTP-pulse time discretization are independent errors.
- Extend 2PI first to K=1, checking Fourier sums and all 16 at the same regulator and pulse. For BS at nonzero total pair momentum the bilocal response carries off-diagonal momentum even though the background G is diagonal. Use fixed total-pair-momentum blocks, with exact flat-coordinate comparison on a small sample.
- Extend 4PI using those conservation blocks and normalized symmetric pair indices, not a dense coordinate-time tensor. For each channel combination, retain the crossing permutations and full overlap constraints. Naive pair dimension scales as (O(N_t^2 N_p)) per fixed momentum block; dense contraction has severe (N_t^6) scaling, so short time windows, separable bare pieces, and validated low-rank remainder compression or similarly controlled matrix-free methods are needed before K=5.
- Carry forward weak-order, contour conjugation, largest-time, matched Fock and time-mesh gates. First establish raw residual/systematic accounts at K≥1 and two λ0 values. Then project to (N_1,N_2,\kappa_4) in generic, opposite, and identical modes, with the identical-mode factorial moment. Compare to direct Hamiltonian out-mode expectation, not a single unequal-time coordinate four-point sample.
- Vary physical momentum cutoff at fixed L, and L at approximately fixed physical cutoff; test changing t_init separately. Claim only the regime supported by these checks, and report any unresolved discrepancy relative to benchmark or original Phase VIII.

## Observed results

- At K=5, λ0=.1, the largest branch error of the independent strict-order value relative to Nmax=7 falls from **9.4393e-3** (Born) to **1.3615e-4** (through one loop). Changing Gauss orders 40→48→56 changes the connected through-one-loop values by at most **2.8e-19 absolute**. This is a pre-closure validation result.
- At K=0 the time-quadrature-corrected h² extrapolation gives the following maximum relative differences across 16 branches. "Finest-pair shift" is the maximum difference between the two finest corrected meshes divided by the finest value, not a rigorous truncation bound.

| lambda0 | 2PI/BS vs Fock | 4PI vs Fock | 2PI finest-pair shift | 4PI finest-pair shift |
| --- | ---: | ---: | ---: | ---: |
| 0.1 | 3.80e-6 | 9.60e-7 | 7.85e-8 | 2.19e-7 |
| 0.5 | 9.40e-5 | 2.37e-5 | 1.96e-6 | 5.40e-6 |
| 2 | 1.45e-3 | 3.66e-4 | 1.10e-5 | 2.91e-5 |

At lambda0=2, the plain h² extrapolation from h=1/6 and h=1/8 gives 1.48e-3 and 3.46e-4, respectively, supporting the approximate factor-four improvement without the quadrature correction. For lambda0=.1 or .5 the plain extrapolation retains a larger common quadrature systematic, so use the explicit correction and keep the shift diagnostic. The raw grids show a reversal at some spacings; do not infer closure accuracy directly from coarse raw quadrature. The maximum oscillator Fock-reference changes at lambda0=2 are 1.55e-9 (Nmax 12→16) and 1.54e-12 (integrator).
