# Lowest-Order Perturbative Benchmark: Work Plan

## 1. First-order vacuum-to-four-particle amplitude

With the interaction-picture convention

\[
H_I(t)=\frac{\lambda(t)}{4!}\int_0^L dx\,:\!\phi_I^4(t,x)\!:\,_m,
\]

the first-order amplitude for four distinct normalized momentum modes is

\[
\mathcal M_{0\to k_1k_2k_3k_4}^{(1)}
=-i\int dt\,\langle k_1k_2k_3k_4|H_I(t)|0\rangle.
\]

With the finite-volume mode normalization in the model document, the matrix element reduces to a factor proportional to

\[
\frac{\delta_{k_1+k_2+k_3+k_4,0}}{L\sqrt{(2\omega_1)(2\omega_2)(2\omega_3)(2\omega_4)}}\,\widetilde\lambda(\omega_1+\omega_2+\omega_3+\omega_4).
\]

The exact combinatorial factor depends on whether the four momenta are distinct and on the normalized Fock-state convention. It must be derived and unit-tested rather than copied into code from a shorthand formula.

This amplitude shows the physical mechanism: a spatially homogeneous time-dependent coupling preserves total momentum, while the external time dependence can supply energy. A finite pulse has nonzero Fourier support at positive total final energy and can therefore create four particles from the free vacuum.

## 2. Expected order counting for the two-particle inclusive observable

For a free vacuum and a normal-ordered quartic pulse, the first-order state correction contains a four-particle component. The operator `a_k^dagger a_q^dagger a_q a_k` has a nonzero diagonal matrix element in that component. Therefore its leading contribution can occur at order `lambda_0^2` for suitable mode pairs. The final two-point occupations are also generated at order `lambda_0^2`, so the Gaussian products built from them are generally of higher order for distinct modes.

This order-counting statement is a prediction to verify with explicit Dyson expansion, including all relevant final states and mode degeneracies. It is not yet a computed result.

## 3. Required derivation

1. Fix normalized finite-volume Fock states, including repeated-momentum cases.
2. Derive the first-order `0 -> 4` amplitude with all factorials and Kronecker deltas.
3. Construct the outgoing state through the order required for `N1` and `N2`.
4. Evaluate `N1(k)`, `N2(k,q)`, and every Wick contraction independently.
5. Extract `kappa4(k,q)` and determine its lowest nonzero order.
6. Compare with the closed-time-path calculation using the same pulse and normalization.

## 4. Validation checks

- `lambda0 -> 0` gives zero out-particle production.
- A spatially homogeneous pulse enforces total momentum conservation.
- The amplitude depends on the pulse Fourier transform at the total final energy.
- The result is invariant under an overall time translation of the pulse up to the expected phase; probabilities are invariant.
- Repeated and distinct momentum modes give consistent limits after correct Bose factors.
- The perturbative result agrees with direct finite-volume Dyson expansion at sufficiently small coupling.