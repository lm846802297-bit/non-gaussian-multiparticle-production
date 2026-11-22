# Finite-volume out-mode cumulant benchmark

This is a **finite box, hard mode cutoff, finite Fock truncation** control for the non-Gaussian particle-production project. It never identifies finite-time contour-field cumulants with out-mode number cumulants. (m=T=1), in-vacuum at (t=-12), and occupations read at (t=12), where the sech-squared pulse is already suppressed. Evolution is archived DOP853 with the exact sparse normal-ordered Hamiltonian in the zero-total-momentum, even-particle-number sector. Source: `eab1f5b3198d03eafc3f162299c9038956beaa0b`, `src/real_time_fock_dynamics.py` SHA-256 `205c2ebcecb53d76e2ca596e27047a43561992a9fdad46771f1fa091795b12c1`.

## Measurements

For each unordered pair define (Q_{kq}=langle a_k^dagger a_q^dagger a_q a_kangle). The (k=q) case is a factorial moment, not (langle n_k^2angle). The connected part is (Q_{kq}-Q^{m G}_{kq}), where (Q^{m G}_{kq}=langle n_kanglelangle n_qangle+|langle a_k^dagger a_qangle|^2+|langle a_k a_qangle|^2). The anomalous term is only present for (k+q=0) by translation invariance. (R_2) uses the Gaussian part as denominator; near zero density it can be very large even if the absolute correlator is tiny and should not be used as a standalone convergence measure. Same-pair and opposite-pair modes remain distinguished.

There are two comparisons. At fixed box (L=20), change (K=5	o6) at each (N_{max}), retaining all common **integer** modes. To match the physical cutoff (Lambda=2pi K/L=pi/2), compare ((K,L)=(5,20)) against ((6,24)) only at common **physical** momenta (0,pmpi/2). Raw cumulants refer to unit-normalized discrete box modes. Multiply by (L) only as a sensitivity diagnostic for the anticipated leading mode-spacing scaling, not as a proof of a continuum, volume-normalized, or coarse-grained measurement. Tests near the cutoff are regulator-dependent.

Each native JSON includes the oscillatory solver settings, sampled norm deviation, final particle-number probabilities, all pair moments, and two algebraically distinct forms of the leading order (Q_{kq}): the analytic Fourier mode sum and the weighted norm of the Fock Born state from (V|0angle) with the analytic pulse transform. The time evolution is performed on a finite window, while the analytic Born control uses the full-line pulse. The integrated (Q_{kq}), weighting off-diagonal pairs twice, is checked against the independent particle-number factorial moment. The weaker-coupling cases check approach to the Born result, not full interacting equivalence.

## Replay

Use `ARCHIVED_SRC=/absolute/path/to/src ./replay.sh` to reproduce native JSON, log, summary, and checksum; runtimes can vary. Run `python validate_out_modes.py` to regenerate the audit from the packaged native JSON without evolving again. This phase does not compare 2PI/4PI or extract a continuum limit.

## Why the exact-cutoff same-sign pair behaves differently

At the same-sign edge ((k,q)=(K,K)) the Born channel enforces (r+s=-2K), so the hard-cutoff sum contains just (r=s=-K). At fixed physical (Lambda), the one term is identical in dimensionless pulse/energy factors and the analytic leading discrete result scales exactly as (1/L^2), not (1/L). In the finite-coupling N8 data the two-cutoff-matched values differ by only about 0.144% after multiplying by (L^2), versus 16.55% after multiplying by (L). This is a special exact-edge phase-space effect, not evidence that (L^2) scaling applies throughout the spectrum. At fixed (L=20), increasing K from 5 to 6 for the fixed pair ((5,5)) opens ((r,s)=(-6,-4),(-4,-6)) in addition to ((-5,-5)), explaining the order-unity shift near that edge. The soft and opposite pairs admit many partner modes, rather than just this one term.
