# Multi-Configuration Validation of the Connected Four-Point Schwinger-Dyson Equation

Status: completed finite-volume validation for four noncoincident all-forward configurations. No continuum extrapolation is claimed.

## 1. Objective and model

We repeat the connected four-point equation test

\[
(\Box_1+m^2)F_{4,c}(1,2,3,4)+\frac{\lambda(t_1)}{6}\mathcal S_{4,c}(1;2,3,4)=0,
\]

where

\[
\mathcal S_{4,c}=M_{33}-M_{31}(1;2)F_2(3,4)-M_{31}(1;3)F_2(2,4)-M_{31}(1;4)F_2(2,3),
\]

and the composite argument is \(A_1=:\phi_1^3:_m\). The cubic source is projected with the same \(P_K\) as the finite-mode Hamiltonian.

Parameters: \(m=1\), \(L=20\), \(K=1\), \(\lambda_0=0.1\), \(T=1\), \(N_{\max}=8\), DOP853, \(rtol=2\times10^{-9}\), \(atol=2\times10^{-11}\), maximum step 0.02. The primary time stencil is \(h=0.003125\), spatial displacement used for the cross-check is \(a=0.1\), and the spatial derivative in the primary residual is spectral.

## 2. Configurations

All events are mutually noncoincident. Event 1 is the latest-time differentiated coordinate, so the small \(t_1\pm h\) stencil does not cross another insertion time.

| Label | Events \((t_i,x_i)\) |
|---|---|
| A reference | \((1,0),(0.5,1),(-0.5,2),(-1,3)\) |
| B spatial translation | \((1,2.3),(0.5,3.3),(-0.5,4.3),(-1,5.3)\) |
| C changed times/separations | \((1.2,0),(0.4,0.7),(-0.3,1.7),(-1.1,3.2)\) |
| D earlier pulse region | \((0.4,0),(-0.1,0.6),(-0.8,1.4),(-1.5,2.1)\) |

## 3. Results

| Configuration | Relative SD residual | Absolute residual |
|---|---:|---:|
| A reference | \(6.12712\times10^{-6}\) | \(1.12818\times10^{-10}\) |
| B spatial translation | \(6.13134\times10^{-6}\) | \(1.12895\times10^{-10}\) |
| C changed times/separations | \(6.10741\times10^{-6}\) | \(7.88955\times10^{-11}\) |
| D earlier pulse region | \(5.37600\times10^{-6}\) | \(2.24238\times10^{-10}\) |

All four configurations satisfy the finite-difference residual criterion at the selected stencil. The spread in residuals is small and consistent with the same second-order time-stencil error scale.

## 4. Independent translation check

Configurations A and B differ by a common spatial translation of 2.3. The relative differences are:

- \(F_{4,c}\): \(2.06\times10^{-13}\);
- \(M_{33}\): \(3.40\times10^{-16}\);
- connected composite source: \(3.40\times10^{-16}\);
- source term: \(3.79\times10^{-16}\).

This verifies spatial translation invariance to numerical precision for the correlators and composite source. The residuals differ slightly because the second time derivative is evaluated with finite differences.

## 5. Refinement checks

For the reference configuration:

- Halving \(h\) from 0.003125 to 0.0015625 reduces the relative residual from \(6.12712\times10^{-6}\) to \(1.89285\times10^{-6}\).
- Increasing \(N_{\max}\) from 8 to 10 changes the reference \(F_{4,c}\), \(M_{33}\), and connected source by less than \(3.0\times10^{-9}\) relative.
- Five automated multi-configuration tests pass, including translation invariance, residual cancellation, time-step refinement, and Fock-cutoff refinement.

## 6. Interpretation and limitations

The six-field composite source and the differentiated four-point term cancel to the expected finite-difference accuracy at four distinct noncoincident configurations. The spatial translation test provides an independent symmetry check.

This is a finite-volume, finite-mode result at \(K=1\), not a continuum result. The calculation covers all-forward time ordering only. It does not yet validate all CTP branch assignments for the six-field source, nor establish convergence in the momentum-mode cutoff. Those remain separate tasks.