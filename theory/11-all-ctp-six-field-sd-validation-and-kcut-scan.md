# All-CTP Six-Field Source Validation and Momentum-Cutoff Scan

Status: completed finite-volume, finite-mode validation for all 16 CTP components at two noncoincident spacetime configurations. A preliminary K=1,2,3 all-forward cutoff scan is also archived. No continuum extrapolation is claimed.

## 1. CTP equation and branch convention

For four distinct spacetime events, define

\[
F_{4,c}^{a_1a_2a_3a_4}(1,2,3,4),\qquad a_i\in\{+,-\},
\]

with the + branch represented by chronologically ordered ket insertions and the - branch by the complementary anti-chronological bra. The composite operator at event 1 is

\[
A_1=:\!\phi_1^3\!:\,_m.
\]

The connected source is the fourth cumulant with \(A_1\) treated as one composite argument:

\[
\mathcal S_{4,c}^{a_1a_2a_3a_4}
=M_{33}^{a_1a_2a_3a_4}
-\sum_{(1j)(kl)}M_{31}^{a_1a_j}(1;j)F_2^{a_ka_l}(k,l).
\]

Away from contour contact surfaces, the branch metric from the CTP action cancels between the equation-of-motion insertion and the branch-resolved correlator identity. Thus the noncontact equation has the same form for each branch assignment:

\[
\boxed{
(\Box_1+m^2)F_{4,c}^{a_1a_2a_3a_4}
+\frac{\lambda(t_1)}6\mathcal S_{4,c}^{a_1a_2a_3a_4}=0.
}
\]

The differentiated event is the latest physical-time event and all four times are distinct. No contact terms are included in this test. At coincident times, the contour delta-function terms must be restored.

## 2. Numerical implementation

The implementation uses the finite-volume Hamiltonian and the same projected normal-ordered cubic operator as the all-forward SD calculation. For each subset of events, a shared-prefix chronological state is built. The CTP component is the overlap of the complementary branch states. The same construction is repeated with event 1 replaced by the projected cubic insertion, and by the spectral second-derivative field insertion.

The finite-mode model is 1+1D, \(m=1\), \(L=20\), \(K=1\), \(\lambda(t)=0.1\operatorname{sech}^2(t)\), and the initial state is the free Gaussian vacuum. The evolution interval is \([-12,12]\), DOP853 uses \(rtol=2\times10^{-9}\), \(atol=2\times10^{-11}\), and \(max\_step=0.02\). The primary time stencil is \(h=0.003125\).

## 3. All 16 components at configuration A

\[
(1,0),\quad(0.5,1),\quad(-0.5,2),\quad(-1,3).
\]

At \(N_{\max}=8\), all 16 components have relative SD residuals below \(6.13\times10^{-6}\). At \(N_{\max}=10\), the maximum is \(6.09\times10^{-6}\).

| CTP component | Relative residual, \(N_{\max}=8\) |
|---|---:|
| ++++ | \(6.12712\times10^{-6}\) |
| +++- | \(1.78479\times10^{-6}\) |
| ++-+ | \(1.83856\times10^{-6}\) |
| ++-- | \(2.01978\times10^{-6}\) |
| +-++ | \(1.98173\times10^{-6}\) |
| +-+- | \(1.64648\times10^{-6}\) |
| +--+ | \(1.81643\times10^{-6}\) |
| +--- | \(6.09629\times10^{-6}\) |
| -+++ | \(6.09629\times10^{-6}\) |
| -++- | \(1.81643\times10^{-6}\) |
| -+-+ | \(1.64648\times10^{-6}\) |
| -+-- | \(1.98173\times10^{-6}\) |
| --++ | \(2.01978\times10^{-6}\) |
| --+- | \(1.83856\times10^{-6}\) |
| ---+ | \(1.78479\times10^{-6}\) |
| ---- | \(6.12712\times10^{-6}\) |

The all-forward ++++ result agrees with the independent time-ordered implementation: \(F_{4,c}\), \(M_{33}\), and the connected source agree at stored precision; the residual difference is at the \(10^{-11}\) relative level due to independent finite-difference evaluation.

## 4. Second spacetime configuration

\[
(1.2,0),\quad(0.4,0.7),\quad(-0.3,1.7),\quad(-1.1,3.2).
\]

All 16 components were computed at \(N_{\max}=8\). The maximum relative SD residual is \(6.29661\times10^{-6}\), and the minimum is \(2.07802\times10^{-6}\).

## 5. CTP identities and Fock convergence

For configuration A:

- Maximum largest-time branch-independence residual of \(F_{4,c}\): \(3.78\times10^{-17}\) absolute.
- Complement-branch conjugation residual of \(F_{4,c}\): zero at stored floating-point precision.
- Complement-branch conjugation residual of the SD residual: zero at stored floating-point precision.
- Largest-time branch difference of the finite-difference SD residual: \(4.05\times10^{-12}\) absolute at \(h=0.003125\).
- Maximum relative change from \(N_{\max}=8\) to 10 across all 16 components: \(1.56\times10^{-8}\) for \(F_{4,c}\), \(1.67\times10^{-9}\) for \(M_{33}\), and \(1.66\times10^{-9}\) for the connected source.

The residual's largest-time branch difference is an absolute finite-difference/integration diagnostic; it is not used as a relative residual because the equation residual itself is already small.

The free-vacuum test gives \(|F_{4,c}|<2\times10^{-17}\) for all branches. The absolute SD residual is below \(2.7\times10^{-12}\); a relative residual is not meaningful when both terms vanish with \(\lambda_0=0\).

## 6. Momentum-cutoff scan

An all-forward scan was performed at fixed \(N_{\max}=5\), \(m=1\), \(L=20\), \(\lambda_0=0.1\), \(T=1\), and \(h=0.00625\). The table reports complex magnitudes; the raw complex values are retained in JSON.

| \(K\) | Basis dimension | \(|F_{4,c}|\) | \(|M_{33}|\) | \(|\mathcal S_{4,c}|\) | Relative SD residual |
|---:|---:|---:|---:|---:|---:|
| 1 | 56 | \(1.51968\times10^{-5}\) | \(1.31476\times10^{-3}\) | \(1.31528\times10^{-3}\) | \(2.04173\times10^{-5}\) |
| 2 | 252 | \(3.74294\times10^{-5}\) | \(3.32906\times10^{-3}\) | \(3.33007\times10^{-3}\) | \(1.77474\times10^{-5}\) |
| 3 | 792 | \(4.49070\times10^{-5}\) | \(3.80564\times10^{-3}\) | \(3.80412\times10^{-3}\) | \(1.29344\times10^{-5}\) |

The observables change substantially as \(K\) increases: relative to \(K=1\), the \(K=3\) values change by about 196% for \(|F_{4,c}|\) and 192% for \(|\mathcal S_{4,c}|\). Therefore the momentum-mode cutoff is not converged and these data do not support a continuum claim.

Fock-cutoff checks show why matched \(N_{\max}\) is essential. At \(K=2\), changing \(N_{\max}=5\to6\) changes \(F_{4,c}\) by about \(0.108\%\) and the source by \(0.00172\%\). At \(K=3\), changing \(N_{\max}=4\to5\) changes \(F_{4,c}\) by about \(0.248\%\) and the source by about \(0.105\%\). A \(K=3,N_{\max}=6\) run is not yet available, so the \(K=3\) Fock convergence remains less controlled than \(K=1,2\).

## 7. Automated validation

The combined repository test suite passes 20/20 tests, covering:

- all 16 branch-resolved SD residuals at two event configurations;
- largest-time identity and complement-branch conjugation;
- independent ++++ comparison;
- \(N_{\max}=8\to10\) convergence for all CTP components;
- free-vacuum connected four-point limit;
- finite and archived \(K=1,2,3\) all-forward scan results;
- previous Phase-IV/V operator, spatial-translation, time-stencil, and Fock-cutoff checks.

## 8. Scope and next work

The all-16 CTP source/residual check is complete at \(K=1\) for two noncoincident event configurations. The all-forward momentum scan reaches \(K=3\), but shows strong cutoff dependence. The next numerical priority is a controlled matched-\(N_{\max}\) and \(K\) scan, especially \(K=3,N_{\max}\ge6\), followed by additional CTP configurations at larger \(K\). No continuum or infinite-volume extrapolation is claimed.