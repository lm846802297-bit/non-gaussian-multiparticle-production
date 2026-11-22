# Exact CTP Four-Point Hierarchy and Finite-Volume Real-Time Dynamics

Status: exact equations established; finite-dimensional real-time Hamiltonian benchmark completed. The calculation is not a continuum-limit result.

## 1. CTP conventions and regulated model

The model is a real scalar field in 1+1 dimensions, in natural units, with

\[
\mathcal L=\frac12(\partial\phi)^2-\frac12m^2\phi^2
-\frac{\lambda(t)}{4!}:\!\phi^4\!:\,_m,
\qquad
\lambda(t)=\lambda_0\operatorname{sech}^2(t/T).
\]

Normal ordering is relative to the free mass-\(m\) vacuum with the same finite momentum-mode regulator used in the Hamiltonian calculation. The initial state is that free Gaussian vacuum. The coupling tends to zero at both temporal ends, so free in/out operators exist.

Let contour labels be implicit in the arguments, e.g. \(1=(x_1,t_1,a_1)\), and define contour-ordered moments without an extra factor of \(-i\):

\[
F_2(1,2)=\langle T_C\phi_1\phi_2\rangle,
\qquad
F_4(1,2,3,4)=\langle T_C\phi_1\phi_2\phi_3\phi_4\rangle.
\]

For the \(\mathbb Z_2\)-symmetric zero-mean state,

\[
F_{4,c}=F_4-F_2(1,2)F_2(3,4)
-F_2(1,3)F_2(2,4)-F_2(1,4)F_2(2,3).
\]

The contour delta is defined with the standard closed-time-path orientation. With the above propagator convention, the free equation is \((\Box+m^2)F_2^{(0)}=-i\delta_C\).

## 2. Exact two-point Schwinger-Dyson equation

The operator equation of motion is

\[
(\Box_x+m^2)\phi(x)+\frac{\lambda(t)}{3!}:\!\phi^3(x)\!:\,_m=0.
\]

Applying it to the first leg of the contour-ordered two-point function gives

\[
(\Box_1+m^2)F_2(1,2)
+\frac{\lambda_1}{6}M_{3;1}(1;2)
=-i\delta_C(1-2),
\]

where

\[
M_{3;1}(1;2)=\langle T_C:\!\phi_1^3\!:\,_m\phi_2\rangle.
\]

Since \(:\!\phi^3\!:\,_m=\phi^3-3\Delta_m(0)\phi\), and the state has zero one-point function,

\[
M_{3;1}(1;2)=F_{4,c}(1,1,1,2)
+3\big[F_2(1,1)-\Delta_m(0)\big]F_2(1,2).
\]

Thus the exact two-point equation already depends on the connected four-point function. Hartree/Gaussian closure amounts to discarding the first term and replacing the coincident two-point function by its self-consistent local value.

## 3. Exact connected four-point equation

The full four-point Schwinger-Dyson equation is

\[
\begin{aligned}
(\Box_1+m^2)F_4(1,2,3,4)
+\frac{\lambda_1}{6}M_{3;3}(1;2,3,4)
={}&-i\delta_C(1-2)F_2(3,4)\\
&-i\delta_C(1-3)F_2(2,4)\\
&-i\delta_C(1-4)F_2(2,3),
\end{aligned}
\]

with

\[
M_{3;3}(1;2,3,4)=
\langle T_C:\!\phi_1^3\!:\,_m\phi_2\phi_3\phi_4\rangle.
\]

Subtracting the three disconnected pair products and using the two-point equation cancels the contour contact terms. The exact connected equation is therefore

\[
\boxed{
\begin{aligned}
(\Box_1+m^2)F_{4,c}(1,2,3,4)
+\frac{\lambda_1}{6}\Big[&M_{3;3}(1;2,3,4)\\
&-M_{3;1}(1;2)F_2(3,4)\\
&-M_{3;1}(1;3)F_2(2,4)\\
&-M_{3;1}(1;4)F_2(2,3)\Big]=0.
\end{aligned}
}
\]

This is an exact CTP identity, not a closure. The composite six-field moment contains the connected six-point function and lower cumulants. In a zero-mean \(\mathbb Z_2\)-symmetric state its ordinary six-field moment decomposes into the connected six-point cumulant, the 15 partitions of type \(4+2\), and the 15 pair partitions of type \(2+2+2\). The normal-ordering subtraction \(-3\Delta_m(0)F_4\) must then be applied to the full moment.

Consequently, the connected four-point equation couples to six-point information. A Gaussian ansatz sets all connected cumulants above second order to zero; a four-point truncation must specify how the six-point sector is approximated. A 2PI calculation with a Bethe-Salpeter response equation and a 4PI truncation are distinct approximations and must not be conflated.

## 4. Finite-volume Hamiltonian realization

For periodic length \(L\), retain modes \(p_n=2\pi n/L\), \(|n|\le K\), with \(\omega_n=\sqrt{p_n^2+m^2}\). In the free Fock basis,

\[
H(t)=H_0+\lambda(t)V,
\qquad
V=\frac1{4!}\int_0^L dx:\!\phi^4(x)\!:\,_m.
\]

The interaction matrix is built by explicitly summing all normal-ordered creation/annihilation monomials. Spatial homogeneity makes total momentum an exact block label; the vacuum evolves only in the zero-total-momentum, even-\(\mathbb Z_2\) sector. The numerical state obeys

\[
i\partial_t|\Psi(t)\rangle=H(t)|\Psi(t)\rangle.
\]

This is a nonperturbative real-time evolution inside the chosen finite momentum and Fock-space truncations. It retains all powers of \(\lambda_0\) represented within that projected Hilbert space. The Fock cutoff \(N_{\max}\) is checked independently from the momentum-mode cutoff \(K\).

## 5. Out-mode observables and Gaussian reconstruction

After the pulse, free out operators are exact. For each measured channel,

\[
N_1(k)=\langle n_k\rangle,
\qquad
N_2(k,q)=\begin{cases}
\langle n_kn_q\rangle,&k\ne q,\\
\langle n_k(n_k-1)\rangle,&k=q.
\end{cases}
\]

The Gaussian reconstruction uses all final two-point contractions:

\[
N_2^G(k,q)=\begin{cases}
N_1(k)N_1(q)+|\langle a_k^\dagger a_q\rangle|^2
+|\langle a_ka_q\rangle|^2,&k\ne q,\\
2N_1(k)^2+|\langle a_ka_k\rangle|^2,&k=q.
\end{cases}
\]

The connected number cumulant is \(\kappa_4=N_2-N_2^G\), and \(\mathcal R_2=\kappa_4/N_2^G\) where the denominator is nonzero. The reported values are evaluated from the full evolved state, not from a Gaussian approximation.

## 6. Equal-time spatial four-point function

In addition to the out-mode projection, the finite-volume state is used to compute equal-time spatial slices

\[
F_{4,c}(t;0,0,0,r)=
\langle\phi(t,0)^3\phi(t,r)\rangle
-3\langle\phi(t,0)^2\rangle\langle\phi(t,0)\phi(t,r)\rangle,
\]

where the three equivalent pairings coincide for the selected geometry and translation-invariant state. More generally the implementation subtracts all three pairings explicitly. The field-operator evaluation basis is extended by four particles so that the four successive field insertions are not cut off at the state-vector boundary.

These are equal-time spatial slices of the CTP connected four-point function. They are not yet the full independent-four-time object \(F_{4,c}(x_1,x_2;y_1,y_2)\); that requires unequal-time contour propagation and is a subsequent extension.

## 7. Numerical findings for the stated finite model

Parameters: \(m=1\), \(L=20\), \(T=1\), \(K=2\), \(N_{\max}=10\), unless stated otherwise.

### Coupling scan, generic channel \((n_k,n_q)=(1,2)\)

| \(\lambda_0/m^2\) | Exact \(N_2\) | \(N_2/N_2^{(2)}\) | \(\mathcal R_2\) |
|---:|---:|---:|---:|
| 0.02 | \(2.61648\times10^{-11}\) | 0.98954 | \(2.2070\times10^9\) |
| 0.05 | \(1.60974\times10^{-10}\) | 0.97407 | \(3.5878\times10^8\) |
| 0.10 | \(6.27249\times10^{-10}\) | 0.94889 | \(9.2044\times10^7\) |
| 0.20 | \(2.38144\times10^{-9}\) | 0.90064 | \(2.4172\times10^7\) |
| 0.30 | \(5.08713\times10^{-9}\) | 0.85508 | \(1.1245\times10^7\) |
| 0.50 | \(1.27459\times10^{-8}\) | 0.77127 | \(4.3832\times10^6\) |

The leading perturbative result overestimates the exact finite-Hilbert-space result by about 1.0%, 2.6%, 5.1%, 9.9%, 14.5%, and 22.9% across this scan. For this vacuum-pulse setup, the Gaussian reconstruction is parametrically much smaller than the full pair spectrum: \(N_2=O(\lambda_0^2)\), while the products of final two-point occupations begin at higher order. The opposite-momentum channel also has an anomalous two-point contraction, which is included in \(N_2^G\).

### Fock-space cutoff convergence

At \(\lambda_0/m^2=0.5\), \(K=2\), the generic-channel \(N_2\) values for \(N_{\max}=4,6,8,10,12\) are respectively

\[
1.25758086,\ 1.27288883,\ 1.27457172,\ 1.27458764,\ 1.27458779
\quad\times10^{-8}.
\]

The change from \(N_{\max}=10\) to 12 is about \(1.1\times10^{-7}\) relative. This demonstrates convergence in the particle-number truncation for this channel and parameter point.

### Momentum-mode cutoff

At \(\lambda_0/m^2=0.2\), \(N_{\max}=6\), the generic-channel values for \(K=2,3,4,5,6,7\) are

\[
2.38085,\ 3.69849,\ 4.11780,\ 4.17485,\ 4.15220,\ 4.12533
\quad\times10^{-9}.
\]

The K=6 to K=7 change is about 0.65%; the values are not monotonic and the remaining variation is not negligible enough to claim momentum-cutoff convergence. This is an explicit outstanding systematic uncertainty; no continuum result is claimed.

### Pulse-duration dependence

At \(\lambda_0/m^2=0.2\), \(K=2\), \(N_{\max}=10\), the generic-channel \(N_2/N_2^{(2)}\) ratio decreases from 0.9748 at \(mT=0.5\) to 0.7871 at \(mT=1.5\). The pulse duration changes both the absolute production and the relative importance of higher-order processes.

### Time-resolved connected correlations

For \(\lambda_0/m^2=0.2\), \(mT=1\), \(K=2\), \(N_{\max}=8\), the generic-channel reference-mode cumulant \(\kappa_4(t)\) is approximately \(2.66\times10^{-7}\) at the pulse maximum \(t=0\), and evolves to \(2.38\times10^{-9}\) at \(t=12T\). The equal-time spatial connected correlator \(F_{4,c}(t;0,0,0,2)\) is approximately \(-3.6571\times10^{-5}\) at \(t=0\), changes sign during the post-pulse evolution, and is approximately \(-3.5918\times10^{-7}\) at \(t=12T\). This displays nontrivial real-time and spatial structure, even though the asymptotic out-mode spectrum is much smaller than the peak instantaneous reference-mode cumulant.

## 8. Validation and limits

- The projected interaction matrix is Hermitian to numerical precision.
- State norm is conserved to approximately \(10^{-12}\) or better in the reported runs.
- The time-step scan with step sizes reduced by factors of two and four gives indistinguishable displayed \(N_2\) values.
- Five unit tests pass, including the normalized four-particle vertex coefficient, zero-coupling limit, weak-coupling agreement with the leading perturbative result, Hermiticity, and coincident-mode factorial moment.
- Fock-space cutoff convergence is demonstrated at \(K=2\); momentum-mode cutoff convergence remains open.
- The calculation is a finite-mode, finite-Fock-space Hamiltonian simulation. It does not yet include a continuum extrapolation, infinite-volume extrapolation, or full unequal-time CTP four-point reconstruction.