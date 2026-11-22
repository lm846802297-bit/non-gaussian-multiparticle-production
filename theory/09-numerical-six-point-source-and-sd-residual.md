# Numerical six-field source and connected four-point Schwinger-Dyson test

Status: Phase-IV numerical validation completed for one noncoincident all-forward configuration in the finite-mode Hamiltonian model. This is a regulated finite-volume result, not a continuum-limit result.

## 1. Model and event configuration

We use the normal-ordered real scalar theory in 1+1 dimensions,

\[
\mathcal L=\frac12(\partial\phi)^2-\frac12m^2\phi^2
-\frac{\lambda(t)}{4!}:\!\phi^4\!:\,_m,
\qquad \lambda(t)=\lambda_0\operatorname{sech}^2(t/T).
\]

The parameters are

\[
m=1,\quad L=20,\quad K=1,\quad \lambda_0=0.1,\quad T=1,
\]

with momentum modes \(n=-1,0,1\), periodic boundary conditions, and a free Gaussian vacuum at \(t_0=-12T\). The final reference time is \(t_f=12T\). The four spacetime events are

\[
1=(1,0),\quad 2=(0.5,1),\quad 3=(-0.5,2),\quad 4=(-1,3).
\]

Event 1 is the latest physical-time event. All four times are distinct, so the time-ordered equation below has no coincident-time contact terms.

## 2. Exact equation tested

Define

\[
F_2(i,j)=\langle T\phi_i\phi_j\rangle,
\]

\[
F_{4,c}=\langle T\phi_1\phi_2\phi_3\phi_4\rangle
-F_2(1,2)F_2(3,4)-F_2(1,3)F_2(2,4)-F_2(1,4)F_2(2,3).
\]

For the normal-ordered interaction, define the composite moments

\[
M_{31}(1;j)=\langle T:\!\phi_1^3\!:\,_m\phi_j\rangle,
\]

\[
M_{33}(1;2,3,4)=\langle T:\!\phi_1^3\!:\,_m\phi_2\phi_3\phi_4\rangle.
\]

Then, away from contact surfaces,

\[
(\Box_1+m^2)F_{4,c}+\frac{\lambda(t_1)}6\,\mathcal S_{4,c}=0,
\]

where

\[
\begin{aligned}
\mathcal S_{4,c}={}&M_{33}(1;2,3,4)\\
&-M_{31}(1;2)F_2(3,4)\\
&-M_{31}(1;3)F_2(2,4)\\
&-M_{31}(1;4)F_2(2,3).
\end{aligned}
\]

The numerical residual is

\[
\mathscr E_{4,c}=(\Box_1+m^2)F_{4,c}+\frac{\lambda(t_1)}6\mathcal S_{4,c}.
\]

The source \(\mathcal S_{4,c}\) is the fourth cumulant with the composite operator \(A_1=:\!\phi_1^3\!:\,_m\) treated as one argument and \(\phi_2,\phi_3,\phi_4\) as the other three arguments. It is not the elementary connected six-field cumulant. In a Gaussian vacuum,

\[
\mathcal S_{4,c}=M_{33}=3!\,F_2(1,2)F_2(1,3)F_2(1,4),
\]

because \(M_{31}(1;j)=0\) while all three fields in the normal-ordered cubic insertion can contract with the three external fields. Thus \(M_{33}\) need not vanish in a Gaussian state; the Schwinger-Dyson source is multiplied by \(\lambda(t_1)\), and it is the composite cumulant rather than an elementary six-point connected function.

## 3. Finite-mode projection and normal ordering

A finite momentum cutoff changes the local equation into the projected equation

\[
(\Box+m^2)\phi_K(x)=-\frac{\lambda(t)}6P_K:\!\phi_K^3(x)\!:\,_m.
\]

Here \(P_K\) projects the cubic composite onto the retained Fourier modes. The calculation therefore uses the projected normal-ordered cubic operator, not an unprojected pointwise cube. This is required for consistency with the finite-mode Hamiltonian.

The finite-regulator tadpole is

\[
\Delta_m(0)=\frac1L\sum_{|n|\le K}\frac1{2\omega_n}
=0.07270141081892326.
\]

An independent Fourier-projection check verifies

\[
P_K\big[\phi_K^3-3\Delta_m(0)\phi_K\big]
=P_K:\!\phi_K^3\!:\,_m
\]

on the interior Fock subspace. A second independent check verifies the canonical equation-of-motion operator identity

\[
-6i[V,\pi_K(x)]=P_K:\!\phi_K^3(x)\!:\,_m
\]

for matrix elements that do not touch the Fock cutoff boundary.

## 4. Numerical method

The full finite-mode Fock basis is used for the field insertions. The projected Hamiltonian is

\[
H(t)=H_0+\lambda(t)V,
\]

and states are propagated with DOP853. Time-ordered moments are computed by inserting the relevant field or projected cubic operator at the event times along a forward-evolved Schrödinger state and taking its overlap with the uninserted evolved state.

The six-field source is evaluated directly as one projected normal-ordered cubic insertion at event 1 and three ordinary field insertions at events 2, 3, and 4. The three \(M_{31}F_2\) subtraction terms are computed independently from their corresponding moments.

For the time derivative, a centered second-order stencil is used:

\[
\partial_{t_1}^2F_{4,c}\simeq
\frac{F_{4,c}(t_1+h)-2F_{4,c}(t_1)+F_{4,c}(t_1-h)}{h^2}.
\]

The spatial second derivative is evaluated spectrally from the retained Fourier modes. Finite-difference spatial stencils are independently compared with this spectral derivative for \(a=0.05,0.1,0.2\).

## 5. Main numerical result

At \(N_{\max}=8\), \(h=0.00625\), and \(a=0.2\),

\[
F_{4,c}=(3.6008526813\times10^{-6})
+i(1.4752437347\times10^{-5}),
\]

\[
M_{33}= -7.3597883011\times10^{-4}
+i(1.0894708180\times10^{-3}),
\]

\[
\mathcal S_{4,c}= -7.3682396732\times10^{-4}
+i(1.0895287487\times10^{-3}).
\]

The two terms in the Schwinger-Dyson equation are

\[
(\Box_1+m^2)F_{4,c}
=5.1575473292\times10^{-6}
-i(7.6257972911\times10^{-6}),
\]

\[
\frac{\lambda(t_1)}6\mathcal S_{4,c}
=-5.1574526760\times10^{-6}
+i(7.6262353148\times10^{-6}).
\]

Their sum is

\[
\mathscr E_{4,c}=(9.46532275\times10^{-11})
+i(4.38023756\times10^{-10}),
\]

with

\[
|\mathscr E_{4,c}|=4.48134\times10^{-10},
\qquad
\epsilon_{SD}=2.43384\times10^{-5}.
\]

Reducing the time stencil from \(h=0.00625\) to \(h=0.003125\) gives

\[
\epsilon_{SD}=6.12712\times10^{-6}.
\]

The residual follows the expected \(h^2\) behavior over the tested range. A linear fit of the complex residual against \(h^2\) gives a zero-step intercept with magnitude \(8.86\times10^{-13}\). This extrapolated value is a discretization diagnostic, not a separate direct simulation result.

## 6. Convergence and independent checks

### Fock cutoff

At \(h=0.00625\), \(a=0.2\), comparing \(N_{\max}=8\) and 10 gives relative changes:

| Quantity | Relative change |
|---|---:|
| \(F_{4,c}\) | \(2.96\times10^{-9}\) |
| \(M_{33}\) | \(6.66\times10^{-10}\) |
| \(\mathcal S_{4,c}\) | \(6.58\times10^{-10}\) |

The SD relative residual is \(2.43384\times10^{-5}\) at \(N_{\max}=8\) and \(2.43110\times10^{-5}\) at \(N_{\max}=10\).

### Time stencil

| \(h\) | Relative SD residual |
|---:|---:|
| 0.025 | \(3.89300\times10^{-4}\) |
| 0.0125 | \(9.73241\times10^{-5}\) |
| 0.00625 | \(2.43384\times10^{-5}\) |
| 0.003125 | \(6.12712\times10^{-6}\) |

### Spatial stencil

The finite-difference spatial derivative differs from the spectral derivative by:

| \(a\) | Relative difference |
|---:|---:|
| 0.05 | \(2.05666\times10^{-5}\) |
| 0.1 | \(8.22456\times10^{-5}\) |
| 0.2 | \(3.28944\times10^{-4}\) |

### Free-vacuum check

At \(\lambda_0=0\),

\[
|F_{4,c}|=3.88\times10^{-18},\qquad
|\mathscr E_{4,c}|=2.30\times10^{-15}.
\]

The normal-ordered six-field moment itself need not vanish in a Gaussian vacuum: it contains contractions of the three external fields with the three fields in the composite insertion. The interaction coefficient is zero in this check, so the equation residual vanishes. The validation therefore tests the connected four-point cancellation and the equation residual, not an incorrect assertion that \(M_{33}=0\).

### Operator and test checks

- Projected Hamiltonian Hermiticity residual: below \(6\times10^{-17}\).
- Out-state norm: within \(2\times10^{-15}\) of unity.
- Projected cubic versus commutator identity: relative matrix-element error below \(4\times10^{-16}\) on the interior Fock subspace.
- Fourier-projected \(\phi^3-3\Delta_m(0)\phi\) versus directly normal-ordered cubic: relative error below \(2\times10^{-12}\).
- Automated tests: 6/6 passed.

## 7. Scope and remaining work

This is a successful numerical test of one noncoincident all-forward connected four-point equation in a finite-volume, finite-momentum-mode Hamiltonian model. The Fock cutoff is controlled between 8 and 10 for the reported observables. The momentum cutoff remains \(K=1\); no continuum or infinite-volume extrapolation is claimed.

The next work is to repeat the source/residual test for other event configurations and contour assignments, increase \(K\), and then compare controlled 2PI/Bethe-Salpeter or 4PI closures against the finite-Hilbert benchmark.