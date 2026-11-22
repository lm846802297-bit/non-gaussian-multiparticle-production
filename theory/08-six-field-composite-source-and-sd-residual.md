# Composite Six-Field Source and a Finite-Projection CTP Schwinger-Dyson Residual

Status: the composite source operator and its connected subtraction are implemented and independently checked. A finite-difference residual diagnostic has been run in several finite bases; it is not yet a converged test of the continuum connected four-point equation.

## 1. Exact connected four-point equation

For the normal-ordered interaction

\[
\mathcal L_{\rm int}=-\frac{\lambda(t)}{4!}:\!\phi^4\!:\,_m,
\]

the field equation is

\[
(\Box_1+m^2)\phi_1+\frac{\lambda_1}{3!}:\!\phi_1^3\!:\,_m=0.
\]

Applying this to the four-point contour correlator and subtracting its three two-point pairings gives, at noncoincident external points,

\[
\boxed{
(\Box_1+m^2)F_{4,c}(1,2,3,4)
+\frac{\lambda_1}{3!}\,\mathcal K_{3;3}(1;2,3,4)=0,
}
\]

where

\[
\begin{aligned}
\mathcal K_{3;3}(1;2,3,4)={}&
\langle T_C :\!\phi_1^3\!:\,\phi_2\phi_3\phi_4\rangle\\
&-\langle T_C :\!\phi_1^3\!:\phi_2\rangle F_2(3,4)\\
&-\langle T_C :\!\phi_1^3\!:\phi_3\rangle F_2(2,4)\\
&-\langle T_C :\!\phi_1^3\!:\phi_4\rangle F_2(2,3).
\end{aligned}
\]

For a \(\mathbb Z_2\)-symmetric state, the odd one-point and three-point functions vanish, so this is the fourth cumulant with the *composite operator* \(A_1=:\phi_1^3:\) treated as one argument and \(\phi_2,\phi_3,\phi_4\) as the other three.

A crucial distinction: \(\mathcal K_{3;3}\) is not simply the elementary connected six-field correlator \(F_{6,c}\). The three constituent fields inside \(:\phi_1^3:\) can each connect to a different external field. Thus even in a Gaussian state, \(\mathcal K_{3;3}\) need not vanish.

## 2. Gaussian identity for the composite source

For a free Gaussian vacuum and normal ordering with respect to that same free vacuum,

\[
\langle :\!\phi_1^3\!:\,\phi_2\phi_3\phi_4\rangle
=3!\,F_2(1,2)F_2(1,3)F_2(1,4),
\]

while

\[
\langle :\!\phi_1^3\!: \phi_j\rangle=0.
\]

Therefore

\[
\boxed{\mathcal K_{3;3}^{\rm Gaussian}=6F_{12}F_{13}F_{14}.}
\]

This is a connected cumulant of the composite variable \(A_1\), even though the elementary six-point cumulant of a Gaussian field vanishes. This distinction is essential when constructing the source term for the four-point hierarchy.

## 3. Finite-basis composite operator

The code constructs the regulated normal-ordered operator directly from the retained free modes:

\[
:\!\Phi(x)^3\!:\,_m
=\sum_{c+a=3}\frac{3!}{\prod_p c_p!\prod_p a_p!}
\left[\prod_p u_p^{c_p+a_p}\right]
 e^{i(\sum_p a_pp-\sum_p c_pp)x}
\prod_p(a_p^\dagger)^{c_p}\prod_p a_p^{a},
\]

with \(u_p=(2\omega_pL)^{-1/2}\). An independent matrix identity check verifies, on interior Fock columns,

\[
:\!\Phi^3\!:\,_m=\Phi^3-3\Delta_m(0)\Phi,
\qquad
\Delta_m(0)=\sum_{|p|\le K_{\rm cut}}\frac{1}{2\omega_pL}.
\]

The test residual is below \(2\times10^{-14}\) for the tested finite bases, and the composite matrix is Hermitian to numerical precision.

## 4. Numerical CTP residual diagnostic

For the all-forward component, the code evaluates \(F_{4,c}^{++++}\) and \(\mathcal K_{3;3}^{++++}\) by the same branch-state overlap method used for the all-16 CTP calculation. It approximates

\[
\Box_1F_{4,c}\simeq
\frac{F(t_1+h)-2F(t_1)+F(t_1-h)}{h^2}
-
\frac{F(x_1+h)-2F(x_1)+F(x_1-h)}{h^2}.
\]

The diagnostic is

\[
\mathcal E_4=(\Box_1+m^2)F_{4,c}+\frac{\lambda(t_1)}{6}\mathcal K_{3;3}.
\]

At \(m=1,L=20,K_{\rm cut}=1,\lambda_0=0.2,T=1\), with event points

\[
((1,0),(0.5,1),(-0.5,2),(-1,3)),
\]

the relative residual (normalized to the larger of the two terms in \(\mathcal E_4\)) is approximately 3.1–3.3% for \(N_{\max}=4,6,8\) and \(h\) between 0.0125 and 0.1. The \(N_{\max}=6\) residual changes only mildly over this finite-difference scan.

A low-cutoff \(K_{\rm cut}=2,N_{\max}=4\) run gives a much larger residual (about 18%), showing that this diagnostic is sensitive to the momentum/Fock regulator combination. This point is not Fock-converged and is not used as a physical conclusion.

## 5. Interpretation

The composite source is now explicitly constructed and its Gaussian Wick structure is validated. The finite-basis equation residual is not yet zero and is not being reported as a successful exact-equation verification. The likely issue to resolve is the relation between the continuum Schwinger-Dyson equation and the projected finite-Fock operator algebra: truncating the Hilbert space modifies canonical commutators and introduces boundary terms in the projected equation. Momentum-mode truncation and finite-difference errors also need to be separated.

The next calculation should derive the *projected* Heisenberg equation, including Fock-boundary terms, and then compare it with the continuum equation as \(N_{\max}\) and \(K_{\rm cut}\) increase. A direct residual should only be interpreted as a continuum check after these regulator effects are controlled.

No continuum result or final novelty claim is made.
