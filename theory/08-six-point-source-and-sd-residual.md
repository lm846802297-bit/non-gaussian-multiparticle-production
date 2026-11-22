# Six-Field Source and Numerical Schwinger-Dyson Residual

Status: Phase-IV derivation and implementation specification. No numerical residual is claimed in this document.

## 1. Objective

The previous phase computed all 16 CTP components of the connected four-field correlator at selected spacetime configurations. The next check is to test the exact connected four-point Schwinger-Dyson equation against the finite-dimensional Hamiltonian evolution, and to isolate the six-field composite source that drives the hierarchy.

The first implementation target is the all-forward (time-ordered) component, at configurations where the differentiated event is not coincident in time with the other three events. This avoids equal-time contact terms in the first residual test. The other 15 contour assignments can then be added using the same branch-state construction.

## 2. Regulated equation of motion

For the normal-ordered theory

\[
\mathcal L=\frac12(\partial\phi)^2-\frac12m^2\phi^2
-\frac{\lambda(t)}{4!}:\!\phi^4\!:\,_m,
\]

the Heisenberg equation is

\[
(\Box_1+m^2)\phi_1
=-\frac{\lambda(t_1)}{3!}:\!\phi_1^3\!:\,_m.
\]

Define the time-ordered two-point function

\[
F_2(1,2)=\langle T\phi_1\phi_2\rangle
\]

and the four-point function

\[
F_4(1,2,3,4)=\langle T\phi_1\phi_2\phi_3\phi_4\rangle.
\]

For mutually noncoincident time arguments, differentiating the time-ordered products produces no delta-function contact terms. The two-point equation is

\[
(\Box_1+m^2)F_2(1,2)
=-\frac{\lambda_1}{6}M_{31}(1;2),
\]

where

\[
M_{31}(1;2)=
\langle T:\!\phi_1^3\!:\,_m\phi_2\rangle.
\]

Similarly,

\[
(\Box_1+m^2)F_4(1,2,3,4)
=-\frac{\lambda_1}{6}M_{33}(1;2,3,4),
\]

with

\[
M_{33}(1;2,3,4)=
\langle T:\!\phi_1^3\!:\,_m\phi_2\phi_3\phi_4\rangle.
\]

## 3. Exact connected four-point equation

Subtract the three disconnected pairings,

\[
\begin{aligned}
F_{4,c}={}&F_4-F_2(1,2)F_2(3,4)\\
&-F_2(1,3)F_2(2,4)
-F_2(1,4)F_2(2,3).
\end{aligned}
\]

Applying \(\Box_1+m^2\) and using the two-point equation gives

\[
\boxed{
\begin{aligned}
(\Box_1+m^2)F_{4,c}(1,2,3,4)
+\frac{\lambda_1}{6}\,\mathcal S_{4,c}(1;2,3,4)=0,
\end{aligned}
}
\]

where the exact connected source is

\[
\boxed{
\begin{aligned}
\mathcal S_{4,c}={}&M_{33}(1;2,3,4)\\
&-M_{31}(1;2)F_2(3,4)\\
&-M_{31}(1;3)F_2(2,4)\\
&-M_{31}(1;4)F_2(2,3).
\end{aligned}
}
\]

This is not a closure: \(\mathcal S_{4,c}\) contains six-field information. The equation is exact within the regulated Hamiltonian theory, away from the contact surfaces. At coincident times, the contour/time-ordering delta terms must be restored explicitly.

## 4. Normal-ordering subtraction in the finite-mode regulator

For the same finite momentum regulator used in the Hamiltonian, define

\[
\Delta_m(0)=\langle0_m|\phi^2(x)|0_m\rangle
=\frac{1}{L}\sum_{|n|\le K}\frac{1}{2\omega_n}.
\]

Then

\[
:\!\phi^3\!:\,_m=\phi^3-3\Delta_m(0)\phi.
\]

Consequently, the six-field source must be evaluated as

\[
\begin{aligned}
M_{33}(1;2,3,4)
={}&\langle T\phi_1\phi_1\phi_1\phi_2\phi_3\phi_4\rangle\\
&-3\Delta_m(0)
\langle T\phi_1\phi_2\phi_3\phi_4\rangle.
\end{aligned}
\]

Likewise,

\[
M_{31}(1;2)=
\langle T\phi_1\phi_1\phi_1\phi_2\rangle
-3\Delta_m(0)F_2(1,2).
\]

The same finite-mode regulator and the same mass-\(m\) vacuum subtraction must be used in the Hamiltonian, the composite operators, and the correlators. Mixing regulators would invalidate the residual test.

## 5. Numerical residual

Choose event 1 to be the latest physical-time event and keep \(t_1\) separated from \(t_2,t_3,t_4\). For the all-forward component, define

\[
\mathscr E_{4,c}(1;2,3,4)
=(\Box_1+m^2)F_{4,c}
+\frac{\lambda_1}{6}\mathcal S_{4,c}.
\]

The target is

\[
\mathscr E_{4,c}\to0
\]

as the time-derivative and spatial-derivative discretizations are refined and the finite Hilbert-space calculation is converged.

A first finite-difference estimator is

\[
\begin{aligned}
\Box_1F_{4,c}\simeq{}&
\frac{F_{4,c}(t_1+h,x_1)-2F_{4,c}(t_1,x_1)
+F_{4,c}(t_1-h,x_1)}{h^2}\\
&-\frac{F_{4,c}(t_1,x_1+a)-2F_{4,c}(t_1,x_1)
+F_{4,c}(t_1,x_1-a)}{a^2}.
\end{aligned}
\]

The implementation must scan \(h\) and \(a\), compare second- and fourth-order stencils where possible, and report the residual relative to the magnitudes of its two terms:

\[
\epsilon_{\rm SD}=
\frac{|\mathscr E_{4,c}|}
{|(\Box_1+m^2)F_{4,c}|+|\lambda_1\mathcal S_{4,c}|/6+\epsilon_{\rm floor}}.
\]

A residual evaluated at only one stencil spacing is not a validation.

## 6. Branch-state evaluation of the six-field term

The existing branch-state method generalizes from four to six field insertions. For a set of events \(e_i=(t_i,x_i)\), a branch subset \(S\) defines a forward-evolved ket with the selected fields inserted in chronological order. An all-forward time-ordered correlator is obtained from the overlap of the fully inserted branch ket with the uninserted evolved state.

For \(M_{33}\), the event list contains three copies of \((t_1,x_1)\), followed by \((t_2,x_2),(t_3,x_3),(t_4,x_4)\). The repeated insertion is an operator product at the same spacetime point; the normal-ordering subtraction above must be applied after evaluating the regulated six-field expectation.

For a full CTP source, branch assignments must be retained for each insertion, including the three coincident fields at point 1. The first validation milestone is the all-forward equation; only after it passes should the complete contour-index source be assembled.

## 7. Required checks

1. Free-vacuum limit: \(F_{4,c}=0\), \(M_{31}\) and \(M_{33}\) cancel under the normal-ordering prescription, and the residual vanishes.
2. Interaction-generated source: verify \(\mathcal S_{4,c}\) is nonzero when \(\lambda(t_1)\ne0\).
3. Normal-ordering consistency: independently evaluate \(\Delta_m(0)\) from the regulated mode sum and from the free-vacuum field matrix.
4. Fock-space convergence: compare at least \(N_{\max}=12,16\), and extend if needed.
5. Momentum-mode convergence: repeat for increasing \(K\); the previous phase showed that this is not yet controlled.
6. Derivative convergence: scan \(h\), \(a\), and stencil order.
7. Independent identity: compare the source equation against direct finite-dimensional Heisenberg evolution of the field operator.
8. Contact-term check: repeat with event times approaching one another only after deriving and including the corresponding delta-function terms.
9. Provenance: archive source hashes, raw correlators, all parameters, runtime/hardware, and numerical residual scans.

## 8. Scope and status

This document specifies the exact equation and the numerical validation protocol. It does not claim that the six-field source or the Schwinger-Dyson residual has already been numerically evaluated. The next calculation is the all-forward source/residual at one noncoincident event configuration, followed by regulator and stencil scans.
