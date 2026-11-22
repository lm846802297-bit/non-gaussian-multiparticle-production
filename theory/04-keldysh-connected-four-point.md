# Lowest-Order Keldysh Derivation of the Dynamically Generated Out-Mode Cumulant

Status: analytic derivation / benchmark in progress. No numerical result is claimed.

## 1. Scope and conventions

We work in \(1+1\) dimensions, in a periodic box of length \(L\), with \(\hbar=c=1\). The initial state at \(t\to-\infty\) is the free Gaussian vacuum of mass \(m>0\). The interaction is a real, spatially homogeneous pulse that vanishes at early and late times:

\[
\mathcal L_{\rm int}(t,x)=-\frac{\lambda(t)}{4!}:\!\phi^4(t,x)\!:\,_m,
\qquad \lambda(t)\to0\quad(t\to\pm\infty).
\]

Normal ordering is with respect to the free mass-\(m\) vacuum. The pulse may be the smooth profile in the model-and-asymptotic-modes document; the derivation below only assumes that its Fourier transform exists:

\[
\widetilde\lambda(\Omega)=\int_{-\infty}^{\infty}dt\,\lambda(t)e^{i\Omega t}.
\]

The free interaction-picture field is expanded as

\[
\phi_I(t,x)=\frac{1}{\sqrt L}\sum_p\frac{1}{\sqrt{2\omega_p}}
\left(a_p e^{-i\omega_pt+ipx}+a_p^\dagger e^{i\omega_pt-ipx}\right),
\qquad \omega_p=\sqrt{p^2+m^2},
\]

with \([a_p,a_q^\dagger]=\delta_{pq}\). All sums are over the discrete momenta \(p=2\pi n/L\).

The interaction-picture Hamiltonian is

\[
H_I(t)=\frac{\lambda(t)}{4!}\int_0^L dx:\!\phi_I^4(t,x)\!:\,_m.
\]

The exact out-mode operator is obtained by projecting the field in the free out region, as given in the model-and-asymptotic-modes document. In the interaction picture, the calculation below is equivalently an in-in expectation with the out number operator inserted after the pulse.

## 2. The observable and the relevant cumulant

For distinct modes \(k\ne q\), define

\[
N_2(k,q)=
\langle 0_{\rm in}|U^\dagger
a_k^\dagger a_q^\dagger a_q a_k
U|0_{\rm in}\rangle,
\qquad
U=T\exp\!\left[-i\int dt\,H_I(t)\right].
\]

The Gaussian reconstruction from the final two-point data is

\[
N_2^{\rm G}(k,q)=
\langle a_k^\dagger a_k\rangle\langle a_q^\dagger a_q\rangle
+\langle a_k^\dagger a_q\rangle\langle a_q^\dagger a_k\rangle
+\langle a_k^\dagger a_q^\dagger\rangle\langle a_q a_k\rangle.
\]

The number-type connected out-mode cumulant is

\[
\kappa_4(k,q)=N_2(k,q)-N_2^{\rm G}(k,q).
\]

This is an ordered cumulant of the operators \(a_k^\dagger,a_q^\dagger,a_q,a_k\). It is not interchangeable with the time-ordered four-field Green function or with the anomalous four-annihilator correlator \(\langle a_{p_1}a_{p_2}a_{p_3}a_{p_4}\rangle_c\). The latter can already be \(O(\lambda)\); the number-type cumulant studied here starts at \(O(\lambda^2)\) for the normal-ordered vacuum-pulse setup.

## 3. Keldysh/in-in expansion and the first nonzero order

Write \(O_{kq}=a_k^\dagger a_q^\dagger a_q a_k\). Since \(O_{kq}|0\rangle=0\) and \(\langle0|O_{kq}=0\), both the zeroth-order and first-order contributions vanish:

\[
N_2^{(0)}=0,
\qquad
N_2^{(1)}=i\int dt\,\langle0|[H_I(t),O_{kq}]|0\rangle=0.
\]

Expanding the forward and backward evolution operators,

\[
U=1-i\int dt\,H_I(t)+O(\lambda^2),
\qquad
U^\dagger=1+i\int dt\,H_I(t)+O(\lambda^2),
\]

the only second-order term that survives between the vacuum and this normal-ordered number operator is the cross-branch term:

\[
\boxed{
N_2^{(2)}(k,q)=
\int dt\,dt'\,
\langle0|H_I(t)\,O_{kq}\,H_I(t')|0\rangle.
}
\]

The two interaction vertices lie on opposite branches of the closed-time path: one on the backward branch and one on the forward branch. This is the cut, positive-probability contribution. Terms with both vertices on the same branch multiply \(O_{kq}|0\rangle\) or \(\langle0|O_{kq}\) and vanish at this order.

Equivalently, inserting a complete set of free out Fock states between the two vertices gives

\[
N_2^{(2)}(k,q)=
\sum_F
\langle0|H_I(t)|F\rangle
\langle F|O_{kq}|F\rangle
\langle F|H_I(t')|0\rangle
\quad\text{integrated over }t,t'.
\]

At first order in \(H_I\), a normal-ordered quartic vertex acting on the vacuum creates exactly four particles. Therefore only four-particle states contribute to this leading term.

## 4. First-order vacuum-to-four-particle amplitude

Let a normalized four-particle occupation state be

\[
|\{n_p\}\rangle=
\prod_p\frac{(a_p^\dagger)^{n_p}}{\sqrt{n_p!}}|0\rangle,
\qquad \sum_p n_p=4.
\]

The four-creation part of \(H_I(t)\), with the stated Fourier convention, gives

\[
\langle\{n_p\}|H_I(t)|0\rangle
=
\frac{\lambda(t)}{L}\,
\frac{\delta_{\sum_p n_p p,\,0}\,
e^{it\sum_p n_p\omega_p}}
{\sqrt{\prod_p n_p!}\,\prod_p(2\omega_p)^{n_p/2}}.
\]

The Kronecker delta enforces total spatial momentum conservation. Integrating over time gives the first-order Dyson amplitude

\[
\boxed{
\mathcal A_{\{n_p\}}^{(1)}
=
-\frac{i}{L}\,
\frac{\delta_{\sum_p n_p p,\,0}\,
\widetilde\lambda\!\left(\sum_p n_p\omega_p\right)}
{\sqrt{\prod_p n_p!}\,\prod_p(2\omega_p)^{n_p/2}},
\qquad \sum_p n_p=4.
}
\]

For four distinct momenta \(p_1,p_2,p_3,p_4\), this reduces to

\[
\mathcal A_{0\to p_1p_2p_3p_4}^{(1)}
=
-\frac{i}{L}\,
\frac{\delta_{p_1+p_2+p_3+p_4,\,0}\,
\widetilde\lambda(\omega_1+\omega_2+\omega_3+\omega_4)}
{\sqrt{(2\omega_1)(2\omega_2)(2\omega_3)(2\omega_4)}}.
\]

The formula for occupation states with repeated momenta includes the Bose factors explicitly through \(\prod_p n_p!\); this avoids applying the distinct-mode expression to degenerate configurations.

The time-dependent pulse supplies energy, so the total final energy appears as the argument of \(\widetilde\lambda\). Energy is not constrained by a delta function because the Hamiltonian is explicitly time dependent.

## 5. Momentum-resolved inclusive two-particle spectrum

For \(k\ne q\), the leading contribution is the sum over all normalized four-particle final states containing at least one particle in each measured mode:

\[
N_2^{(2)}(k,q)=
\sum_{\{n_p\}:\,\sum n_p=4}
n_k n_q\,\left|\mathcal A_{\{n_p\}}^{(1)}\right|^2.
\]

Substituting the amplitude gives the finite-volume expression

\[
\boxed{
N_2^{(2)}(k,q)=
\frac{1}{L^2}
\sum_{\substack{\{n_p\}:\,\sum n_p=4\\n_k\ge1,\ n_q\ge1}}
\frac{n_kn_q\,
\delta_{\sum_p n_p p,\,0}\,
\left|\widetilde\lambda\!\left(\sum_p n_p\omega_p\right)\right|^2}
{\left(\prod_p n_p!\right)\left(\prod_p(2\omega_p)^{n_p}\right)}.
}
\]

The sum is over distinct occupation-number configurations, not ordered momentum tuples. This is the expression to implement for a finite-volume perturbative benchmark. It automatically includes configurations in which an unobserved particle shares momentum with \(k\) or \(q\), provided the occupation-number factorials are retained.

For four distinct final momenta \(\{k,q,r,s\}\), the corresponding contribution is simply

\[
\frac{1}{L^2}
\frac{\delta_{k+q+r+s,\,0}\,
|\widetilde\lambda(\omega_k+\omega_q+\omega_r+\omega_s)|^2}
{(2\omega_k)(2\omega_q)(2\omega_r)(2\omega_s)},
\]

with each unordered final Fock configuration counted once.

The analogous leading single-particle spectrum is

\[
N_1^{(2)}(k)=
\frac{1}{L^2}
\sum_{\substack{\{n_p\}:\,\sum n_p=4\\n_k\ge1}}
\frac{n_k\,
\delta_{\sum_p n_p p,\,0}\,
|\widetilde\lambda(\sum_p n_p\omega_p)|^2}
{\left(\prod_p n_p!\right)\left(\prod_p(2\omega_p)^{n_p}\right)}.
\]

## 6. Why the leading contribution is non-Gaussian

Normal ordering removes the one-vertex vacuum-to-two-particle matrix element. Thus, for the free vacuum and this pulse,

\[
\langle a_k^\dagger a_k\rangle=O(\lambda_0^2),
\qquad
\langle a_k a_q\rangle=O(\lambda_0^2)
\]

or smaller. Translation invariance gives \(\langle a_k^\dagger a_q\rangle=0\) for \(k\ne q\), and \(\langle a_k a_q\rangle\) can be nonzero only in the opposite-momentum channel \(k+q=0\).

Consequently, for distinct modes, every term in the Gaussian reconstruction is \(O(\lambda_0^4)\) or smaller:

\[
N_2^{\rm G}(k,q)=O(\lambda_0^4).
\]

By contrast, the four-particle production probability gives

\[
N_2(k,q)=N_2^{(2)}(k,q)+O(\lambda_0^3)
\]

in formal coupling order, with \(N_2^{(2)}=O(\lambda_0^2)\). Therefore,

\[
\boxed{
\kappa_4(k,q)=N_2^{(2)}(k,q)+O(\lambda_0^3)
}
\]

for distinct modes, subject to the stated normal-ordering and initial-vacuum assumptions. The explicit coefficient at order \(\lambda_0^2\) is the positive occupation-state sum in Section 5. Higher orders must be calculated before assigning the next nonzero power; the displayed remainder is a conservative order symbol, not a claim that an \(O(\lambda_0^3)\) term is nonzero.

This order counting is the key result of the present benchmark: the leading pair-inclusive signal comes from a correlated four-particle production event, whereas reconstruction from final two-point functions starts at higher order.

## 7. Relation to the connected field correlator

The number cumulant is obtained by applying the exact out-mode projectors to the equal-time, operator-ordered connected four-field correlator, including the required time derivatives because \(a_k^{out}\) is linear in both \(\phi_k\) and \(\pi_k\). The ordering must remain that of

\[
a_k^\dagger a_q^\dagger a_q a_k.
\]

A different projection, with four annihilation operators, measures the anomalous coherence \(\langle a_{p_1}a_{p_2}a_{p_3}a_{p_4}\rangle_c\), which is already \(O(\lambda)\) and is proportional to the amplitude in Section 4. Its modulus squared contributes to probabilities at \(O(\lambda^2)\), but it is not itself the number-type cumulant.

In Keldysh language, the \(O(\lambda)\) four-annihilator coherence is represented by a single-branch quartic vertex. The \(O(\lambda^2)\) number-type cumulant is represented by the two-vertex cut diagram with one vertex on each contour branch. This distinction prevents confusing a tree-level connected four-field Green function with the inclusive out-state number cumulant.

## 8. Checks required before numerical implementation

1. **Finite-volume combinatorics:** verify the amplitude for four distinct modes and for repeated-mode configurations by direct operator algebra.
2. **Momentum conservation:** every contributing configuration must satisfy \(\sum_p n_p p=0\).
3. **Pulse dependence:** probabilities depend on \(|\widetilde\lambda(\Omega)|^2\), and are invariant under a common time translation of the pulse.
4. **Free limit:** all production vanishes as \(\lambda_0\to0\).
5. **Gaussian subtraction:** compute all two-point contractions independently and verify that they do not contaminate the \(O(\lambda_0^2)\) coefficient for distinct modes.
6. **Opposite-momentum channel:** test \(q=-k\) separately, retaining the anomalous two-point contraction.
7. **Coincident modes:** derive \(N_2(k,k)=\langle n_k(n_k-1)\rangle\) separately; do not use the distinct-mode formula.
8. **Keldysh/Dyson agreement:** evaluate the cross-branch expression and compare it with the normalized Fock-state sum.
9. **Regulator:** keep finite \(L\) and a specified momentum cutoff for the first benchmark; test cutoff convergence for the smooth pulse before making a continuum claim.

## 9. Current status

Established in this derivation:
- the first nonzero order of the number-type two-particle observable is \(O(\lambda_0^2)\);
- the surviving second-order in-in term is the cross-branch Keldysh contribution;
- the normalized finite-volume \(0\to4\) amplitude, including repeated-mode factors;
- the finite-volume occupation-state sum for the leading inclusive \(N_2(k,q)\);
- the leading Gaussian/non-Gaussian order separation for distinct modes.

Not yet completed:
- an independent symbolic/combinatorial implementation of the repeated-mode formula;
- explicit evaluation for the chosen smooth pulse;
- a complete next-order calculation;
- a numerical comparison with exact finite-volume Hamiltonian evolution;
- continuum and cutoff convergence.

No numerical result or final novelty claim is made here.
