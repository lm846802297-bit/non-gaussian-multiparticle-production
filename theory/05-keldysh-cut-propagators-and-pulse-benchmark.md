# Four-Cut Keldysh Verification and Explicit Pulse Benchmark

Status: leading-order perturbative benchmark completed for the stated finite-volume parameters. This is not a nonperturbative real-time simulation.

## 1. Model and conventions

We use the normal-ordered interaction

\[
H_I(t)=\frac{\lambda(t)}{4!}\int_0^L dx\,:\!\phi_I^4(t,x)\!:\,_m,
\]

in 1+1 dimensions, with periodic length \(L\), discrete momenta \(p=2\pi n/L\), \(\omega_p=\sqrt{p^2+m^2}\), and

\[
\phi_I(t,x)=\frac{1}{\sqrt L}\sum_p\frac{a_p e^{-i\omega_pt+ipx}+a_p^\dagger e^{i\omega_pt-ipx}}{\sqrt{2\omega_p}},
\qquad [a_p,a_q^\dagger]=\delta_{pq}.
\]

The initial state is the free Gaussian vacuum. The interaction is spatially homogeneous and tends to zero at early and late times. The measured observable is

\[
N_2(k,q)=\langle a_k^\dagger a_q^\dagger a_q a_k\rangle_{\rm out}.
\]

For \(k\ne q\), its occupation weight in a final Fock state is \(n_kn_q\). For \(k=q\), it is \(n_k(n_k-1)\), not \(n_k^2\).

## 2. Independent closed-time-path cut calculation

At the first nonzero order, there is one quartic vertex on the forward branch and one on the backward branch. Their vertex factors are \((-i\lambda(t))(+i\lambda(t'))/(4!)^2\), so their product has positive sign. All four lines connecting the two vertices are cross-branch Wightman (cut) propagators.

Define the positive-frequency cut line, with no extra \(i\) convention, by

\[
\Delta_p^{\rm cut}(x,y)=\frac{1}{2\omega_p L}
 e^{-i\omega_p(t_x-t_y)+ip(x-y)},
\]

and the free Wightman function

\[
D^>(x,y)=\langle0|\phi_I(x)\phi_I(y)|0\rangle
=\sum_p\Delta_p^{\rm cut}(x,y).
\]

Two of the four cut lines are fixed by the measured out-mode insertions \(k,q\); the other two are unobserved and are summed by \([D^>(x,y)]^2\).

The contraction count is explicit. At each vertex, the two ordered measured insertions can be attached to two distinct legs in \(4\times3\) ways. The two remaining legs are joined across the cut in \(2!\) ways. Dividing by the two vertex factors \((4!)^2\), the net factor is

\[
\frac{(4\cdot3)^2\,2!}{(4!)^2}=\frac12.
\]

Therefore the independent cut-diagram expression is

\[
\boxed{
N_2^{(2)}(k,q)=\frac12\int dt_xdt_y\,d x d y\,
\lambda(t_x)\lambda(t_y)
\Delta_k^{\rm cut}(x,y)\Delta_q^{\rm cut}(x,y)
[D^>(x,y)]^2.
}
\]

This expression applies to both distinct and coincident measured modes, provided the coincident observable is the factorial moment \(\langle n_k(n_k-1)\rangle\). The two ordered identical-mode insertions are precisely why the same \(4\times3\) assignment count applies.

Expanding the two unobserved cut lines gives

\[
[D^>(x,y)]^2=\sum_{r,s}\Delta_r^{\rm cut}(x,y)\Delta_s^{\rm cut}(x,y).
\]

The spatial integrals produce

\[
\int_0^L dx\int_0^Ldy\,e^{i(k+q+r+s)(x-y)}
=L^2\delta_{k+q+r+s,0}.
\]

The time integrals produce the squared pulse transform. Thus

\[
\boxed{
N_2^{(2)}(k,q)=\frac{1}{2L^2}\sum_{r,s}
\frac{\delta_{k+q+r+s,0}\,
|\widetilde\lambda(\omega_k+\omega_q+\omega_r+\omega_s)|^2}
{(2\omega_k)(2\omega_q)(2\omega_r)(2\omega_s)}.
}
\]

The factor \(1/2\) is not optional: it is the residual symmetry factor after the measured legs have been marked and the two unobserved cut lines are summed as ordered momenta.

## 3. Independent Fock-state sum

The normalized occupation state with four particles is

\[
|\{n_p\}\rangle=\prod_p\frac{(a_p^\dagger)^{n_p}}{\sqrt{n_p!}}|0\rangle,
\qquad \sum_p n_p=4.
\]

The first-order vacuum-to-four-particle amplitude is

\[
\mathcal A_{\{n_p\}}^{(1)}=-\frac{i}{L}
\frac{\delta_{\sum_p n_pp,0}\,
\widetilde\lambda(\sum_p n_p\omega_p)}
{\sqrt{\prod_p n_p!}\,\prod_p(2\omega_p)^{n_p/2}}.
\]

Consequently, the independent Fock-state expression is

\[
N_2^{(2)}(k,q)=\frac{1}{L^2}
\sum_{\substack{\{n_p\}:\sum_p n_p=4}}
\frac{w_{kq}(\{n_p\})\,\delta_{\sum_p n_pp,0}
|\widetilde\lambda(\sum_p n_p\omega_p)|^2}
{\left(\prod_p n_p!\right)\left(\prod_p(2\omega_p)^{n_p}\right)},
\]

where

\[
w_{kq}(\{n_p\})=
\begin{cases}
n_kn_q,&k\ne q,\\
n_k(n_k-1),&k=q.
\end{cases}
\]

The Fock-state sum enumerates unordered occupation configurations. The cut sum enumerates the two unobserved momenta as an ordered pair and carries the factor \(1/2\). The two forms are independent implementations of the same leading-order result.

## 4. Explicit smooth pulse

Choose

\[
\boxed{\lambda(t)=\lambda_0\,\operatorname{sech}^2(t/T)}.
\]

This pulse is smooth and tends exponentially to zero as \(t\to\pm\infty\), so free in/out modes are well defined. With

\[
\widetilde\lambda(\Omega)=\int_{-\infty}^{\infty}dt\,\lambda(t)e^{i\Omega t},
\]

the Fourier transform is

\[
\boxed{
\widetilde\lambda(\Omega)=
\lambda_0\frac{\pi T^2\Omega}{\sinh(\pi T\Omega/2)},
\qquad \widetilde\lambda(0)=2\lambda_0T.
}
\]

The transform is real and even for this pulse. The pair of Keldysh time integrals is

\[
\int dt_xdt_y\,\lambda(t_x)\lambda(t_y)e^{-i\Omega(t_x-t_y)}
=|\widetilde\lambda(\Omega)|^2.
\]

## 5. Numerical benchmark parameters

The finite-volume momentum sum was evaluated with

- \(m=1\);
- \(L=20\);
- \(T=1\);
- \(\lambda_0/m^2=0.1\);
- \(\Delta k=2\pi/L=0.3141592653589793\);
- symmetric momentum-index cutoff \(|n|\le N_{\rm cut}\).

Three measured channels were checked:

| Channel | Mode indices | \((k/m,q/m)\) |
|---|---:|---:|
| Generic | \((1,2)\) | \((0.3141593,0.6283185)\) |
| Opposite momenta | \((1,-1)\) | \((0.3141593,-0.3141593)\) |
| Coincident | \((1,1)\) | \((0.3141593,0.3141593)\) |

At \(N_{\rm cut}=24\), the results are:

| Channel | \(N_2^{(2)}\) |
|---|---:|
| Generic | \(1.2675866570249926\times10^{-9}\) |
| Opposite momenta | \(3.7983391362711916\times10^{-9}\) |
| Coincident | \(2.8617293918174372\times10^{-9}\) |

For all three channels, the independent cut sum and Fock-state sum agree to relative differences below \(5\times10^{-16}\) at the displayed cutoffs. The coincident-mode calculation uses \(n_k(n_k-1)\); replacing it with \(n_k^2\) gives an incorrect result.

The leading single-particle result for mode index \(n_k=1\) is

\[
N_1^{(2)}(k)=5.397479856651738\times10^{-9}.
\]

## 6. Momentum-cutoff convergence

The following values show convergence of the finite-volume momentum sum:

| Channel | \(N_{\rm cut}=16\) | \(20\) | \(24\) | \(32\) |
|---|---:|---:|---:|---:|
| Generic | \(1.2675866570163574\times10^{-9}\) | \(1.2675866570249889\times10^{-9}\) | \(1.2675866570249926\times10^{-9}\) | same as 24 |
| Opposite | \(3.7983391362704273\times10^{-9}\) | \(3.7983391362711916\times10^{-9}\) | same as 20 | same |
| Coincident | \(2.8617293918120129\times10^{-9}\) | \(2.8617293918174351\times10^{-9}\) | \(2.8617293918174372\times10^{-9}\) | same as 24 |

The relative change from cutoff 20 to 24 is approximately \(2.94\times10^{-15}\) (generic), zero at double precision (opposite), and \(7.23\times10^{-16}\) (coincident). This demonstrates convergence of the unobserved-pair momentum sum for this smooth pulse and these fixed external modes. It is not a continuum or infinite-volume extrapolation.

## 7. Independent Fourier-transform check

The analytic transform was compared with composite Simpson integration over \([-20T,20T]\), using 200,000 subintervals:

| \(\Omega/m\) | Analytic \(\widetilde\lambda/m\) | Relative error |
|---:|---:|---:|
| 2 | \(5.4405810996426635\times10^{-2}\) | \(2.08\times10^{-13}\) |
| 4 | \(4.6934118610807577\times10^{-3}\) | \(2.44\times10^{-14}\) |
| 6 | \(3.0423001583780265\times10^{-4}\) | \(3.06\times10^{-14}\) |
| 8 | \(1.7529294603123114\times10^{-5}\) | \(5.42\times10^{-12}\) |

The numerical integral reproduces the analytic Fourier transform at the stated accuracy.

## 8. Physical interpretation and scope

The homogeneous pulse enforces \(k+q+r+s=0\) but does not enforce energy conservation. Instead, the pulse supplies the total final energy through \(|\widetilde\lambda(\omega_k+\omega_q+\omega_r+\omega_s)|^2\). The opposite-momentum channel is enhanced for the chosen parameters because its measured pair has zero total momentum, leaving a larger set of unobserved pairs satisfying the momentum constraint. This is a property of this finite-volume benchmark, not a universal statement about all pulse parameters.

The same formula applies to coincident modes only when the measured observable is the factorial moment \(\langle n_k(n_k-1)\rangle\). This is the required Bose-statistics correction.

The calculation is leading-order perturbation theory, \(O(\lambda_0^2)\). It does not include higher-order rescattering, loop corrections to the state, a nonperturbative time evolution, continuum extrapolation, or a comparison with an exact Hamiltonian solver.

## 9. Reproducibility

- Source: `src/keldysh_cut_benchmark.py`.
- Native output: `results/keldysh-cut-001/result.json`.
- Run metadata and SHA-256 hashes: `results/keldysh-cut-001/run-manifest.json`.
- The script independently computes the cut-Wightman sum and the unordered normalized-Fock-state sum, asserts their agreement, and checks the Fourier transform against numerical Simpson integration.
