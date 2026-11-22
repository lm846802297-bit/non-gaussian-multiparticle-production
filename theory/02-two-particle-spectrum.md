# Outgoing Two-Particle Spectrum and Connected Cumulant

## 1. Inclusive spectra

For out-mode operators defined in `theory/01-model-and-asymptotic-modes.md`,

\[
N_1(k)=\langle a_k^\dagger a_k\rangle_{out},
\]

and for distinct modes `k != q`,

\[
N_2(k,q)=\langle a_k^\dagger a_q^\dagger a_q a_k\rangle_{out}.
\]

This is an inclusive factorial two-particle spectrum: it counts ordered mode occupations in the state, not the exclusive probability of producing exactly two particles. The coincident-mode observable is `N_2(k,k)=<a_k^dagger a_k^dagger a_k a_k>` and must be handled separately.

## 2. Gaussian reconstruction

For a zero-mean Gaussian state, Wick's theorem gives, for distinct modes,

\[
\begin{aligned}
N_2^G(k,q)={}&\langle a_k^\dagger a_k\rangle\langle a_q^\dagger a_q\rangle\\
&+\langle a_k^\dagger a_q\rangle\langle a_q^\dagger a_k\rangle\\
&+\langle a_k^\dagger a_q^\dagger\rangle\langle a_q a_k\rangle.
\end{aligned}
\]

Equivalently,

\[
N_2^G(k,q)=N_1(k)N_1(q)+|\langle a_k^\dagger a_q\rangle|^2+|\langle a_k a_q\rangle|^2.
\]

The anomalous term is important in channels such as `q=-k`; it must not be dropped merely because `k != q`.

## 3. Connected four-operator cumulant

Define

\[
\boxed{\kappa_4(k,q)=N_2(k,q)-N_2^G(k,q).}
\]

For a zero-mean state this is the connected cumulant associated with the operator ordering `a_k^dagger, a_q^dagger, a_q, a_k`. It is not automatically identical to an arbitrary time-ordered field-theory four-point function. The exact relation is obtained by substituting the linear out-mode projection for each operator and preserving the specified Wightman/operator ordering.

Thus the exact target is

\[
\boxed{N_2(k,q)=N_2^G(k,q)+\kappa_4(k,q).}
\]

A nonzero `kappa4` proves that the two-particle inclusive spectrum is not reconstructible from the final two-point data alone.

## 4. Single-particle spectrum from equal-time data

The out-mode projection also yields

\[
\begin{aligned}
N_1(k)=\frac{1}{2\omega_k}\big[&\langle\pi_k^\dagger\pi_k\rangle+\omega_k^2\langle\phi_k^\dagger\phi_k\rangle\\
&+i\omega_k\langle\phi_k^\dagger\pi_k-\pi_k^\dagger\phi_k\rangle\big]_{out}.
\end{aligned}
\]

This identity is exact in the free out region, with the Fourier convention stated in the model document.

## 5. Research diagnostic

A useful relative diagnostic is

\[
\mathcal R_2(k,q)=\frac{\kappa_4(k,q)}{N_2^G(k,q)},
\]

where the denominator is nonzero. Report `kappa4` and `N2_G` separately as well; a ratio alone can become large when the Gaussian denominator is small.

For translation-invariant states, the normal two-point contraction is diagonal in momentum and the anomalous contraction pairs opposite momenta. Report generic `k,q` and the `q=-k` channel separately.