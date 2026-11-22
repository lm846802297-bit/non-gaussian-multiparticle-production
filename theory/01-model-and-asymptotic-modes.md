# Model, Conventions, and Asymptotic Particle Operators

## 1. Model

Work in 1+1 dimensions and natural units `hbar = c = 1`:

\[
\mathcal L(t,x)=\frac12(\partial_t\phi)^2-\frac12(\partial_x\phi)^2-
\frac12m^2\phi^2-\frac{\lambda(t)}{4!}:\!\phi^4\!:\,_m.
\]

The normal ordering is with respect to the free mass-`m` vacuum. This is a proposed renormalization convention for the first implementation; the required vacuum-energy subtraction and the regulator-to-continuum prescription must be documented and tested explicitly.

The canonical momentum is `pi = partial_t phi`, and the Hamiltonian is

\[
H(t)=\int_0^L dx\left[\frac12\pi^2+\frac12(\partial_x\phi)^2+\frac12m^2\phi^2+\frac{\lambda(t)}{4!}:\!\phi^4\!:\,_m\right].
\]

The Heisenberg equations are

\[
\dot\phi=\pi,\qquad
\dot\pi=\partial_x^2\phi-m^2\phi-\frac{\lambda(t)}{3!}:\!\phi^3\!:\,_m.
\]

A spatially homogeneous pulse preserves total momentum. Since `H(t)` depends explicitly on time, energy is not conserved during the pulse:

\[
\frac{d}{dt}\langle H(t)\rangle=\left\langle\frac{\partial H}{\partial t}\right\rangle
=\frac{\dot\lambda(t)}{4!}\int_0^L dx\,\langle:\!\phi^4\!:\,_m\rangle.
\]

After the pulse, `lambda(t)=0` and the free Hamiltonian is time independent.

## 2. Finite-volume Fourier convention

Use a periodic box of length `L`, with momenta `k = 2 pi n/L`. Define

\[
\phi_k(t)=\frac{1}{\sqrt L}\int_0^L dx\,e^{-ikx}\phi(t,x),\qquad
\pi_k(t)=\frac{1}{\sqrt L}\int_0^L dx\,e^{-ikx}\pi(t,x).
\]

For a free out region,

\[
\phi_k(t)=\frac{1}{\sqrt{2\omega_k}}\left(a_k^{out}e^{-i\omega_kt}+a_{-k}^{out\dagger}e^{i\omega_kt}\right),
\]

\[
\pi_k(t)=-i\sqrt{\frac{\omega_k}{2}}a_k^{out}e^{-i\omega_kt}+i\sqrt{\frac{\omega_k}{2}}a_{-k}^{out\dagger}e^{i\omega_kt},
\qquad \omega_k=\sqrt{k^2+m^2}.
\]

Combining these equations gives the exact out-mode projection

\[
\boxed{a_k^{out}=e^{i\omega_kt}\left[\sqrt{\frac{\omega_k}{2}}\phi_k(t)+\frac{i}{\sqrt{2\omega_k}}\pi_k(t)\right]}
\]

for any time in the free out region. The formula fixes the finite-volume normalization used throughout the project.

## 3. Pulse protocol

Use a smooth finite-duration pulse, for example

\[
\lambda(t)=\frac{\lambda_0}{4}
\left[1+\tanh\left(\frac{t-t_i}{\tau}\right)\right]
\left[1-\tanh\left(\frac{t-t_f}{\tau}\right)\right].
\]

The pulse is spatially homogeneous. It injects energy but does not inject net spatial momentum. The Fourier transform convention is

\[
\widetilde\lambda(\Omega)=\int_{-\infty}^{\infty}dt\,\lambda(t)e^{i\Omega t}.
\]

The pulse parameters and the treatment of the small nonzero tails at finite simulation start/end times must be recorded in each run manifest.

## 4. Initial state and symmetry

The baseline initial state is the free Gaussian vacuum of mass `m`, with zero one-point function. The theory and pulse preserve `Z_2: phi -> -phi`; odd field correlators vanish if the initial state respects this symmetry. Translation invariance enforces total-momentum selection rules.

## 5. Regulator and renormalization are mandatory inputs

The numerical implementation will use a spatial lattice with spacing `a`, `N_x` sites, and `L=N_x a`. The momentum cutoff is of order `pi/a`. A continuum claim requires a sequence of `a -> 0` calculations at fixed physical volume and renormalized parameters. Normal ordering relative to the free mass is a chosen prescription, not a substitute for documenting the complete counterterm and vacuum-energy convention.