# Project Charter

## Working title

Non-Gaussian Multiparticle Production in Nonequilibrium Quantum Field Theory

## Central question

For a specified initial quantum state and interacting field theory, can the dynamically generated connected four-point correlation be evolved and projected onto the outgoing two-particle inclusive spectrum, and when does it invalidate a reconstruction based only on two-point functions?

## First model

A real scalar field in 1+1 dimensions with a spatially homogeneous, finite-duration quartic interaction pulse:

\[
\mathcal L(t,x)=\frac12(\partial_\mu\phi)(\partial^\mu\phi)-\frac12m^2\phi^2-
\frac{\lambda(t)}{4!}:\!\phi^4\!:\,_m.
\]

The normal-ordering reference mass, counterterms, regulator, and continuum prescription must be fixed before quantitative calculations. The initial state is the free vacuum associated with mass `m`. The interaction is switched off at late times, defining free out modes.

## Main deliverables

1. Exact finite-volume projection from late-time field and momentum modes to `a_k^out`.
2. Exact decomposition of the two-particle inclusive spectrum into Gaussian Wick contractions plus a connected four-operator cumulant.
3. Closed-time-path or equivalent Hamiltonian evolution equations sufficient to calculate the relevant connected four-point function, with a stated truncation if the hierarchy is not exact.
4. Lowest nonzero perturbative benchmark for a finite interaction pulse.
5. Full quantum finite-volume computation of `N1(k)`, `N2(k,q)`, and `kappa4(k,q)`.
6. Convergence, uncertainty, and independent physical validation.

## Scope exclusions for the first paper

- No full QCD implementation in Phase I.
- No claim to solve the complete infinite hierarchy.
- No replacement of the full connected four-point function by a Gaussian closure.
- No use of a finite-time quasiparticle number as a substitute for the asymptotic out-particle number when the interaction is still active.

## Scientific hypothesis to test

A spatially homogeneous finite-duration quartic interaction can produce correlated multiparticle out states. For an initial free vacuum, the first-order amplitude for vacuum-to-four-particle production is nonzero when the pulse has Fourier support at the total final energy. The resulting two-particle inclusive spectrum can therefore receive a connected contribution that is not reconstructible from the final two-point function alone.

This is a testable hypothesis, not yet a completed result. Its precise perturbative order and normalization must be verified in the chosen finite-volume convention.

## Decision criteria

The project advances beyond a reproduction only if it produces at least one of:

- a new controlled calculation of the momentum-resolved out-state connected two-particle spectrum;
- a quantitative domain-of-validity criterion for Gaussian reconstruction;
- a benchmark showing a discrepancy between exact finite-volume evolution and a specified Gaussian/2PI/kinetic closure, with controlled errors.