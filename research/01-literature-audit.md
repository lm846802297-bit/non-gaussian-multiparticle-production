# Literature Audit and Provisional Novelty Assessment

Audit date: 2026-11-22. This is a targeted first-pass audit, not an exhaustive systematic review. The records below were checked against DOI/arXiv metadata and, where available, full text.

## Directly relevant prior work

### 1. Non-Gaussian initial correlations in Kadanoff-Baym dynamics

M. Garny and M. M. Müller, *Kadanoff-Baym Equations with Non-Gaussian Initial Conditions: The Equilibrium Limit*, Phys. Rev. D 80, 085011 (2009), DOI: [10.1103/PhysRevD.80.085011](https://doi.org/10.1103/PhysRevD.80.085011), arXiv: [0904.3600](https://arxiv.org/abs/0904.3600).

They derive Kadanoff-Baym equations for non-Gaussian correlated initial states from the 2PI effective action and numerically study a real scalar `Phi^4` theory with the thermal initial four-point correlation as the leading non-Gaussian correction. This rules out novelty claims based merely on adding non-Gaussian initial four-point data to Kadanoff-Baym evolution. Their focus is initial-state correlations and the equilibrium limit, not the proposed out-state momentum-resolved `kappa4(k,q)` after a finite interaction pulse.

### 2. Far-from-equilibrium two-point and four-vertex dynamics

L. Shen and J. Berges, *Spectral, Statistical, and Vertex Functions in Scalar Quantum Field Theory Far from Equilibrium*, Phys. Rev. D 101, 056009 (2020), DOI: [10.1103/PhysRevD.101.056009](https://doi.org/10.1103/PhysRevD.101.056009), arXiv: [1912.07565](https://arxiv.org/abs/1912.07565).

They compute far-from-equilibrium dynamics in an O(N) scalar theory using a next-to-leading-order 1/N expansion of the 2PI effective action, including two-point statistical/spectral functions and the momentum-dependent four-vertex. This rules out novelty claims based simply on real-time evolution of a four-vertex or far-from-equilibrium scalar-field correlators. Their reported target is scaling of two-point functions and the four-vertex, not a momentum-resolved outgoing connected two-particle spectrum for the finite-pulse setup.

### 3. Field-theory mapping to particle spectra and interferometry

F. Rennecke, R. D. Pisarski, and D. H. Rischke, *Particle Interferometry in a Moat Regime*, Phys. Rev. D 107, 116011 (2023), DOI: [10.1103/PhysRevD.107.116011](https://doi.org/10.1103/PhysRevD.107.116011), arXiv: [2301.11484](https://arxiv.org/abs/2301.11484).

They develop a field-theoretic relation between particle spectra and in-medium real-time correlation functions on curved spacetime hypersurfaces, and study two-particle Bose-Einstein correlations in a moat regime. Gaussian/Wick reconstruction is central to the tractable two-particle result. This rules out novelty claims based on the general `2n`-point-to-`n`-particle-spectrum mapping or Gaussian HBT formulas. The proposed project must focus on a dynamically generated connected non-Gaussian out-state cumulant beyond that reconstruction.

### 4. Kinetic theory beyond the Boltzmann equation

X.-Y. Hu and V. Rosenhaus, *Beyond the Boltzmann Equation for Weakly Coupled Quantum Fields*, JHEP 10 (2025) 070, DOI: [10.1007/JHEP10(2025)070](https://doi.org/10.1007/JHEP10(2025)070), arXiv: [2503.09932](https://arxiv.org/abs/2503.09932).

They derive a perturbative kinetic equation for occupation-number evolution in weakly coupled, nearly homogeneous, nearly stationary, nearly Gaussian states. The derivation uses equal-time four-point functions and includes higher-order/multiparticle scattering corrections; Appendix C.1 treats `lambda phi^4`. This rules out novelty claims based merely on using a four-point correlator to derive a correction to `dn_k/dt` or on weak-coupling higher-order kinetic terms. Their target is a kinetic equation for occupation numbers under stated near-stationary assumptions, not the full momentum-resolved outgoing `N2(k,q)` and its connected cumulant for a finite-duration pulse.

### 5. Non-Gaussian initial cumulants in Keldysh theory

A. G. Mikhaylenko and A. G. Semenov, *Keldysh Diagram Technique with Non-Gaussian Initial Correlations*, Theoretical and Mathematical Physics 227 (2026), 680–704, DOI: [10.1134/S0040577926040094](https://doi.org/10.1134/S0040577926040094).

They formulate Keldysh diagrammatics with non-Gaussian initial correlations through initial cumulant vertices. This rules out novelty claims based merely on introducing initial non-Gaussian cumulants into the Keldysh action. The proposed target is instead dynamical generation from a Gaussian in-state and its projection onto an out-state spectrum.

### 6. Nonperturbative real-time particle production with tensor networks

E. Budd, A. Florio, D. Frenklakh, and S. Mukherjee, *Quantum Dynamics of Cosmological Particle Production: Interacting Quantum Field Theories with Matrix Product States*, JHEP 04 (2026) 183, DOI: [10.1007/JHEP04(2026)183](https://doi.org/10.1007/JHEP04(2026)183).

They use matrix product states to study real-time dynamics and particle production in interacting scalar and gauge theories in 1+1 dimensions, including scalar `lambda phi^4`. This rules out novelty claims based merely on MPS real-time scalar-field evolution or single-particle production. The specific connected outgoing two-particle spectrum is not identified in the metadata and abstract examined in this first-pass audit; the full paper should be checked before final novelty claims.

### 7. Real-time two-particle scattering in `phi^4`

B. Al Sayegh and W. Chemissany, *Real-Time Scattering in `phi^4` Theory Using Matrix Product States*, Phys. Rev. Research 8, 023110 (2026), DOI: [10.1103/5v2j-12ps](https://doi.org/10.1103/5v2j-12ps), arXiv: [2511.15697](https://arxiv.org/abs/2511.15697).

They use uniform MPS and TDVP to study real-time two-particle scattering in 1+1D `phi^4`, extracting elastic probabilities and Wigner time delays. This rules out novelty claims based merely on real-time two-particle scattering with MPS. The proposed observable is an inclusive momentum-resolved two-particle spectrum and its connected cumulant after a finite-duration interaction pulse, not only the elastic `1+1 -> 1+1` channel.

## Provisional gap

The targeted search did not identify a paper that simultaneously performs all of the following in the same controlled setup:

1. starts from a specified Gaussian in-state;
2. evolves a finite-duration interacting `lambda phi^4` theory in real time;
3. computes the dynamically generated connected four-operator out-mode cumulant;
4. obtains the momentum-resolved inclusive `N2(k,q)` from the full out-state;
5. compares it quantitatively with the complete Gaussian/Wick reconstruction from the same final two-point data;
6. supplies finite-volume, regulator, truncation, and independent validation errors.

This is a provisional novelty gap, not proof of priority. Before submission, trace citations and references of the papers above and search for exact phrases around out-particle number correlations, factorial moments, connected four-point functions, and real-time quenches.

## Explicit no-go list

Do not claim as new: Wigner transforms; the general relation between `2n`-point functions and `n`-particle inclusive spectra; Gaussian Wick reconstruction of two-particle spectra; non-Gaussian initial-state Kadanoff-Baym equations; generic connected-correlation hierarchies; weak-coupling corrections to Boltzmann occupation-number evolution; generic MPS real-time evolution of `phi^4`.