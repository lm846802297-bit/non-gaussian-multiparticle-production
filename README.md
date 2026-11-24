# Non-Gaussian Multiparticle Production in Nonequilibrium Quantum Field Theory

**Research governance:** This repository adopts the shared standards pinned in [research/standards-lock.yaml](research/standards-lock.yaml). See the [Central Research Standards Index](https://github.com/lm846802297-bit/high-energy-theory-research/blob/14155487afc075d301fa7f6cf452c7bff137e869/docs/governance/CENTRAL_STANDARDS_INDEX.md). The project-specific QFT model, CTP conventions, out-mode observables, numerical methods, and physical validation remain local.

Research archive for the study of dynamically generated connected four-point correlations in a time-dependent interacting quantum field theory and their contribution to outgoing two-particle spectra.

**Repository name:** `non-gaussian-multiparticle-production`  
**Visibility:** public  
**Status:** research bootstrap / theory development; no numerical result is claimed.

## Research question

Starting from a specified initial density operator and an interacting quantum field theory, how are non-Gaussian connected correlations generated in real time, and how do they contribute to measurable inclusive multiparticle spectra?

The first controlled model is a real scalar field in 1+1 dimensions with a finite-duration, spatially homogeneous quartic interaction pulse,

\[
\mathcal L=\frac12(\partial\phi)^2-\frac12m^2\phi^2-
\frac{\lambda(t)}{4!}:\!\phi^4\!:\,.
\]

The initial state is the free Gaussian vacuum. The coupling is switched on and off smoothly so that free in/out particle operators are unambiguously defined.

## Primary observables

\[
N_1(k)=\langle a_k^\dagger a_k\rangle_{\rm out},
\qquad
N_2(k,q)=\langle a_k^\dagger a_q^\dagger a_q a_k\rangle_{\rm out}.
\]

For distinct modes, decompose

\[
N_2(k,q)=N_2^{\rm G}(k,q)+\kappa_4(k,q),
\]

where `N2_G` is the Wick reconstruction from all two-point contractions and `kappa4` is the connected four-operator cumulant. The central target is a controlled calculation of `kappa4` and its effect on the outgoing two-particle spectrum.

## Start here

1. [Project charter and scope](research/00-project-charter.md)
2. [Literature audit and novelty assessment](research/01-literature-audit.md)
3. [Model and asymptotic particle operators](theory/01-model-and-asymptotic-modes.md)
4. [Two-particle spectrum and cumulant decomposition](theory/02-two-particle-spectrum.md)
5. [Perturbative benchmark plan](theory/03-perturbative-benchmark.md)
6. [Research status](STATUS.md)
7. [Provenance and archival policy](provenance/README.md)
8. [Bibliography](references/references.bib)

## Research workflow

Scientific question -> literature audit -> derivation -> analytic benchmark -> numerical implementation -> convergence and error analysis -> independent physical validation -> claim/evidence review -> pull request -> review -> merge.

No result is considered complete merely because code runs. Every numerical result must be tied to an exact source revision, input hashes, runtime environment, hardware, parameters, raw outputs, convergence evidence, and an independent physics check.

## Repository structure

- `research/`: scope, literature audit, novelty assessment, decisions, and research log.
- `theory/`: derivations and mathematical conventions.
- `references/`: BibTeX records and literature notes.
- `src/`: numerical and symbolic source code.
- `tests/`: analytic limits, regression tests, and independent checks.
- `data/raw/`: immutable input data; never overwrite original inputs.
- `data/processed/`: derived datasets with provenance.
- `results/`: native numerical outputs and run manifests.
- `figures/`: generated figures with source scripts and metadata.
- `provenance/`: environment snapshots, hashes, and archival rules.

## Scientific status boundary

This bootstrap contains a research plan, a preliminary literature map, and derivation notes. It contains no completed numerical calculation, no validated non-Gaussian spectrum, and no claim that the proposed novelty gap is conclusively established. The literature audit must be extended before manuscript submission.