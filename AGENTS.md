# Research Repository Instructions

## Scientific integrity

- Distinguish established identities, derived results, approximations, conjectures, and pending calculations.
- Never label a calculation as completed unless the native output exists and can be tied to a source revision and run manifest.
- Preserve exact equations, conventions, normalization factors, and operator ordering.
- When correcting a derivation, record the correction and its consequences; do not silently replace earlier reasoning.
- Do not infer physical conclusions from a numerical run before convergence and independent validation.

## Reproducibility requirements

For every substantive calculation, archive:

1. Exact source commit SHA and a clean/dirty working-tree status.
2. SHA-256 hashes of all input files and configuration files.
3. Runtime environment, compiler/interpreter and library versions.
4. Hardware and execution parameters.
5. Native, uncompressed numerical outputs.
6. Numerical convergence and uncertainty analysis.
7. Independent physical validation and analytic limiting cases.
8. A run manifest linking all artifacts.

Do not commit credentials, tokens, private keys, personal data, or machine-specific secrets.

## Git workflow

- Work on a named research branch; do not develop directly on `main`.
- Open a pull request for substantive research changes.
- Do not merge without review and successful, relevant validation.
- Preserve raw data and historical results. New runs create new immutable result directories.
- Record the exact base and head SHAs in the PR description.

## Physics-specific rules

- Use natural units `hbar = c = 1` unless stated otherwise.
- State Fourier-transform conventions and finite-volume normalization explicitly.
- Distinguish equal-time correlators from contour-ordered, Wightman, retarded, and statistical correlators.
- Distinguish a field-theory connected four-point function from a connected cumulant of out-mode creation/annihilation operators; state the projection and ordering connecting them.
- For `k != q`, specify whether `q = -k` is included; anomalous contractions can contribute in that channel.
- Treat the coincident-mode case `k = q` separately because of bosonic commutators and factorial moments.
- State the renormalization/normal-ordering prescription and regulator before comparing results.
- A time-dependent coupling injects energy. Momentum is conserved for a spatially homogeneous pulse, while energy need not be conserved during the pulse.