# Research Repository Instructions

The shared minimum requirements for research governance, compute execution, task handoff, scientific-change tracking, quality gates, result archival, structured records, document review, and Git integration are maintained by the central research infrastructure repository. The adopted revision is pinned in [research/standards-lock.yaml](research/standards-lock.yaml). Start from the [Central Research Standards Index](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/governance/CENTRAL_STANDARDS_INDEX.md).

The central standards define common research-process requirements. This repository remains authoritative for the physical model, field-theory conventions, finite-volume normalization, CTP/operator ordering, out-mode observables, numerical implementation, and project-specific scientific acceptance criteria.

## 1. Preflight and task handoff

Before substantial work:

1. Read this file, the project README, STATUS.md, and research/standards-lock.yaml.
2. Read the pinned Central Research Standards Index.
3. For a new derivation, numerical run, parameter scan, implementation change, or multi-agent handoff, prepare a concise Research Task Brief using the [central Task Handoff Standard](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/RESEARCH_TASK_HANDOFF_STANDARD.md) and [template](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/templates/research-task-brief.md). Use the existing research/task-record layout; do not restructure the project solely to host a brief. Historical work does not need retroactive briefs.
4. Select the task-appropriate [central Research Quality Gate](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/governance/RESEARCH_QUALITY_GATE.md) profile and define relevant acceptance checks before interpreting new results.
5. Inspect existing source, run manifests, raw outputs, parent runs, and provenance before starting. Reuse verified prior results rather than repeating calculations without scientific need.

## 2. Project-specific physics conventions

The following local scientific requirements remain authoritative:

- Work in natural units hbar = c = 1 unless stated otherwise.
- State the metric, Fourier-transform convention, finite-volume normalization, momentum labels, and normal-ordering/renormalization prescription.
- Distinguish equal-time, time-ordered, contour-ordered, Wightman, retarded, and statistical correlators.
- Distinguish a field-theory connected four-point function from the connected cumulant of out-mode creation/annihilation operators. State the projection and operator ordering connecting them.
- For k != q, state whether q = -k is included; anomalous contractions can contribute in that channel.
- Treat k = q separately because of bosonic commutators and factorial moments.
- A spatially homogeneous time-dependent coupling conserves total momentum but can inject energy during the pulse.
- State the regulator, finite-volume/Fock-space truncation, counterterms, vacuum-energy convention, and continuum prescription before making claims about regulator or continuum limits.
- Do not replace the connected four-point sector by a Gaussian closure without explicitly defining and validating that approximation.

## 3. Scientific integrity and validation

For every substantive calculation, preserve:

1. Exact source commit SHA and clean/dirty working-tree state.
2. SHA-256 hashes of all inputs, source files, and configuration files.
3. Actual runtime environment, interpreter/compiler and dependency versions.
4. Hardware and execution parameters.
5. Native, uncompressed numerical outputs.
6. Numerical convergence, truncation dependence, and uncertainty analysis.
7. Independent physical validation and analytic limiting cases.
8. A run manifest linking the source, inputs, environment, outputs, and validation evidence.

Do not label a calculation complete unless its native output exists and is tied to the actual source revision and run manifest. Do not infer physical conclusions from a passing software test alone. Distinguish computed, archived, numerically validated, and physically validated status.

Preserve exact equations, signs, normalization factors, operator ordering, and historical results. When correcting a derivation or implementation, identify the change and its effect; do not silently replace earlier reasoning or overwrite previous runs.

## 4. Repository and Git workflow

- Work on a named research branch; do not develop directly on main.
- Open a focused PR for substantive research or governance changes.
- Follow the pinned [central Git merge policy](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/governance/GIT_MERGE_POLICY.md). Substantive research defaults to Create a merge commit.
- Preserve raw inputs, native outputs, and historical results. New runs create new immutable result directories.
- Record exact base/head SHAs and relevant acceptance checks in the PR.
- After merge, verify the integrated main state and remote files separately from scientific validation.

## 5. Reporting

Follow the central [Research Task Handoff Standard](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/RESEARCH_TASK_HANDOFF_STANDARD.md) and [Scientific Diff reporting standard](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/scientific-diff-workflow.md#11-research-reporting-and-communication-standard).

For this project, report the exact regulator, volume, Fock/momentum cutoffs, solver settings, source/input fingerprints, raw outputs, convergence data, and independent checks. Keep finite-volume/Fock-space results distinct from continuum or infinite-volume claims.
