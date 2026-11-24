# Provenance and Archival Policy

**Status:** Project-specific extension to the central Calculation-result Archival Policy.  
**Central source commit:** 14155487afc075d301fa7f6cf452c7bff137e869

The shared requirements for native-file archival, extract-before-upload, run provenance, checksums, remote verification, and separation of archival from scientific validation are defined by the [central Calculation-result Archival Policy](https://github.com/lm846802297-bit/high-energy-theory-research/blob/14155487afc075d301fa7f6cf452c7bff137e869/docs/CALCULATION_RESULT_ARCHIVAL_POLICY.md). The common run-manifest, result-summary, validation-report, and 14-field source-control contract is defined by [central RCRPS v1.0.1](https://github.com/lm846802297-bit/high-energy-theory-research/blob/14155487afc075d301fa7f6cf452c7bff137e869/docs/STRUCTURED_RESULT_DATA_STANDARD_v1.0.1.md). This document retains only the project's immutable run-bundle and QFT-specific provenance requirements.

## 1. Immutable run bundle

Every substantial numerical run receives a unique immutable directory under results/<run_id>/. Do not overwrite a previous run. For an incremental run, record its parent and create a new run ID.

Keep the following roles separate:

- run_manifest.json: exact source revision, working-tree/source-control context, inputs, configuration, environment, command, outputs, and archive state;
- result_summary.json: physical/numerical results, units, benchmark, uncertainty, scope, and scientific status;
- validation_report.json: checks actually executed, methods, observed values, tolerances, and evidence;
- environment.txt: OS, kernel, CPU/GPU, RAM, interpreter/compiler, dependency versions, and relevant non-secret environment information;
- input_sha256.txt: SHA-256 hashes of all input and configuration files;
- native output files without lossy conversion;
- convergence tables, uncertainty estimates, and independent validation output;
- RESULTS.md: raw observations, derived quantities, interpretation, and unresolved checks.

The final run records must conform to central RCRPS v1.0.1, including all 14 fields in execution.source_control. The manifest records the actual checkout used by the calculation, not the later commit that archives its result.

## 1.1 Historical run records

Existing records created before adoption remain immutable historical records. In particular, the existing `run-manifest.json` files and the Phase-VII `summary.json`, `phase7-environment.txt`, and `phase7-source-manifest.txt` are legacy project records, not retroactive RCRPS v1.0.1 records. They do not contain the full central 14-field `execution.source_control` object and must not be rewritten or relabelled as v1.0.1.

For new substantial runs after adoption, create the central RCRPS v1.0.1 `run_manifest.json`, `result_summary.json`, and `validation_report.json` alongside the existing project-specific progress and native-output records. Where useful, add an index or cross-reference from a new RCRPS record to an older run; do not mutate the old run bundle or fabricate missing source-control fields.


## 2. QFT-specific run metadata

In addition to central RCRPS fields, record the parameters that define this project's calculation:

- spacetime dimension, metric, natural-unit convention, and Fourier-transform convention;
- spatial volume/length, boundary conditions, lattice spacing or momentum-mode cutoff, and retained momentum modes;
- Fock-space symmetry sector, basis dimension, Fock cutoff, and any projection/truncation rule;
- mass, coupling amplitude and time profile, pulse start/end/tail treatment, and time-origin convention;
- normal-ordering reference mass, counterterms, vacuum-energy subtraction, regulator, and renormalization prescription;
- initial state, in/out definitions, observation times, and exact out-mode projection;
- solver, integration method, tolerances, maximum step, time stencil, precision, and random seeds if any;
- operator ordering and branch assignment for each CTP/Wightman component;
- output correlators, cumulants, residual definitions, units, normalization, and data shapes.

Keep finite-volume/Fock-space results distinct from regulator extrapolations and continuum/infinite-volume claims.

## 3. Native outputs and derived data

Preserve native JSON/CSV tables, complex-valued correlator components, Fock-space or Hamiltonian data needed for independent checks, solver diagnostics, and plotting inputs. Keep raw simulation outputs, processed observables, derived cumulants, and figures distinguishable.

Record the source script, exact command, byte size, row/record count, units, and SHA-256 for each output. Preserve sufficient precision for independent recomputation. Do not use a figure or a prose summary as a substitute for its numerical source data.

### 3.1 Three-tier storage selection

Apply the central three-tier policy according to data role, size, access pattern, durability, licensing, and project budget. Do not migrate historical run bundles merely to normalize storage.

- **Ordinary Git:** retain source code, theory notes, input definitions, run manifests, provenance, validation reports, and reasonably sized native CSV/JSON numerical outputs. Preserve exact native precision and immutable run directories.
- **Git LFS:** consider prospectively for a large binary numerical artifact that must be versioned with the source and remains practical to retrieve with the project working set. Record the LFS object identity/OID, SHA-256, byte size, format, run ID, and referencing commit; verify the actual downloaded bytes.
- **External object storage:** use for bulk raw simulation outputs or long-term datasets whose volume or access pattern makes Git/LFS unsuitable. Keep an immutable object/version identifier, exact URI, byte size, source run ID, SHA-256, retrieval instructions, and retention responsibility in a repository-side manifest.

After upload, retrieve each LFS/external object and verify its exact bytes against the recorded SHA-256. Never rely on a temporary runtime workspace or expiring CI artifact as the sole permanent copy. Do not introduce LFS or external storage without a concrete need.

## 4. Validation and scientific status

Validation should be selected according to the claim and may include:

- free-vacuum and zero-coupling limits;
- canonical commutators, Hermiticity, norm conservation, and symmetry-sector checks;
- momentum conservation and spatial translation invariance;
- CTP branch conjugation and largest-time identities;
- independent operator construction and normal-ordering identities;
- time-step/stencil refinement, solver tolerance, Fock-cutoff and momentum-cutoff dependence;
- independent perturbative/Dyson expansion or analytic pulse Fourier-transform checks;
- comparison of independently implemented observables and Wick contractions.

Record only checks actually executed. A schema test, unit test, successful solver run, or checksum match does not by itself establish physical correctness or continuum convergence.

## 5. Archive verification

After pushing, retrieve the remote files and compare paths, byte sizes, row/record counts, and SHA-256. Report one of the actual archive states: local only, committed to PR branch, or merged and verified on main.

Keep archive state separate from execution, numerical validation, and physical validation. The bootstrap environment snapshot describes repository initialization only; it is not the environment of a scientific run.
