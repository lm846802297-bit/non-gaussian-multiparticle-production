# Provenance and Archival Policy

Every numerical run receives a unique immutable directory under `results/` containing:

- `run_manifest.json`: run ID, timestamp, git commit, dirty-tree flag, command, model parameters, lattice/volume parameters, solver settings, tolerances, random seeds, and output file list;
- `environment.txt`: OS, kernel, CPU/GPU, RAM, compiler/interpreter, package versions, and relevant environment variables (never secrets);
- `input_sha256.txt`: SHA-256 hashes of all input and configuration files;
- native output files without lossy conversion;
- convergence tables and error estimates;
- independent validation output;
- a short `RESULTS.md` separating raw observation, derived quantity, interpretation, and unresolved checks.

Rules:

1. Never overwrite raw inputs or prior run outputs.
2. Preserve exact source revisions and input hashes.
3. Record actual hardware and execution parameters, not assumed defaults.
4. Keep generated plots traceable to raw data and plotting scripts.
5. Do not treat a passing software test as physical validation.
6. Do not merge a scientific result into `main` until the corresponding PR has been reviewed.
7. Keep secrets and private credentials out of the repository.

The bootstrap environment snapshot records the environment in which this archive was initialized; it is not a scientific run environment.