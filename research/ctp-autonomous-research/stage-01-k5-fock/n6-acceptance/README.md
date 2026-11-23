# Integrated K=5, Nmax=6 numerical acceptance, configuration A

This continues the restored Phase VII baseline on the existing research
branch. It does not claim to have recovered the later historical runner.
The CF4/Chebyshev adapter was reconstructed around the recovered operators;
the existing DOP853 acceptance snapshot was not regenerated.

Physical setup: m=1, L=20, lambda0=.1, T=1, t_init=-12, 12,376 Fock states,
events (1,0), (.5,1), (-.5,2), (-1,3). All 16 signed complex contour
components are retained. There are no coincident-time contact terms.

## Independent paths

- CF4 Gauss-node propagation with Chebyshev/Bessel action and a shared
  chronological prefix tree, advancing both branches only up to the latest
  insertion. The two-mesh dataset records its block scheduling explicitly.
- DOP853 direct insertion paths use the pinned reference builders and do
  not share the CF4 state or the shared-prefix contraction. Independently
  log signed M33, three M31, three F2; reconstruct each disconnected
  product and the connected source only after those values are preserved.
- Distinguish the finite-difference time step h=0.00625, 0.003125,
  0.0015625 from the CF4 internal maximum step dt=0.10 or 0.05.
  Compare the signed SD residual to the finite-matrix analytic insertion
  (i[\dot H,\phi]-[H,[H,\phi]]), and do not estimate stencil order from
  ratios of maxima of the normalized SD residual. The residual has a
  nonzero finite-Fock contribution.
- A third check reads, but does not regenerate, the Stage 01 DOP853
  snapshot. The acceptance thresholds below concern numerical agreement
  at fixed Nmax. The Nmax=5→6 physical Fock shift is recorded separately.

## Replay and artifacts

From the root of the pinned research branch:
```sh
root="$PWD"
here="$root/research/ctp-autonomous-research/stage-01-k5-fock/n6-acceptance"
phase7="$root/research/ctp-autonomous-research/stage-01-k5-fock/phase7-reconciled/src"
baseline="$root/research/ctp-autonomous-research/stage-01-k5-fock/results/k5-n6.json"
export PYTHONPATH="$phase7:$root/src:$here/cf4:$here/direct:$here/audit"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
```
The component run commands follow in their README files. Run any
reproduction in a fresh scratch directory; committed native files should
never be overwritten.

The native commands are:
```sh
python "$here/cf4/gate_cf4.py" --K 5 --N 6 --steps .10 .05 --baseline "$baseline" --output k5-n6-cf4.json --checkpoints
python "$here/cf4/verify_vector_mode.py" k5-n6-cf4.json k5-n6-cf4-dt0.05-states.npz k5-n6-cf4-single-vector-bridge.json
python "$here/direct/direct_moments.py" k5-n6-direct-moments.json --archived "$(dirname "$baseline")/k5-n6-segmented.json"
```
For replay of the signed audit with the original archive hash, supply the
byte-identical Phase VII tar.xz (digest in the reconciled baseline notes)
and use the file paths accepted by `audit/gate_audit.py`. The portable
replay helper is explained in the audit notes. The two large NPZ checkpoints
are preserved separately; the native CF4 JSON records each size and
SHA-256. They can also be recreated with the above `--checkpoints` option.
The digest of the preserved checkpoint ZIP is
`b1d940d5b7a8292b37e0f991c8a3cac56a321dd03fd124605c52c2bf6336932c`.

## Results and interpretation

- At fixed K=5, Nmax=6 the fine CF4 versus independently propagated DOP853
  max relative complex differences are 3.42e-8 for F4c and 3.40e-10
  for the connected source. Both are below 1e-3. Reducing maximum
  propagation step .10→.05 changes the same quantities by at most
  4.96e-7 and 5.54e-9. The seven signed primitive moments agree to
  9.90e-12 maximum absolute and 2.18e-7 maximum relative across all
  112 entries. Products and connected subtraction are separately audited.
- A true single-vector CF4 shared-prefix replay at dt=.05 reproduces all
  24 stored central prefix/field/source states (L2 differences zero at
  binary64) and the raw moment values to 2.12e-16 absolute. This bridge
  is central-only; the three stencils were collected in the block run.
- At dt=.05 the max normalized finite-difference SD residuals are
  4.46e-5, 1.36e-5, 8.16e-6; the no-stencil maximum is 7.99e-6.
  The signed analytic residual comparison to archived DOP853 has a
  6.14e-16 maximum absolute difference. The finest stencil shows
  cancellation noise and does not have uniform branchwise second order.
  No physical cutoff error threshold is imposed on that residual.
- Comparing archived Nmax=5→6 in the all-forward component at the same
  physical point gives 2.24e-3 for F4c and 3.69e-5 for the source.
  Fixed-cutoff numerical acceptance is not evidence of convergence
  with respect to the Fock cutoff or continuum/volume limits.

To rerun the signed audit and its 11 failure-path controls in a fresh
temporary fixture, use:
```sh
bash "$here/audit/run-from-branch.sh" /path/to/phase7-cutoff-scan-001.tar.xz /tmp/n6-acceptance-replay.json
cmp /tmp/n6-acceptance-replay.json "$here/audit/gate_validation.json"
```
