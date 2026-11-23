# Adopted Phase VII baseline and current acceptance gate

This is a provenance and integration checkpoint on the existing research
branch, not a separate physics run or a new PR. Earlier results remain
unchanged and do not close the current acceptance gate.

## Verified recovered baseline

- Exact archive: `phase7-cutoff-scan-001.tar.xz`; SHA-256
  `e8e720bcb253ddd9bfd657e42a0902226a5a2dbb07a81c3b39560f1bd78a8b87`,
  matching the digest recorded in the Phase VII PR provenance.
- All 20 entries of the native output manifest match SHA-256. This includes
  manifests, the archive's test log and summary, as well as native cases.
- All three optimized source files match both the archive manifest and
  the PR provenance. The Phase VI dependencies used in replay match the
  two recorded Git blob IDs. The exact archived suite passed 13/13 tests
  on replay; those tests are not a K=5 N=6 propagation result.
- Test identity note: the archived test file hashes to
  `3dffa0f668dc400b4fc76b39423c8d5f25773346b6a926db47dfc2059d02690c`
  as specified by the manifest inside the archive. The shorter PR
  provenance lists `06243d6d1df12efc970aa07391104621aa0c8c1b09e3cd45f7acb7b00035af20`
  instead. The archive SHA pins exactly what was replayed; these test IDs
  must not be treated as identical.

## Evidence that carries forward

The registered Phase VII record includes K=3 N=6 DOP853 and a K=5 N=5
forward point. The restored base is the dependency for subsequent code.
The existing Stage 01 K=5 N=6 DOP853 run captures 16 signed components
at three finite-difference steps and independently validates the cubic
matrix. It does not contain the later compact CF4/Chebyshev state's bytes
or the seven independent propagated moments. The Stage 02 CF4 calculation
is at N=7 and uses Taylor action; it cannot substitute for the N=6
CF4/Chebyshev integrator-step test.

`operator_bridge.py` is a stationary operator compatibility check only:
it compares the recovered optimized builders to the pinned Hamiltonian
and independent sparse Wick product at K=5, N=6. It does not propagate
states or establish independent composite-moment agreement.

## Active gate

1. Recover the compact shared-prefix CF4/Chebyshev runner, its native
   configuration/manifest and any available K=5 N=6 state files. The
   Phase VII archive does not contain these later files.
2. Use those same inputs and states (or continue that exact implementation
   from the last checkpoint). Retain provenance for the four-point
   function and the connected source. Save the shared-prefix contraction
   and independently propagated M33, three M31, and three F2 terms so
   their signed complex values can be compared before subtraction.
3. Separate time-derivative stencil h=0.00625, 0.003125, 0.0015625 from
   CF4 propagation max-step refinement. Report complex residuals by
   branch, minimum comparison scale/absolute errors, and convergence of
   F4c and source with the 1e-3 relative acceptance limits.
4. Do not promote existing higher-cutoff or closure pilots, or update
   PR #8 scientific conclusions, before this gate is closed.

Replay the original Phase VII test archive by extracting it unchanged
to a temporary directory, setting PYTHONPATH to its src and the pinned
repository src, and running `python -m unittest discover -s tests -p test_phase7.py -v`
from the extracted root. The bridge uses the same PYTHONPATH plus the
Stage 01 scripts directory, which supplies the independent Wick builder;
run it from the extracted archive root.

## K=5 N=6 operator bridge result

For D=12,376 at the reference parameters, the maximum absolute element
differences are 2.78e-17 for the interaction matrix and 1.67e-16 for the
cubic matrix against the independent Wick construction. No matrix entry
differs by more than 1e-14. These are operator checks, not propagated
moments. The three exact restored builders are bundled for reuse.
