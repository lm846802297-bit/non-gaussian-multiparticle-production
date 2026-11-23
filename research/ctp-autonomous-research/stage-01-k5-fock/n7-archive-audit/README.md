# K=5, Nmax=7, configuration A: archival audit only

This complements the fixed-Nmax=6 numerical acceptance on the same research
branch. It performs no time evolution and checks only byte-identified records
already present before this audit. Source inputs, Git blobs and SHA-256
are recorded in the JSON. A matched Fock change is not the same as a
fixed-cutoff solver comparison.

- Archived DOP853: 16 signed branches and all three FD steps, with
  complement conjugation, latest-time and residual-normalization checks
  passing. The worst relative no-stencil SD residual is 1.809e-6. The
  finite-difference maxima are 4.156e-5, 1.067e-5, and 3.298e-6. Signed
  branchwise slopes against the existing analytic residual vary,
  particularly at the finest stencil; no uniform second-order claim.
- Independent CF4/Taylor includes two propagation max steps, .10 and .05,
  but only one FD h=.003125 and the all-forward (++++) component. At .05,
  relative differences vs DOP853 are 1.949e-8 for F4c and 3.400e-10 for
  the connected source (both below 1e-3 in ++++ only). Its complex
  analytic SD residual differs from DOP853 by 2.03e-16 absolute.
- The Wick-source matrix check has max absolute discrepancy 1.665e-16
  with zero entries above 1e-14. It is not independent propagated-moment
  evidence.
- The matched A-case Nmax=6→7 maximum relative change across the 16
  branches (reference Nmax=6) is 3.030e-5 for F4c and 3.727e-6 for
  the connected source. These are adjacent finite-cutoff displacements
  only, not error bounds or evidence of continuum convergence.

**Open before acceptance at the N6 standard:** independent seven signed
primitive moments for every branch (the existing N7 CF4 files retain the
three products, not the separate factors); independent all-branch CF4
stencils and single-vector bridge at N7. The DOP853 N7 archive does not
contain the shifted correlators needed to reconstruct d2F4c independently
from its residual. Preserve the current data; do not rerun its existing
DOP853 sweep merely to repeat the reported numbers.

Replay from a clone of the research branch, with the checkout root as the
first argument and an output outside the repository:
```sh
python research/ctp-autonomous-research/stage-01-k5-fock/n7-archive-audit/audit_n7.py "$PWD" /tmp/n7-audit-replay.json
cmp /tmp/n7-audit-replay.json research/ctp-autonomous-research/stage-01-k5-fock/n7-archive-audit/audit.json
```
