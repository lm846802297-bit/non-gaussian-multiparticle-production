# Stage 02: independent propagation and regulator sensitivity

All data below are finite-regulator results with the Phase VII Hamiltonian
at commit eab1f5b3198d03eafc3f162299c9038956beaa0b. Previous raw files are
immutable. No high-cutoff convergence, continuum extrapolation, or higher-K
self-consistent nPI calculation is asserted.

- **Independent propagator:** for K=5 N=7, A, `++++`, fourth-order CF4 at
  maximum integration step .05 matches DOP853 by 1.95e-8 (F4c), 4.47e-10
  (M33), and 3.40e-10 (connected source) relatively. No-stencil SD absolute
  difference is 2.03e-16. Other branches are DOP853 only; CF4 step is
  distinct from the derivative-stencil h.
- **Configuration C:** (1.2,0), (.4,.7), (-.3,1.7), (-1.1,3.2).
  The no-stencil maximum relative SD residual falls from 2.58e-6 to 5.67e-7
  across N=6→7. N=7 largest relative changes of F4c and source against N=6
  are 2.12e-5 and 5.31e-6 respectively, across 16 branches.
- **Exploratory CTP regulator check (N=6; not Fock-extrapolated at K=6):**
  at fixed L=20, K=5→6 can change F4c and source by up to 35.7% and 47.4%
  across branches. At matched physical cutoff K/L=1/4, L=20→24
  changes them by up to 4.18% and 5.69%. These two cases are not yet a
  cutoff or volume extrapolation.
- **Out-mode observables:** eleven independent native runs in the physical
  even-parity zero-momentum sector. N=6 and N=7 then coincide; the real
  truncation test is N=6→8 (max κ4 drift 5.84e-4 relative on this grid).
  At fixed cutoff the zero and opposite edge pairs admit many partner
  channels and Lκ4 varies by about .04% and .13%; a same-sign edge pair
  has one Born partner channel, so it instead carries L^-2 scaling.
  These numerical patterns are scoped to the reported boxes and mode pairs.
- **Pre-closure perturbative gate:** explicit A contour integration with all
  three bubble channels decreases the worst relative F4c discrepancy
  vs K=5 N=7 Fock from 9.44e-3 (Born) to 1.36e-4 (through O(λ²)).
  This is not a 2PI or 4PI numerical solution.

## Replay
From the repository root, let BASE be the absolute path to
research/ctp-autonomous-research/stage-01-k5-fock and let S2 be the absolute
path to this directory. Keep native outputs immutable by using fresh output
paths. Set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1.

```sh
export PYTHONPATH="$PWD/src:$BASE/scripts"
python "$S2/scan/next_scan.py" 5 7 /tmp/replay-ctp-C.json --config C
python "$S2/scan/next_scan.py" 6 6 /tmp/replay-ctp-L24.json --length 24
python "$S2/cf4/cf4.py" --baseline "$BASE/results/k5-n7.json" --steps .1 .05 --output /tmp/replay-ctp-cf4.json
python "$S2/perturbation/born_bubble.py" --fock "$BASE/results/k5-n7.json" --orders 40 48 56 /tmp/replay-ctp-bubble.json
```
For `out/replay.sh`, copy the out directory to a scratch folder before running
and set ARCHIVED_SRC to the repository src path. The stage manifest covers all
files below, and the component manifests preserve the original sub-run digests.
