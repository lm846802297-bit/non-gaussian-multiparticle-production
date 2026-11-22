# Non-Gaussian CTP dynamics: autonomous research checkpoint

Baseline: Phase VII draft PR #8, commit
`eab1f5b3198d03eafc3f162299c9038956beaa0b`.
Results here are isolated on the research branch; no change is made to main.

## Published stages

1. [K=5 finite-Fock/time-derivative audit](stage-01-k5-fock/README.md)
   ([replay](stage-01-k5-fock/REPLAY.md),
   [SHA-256](stage-01-k5-fock/SHA256SUMS.txt)).
2. [Independent CF4, configuration C, and regulator checks](stage-02-ctp-regulators/README.md)
   ([SHA-256](stage-02-ctp-regulators/SHA256SUMS.txt)).
3. [Single-mode nPI comparison pilot](stage-03-single-mode-npi/README.md)
   ([SHA-256](stage-03-single-mode-npi/SHA256SUMS.txt)).
4. [Nonzero-mode K=1 three-loop 2PI/BS comparison](stage-04-nonzero-mode-2pi/README.md)
   ([SHA-256](stage-04-nonzero-mode-2pi/SHA256SUMS.txt)).
5. [K=6 Fock-truncation and two-volume audit](stage-05-k6-fock-volume/README.md)
   ([SHA-256](stage-05-k6-fock-volume/SHA256SUMS.txt)).

## Remaining gates

- Stage 05 verifies K=6 N=6→7 at L=20 and L=24 with signed per-branch
  residuals. An independent K=6 N=7 propagator, higher momentum cutoffs,
  and further matched-cutoff volumes remain open. The existing two-volume
  difference is not an infinite-volume extrapolation; retain the endpoint
  prescription and mode normalization when matching K/L.
- Cross-check the closure implementation on a nonzero-mode case; exploit
  total pair momentum so that solving a full dense (branch,time,mode)^4
  vertex is not the first computational route. Compare regulated F4c at
  common events and apply the same out-mode projection before making
  statements about κ4.
- Keep numerical convergence, regulator dependence, and closure error
  distinct. Archive native datasets, inputs, replay, and digests as separate
  checkpoints. Mark any conclusion that cannot be checked as pending.
