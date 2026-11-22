# Stage 05: K=6 finite-Fock and two-volume checks

Two completed native DOP853 calculations for configuration A at
`(K,Nmax,L)=(6,7,20)` and `(6,7,24)` are archived under [run/](run/).
Parameters are m=T=1 and λ₀=0.1, initial free vacuum imposed at
t=-12, fixed physical insertion points
`((1,0),(.5,1),(-.5,2),(-1,3))`, and all 16 CTP assignments.
Each Fock basis has 77,520 states and 3,995,140 nonzero interaction entries.
Both use DOP853 with rtol 2e-10, atol 2e-12 and max step 0.1.
Native wall times are 850.9 and 840.5 s, respectively; replay sequentially
unless the available memory is known to suffice for two concurrent runs.
The full Hamiltonian/source implementation is inherited from exact Phase VII
baseline `eab1f5b3198d03eafc3f162299c9038956beaa0b`.

## Finite-Fock check

Values are maximum across CTP assignments, unless noted. Complex signed
residual data and the separate audit are available in
[high_signed_audit.json](high_signed_audit.json). The maximum normalized
analytic residual uses `|R|/(|R-S|+|S|)`, where
`S = λ(t₁) S₄,c / 6` and `R` is the signed equation
residual. The analytic temporal derivative inserts
`−[H,[H,φ]] + i[Ḣ,φ]`, avoiding the temporal stencil but sharing
the propagated DOP853 states. It is not an independent K=6/N=7 propagation
method. The auditor recomputes each branch's complex Richardson combination,
not an extrapolation of unsigned maxima.

| L at K=6 | Analytic residual N=6 | Analytic residual N=7 | Δ F4c N=6→7 | Δ connected composite source N=6→7 |
|---:|---:|---:|---:|---:|
| 20 | 9.616e-6 | 1.774e-6 | 3.000e-5 | 3.828e-6 |
| 24 | 8.111e-6 | 1.908e-6 | 3.056e-5 | 3.901e-6 |

The N=7 stencil runs at `h=(.00625,.003125,.0015625)` give maximum
normalized residuals `(5.112e-5,1.315e-5,3.894e-6)` at L=20 and
`(4.041e-5,1.045e-5,3.882e-6)` at L=24. The signed independent audit
checks metadata, all branches, complement conjugation and the latest-time
relation, unnormalized versus normalized residuals, and complex
Richardson comparisons. Small fixed-K Fock differences are not a bound on
momentum-cutoff or volume dependence.

## Separate regulator comparisons at N=7

Here `Λ=2πK/L` labels the largest retained mode, with the
endpoints `|n|≤K` included. The physical spacetime coordinates of the
insertions are held fixed. Relative differences use the first point as
denominator and are maximized across branches.

| Comparison (K,L), first → second | Regulator change | Max Δ F4c | Max Δ connected composite source |
|---|---|---:|---:|
| (5,20) → (6,20) | fixed L; changes Λ | 35.73% | 47.35% |
| (5,20) → (6,24) | matched Λ=1.5707963 | 4.18% | 5.69% |
| (6,20) → (6,24) | fixed integer K; changes both L and Λ | 34.28% | 41.75% |

The matched ratio controls the highest retained physical momentum, but
changes mode spacing and the number of retained modes. No continuum or
infinite-volume limit follows from these points, and the fixed-K row is
not an isolated volume change. Increasing K and adding matched-cutoff
volume points remain open. No K=6/N=7 CF4 or other independent time
propagator run is claimed here.

## Replay and provenance

Verify bytes with `sha256sum -c SHA256SUMS.txt`. The two existing imports
`src/sd_residual.py` and `src/ctp_sd_all_branches.py` match the baseline
Git blobs `29139f6f8a35d3d874c401d240318ed2c8cadc40` and
`4a829bec5cf295032808ab28f13684d86bea5bea`. Earlier inputs are linked
rather than recopied:

- [Stage 01 K=5 N=6 A](../stage-01-k5-fock/results/k5-n6-segmented.json) and
  [K=5 N=7 A](../stage-01-k5-fock/results/k5-n7.json).
- [Stage 02 K=6 N=6 L=20](../stage-02-ctp-regulators/scan/run/k6-n6-A.json)
  and [K=6 N=6 L=24](../stage-02-ctp-regulators/scan/run/k6-n6-A-L24.json).
- [Stage 02 configuration C inputs](../stage-02-ctp-regulators/scan/run/)
  are used by the general-purpose `verify_high.py`.

Run `bash replay-audits.sh /path/to/new-output-dir` from the checkout;
it links the inputs read-only into a temporary directory and byte-compares
the regenerated audits against the included originals.

For a full numerical rerun, start at the repository root and choose a new
results directory to preserve the archived files:

```sh
p=research/ctp-autonomous-research/stage-05-k6-fock-volume
export PYTHONPATH="$PWD/src"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p /tmp/k6-stage05-replay
python "$p/next_scan.py" 6 7 /tmp/k6-stage05-replay/k6-n7-A.json --length 20 --config A --rtol 2e-10 --atol 2e-12 --step .1
python "$p/next_scan.py" 6 7 /tmp/k6-stage05-replay/k6-n7-A-L24.json --length 24 --config A --rtol 2e-10 --atol 2e-12 --step .1
```

Native runs used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0; their recorded
wall times prevent byte-identical JSON comparisons after a full
recomputation. Replay of the existing signed audit, however, is byte exact.
