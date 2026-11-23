# Independent direct CTP signed-moment check

Case: exact Phase VII pinned Hamiltonian, K=5, Nmax=6, m=T=1, L=20,
λ₀=.1, vacuum at t=-12, configuration A
`((1,0),(.5,1),(-.5,2),(-1,3))`. Operator module is
`src/sd_residual.py` from commit `eab1f5b3198d03eafc3f162299c9038956beaa0b`,
exact Git blob `29139f6f8a35d3d874c401d240318ed2c8cadc40`.
The source is the pinned projected normal-ordered cubic at event 0.

This DOP853 runner independently develops each complete chronological
insertion path, with no shared-prefix tree and without using any CF4
checkpoint state. All paths begin in the same physical state at the
earliest event (pre-evolved from the vacuum); only identical completed
paths are cached, not intermediate states between insertions. The bra and
ket are evaluated up to the latest time of the required observable,
then contracted in the CTP convention; no conjugate label is filled by
symmetry. Conjugation and the latest relevant time relation are checked
afterwards. The Hamiltonian and field/source operators are identical to
the pinned model, not rederived here.

For each of the 16 CTP assignments the raw signed complex observables
are `M33`, `M31[01,02,03]` (event 0 cubic and event j), and
`F2[23,13,12]` (complementary pair).
Then derive the three ordered M31×F2 terms, connected composite source
and F4c independently; raw F4 and the other F2 pairs are included.
JSON complex values are [real,imag]. This is an independent
propagation/contraction check against CF4, not a standalone test of
the discrete commutator expression or the unrecovered Phase VIII code.

## Reproduce

Run in this directory; `pin` is the root of the exact archived operator
source checkout, and `reference` is `k5-n6-segmented.json` from stage 01.

```sh
export PYTHONPATH="$pin/src"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -m unittest -v test_direct
python direct_moments.py raw.json --archived "$reference"
python verify_direct.py raw.json "$pin/src/sd_residual.py" "$reference" --out direct-validation.json
```

The weak/free controls compare a small interacting case against the shared
archive, and compare the free projected cubic contraction against an
independent triple-momentum Wick oracle with all CTP labels. A complete
native replay records runtime and will not be byte-identical by design.

## Signed comparison with independent CF4

The companion native CF4 run is `k5-n6-cf4.json` in the separate CF4
archive. It is not recopied here. Compare both integrations with:

```sh
python compare_cf4.py k5-n6-direct-moments.json /path/to/k5-n6-cf4.json direct-cf4-signed-comparison.json
```

The script checks that both calculations use the same case, spacetime
insertions, initial time, archived Git pin and pinned operator module.
It archives the real and imaginary differences separately before
reporting norms. At fine CF4 step .05, the largest of 112 raw complex
discrepancies is `9.896e-12` for M31:01 (relative `2.175e-7`).
The maximum differences for M33, the connected source, and F4c are
`2.008e-12`, `1.250e-12`, and `1.466e-12`, respectively.
The largest raw difference decreases by approximately 16 times when
the CF4 step is halved from .1 to .05. This does not reproduce the
unrecovered Phase VIII code or serve as an independent analytic
equation-of-motion implementation.
