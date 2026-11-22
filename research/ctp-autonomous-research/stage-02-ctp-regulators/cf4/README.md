# Independent CF4 cross-check, K=5 / Nmax=7 / configuration A

This is a research-side calculation, not a modification to main or the pinned
upstream source. The Hamiltonian, one-field operator and exact spatial
derivative are evaluated with the same pinned finite-volume Hamiltonian and
operator conventions as the DOP853 benchmark. The cubic source uses the
independent sparse Wick construction. All-forward connected moments are
evaluated explicitly from the three pair subtractions, independently from
the production CTP subtraction functions.

The independent propagator uses a symmetric two-exponential fourth-order
commutator-free Magnus scheme at the two Gauss-Legendre nodes. Each sparse
exponential action uses a Taylor series with an a priori Hermitian spectral
norm bound, not the reference DOP853 propagation. Matrix-vector arithmetic is
complex binary64; the tolerance bounds Taylor truncation, not roundoff.
Only the common physical-time surface t=1 is needed for the central moments.
All reported h values for the second-derivative stencil are distinct from
the maximum integration time step; the CF4 run refines the latter.

From this directory, with an exact copy of the baseline source and the
independent_source.py dependency under /tmp/ctpproject:
    
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
    export PYTHONPATH=/tmp/ctpproject/src:/tmp/ctpproject
    python -m unittest -v test_cf4
    python -u cf4.py --baseline /tmp/ctpproject/run/k5-n7.json \
        --steps .1 .05 --output cf4_k5_n7_A.json

The benchmark JSON supplies an immutable SHA-256 in the CF4 JSON, rather than
requiring GitHub network access during replay. No continuum or infinite-volume
limit has been tested by this comparison.
