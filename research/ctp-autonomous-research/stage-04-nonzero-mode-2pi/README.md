# K=1 nonzero-mode 2PI + Bethe–Salpeter benchmark

An actual real-time three-loop 2PI solution and its consistent, differentiated
Bethe–Salpeter (BS) response, compared to the same projected Fock Hamiltonian.
This is **not** 4PI at K=1 and does not estimate a continuum limit or inclusive
out-particle cumulants. The spatial regulator is only K=1, with momenta
0, ±2π/L. Initial vacuum is imposed at (t=-4), not (t=-12); absolute
time-tail systematics are not yet determined. Both methods share this initial
condition. Units (m=T=1,L=20).

Configurations (time, position):
- A `((1,0),(.5,1),(-.5,2),(-1,3))`.
- C `((1.2,0),(.4,.7),(-.3,1.7),(-1.1,3.2))`.

## Methods and nonzero pair-momentum check

Use `phi(x)=L^(-1/2) sum_k exp(i*p_k*x) q_k`, (g(t)=λ(t)/L),
(G_k(z,w)=\langle T_C q_k(z)q_{-k}(w)\rangle) with no external factor of
(i), and (int_C=int_+-int_-).
The background (G) is diagonal in the translational sense. Its Dyson
insertion includes both the time-local tadpole, (-ig(z)\sum_p[G_p(z,z)-1/(2ω_p)]/2),
with the contour delta, and the nonlocal sunset
(-g(z)g(w)\sum_{a+b+c=k}G_aG_bG_c/6). This uses finite-mode normal
ordering before comparison; no fitted mass counterterm is used.

For a bilocal source, let (X_{kl}) denote the complete normalized response.
(Q=k+l) stays fixed, but (X_{kl}) is not artificially restricted to
(l=-k). For each Q, the functional derivative of the local insertion
is proportional to (sum_{r+s=Q}X_{rs}(z,z)). Its nonlocal term is
(-g(z)g(w)\sum_{a,b}G_aG_b X_{k-a-b,l+a+b}/2), retaining only the mode indices
within cutoff. Both contour legs carry their branch orientation.
Iterate (X_{kl}=X_{0,kl}+G_k K'_{kl}[X]G_l); (X_0) is the sum
of the two cross-side Gaussian contractions with dressed (G). Connected
response is (X-X_0), Fourier transformed and divided by (L^2) for four
physical external fields. All 5 Q blocks and 16 contour labels are retained.
See [Berges (2004)](https://arxiv.org/abs/hep-ph/0409233) and
[Carrington et al., arXiv:1310.4352, Eqs. 14–19](https://arxiv.org/abs/1310.4352)
for the bilocal-response framework.

The independent Born/one-loop computation uses explicit spatial quadrature
for the weak checks, with all three channels. Holding the right pair fixed
and translating both left insertions projects Q by a discrete Fourier
transform. Thus each Q is checked against an independent coordinate-space
calculation, not merely against an index identity in the 2PI code.
In addition, `reference-A-Q2.json` contains a separate Fock DFT extraction
for A and λ0=2: 5 versus 7 shifts agree to (2.05×10^{-18}) absolute, and
Q reconstruction agrees with the full Fock value to (1.40×10^{-18}).
At the same coupling the corrected maximum relative differences for Q=0,
|Q|=1, |Q|=2 are 6.020e-3, 7.115e-3, 8.194e-3. Their maximum absolute
errors are 9.83e-7, 5.60e-7, 1.46e-7 respectively: the largest relative error
is not the dominant absolute error.

## Results

The error is the maximum of the 16 complex relative norm errors compared to
Fock Nmax=14 with tighter DOP853 tolerances. All raw meshes are archived; the
central estimates replace discrete Born+one-loop terms by independently
converged Gauss/Duffy integrals and extrapolate the remaining part in h².
"Fine-pair shift" denotes the difference of the two finest corrected meshes
relative to the finest value, not a rigorous truncation bound.

| Geometry | lambda0 | One-loop fixed order | 2PI+BS | Fine-pair shift |
|---|---:|---:|---:|---:|
| A | 0.1 | 2.925e-5 | 1.798e-5 | 3.71e-8 |
| A | 0.5 | 7.245e-4 | 4.391e-4 | 9.25e-7 |
| A | 2 | 1.122e-2 | 6.560e-3 | 1.47e-5 |
| A | 5 | 6.629e-2 | 3.668e-2 | 9.05e-5 |
| C | 0.1 | 3.106e-5 | 1.890e-5 | 2.42e-8 |
| C | 0.5 | 7.640e-4 | 4.603e-4 | 6.01e-7 |
| C | 2 | 1.178e-2 | 6.722e-3 | 9.45e-6 |

For λ0=2 the Nmax=12→14 relative Fock-reference variation is about 1.1e-8
for both configurations; it is below 1.8e-6 even at λ0=5 in A. All-branch
CTP conjugation and latest-time identities pass on the stored runs.
2PI+BS improves over strict one loop across the cases tested, but the
remaining error grows with coupling. These statements are about the
connected four-point cumulant at the stated finite regulator and samples.

Next tests: independent t_init extension at K=1; consistent K=1 4PI for both
geometries; higher K and separate L tests; and an out-mode projection
validated against direct Hamiltonian expectations before comparing
inclusive non-Gaussian spectra.

## Replay

Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 were used. Work in this directory.
`sd_residual.py` is the unchanged pinned Fock builder from project baseline
`eab1f5b3198d03eafc3f162299c9038956beaa0b`. Its original Git blob identity
can be checked with `git hash-object sd_residual.py` against the baseline.
`born_bubble.py` and `two_pi_k0.py` are copies of the previous validated
closures package. The independent zero-mode implementation is for test
reduction, not used in nonzero-mode solves.

```sh
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
python -m unittest -v test_k1
python reference.py reference-A-high.json --case A
python reference.py reference-C-high.json --case C --couplings .1 .5 2
python pilot.py reference-A-high.json k1-2pi-A.json --h .25 .125 .0625 .03125
python pilot.py reference-C-high.json k1-2pi-C.json --h .1 .05 .025 --couplings .1 .5 2
python exact_q.py reference-A-Q2.json
python q_summary.py k1-2pi-A.json reference-A-Q2.json k1-2pi-A-Q2.json
```

Native JSON retains complex values as `[real,imag]`. Wall times and runtime
versions are recorded where applicable. `SHA256SUMS.txt` authenticates the
files in the supplied execution, not byte-identical replay with timing.
