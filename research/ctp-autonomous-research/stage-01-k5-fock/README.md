# Autonomous research record: non-Gaussian CTP four-point dynamics

This development branch builds from Phase VII draft PR #8, exact commit
`eab1f5b3198d03eafc3f162299c9038956beaa0b`. The pre-existing default branch,
historical data, and draft PR are unchanged. Each stage has raw JSON and a
file digest manifest. Scientific status is stage-specific; finite-cutoff
tests are not a continuum limit or an nPI comparison.

## Stage 01: K=5 finite-Fock check

The comparison holds m=1, L=20, T=1, λ₀=0.1, the free-vacuum normal-ordering
prescription and periodic mode set fixed. Configuration A events (t,x) are
(1,0), (0.5,1), (-0.5,2), (-1,3). The contour insertion is differentiable
near the latest event; no contact terms are included.

- Independent finite-matrix time derivative i[Ḣ,φ] − [H,[H,φ]] separates
  stencil error from the residual in the projected Fock space.
- At h=.003125, K=5 N=6 gives max relative SD residual 1.36e−5 at both
  recorded DOP853 tolerance choices. The prior reported ~1.10e−3 was not
  reproduced by this archived baseline. The unarchived runner
  and its native output were unavailable, so the discrepancy is not assigned
  to a cause.
- Analytic derivative max relative residual drops from 7.99e−6 to 1.81e−6
  from N=6 to 7. Largest relative changes across branches are 3.03e−5
  (F4c) and 3.73e−6 (connected six-field source). These observations are
  limited to the stated test configuration.
- Independent projected sparse Wick matrices agree at both N=6,7 with
  the upstream normal-ordered monomial builder to 1.67e−16 maximum absolute
  element difference.
- Native data includes 3 stencil sizes, complex signed residuals, and all
  16 branches; the archived stage passed four supplemental and eight prior
  targeted tests. A smaller step is not an independent integrator.
