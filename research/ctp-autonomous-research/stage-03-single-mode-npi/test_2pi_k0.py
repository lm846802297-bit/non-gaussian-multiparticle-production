import unittest
import numpy as np
from two_pi_k0 import solve, EV
from born_bubble import compute

class WithBilocalSource(unittest.TestCase):
    def test_free(self):
        r=solve(0.,.25)
        self.assertLess(max(map(abs,r['F4c'].values())),1e-14)
        self.assertLess(r['G0_delta'],1e-14)
    def test_weak_matches_diagrams(self):
        t=[(v,0.) for v in EV]; split=compute(0,20.,1.,1.,-12.,t,42)
        l=.5
        p={m:l*a+l*l*b for m,(a,b) in split.items()}
        r=solve(l,.125)
        for k,v in p.items():
            np.testing.assert_allclose(r['F4c'][k],v,rtol=3e-2,atol=2e-9)
            tr=str.maketrans('+-','-+')
            np.testing.assert_allclose(r['F4c'][k],np.conj(r['F4c'][k.translate(tr)]),rtol=1e-10,atol=1e-12)
            np.testing.assert_allclose(r['F4c'][k],r['F4c'][k[0].translate(tr)+k[1:]],rtol=1e-10,atol=1e-12)
