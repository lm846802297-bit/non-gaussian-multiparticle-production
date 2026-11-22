import unittest
import numpy as np
from four_pi_k0 import solve
from born_bubble import compute

class VertexChecks(unittest.TestCase):
    def test_free(self):
        f=solve(0.,.5,tmin=-4.)
        self.assertLess(max(map(abs,f['F4c'].values())),1e-14)
        self.assertLess(f['G0_delta'],1e-14)
    def test_weak_and_identities(self):
        lam=.5
        f=solve(lam,.25,tmin=-4.)
        p=compute(0,20.,1.,1.,-4.,[(t,0.) for t in (1.,.5,-.5,-1.)],42)
        for k,(a,b) in p.items():
            v=lam*a+lam*lam*b
            np.testing.assert_allclose(f['F4c'][k],v,rtol=.09,atol=2e-9)
            tr=str.maketrans('+-','-+')
            np.testing.assert_allclose(f['F4c'][k],np.conj(f['F4c'][k.translate(tr)]),rtol=1e-10,atol=1e-12)
            np.testing.assert_allclose(f['F4c'][k],f['F4c'][k[0].translate(tr)+k[1:]],rtol=1e-10,atol=1e-12)
