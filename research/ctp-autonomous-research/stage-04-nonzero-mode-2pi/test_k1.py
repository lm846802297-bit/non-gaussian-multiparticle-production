import unittest
import numpy as np
from k1_bs import solve
from two_pi_k0 import solve as zero
from explicit_wick import wick_by_q
from reference import A

class FourPointResponse(unittest.TestCase):
    def test_zero_mode_reduction(self):
        a=solve(.5,.5,A,K=0)
        b=zero(.5,.5,tmin=-4.)
        for tag in a['F4c']:
            np.testing.assert_allclose(a['F4c'][tag],b['F4c'][tag],rtol=1e-10,atol=1e-13)
    def test_free(self):
        a=solve(0.,.5,A)
        self.assertEqual(max(map(abs,a['F4c'].values())),0.)
    def test_independent_spatial_wick_in_all_Q(self):
        h=.5; lam=.01
        a,b=[solve(s*lam,h,A) for s in (+1,-1)]
        p=wick_by_q(A,h)
        for Q,P in p.items():
            for tag,(c1,c2) in P.items():
                u,v=a['F4c_by_Q'][Q][tag],b['F4c_by_Q'][Q][tag]
                np.testing.assert_allclose((u-v)/(2*lam),c1,rtol=1e-6,atol=1e-11)
                np.testing.assert_allclose((u+v)/(2*lam*lam),c2,rtol=1e-4,atol=2e-11)
    def test_shift_and_contour_identities(self):
        a=solve(.5,.5,A)
        s=solve(.5,.5,[(t,x+3.21) for t,x in A])
        tr=str.maketrans('+-','-+')
        for tag,v in a['F4c'].items():
            np.testing.assert_allclose(v,s['F4c'][tag],rtol=1e-10,atol=1e-13)
            np.testing.assert_allclose(v,np.conj(a['F4c'][tag.translate(tr)]),rtol=1e-10,atol=1e-13)
            np.testing.assert_allclose(v,a['F4c'][tag[0].translate(tr)+tag[1:]],rtol=1e-10,atol=1e-13)
