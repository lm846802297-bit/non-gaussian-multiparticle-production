import unittest
import numpy as np
from independent_source import setup, wick_matrix, summarize
from sd_residual import projected_normal_ordered_cubic
class Check(unittest.TestCase):
    def test_against_monomials_at_nonzero_x(self):
        b, p, i, m = setup(2, 4)
        a = wick_matrix(b, p, i, m, 2, 4, 1, 20, 1.7)
        ref = projected_normal_ordered_cubic(tuple(m),list(map(tuple,b[:a.shape[0]])),{tuple(s):j for j,s in enumerate(b[:a.shape[0]])},1,20,1.7)
        self.assertLess(summarize(a-ref)['max_abs'], 7.e-16)
    def test_vacuum_no_one_particle_matrix_element(self):
        b,p,i,m=setup(2,4)
        a=wick_matrix(b,p,i,m,2,4,1,20,1.7)
        self.assertLess(np.max(np.abs(a[:6,0].toarray())),7.e-16)
if __name__=='__main__': unittest.main()
