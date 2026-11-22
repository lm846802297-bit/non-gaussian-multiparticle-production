import unittest
import numpy as np
import scan
from ctp_sd_all_branches import evaluate_all
class ScanTest(unittest.TestCase):
 def test_nonzero_coupling(self):
  a=scan.run(1, 5, [0.003125], .08, 2e-10, 2e-12)["rows"][0]["components"]
  b=evaluate_all(kcut=1,nmax=5,h=.003125,rtol=2e-10,atol=2e-12,max_step=.08)["components"]
  for k,v in a.items():
   for q in ("F4c","M33","connected_composite_source","SD_residual"):
    np.testing.assert_allclose(v[q],b[k][q],atol=7e-10,rtol=2e-8)
 def test_free_derivative(self):
  a=scan.run(1, 5, [.003125], .1, 2e-10, 2e-12, lam=0)
  self.assertLess(a["max_analytic_abs"],1e-10)
