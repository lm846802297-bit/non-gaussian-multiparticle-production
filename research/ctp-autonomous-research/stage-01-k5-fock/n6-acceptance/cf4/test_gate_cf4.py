import unittest
import numpy as np
from scipy.sparse import csr_matrix
from scipy.integrate import solve_ivp
from gate_cf4 import ChebCF4, build, measure, EV, PHYS
from ctp_sd_all_branches import evaluate_all
class TestCF4(unittest.TestCase):
 def test_cheb_and_cf4_order(self):
  E=np.array([.2,2.7]);V=csr_matrix([[.7,-.84],[-.84,.3]]);v=np.array([.8,.6j])
  profile=lambda t:.4+np.sin(1.5*t)
  P=ChebCF4(E,V,profile)
  ref=solve_ivp(lambda t,w:-1j*(E*w+profile(t)*(V@w)),(-1.,1.4),v,method='DOP853',rtol=3e-13,atol=3e-15).y[:,-1]
  a=P(v,-1.,1.4,.12);b=P(v,-1.,1.4,.06)
  ratio=np.linalg.norm(a-ref)/np.linalg.norm(b-ref)
  self.assertTrue(15.5<ratio<16.5,ratio)
  self.assertLess(abs(np.linalg.norm(b)-1),2e-14)
  # Time reversal follows from symmetric CF4 and Hermitian action.
  np.testing.assert_allclose(P(b,1.4,-1.,.06),v,rtol=0.,atol=1e-14)
 def test_all_branches_small_exact_dop853(self):
  c=build(1,5,PHYS,EV)
  obs=measure(c,.04,None,hs=[.00625],verbose=False)
  ref=evaluate_all(kcut=1,nmax=5,rtol=2e-11,atol=2e-13,max_step=.06,h=.00625)
  got=obs['central_observables']
  for tag,b in ref['components'].items():
   a=got[tag]
   for key in ('F4c','M33','connected_composite_source'):
    self.assertLess(abs(complex(*a[key])-complex(*b[key])),2e-10)
   for j in range(3):
    self.assertLess(abs(complex(*a['M31_F2_terms'][j])-complex(*b['M31_F2_terms'][j])),2e-10)
  for x in obs['central_checks'].values():self.assertLess(x,1e-14)
 def test_zero_pulse_exact_free_four_field(self):
  c=build(1,5,{**PHYS,'lambda0':0.},EV)
  g=measure(c,.12,None,hs=[.00625],verbose=False)
  self.assertLess(max(abs(complex(*r['F4c'])) for r in g['central_observables'].values()),2e-16)
  self.assertLess(max(abs(complex(*r['analytic_residual'])) for r in g['central_observables'].values()),2e-16)
if __name__=='__main__': unittest.main()
