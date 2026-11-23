import unittest
import numpy as np
from ctp_sd_all_branches import evaluate_all
from direct_moments import run

class Direct(unittest.TestCase):
 def test_against_shared_archive_at_small_K(self):
  a=run(1,5,rtol=2e-10,atol=2e-12,max_step=.08)
  expected=evaluate_all(kcut=1,nmax=5,h=.003125,rtol=2e-10,atol=2e-12,max_step=.08)['components']
  for label,c in a['components'].items():
   f=expected[label]
   for name in ('F4c','M33','connected_composite_source','M31_F2_terms'):
    np.testing.assert_allclose(c[name],f[name],rtol=5e-8,atol=5e-11)
 def test_free_gaussian(self):
  a=run(1,5,lam=0,rtol=2e-10,atol=2e-12,max_step=.08)
  for c in a['components'].values():
   self.assertLess(abs(complex(*c['F4c'])),2e-10)
   self.assertLess(max(abs(complex(*v)) for v in c['M31'].values()),2e-10)
 def test_free_projected_wick_oracle(self):
  a=run(1,5,lam=0,rtol=2e-10,atol=2e-12,max_step=.08)
  k=np.arange(-1,2);p=2*np.pi*k/20.;w=np.hypot(p,1);times=[t for t,_ in a['events']];x=[x for _,x in a['events']]
  keep=abs(k[:,None,None]+k[None,:,None]+k[None,None,:])<=1
  def g(i,j,tag):
   d=times[i]-times[j];left,right=tag[i],tag[j]
   tau=(abs(d) if left=='+' else -abs(d)) if left==right else (-d if left=='+' else d)
   return np.exp(-1j*w*tau+1j*p*(x[i]-x[j]))/(40*w)
  for tag,c in a['components'].items():
   a1,a2,a3=[g(0,j,tag) for j in (1,2,3)]
   analytic_M33=6*np.sum((a1[:,None,None]*a2[None,:,None]*a3[None,None,:])[keep])
   np.testing.assert_allclose(complex(*c['M33']),analytic_M33,rtol=1e-8,atol=1e-10)
   for ij,v in c['F2'].items():
    np.testing.assert_allclose(complex(*v),g(int(ij[0]),int(ij[1]),tag).sum(),rtol=1e-8,atol=1e-10)
