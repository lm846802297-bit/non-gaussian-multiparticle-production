import unittest
from itertools import product
import numpy as np
from scipy.sparse import csr_matrix
from scipy.integrate import solve_ivp
from cf4 import Propagator, setup_ops, evolve_all_forward, evaluate
class Check(unittest.TestCase):
    def data(self):
        return np.array([.2,2.7]),csr_matrix([[.7,-.84],[-.84,.3]]),np.array([.8,.6j])
    def test_noncommuting_order_four(self):
        E,V,v=self.data()
        l=lambda t: .4+np.sin(1.5*t)
        ref=solve_ivp(lambda t,w:-1j*(E*w+l(t)*(V@w)),(-1.,1.4),v,method='DOP853',rtol=3e-13,atol=3e-15).y[:,-1]
        P=Propagator(E,V,l)
        v1=P(v,-1.,1.4,.12); v2=P(v,-1.,1.4,.06)
        e1=np.linalg.norm(v1-ref);e2=np.linalg.norm(v2-ref)
        self.assertGreater(e1/e2,15.5)
        self.assertLess(e1/e2,16.5)
        self.assertLess(e2,5e-7)
        self.assertLess(abs(np.linalg.norm(v2)-1),2e-14)
    def test_constant_energy_exact_and_reversal(self):
        E,V,v=self.data()
        P=Propagator(E,V,lambda t:0.)
        ans=P(v,-1.,1.4,.15)
        np.testing.assert_allclose(ans,v*np.exp(-2.4j*E),rtol=0.,atol=3e-15)
        np.testing.assert_allclose(P(ans,1.4,-1.,.15),v,rtol=0.,atol=3e-15)
    def test_free_all_forward_matches_wick(self):
        K,N,m,L=1,5,1.,20.
        V,E,X,C,DX,_=setup_ops(K,N,m,L)
        prop=Propagator(E,V,lambda t:0.)
        states=evolve_all_forward(prop,X,.07)
        response=evaluate(prop,states,E,V,X,C,DX,0.,.07,.003125)
        q=2*np.pi/L*np.arange(-K,K+1)
        w=np.sqrt(m*m+q*q)
        data=[np.exp(-1j*w*(1.-t)-1j*q*x)/(2*L*w)
              for t,x in [(.5,1.),(-.5,2.),(-1.,3.)]]
        exact=6*sum(data[0][i]*data[1][j]*data[2][k]
                    for i,j,k in product(range(2*K+1),repeat=3)
                    if abs(i+j+k-3*K)<=K)
        self.assertLess(abs(complex(*response['F4c'])),7e-16)
        self.assertLess(abs(complex(*response['connected_composite_source'])-exact),7e-16)
        self.assertLess(abs(complex(*response['analytic_residual'])),7e-16)
if __name__=='__main__': unittest.main()
