"""On same quadrature grid verify both resummations' actual lambda terms."""
import itertools
import unittest
import numpy as np
from scipy.special import roots_legendre
from two_pi_k0 import solve as two, EV
from four_pi_k0 import solve as four

def wick_direct(h=.5,tmin=-4.,L=20.):
    t=np.linspace(tmin,max(EV),round((max(EV)-tmin)/h)+1)
    xnodes=np.tile(t,2)
    sign=np.repeat([1,-1],len(t))
    w=np.full(len(t),h);w[[0,-1]]*=.5
    weight=sign*np.tile(w,2)/(L*np.cosh(xnodes)**2)
    def pair(s,u,a,b):
        tau=s-u
        if a==b:tau=np.abs(tau)*a
        elif a==1:tau=-tau
        return np.exp(-1j*tau)/2
    G=pair(xnodes[:,None],xnodes[None,:],1,1)
    for a,b in itertools.product((+1,-1),repeat=2):
        ai,bi=np.where(sign==a)[0],np.where(sign==b)[0]
        G[np.ix_(ai,bi)]=pair(xnodes[ai,None],xnodes[None,bi],a,b)
    out={}
    for branch in itertools.product((+1,-1),repeat=4):
        f=np.array([np.array([pair(ev,n,sg,s) for n,s in zip(xnodes,sign)]) for ev,sg in zip(EV,branch)])
        c1=-1j*np.dot(weight,np.prod(f,axis=0))/L**2
        c2=-sum((weight*f[0]*f[j])@(G**2)@(weight*f[k]*f[l]) for j,(k,l) in ((1,(2,3)),(2,(1,3)),(3,(1,2))))/(2*L**2)
        out[''.join('+' if a==1 else '-' for a in branch)]=(c1,c2)
    return out

class WeakOrder(unittest.TestCase):
    def test_explicit_wick_coefficients_both_solvers(self):
        p=wick_direct()
        lam=.01
        for solver in (two,four):
            a,b=[solver(s*lam,.5,tmin=-4.)['F4c'] for s in (+1,-1)]
            for m,(c1,c2) in p.items():
                np.testing.assert_allclose((a[m]-b[m])/(2*lam),c1,rtol=2e-7,atol=1e-11)
                np.testing.assert_allclose((a[m]+b[m])/(2*lam**2),c2,rtol=2e-4,atol=1e-11)
