import unittest
import numpy as np
from scipy.integrate import solve_ivp
from born_bubble import compute

class BubbleCheck(unittest.TestCase):
    def test_oscillator_fock_coefficients(self):
        t=.3; lam=.005
        vals=compute(0,20.,1.,1.,-12.,[(t,0.)]*4,48)['++++']
        D=16; a=np.diag(np.sqrt(np.arange(1,D)),1); ad=a.T
        u=1/np.sqrt(40.); phi=u*(a+ad)
        W=20*u**4/24*(ad@ad@ad@ad+4*ad@ad@ad@a+6*ad@ad@a@a+4*ad@a@a@a+a@a@a@a)
        p2=phi@phi; p4=p2@p2
        n=np.arange(D)
        ini=np.zeros(D,complex); ini[0]=1
        def cumulant(lam):
            sol=solve_ivp(lambda s,v:-1j*(n*v+lam/np.cosh(s)**2*(W@v)),(-12.,t),ini,method='DOP853',rtol=2e-13,atol=5e-16,max_step=.06)
            v=sol.y[:,-1]
            return (np.vdot(v,p4@v)-3*np.vdot(v,p2@v)**2).real
        pm,pp=cumulant(-lam),cumulant(lam)
        np.testing.assert_allclose((pp-pm)/(2*lam),vals[0],rtol=1e-5,atol=1e-10)
        np.testing.assert_allclose((pp+pm)/(2*lam**2),vals[1],rtol=5e-4,atol=5e-10)
    def test_ctp_identities(self):
        out=compute(1,20.,1.,1.,-12.,[(1.,0.),(.5,1.),(-.5,2.),(-1.,3.)],22)
        for k,v in out.items():
            f=k.translate(str.maketrans('+-','-+'))
            np.testing.assert_allclose(v,np.conj(out[f]),rtol=1e-10,atol=1e-12)
            np.testing.assert_allclose(v,out[f[0]+k[1:]],rtol=1e-10,atol=1e-12)
