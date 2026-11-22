"""Independent all-forward CF4 check; exact model operators from pinned baseline, source via Wick identity."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from sd_residual import make_hamiltonian, field_matrix, Correlator, pulse
from independent_source import setup, wick_matrix

class Propagator:
    """Two Gauss nodes/fourth-order commutator-free Magnus, manual Taylor action."""
    def __init__(self,E,V,lam,tol=2.e-17):
        self.center=.5*(max(E)+min(E))
        self.E=np.asarray(E)-self.center
        self.V=V.astype(complex)
        self.Enorm=float(np.max(np.abs(self.E)))
        self.Vnorm=float(np.max(np.asarray(abs(V).sum(axis=1))))
        self.lam=lam
        self.tol=tol
        self.count=0
        self.max_order=0
        self.max_radius=0.
    def exponential(self,w,dt,q):
        shp=(-1,)+(1,)*(w.ndim-1)
        if q==0.:
            return np.exp(-.5j*dt*(self.E+self.center)).reshape(shp)*w
        radius=abs(dt)*(.5*self.Enorm+abs(q)*self.Vnorm)
        self.max_radius=max(self.max_radius,radius)
        order=0
        remainder=radius
        while remainder>self.tol:
            order+=1
            remainder*=radius/(order+1)
        self.count+=1
        self.max_order=max(self.max_order,order)
        diagonal=(-.5j*dt*self.E).reshape(shp)
        coeff=-1j*dt*q
        term=w
        v=w.copy()
        for n in range(1,order+1):
            term=(diagonal*term+coeff*(self.V@term))/n
            v+=term
        return np.exp(-.5j*dt*self.center)*v
    def __call__(self,w,t0,t1,hmax):
        if t0==t1: return w.copy()
        if hmax<=0: raise ValueError("step must be positive")
        n=math.ceil(abs(t1-t0)/hmax-1.e-13)
        n=max(n,1)
        dt=(t1-t0)/n
        g=math.sqrt(3)/6
        a=(3-2*math.sqrt(3))/12
        b=(3+2*math.sqrt(3))/12
        for i in range(n):
            t=t0+(i+.5)*dt
            l,r=self.lam(t-g*dt),self.lam(t+g*dt)
            w=self.exponential(w,dt,b*l+a*r)
            w=self.exponential(w,dt,a*l+b*r)
        return w

def zz(z):
    return [float(z.real),float(z.imag)]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def setup_ops(K,N,m,L):
    modes,basis,ind,V,E,ws,herm=make_hamiltonian(K,N,m,L)
    c=Correlator.__new__(Correlator)
    c.modes,c.basis,c.index,c.mass,c.L,c.dxxphis=modes,basis,ind,m,L,{}
    # Matrix via independent sparse Wick construction, not the baseline monomial routine.
    b,p,occ,mm=setup(K,N)
    C=wick_matrix(b,p,occ,mm,K,N,m,L,0.)
    X=[field_matrix(modes,basis,ind,m,L,x) for x in range(4)]
    return V,E,X,C,c.dxxphi(1.,0.),herm

def evolve_all_forward(P,X,hstep):
    """Create all early-time ordered ket products once; evaluate with the vacuum at t1."""
    vacuum=np.zeros(len(P.E),complex)
    vacuum[0]=1.
    cols=P(vacuum,-12.,-1.,hstep)[:,None]
    prev=-1.
    for k,t in [(3,-1.),(2,-.5),(1,.5)]:
        cols=P(cols,prev,t,hstep)
        cols=np.hstack((cols,X[k]@cols))
        prev=t
    return P(cols,prev,1.,hstep)

def evaluate(P,cols,E,V,X,C,DX,la,hstep,h):
    # Columns correspond to subsets of the early fields ordered by binary labels
    # x=3, 2, 1. Multipliers enumerate the exact three pair subtractions.
    p=((1,6),(2,5),(4,3))
    def connected(a,endop):
        o=endop(a[:,[7,1,2,4]])
        bra=a[:,0].conj()
        raw,m1,m2,m4=bra@o
        terms=[z*np.vdot(a[:,0],a[:,j]) for (_,j),z in zip(p,[m1,m2,m4])]
        return raw-sum(terms),raw,terms
    f,raw4,pair=connected(cols,lambda v:X[0]@v)
    s,raw6,m31=connected(cols,lambda v:C@v)
    dxx,_,_=connected(cols,lambda v:DX@v)
    def H(w):
        return E[:,None]*w+la*(V@w)
    def comm(w):
        return H(X[0]@w)-X[0]@H(w)
    def dtime(w):
        vd=V@(X[0]@w)-X[0]@(V@w)
        return -H(comm(w))+comm(H(w))-2j*np.tanh(1.)*la*vd
    tt,_,_=connected(cols,dtime)
    fminus,*_=connected(P(cols,1.,1.-h,hstep),lambda v:X[0]@v)
    fplus,*_=connected(P(cols,1.,1.+h,hstep),lambda v:X[0]@v)
    fd=((fplus-f)-(f-fminus))/(h*h)
    B=tt-dxx+f
    B2=fd-dxx+f
    R=B+la*s/6
    R2=B2+la*s/6
    return {'F4':zz(raw4),'F4c':zz(f),'F2_F2':list(map(zz,pair)),'M33':zz(raw6),
            'M31_F2':list(map(zz,m31)),'connected_composite_source':zz(s),
            'analytic_dtt':zz(tt),'dxx':zz(dxx),'analytic_residual':zz(R),
            'analytic_relative_residual':float(abs(R)/(abs(B)+abs(la*s/6)+1e-30)),
            'h':h,'SD_residual':zz(R2),'SD_relative_residual':float(abs(R2)/(abs(B2)+abs(la*s/6)+1e-30)),
            'vacuum_norm':float(np.vdot(cols[:,0],cols[:,0]).real)}

def compare(a,b):
    out={}
    for key in ('F4c','M33','connected_composite_source','analytic_residual','SD_residual'):
        x,y=complex(*a[key]),complex(*b[key])
        out[key]={'absolute':abs(x-y),'relative_to_reference':abs(x-y)/max(abs(y),1e-30)}
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--K',type=int,default=5); ap.add_argument('--N',type=int,default=7)
    ap.add_argument('--steps',type=float,nargs='+',default=[.1,.05])
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--baseline',type=Path)
    args=ap.parse_args()
    tic=time.monotonic()
    V,E,X,C,DX,herm=setup_ops(args.K,args.N,1.,20.)
    print(f'constructed dimension {len(E)} in {time.monotonic()-tic:.2f}s',flush=True)
    out={'K':args.K,'Nmax':args.N,'basis_dimension':len(E),'hamiltonian_nnz':V.nnz,'hermiticity_error':herm,
         'events':[[1.,0.],[.5,1.],[-.5,2.],[-1.,3.]],'physical':{'m':1.,'L':20.,'lambda0':.1,'T':1.,'t_initial':-12.},
         'source_commit':'eab1f5b3198d03eafc3f162299c9038956beaa0b',
         'runner_sha256':sha(__file__),'source_module_sha256':sha(__import__('sd_residual').__file__),
         'wick_module_sha256':sha(__import__('independent_source').__file__),
         'machine':platform.machine(),'platform':platform.platform(),'python':platform.python_version(),
         'numpy':np.__version__,'scipy':scipy.__version__,'blas_threads':os.getenv('OPENBLAS_NUM_THREADS'),
         'steps':[]}
    ref=None
    if args.baseline:
        out['baseline_sha256']=sha(args.baseline)
        base=json.loads(args.baseline.read_text())
        assert (base['K'],base['Nmax'],base['events'])==(args.K,args.N,out['events'])
        assert all(out['physical'][k]==base[k] for k in ('m','L','lambda0','T'))
        ref=next(r for r in base['rows'] if r['h']==.003125)['components']['++++']
    for hstep in args.steps:
        step_tic=time.monotonic()
        P=Propagator(E,V,lambda t:pulse(t,.1,1.))
        cols=evolve_all_forward(P,X,hstep)
        data=evaluate(P,cols,E,V,X,C,DX,pulse(1.,.1,1.),hstep,.003125)
        if ref is not None: data['comparison_DOP853']=compare(data,ref)
        if out['steps']: data['comparison_previous_CF4']=compare(data,out['steps'][-1])
        data.update({'max_step':hstep,'max_taylor_order':P.max_order,'max_taylor_radius_bound':P.max_radius,
                     'taylor_action_tolerance':P.tol,'exp_actions':P.count,'wall_seconds':time.monotonic()-step_tic})
        out['steps'].append(data)
        out['elapsed_wall_seconds']=time.monotonic()-tic
        args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
        print(json.dumps(data,sort_keys=True),flush=True)

if __name__=='__main__':
    main()
