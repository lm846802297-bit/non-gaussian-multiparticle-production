"""Three-loop 2PI with consistent bilocal Bethe–Salpeter response at K=0.

The canonical zero-mode coordinate q is sqrt(L)*phi. G=<T_C q q>.
Trapezoid branch weights include eta. The variational G, Hartree vacuum
subtraction, and derivative of the sunset kernel are all included.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import time

import numpy as np
from born_bubble import ctp_mode
EV=(1.,.5,-.5,-1.)

def solve(lam,h,L=20.,m=1.,T=1.,tmin=-12.,tmax=1.,tol=2e-13,maxiter=100):
    n=round((tmax-tmin)/h)
    assert abs(n*h-tmax+tmin)<1e-10
    tt=np.linspace(tmin,tmax,n+1)
    idx=np.array([round((t-tmin)/h) for t in EV])
    assert np.max(np.abs(tt[idx]-EV))<1e-10
    w=np.full(n+1,h); w[[0,-1]]*=.5
    eta_g=np.concatenate((w,-w))*np.tile(lam/(L*np.cosh(tt/T)**2),2)
    eta_gg=np.outer(eta_g,eta_g)
    g0=np.block([[ctp_mode(tt[:,None],tt[None,:],a,b,np.array([m]))[...,0] for b in (0,1)] for a in (0,1)])
    def local(g):
        return -.5j*eta_g*(np.diag(g)-.5/m)
    def nonlocal_kernel(g):
        return -(1/6)*eta_gg*(g**3)
    g=g0.copy()
    for it in range(maxiter):
        d=g0@(local(g)[:,None]*g+nonlocal_kernel(g)@g)
        gnew=g0+d
        residual=float(np.max(np.abs(gnew-g)))
        g=gnew
        if residual<tol:break
    else:raise RuntimeError(f'2PI Dyson residual {residual:g} not converged.')
    def kernel_derivative(x):
        return -.5j*eta_g*np.diag(x), -.5*eta_gg*g*g*x
    dres={}
    c4={}
    for b2,b3 in itertools.product((0,1),repeat=2):
        k,l=idx[2]+b2*(n+1),idx[3]+b3*(n+1)
        x0=np.outer(g[:,k],g[l,:])+np.outer(g[:,l],g[k,:])
        c=np.zeros_like(g)
        for nt in range(maxiter):
            kd,knd=kernel_derivative(x0+c)
            nxt=g@(kd[:,None]*g+knd@g)
            check=float(np.max(np.abs(nxt-c)))
            c=nxt
            if check<tol:break
        else:raise RuntimeError(f'BS residual {check:g} not converged.')
        dres[('+' if not b2 else '-')+('+' if not b3 else '-')]=[nt+1,check]
        for b0,b1 in itertools.product((0,1),repeat=2):
            name=''.join('+' if x==0 else '-' for x in (b0,b1,b2,b3))
            c4[name]=c[idx[0]+b0*(n+1),idx[1]+b1*(n+1)]/L**2
    return dict(G0_delta=float(np.max(np.abs(g-g0))),g_conj=float(np.max(np.abs(g-np.conj(np.roll(np.roll(g,n+1,axis=0),n+1,axis=1))))),
          F4c=c4,G_last=float(abs(g[n,n]-g[-1,-1])),dyson=[it+1,residual],bs=dres,nt_phys=n+1)

def cli():
    p=argparse.ArgumentParser();p.add_argument('--h',type=float,nargs='+',default=[.25,.125,.0625])
    p.add_argument('--couplings',type=float,nargs='+',default=[.1,.5,2.])
    p.add_argument('--include-reference-dir',type=Path)
    p.add_argument('out',type=Path); a=p.parse_args()
    rr={}
    start=time.monotonic()
    for l in a.couplings:
        key=f'{l:g}';rr[key]=[]
        refs=None
        if a.include_reference_dir is not None:
            refs=[json.loads((a.include_reference_dir/f'k0-exact-l{key}-n{n}.json').read_text()) for n in (12,16)]
            refa,refb=[{b:complex(*row['F4c']) for b,row in r['rows'][0]['components'].items()} for r in refs]
        for h in a.h:
            r=solve(l,h); r['h']=h
            if refs is not None:
                r['max_exact_fock_convergence_rel']=max(abs(refa[b]-refb[b])/abs(refb[b]) for b in refb)
                r['exact']={b:[z.real,z.imag] for b,z in refb.items()}
                r['rel_diff_exact']={b:abs(v-refb[b])/abs(refb[b]) for b,v in r['F4c'].items()}
            r['F4c']={k:[z.real,z.imag] for k,z in r['F4c'].items()}
            rr[key].append(r)
    out={'description':__doc__,'K':0,'L':20.,'m':1.,'T':1.,'tmin':-12.,'events':EV,
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'runtime_seconds':time.monotonic()-start,
         'python':platform.python_version(),'numpy':np.__version__,'runs':rr}
    a.out.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':cli()
