"""Matching finite Fock-Hamiltonian reference for the K=1 closure pilot."""
import argparse
import itertools
import json
import time
import hashlib
from pathlib import Path
import numpy as np
import scipy
from sd_residual import make_hamiltonian, field_matrix, evolve

A=((1.,0.),(.5,1.),(-.5,2.),(-1.,3.))
B=((1.5,0.),(.75,1.5),(-.75,3.),(-1.5,4.5))
C=((1.2,0.),(.4,.7),(-.3,1.7),(-1.1,3.2))

def fock(lam,N,events=A,K=1,L=20.,m=1.,T=1.,tmin=-4.,rt=2e-12,at=2e-14,step=.08):
    ta=max(t for t,_ in events)
    pm,basis,index,mat,e,omegas,herm=make_hamiltonian(K,N,m,L)
    ins=[field_matrix(pm,basis,index,m,L,x) for _,x in events]
    ini=np.zeros(len(basis),complex);ini[0]=1
    def prop(x,a,b):return evolve(x,a,b,e,mat,lam,T,rt,at,step)
    rooted={0:(tmin,ini)}
    states={}
    for mask in range(16):
        if mask:
            j=max([q for q in range(4) if mask>>q&1],key=lambda q:events[q][0])
            ti,v=rooted[mask^(1<<j)]; te=events[j][0]
            rooted[mask]=(te,ins[j]@prop(v,ti,te))
        ti,v=rooted[mask]
        states[mask]=prop(v,ti,ta)
    def moment(b,s):
        return np.vdot(states[s^(s&b)],states[s&b])
    f4={}
    for mask in range(16):
        v=moment(mask,15)-sum(moment(mask,s)*moment(mask,s^15) for s in (3,5,9))
        label=''.join('+' if mask&(1<<i) else '-' for i in range(4))
        f4[label]=v
    return f4,len(basis),herm

def cli():
    pp=argparse.ArgumentParser(); pp.add_argument('out',type=Path); pp.add_argument('--case',choices='ABC',default='A')
    pp.add_argument('--couplings',nargs='+',type=float,default=[.1,.5,2.,5.]);args=pp.parse_args()
    ev={'A':A,'B':B,'C':C}[args.case]
    st=time.monotonic();runs={}
    for l in args.couplings:
        k=f'{l:g}'
        runs[k]={}
        for N,tol in ((8,'tight'),(10,'tight'),(12,'tight'),(14,'tight'),(14,'tighter')):
            p={} if tol=='tight' else {'rt':2e-13,'at':2e-15,'step':.04}
            f,d,herm=fock(l,N,ev,**p)
            runs[k][f'N{N}_{tol}']={'D':d,'H_herm':herm,'F4c':{b:[v.real,v.imag] for b,v in f.items()}}
        print(k, runs[k]['N14_tighter']['F4c']['++++'],flush=True)
    scripts=['reference.py','sd_residual.py']
    out={'K':1,'m':1.,'L':20.,'T':1.,'tmin':-4.,'events':ev,'runs':runs,'runtime_seconds':time.monotonic()-st,'numpy':np.__version__,'scipy':scipy.__version__,
         'Fock_builder_reference_commit':'eab1f5b3198d03eafc3f162299c9038956beaa0b',
         'sha256':{s:hashlib.sha256(Path(s).read_bytes()).hexdigest() for s in scripts}}
    args.out.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':cli()
