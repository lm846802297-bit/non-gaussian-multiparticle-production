"""Finite-mode CTP vacuum connected four-point Born and one-loop terms.

Returns coefficients c_1,c_2 of lambda0*c_1+lambda0**2*c_2;
fixed-order diagrams, NOT a self-consistent 2PI/4PI solver.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import time

import numpy as np
from scipy.special import roots_legendre

PAIRINGS = (((0,1),(2,3)), ((0,2),(1,3)), ((0,3),(1,2)))
BRANCHES = np.array(list(itertools.product([0,1], repeat=4)))
LABELS = [''.join('+' if a==0 else '-' for a in row) for row in BRANCHES]


def ctp_mode(t, s, a, b, omega):
    """<T_C phi_p(t)_a phi_{-p}(s)_b> with mode normalization 1/(2 omega)."""
    dt=np.asarray(t)-np.asarray(s)
    if a==b:
        d=np.abs(dt)*(1 if a==0 else -1)
    else:
        d=dt*(1 if a==1 else -1)
    return np.exp(-1j*d[...,None]*omega)/(2*omega)


def compute(K,L,m,T,tmin,events,order):
    n=np.arange(-K,K+1)
    ps=2*np.pi*n/L
    omega=np.hypot(ps,m)
    Qs=np.arange(-2*K,2*K+1)
    pairs=[(np.nonzero((n[0]<=Q-n)&(Q-n<=n[-1]))[0],Q) for Q in Qs]

    def pair_reduce(z,w):
        # z and w have common trailing dimensions [sample, mode].
        out=[]
        for p,Q in pairs:
            q=Q-n[p]+K
            out.append(np.sum(z[...,p]*w[...,q],axis=-1))
        return np.stack(out,axis=-1)

    tx=np.array(events)
    phases=np.exp(1j*tx[:,1,None]*ps)

    def external(t,vertex):
        # shape [event, external branch, sample, mode].
        return np.stack([
            np.stack([ctp_mode(ti,t,a,vertex,omega)*phases[i] for a in (0,1)])
            for i,ti in enumerate(tx[:,0])
        ])

    def paired(ext,i,j):
        return pair_reduce(ext[i,:,None],ext[j,None,:]).reshape((4,)+ext.shape[2:-1]+(len(Qs),))

    if np.any(tx[:,0]<tmin): raise ValueError('called before initial time')
    edges=sorted(set([tmin]+tx[:,0].tolist()))
    nodes,w=roots_legendre(order)
    u,w=(nodes+1)/2,w/2
    intervals=[(lo,hi-lo,lo+(hi-lo)*u,(hi-lo)*w) for lo,hi in zip(edges,edges[1:])]
    born=np.zeros(16,complex)
    bubble=np.zeros(16,complex)
    for _,_,t,w1 in intervals:
        for a in (0,1):
            e=external(t,a)
            i,j,k,l=[e[v,BRANCHES[:,v]] for v in range(4)]
            ij=pair_reduce(i,j)
            kl=pair_reduce(k,l)
            born+=(-1j*(-1)**a/L**3)*np.einsum('ctq,ctq,t->c',ij,kl[...,::-1],w1/np.cosh(t/T)**2)
    for i_int,(lo,h,t,w1) in enumerate(intervals):
        for j_int,(_,_,s,w2) in enumerate(intervals):
            if i_int!=j_int:
                rules=[(np.repeat(t,order),np.tile(s,order),np.outer(w1,w2).ravel())]
            else:
                # Split at the |t-s| cusp using a Duffy map and its reflection.
                tt=np.repeat(t,order)
                ss=(lo+h*u[:,None]*u[None,:]).ravel()
                ww=(h*h*u[:,None]*np.outer(w,w)).ravel()
                rules=[(tt,ss,ww),(ss,tt,ww)]
            for tvec,svec,ww in rules:
                ww=ww/(np.cosh(tvec/T)*np.cosh(svec/T))**2
                es=[external(tvec,a) for a in (0,1)]
                fs=[external(svec,b) for b in (0,1)]
                for a,b in itertools.product((0,1),repeat=2):
                    f=ctp_mode(tvec,svec,a,b,omega)
                    internal=pair_reduce(f,f)*ww[:,None]
                    for (i,j),(k,l) in PAIRINGS:
                        left=paired(es[a],i,j)
                        # Momentum conservation routes Q from external left
                        # into the internal pair, -Q into external right.
                        right=paired(fs[b],k,l)[...,::-1]
                        mtx=np.einsum('atq,btq,tq->ab',left,right,internal,optimize=False)
                        bubble+=(-.5*(-1)**(a+b)/L**4)*mtx[2*BRANCHES[:,i]+BRANCHES[:,j],2*BRANCHES[:,k]+BRANCHES[:,l]]
    return {name:(c1,c2) for name,c1,c2 in zip(LABELS,born,bubble)}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--K',type=int,default=5)
    parser.add_argument('--L',type=float,default=20.)
    parser.add_argument('--m',type=float,default=1.)
    parser.add_argument('--T',type=float,default=1.)
    parser.add_argument('--tmin',type=float,default=-12.)
    parser.add_argument('--lam',type=float,default=.1)
    parser.add_argument('--orders',type=int,nargs='+',default=[30,40,52])
    parser.add_argument('--fock',type=Path)
    parser.add_argument('--events',type=json.loads)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    ref=None
    if args.fock:
        ref=json.loads(args.fock.read_text())
        assert ref['K']==args.K and ref['L']==args.L and ref['m']==args.m and ref['T']==args.T and ref['lambda0']==args.lam
        ref_val={name:complex(*array['F4c']) for name,array in ref['rows'][0]['components'].items()}
    events=args.events or (ref['events'] if ref is not None else [[1.,0.],[.5,1.],[-.5,2.],[-1.,3.]])
    if ref is not None: assert ref['events']==events
    start=time.monotonic()
    runs=[]
    z=lambda v:[v.real,v.imag]
    for order in args.orders:
        coeff=compute(args.K,args.L,args.m,args.T,args.tmin,events,order)
        res={}
        for name,(c1,c2) in coeff.items():
            born=args.lam*c1
            two=born+args.lam**2*c2
            row={'c1':z(c1),'c2':z(c2),'O_lambda':z(born),'through_O_lambda2':z(two)}
            if ref is not None:
                v=ref_val[name]
                row.update(Fock=z(v),relative_tree_error=abs(born-v)/abs(v),relative_through_lambda2_error=abs(two-v)/abs(v))
            res[name]=row
        runs.append({'order':order,'components':res})
    out=dict(K=args.K,L=args.L,m=args.m,T=args.T,tmin=args.tmin,lambda0=args.lam,events=events,runs=runs,walltime_seconds=time.monotonic()-start,
        implementation='Strictly fixed-order contour Born and one-loop bubbles in all three pair partitions; NOT a 2PI/4PI resummation.',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        fock_input_sha256=hashlib.sha256(args.fock.read_bytes()).hexdigest() if args.fock else None,
        python=platform.python_version(),numpy=np.__version__,scipy=__import__('scipy').__version__)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
