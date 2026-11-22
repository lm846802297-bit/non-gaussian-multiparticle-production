#!/usr/bin/env python3
"""Selected unequal-time CTP/Feynman four-point function in a finite Fock space.

All four fields are on the forward branch and their times are supplied in
non-increasing operator order. The calculation evaluates the time-ordered
four-field moment and subtracts all three two-point pairings. It is a finite
mode/Fock-space benchmark, not a continuum result.
"""
from __future__ import annotations
import itertools, math, json, sys, os, platform, time
from collections import Counter
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix


def full_basis(kcut,nmax):
    modes=tuple(range(-kcut,kcut+1)); M=len(modes); basis=[]
    def rec(i,rem,occ):
        if i==M:
            basis.append(tuple(occ)); return
        for n in range(rem+1):
            occ.append(n); rec(i+1,rem-n,occ); occ.pop()
    rec(0,nmax,[])
    basis.sort(key=lambda s:(sum(s),s))
    return modes,basis,{s:i for i,s in enumerate(basis)}


def build_full_hamiltonian(kcut,nmax,mass,length):
    modes,basis,index=full_basis(kcut,nmax); M=len(modes)
    ws=[math.sqrt((2*math.pi*p/length)**2+mass**2) for p in modes]
    u=[1/math.sqrt(2*w*length) for w in ws]
    mon=[]
    for nc in range(5):
        for ct in itertools.combinations_with_replacement(range(M),nc):
            cc=Counter(ct); cp=sum(modes[i]*n for i,n in cc.items())
            cf=math.prod(u[i]**n/math.factorial(n) for i,n in cc.items())
            for at in itertools.combinations_with_replacement(range(M),4-nc):
                ac=Counter(at); ap=sum(modes[i]*n for i,n in ac.items())
                if ap!=cp: continue
                af=math.prod(u[i]**n/math.factorial(n) for i,n in ac.items())
                cv=tuple(cc.get(i,0) for i in range(M)); av=tuple(ac.get(i,0) for i in range(M))
                mon.append((cv,av,length*cf*af))
    rows=[]; cols=[]; vals=[]
    for j,state in enumerate(basis):
        for cv,av,coef in mon:
            new=list(state); amp=coef; ok=True
            for i,nann in enumerate(av):
                if new[i]<nann: ok=False; break
                amp*=math.sqrt(math.factorial(new[i])/math.factorial(new[i]-nann)); new[i]-=nann
            if not ok: continue
            for i,ncre in enumerate(cv):
                amp*=math.sqrt(math.factorial(new[i]+ncre)/math.factorial(new[i])); new[i]+=ncre
            row=index.get(tuple(new))
            if row is not None and amp:
                rows.append(row); cols.append(j); vals.append(amp)
    V=coo_matrix((vals,(rows,cols)),shape=(len(basis),len(basis)),dtype=float).tocsr(); V.sum_duplicates()
    anti=V-V.T; herm=float(np.max(np.abs(anti.data))) if anti.nnz else 0.0
    E=np.array([sum(n*w for n,w in zip(s,ws)) for s in basis])
    return modes,basis,index,V,E,ws,herm


def pulse(t,lam,T):
    z=t/T
    if abs(z)>350:return 0.0
    c=math.cosh(z); return lam/(c*c)


def propagate(v,t0,t1,E,V,lam,T,rtol=2e-9,atol=2e-11):
    if abs(t1-t0)<1e-15:return v.copy()
    def rhs(t,y): return -1j*(E*y+pulse(t,lam,T)*(V@y))
    max_step=min(T/24,0.02)
    sol=solve_ivp(rhs,(t0,t1),v,method='DOP853',rtol=rtol,atol=atol,max_step=max_step)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.y[:,-1]


def field_matrix(modes,basis,index,mass,length,x):
    rows=[];cols=[];vals=[]
    for j,s in enumerate(basis):
        for i,nmode in enumerate(modes):
            p=2*math.pi*nmode/length; w=math.sqrt(p*p+mass*mass); u=1/math.sqrt(2*w*length)
            if s[i]:
                ns=list(s);ns[i]-=1;r=index.get(tuple(ns))
                if r is not None: rows.append(r);cols.append(j);vals.append(u*np.exp(1j*p*x)*math.sqrt(s[i]))
            ns=list(s);ns[i]+=1;r=index.get(tuple(ns))
            if r is not None: rows.append(r);cols.append(j);vals.append(u*np.exp(-1j*p*x)*math.sqrt(s[i]+1))
    return coo_matrix((vals,(rows,cols)),shape=(len(basis),len(basis)),dtype=np.complex128).tocsr()


def ordered_moment(events,psi0,t0,operators,E,V,lam,T):
    """events [(t,x),...] in descending time and matching operator order."""
    times=[e[0] for e in events]
    if any(times[i]<times[i+1] for i in range(len(times)-1)):
        raise ValueError('events must be ordered with non-increasing times')
    v=propagate(psi0,t0,times[-1],E,V,lam,T)
    for j in range(len(events)-1,-1,-1):
        t,x=events[j]
        v=operators[x]@v
        if j>0:
            v=propagate(v,times[j],times[j-1],E,V,lam,T)
    bra=propagate(psi0,t0,times[0],E,V,lam,T)
    return np.vdot(bra,v)


def connected_four(events,psi0,t0,operators,E,V,lam,T):
    f4=ordered_moment(events,psi0,t0,operators,E,V,lam,T)
    pairs=[]
    for i,j,k,l in [(0,1,2,3),(0,2,1,3),(0,3,1,2)]:
        g1=ordered_moment([events[i],events[j]],psi0,t0,operators,E,V,lam,T)
        g2=ordered_moment([events[k],events[l]],psi0,t0,operators,E,V,lam,T)
        pairs.append(g1*g2)
    return {"F4": [float(f4.real),float(f4.imag)],
            "pair_products": [[float(x.real),float(x.imag)] for x in pairs],
            "F4_connected": [float((f4-sum(pairs)).real),float((f4-sum(pairs)).imag)]}


def calculate(nmax=16, lam=.2, T=1.0, kcut=1, mass=1.0, length=20.0):
    t0=-12*T
    modes,basis,index,V,E,ws,herm=build_full_hamiltonian(kcut,nmax,mass,length)
    psi0=np.zeros(len(basis),complex);psi0[index[tuple(0 for _ in modes)]]=1
    positions=[0.0,1.0,2.0,3.0]
    operators={x:field_matrix(modes,basis,index,mass,length,x) for x in positions}
    configs={
      "equal_time": [ (0.0,0.0),(0.0,0.0),(0.0,0.0),(0.0,2.0) ],
      "unequal_time_A": [ (1.0,0.0),(0.5,1.0),(-0.5,2.0),(-1.0,3.0) ],
      "unequal_time_B": [ (2.0,0.0),(1.0,1.0),(-1.0,2.0),(-2.0,3.0) ]
    }
    out={"model":"1+1D normal-ordered lambda(t) phi^4, finite full Fock basis",
         "parameters":{"Kcut":kcut,"Nmax":nmax,"m":mass,"L":length,"lambda0":lam,"T":T,"t0":t0},
         "basis_dimension":len(basis),"interaction_hermiticity_error":herm,"results":{}}
    for name,events in configs.items():
        out["results"][name]={"events":[[t,x] for t,x in events],**connected_four(events,psi0,t0,operators,E,V,lam,T)}
    out["runtime"]={"python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__,"platform":platform.platform(),"machine":platform.machine(),"cpu_count":os.cpu_count()}
    return out


def main():
    import argparse
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser()
    ap.add_argument('--nmax',type=int,default=16)
    ap.add_argument('--lambda0',type=float,default=.2)
    ap.add_argument('--duration',type=float,default=1.0)
    ap.add_argument('--output',type=Path,default=root/'results/phase2-dynamics-001/unequal-time-ctp-four-point.json')
    args=ap.parse_args()
    out=calculate(nmax=args.nmax,lam=args.lambda0,T=args.duration)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=='__main__': main()