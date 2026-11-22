#!/usr/bin/env python3
"""Finite-volume test of the connected four-point Schwinger-Dyson equation.

1+1D normal-ordered real scalar lambda(t) phi^4, finite momentum projection.
The cubic source is projected onto the retained momentum modes, matching the
projected Hamiltonian exactly. This script computes selected all-forward
(time-ordered) moments and a finite-difference SD residual.
"""
from __future__ import annotations
import argparse, itertools, json, math, os, platform, sys, time
from collections import Counter
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix


def make_basis(kcut, nmax):
    modes=tuple(range(-kcut,kcut+1)); basis=[]
    def rec(i,left,occ):
        if i==len(modes): basis.append(tuple(occ)); return
        for n in range(left+1):
            occ.append(n); rec(i+1,left-n,occ); occ.pop()
    rec(0,nmax,[]); basis.sort(key=lambda s:(sum(s),s))
    return modes,basis,{s:i for i,s in enumerate(basis)}

def make_hamiltonian(kcut,nmax,mass,length):
    modes,basis,index=make_basis(kcut,nmax); nm=len(modes)
    ws=[math.sqrt((2*math.pi*p/length)**2+mass**2) for p in modes]
    u=[1/math.sqrt(2*w*length) for w in ws]
    mons=[]
    for nc in range(5):
        for ct in itertools.combinations_with_replacement(range(nm),nc):
            cc=Counter(ct); pc=sum(modes[i]*n for i,n in cc.items())
            fc=math.prod(u[i]**n/math.factorial(n) for i,n in cc.items())
            for at in itertools.combinations_with_replacement(range(nm),4-nc):
                ac=Counter(at); pa=sum(modes[i]*n for i,n in ac.items())
                if pa!=pc: continue
                fa=math.prod(u[i]**n/math.factorial(n) for i,n in ac.items())
                cv=tuple(cc.get(i,0) for i in range(nm)); av=tuple(ac.get(i,0) for i in range(nm))
                mons.append((cv,av,length*fc*fa))
    rows=[]; cols=[]; vals=[]
    for j,state in enumerate(basis):
        for cv,av,coef in mons:
            new=list(state); amp=coef; ok=True
            for i,n in enumerate(av):
                if new[i]<n: ok=False; break
                amp*=math.sqrt(math.factorial(new[i])/math.factorial(new[i]-n)); new[i]-=n
            if not ok: continue
            for i,n in enumerate(cv):
                amp*=math.sqrt(math.factorial(new[i]+n)/math.factorial(new[i])); new[i]+=n
            row=index.get(tuple(new))
            if row is not None and amp:
                rows.append(row); cols.append(j); vals.append(amp)
    V=coo_matrix((vals,(rows,cols)),shape=(len(basis),len(basis)),dtype=float).tocsr(); V.sum_duplicates()
    herm=float(np.max(np.abs((V-V.T).data))) if (V-V.T).nnz else 0.0
    energies=np.array([sum(n*w for n,w in zip(s,ws)) for s in basis])
    return modes,basis,index,V,energies,ws,herm

def pulse(t,lam,T):
    z=t/T
    if abs(z)>350: return 0.0
    c=math.cosh(z); return lam/(c*c)

def evolve(v,t0,t1,E,V,lam,T,rtol=2e-9,atol=2e-11,max_step=0.02):
    if abs(t1-t0)<1e-14: return v.copy()
    def rhs(t,y): return -1j*(E*y+pulse(t,lam,T)*(V@y))
    sol=solve_ivp(rhs,(t0,t1),v,method='DOP853',rtol=rtol,atol=atol,max_step=max_step)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.y[:,-1]

def field_matrix(modes,basis,index,mass,L,x):
    rows=[]; cols=[]; vals=[]
    for j,state in enumerate(basis):
        for i,nmode in enumerate(modes):
            p=2*math.pi*nmode/L; w=math.sqrt(p*p+mass*mass); u=1/math.sqrt(2*w*L)
            if state[i]:
                ns=list(state); ns[i]-=1; row=index.get(tuple(ns))
                if row is not None:
                    rows.append(row); cols.append(j); vals.append(u*np.exp(1j*p*x)*math.sqrt(state[i]))
            ns=list(state); ns[i]+=1; row=index.get(tuple(ns))
            if row is not None:
                rows.append(row); cols.append(j); vals.append(u*np.exp(-1j*p*x)*math.sqrt(state[i]+1))
    return coo_matrix((vals,(rows,cols)),shape=(len(basis),len(basis)),dtype=np.complex128).tocsr()

def projected_normal_ordered_cubic(modes,basis,index,mass,L,x):
    """P_K :phi_K(x)^3:_m in the same finite Fock representation.

    For three normal-ordered factors, group creation/annihilation multisets.
    Keep only monomials with |sum(p_ann)-sum(p_cre)| <= K, i.e. retained
    Fourier momentum of the composite. No tadpole contractions are included.
    """
    nm=len(modes); ws=[math.sqrt((2*math.pi*p/L)**2+mass**2) for p in modes]
    u=[1/math.sqrt(2*w*L) for w in ws]; K=max(abs(p) for p in modes)
    mons=[]
    for nc in range(4):
        na=3-nc
        for ct in itertools.combinations_with_replacement(range(nm),nc):
            cc=Counter(ct); pc=sum(modes[i]*n for i,n in cc.items())
            fc=math.prod(u[i]**n for i,n in cc.items())
            for at in itertools.combinations_with_replacement(range(nm),na):
                ac=Counter(at); pa=sum(modes[i]*n for i,n in ac.items())
                Q=pa-pc
                if abs(Q)>K: continue
                fa=math.prod(u[i]**n for i,n in ac.items())
                cv=tuple(cc.get(i,0) for i in range(nm)); av=tuple(ac.get(i,0) for i in range(nm))
                # Ordered triple expansion has 3!/prod(n!) multiplicity.
                mult=math.factorial(3)/(math.prod(math.factorial(n) for n in cc.values())*math.prod(math.factorial(n) for n in ac.values()))
                phase=np.exp(1j*(2*math.pi/L)*Q*x)
                mons.append((cv,av,mult*fc*fa*phase))
    rows=[]; cols=[]; vals=[]
    for j,state in enumerate(basis):
        for cv,av,coef in mons:
            new=list(state); amp=coef; ok=True
            for i,n in enumerate(av):
                if new[i]<n: ok=False; break
                amp*=math.sqrt(math.factorial(new[i])/math.factorial(new[i]-n)); new[i]-=n
            if not ok: continue
            for i,n in enumerate(cv):
                amp*=math.sqrt(math.factorial(new[i]+n)/math.factorial(new[i])); new[i]+=n
            row=index.get(tuple(new))
            if row is not None and amp:
                rows.append(row); cols.append(j); vals.append(amp)
    O=coo_matrix((vals,(rows,cols)),shape=(len(basis),len(basis)),dtype=np.complex128).tocsr(); O.sum_duplicates()
    return O

class Correlator:
    def __init__(self,kcut,nmax,mass,L,lam,T,rtol=2e-9,atol=2e-11,max_step=0.02):
        self.modes,self.basis,self.index,self.V,self.E,self.ws,self.herm=make_hamiltonian(kcut,nmax,mass,L)
        self.mass=mass; self.L=L; self.lam=lam; self.T=T; self.kcut=kcut
        self.t0=-12*T; self.tf=12*T; self.rtol=rtol; self.atol=atol; self.max_step=max_step
        self.psi0=np.zeros(len(self.basis),complex); self.psi0[self.index[tuple(0 for _ in self.modes)]]=1
        self.vac_out=evolve(self.psi0,self.t0,self.tf,self.E,self.V,lam,T,rtol,atol,max_step)
        self.phis={}; self.cubics={}; self.dxxphis={}; self.cache={}; self.prefix_cache={(): (self.psi0.copy(), self.t0)}
        self.delta=sum(1/(2*w*self.L) for w in self.ws)
    def phi(self,t,x):
        key=float(x)
        if key not in self.phis: self.phis[key]=field_matrix(self.modes,self.basis,self.index,self.mass,self.L,key)
        return self.phis[key]
    def cubic(self,t,x):
        key=float(x)
        if key not in self.cubics: self.cubics[key]=projected_normal_ordered_cubic(self.modes,self.basis,self.index,self.mass,self.L,key)
        return self.cubics[key]
    def dxxphi(self,t,x):
        key=float(x)
        if key not in self.dxxphis:
            rows=[]; cols=[]; vals=[]
            for j,state in enumerate(self.basis):
                for i,nmode in enumerate(self.modes):
                    p=2*math.pi*nmode/self.L; w=math.sqrt(p*p+self.mass*self.mass); u=1/math.sqrt(2*w*self.L)
                    coeff=-p*p*u
                    if state[i]:
                        ns=list(state); ns[i]-=1; row=self.index.get(tuple(ns))
                        if row is not None:
                            rows.append(row); cols.append(j); vals.append(coeff*np.exp(1j*p*key)*math.sqrt(state[i]))
                    ns=list(state); ns[i]+=1; row=self.index.get(tuple(ns))
                    if row is not None:
                        rows.append(row); cols.append(j); vals.append(coeff*np.exp(-1j*p*key)*math.sqrt(state[i]+1))
            self.dxxphis[key]=coo_matrix((vals,(rows,cols)),shape=(len(self.basis),len(self.basis)),dtype=np.complex128).tocsr()
        return self.dxxphis[key]
    def moment(self,events):
        # Prefix-state memoization: sequences sharing an identical early insertion
        # history reuse the already propagated ket. This changes no physics or
        # integration tolerances; it only avoids repeated DOP853 integrations.
        ev=tuple(sorted([(float(t),float(x),str(k)) for t,x,k in events],key=lambda z:z[0]))
        if ev in self.cache: return self.cache[ev]
        if ev in self.prefix_cache:
            state, now = self.prefix_cache[ev]
        else:
            prefix = ev
            while prefix not in self.prefix_cache:
                prefix = prefix[:-1]
            state, now = self.prefix_cache[prefix]
            for i in range(len(prefix), len(ev)):
                t,x,kind=ev[i]
                state=evolve(state,now,t,self.E,self.V,self.lam,self.T,self.rtol,self.atol,self.max_step)
                if kind=='phi': op=self.phi(t,x)
                elif kind=='cubic': op=self.cubic(t,x)
                elif kind=='dxx': op=self.dxxphi(t,x)
                else: raise ValueError(kind)
                state=op@state
                now=t
                self.prefix_cache[ev[:i+1]]=(state.copy(),now)
        final=evolve(state,now,self.tf,self.E,self.V,self.lam,self.T,self.rtol,self.atol,self.max_step)
        val=np.vdot(self.vac_out,final)
        self.cache[ev]=val
        return val
    def F2(self,a,b): return self.moment([(a[0],a[1],'phi'),(b[0],b[1],'phi')])
    def F4(self,events):
        e=[(t,x,'phi') for t,x in events]
        f4=self.moment(e)
        pairs=((0,1,2,3),(0,2,1,3),(0,3,1,2))
        sub=0j
        for i,j,k,l in pairs: sub+=self.F2(events[i],events[j])*self.F2(events[k],events[l])
        return f4-sub, f4
    def M31(self,a,b): return self.moment([(a[0],a[1],'cubic'),(b[0],b[1],'phi')])
    def F4_dxx_first(self,events):
        a,b,c,d=events
        raw=self.moment([(a[0],a[1],'dxx'),(b[0],b[1],'phi'),(c[0],c[1],'phi'),(d[0],d[1],'phi')])
        sub=(self.moment([(a[0],a[1],'dxx'),(b[0],b[1],'phi')])*self.F2(c,d)
             +self.moment([(a[0],a[1],'dxx'),(c[0],c[1],'phi')])*self.F2(b,d)
             +self.moment([(a[0],a[1],'dxx'),(d[0],d[1],'phi')])*self.F2(b,c))
        return raw-sub
    def M33(self,a,b,c,d): return self.moment([(a[0],a[1],'cubic'),(b[0],b[1],'phi'),(c[0],c[1],'phi'),(d[0],d[1],'phi')])
    def source(self,events):
        a,b,c,d=events
        m33=self.M33(a,b,c,d)
        terms=[self.M31(a,b)*self.F2(c,d),self.M31(a,c)*self.F2(b,d),self.M31(a,d)*self.F2(b,c)]
        return m33-sum(terms),m33,terms

def evaluate(kcut,nmax,mass,L,lam,T,h,a,events,rtol,atol,max_step):
    calc=Correlator(kcut,nmax,mass,L,lam,T,rtol,atol,max_step)
    e0=[tuple(x) for x in events]
    F0,_=calc.F4(e0)
    # Event 0 is the differentiated, latest-time coordinate.
    ep=[(e0[0][0]+h,e0[0][1]),*e0[1:]]
    em=[(e0[0][0]-h,e0[0][1]),*e0[1:]]
    xp=[(e0[0][0],e0[0][1]+a),*e0[1:]]
    xm=[(e0[0][0],e0[0][1]-a),*e0[1:]]
    Ftp,_=calc.F4(ep); Ftm,_=calc.F4(em); Fxp,_=calc.F4(xp); Fxm,_=calc.F4(xm)
    dtt=(Ftp-2*F0+Ftm)/(h*h)
    dxx_fd=(Fxp-2*F0+Fxm)/(a*a)
    dxx_spectral=calc.F4_dxx_first(e0)
    box=dtt-dxx_spectral
    src,m33,terms=calc.source(e0)
    lam1=pulse(e0[0][0],lam,T)
    eom=(box+mass*mass*F0)+(lam1/6)*src
    denom=abs(box+mass*mass*F0)+abs(lam1*src/6)+1e-30
    return {
      'parameters':{'Kcut':kcut,'Nmax':nmax,'m':mass,'L':L,'lambda0':lam,'T':T,'h':h,'a':a,'events':events,'rtol':rtol,'atol':atol,'max_step':max_step},
      'basis_dimension':len(calc.basis),'interaction_hermiticity_error':calc.herm,'vacuum_norm_out':float(np.vdot(calc.vac_out,calc.vac_out).real),'normal_order_delta':calc.delta,
      'F4c':[F0.real,F0.imag],'dtt_F4c':[dtt.real,dtt.imag],'dxx_F4c_fd':[dxx_fd.real,dxx_fd.imag], 'dxx_F4c_spectral':[dxx_spectral.real,dxx_spectral.imag],
      'box_plus_mass_F4c':[(box+mass*mass*F0).real,(box+mass*mass*F0).imag],
      'M33_projected':[m33.real,m33.imag], 'M31_F2_terms':[[z.real,z.imag] for z in terms], 'connected_source':[src.real,src.imag],
      'lambda_at_t1':lam1,'source_term':[ (lam1*src/6).real,(lam1*src/6).imag],
      'SD_residual':[eom.real,eom.imag],'SD_residual_abs':abs(eom),'SD_relative_residual':abs(eom)/denom,
      'moment_cache_entries':len(calc.cache),
      'runtime':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'machine':platform.machine(),'processor':platform.processor(),'cpu_count':os.cpu_count()}
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--kcut',type=int,default=1); ap.add_argument('--nmax',type=int,default=8)
    ap.add_argument('--mass',type=float,default=1.0); ap.add_argument('--length',type=float,default=20.0)
    ap.add_argument('--lambda0',type=float,default=0.2); ap.add_argument('--duration',type=float,default=1.0)
    ap.add_argument('--h',type=float,default=0.05); ap.add_argument('--a',type=float,default=0.2)
    ap.add_argument('--rtol',type=float,default=2e-9); ap.add_argument('--atol',type=float,default=2e-11); ap.add_argument('--max-step',type=float,default=0.02)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    events=[(1.0,0.0),(0.5,1.0),(-0.5,2.0),(-1.0,3.0)]
    result=evaluate(args.kcut,args.nmax,args.mass,args.length,args.lambda0,args.duration,args.h,args.a,events,args.rtol,args.atol,args.max_step)
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__': main()