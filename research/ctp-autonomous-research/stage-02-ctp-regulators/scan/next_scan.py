"""Supplemental Phase VII baseline derivative audit; does not require Phase VIII source."""
import argparse, json, platform, time, hashlib
from pathlib import Path
import numpy as np
import scipy
from sd_residual import Correlator, make_hamiltonian, evolve, pulse
from ctp_sd_all_branches import ctp_f4c, ctp_source, derivative_dxx_f4c

BASE="eab1f5b3198d03eafc3f162299c9038956beaa0b"
EV=[(1.,0.),(.5,1.),(-.5,2.),(-1.,3.)]
def z(q):return [q.real,q.imag]
def label(m):return ''.join('+' if m>>i&1 else '-' for i in range(4))

def run(K,N,hs,ms,rt,at,lam=.1,events=EV,length=20.):
 EV=events
 tic=time.monotonic()
 # Same Hamiltonian/operators as the pin, without its redundant t=+12 evolution.
 c=Correlator.__new__(Correlator)
 c.modes,c.basis,c.index,c.V,c.E,c.ws,c.herm=make_hamiltonian(K,N,1.,length)
 c.mass=1.; c.L=length; c.kcut=K
 c.phis={}; c.cubics={};c.dxxphis={}
 def prop(v,a,b):return evolve(v,a,b,c.E,c.V,lam,1.,rt,at,ms)
 print(f'built H: {len(c.basis)} {c.V.nnz} {time.monotonic()-tic:.2f}s',flush=True)
 vacuum=np.zeros(len(c.basis),complex); vacuum[c.index[(0,)*len(c.modes)]]=1
 after={0:(vacuum,-12.)}
 for m in range(2,16,2):
  j=min(i for i in range(1,4) if m>>i&1)
  t,x=EV[j]; v,s=after[m^(1<<j)]
  after[m]=(c.phi(t,x)@prop(v,s,t),t)
 # For any given evaluation time, the overlap after the last insertion is invariant to further propagation.
 def base(t):return {m:prop(v,s,t) for m,(v,s) in after.items()}
 def insert(even,apply):return even|{m|1:apply(v) for m,v in even.items()}
 t0,x0=EV[0]; phi=c.phi(t0,x0); dxx=c.dxxphi(t0,x0)
 bb=base(t0); s0=insert(bb, lambda v:phi@v)
 C=c.cubic(t0,x0); sc=insert(bb,lambda v:C@v); sx=insert(bb,lambda v:dxx@v)
 lt=pulse(t0,lam,1.)
 def H(v):return c.E*v+lt*(c.V@v)
 def comm(v):return H(phi@v)-phi@H(v)
 def d2(v):
  return -H(comm(v))+comm(H(v))-2j*np.tanh(t0)*lt*(c.V@(phi@v)-phi@(c.V@v))
 s2=insert(bb,d2)
 cc={m:ctp_f4c(s0,m)[0] for m in range(16)}
 src={m:ctp_source(s0,sc,m) for m in range(16)}
 dx={m:derivative_dxx_f4c(sx,s0,m) for m in range(16)}
 an={m:derivative_dxx_f4c(s2,s0,m)-dx[m]+cc[m]+lt*src[m][0]/6 for m in range(16)}
 rows=[]
 for h in hs:
  sp,sm=[insert({m:prop(v,t0,t0+h*sign) for m,v in bb.items()},lambda v:phi@v) for sign in [1,-1]]
  out={}
  for m in range(16):
   dtt=(ctp_f4c(sp,m)[0]-2*cc[m]+ctp_f4c(sm,m)[0])/h**2
   b=dtt-dx[m]+cc[m]; source=lt*src[m][0]/6; r=b+source
   out[label(m)]={"F4c":z(cc[m]),"M33":z(src[m][1]),"connected_composite_source":z(src[m][0]),"SD_residual":z(r),"SD_relative_residual":abs(r)/(abs(b)+abs(source)+1e-30),"analytic_residual":z(an[m]),"analytic_relative_residual":abs(an[m])/(abs(an[m]-source)+abs(source)+1e-30),"fd_error_vs_analytic":abs(r-an[m])}
  rows.append({"h":h,"components":out,"max_relative":max(a["SD_relative_residual"] for a in out.values()),"max_fd_error":max(a["fd_error_vs_analytic"] for a in out.values())})
 return {"source_commit":BASE,"runner_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"basis_dimension":len(c.basis),"V_nnz":c.V.nnz,"hermiticity_error":c.herm,"K":K,"Nmax":N,"m":1.,"L":length,"T":1.,"lambda0":lam,"events":EV,"rtol":rt,"atol":at,"max_step":ms,"max_analytic_abs":max(map(abs,an.values())),"rows":rows,"runtime":{"python":platform.python_version(),"scipy":scipy.__version__,"numpy":np.__version__,"platform":platform.platform(),"wall_seconds":time.monotonic()-tic}}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('K',type=int);p.add_argument('N',type=int);p.add_argument('out',type=Path)
 p.add_argument('--rtol',type=float,default=2e-10);p.add_argument('--atol',type=float,default=2e-12);p.add_argument('--step',type=float,default=.1)
 p.add_argument('--config', choices=['A','C'], default='A');p.add_argument('--length',type=float,default=20.)
 a=p.parse_args()
 result=run(a.K,a.N,[.00625,.003125,.0015625],a.step,a.rtol,a.atol,events=EV if a.config=='A' else [(1.2,0.),(.4,.7),(-.3,1.7),(-1.1,3.2)],length=a.length)
 a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 print(result["max_analytic_abs"],[(r["h"],r["max_relative"],r["max_fd_error"]) for r in result["rows"]],flush=True)
