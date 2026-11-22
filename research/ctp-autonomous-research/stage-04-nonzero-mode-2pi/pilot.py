"""Replayable K=1 three-loop real-time 2PI+BS/Fock comparison."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time
import numpy as np
import scipy
from born_bubble import compute
from explicit_wick import wick
from k1_bs import solve

def z(d):return {k:[float(v.real),float(v.imag)] for k,v in d.items()}
def cz(d):return {k:complex(*v) for k,v in d.items()}
def diff(a,b):return {k:abs(a[k]-v)/abs(v) for k,v in b.items()}
def rich(a,b,h1,h2):return {k:(h1*h1*b[k]-h2*h2*a[k])/(h1*h1-h2*h2) for k in a}
def checks(v,latest):
    tr=str.maketrans('+-','-+')
    return dict(max_complement_abs=max(abs(a-np.conj(v[k.translate(tr)])) for k,a in v.items()),
          max_largest_time_abs=max(abs(a-v[k[:latest]+k[latest].translate(tr)+k[latest+1:]]) for k,a in v.items()))
def cli():
    p=argparse.ArgumentParser();p.add_argument('reference',type=Path);p.add_argument('out',type=Path)
    p.add_argument('--h',nargs='+',type=float,default=[.25,.125,.0625]);p.add_argument('--couplings',nargs='+',type=float,default=[.1,.5,2.,5.])
    a=p.parse_args()
    K=1; L=20.; m=1.; T=1.; tmin=-4.
    raw=json.loads(a.reference.read_text()); ev=raw['events']; s=raw['runs']
    assert all(raw[k]==v for k,v in dict(K=K,L=L,m=m,T=T,tmin=tmin).items())
    assert all(x>y for x,y in zip(a.h,a.h[1:])) and len(a.h)>=2
    end=max(range(4),key=lambda i:ev[i][0])
    start=time.monotonic()
    quadrature=[compute(K,L,m,T,tmin,ev,z) for z in (40,48,56)]
    p=quadrature[-1]
    qdelta=max(abs(v-p[k][i]) for q in quadrature[:-1] for k,vs in q.items() for i,v in enumerate(vs))
    pw={h:wick(ev,h,K,L,m,T,tmin) for h in a.h}
    runs={}
    for l in a.couplings:
        f=s[f'{l:g}']; r=cz(f['N14_tighter']['F4c'])
        ref_errors={lab:max(diff(cz(v['F4c']),r).values()) for lab,v in f.items() if lab!='N14_tighter'}
        ep={k:l*c1+l*l*c2 for k,(c1,c2) in p.items()}
        born={k:l*c1 for k,(c1,c2) in p.items()}
        rows=[]
        for h in a.h:
            data=solve(l,h,ev,K,L,m,T,tmin)
            v=data.pop('F4c'); by=data.pop('F4c_by_Q')
            correction={k:ep[k]-(l*c1+l*l*c2) for k,(c1,c2) in pw[h].items()}
            vc={k:u+correction[k] for k,u in v.items()}
            data.update(F4c=z(v),F4c_by_Q={Q:z(q) for Q,q in by.items()},through_lambda2_correction=z(correction),
                        corrected=z(vc),raw_err=diff(v,r),corrected_err=diff(vc,r),**checks(v,end))
            rows.append(data)
        before,last=rows[-2:]; h1,h2=before['step'],last['step']
        plain=rich(cz(before['F4c']),cz(last['F4c']),h1,h2)
        converged=rich(cz(before['corrected']),cz(last['corrected']),h1,h2)
        summary={'raw_Richardson':z(plain),'corrected_Richardson':z(converged),'plain_errors':diff(plain,r),
                 'corrected_errors':diff(converged,r),
                 'corrected_pair_shift':diff(cz(before['corrected']),cz(last['corrected']))}
        runs[f'{l:g}']=dict(Fock=z(r),Fock_error_max=ref_errors,Born=z(born),through_lambda2=z(ep),
                 Born_err=diff(born,r),through_lambda2_err=diff(ep,r),grids=rows,richardson=summary)
        print(l,max(summary['corrected_errors'].values()),max(summary['corrected_pair_shift'].values()),ref_errors,flush=True)
    scripts=['pilot.py','k1_bs.py','explicit_wick.py','born_bubble.py']
    sig=lambda s:hashlib.sha256(Path(s).read_bytes()).hexdigest()
    out=dict(K=K,L=L,m=m,T=T,tmin=tmin,events=ev,program='bilocal BS from full three-loop 2PI kernel derivative',
             quadrature_orders=[40,48,56],quad_coeff_max_delta=qdelta,
             source_hashes={s:sig(s) for s in scripts},ref_hash=sig(a.reference),runs=runs,
             elapsed_seconds=time.monotonic()-start,numpy=np.__version__,scipy=scipy.__version__,python=platform.python_version())
    a.out.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':cli()
