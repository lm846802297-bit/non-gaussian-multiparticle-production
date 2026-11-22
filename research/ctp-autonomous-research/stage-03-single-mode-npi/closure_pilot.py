"""Matched K=0 scalar-pulse CTP reference: oscillator vs 2PI/BS and 4PI."""
import argparse
from itertools import combinations, product
import json
from pathlib import Path
import hashlib
import platform
import time
import numpy as np
from scipy.integrate import solve_ivp
from two_pi_k0 import solve as two, EV
from four_pi_k0 import solve as four
from born_bubble import compute

L=20.
TMIN=-4.

def fock(lam, D, rt=2e-12, at=2e-14, step=.08):
    a=np.diag(np.sqrt(np.arange(1,D)),1)
    ad=a.T
    q=(a+ad)/np.sqrt(2.)
    w=(ad@ad@ad@ad+4*ad@ad@ad@a+6*ad@ad@a@a+4*ad@a@a@a+a@a@a@a)/(96.*L)
    num=np.arange(D)
    ini=np.zeros(D,complex);ini[0]=1
    ts=sorted(EV)
    def fun(t,s):return -1j*(num*s+lam/np.cosh(t)**2*(w@s))
    rooted={0:(TMIN,ini)}
    states={}
    for sub in range(16):
        if sub:
            j=min(i for i in range(4) if sub&(1<<i))
            ti,z=rooted[sub^(1<<j)]; te=EV[j]
            r=solve_ivp(fun,(ti,te),z,method='DOP853',rtol=rt,atol=at,max_step=step)
            assert r.success
            rooted[sub]=(te,q@r.y[:,-1])
        ti,z=rooted[sub]
        if ti!=ts[-1]:
            r=solve_ivp(fun,(ti,ts[-1]),z,method='DOP853',rtol=rt,atol=at,max_step=step)
            assert r.success
            z=r.y[:,-1]
        states[sub]=z
    def moment(b,s):
        return np.vdot(states[s^(s&b)],states[s&b])
    f4={}
    for mask in range(16):
        v=moment(mask,15)-sum(moment(mask,s)*moment(mask,s^15) for s in (3,5,9))
        label=''.join('+' if mask&(1<<i) else '-' for i in range(4))
        f4[label]=v/L**2
    return f4

def main():
    p=argparse.ArgumentParser();p.add_argument('--h',type=float,nargs='+',default=[.5,.25,1/6])
    p.add_argument('--couplings',type=float,nargs='+',default=[.1,.5,2.])
    p.add_argument('out',type=Path);a=p.parse_args()
    rr={};start=time.monotonic()
    z=lambda c:[c.real,c.imag]
    def diff(a,b):return {tag:abs(v-b[tag])/abs(b[tag]) for tag,v in a.items()}
    f0=compute(0,L,1.,1.,TMIN,[(t,0.) for t in EV],52)
    for l in a.couplings:
        key=f'{l:g}';rr[key]=dict(rows=[])
        ed12,ed16,ed16tight=[fock(l,n,rt,at,step) for n,rt,at,step in ((13,2e-12,2e-14,.08),(17,2e-12,2e-14,.08),(17,2e-13,2e-15,.04))]
        rr[key].update(Fock_N12={b:z(v) for b,v in ed12.items()},Fock_N16={b:z(v) for b,v in ed16.items()},Fock_N16_tight={b:z(v) for b,v in ed16tight.items()},
                       Fock_truncation_error=diff(ed12,ed16tight),Fock_integrator_error=diff(ed16,ed16tight),
                       Born={b:z(l*c1) for b,(c1,_) in f0.items()},through_O_lambda2={b:z(l*c1+l*l*c2) for b,(c1,c2) in f0.items()})
        for h in a.h:
            tt=two(l,h,tmin=TMIN)
            ff=four(l,h,tmin=TMIN)
            tv,fv=tt.pop('F4c'),ff.pop('F4c')
            def checks(v):
                tr=str.maketrans('+-','-+')
                return dict(max_largest_time_abs=max(abs(v[b]-v[b[0].translate(tr)+b[1:]]) for b in v),
                    max_complement_abs=max(abs(v[b]-np.conj(v[b.translate(tr)])) for b in v))
            row=dict(h=h,two_pi=dict(F4c={b:z(v) for b,v in tv.items()},rel_diff_Fock=diff(tv,ed16tight),**checks(tv),**tt),
                          four_pi=dict(F4c={b:z(v) for b,v in fv.items()},rel_diff_Fock=diff(fv,ed16tight),**checks(fv),**ff),
                          four_vs_two_rel=diff(fv,tv))
            rr[key]['rows'].append(row)
            print(l,h,max(row['two_pi']['rel_diff_Fock'].values()),max(row['four_pi']['rel_diff_Fock'].values()),max(row['four_vs_two_rel'].values()),flush=True)
        if len(a.h)>1:
            sa,sb=rr[key]['rows'][-2:]
            ha,hb=sa['h'],sb['h']; assert ha>hb
            out={}
            for f in ['two_pi','four_pi']:
                va,vb=[{k:complex(*v) for k,v in x[f]['F4c'].items()} for x in (sa,sb)]
                ex={b:(ha*ha*vb[b]-hb*hb*va[b])/(ha*ha-hb*hb) for b in va}
                out[f]=dict(F4c={b:z(v) for b,v in ex.items()},rel_diff_Fock=diff(ex,ed16tight),rel_difference_vs_finest=diff(ex,vb))
            rr[key]['second_order_extrapolation']=out
    scripts=['born_bubble.py','two_pi_k0.py','four_pi_k0.py',Path(__file__).name]
    out={'description':__doc__,'K':0,'L':L,'m':1.,'T':1.,'tmin':TMIN,'events':[[t,0.] for t in EV],
         'source_sha256':{s:hashlib.sha256(Path(s).read_bytes()).hexdigest() for s in scripts},'runtime_seconds':time.monotonic()-start,
         'python':platform.python_version(),'numpy':np.__version__,'scipy':__import__('scipy').__version__,'runs':rr}
    a.out.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
