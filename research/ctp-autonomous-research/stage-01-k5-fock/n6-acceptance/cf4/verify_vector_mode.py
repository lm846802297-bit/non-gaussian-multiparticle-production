"""Actual single-vector prefix path; compare against the saved block checkpoint."""
from __future__ import annotations
import argparse, hashlib, json, platform, time
from pathlib import Path
import numpy as np
import sd_residual as pinned
from gate_cf4 import build, ChebCF4, EV, PHYS, moment, lab, conn, z, sha

def build_vector_prefixes(c,p,dt):
    work={0:c.psi0.copy()}
    current=c.t0
    for j in sorted(range(1,4),key=lambda j:c.events[j][0]):
        now=c.events[j][0]
        nxt={}
        for mask,vec in work.items():
            assert vec.ndim==1
            w=p(vec,current,now,dt)
            nxt[mask]=w
            nxt[mask|(1<<j)]=c.X[j]@w
        current=now;work=nxt
    return {mask:p(vec,current,c.events[0][0],dt) for mask,vec in sorted(work.items())}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('native',type=Path); parser.add_argument('checkpoint',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    tic=time.monotonic()
    nrp=json.loads(args.native.read_text())
    assert nrp['events']==[list(a) for a in EV]
    assert all(nrp[key]==val for key,val in PHYS.items())
    with np.load(args.checkpoint) as a:
        dt=float(a['integrator_max_step'])
        run=next(r for r in nrp['runs'] if r['integrator_max_step']==dt)
        assert run['checkpoint']['sha256']==sha(args.checkpoint)
        assert run['checkpoint']['filename']==args.checkpoint.name
        c=build(nrp['K'],nrp['Nmax'],PHYS,EV)
        actions={'vector':0}
        class OnlyVector(ChebCF4):
            def __call__(self,w,*args,**kw):
                assert w.ndim==1
                actions['vector']+=1
                return super().__call__(w,*args,**kw)
        p=OnlyVector(c.E,c.V,lambda t:pinned.pulse(t,c.lam,c.T))
        bb=build_vector_prefixes(c,p,dt)
        assert tuple(bb)==tuple(a['early_masks'])
        pt={mask|1:c.X[0]@w for mask,w in bb.items()}
        ct={mask|1:c.C@w for mask,w in bb.items()}
        states={}
        for key,tab in [('early_fields_at_t1',bb),('phi_at_t1',pt),('cubic_at_t1',ct)]:
            states[key]={lab(mask):{'max_abs':float(np.max(np.abs(w-x))), 'l2_error':float(np.linalg.norm(w-x)),
                'l2_scale':float(np.linalg.norm(x))} for (mask,w),x in zip(tab.items(),a[key].T)}
    terms={}
    for m in range(16):
        P=bb|pt;Q=bb|ct
        S,M33,left,right=conn(Q,P,m)
        b=run['central_observables'][lab(m)]
        got={'F4c':conn(P,P,m)[0],'connected_composite_source':S,'M33':M33,
          **{f'M31_j{j}':x for j,x in enumerate(left,1)},
          **{f'F2_kl_{kl}':x for kl,x in zip(('23','13','12'),right)}}
        terms[lab(m)]={k:{'vector':z(x),'block':b[k],
             'abs':float(abs(x-complex(*b[k])))} for k,x in got.items()}
    max_state_error=max(x['l2_error'] for values in states.values() for x in values.values())
    max_moment_error=max(x['abs'] for vs in terms.values() for x in vs.values())
    assert max_state_error<1e-12 and max_moment_error<1e-12
    assert actions=={'vector':15}
    d={'schedule':'separately propagated single vectors, with prefixes shared before branching; central time only',
       'propagation_calls':actions,'K':c.K,'Nmax':c.N,'integrator_max_step':dt,
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'main_result_sha256':sha(args.native),
       'checkpoint_sha256':sha(args.checkpoint),'states':states,'raw_moment_comparisons':terms,
       'max_state_l2_error':max_state_error,'max_moment_abs_error':max_moment_error,
       'python':platform.python_version(),'seconds':time.monotonic()-tic}
    args.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print('PASS',max_state_error,max_moment_error,'vector calls',actions['vector'],flush=True)
if __name__=='__main__':main()
