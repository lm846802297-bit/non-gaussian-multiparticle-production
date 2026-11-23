"""Direct CTP signed raw moments, with independently propagated insertion paths.

This uses archived operators, not the shared-prefix branch construction.
Each complete chronological insertion path is solved from the same vacuum
state at the earliest physical event. Only completed identical paths are cached.
"""
import argparse
from functools import cache
from itertools import product
from pathlib import Path
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import sd_residual as pin

SOURCE_PIN='eab1f5b3198d03eafc3f162299c9038956beaa0b'
EVENTS=((1.,0.),(.5,1.),(-.5,2.),(-1.,3.))
LABELS=[''.join(a) for a in product('+-',repeat=4)]
INV=str.maketrans('+-','-+')
PAIRS=((1,2,3),(2,1,3),(3,1,2))
KEYMAP={'M33':(0,1,2,3),'01':(0,1),'02':(0,2),'03':(0,3),'23':(2,3),'13':(1,3),'12':(1,2)}

def pack(z):return [float(z.real),float(z.imag)]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def run(K=5,N=6,L=20.,m=1.,T=1.,lam=.1,events=EVENTS,t_initial=-12.,rtol=2e-10,atol=2e-12,max_step=.1,progress=None):
    start=time.monotonic()
    modes,basis,index,V,E,omega,H_error=pin.make_hamiltonian(K,N,m,L)
    VC=V.astype(complex)
    phi=[pin.field_matrix(modes,basis,index,m,L,x) for _,x in events]
    cube=pin.projected_normal_ordered_cubic(modes,basis,index,m,L,events[0][1])
    t=[e[0] for e in events]
    assert len(set(t))==len(t)==4 and t[0]==max(t) and min(t)>t_initial
    o0=np.zeros(len(basis),complex);o0[index[(0,)*len(modes)]]=1.
    counter={'segments':0,'RHS':0}
    def evolve(v,a,b):
        if a==b:return v.copy()
        assert a<b
        def f(s,q):return -1j*(E*q+lam/np.cosh(s/T)**2*(VC@q))
        sol=solve_ivp(f,(a,b),v,method='DOP853',rtol=rtol,atol=atol,max_step=max_step,t_eval=[b])
        if not sol.success or sol.y.shape[-1]!=1:raise RuntimeError(sol.message)
        counter['segments']+=1
        counter['RHS']+=sol.nfev
        return sol.y[:,-1]
    t_earliest=min(t)
    initial=evolve(o0,t_initial,t_earliest)
    @cache
    def path(seq,end):
        v=initial.copy();now=t_earliest
        for i,is_cube in seq:
            v=evolve(v,now,t[i])
            v=(cube if is_cubic_event(i,is_cubic_flag=is_cube) else phi[i])@v
            now=t[i]
        return evolve(v,now,end)
    def measured(label,event_ids,cubic_flag=False):
        end=max(t[i] for i in event_ids)
        ordered=sorted(event_ids,key=t.__getitem__)
        def seq(sign):
            return tuple((i,cubic_flag and i==0) for i in ordered if label[i]==sign)
        return np.vdot(path(seq('-'),end),path(seq('+'),end))
    if progress is not None:progress('Hamiltonian/projected cubic built; independent paths started',time.monotonic()-start)
    rr={}
    for count,lab in enumerate(LABELS):
        M33=measured(lab,KEYMAP['M33'],True)
        M31={f'0{j}':measured(lab,(0,j),True) for j,_,_ in PAIRS}
        F2={f'{k}{l}':measured(lab,(k,l)) for _,k,l in PAIRS}
        F2more={f'0{j}':measured(lab,(0,j)) for j,_,_ in PAIRS}
        mm=[M31[f'0{j}']*F2[f'{k}{l}'] for j,k,l in PAIRS]
        F4raw=measured(lab,(0,1,2,3))
        F4c=F4raw-sum(F2more[f'0{j}']*F2[f'{k}{l}'] for j,k,l in PAIRS)
        rr[lab]={'M33':pack(M33),'M31':{k:pack(v) for k,v in M31.items()},
                 'F2':{k:pack(v) for k,v in F2.items()},'M31_F2_terms':[pack(v) for v in mm],
                 'connected_composite_source':pack(M33-sum(mm)),'F4c':pack(F4c),'F4raw':pack(F4raw),
                 'F2_0j_for_F4':{k:pack(v) for k,v in F2more.items()}}
        if progress is not None and count in (3,7,11,15):progress(f'{count+1}/16 signed assignments',time.monotonic()-start)
    # Every label's bra and ket are evaluated from their paths before symmetry checks.
    errors={}
    for name in ('M33','01','02','03','23','13','12','connected_composite_source','F4c'):
        def v(l):
            c=rr[l]
            if name in c:return complex(*c[name])
            return complex(*(c['M31'] if name[0]=='0' else c['F2'])[name])
        latest=min(KEYMAP.get(name,(0,)),key=lambda i:-t[i])
        errors[name]={
            'max_conjugation_abs':max(abs(v(a)-v(a.translate(INV)).conjugate()) for a in LABELS),
            'max_latest_time_abs':max(abs(v(a)-v(a[:latest]+a[latest].translate(INV)+a[latest+1:])) for a in LABELS)
        }
    runsrc=Path(__file__).resolve()
    pinpath=Path(pin.__file__).resolve()
    pinbytes=pinpath.read_bytes()
    pinned_blob=hashlib.sha1(b'blob '+str(len(pinbytes)).encode()+b'\0'+pinbytes).hexdigest()
    assert pinned_blob=='29139f6f8a35d3d874c401d240318ed2c8cadc40'
    return {'source_commit':SOURCE_PIN,'K':K,'Nmax':N,'m':m,'L':L,'T':T,'lambda0':lam,'events':events,'t_initial':t_initial,
            'rtol':rtol,'atol':atol,'max_step':max_step,'solver':'DOP853, independently propagated per chronological insertion path',
            'operator_git_blob':pinned_blob,'operator_sha256':sha(pinpath),'runner_sha256':sha(runsrc),
            'basis_dimension':len(basis),'V_nnz':V.nnz,'hamiltonian_hermiticity_error':H_error,
            'source_order':['M33','M31:01','M31:02','M31:03','F2:23','F2:13','F2:12'],
            'components':rr,'independently_evaluated_contour_labels':LABELS,'checks':errors,
            'zero_insertion_vacuum_norm_abs_error':float(abs(np.vdot(initial,initial)-1)),
            'counts':{**counter,'complete_unique_paths':path.cache_info().currsize,'exact_reuse_of_complete_paths':path.cache_info().hits},
            'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'wall_seconds':time.monotonic()-start}}

def is_cubic_event(i,is_cubic_flag=False):
    if is_cubic_flag:
        assert i==0
    return is_cubic_flag

def archived_compare(o,archived):
    p=Path(archived);v=json.loads(p.read_text())
    for field in ('K','Nmax','m','L','T','lambda0'):
        assert o[field]==v[field],field
    assert np.array_equal(o['events'],v['events'])
    cc=v['rows'][1]['components'];rr={}
    for key in ('F4c','M33','connected_composite_source'):
        errors={}
        for tag,a in o['components'].items():
            ref=complex(*cc[tag][key]);ours=complex(*a[key]);delta=ours-ref
            errors[tag]={'complex_delta':pack(delta),'absolute':abs(delta),'relative_to_archived':abs(delta)/max(abs(ref),1e-30)}
        rr[key]={'max_absolute':max(x['absolute'] for x in errors.values()),
                 'max_relative':max(x['relative_to_archived'] for x in errors.values()),'by_branch':errors}
    return {'archived_file':p.name,'archived_sha256':sha(p),'check':rr}

def cli():
    p=argparse.ArgumentParser()
    p.add_argument('out',type=Path)
    p.add_argument('--K',type=int,default=5);p.add_argument('--N',type=int,default=6)
    p.add_argument('--rtol',type=float,default=2e-10);p.add_argument('--atol',type=float,default=2e-12);p.add_argument('--max-step',type=float,default=.1)
    p.add_argument('--archived',type=Path);a=p.parse_args()
    o=run(a.K,a.N,rtol=a.rtol,atol=a.atol,max_step=a.max_step,progress=lambda what,secs:print(f'[{secs:.2f}s] {what}',flush=True))
    if a.archived is not None:o['archived_comparison']=archived_compare(o,a.archived)
    a.out.write_text(json.dumps(o,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(o['counts'],o['runtime']['wall_seconds'],flush=True)
if __name__=='__main__':cli()
