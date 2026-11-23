"""Reconstructed Phase-VII-compatible CF4/Chebyshev runner (not original recovered code)."""
from __future__ import annotations
import argparse, hashlib, math, os, platform, time
import json
from pathlib import Path
import numpy as np
from scipy.special import jv
import scipy
import sd_residual as pinned
import phase7_fast_ops as restored

BASE='eab1f5b3198d03eafc3f162299c9038956beaa0b'
RESTORE='e8e720bcb253ddd9bfd657e42a0902226a5a2dbb07a81c3b39560f1bd78a8b87'
EV=[(1.,0.),(.5,1.),(-.5,2.),(-1.,3.)]
PHYS={'m':1.,'L':20.,'lambda0':.1,'T':1.,'t_init':-12.}
HS=[.00625,.003125,.0015625]
PAIRS=[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]
REPORT_KEYS=('F4c','M33','connected_composite_source','analytic_residual','dxx_F4c')
def lab(m):return ''.join('+' if m>>i&1 else '-' for i in range(4))
def z(q):return [float(q.real),float(q.imag)]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

class ChebCF4:
    """Two Gauss nodes, two time-ordered exponentials, Chebyshev/Bessel action."""
    def __init__(self,E,V,profile,tol=2e-17):
        self.center=float(.5*(max(E)+min(E)))
        self.E=E-self.center;self.enorm=float(np.max(np.abs(self.E)))
        self.V=V.astype(complex);self.vnorm=float(np.max(np.asarray(abs(V).sum(axis=1))))
        self.profile=profile;self.tol=tol;self.max_order=0;self.actions=0;self.max_scaled_time=0.
    def action(self, w, dt, q):
        self.actions+=1
        sh=(-1,)+(1,)*(w.ndim-1)
        if q==0.:
            return np.exp(-.5j*dt*(self.E+self.center)).reshape(sh)*w
        radius=.5*self.enorm+abs(q)*self.vnorm
        if radius==0.:return np.exp(-.5j*dt*self.center)*w
        ab=abs(dt)*radius;self.max_scaled_time=max(self.max_scaled_time,ab)
        x=ab/2.
        # |J_n(ab)| <= x**n/n! for real ab, using the normalized Poisson integral.
        # A Hermitian rescaled generator has ||T_n||<=1; bound its discarded tail.
        n=0
        term=x
        while x/(n+2)>=1. or 2*term/(1-x/(n+2))>self.tol:
            n+=1;term*=x/(n+1)
        self.max_order=max(n,self.max_order)
        coeff=jv(np.arange(n+1),ab).astype(complex)
        coeff[1:]*=2.*((-1j if dt>=0 else 1j)**np.arange(1,n+1))
        d=(.5*self.E/radius).reshape(sh)
        v=jv(0,ab)*w
        if n:
            def a(v):return d*v+(q/radius)*(self.V@v)
            t0,t1=w,a(w)
            v+=coeff[1]*t1
            for j in range(2,n+1):
                t0,t1=t1,2*a(t1)-t0
                v+=coeff[j]*t1
        return np.exp(-.5j*dt*self.center)*v
    def __call__(self,w,t0,t1,hmax):
        if t0==t1:return w.copy()
        if hmax<=0:raise ValueError('integrator maximum step must be positive')
        n=max(1,math.ceil(abs(t1-t0)/hmax-1e-13))
        dt=(t1-t0)/n
        g=math.sqrt(3.)/6.;a=(3-2*math.sqrt(3.))/12.;b=(3+2*math.sqrt(3.))/12.
        for j in range(n):
            mid=t0+(j+.5)*dt
            left,right=self.profile(mid-g*dt),self.profile(mid+g*dt)
            w=self.action(w,dt,b*left+a*right)
            w=self.action(w,dt,a*left+b*right)
        return w

def build(K,N,phys,events):
    if len(events)!=4 or len({t for t,x in events})!=4 or events[0][0]<=max(t for t,x in events[1:]):
        raise ValueError('the four times must be distinct, and event 0 must be latest')
    c=pinned.Correlator.__new__(pinned.Correlator)
    c.K,c.N=K,N;c.mass,c.L,c.lam,c.T,c.t0=(phys[k] for k in ('m','L','lambda0','T','t_init'))
    c.events=events
    c.modes,c.basis,c.index,c.V,c.E,c.ws,c.herm=restored.make_hamiltonian_fast(K,N,c.mass,c.L)
    c.phis={};c.dxxphis={}
    c.X=[c.phi(t,x) for t,x in events]
    t,x=events[0]
    c.C=restored.projected_normal_ordered_cubic_fast(c.modes,c.basis,c.index,c.mass,c.L,x)
    c.DX=c.dxxphi(t,x)
    c.psi0=np.zeros(len(c.basis),complex);c.psi0[c.index[(0,)*len(c.modes)]]=1.
    return c

def early(c,P,dt):
    masks=[0];u=c.psi0[:,None];last=c.t0
    for i in sorted(range(1,4),key=lambda i:c.events[i][0]):
        now=c.events[i][0]
        u=P(u,last,now,dt);u=np.hstack((u,c.X[i]@u));masks+=[m|(1<<i) for m in masks]
        last=now
    u=P(u,last,c.events[0][0],dt)
    return u[:,np.argsort(masks)].copy()

def states(even,odd):
    return {m: (odd if m&1 else even)[:,m>>1] for m in range(16)}

def msk(ij):
    return sum(1<<v for v in ij)
def rev(tag):return tag.translate(str.maketrans('+-','-+'))

def moment(table, m, sub):
    ket=m&sub;bra=sub^ket
    return np.vdot(table[bra],table[ket])

def conn(table, orig, m):
    raw=moment(table,m,15)
    left3=[moment(table,m,msk(a)) for a,b in PAIRS]
    right3=[moment(orig,m,msk(b)) for a,b in PAIRS]
    return raw-sum(a*b for a,b in zip(left3,right3)),raw,left3,right3

def central_checks(ob):
    def mm(key):return {b:complex(*a[key]) for b,a in ob.items()}
    return {key+'_'+suffix:max(abs(x-d[other(b)].conjugate() if suffix=='complement' else x-d[other(b)]) for b,x in d.items())
            for key in ('F4c','M33','connected_composite_source','analytic_residual','M31_j1','M31_j2','M31_j3','F2_kl_23','F2_kl_13','F2_kl_12')
            for d in [mm(key)] for suffix,other in [('complement',rev),('largest_time',lambda tag:rev(tag[0])+tag[1:])]}

def measure(c,dt,checkpoint=None,hs=HS,verbose=True):
    tic=time.monotonic();t=c.events[0][0]
    P=ChebCF4(c.E,c.V,lambda t:pinned.pulse(t,c.lam,c.T))
    bb=early(c,P,dt)
    lt=pinned.pulse(t,c.lam,c.T)
    phi=c.X[0]
    def H(w):return c.E[:,None]*w+lt*(c.V@w)
    def comm(w):return H(phi@w)-phi@H(w)
    def D2(w):
        return -H(comm(w))+comm(H(w))-(2j/c.T)*np.tanh(t/c.T)*lt*(c.V@(phi@w)-phi@(c.V@w))
    yy=phi@bb; cc=c.C@bb; dd=c.DX@bb; aa=D2(bb)
    s0,sc,sx,s2=[states(bb,u) for u in (yy,cc,dd,aa)]
    ob={}
    for m in range(16):
        F,F4,_,_=conn(s0,s0,m)
        S,M33,left,right=conn(sc,s0,m)
        Dxx=conn(sx,s0,m)[0]
        Dtt=conn(s2,s0,m)[0];K=Dtt-Dxx+c.mass*c.mass*F
        S6=lt*S/6.
        R=K+S6;scale=abs(K)+abs(S6)+1e-30
        ob[lab(m)]={'F4c':z(F),'F4':z(F4),'M33':z(M33),'M31':dict(zip(('01','02','03'),map(z,left))),
            'F2':dict(zip(('23','13','12'),map(z,right))),
            **{f'M31_j{j}':z(x) for j,x in enumerate(left,1)},
            **{f'F2_kl_{ij}':z(x) for ij,x in zip(('23','13','12'),right)},
            'M31_F2_terms':[z(a*b) for a,b in zip(left,right)],
            'connected_composite_source':z(S),'source_term':z(S6),'analytic_dtt_F4c':z(Dtt),'dxx_F4c':z(Dxx),
            'analytic_box_plus_mass_F4c':z(K),'analytic_residual':z(R),'analytic_scale':float(scale),
            'analytic_relative_residual':float(abs(R)/scale)}
    rows=[];mats={}
    for h in hs:
        by_step=[]
        for sign in (1,-1):
            moved=P(bb,t,t+sign*h,dt);mats[f'stencil_{h:g}_{sign:+}']=moved
            table=states(moved,phi@moved)
            by_step.append([conn(table,table,m)[0] for m in range(16)])
        rr={}
        for m in range(16):
            tag=lab(m);center=ob[tag];F=complex(*center['F4c'])
            fp,fm=(u[m] for u in by_step)
            d=((fp-F)-(F-fm))/(h*h)
            B=d-complex(*center['dxx_F4c'])+c.mass*c.mass*F
            S6=complex(*center['source_term'])
            R=B+S6;scale=abs(B)+abs(S6)+1e-30
            rr[tag]={**center,'F4c_t_plus_h':z(fp),'F4c_t_minus_h':z(fm),'dtt_F4c':z(d),
               'box_plus_mass_F4c':z(B),'SD_residual':z(R),'SD_residual_abs':float(abs(R)),'SD_scale':float(scale),
               'SD_relative_residual':float(abs(R)/scale),
               'fd_error_vs_analytic':float(abs(R-complex(*center['analytic_residual'])))}
        rows.append({'h':h,'components':rr,'max_relative':max(a['SD_relative_residual'] for a in rr.values()),
                     'max_fd_error':max(a['fd_error_vs_analytic'] for a in rr.values()),
                     'min_SD_scale':min(a['SD_scale'] for a in rr.values())})
    rich=[]
    for a,b in zip(rows,rows[1:]):
        assert a['h']==2*b['h']
        d={}
        for tag in ob:
            r1=complex(*a['components'][tag]['SD_residual']);r2=complex(*b['components'][tag]['SD_residual'])
            e=(4*r2-r1)/3
            d[tag]={'signed_R0':z(e),'abs_minus_analytic':abs(e-complex(*ob[tag]['analytic_residual']))}
        rich.append({'coarse_h':a['h'],'fine_h':b['h'],'branches':d,'max_abs_minus_analytic':max(q['abs_minus_analytic'] for q in d.values())})
    cp=None
    if checkpoint:
        np.savez(checkpoint,t_center=np.array(t),integrator_max_step=np.array(dt),
                 early_masks=np.arange(0,16,2),early_fields_at_t1=bb,phi_at_t1=yy,
                 cubic_at_t1=cc,dxx_at_t1=dd,analytic_dtt_at_t1=aa,**mats)
        cp={'filename':checkpoint.name,'size_bytes':checkpoint.stat().st_size,'sha256':sha(checkpoint)}
    out={'integrator_max_step':dt,'central_observables':ob,'central_checks':central_checks(ob),'rows':rows,
         'richardson_complex':rich,'checkpoint':cp,'vacuum_norm_t1':float(np.vdot(bb[:,0],bb[:,0]).real),
         'max_analytic_relative_residual':max(q['analytic_relative_residual'] for q in ob.values()),
         'max_analytic_abs_residual':max(abs(complex(*q['analytic_residual'])) for q in ob.values()),
         'engine':{'order':4,'maximum_chebyshev_order':P.max_order,'maximum_scaled_time_bound':P.max_scaled_time,
                   'exp_actions':P.actions,'action_truncation_tol_in_exact_arithmetic':P.tol},
         'wall_seconds':time.monotonic()-tic}
    if verbose: print('dt',dt,'analytic_max',out['max_analytic_relative_residual'],'FD_max',[(r['h'],r['max_relative'],r['max_fd_error']) for r in rows],'elapsed',out['wall_seconds'],flush=True)
    return out

def compare(a,b,keys):
    result={}
    for key in keys:
        rr={}
        for tag,ac in a.items():
            d,e=complex(*ac[key]),complex(*b[tag][key])
            scale=abs(e)
            rr[tag]={'abs':abs(d-e),'reference_abs':scale,'relative':(abs(d-e)/scale if scale else None)}
        result[key]={'branches':rr,'min_reference_abs':min(d['reference_abs'] for d in rr.values()),
                     'max_abs':max(d['abs'] for d in rr.values()),
                     'max_relative':max([d['relative'] for d in rr.values() if d['relative'] is not None],default=None)}
    return result

def main():
    a=argparse.ArgumentParser()
    a.add_argument('--K',type=int,default=5);a.add_argument('--N',type=int,default=6)
    a.add_argument('--steps',type=float,nargs='+',default=[.1,.05])
    a.add_argument('--baseline',type=Path);a.add_argument('--output',type=Path,required=True)
    a.add_argument('--checkpoints',action='store_true')
    args=a.parse_args()
    tic=time.monotonic();c=build(args.K,args.N,PHYS,EV)
    print('operator dimension',len(c.basis),'wall seconds',time.monotonic()-tic,flush=True)
    ref=None
    if args.baseline:
        ref=json.loads(args.baseline.read_text())
        assert(ref['K'],ref['Nmax'],ref['events'])==(args.K,args.N,[list(a) for a in EV])
        assert all(ref[k]==PHYS[k] for k in ('m','L','lambda0','T'))
        assert ref['basis_dimension']==len(c.basis)
    out={'K':args.K,'Nmax':args.N,'events':EV,**PHYS,'basis_dimension':len(c.basis),'V_nnz':c.V.nnz,
         'C_nnz':c.C.nnz,'interaction_hermiticity_error':c.herm,
         'branch_convention':'+ chronological ket, - complementary anti-ordered bra; event order 0,1,2,3',
         'adapter':'Reconstructed CF4/Chebyshev; exact restored Phase VII optimized H and cubic; not original later runner',
         'phase7_git_base':BASE,'phase7_original_archive_sha256':RESTORE,'native_source_sha256':{
              'gate_cf4.py':sha(__file__),'sd_residual.py':sha(pinned.__file__),'phase7_fast_ops.py':sha(restored.__file__)},
         'baseline_dop853_sha256':sha(args.baseline) if args.baseline else None,
         'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),
               'machine':platform.machine(),'cpu_count':os.cpu_count(),'OPENBLAS_NUM_THREADS':os.getenv('OPENBLAS_NUM_THREADS'),
               'OMP_NUM_THREADS':os.getenv('OMP_NUM_THREADS')},'runs':[]}
    for dt in args.steps:
        cp=args.output.parent/f'k{args.K}-n{args.N}-cf4-dt{dt:g}-states.npz' if args.checkpoints else None
        run=measure(c,dt,cp)
        if ref:
            run['DOP853_central_comparison']=compare(run['central_observables'],ref['rows'][0]['components'],REPORT_KEYS[:-1])
            for b_row,r_row in zip(ref['rows'],run['rows']):
                assert b_row['h']==r_row['h']
                r_row['DOP853_stencil_comparison']=compare(r_row['components'],b_row['components'],('SD_residual',))
        if out['runs']:
            previous=out['runs'][-1]
            run['previous_max_step_comparison']=compare(previous['central_observables'],run['central_observables'],REPORT_KEYS)
            run['previous_max_step']=previous['integrator_max_step']
        out['runs'].append(run);out['wall_seconds']=time.monotonic()-tic
        args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
