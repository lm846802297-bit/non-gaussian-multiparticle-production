"""Read-only audit of existing matched K=5 archives; no evolution."""
import argparse, hashlib, json, math
from pathlib import Path
from itertools import product

def require(b,msg):
    if not b: raise ValueError(msg)
def sha(b): return hashlib.sha256(b).hexdigest()
def z(v):
    require(isinstance(v,list) and len(v)==2 and all(type(x) in (int,float) and math.isfinite(x) for x in v), 'invalid complex')
    return complex(*v)
def zz(v): return [v.real,v.imag]
def close(a,b,e): return abs(a-b)<=e
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('source',type=Path);p.add_argument('output',type=Path)
p.add_argument('--flat',action='store_true',help='use logical input names in source/')
a=p.parse_args()
d={};prov=json.loads(Path(__file__).with_name('inputs.json').read_text())
for k,v in prov.items():
    b=(a.source/(k if a.flat else v['path'])).read_bytes()
    require(len(b)==v['size'] and hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==v['sha'], 'git input mismatch '+k)
    v['SHA256']=sha(b)
    if k in ('n7','n6','cf','wick'):d[k]=json.loads(b)
n6,n7,cf,w=(d[x] for x in ('n6','n7','cf','wick'))
keys=('K','m','L','T','lambda0','events')
for k in keys:require(n6[k]==n7[k], 'DOP parameters mismatch '+k)
for k in ('m','L','T','lambda0'):require(cf['physical'][k]==n7[k], 'CF parameters mismatch '+k)
require(cf['physical']['t_initial']==-12. and cf['events']==n7['events'] and n6['K']==n7['K']==cf['K']==5 and (n6['Nmax'],n7['Nmax'],cf['Nmax'])==(6,7,7),'configuration mismatch')
require(n7['basis_dimension']==cf['basis_dimension']==31824 and n6['basis_dimension']==12376,'basis dimension')
require(n7['source_commit']==n6['source_commit']==cf['source_commit']==w['source_commit']=='eab1f5b3198d03eafc3f162299c9038956beaa0b','base source pin')
require(cf['baseline_sha256']==prov['n7']['SHA256'] and cf['source_module_sha256']==prov['ops']['SHA256'] and cf['wick_module_sha256']==prov['wickrunner']['SHA256'] and cf['runner_sha256']==prov['cfrunner']['SHA256'],'CF runner/reference pin')
require(n7['runner_sha256']==n6['runner_sha256']==prov['n7runner']['SHA256'],'archived DOP runner pin')
labels={''.join(b) for b in product('+-',repeat=4)}
flip=str.maketrans('+-','-+')
hs=[.00625,.003125,.0015625]
F='F4c';S='connected_composite_source';M='M33';A='analytic_residual';R='SD_residual'
metrics=(F,S,M,A,R)
rows=n7['rows'];base=rows[1]['components']
require([r['h'] for r in rows]==[r['h'] for r in n6['rows']]==hs,'stencils incomplete')
g=n7['lambda0']/math.cosh(n7['events'][0][0]/n7['T'])**2/6
def summary(r):
    h=r['h']; c=r['components']
    require(set(c)==labels,'not all branches')
    sym={k:{'complement':max(abs(z(c[l][k])-z(c[l.translate(flip)][k]).conjugate()) for l in labels),
            'latest':max(abs(z(c[l][k])-z(c[l[0].translate(flip)+l[1:]][k])) for l in labels)} for k in metrics}
    require(max(x['complement'] for x in sym.values())<1e-14 and all(sym[k]['latest']<(2e-15/h**2 if k==R else 1e-12) for k in metrics),'CTP symmetry mismatch')
    by={}
    for l,q in sorted(c.items()):
        require(all(q[k]==base[l][k] for k in (F,S,M,A)),'center drift')
        s=g*z(q[S]); r0=z(q[A]); rf=z(q[R])
        rrel=abs(rf)/(abs(rf-s)+abs(s)+1e-30)
        arel=abs(r0)/(abs(r0-s)+abs(s)+1e-30)
        require(close(rrel,q['SD_relative_residual'],1e-16) and close(arel,q['analytic_relative_residual'],1e-16) and close(abs(rf-r0),q['fd_error_vs_analytic'],1e-20),'signed normalization mismatch')
        by[l]={'stencil_minus_analytic':zz(rf-r0),'abs_error':abs(rf-r0),'normalized_residual':rrel}
    require(close(max(q['normalized_residual'] for q in by.values()),r['max_relative'],1e-16),'reported worst branch mismatch')
    return {'h':h,'symmetry_maxima_abs':sym,'max_normalized_residual':r['max_relative'],'max_signed_error_abs':max(q['abs_error'] for q in by.values()),'by_branch':by}
ss=[summary(r) for r in rows]
orders={}
for u,v in zip(ss,ss[1:]):
    o={l:math.log2(u['by_branch'][l]['abs_error']/v['by_branch'][l]['abs_error']) for l in sorted(labels) if min(u['by_branch'][l]['abs_error'],v['by_branch'][l]['abs_error'])>0}
    orders[str(u['h'])]={'n_nonzero':len(o),'min':min(o.values()),'max':max(o.values()),'by_branch':o}
def delta(a,b,keys):
    return {k:{'complex_difference':zz(z(b[k])-z(a[k])),'abs':abs(z(b[k])-z(a[k])),'relative_to_first':abs(z(b[k])-z(a[k]))/abs(z(a[k]))} for k in keys}
require([q['max_step'] for q in cf['steps']]==[.1,.05] and all(q['h']==hs[1] for q in cf['steps']),'CF coverage mismatch')
cf_comparisons={}
for q in cf['steps']:
    require(close(z(q[F]),z(q['F4'])-sum(z(x) for x in q['F2_F2']),1e-17),'CF4 connected F reconstruction')
    require(close(z(q[S]),z(q[M])-sum(z(x) for x in q['M31_F2']),1e-17),'CF4 connected source reconstruction')
    cc=delta(base['++++'],q,metrics)
    require(all(close(cc[k]['abs'],q['comparison_DOP853'][k]['absolute'],1e-20) for k in metrics),'CF/DOP cross report')
    cf_comparisons[str(q['max_step'])]=cc
cf_ref=delta(*cf['steps'],(F,S,M,A))
cut={l:delta(n6['rows'][1]['components'][l],base[l],(F,S,M)) for l in sorted(labels)}
require(w['K']==5 and w['Nmax']==7 and w['comparison']['shape']==[31824,31824] and w['comparison']['nnz_over_tolerance']==0,'Wick metadata mismatch')
out={'study':'archived A case only; no new propagation','physical':{x:n7[x] for x in keys},'inputs':prov,
     'dop853':{'branches':16,'rows':ss,'max_analytic_normalized_residual':max(q['analytic_relative_residual'] for q in base.values()),'signed_stencil_orders':orders},
     'CF4_Taylor_all_forward_only':{'vs_DOP853':cf_comparisons,'refinement':cf_ref,'fixed_F_S_within_1e-3_++++':all(cf_comparisons['0.05'][k]['relative_to_first']<1e-3 for k in (F,S))},
     'n6_to_n7_physical_Fock_sensitivity':{'max_relative_to_N6':{k:max(v[k]['relative_to_first'] for v in cut.values()) for k in (F,S,M)},'by_branch':cut},
     'wick_operator_only':w['comparison'],
     'missing_for_full_N6_like_acceptance':['independent 16x7 signed primitive moments (CF4 archive stores only three products, not separated factors)','independent CF4 all-branch three-stencil and single-vector bridge at N7','N7 DOP853 archive does not retain the shifted F4c values, so signed residual normalization can be checked but the second derivative cannot be rebuilt from raw correlators']}
a.output.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+'\n')
print({k:v for k,v in out.items() if k in ('missing_for_full_N6_like_acceptance','wick_operator_only')})
