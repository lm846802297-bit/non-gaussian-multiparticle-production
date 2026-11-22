from pathlib import Path
import json,math,hashlib
from itertools import product
def z(x,Q):return complex(*x[Q])
keys={''.join(x) for x in product('+-',repeat=4)}
tr=str.maketrans('+-','-+')
expected={'k5-n6-C':(5,6,20.,'C'),'k5-n7-C':(5,7,20.,'C'),'k6-n6-A':(6,6,20.,'A'),'k6-n6-A-L24':(6,6,24.,'A')}
evs={'A':[[1.,0.],[.5,1.],[-.5,2.],[-1.,3.]],'C':[[1.2,0.],[.4,.7],[-.3,1.7],[-1.1,3.2]]}
dataset={}
info={}
for key,p in [*( (k,Path(f'run/{k}.json')) for k in expected), *((f'k5-n{n}-A',Path(f'/tmp/ctpproject/run/k5-n{n}'+('-segmented' if n==6 else '')+'.json')) for n in (6,7))]:
 if not p.exists():continue
 d=json.loads(p.read_text())
 if key in expected:
  K,N,L,e=expected[key]
  assert d['runner_sha256']==hashlib.sha256(Path('next_scan.py').read_bytes()).hexdigest()
  assert (d['K'],d['Nmax'],d['L'])==(K,N,L) and d['events']==evs[e]
  assert [x['h'] for x in d['rows']]==[.00625,.003125,.0015625]
 assert d['basis_dimension']==math.comb(2*d['K']+1+d['Nmax'],d['Nmax'])
 assert d['source_commit']=='eab1f5b3198d03eafc3f162299c9038956beaa0b'
 assert math.isfinite(d['max_analytic_abs']) and d['hermiticity_error']<1e-14
 dataset[key]=(d, d['rows'][1]['components'])
 largest=conj=fd_largest=0.
 for row in d['rows']:
  c=row['components'];assert set(c)==keys
  for k in keys:
   for u in ('F4c','M33','connected_composite_source','analytic_residual','SD_residual'):
    v=z(c[k],u);assert math.isfinite(abs(v))
    conj=max(conj,abs(v-z(c[k.translate(tr)],u).conjugate()))
    err=abs(v-z(c[k[0].translate(tr)+k[1:]],u))
    if u=='SD_residual':
     fd_largest=max(fd_largest,err);assert err<2e-15/row['h']**2
    else: largest=max(largest,err)
 assert largest<1e-12 and conj<1e-12
 a=dataset[key][1]
 info[key]={'dimension':d['basis_dimension'],'physical_cutoff':2*math.pi*d['K']/d['L'],'max_relative':[r['max_relative'] for r in d['rows']],'max_fd_error':[r['max_fd_error'] for r in d['rows']],'max_analytic_relative':max(c['analytic_relative_residual'] for c in a.values()),'max_analytic_abs':d['max_analytic_abs'],'largest_time_max_abs_no_stencil':largest,'largest_time_max_abs_fd':fd_largest,'complement_max_abs':conj,'min_source_abs':min(abs(z(x,'connected_composite_source')) for x in a.values()),'wall_seconds':d['runtime']['wall_seconds']}
comparisons={}
for name,i,j in [
 ('C_Fock','k5-n6-C','k5-n7-C'),
 ('A_Fock','k5-n6-A','k5-n7-A'),
 ('A_K_at_fixed_L_N6','k5-n6-A','k6-n6-A'),
 ('A_L_K_over_L_matched_N6','k5-n6-A','k6-n6-A-L24'),
 ('A_L_at_fixed_K_N6','k6-n6-A','k6-n6-A-L24')]:
 if i not in dataset or j not in dataset:continue
 x,y=(dataset[k][1] for k in (i,j))
 comparisons[name]={q:{'max_absolute':max(abs(z(y[k],q)-z(x[k],q)) for k in keys),'max_relative_to_first':max(abs(z(y[k],q)-z(x[k],q))/max(abs(z(x[k],q)),1e-30) for k in keys),'all_forward_relative_to_first':abs(z(y['++++'],q)-z(x['++++'],q))/max(abs(z(x['++++'],q)),1e-30)} for q in ('F4c','M33','connected_composite_source')}
print(json.dumps({'cases':info,'comparisons':comparisons},indent=2,sort_keys=True))
