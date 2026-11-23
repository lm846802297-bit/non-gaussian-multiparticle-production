"""Cross-check immutable inputs, both propagation meshes, and signed independent moments."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from gate_cf4 import compare,states,conn,lab,sha

PAIRS=(('01','23'),('02','13'),('03','12'))
DKEYS=['M33',*[f'M31_j{i}' for i in (1,2,3)],*[f'F2_kl_{ij}' for a,ij in PAIRS]]
def load(p):return json.loads(Path(p).read_text())

def main():
 a=argparse.ArgumentParser()
 a.add_argument('--native',default='k5-n6-cf4.json');a.add_argument('--direct',default='references/k5-n6-direct-moments.json')
 a.add_argument('--archive',default='references/archived-dop-k5-n6.json')
 a.add_argument('--output',default='verification.json')
 a.add_argument('--vector',default='k5-n6-cf4-single-vector-bridge.json')
 x=a.parse_args()
 d,r,v=map(load,(x.native,x.direct,x.vector))
 assert sha(x.archive)==d['baseline_dop853_sha256']
 assert d['native_source_sha256']['gate_cf4.py']==sha('gate_cf4.py')
 assert d['native_source_sha256']['sd_residual.py']==sha('dependencies/sd_residual.py')
 assert d['native_source_sha256']['phase7_fast_ops.py']==sha('dependencies/phase7_fast_ops.py')
 assert v['main_result_sha256']==sha(x.native)
 assert v['source_sha256']==sha('verify_vector_mode.py')
 assert v['propagation_calls']=={'vector':15}
 for k in ('K','Nmax','L','m','T','lambda0','basis_dimension','events'):assert d[k]==r[k]
 assert d['t_init']==r['t_initial']==-12
 assert d['K']==5 and d['Nmax']==6 and d['basis_dimension']==12376
 assert len(r['components'])==16 and len(r['independently_evaluated_contour_labels'])==16
 steps=[s['integrator_max_step'] for s in d['runs']]
 assert steps==[.1,.05]
 direct={}
 for tag,o in r['components'].items():
  direct[tag]={**o,**{f'M31_j{j}':o['M31'][ij] for j,(ij,_) in enumerate(PAIRS,1)},
                   **{f'F2_kl_{ij}':o['F2'][ij] for _,ij in PAIRS}}
 successful=[]
 for q in d['runs']:
  p=Path(q['checkpoint']['filename'])
  assert p.stat().st_size==q['checkpoint']['size_bytes'] and sha(p)==q['checkpoint']['sha256']
  obs=q['central_observables']
  assert set(obs)==set(direct)
  assert len(q['rows'])==3 and [b['h'] for b in q['rows']]==[.00625,.003125,.0015625]
  assert max(q['central_checks'].values())<1e-14
  assert abs(q['vacuum_norm_t1']-1.)<1e-12
  with np.load(p) as f:
   assert np.array_equal(f['early_masks'],np.arange(0,16,2))
   assert f['early_fields_at_t1'].shape==(12376,8)
   for key in ['phi_at_t1','cubic_at_t1','dxx_at_t1','analytic_dtt_at_t1']:
    assert f[key].shape==(12376,8)
   assert float(f['t_center'])==1 and float(f['integrator_max_step'])==q['integrator_max_step']
   s0=states(f['early_fields_at_t1'],f['phi_at_t1'])
   sc=states(f['early_fields_at_t1'],f['cubic_at_t1'])
   for m in range(16):
    o=obs[lab(m)]
    S,M,L,R=conn(sc,s0,m)
    for value,expected in [(S,o['connected_composite_source']),(M,o['M33']),*zip(L,o['M31'].values()),*zip(R,[o['F2'][p] for p in ('23','13','12')])]:
     assert abs(value-complex(*expected))<1e-15
  cmp=compare(obs,direct,[*DKEYS,'F4c','connected_composite_source'])
  assert max(x['max_relative'] for x in cmp.values())<1e-3
  ratios=[m['max_fd_error']/n['max_fd_error'] for m,n in zip(q['rows'],q['rows'][1:])]
  p2_orders=list(map(np.log2,ratios))
  successful.append({'dt':q['integrator_max_step'],'independent_direct_comparison':cmp,'FD_max_error_halving_ratios':ratios,'FD_max_error_effective_orders':p2_orders})
 finest=d['runs'][-1]
 assert min(successful[-1]['FD_max_error_halving_ratios'])>3.5
 assert finest['previous_max_step_comparison']['F4c']['max_relative']<1e-3
 assert finest['previous_max_step_comparison']['connected_composite_source']['max_relative']<1e-3
 assert finest['DOP853_central_comparison']['analytic_residual']['max_abs']<1e-13
 assert v['max_state_l2_error']<1e-12 and v['max_moment_abs_error']<1e-12
 summ={'status':'PASS','native_sha256':sha(x.native),'direct_sha256':sha(x.direct),'vector_bridge_sha256':sha(x.vector),'runs':successful}
 Path(x.output).write_text(json.dumps(summ,indent=2,sort_keys=True)+'\n')
 print('PASS: all 16 signed branches × seven primitive moments, three h at two independent dt, saved states, vector bridge')
 for q in successful:
  print('dt',q['dt'],'O(h^2) ratio',q['FD_max_error_halving_ratios'])
  for key,v in q['independent_direct_comparison'].items():print(key,v['max_relative'],v['max_abs'],v['min_reference_abs'])
if __name__=='__main__':main()
