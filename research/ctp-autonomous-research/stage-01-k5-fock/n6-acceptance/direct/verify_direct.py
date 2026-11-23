"""Replayable verification of the independent direct-moment native output."""
import argparse
from pathlib import Path
import hashlib
import json
import math
from itertools import product

LABELS={''.join(a) for a in product('+-',repeat=4)}
PAIRS=(('01','23'),('02','13'),('03','12'))
SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def Z(x):
 a=complex(*x);assert math.isfinite(abs(a));return a
def main():
 p=argparse.ArgumentParser()
 p.add_argument('native',type=Path);p.add_argument('operator',type=Path);p.add_argument('archived',type=Path)
 p.add_argument('--out',type=Path)
 a=p.parse_args()
 d=json.loads(a.native.read_text());components=d['components']
 assert set(components)==set(d['independently_evaluated_contour_labels'])==LABELS
 assert (d['K'],d['Nmax'],d['L'],d['m'],d['T'],d['lambda0'],d['t_initial'])==(5,6,20.,1.,1.,.1,-12.)
 assert d['events']==[[1.,0.],[.5,1.],[-.5,2.],[-1.,3.]]
 assert d['basis_dimension']==12376 and d['V_nnz']==379654
 assert d['runner_sha256']==SHA(Path(__file__).with_name('direct_moments.py'))
 assert d['operator_sha256']==SHA(a.operator)
 assert d['archived_comparison']['archived_sha256']==SHA(a.archived)
 c_ref=json.loads(a.archived.read_text())['rows'][1]['components']
 checks={name:0. for name in ('source_reconstruction','F4_reconstruction','raw_product_reconstruction')}
 comp={name:{'max_abs':0.,'max_rel':0.} for name in ('M33','F4c','connected_composite_source')}
 for label,c in components.items():
  assert set(c['M31'])=={'01','02','03'} and set(c['F2'])=={'12','13','23'}
  terms=[Z(c['M31'][i])*Z(c['F2'][j]) for i,j in PAIRS]
  checks['raw_product_reconstruction']=max(checks['raw_product_reconstruction'],*(abs(v-Z(z)) for v,z in zip(terms,c['M31_F2_terms'])))
  checks['source_reconstruction']=max(checks['source_reconstruction'],abs(Z(c['M33'])-sum(terms)-Z(c['connected_composite_source'])))
  F4_disconnected=sum(Z(c['F2_0j_for_F4'][i])*Z(c['F2'][j]) for i,j in PAIRS)
  checks['F4_reconstruction']=max(checks['F4_reconstruction'],abs(Z(c['F4raw'])-F4_disconnected-Z(c['F4c'])))
  for name,v in comp.items():
   x,y=Z(c[name]),Z(c_ref[label][name])
   v['max_abs']=max(v['max_abs'],abs(x-y));v['max_rel']=max(v['max_rel'],abs(x-y)/abs(y))
 for value in checks.values():assert value<2e-16
 for v in comp.values():
  assert v['max_abs']<5e-10 and v['max_rel']<1e-6
 for v in d['checks'].values():assert max(v.values())<5e-10
 summary={'native_sha256':SHA(a.native),'operator_sha256':SHA(a.operator),'archived_sha256':SHA(a.archived),
          'checked_ctp_labels':16,'checked_raw_moments_per_label':7,'reconstruction':checks,'archived':comp,
          'branch_symmetries':d['checks'],'counts':d['counts']}
 if a.out:a.out.write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
 print(json.dumps(summary,indent=2,sort_keys=True))
if __name__=='__main__':main()
