"""Direct DOP853 vs separately reconstructed shared-prefix CF4, signed and per moment."""
import argparse
from hashlib import sha256
from pathlib import Path
import json, math

FIELDS=(('M33',),('M31','01'),('M31','02'),('M31','03'),
        ('F2','23'),('F2','13'),('F2','12'),('connected_composite_source',),('F4c',))
def read(d,field):
 for k in field:d=d[k]
 z=complex(*d);assert math.isfinite(abs(z))
 return z
def pack(z):return [float(z.real),float(z.imag)]
def cmp(d,c):
 assert set(d)==set(c) and len(c)==16
 reports={}
 for field in FIELDS:
  values={}
  for key in sorted(c):
   a,b=read(d[key],field),read(c[key],field)
   diff=b-a
   values[key]={'DOP853':pack(a),'CF4':pack(b),'CF4_minus_DOP853':pack(diff),
       'absolute':abs(diff),'relative_to_DOP853':abs(diff)/max(abs(a),1e-30)}
  reports[':'.join(field)]={'worst_abs':max(x['absolute'] for x in values.values()),
                           'worst_rel':max(x['relative_to_DOP853'] for x in values.values()),
                           'branches':values}
 return reports
def main():
 p=argparse.ArgumentParser();p.add_argument('direct',type=Path);p.add_argument('cf4',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
 d,c=[json.loads(f.read_text()) for f in (a.direct,a.cf4)]
 for k in ('K','Nmax','L','m','T','lambda0','events','basis_dimension','V_nnz'):
  assert d[k]==c[k],k
 assert d['t_initial']==c['t_init']
 assert d['source_commit']==c['phase7_git_base']
 assert d['operator_sha256']==c['native_source_sha256']['sd_residual.py']
 v={}
 for r in c['runs']:
  h=r['integrator_max_step'];v[str(h)]=cmp(d['components'],r['central_observables'])
  print(h,[(k,x['worst_abs'],x['worst_rel']) for k,x in v[str(h)].items()])
 assert len(v)==2 and sorted(map(float,v))==[.05,.1]
 for col,rec in v['0.05'].items():
  assert rec['worst_abs']<1e-9 and rec['worst_rel']<1e-4,(col,rec['worst_abs'],rec['worst_rel'])
 sources={x.name:sha256(x.read_bytes()).hexdigest() for x in (a.direct,a.cf4,Path(__file__))}
 a.out.write_text(json.dumps({'max_absolute_gate_fine':1e-9,'max_relative_gate_fine':1e-4,
     'checked_ctp_labels':16,'raw_moments':7,'sha256':sources,'by_CF4_max_dt':v},indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
