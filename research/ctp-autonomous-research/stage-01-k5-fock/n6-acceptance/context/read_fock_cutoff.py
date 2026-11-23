"""Read archived ++++ observables. This does not propagate any states."""
import argparse,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('n5',type=Path);p.add_argument('n6',type=Path)
args=p.parse_args()
a,b=(json.loads(f.read_text()) for f in (args.n5,args.n6))
q=a['parameters']
assert (q['Kcut'],q['Nmax'],b['K'],b['Nmax'])==(5,5,5,6)
for key in ('events','L','T','m','lambda0'): assert q[key]==b[key]
a_h=q['h']; b_h=next(row for row in b['rows'] if row['h']==a_h)
row=b_h['components']['++++']
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
dat={}
for x,y in [('F4c','F4c'),('M33_projected','M33'),('connected_source','connected_composite_source')]:
    u,v=complex(*a[x]),complex(*row[y])
    dat[y]={'N5':a[x],'N6':row[y],'complex_absolute_delta':abs(u-v),'relative_to_N6':abs(u-v)/abs(v)}
print(json.dumps({'scope':'physical Fock cutoff sensitivity at fixed K=5, A, all-forward; not fixed-Nmax numerical acceptance','common_parameters':{key:q[key] for key in ('events','L','T','m','lambda0','h')},'N5_sha256':sha(args.n5),'N6_sha256':sha(args.n6),'observables':dat},indent=2))
