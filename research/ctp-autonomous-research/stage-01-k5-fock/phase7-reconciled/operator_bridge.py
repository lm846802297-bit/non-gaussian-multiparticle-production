"""Stationary matrix compatibility check only; no time propagation."""
from pathlib import Path
from hashlib import sha256
import json,numpy as np
from sd_residual import make_hamiltonian
from phase7_fast_ops import make_hamiltonian_fast, projected_normal_ordered_cubic_fast
from independent_source import setup,wick_matrix
def diff(a,b):
 d=(a-b).tocsr()
 return {'max_abs':float(np.max(np.abs(d.data),initial=0)),
         'relative_frobenius':float(np.linalg.norm(d.data)/np.linalg.norm(b.data)),
         'num_above_1e-14':int(np.count_nonzero(np.abs(d.data)>1e-14))}
K,N,m,L,x=5,6,1.,20.,0.
p=make_hamiltonian_fast(K,N,m,L)
r=make_hamiltonian(K,N,m,L)
assert p[0]==r[0] and p[1]==r[1]
np.testing.assert_array_equal(p[4],r[4])
S=projected_normal_ordered_cubic_fast(p[0],p[1],p[2],m,L,x)
b,ps,occ,ns=setup(K,N)
assert np.array_equal(b[:len(p[1])],np.asarray(p[1]))
W=wick_matrix(b,ps,occ,ns,K,N,m,L,x)
d={'K':K,'Nmax':N,'L':L,'m':m,'x':x,'dimension':len(p[1]),'base_commit':'eab1f5b3198d03eafc3f162299c9038956beaa0b','hamiltonian':diff(p[3],r[3]),'cubic_vs_independent_wick':diff(S,W),'source_sha256':{p:sha256(Path(p).read_bytes()).hexdigest() for p in ['src/phase7_fast_ops.py','operator_bridge.py']}}
assert d['hamiltonian']['num_above_1e-14']==0 and d['cubic_vs_independent_wick']['num_above_1e-14']==0
Path('operator_bridge.json').write_text(json.dumps(d,indent=2)+'\n')
print(json.dumps(d,indent=2))
