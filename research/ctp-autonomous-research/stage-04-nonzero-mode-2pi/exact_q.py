"""Independent Q-sector Fock check by translating both left external insertions."""
import argparse
import json
from pathlib import Path
from hashlib import sha256
import numpy as np
from reference import fock,A

def cli():
    p=argparse.ArgumentParser();p.add_argument('out',type=Path);p.add_argument('--lam',type=float,default=2.)
    a=p.parse_args();lam=a.lam; L=20.; Qs=range(-2,3)
    def calc(S):
        shifts=np.arange(S)*L/S
        rows=[fock(lam,14,[(t,x+(d if j<2 else 0.)) for j,(t,x) in enumerate(A)],rt=2e-13,at=2e-15,step=.04)[0] for d in shifts]
        out={Q:{tag:np.dot([r[tag] for r in rows],np.exp(-2j*np.pi*Q*shifts/L))/S for tag in rows[0]} for Q in Qs}
        return rows,out
    enc=lambda z:{k:[float(v.real),float(v.imag)] for k,v in z.items()}
    rs,Qs5=calc(5);rr,Qs7=calc(7)
    err=max(abs(v-Qs7[Q][b]) for Q,p in Qs5.items() for b,v in p.items())
    assert err<5e-13,err
    reconstruct=max(abs(sum(Qs7[Q][b] for Q in Qs)-v) for b,v in rr[0].items())
    assert reconstruct<5e-13,reconstruct
    out=dict(K=1,Nmax=14,l=lam,L=L,m=1.,T=1.,tmin=-4.,events=A,
        scripts={s:sha256(Path(s).read_bytes()).hexdigest() for s in ('reference.py','sd_residual.py','exact_q.py')},
        shift_raw_5=[enc(v) for v in rs],shift_raw_7=[enc(v) for v in rr],
        Q5={Q:enc(v) for Q,v in Qs5.items()},Q7={Q:enc(v) for Q,v in Qs7.items()},
        max_abs_5_vs_7=err,max_reconstruct_7=reconstruct)
    a.out.write_text(json.dumps(out,indent=2)+'\n')
    print('max Q sector 5-vs-7',err,'max reconstruct error',reconstruct)
if __name__=='__main__':cli()
