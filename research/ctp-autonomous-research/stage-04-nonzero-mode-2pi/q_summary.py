"""Compare separately projected Fock and 2PI/BS total pair-momentum sectors."""
import argparse
import json
import hashlib
from pathlib import Path
import numpy as np
from born_bubble import compute
from explicit_wick import wick_by_q
from pilot import cz,z,diff,rich

def cli():
    p=argparse.ArgumentParser()
    p.add_argument('pilot',type=Path);p.add_argument('qref',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
    pil=json.loads(a.pilot.read_text()); qref=json.loads(a.qref.read_text())
    assert all(pil[n]==qref[n] for n in ('K','L','m','T','tmin','events'))
    K,L,m,T,tmin,ev=[pil[n] for n in ('K','L','m','T','tmin','events')]; lam=qref['l']
    ref={int(Q):cz(v) for Q,v in qref['Q7'].items()}; Qs=sorted(ref)
    run=pil['runs'][f'{lam:g}']; tags=list(run['Fock']); intervals=run['grids'][-2:]
    assert max(abs(sum(v[b] for v in ref.values())-complex(*run['Fock'][b])) for b in tags)<5e-13
    shifts=np.arange(4*K+1)*L/(4*K+1)
    born=[]
    for sh in shifts:
        born.append(compute(K,L,m,T,tmin,[(t,x+sh if j<2 else x) for j,(t,x) in enumerate(ev)],48))
    cont={Q:{b:sum((lam*v[b][0]+lam*lam*v[b][1])*np.exp(-2j*np.pi*Q*sh/L) for sh,v in zip(shifts,born))/len(shifts) for b in tags} for Q in Qs}
    corrected=[]
    for gr in intervals:
        h=gr['step']
        wick=wick_by_q(ev,h,K,L,m,T,tmin)
        d={Q:{b:complex(*gr['F4c_by_Q'][str(Q)][b])+cont[Q][b]-(lam*c1+lam*lam*c2) for b,(c1,c2) in q.items()} for Q,q in wick.items()}
        assert max(abs(sum(v[b] for v in d.values())-complex(*gr['corrected'][b])) for b in tags)<5e-13
        corrected.append(d)
    out={}
    for Q in Qs:
        ex=rich(corrected[0][Q],corrected[1][Q],intervals[0]['step'],intervals[1]['step'])
        out[Q]=dict(Fock=z(ref[Q]),extrapolated=z(ex),finest_corrected=z(corrected[1][Q]),
                    rel_errors=diff(ex,ref[Q]),abs_errors={b:abs(ex[b]-ref[Q][b]) for b in tags},
                    finest_pair_shift=diff(corrected[0][Q],corrected[1][Q]))
        print(Q,'max rel error',max(out[Q]['rel_errors'].values()),'max abs',max(out[Q]['abs_errors'].values()))
    sig=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    a.out.write_text(json.dumps(dict(K=K,L=L,m=m,T=T,tmin=tmin,events=ev,lam=lam,grids=[g['step'] for g in intervals],
        Q_fock_5_vs_7_shift_abs=qref['max_abs_5_vs_7'],SHA256={'pilot':sig(a.pilot),'qref':sig(a.qref),'source':sig(__file__)},compared=out),indent=2)+'\n')
if __name__=='__main__':cli()
