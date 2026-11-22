"""Remove known fixed-order quadrature error, then diagnose remainder mesh error."""
import argparse
import json
from pathlib import Path
from test_discrete_weak_order import wick_direct

def z(c):return [c.real,c.imag]
def vals(k):return {i:complex(*x) for i,x in k.items()}
def rr(m,n):return {b:abs(m[b]-n[b])/abs(n[b]) for b in m}
def extr(a,b,ha,hb):return {k:(ha*ha*b[k]-hb*hb*a[k])/(ha*ha-hb*hb) for k in a}

def main():
    p=argparse.ArgumentParser();p.add_argument('source',type=Path,nargs='+');p.add_argument('out',type=Path);a=p.parse_args()
    q={}
    for src in a.source:
        src=json.loads(src.read_text())
        for l,v in src['runs'].items():
            if l in q:
                assert v['Fock_N16_tight']==q[l]['Fock_N16_tight']
                q[l]['rows']+=v['rows']
            else:q[l]=v
    out={}
    for l,d in q.items():
        lam=float(l);corr={}
        ref=vals(d['Fock_N16_tight']);continuum=vals(d['through_O_lambda2'])
        rows=sorted(d['rows'],key=lambda x:-x['h'])
        for r in rows:
            h=r['h']
            wk=wick_direct(h)
            defect={b:continuum[b]-(lam*u+lam*lam*v) for b,(u,v) in wk.items()}
            corr[h]={}
            for m in ('two_pi','four_pi'):
                mesh=vals(r[m]['F4c'])
                corr[h][m]={b:v+defect[b] for b,v in mesh.items()}
        ha,hb=[r['h'] for r in rows[-2:]]
        out[l]={'Fock_truncation_max':max(d['Fock_truncation_error'].values()),
                'Fock_integrator_max':max(d['Fock_integrator_error'].values()),'meshes':list(corr)}
        for m in ('two_pi','four_pi'):
            ext=extr(corr[ha][m],corr[hb][m],ha,hb)
            ex0=extr(vals(rows[-2][m]['F4c']),vals(rows[-1][m]['F4c']),ha,hb)
            out[l][m]=dict(finest_defect_removed_rel_err_max=max(rr(corr[hb][m],ref).values()),
                step_uncertainty_max=max(rr(corr[ha][m],corr[hb][m]).values()),
                plain_h2_extrap_rel_err_max=max(rr(ex0,ref).values()),
                defect_removed_h2_extrap_rel_err_max=max(rr(ext,ref).values()),
                all_16_relative_errors=rr(ext,ref),F4c_extrapolated={b:z(v) for b,v in ext.items()},
                branch4=ext['++++'].__repr__())
    a.out.write_text(json.dumps(out,indent=2)+'\n')
    for l,v in out.items():
        print(l,[(m,d['finest_defect_removed_rel_err_max'],d['step_uncertainty_max'],d['defect_removed_h2_extrap_rel_err_max']) for m,d in v.items() if m in ('two_pi','four_pi')])
if __name__=='__main__':main()
