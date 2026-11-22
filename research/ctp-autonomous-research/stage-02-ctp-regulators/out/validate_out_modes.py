"""Audit native JSON, integrated factorial sum rule, and matched regulator comparisons."""
import hashlib, itertools, json, math, sys
from pathlib import Path

def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def rel(x,y):
    return (y-x)/x if x!=0 else (0. if y==0 else None)

def values(rows):
    return {tuple(r['modes']):r for r in rows}

def selected_pair(a,b,p,q):
    a,b=a['table'][tuple(sorted(p))],b['table'][tuple(sorted(q))]
    assert all(abs(v-w)<1e-14 for v,w in zip(a['momenta'],b['momenta']))
    av,bv=a['out']['kappa4'],b['out']['kappa4']
    return {'physical_p':a['momenta'],'kappa4_first':av,'kappa4_second':bv,
            'frac_delta_kappa4':rel(av,bv),
            'R2_first':a['out']['R2'],'R2_second':b['out']['R2']}

def full_pair_delta(a,b,fields=('N1_k','N1_q','N2','N2_G','kappa4')):
    rows=[(k,a['table'][k]['out'],b['table'][k]['out']) for k in a['table']]
    result={}
    for f in fields:
        rk,rx,ry=max(rows,key=lambda r:abs(rel(r[1][f],r[2][f])))
        ak,ax,ay=max(rows,key=lambda r:abs(r[1][f]-r[2][f]))
        result[f]={'max_rel':abs(rel(rx[f],ry[f])),'max_rel_modes':rk,
                   'max_rel_first':rx[f],'max_rel_second':ry[f],
                   'max_abs':abs(ax[f]-ay[f]),'max_abs_modes':ak}
    return result

def main():
    outputs={}
    for fn in Path('out').glob('*.json'):
        d=json.loads(fn.read_text()); c=d['input']; K=c['K']
        table=values(d['moments'])
        assert len(table) == len(d['moments']) == (2*K+1)*(K+1)
        assert d['method']['runner_sha256']==digest('run_out_modes.py')
        assert abs(d['diagnostics']['norm']-1)<1e-11
        assert d['diagnostics']['max_sampled_norm_error']<1e-11
        assert d['diagnostics']['hermiticity_error']<1e-15
        assert d['diagnostics']['max_radius_prob_Born_check']<4e-14
        num=gauss=conn=0.
        for pair,r in table.items():
            assert all(math.isfinite(s) for k in ('C_kq','A_kq') for s in r['out'][k])
            assert -K<=pair[0]<=pair[1]<=K
            out=r['out']
            assert all(out[f]>=0 for f in ('N1_k','N1_q','N2','N2_G'))
            assert out['kappa4']==out['N2']-out['N2_G']
            if pair[0]+pair[1]!=0:
                assert out['A_kq']==[0.,0.]
            if pair[0]!=pair[1]:
                assert out['C_kq']==[0.,0.]
            factor=1 if pair[0]==pair[1] else 2
            num+=factor*out['N2']; gauss+=factor*out['N2_G']; conn+=factor*out['kappa4']
        prob=d['diagnostics']['particle_number_probability']
        exp=sum(int(n)*(int(n)-1)*v for n,v in prob.items())
        assert abs(exp-num)<1e-13*abs(exp)+1e-30
        if not c['lambda0']:
            assert num==gauss==conn==0. and all(r['out']['R2'] is None for r in table.values())
        outputs[fn.stem]={'input':c,'table':table,'diag':d['diagnostics'],
                          'Nsummed':{'N2':num,'N2_G':gauss,'kappa4':conn,'sum_rule_error':abs(exp-num)}}
    base='k5-n6-L20-lam0.1'
    tests={'vacuum': 'k5-n6-L20-lam0' in outputs}
    comparisons={}
    a=outputs[base]
    for n,aa,bb in (
        ('K5_N6_to_N8_L20',base,'k5-n8-L20-lam0.1'),
        ('K6_N6_to_N8_L20','k6-n6-L20-lam0.1','k6-n8-L20-lam0.1'),
        ('K6_N6_to_N8_L24','k6-n6-L24-lam0.1','k6-n8-L24-lam0.1'),
        ('K5_N6_to_N7_even_parity',base,'k5-n7-L20-lam0.1'),
        ('K5_N6_loose_tight',base,base+'-tight'),
        ('fixed_L20_K5_to_K6_N6',base,'k6-n6-L20-lam0.1'),
        ('fixed_L20_K5_to_K6_N8','k5-n8-L20-lam0.1','k6-n8-L20-lam0.1'),
    ):
        if aa in outputs and bb in outputs:
            comparisons[n]=full_pair_delta(outputs[aa],outputs[bb])
    if 'K5_N6_to_N7_even_parity' in comparisons:
        assert all(x['max_abs']==0 for x in comparisons['K5_N6_to_N7_even_parity'].values())
        tests['N6_N7_same_sector']=True
    if 'K5_N6_loose_tight' in comparisons:
        assert comparisons['K5_N6_loose_tight']['kappa4']['max_rel']<2e-6
        tests['time_integration_tightening']=True
    for N in (6,8):
        aa,bb=[f'k{K}-n{N}-L{L}-lam0.1' for K,L in ((5,20),(6,24))]
        if aa in outputs and bb in outputs:
            x,y=outputs[aa],outputs[bb]
            assert abs(x['input']['cutoff_momentum']-y['input']['cutoff_momentum'])<1e-15
            volume=[selected_pair(x,y,p,q) for p,q in zip(
                        list(itertools.combinations_with_replacement((-5,0,5),2)),
                        list(itertools.combinations_with_replacement((-6,0,6),2)))]
            for v in volume:
                v['frac_delta_L_times_kappa4']=rel(20*v['kappa4_first'],24*v['kappa4_second'])
                v['frac_delta_L2_times_kappa4']=rel(20**2*v['kappa4_first'],24**2*v['kappa4_second'])
            comparisons[f'fixed_Lambda_N{N}']=volume
    weak={}
    for v in (0.1,0.025,0.0125):
        key=f'k5-n6-L20-lam{v}'
        if key in outputs:
            rows=outputs[key]['table']
            weak[str(v)]={'max_N2_relative_deviation_from_Born':max(
                abs(rel(r['leading_N2_mode_sum'],r['out']['N2'])) for r in rows.values()),
             'reference_pairs':{str(k):{'N2_over_lam2':rows[k]['out']['N2']/v**2,
                'Born_N2_over_lam2':rows[k]['leading_N2_mode_sum']/v**2,
                'kappa4_over_lam2':rows[k]['out']['kappa4']/v**2,
                'Gaussian_over_lam4':rows[k]['out']['N2_G']/v**4} for k in ((0,0),(-5,5),(1,2))}}
    if len(weak)==3:
        assert weak['0.0125']['max_N2_relative_deviation_from_Born'] < .6*weak['0.025']['max_N2_relative_deviation_from_Born']
        tests['weak_deviation_declines_linearly']=True
    result={'available_native_runs':len(outputs),'tests':tests,'per_run':{k:{c:v[c] for c in ('input','diag','Nsummed')} for k,v in sorted(outputs.items())},
       'pairs':{'fixed_L20_N6_examples': [selected_pair(a,outputs['k6-n6-L20-lam0.1'],k,k) for k in ((0,0),(-5,5),(1,2))]},
       'audit':comparisons,'weak_coupling':weak}
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=='__main__':
    main()
