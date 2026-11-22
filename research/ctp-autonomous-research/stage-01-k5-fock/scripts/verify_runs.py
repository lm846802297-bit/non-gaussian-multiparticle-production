"""Validate run shape, source identity, and CTP symmetry; report physical metrics."""
import hashlib, json, math,sys
from pathlib import Path

def val(x,q): return complex(*x[q])
def flip(mask): return ''.join('+' if a=='-' else '-' for a in mask)
labels={f'{i:04b}'.translate(str.maketrans('01','-+')) for i in range(16)}
results={}
for p in sorted(Path("run").glob("k5-*.json")):
 d=json.loads(p.read_text())
 assert d['source_commit']=='eab1f5b3198d03eafc3f162299c9038956beaa0b'
 runner=Path('scan_segmented.py' if ('segmented' in p.name or d['Nmax']==7) else 'scan.py')
 assert hashlib.sha256(runner.read_bytes()).hexdigest()==d['runner_sha256']
 assert d['basis_dimension']==math.comb(11+d['Nmax'],11)
 assert [r['h'] for r in d['rows']]==[.00625,.003125,.0015625]
 assert all(set(r['components'])==labels for r in d['rows'])
 center=d['rows'][1]['components']
 for r in d['rows']:
  c=r['components']
  assert all(c[m]['F4c']==center[m]['F4c'] and c[m]['M33']==center[m]['M33'] for m in labels)
  assert all(math.isfinite(c[m]['SD_relative_residual']) for m in labels)
  for q in ['F4c','M33','connected_composite_source','SD_residual','analytic_residual']:
   assert max(abs(val(c[m],q)-val(c[flip(m)],q).conjugate()) for m in labels)<2e-12
   assert max(abs(val(c[m],q)-val(c[flip(m[0])+m[1:]],q)) for m in labels)<(4e-16/r['h']**2 if q=='SD_residual' else 2e-12)
 a=max(center,key=lambda m:center[m]['analytic_relative_residual'])
 rich=max(abs((4*val(d['rows'][2]['components'][m],'SD_residual')-val(center[m],'SD_residual'))/3-val(center[m],'analytic_residual')) for m in labels)
 results[p.stem]={'basis':d['basis_dimension'],'max_rel': [r['max_relative'] for r in d['rows']], 'max_fd_error':[r['max_fd_error'] for r in d['rows']], 'max_analytic_relative':center[a]['analytic_relative_residual'], 'max_analytic_branch':a, 'max_analytic_abs':d['max_analytic_abs'], 'max_richardson_disagreement':rich,'++++':{q:center['++++'][q] for q in ['F4c','M33','connected_composite_source']}}
if 'k5-n7' in results and 'k5-n6-segmented' in results:
 a,b=(json.loads(Path('run/'+n+'.json').read_text())['rows'][1]['components'] for n in ('k5-n6-segmented','k5-n7'))
 results['n6_to_n7']={}
 for q in ['F4c','M33','connected_composite_source']:
  changes={m:abs(val(b[m],q)-val(a[m],q))/abs(val(a[m],q)) for m in labels}
  results['n6_to_n7'][q]={'maximum':max(changes.values()),'++++':changes['++++']}
print(json.dumps(results,indent=2))
