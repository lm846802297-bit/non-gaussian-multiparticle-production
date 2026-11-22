"""Validate complete K5,N7 CF4 step refinement and DOP853 agreement outputs."""
from pathlib import Path
import json, hashlib, math

p=Path('cf4_k5_n7_A.json')
j=json.loads(p.read_text())
assert j['basis_dimension']==31824 and j['K']==5 and j['Nmax']==7
assert j['source_commit']=='eab1f5b3198d03eafc3f162299c9038956beaa0b'
assert j['runner_sha256']==hashlib.sha256(Path('cf4.py').read_bytes()).hexdigest()
assert j['baseline_sha256']==hashlib.sha256(Path('/tmp/ctpproject/run/k5-n7.json').read_bytes()).hexdigest()
assert j['source_module_sha256']==hashlib.sha256(Path('/tmp/ctpproject/src/sd_residual.py').read_bytes()).hexdigest()
assert j['wick_module_sha256']==hashlib.sha256(Path('/tmp/ctpproject/independent_source.py').read_bytes()).hexdigest()
assert [x['max_step'] for x in j['steps']]==[.1,.05]
coarse,fine=j['steps']
for row in (coarse,fine):
    assert abs(row['vacuum_norm']-1)<5e-13
    assert row['h']==.003125
    assert row['comparison_DOP853']['analytic_residual']['absolute']<5e-13
    assert row['comparison_DOP853']['SD_residual']['absolute']<1e-10
    for key in ('F4c','M33','connected_composite_source'):
        assert row['comparison_DOP853'][key]['relative_to_reference']<1e-6
for key in ('F4c','M33','connected_composite_source','analytic_residual','SD_residual'):
    a=coarse['comparison_DOP853'][key];b=fine['comparison_DOP853'][key]
    d=fine['comparison_previous_CF4'][key]
    assert all(math.isfinite(x['relative_to_reference']) for x in (a,b,d))
    print(f"{key:27} vs DOP853 coarse {a['absolute']:.8g} ({a['relative_to_reference']:.8g}) fine {b['absolute']:.8g} ({b['relative_to_reference']:.8g}); |coarse - fine| {d['absolute']:.8g}")
print("PASS: identity, hashes, both steps, unitarity, and DOP853 comparison")
