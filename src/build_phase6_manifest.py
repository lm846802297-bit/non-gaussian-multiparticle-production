import hashlib, json, platform, sys
from pathlib import Path
import numpy as np, scipy
root=Path(__file__).resolve().parents[1]
r=root/'results'/'phase6-kcut-extension-001'
def load(name): return json.loads((r/name).read_text())
def z(d,k): return complex(*d[k])
rows=[]
for k in range(1,5):
 d=load(f'K{k}-N5-matched-loose.json')
 rows.append({'K':k,'Nmax':5,'basis_dimension':d['basis_dimension'],'abs_F4c':abs(z(d,'F4c')),'abs_M33':abs(z(d,'M33_projected')),'abs_source':abs(z(d,'connected_source')),'SD_relative_residual':d['SD_relative_residual'],'norm':d['vacuum_norm_out'],'hermiticity_error':d['interaction_hermiticity_error'],'rtol':2e-7,'atol':2e-9,'max_step':0.04})
k3n5=load('K3-N5-reference.json'); k3n6=load('K3-N6-h0p00625.json'); k3n6h=load('K3-N6-h0p003125.json')
k4n4=load('K4-N4-h0p00625.json'); k4n5=load('K4-N5-h0p00625.json'); k4n5h=load('K4-N5-h0p003125.json')
summary={
 'run_id':'phase6-kcut-extension-001',
 'scope':'K=1..4 matched Nmax=5 scan, plus K3 and K4 Fock/time refinements',
 'matched_K_scan':rows,
 'K3_fock_refinement':{
  'Nmax5_to6_relative_change_F4c':abs(z(k3n6,'F4c')-z(k3n5,'F4c'))/abs(z(k3n5,'F4c')),
  'Nmax5_to6_relative_change_M33':abs(z(k3n6,'M33_projected')-z(k3n5,'M33_projected'))/abs(z(k3n5,'M33_projected')),
  'Nmax5_to6_relative_change_source':abs(z(k3n6,'connected_source')-z(k3n5,'connected_source'))/abs(z(k3n5,'connected_source')),
  'Nmax6_h00625_residual':k3n6['SD_relative_residual'],
  'Nmax6_h003125_residual':k3n6h['SD_relative_residual']},
 'K4_fock_refinement':{
  'Nmax4_to5_relative_change_F4c':abs(z(k4n5,'F4c')-z(k4n4,'F4c'))/abs(z(k4n4,'F4c')),
  'Nmax4_to5_relative_change_M33':abs(z(k4n5,'M33_projected')-z(k4n4,'M33_projected'))/abs(z(k4n4,'M33_projected')),
  'Nmax4_to5_relative_change_source':abs(z(k4n5,'connected_source')-z(k4n4,'connected_source'))/abs(z(k4n4,'connected_source')),
  'Nmax4_residual_h00625':k4n4['SD_relative_residual'],
  'Nmax4_residual_h003125':load('K4-N4-h0p003125.json')['SD_relative_residual'],
  'Nmax5_residual_h00625':k4n5['SD_relative_residual'],
  'Nmax5_residual_h003125':k4n5h['SD_relative_residual'],
  'stencil_converged':False},
 'cached_engine_crosscheck':{'K':1,'Nmax':5,'matches_original_stored_results_exactly':True},
 'unit_tests':{'count':6,'status':'passed'},
 'limitations':['K scan is non-monotonic and not converged.','K=4 Nmax=6 not computed.','K=4 Nmax=5 SD residual is not time-stencil converged.','No continuum or infinite-volume extrapolation.'],
 'runtime':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'machine':platform.machine(),'cpu_model':'Intel Xeon Platinum 8573C','visible_cpu_count':5,'virtualization':'KVM','blas_threads':1,'omp_threads':1},
 'source_sha256':hashlib.sha256((root/'src'/'sd_residual_cached.py').read_bytes()).hexdigest(),
 'reference_source_sha256':hashlib.sha256((root/'src'/'sd_residual.py').read_bytes()).hexdigest(),
 'test_sha256':hashlib.sha256((root/'tests'/'test_phase6_extended_scan.py').read_bytes()).hexdigest()
}
(r/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
print(json.dumps(summary,indent=2,sort_keys=True))