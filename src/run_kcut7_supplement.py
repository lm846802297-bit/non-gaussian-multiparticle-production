#!/usr/bin/env python3
"""Supplemental Kcut=7 convergence point for phase2-dynamics-001."""
import json, os, platform, sys, time
from pathlib import Path
import numpy as np
import scipy
from real_time_fock_dynamics import evolve, mode_observables, leading_n2

ROOT=Path(__file__).resolve().parents[1]
res=evolve(7,6,1.0,20.0,0.2,1.0,rtol=5e-10,atol=5e-12)
obs=mode_observables(res['modes'],res['basis'],res['psi'],1,2)
lo=leading_n2(7,1.0,20.0,0.2,1.0,1,2)
out={
 "parameters":{"m":1.0,"L":20.0,"T":1.0,"lambda0_over_m2":0.2,"Kcut":7,"Nmax":6,"channel":"generic (n_k,n_q)=(1,2)","rtol":5e-10,"atol":5e-12,"time_interval":"[-12T,+12T]"},
 "basis_dimension":res['dimension'],"interaction_hermiticity_error":res['hermiticity_error'],
 "state_norm":res['norm'],"solver_nfev":res['nfev'],"elapsed_seconds":res['elapsed'],
 "N2_exact":obs['N2'],"N2_G":obs['N2_G'],"kappa4":obs['kappa4'],"R2":obs['R2'],
 "N2_LO":lo,"N2_over_LO":obs['N2']/lo,
 "runtime":{"python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__,"platform":platform.platform(),"machine":platform.machine(),"processor":platform.processor(),"cpu_count":os.cpu_count()}
}
p=ROOT/'results/phase2-dynamics-001/kcut7-supplement.json'
p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))