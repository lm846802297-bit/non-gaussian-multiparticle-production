"""Out-mode cumulants, archived finite-box Fock Hamiltonian and solver."""
import argparse, hashlib, itertools, json, math, os, platform, sys, time
from pathlib import Path

import numpy as np
import scipy
import real_time_fock_dynamics as f

def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--K',type=int,required=True)
    p.add_argument('--N',type=int,required=True)
    p.add_argument('--L',type=float,required=True)
    p.add_argument('--lam',type=float,required=True)
    p.add_argument('--step-scale',type=float,default=1.)
    p.add_argument('--rtol',type=float,default=2e-10)
    p.add_argument('--atol',type=float,default=2e-12)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); m=T=1.
    start=time.monotonic()
    e=f.evolve(a.K,a.N,m,a.L,a.lam,T,a.rtol,a.atol,save_times=49,max_step_scale=a.step_scale)
    print('evolved',a.K,a.N,a.L,a.lam,e['dimension'],time.monotonic()-start,flush=True)
    o=np.array(e['basis'],dtype=int)
    vac=int(e['index'][tuple(0 for _ in e['modes'])])
    pulse_fourier=np.array([f.lambda_tilde(E,a.lam,T) for E in e['energies']])
    amp1=e['V'][:,vac].toarray().flatten()*pulse_fourier
    prob1=amp1**2
    moments=[]
    for nk,nq in itertools.combinations_with_replacement(e['modes'],2):
        ik,iq=e['modes'].index(nk),e['modes'].index(nq)
        w=o[:,ik]*(o[:,iq]-int(nk==nq))
        analytic=f.leading_n2(a.K,m,a.L,a.lam,T,nk,nq)
        born_state=float(w@prob1)
        moments.append({'modes':[nk,nq],
               'momenta':[2*math.pi*nk/a.L,2*math.pi*nq/a.L],
               'out':f.mode_observables(e['modes'],e['basis'],e['psi'],nk,nq),
               'leading_N2_mode_sum':analytic,'leading_N2_Fock_Born':born_state})
    norms=np.sum(np.abs(e['states'])**2,axis=0)
    runtime=time.monotonic()-start
    data={'input':{'K':a.K,'Nmax':a.N,'L':a.L,'m':m,'T':T,'lambda0':a.lam,
                   'cutoff_momentum':2*math.pi*a.K/a.L,'step_scale':a.step_scale,
                   'rtol':a.rtol,'atol':a.atol,'max_step':e['max_step'],
                   'time_interval':[-e['tmax'],e['tmax']]},
          'method':{'module_file':Path(f.__file__).name,'module_sha256':digest(f.__file__),
                    'runner_sha256':digest(__file__),'solver':'archived DOP853',
                    'reference_commit':'eab1f5b3198d03eafc3f162299c9038956beaa0b',
                    'preparation':'free vacuum in P=0, even particle number sector',
                    'platform':platform.platform(),'python':platform.python_version(),
                    'numpy':np.__version__,'scipy':scipy.__version__,
                    'threads':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS')}},
          'moments':moments,
          'diagnostics':{'norm':e['norm'],'max_sampled_norm_error':float(np.max(np.abs(norms-1))),
             'hermiticity_error':e['hermiticity_error'],'dimension':e['dimension'],
             'nfev':e['nfev'],'elapsed_s':runtime,'time_sample_points':len(e['times']),
             'particle_number_probability':{str(k):float(sum(abs(e['psi'][o.sum(axis=1)==k])**2))
                                            for k in range(0,a.N+1,2)},
             'max_radius_prob_Born_check':max(abs(x['leading_N2_mode_sum']-x['leading_N2_Fock_Born'])/max(abs(x['leading_N2_mode_sum']),1e-300) for x in moments)}}
    a.out.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':
    main()
