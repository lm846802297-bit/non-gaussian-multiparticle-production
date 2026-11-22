import json, os
import numpy as np
import matplotlib.pyplot as plt
base='/mnt/data/ngmp_phase4/results/phase4-sd-001'
hs=np.array([0.025,0.0125,0.00625,0.003125])
res=[]
for h in hs:
    d=json.load(open(os.path.join(base,f'spectral-n8-h{h}-a0.2.json')))
    res.append(d['SD_relative_residual'])
plt.figure()
plt.loglog(hs,res,'o-',label='SD relative residual')
plt.loglog(hs,res[-1]*(hs/hs[-1])**2,'--',label='reference $h^2$')
plt.xlabel('time stencil spacing h')
plt.ylabel('relative Schwinger-Dyson residual')
plt.legend()
plt.tight_layout()
plt.savefig('/mnt/data/ngmp_phase4/figures/phase4_sd_residual_vs_h.svg')
plt.close()

as_ = np.array([0.05,0.1,0.2])
errs=[]
for a in as_:
    d=json.load(open(os.path.join(base,f'spectral-n8-h0.00625-a{a}.json')))
    fd=np.array(d['dxx_F4c_fd']); sp=np.array(d['dxx_F4c_spectral'])
    errs.append(np.linalg.norm(fd-sp)/np.linalg.norm(sp))
plt.figure()
plt.loglog(as_,errs,'o-',label='finite-difference vs spectral derivative')
plt.xlabel('space stencil spacing a')
plt.ylabel('relative spatial-derivative error')
plt.legend()
plt.tight_layout()
plt.savefig('/mnt/data/ngmp_phase4/figures/phase4_spatial_stencil_error.svg')
plt.close()