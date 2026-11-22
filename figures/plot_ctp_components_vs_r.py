import json
from pathlib import Path
import matplotlib.pyplot as plt

base=Path('results/phase3-ctp-all-001')
files=[base/'spatial-r1-Nmax16.json',base/'spatial-r2-Nmax16.json',base/'ctp-components-Nmax16.json']
rs=[1.0,2.0,3.0]
labels=['++++','+++-','++-+','++--']
values={lab:[] for lab in labels}
for path in files:
    d=json.load(open(path))
    for lab in labels:
        z=complex(*d['components'][lab]['F4_connected'])
        values[lab].append(z)
fig,axes=plt.subplots(2,1,figsize=(7.2,7.2),sharex=True)
for lab in labels:
    axes[0].plot(rs,[z.real for z in values[lab]],marker='o',label=lab)
    axes[1].plot(rs,[z.imag for z in values[lab]],marker='o',label=lab)
axes[0].set_ylabel(r'Re $F_{4,c}^{a_1a_2a_3a_4}$')
axes[1].set_ylabel(r'Im $F_{4,c}^{a_1a_2a_3a_4}$')
axes[1].set_xlabel(r'Spatial extent parameter $r$')
axes[0].set_title(r'All-branch connected four-point components ($t_i=(1,0.5,-0.5,-1)$)')
for ax in axes:
    ax.grid(True,alpha=.25)
    ax.legend(ncol=4,fontsize=8)
fig.tight_layout()
fig.savefig('phase3-ctp-components-vs-r.svg')
fig.savefig('phase3-ctp-components-vs-r.png',dpi=180)
