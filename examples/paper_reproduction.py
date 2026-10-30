"""Published SI Table S10 concentration dependence: run with python examples/paper_reproduction.py."""
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from cstr_convergence import XyloseModel

model=XyloseModel()
concs=np.linspace(0.033,0.34,27)
phi=np.array([model.phi(float(c)) for c in concs])
out=Path('results');out.mkdir(exist_ok=True)
with (out/'xylose_phi.csv').open('w',newline='') as fh:
    wr=csv.writer(fh);wr.writerow(['C_X0_M','Phi']);wr.writerows(zip(concs,phi))
fig,ax=plt.subplots(figsize=(10,8));ax.plot(concs,phi,'o-',color='#800020',lw=3,markersize=8)
ax.set_xlabel('Inlet xylose concentration (mol/L)',fontsize=17,fontweight='bold')
ax.set_ylabel('Optimized finite-stage coefficient Φ',fontsize=17,fontweight='bold')
ax.tick_params(labelsize=13,width=2)
for spine in ax.spines.values():spine.set_linewidth(2)
fig.tight_layout();fig.savefig(out/'xylose_phi.png',dpi=200)
print('Saved results/xylose_phi.csv and results/xylose_phi.png')
for c in [0.033,0.176,0.340]:
    print(f'C={c:.3f} M | Phi={model.phi(c):.4f} | dPhi/dln(gamma)={model.phi_sensitivity(c):.4f}')
