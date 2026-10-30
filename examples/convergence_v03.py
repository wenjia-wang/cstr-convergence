"""Run: python examples/convergence_v03.py (plot extra required)."""
import numpy as np
import matplotlib.pyplot as plt
from cstr_convergence import series_pfr, series_cstr_exact, series_asymptotic_error

ns = np.arange(1, 401)
a,b = 6.,.4
pfr=series_pfr(a,b)[1]
errors=np.array([series_cstr_exact(a,b,int(n))[1]-pfr for n in ns])
fig,ax=plt.subplots(figsize=(10,8))
ax.plot(ns,errors,color='#800020',lw=3,label='Exact finite CSTR')
ax.plot(ns,series_asymptotic_error(a,b)[1]/ns,'--',color='#002FA7',lw=3,label='Leading asymptotic term')
ax.axhline(0,color='black',lw=1.5)
ax.set_xlabel('CSTR stages, N',fontweight='bold',fontsize=17)
ax.set_ylabel('Intermediate B yield error',fontweight='bold',fontsize=17)
ax.tick_params(labelsize=14,width=2)
for tick in ax.get_xticklabels()+ax.get_yticklabels(): tick.set_fontweight('bold')
for spine in ax.spines.values():spine.set_linewidth(2)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig('series_b_convergence.png',dpi=220,facecolor='white')
