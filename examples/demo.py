"""Run: python examples/demo.py  (after pip install -e .)"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from cstr_convergence import first_order_cstr, first_order_pfr, first_order_minimum_stages, series_pfr, series_cstr

out = Path('outputs')
out.mkdir(exist_ok=True)
ns = np.arange(1, 101)
fig, ax = plt.subplots(figsize=(10, 8))
for da, color in [(0.5, '#800020'), (2., '#002FA7'), (5., '#00754A')]:
    errors = [abs(first_order_pfr(da)-first_order_cstr(da, int(n))) for n in ns]
    ax.plot(ns, errors, lw=3, color=color, label=f'Da = {da:g}')
ax.set(xlabel='Number of CSTRs, N', ylabel='Absolute conversion error')
ax.set_yscale('log')
for spine in ax.spines.values(): spine.set_linewidth(2)
ax.tick_params(width=2, labelsize=14)
ax.xaxis.label.set_size(19); ax.yaxis.label.set_size(19)
ax.xaxis.label.set_weight('bold'); ax.yaxis.label.set_weight('bold')
for tick in ax.get_xticklabels()+ax.get_yticklabels(): tick.set_fontweight('bold')
ax.legend(frameon=False, fontsize=14)
fig.tight_layout()
fig.savefig(out/'convergence.png', dpi=200, facecolor='white')
plt.close(fig)

print('Da = 2, absolute conversion tolerance = 0.01')
print('Minimum CSTR count:', first_order_minimum_stages(2., 0.01))
print('A -> B -> C (Da1=2, Da2=1):')
print('PFR:', series_pfr(2., 1.))
print('10 CSTRs:', series_cstr(2., 1., 10))
print('Saved outputs/convergence.png')
