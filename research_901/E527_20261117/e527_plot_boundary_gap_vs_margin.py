#!/usr/bin/env python3
"""Visual diagnostic of source-component gap pixels vs archived-cost margins."""
import csv,pathlib
from decimal import Decimal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
base=pathlib.Path('/mnt/data/901_E527_20261117/results')
rows=list(csv.DictReader((base/'E527_151_exact_source_cost_component_geometry_decomposition.csv').open()))
fig,ax=plt.subplots(figsize=(8.5,5.2))
for gl,marker in [('Pos5_GL35','o'),('Pos12_GL22','x')]:
    rr=[r for r in rows if r['GL']==gl and Decimal(r['delta_exact_source_cost'])<=Decimal('0.001')]
    ax.scatter([int(r['parent_missing_in_daughters_px']) for r in rr],[float(r['delta_exact_source_cost']) for r in rr],marker=marker,label=f'{gl}, n={len(rr)}',alpha=0.75)
ax.set_yscale('log')
ax.set_xlabel('Pixel gap: whole mask minus union of two daughter masks (pixels)')
ax.set_ylabel('Exact additional archived assignment cost (log scale)')
ax.set_title('901 E527: 55 near-optimal 1-frame software-division alternatives')
ax.grid(True,alpha=0.2)
ax.legend(loc='upper left')
for gl,fr,lab in [('Pos5_GL35','960','P5 f960'),('Pos12_GL22','118','P12 f118')]:
    rr=[r for r in rows if r['GL']==gl and r['source_division_t']==fr and r['direction']=='early']
    assert len(rr)==1
    r=rr[0]
    ax.annotate(lab,xy=(int(r['parent_missing_in_daughters_px']),float(r['delta_exact_source_cost'])),xytext=(7,8),textcoords='offset points',fontsize=8)
fig.tight_layout()
p=base/'E527_55_near_optimal_gap_pixels_vs_exact_cost.png';fig.savefig(p,dpi=150)
plt.close(fig)
print('E527_PLOT_PASS',p,p.stat().st_size)
