#!/usr/bin/env python3
"""Make non-inferential original-software cost margin rank plot; no biology labels."""
import csv,pathlib
from decimal import Decimal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=pathlib.Path('/mnt/data/901_E526_20261117/results')
rows=list(csv.DictReader((P/'E526_151_local_moves_exact_integer_ranked.csv').open()))
fig,ax=plt.subplots(figsize=(9,5.2),dpi=160)
for gl in ['Pos5_GL35','Pos12_GL22']:
 vals=sorted(Decimal(r['delta_exact_decimal']) for r in rows if r['GL']==gl)
 ax.plot(range(1,len(vals)+1),[float(v) for v in vals],marker='.',markersize=3,linewidth=1.5,label=f'{gl}: {len(vals)} feasible 5-action shifts')
ax.set_yscale('log')
ax.set_xlabel('Within-growthlane rank among feasible one-frame five-action swaps')
ax.set_ylabel('Exact archived software cost increase (log scale)')
ax.set_title('E526: source-price exact cost margins for 151 local timing alternatives')
ax.legend(loc='upper left')
ax.grid(True,alpha=0.2)
fig.tight_layout()
fig.savefig(P/'E526_151_exact_local_swaps_cost_rank.png')
print(P/'E526_151_exact_local_swaps_cost_rank.png')
