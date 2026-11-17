#!/usr/bin/env python3
"""901 E520 true-source retained lineage vs model selected Div -- not biology."""
from pathlib import Path
import json,csv,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/data/901_E520_20261117/results')
s=json.loads((P/'E520_original_source_solver_vs_CellStats_reconciliation.json').read_text())['comparisons']['Pos12_GL22']
assert s['pruned_main_solver_division_selected_only_n']==14 and s['cost_unbroken_Mapping_consecutive_ASC_count']==160
q=list(csv.DictReader((P/'E520_Pos12_17_original_extra_solver_DIV_vs_CellStats_export.csv').open()))
assert len(q)==17
fig,ax=plt.subplots(figsize=(11.4,4.0),dpi=160)
ax.hlines(3.1,243,607,linewidth=8,color='#26836c',label='Alternative CellStats: cell ID 76, f243–607 (ends exit)')
ax.hlines(2.2,246,606,linewidth=3,color='#26836c',alpha=.6,label='Alternative solver: 361 selected mappings')
ax.hlines(2.2,246,405,linewidth=7,color='#00747d',alpha=.93,label='160 consecutive forced ASC mappings (f246–405)')
ax.hlines(1.3,247,261,linewidth=8,color='#bf693b',label='Standard CellStats: cell ID 79, f247–261 (ends pruned)')
for kind,c,mark,label in [('True','#d9534f','s','Standard solver extra Div, omitted after PR: 14'),('False','#374c7a','o','Standard solver extra Div, exported: 3')]:
    hits=[int(r['frame']) for r in q if str(r['source_main_PR_root_descendant'])==kind]
    ax.scatter(hits,[.44]*len(hits),s=72,marker=mark,edgecolors='white',linewidths=.8,color=c,label=label,zorder=4)
ax.axvline(261,ls='--',lw=1,color='#d9534f')
ax.text(266,1.5,'PR at f261',fontsize=9,color='#b0473e')
ax.axvspan(369,405,color='#d9534f',alpha=.05)
ax.text(386,.75,'12 unexported selected Div\nin this short interval',fontsize=8,ha='left',color='#a3473a')
ax.set_yticks([.44,1.3,2.2,3.1])
ax.set_yticklabels(['Solver Div-only discrepancies','Standard exported track','Alternative solver mappings','Alternative exported track'],fontsize=9)
ax.set_ylim(-.12,4.25);ax.set_xlim(225,625)
ax.set_xlabel('Verified original MoMA zero-based frame')
ax.set_title('Pos12_GL22: archive-only selected divisions versus pruned CellStats lineage',fontsize=11)
ax.grid(axis='x',lw=.4,color='#d6d6d6')
ax.legend(loc='upper center',bbox_to_anchor=(.5,1.02),ncol=2,fontsize=7.5,frameon=False)
ax.spines[['top','right','left']].set_visible(False)
ax.text(.02,.03,'All curves: MoMA source-software decisions; neither independent phase-cell annotation nor biological fate truth',transform=ax.transAxes,fontsize=7.5,color='#6a6a6a')
fig.tight_layout();fig.savefig(P/'E520_Pos12_PR_f261_export_vs_solver_timeline.png',bbox_inches='tight');plt.close(fig)
