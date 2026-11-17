#!/usr/bin/env python3
"""Plot archived source ASC contradictions only; no biological labels or p-values."""
import argparse,csv,pathlib,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
A=argparse.ArgumentParser();A.add_argument('--input',required=True);A.add_argument('--output',required=True);a=A.parse_args()
rows=list(csv.DictReader(open(a.input,encoding='utf8')))
fig,ax=plt.subplots(figsize=(10.5,3.9),dpi=160)
for i,gl in enumerate(['Pos5_GL35','Pos12_GL22']):
 sub=[r for r in rows if r['gl']==gl]
 for class_bool,marker,lab,color in [('True','o','Direct disputed-mother ASC','#246a9f'),('False','s','Other forced ASC','#d98936')]:
  arr=[r for r in sub if str(r['directly_same_mother_dispute'])==class_bool]
  ax.scatter([int(r['ASC_source_frame_zero_based']) for r in arr],[1-i]*len(arr),s=95,marker=marker,c=color,
             edgecolors='white',linewidths=.75,alpha=.92,label=lab if i==0 else None,zorder=3)
 ax.plot([0,650],[1-i,1-i],lw=.6,c='#a5a5a5',alpha=.7,zorder=1)
ax.axvspan(369,405,color='#d98936',alpha=.08,zorder=0)
ax.text(378,1.2,'Pos12: 12 other constraints\ncluster at frames 369–405',fontsize=8,c='#60422e',ha='left')
ax.set_yticks([0,1]);ax.set_yticklabels(['Pos12_GL22','Pos5_GL35']);ax.set_xlim(0,660);ax.set_ylim(-.5,1.65)
ax.set_xlabel('Original source (zero-based) time frame')
ax.set_title('25 archived ASC equalities incompatible with the fixed standard MoMA solution',fontsize=11)
ax.legend(loc='upper left',fontsize=8,frameon=False,ncol=2)
ax.spines[['top','right','left']].set_visible(False)
ax.grid(axis='x',color='#ddd',lw=.5)
ax.text(.01,.02,'Source: E512 official ZIP64 MoMA curation + Gurobi .sol (not biological division truth)',transform=ax.transAxes,fontsize=7,color='#777')
fig.tight_layout();fig.savefig(a.output,bbox_inches='tight');plt.close(fig)
