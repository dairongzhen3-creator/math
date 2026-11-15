#!/usr/bin/env python3
"""E491 numerical QA independent of the main analysis functions.
Original RData bytes were verified against published Git blob by a separate cloud R runner;
this QA parses its exported source matrices using Python stdlib and checks raw source IDs.
The source decoder and processed measurements are shared; no independent wet-lab repeat.
"""
from pathlib import Path
import csv,math,statistics,json,hashlib
from collections import defaultdict,Counter
ROOT=Path(__file__).resolve().parent
meta={}
with (ROOT/'metadata_original.csv').open(newline='') as f:
 for r in csv.DictReader(f):meta[r['cell_id']]={'mother':r['mother_id'],'birth':int(r['cell_birth_timepoint'])}
assert len(meta)==85 and sum(x['mother']=='root' for x in meta.values())==25
prots=['Cdc10','Stb3','CLB5','Whi5','Xbp1','Tup1']
values=defaultdict(dict)
with (ROOT/'time_series_original_long.csv').open(newline='') as f:
 for r in csv.DictReader(f):
  v=r['source_value']
  if v:values[(r['cell_id'],r['protein'])][int(r['time_index'])]=float(v)
assert len(values)==85*6
parents={x['mother'] for x in meta.values() if x['mother']!='root'}
rows=[]
for child,record in meta.items():
 if record['mother']=='root' or record['birth']<=5:continue
 b=record['birth'];parent=record['mother'];prefix=child.split('-')[0]
 possible=[q for q in parents if q.split('-')[0]==prefix and q!=child and meta[q]['birth']<=b-1 and all(b-1 in values[(q,p)] for p in prots)]
 assert parent in possible and len(possible)>=2
 for p in prots:
  target=values[child,p][b]
  dist=abs(values[parent,p][b]-target)
  competitors=[abs(values[q,p][b]-target) for q in possible if q!=parent]
  rank=sum(v>dist for v in competitors)/len(competitors)
  rows.append({'child':child,'p':p,'rank':rank,'strict_closest':dist<min(competitors),'candidate_n':len(possible)})
assert len(rows)==48*6
by={p:[r for r in rows if r['p']==p] for p in prots}
identity_results={p:{'n':len(s),'rank_mean':sum(r['rank'] for r in s)/len(s),'strict_best_n':sum(bool(r['strict_closest']) for r in s)} for p,s in by.items()}
source_results=json.loads((ROOT/'E491_parent_identity_summary.json').read_text())['by_alignment']['coincident_birth_interpolated']
assert all(abs(identity_results[p]['rank_mean']-source_results[p]['mean_fraction_impostors_farther'])<1e-13 for p in prots)
assert all(identity_results[p]['strict_best_n']==round(source_results[p]['matched_strict_closest_fraction']*48) for p in prots)

x=[];date_sample=defaultdict(list);schema=Counter();units=set()
with (ROOT/'E491_frozen_division_sum_audit.csv').open(newline='') as f:
 for r in csv.DictReader(f):
  parent=float(r['parent_last_amp']);daughters=float(r['daughter_sum_amp']);
  r0=daughters/parent
  assert abs(r0-float(r['fluorescence_sum_ratio']))<1e-12
  assert int(r['daughter_birth_frame'])==int(r['parent_end_frame'])+1
  x.append(r0);date_sample[r['date']].append(r0);schema[r['schema']]+=1;units.add(r['lane_file'])
assert len(x)==42658 and len(date_sample)==27 and len(units)==857
ratio_median=statistics.median(x)
within=sum(.95<=v<=1.05 for v in x)/len(x)
abovel=sum(v>1 for v in x)/len(x)
basis=json.loads((ROOT/'E491_fluorescence_ratio_summary.json').read_text())['schema']['all']
assert abs(ratio_median-basis['fluorescence_ratio_median'])<1e-12
assert abs(within-basis['fraction_within_5pct_of_unity'])<1e-12
assert abs(abovel-basis['fraction_daughters_sum_greater_than_parent'])<1e-12
r={
 'stage':'E491','QA_pass':True,
 'identity_source_type':'published RData verified Git blob + its author-structure CSV export',
 'independent_qa_method':'Python stdlib CSV parser and standalone loops; no imported main analysis functions',
 'source_cells':len(meta),'founders':sum(x['mother']=='root' for x in meta.values()),
 'parents_known_in_cell_table':all(p in meta for p in parents),
 'eligible_daughters':48,'measured_proteins':prots,
 'coincident_recomputed_identity_metrics':identity_results,
 'frozen_bacterial_events':len(x),'frozen_bacterial_dates':len(date_sample),'frozen_lanes':len(units),
 'frozen_bacterial_median_ratio':ratio_median,'frozen_bacterial_fraction_within5pct':within,
 'frozen_bacterial_fraction_daughter_sum_gt_mother':abovel,
 'source_scope':'SOURCE_CODE_NUMERICAL_QA_ONLY; NOT photon technical repeat; NOT molecular count calibration',
 'molecular_information_transferred':'NOT_IDENTIFIED'
}
(ROOT/'E491_independent_QA.json').write_text(json.dumps(r,ensure_ascii=False,indent=2))
print(json.dumps({k:v for k,v in r.items() if k not in ('coincident_recomputed_identity_metrics',)},ensure_ascii=False,indent=2))
