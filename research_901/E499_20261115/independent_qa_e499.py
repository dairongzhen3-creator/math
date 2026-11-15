#!/usr/bin/env python3
"""E499 independent scientific QA: pure stdlib read from source-authenticated E498
bleaching gzip, independently reconstruct every per-cell weighted log-intensity
slope, cadence, cross-phase matched IDs and source grouped diagnostics.
Does not import primary E499 code, pandas/numpy or other source-derived results.
"""
from pathlib import Path
import gzip,csv,json,hashlib,statistics,math
from collections import defaultdict,Counter

P=Path(__file__).resolve().parent;B=P.parent/'p901_E498' if (P.parent/'p901_E498').exists() else P/'source_E498'
LOCK={
'myframes_bleaching.feather.csv.gz':'0a9b630b95af875ae6bba5be5bac5430acbf96cfcf4bc767576965e59f8251db',
'E498_bleaching_per_cell_weighted_slopes.csv':'3e2a7c5a620d767ad9f0f2192f72008d52af90e0d0fd30d3eac188c0f0d0b45a'}
for n,h in LOCK.items():assert hashlib.sha256((B/n).read_bytes()).hexdigest()==h
cells=defaultdict(list)
unique_keys=set();nrows=0
with gzip.open(B/'myframes_bleaching.feather.csv.gz','rt',newline='') as fh:
 for x in csv.DictReader(fh):
  key=(x['date'],x['strain'],x['Fluo_sig'],int(x['exp_series']),int(x['idx']),x['cell'])
  datum=(int(x['time_sec']),float(x['signal']))
  assert datum[1]>0
  assert (key,datum[0]) not in unique_keys
  unique_keys.add((key,datum[0]))
  cells[key].append(datum);nrows+=1
assert nrows==91838
assert all((20 if key[-2]==1 else 5)==(20-15*int(key[-2]==2)) for key in cells)
cad=defaultdict(lambda:Counter())
for key,data in cells.items():
 date,strain,fluo,series,phase,cell=key
 ordered=sorted(data)
 assert len({t for t,y in ordered})==len(ordered)
 gaps=[(b[0]-a[0])/60 for a,b in zip(ordered[:-1],ordered[1:])]
 cad[(date,strain,fluo,series,phase)].update(gaps)
assert len(cad)==14
for key,gaps in cad.items():
 mode=gaps.most_common(1)[0][0]
 assert mode==(3 if key[-1]==1 else 12)

# Independently recalculate published source E498 WLS using positive per frame input.
raw={}
for key,obs in cells.items():
 if len(obs)<6:continue
 times=[t/3600 for t,v in obs];logs=[math.log(v) for t,v in obs];weights=[1/v for t,v in obs]
 W=sum(weights)
 mt=sum(w*t for w,t in zip(weights,times))/W
 my=sum(w*y for w,y in zip(weights,logs))/W
 vt=sum(w*(t-mt)**2 for w,t in zip(weights,times))/W
 cv=sum(w*(t-mt)*(y-my) for w,t,y in zip(weights,times,logs))/W
 if vt<=0:continue
 raw[(key[0],key[3],key[4],key[5],key[2])]=-cv/vt
assert len(raw)==1659
source={}
with (B/'E498_bleaching_per_cell_weighted_slopes.csv').open(newline='') as fh:
 for x in csv.DictReader(fh):
  key=(x['date'],int(x['experiment_series']),int(x['phase_idx']),x['cell'],x['Fluo_sig'])
  assert key not in source
  source[key]=float(x['raw_degradation_corrected_signal_rate_per_h'])
assert set(raw)==set(source)
assert all(abs(value-source[key])<1e-10 for key,value in raw.items())

# Separate native grouping of source signals, no primary outputs used.
by=defaultdict(dict);allgroup=defaultdict(list)
for (date,series,phase,cell,fluo),k in raw.items():
 by[(date,fluo,series,cell)][phase]=k
 allgroup[(date,fluo,phase)].append(k)
matched=[(k,v[1],v[2]) for k,v in by.items() if 1 in v and 2 in v]
assert len(matched)==16
r=[]
for (date,fluo,phase),values in sorted(allgroup.items()):
 v=[x for x in values if x>0]
 r.append((date,fluo,phase,statistics.median(v),len(values),len(v)))
assert len(r)==14
n_groups=0;rat=[]
for date,fluo in sorted({(t[0],t[1]) for t in r}):
 p1=next(t for t in r if t[:3]==(date,fluo,1));p2=next(t for t in r if t[:3]==(date,fluo,2))
 ratio=p2[3]/p1[3];n_groups+=1;rat.append(ratio)
 assert ratio<.25

e499=json.loads((P/'E499_results.json').read_text())
assert abs(statistics.median(rat)-e499['new_positive_only_comparison']['group_hourly_starve_over_growth_ratio_min_median_max'][1])<1e-13
assert sum(x>y for k,x,y in matched)==e499['same_exact_cell_cross_phase']['growth_hourly_greater']==14
assert sum(x/20>y/5 for k,x,y in matched)==e499['same_exact_cell_cross_phase']['growth_per_frame_greater']==13
assert sum(x>0 and y>0 for k,x,y in matched)==9
n_grows=sum(sum(v.values()) for k,v in cad.items() if k[-1]==1)
growth_exact=sum(v[3] for k,v in cad.items() if k[-1]==1)
n_starves=sum(sum(v.values()) for k,v in cad.items() if k[-1]==2)
starve_exact=sum(v[12] for k,v in cad.items() if k[-1]==2)
assert abs(growth_exact/n_grows-e499['source_cadence']['fraction_native_gap_growth'])<1e-14
assert abs(starve_exact/n_starves-e499['source_cadence']['fraction_native_gap_starvation'])<1e-14

qa={'project':901,'stage':'E499','PASS':True,'independent_source_gzip_rows':nrows,
 'no_duplicate_true_source_cell_time':True,
 'raw_cell_phase_groups':len(cells),'wls_cell_phase_n':len(raw),
 'every_source_E498_per_cell_WLS_recomputed_from_raw_rows':True,
 'source_nominal_gap_phase_groups':len(cad),
 'source_exact_3min_gaps':growth_exact,'source_all_growth_gaps':n_grows,
 'source_exact_12min_gaps':starve_exact,'source_all_starvation_gaps':n_starves,
 'new_observation_group_pairs':n_groups,
 'positive_only_group_ratio_below_quarter_n':sum(a<.25 for a in rat),
 'group_ratio_median':statistics.median(rat),
 'same_cell_both_phases_with_6plus_frames':len(matched),
 'same_cell_greater_loss_hour_n':sum(x>y for k,x,y in matched),
 'same_cell_greater_loss_per_planned_frame_n':sum(x/20>y/5 for k,x,y in matched),
 'measurement_exposure_phase_collinearity':'rank-deficient; no same-phase cross-frequency treatment',
 'source_scope':'Independent numerical source repeat on same original experimental data, not an independent camera dose or absolute molecules'}
(P/'E499_independent_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2))
print(json.dumps(qa,ensure_ascii=False,indent=2))
