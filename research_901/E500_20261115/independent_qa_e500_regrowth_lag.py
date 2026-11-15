#!/usr/bin/env python3
"""Independent E500 source regrowth-lag QA from 71,755 exact source-exported rows.
Pure stdlib implementation, no pandas, NumPy or reuse of first-cloud-model functions.
The original Feather's exact SHA is authenticated but this alternative decoder uses a source-export CSV.gz.
"""
from pathlib import Path
import hashlib,csv,gzip,json,math,collections,statistics
P=Path(__file__).resolve().parent
orig=P/'myframes_control.feather'
assert orig.stat().st_size==302295912
h=hashlib.sha256()
with orig.open('rb') as f:
 for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
assert h.hexdigest()=='5cb1850cb127a64bf5e3dffa4d4e13e7cc648f73af176f6d03cf4d0e1fa9bb7c'
rows=collections.defaultdict(list);seen=set();N=0
with gzip.open(P/'E500_238_true_mothers_full_GM2_timeseries.csv.gz','rt',newline='') as fh:
 for x in csv.DictReader(fh):
  N+=1
  assert x['c_step']=='GM.2' and int(float(x['b_rank']))==0
  assert float(x['time_start'])<float(x['m_start'])
  key=(x['condition'],x['cell']);t=float(x['time_sec'])
  assert (key,t) not in seen;seen.add((key,t))
  rows[key].append((t,float(x['length_um']),x))
assert N==71755 and len(rows)==238
original={}
with (P/'E500_238_mothers_regrowth_lag_author_11frame_qualified.csv').open(newline='') as fh:
 for r in csv.DictReader(fh):
  key=(r['condition'],r['cell']);assert key not in original;original[key]=r
assert set(original)==set(rows)

def fit(vals):
 n=len(vals);X=[t for t,l,rr in vals[-11:]];Y=[math.log(l) for t,l,rr in vals[-11:]]
 mx=sum(X)/11;my=sum(Y)/11
 slope=sum((t-mx)*(y-my) for t,y in zip(X,Y))/sum((t-mx)**2 for t in X)
 intercept=my-slope*mx
 ss=sum((y-my)**2 for y in Y)
 r2=1-sum((y-(intercept+slope*t))**2 for t,y in zip(X,Y))/ss if ss>0 else 0
 l0=sum(l for t,l,rr in vals[:10])/10
 return slope,intercept,r2,l0

passed=collections.defaultdict(list);reasons=collections.Counter();near=0
for key,vals in rows.items():
 vals.sort(key=lambda v:v[0]);d=original[key]
 assert len(vals)==int(d['n_raw_frames'])
 assert vals[0][-1]['end_type']==d['end_type']
 if len(vals)<=10:reason='N_LE_10'
 elif min(v[1] for v in vals)<=0:reason='NONPOSITIVE_LENGTH'
 else:
  b,a,r2,l0=fit(vals)
  for field,v in [('source_last11_log_length_slope_per_sec',b),('source_last11_log_length_intercept',a),('source_last11_r_squared',r2),('source_initial10_length_um_mean',l0)]:
   assert abs(float(d[field])-v) < 1e-9*max(1,abs(v)),(key,field,d[field],v)
  if r2<=0.8:reason='R2_LE_0p8'
  elif d['end_type'] in ('pruned','touchtop'):reason='PRUNED_OR_TOUCHTOP'
  elif d['end_type']!='div':reason='NO_OBSERVED_DIV'
  elif b==0:reason='ZERO_SLOPE'
  else:
   reason='PASS_TRUE_SURVIVOR'
   lag=((math.log(l0)-a)/b-72*3600)/3600
   assert abs(float(d['source_Rmd_lag_hours'])-lag)<1e-8,(key,lag,d['source_Rmd_lag_hours'])
   passed[key[0]].append(lag)
 assert reason==d['rmd_selection_reason'],(key,reason,d['rmd_selection_reason'])
 reasons[reason]+=1
receipt=json.loads((P/'E500_mother_regrowth_lag_original_receipt.json').read_text())
assert reasons==receipt['regrowth_lag_selection_counts']
assert len(passed['M9_M9zero_60h'])==58 and len(passed['M9_M9zero_60h_NO_LIGHT'])==38
for name,lag in passed.items():
 assert abs(statistics.median(lag)-receipt['group_lag_stats'][name]['median_lag_h'])<1e-8
summary={'project':901,'stage':'E500','QA_PASS':True,'source_Feather_SHA':h.hexdigest(),
         'source_length_rows_reparsed_from_frozen_source_export':N,'source_true_mother_IDs':len(rows),
         'per_source_mother_all_fits_recomputed':238,'recomputed_author_11frames_R2_and_lags_exactly_match_source_first_parser':True,
         'source_fit_rejection_counts':dict(reasons),'passed_lag_by_condition_n':{k:len(v) for k,v in passed.items()},
         'passed_lag_median_by_condition_h':{k:statistics.median(v) for k,v in passed.items()},
         'source_note':'Independent stdlib numerical reconstruction from original Arrow row-preserving export; same source, not biological replicate; source Rmd original fits all cell ranks whereas this scoped analysis targets one mother cohort only'}
(P/'E500_independent_QA_regrowth_lag.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False));print(json.dumps(summary,indent=2,ensure_ascii=False))
