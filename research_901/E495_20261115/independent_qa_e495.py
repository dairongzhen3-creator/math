#!/usr/bin/env python3
"""901 E495 independent stdlib audit of new raw optical identity gate.
No old model retraining, pandas, numpy, or primary E495 function imports.
Not a second experiment and does not identify molecules.
"""
from pathlib import Path
from collections import Counter,defaultdict
import csv,json,hashlib,statistics,math
ROOT=Path(__file__).resolve().parent
SOURCE={
'E494_frozen_true_divisions_size_fluorescence_split.csv':'57cf0104f01fd93ba1ddde9da085c01b25d96851499828a84eca373a04a47c7a',
'E494_leave_date_out_mse.csv':'a74c4419a1efbae1a9b90d425f405051ace1e89833c869b5f97a50fd8bbf4a46'}
for n,h in SOURCE.items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h
counts=Counter();d=defaultdict(list);schema=defaultdict(list);keys=set();lanes=set()
with (ROOT/'E494_frozen_true_divisions_size_fluorescence_split.csv').open(newline='',encoding='utf8') as f:
 for r in csv.DictReader(f):
  k=tuple(r[p] for p in ('date','lane_file','parent_id','parent_end_frame'))
  assert k not in keys;keys.add(k)
  date=r['date'];s=r['schema'];x=float(r['fraction_top_length']);y=float(r['fraction_top_fluorescence'])
  assert 0<x<1 and 0<y<1 and abs((y-x)-float(r['fluorescence_minus_length_fraction']))<1e-13
  counts[s]+=1;d[date].append((x,y));schema[s].append((x,y));lanes.add(r['lane_file'])
assert len(keys)==42658 and len(d)==27 and len(lanes)==857
assert counts==Counter({'dual':26773,'single':15885})

def slope(rows):
 n=len(rows);sx=sy=sxx=sxy=0.0
 for x,y in rows:sx+=x;sy+=y;sxx+=x*x;sxy+=x*y
 return (sxy-sx*sy/n)/(sxx-sx*sx/n)

lock={}
with (ROOT/'E494_leave_date_out_mse.csv').open(newline='',encoding='utf8') as f:
 for r in csv.DictReader(f):lock[r['date']]=r
assert len(lock)==27 and set(lock)==set(d)
by=[];total_ideal=total_base=total_fit=0.0
for date in sorted(d):
 records=d[date]; orig=lock[date]
 assert len(records)==int(orig['events'])
 ideal=sum((y-x)**2 for x,y in records)
 base=float(orig['mse_constant_total']);fit=float(orig['mse_size_total'])
 b=slope(records)
 total_ideal+=ideal;total_base+=base;total_fit+=fit
 by.append((date,len(records),ideal,base,fit,b))
assert sum(n for _,n,*_ in by)==42658
primary=json.loads((ROOT/'E495_results.json').read_text())['new_identity_gate']
assert abs(total_ideal-primary['total_sse_direct_identity'])<1e-9
assert abs(total_base-primary['total_sse_E494_frozen_date_heldout_constant'])<1e-10
assert abs(total_fit-primary['total_sse_E494_frozen_date_heldout_regression'])<1e-10
assert sum(fit<ideal for _,n,ideal,base,fit,b in by)==27
assert sum(ideal<base for _,n,ideal,base,fit,b in by)==12
assert all(0<b<1 for *_,b in by)
for s,records in schema.items():assert abs(slope(records)-primary['slope_by_schema'][s]['observed_slope_OLS'])<1e-11
with (ROOT/'E495_by_date_direct_identity_gate.csv').open(newline='',encoding='utf8') as f:
 out={r['date']:r for r in csv.DictReader(f)}
assert set(out)==set(d)
for date,n,ideal,base,fit,b in by:
 assert abs(float(out[date]['same_cell_identity_qF_equal_qL_SSE'])-ideal)<1e-11
 assert abs(float(out[date]['source_observed_optical_slope'])-b)<1e-11
probe=json.loads((ROOT/'E495_partial_official_tar_index.json').read_text())
assert probe['http_status']==206 and probe['compressed_bytes_observed']==50331648
assert probe['content_range']=='bytes 0-50331647/4430365592'
assert len(probe['members'])==6 and all('/preproc/' not in k['name'] for k in probe['members'])
qa={
 'stage':'E495','pass':True,'source_SHA_checked':True,'original_verified_event_count':len(keys),
 'true_dates':len(d),'true_lanes':len(lanes),'n_heldout_days_fit_better_than_identity':27,
 'n_heldout_days_identity_better_than_baseline':12,
 'independent_stdlib_OLS_slope_by_schema':{s:slope(v) for s,v in schema.items()},
 'independent_global_sse':{'ideal_direct_equality':total_ideal,'frozen_baseline':total_base,'frozen_E494_size_regression':total_fit},
 'independent_ideal_mse':total_ideal/len(keys),'independent_ideal_vs_baseline_increase':total_ideal/total_base-1,
 'recomputed_per_date_slopes_27_all_between_0_and_1':True,
 'verified_source_range_response_206':True,
 'independent_original_357_calibration_fit':'NOT_DONE_NO_357_RAW_PAIR_INPUT',
 'source_scope':'Optical identity proxy validation and compressed range probe, NOT physical molecule transfer or biological replication'}
(ROOT/'E495_independent_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2))
print(json.dumps(qa,indent=2))
