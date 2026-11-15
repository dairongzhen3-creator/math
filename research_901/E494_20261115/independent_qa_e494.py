#!/usr/bin/env python3
"""Independent E494 stdlib row-level QA from *two frozen E488 original CSVs*.
Different decoder and sufficient-statistics formulas; does not import primary pandas/numpy code.
This is an independent computational audit, not optical technical replicate.
"""
from pathlib import Path
from collections import defaultdict,Counter
from hashlib import sha256
from math import sqrt,log2
import csv,json
ROOT=Path(__file__).resolve().parent
P=ROOT/'E488_same_event_size_and_background.csv';Q=ROOT/'E488_cv_predictions.csv'
assert sha256(P.read_bytes()).hexdigest()=='7e4dadde0b754c78e21fa22ab03f9fdcca911a608cecefe7ab7b23ce0a527072'
assert sha256(Q.read_bytes()).hexdigest()=='b06ada09bbeb9b536305d165b08550d62c67a9a39fc7360d0409be2d4438c8a1'
KEYS=('date','lane_file','parent_id','parent_end_frame')
def key(r):return tuple(r[k] for k in KEYS)
with Q.open(newline='',encoding='utf8') as f:
    xs={key(r) for r in csv.DictReader(f)}
assert len(xs)==42658
records=[];all_original_keys=set();source_by_key={}
with P.open(newline='',encoding='utf8') as f:
    for row in csv.DictReader(f):
        k=key(row)
        assert k not in all_original_keys
        all_original_keys.add(k)
        if k not in xs:continue
        top=float(row['daughter_top_amp']);bottom=float(row['daughter_bottom_amp'])
        a=float(row['daughter_top_length_px']);b=float(row['daughter_bottom_length_px'])
        assert min(top,bottom,a,b)>0
        fluor=top/(top+bottom);length=a/(a+b)
        assert abs(float(row['daughter_sum_amp'])-top-bottom)<1e-9
        assert abs(float(row['daughter_sum_length_px'])-a-b)<1e-10
        cat='TOP_LARGER' if a>b else 'BOTTOM_LARGER' if a<b else 'LENGTH_EQUAL'
        records.append((k[0],k[1],row['schema'],fluor,length,cat,int(top>bottom)))
        source_by_key[k]=(row['schema'],fluor,length,cat,int(top>bottom))
assert len(records)==42658
assert len({x[0] for x in records})==27 and len({x[1] for x in records})==857

def r_of(rows):
 n=len(rows); sx=sy=sxx=syy=sxy=0.0
 for x in rows:
  a,b=x[4],x[3];sx+=a;sy+=b;sxx+=a*a;syy+=b*b;sxy+=a*b
 return (sxy-sx*sy/n)/sqrt((sxx-sx*sx/n)*(syy-sy*sy/n))

def lin(rows):
 n=len(rows);sx=sy=sxx=sxy=0
 for z in rows:
  x,y=z[4],z[3];sx+=x;sy+=y;sxx+=x*x;sxy+=x*y
 beta=(sxy-sx*sy/n)/(sxx-sx*sx/n)
 return sy/n-beta*sx/n,beta
r=r_of(records)
no_tie=[x for x in records if x[5]!='LENGTH_EQUAL']
agree=sum(((x[4]>.5)==(x[3]>.5)) for x in no_tie)
bydate=defaultdict(list)
for rec in records:bydate[rec[0]].append(rec)
assert len(bydate)==27
all_names=sorted(bydate)
base_sum=size_sum=0.0;good_days=0
score_sum=0.0;pos_score_dates=0
for d in all_names:
 train=[x for day,items in bydate.items() if day!=d for x in items]
 test=bydate[d]
 trained={};n0=n1=0.0;day_gain=0.0
 for schema in ('single','dual'):
  s=[x for x in train if x[2]==schema]
  trained[schema]=(sum(x[3] for x in s)/len(s),lin(s))
  counts=Counter((x[5],x[6]) for x in s)
  positive=sum(x[6] for x in s);p0=(positive+1)/(len(s)+2)
  for v in (x for x in test if x[2]==schema):
   mid,(a,b)=trained[schema]
   n0+=(v[3]-mid)**2
   n1+=(v[3]-(a+b*v[4]))**2
   c=v[5];n=counts[(c,0)]+counts[(c,1)]
   p1=(counts[(c,1)]+1)/(n+2)
   gain=log2(p1/p0) if v[6] else log2((1-p1)/(1-p0))
   score_sum+=gain;day_gain+=gain
 base_sum+=n0;size_sum+=n1
 good_days+=int(n1<n0)
 pos_score_dates+=int(day_gain>0)
source=json.loads((ROOT/'E494_results.json').read_text())
assert abs(r-source['split_pearson_r'])<1e-12
assert len(no_tie)==source['n_non_tied_daughter_lengths']==35718
assert agree==source['split_size_sign_fluorescence_sign_agreement_n']
reduction=1-size_sum/base_sum
assert abs(reduction-source['out_of_date_linear_prediction']['pooled_relative_mse_reduction'])<1e-11
assert good_days==27
bitgain=score_sum/len(records)
assert abs(bitgain-source['out_of_date_optical_logscore']['bits_gain_per_event'])<1e-12
assert pos_score_dates==26
# Row-level original key + both split fractions cross check, without reading the full source copy from primary.
with (ROOT/'E494_frozen_true_divisions_size_fluorescence_split.csv').open(newline='',encoding='utf8') as f:
 derived=list(csv.DictReader(f))
assert len(derived)==len(records)
derived_by_key={key(r):r for r in derived}
assert len(derived_by_key)==len(records)
assert set(derived_by_key)==set(source_by_key)
for k, v in source_by_key.items():
    rec=derived_by_key[k]
    assert v[0]==rec['schema']
    assert abs(v[1]-float(rec['fraction_top_fluorescence']))<1e-13
    assert abs(v[2]-float(rec['fraction_top_length']))<1e-13
    assert v[3]==rec['size_category']
    assert v[4]==int(rec['top_brighter'])
qa={'stage':'E494','passed':True,'n_source_events':len(records),'n_original_feature_keys':len(all_original_keys),
 'n_heldout_event_keys':len(xs),'n_dates':len(bydate),'source_lanes':len({x[1] for x in records}),
 'n_non_length_ties':len(no_tie),'n_sign_direction_match':agree,
 'independent_direct_pearson':r,'independent_heldout_mse_gain_fraction':reduction,
 'heldout_days_mse_improves':good_days,'independent_heldout_optical_log_score_gain_per_event':bitgain,
 'heldout_days_optical_log_score_positive':pos_score_dates,
 'method':'stdlib CSV streaming and independent sufficient-statistic calculations; no primary pandas/numpy imported; same laboratory observations; no real technical calibration',
 'physical_protein_transfer_bits':'NOT_IDENTIFIED'}
(ROOT/'E494_independent_QA.json').write_text(json.dumps(qa,indent=2,ensure_ascii=False),encoding='utf8')
print(json.dumps(qa,indent=2,ensure_ascii=False))
