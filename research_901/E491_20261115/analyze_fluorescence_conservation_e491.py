#!/usr/bin/env python3
"""E491 independent fluorescence-sum audit of E488 frozen divisions.
We DO NOT estimate protein molecules or fit Rosenfeld et al 2006 binomial calibration.
E487 Julou sustained expression / temporal offset violate no-synthesis pulse-chase premise.
"""
from pathlib import Path
import os
import hashlib, json
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parent
SOURCE_A=Path(os.environ.get('E488_FEATURES_PATH',str(ROOT/'E488_same_event_size_and_background.csv')))
SOURCE_B=Path(os.environ.get('E488_CV_PATH',str(ROOT/'E488_cv_predictions.csv')))
assert SOURCE_A.exists() and SOURCE_B.exists()
sa=hashlib.sha256(SOURCE_A.read_bytes()).hexdigest()
sb=hashlib.sha256(SOURCE_B.read_bytes()).hexdigest()
a=pd.read_csv(SOURCE_A)
b=pd.read_csv(SOURCE_B,usecols=['date','lane_file','parent_id','parent_end_frame'])
key=['date','lane_file','parent_id','parent_end_frame']
assert not a.duplicated(key).any() and not b.duplicated(key).any()
c=a.merge(b,on=key,how='inner',validate='one_to_one')
assert len(c)==42658 and c.date.nunique()==27 and c.lane_file.nunique()==857
assert c['daughter_birth_frame'].eq(c['parent_end_frame']+1).all()
cols=['date','lane_file','schema','parent_id','parent_end_frame','daughter_birth_frame','parent_last_amp','daughter_sum_amp','parent_last_length_px','daughter_sum_length_px']
d=c[cols].copy()
assert (d[['parent_last_amp','daughter_sum_amp','parent_last_length_px','daughter_sum_length_px']]>0).all().all()
d['fluorescence_sum_ratio']=d.daughter_sum_amp/d.parent_last_amp
d['length_sum_ratio']=d.daughter_sum_length_px/d.parent_last_length_px
d.to_csv(ROOT/'E491_frozen_division_sum_audit.csv',index=False,float_format='%.15g')

def stats(x):
 r=x.fluorescence_sum_ratio
 q=x.length_sum_ratio
 return {'n_events':len(x),'n_dates':int(x.date.nunique()),'n_lanes':int(x.lane_file.nunique()),
         'fluorescence_ratio_median':float(r.median()),'fluorescence_ratio_iqr':[float(r.quantile(.25)),float(r.quantile(.75))],
         'fluorescence_ratio_tail_2p5_97p5':[float(r.quantile(.025)),float(r.quantile(.975))],
         'fraction_within_5pct_of_unity':float(r.between(.95,1.05).mean()),
         'fraction_daughters_sum_greater_than_parent':float((r>1).mean()),
         'length_ratio_median':float(q.median()),
         'positive_source_fraction':float((r>0).mean())}
res={'stage':'E491','date':'2026-11-15','status':'OBSERVED_SOURCE_FLUORESCENCE_RATIO_ONLY',
 'source_E488_event_features_sha256':sa,'source_E488_cv_keys_sha256':sb,
 'selection':'Source E488 held-out 42,658 events, 27 author-QC-included dates, 857 lanes',
 'schema':{k:stats(x) for k,x in [('all',d),*list(d.groupby('schema'))]},
 'published_2006_method':'Rosenfeld et al Biophys J DOI 10.1529/biophysj.105.073098',
 'physical_mass_balance_interpretation':'Raw fluorescence measured at different acquisition frames; no matching independent noise calibration; ongoing reporter production, maturation, segmentation and bleaching may influence ratio. Unchecked N_parent=N_daughter1+N_daughter2 assumption is not valid to impose on this source.',
 'molecule_per_fluorescence_parameter':'NOT_IDENTIFIED','physical_transfer_information_bit':'NOT_IDENTIFIED'}

dates=d.groupby('date',sort=True).fluorescence_sum_ratio
res['date_median_min_max']=[float(dates.median().min()),float(dates.median().max())]
date_summary=d.groupby('date').agg(n=('fluorescence_sum_ratio','size'),median_ratio=('fluorescence_sum_ratio','median'),within5pct=('fluorescence_sum_ratio',lambda x:float(x.between(.95,1.05).mean())),aboveunity=('fluorescence_sum_ratio',lambda x:float((x>1).mean()))).reset_index()
date_summary.to_csv(ROOT/'E491_division_sum_by_date.csv',index=False,float_format='%.15g')
# Block-resampling by source dates is descriptive; 27 batches not a population CI.
rng=np.random.default_rng(90149102)
bydate={day:v.fluorescence_sum_ratio.to_numpy() for day,v in d.groupby('date')}
keys=list(bydate)
medians=[]
for _ in range(2000):
 draw=rng.choice(keys,size=len(keys),replace=True)
 medians.append(float(np.median(np.concatenate([bydate[k] for k in draw]))))
res['source_date_block_bootstrap_median_ratio_2000_2p5_50_97p5']=list(map(float,np.quantile(medians,[.025,.5,.975])))
(ROOT/'E491_fluorescence_ratio_summary.json').write_text(json.dumps(res,indent=2,ensure_ascii=False))
print(json.dumps(res,indent=2,ensure_ascii=False))
