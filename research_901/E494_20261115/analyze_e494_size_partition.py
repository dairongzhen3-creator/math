#!/usr/bin/env python3
"""E494 new analysis on E488 *frozen* genuine mother/sister events.
Outcome is uncalibrated cell fluorescence amplitude split. Length split is not volume split.
Scientific NOT_IDENTIFIED: molecular transmission / copy / bits / energy.
"""
from pathlib import Path
from collections import defaultdict
from hashlib import sha256
from math import log2
import json, csv
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
P = ROOT/'E488_same_event_size_and_background.csv'
Q = ROOT/'E488_cv_predictions.csv'
EXP_P='7e4dadde0b754c78e21fa22ab03f9fdcca911a608cecefe7ab7b23ce0a527072'
EXP_Q='b06ada09bbeb9b536305d165b08550d62c67a9a39fc7360d0409be2d4438c8a1'
assert sha256(P.read_bytes()).hexdigest()==EXP_P
assert sha256(Q.read_bytes()).hexdigest()==EXP_Q
keys=['date','lane_file','parent_id','parent_end_frame']
feat=pd.read_csv(P,low_memory=False)
pred=pd.read_csv(Q,usecols=keys,low_memory=False)
assert not feat.duplicated(keys).any() and not pred.duplicated(keys).any()
d=feat.merge(pred,on=keys,how='inner',validate='one_to_one',indicator=True).drop(columns='_merge')
assert len(d)==42658 and d.date.nunique()==27 and d.lane_file.nunique()==857
amp=['daughter_top_amp','daughter_bottom_amp','daughter_sum_amp']
length=['daughter_top_length_px','daughter_bottom_length_px','daughter_sum_length_px']
assert (d[amp+length]>0).all().all()
assert np.max(abs(d.daughter_top_amp+d.daughter_bottom_amp-d.daughter_sum_amp))<1e-9
assert np.max(abs(d.daughter_top_length_px+d.daughter_bottom_length_px-d.daughter_sum_length_px))<1e-10
assert set(d.schema)=={'single','dual'}
d['fraction_top_fluorescence']=d.daughter_top_amp/d.daughter_sum_amp
d['fraction_top_length']=d.daughter_top_length_px/d.daughter_sum_length_px
d['fluorescence_minus_length_fraction']=d.fraction_top_fluorescence-d.fraction_top_length
d['size_category']=np.where(d.fraction_top_length>.5,'TOP_LARGER',np.where(d.fraction_top_length<.5,'BOTTOM_LARGER','LENGTH_EQUAL'))
d['top_brighter']=(d.fraction_top_fluorescence>.5).astype(int)
d['lane_id']=d['lane_file'].astype(str)
assert d.date.nunique()==27

source_n={'n_events':len(d),'n_dates':int(d.date.nunique()),'n_lanes':int(d.lane_id.nunique()),'schema_counts':d.schema.value_counts().to_dict()}
size_sign=d[d.size_category!='LENGTH_EQUAL']
def corr(a,b):
 a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
 aa=a-a.mean();bb=b-b.mean()
 return float(np.dot(aa,bb)/(np.linalg.norm(aa)*np.linalg.norm(bb)))

def slope_ols(x,y):
 x=np.asarray(x);y=np.asarray(y)
 beta=float(np.dot(x-x.mean(),y-y.mean())/np.dot(x-x.mean(),x-x.mean()))
 alpha=float(y.mean()-beta*x.mean())
 return alpha,beta

# Out-of-date prediction tests with no leakage: whole source experiment day held out.
cv_rows=[]
for date in sorted(d.date.unique()):
 tr=d[d.date!=date];te=d[d.date==date]
 actual=te.fraction_top_fluorescence.to_numpy()
 baseline=np.empty(len(te));linear=np.empty(len(te))
 for schema in ['dual','single']:
  s=tr[tr.schema==schema]
  z=te[te.schema==schema]
  if not len(z):continue
  baseline[te.schema.to_numpy()==schema]=float(s.fraction_top_fluorescence.mean())
  alpha,beta=slope_ols(s.fraction_top_length.to_numpy(),s.fraction_top_fluorescence.to_numpy())
  linear[te.schema.to_numpy()==schema]=alpha+beta*z.fraction_top_length.to_numpy()
 mse0=float(np.sum((actual-baseline)**2));mse1=float(np.sum((actual-linear)**2))
 cv_rows.append({'date':str(date),'events':len(te),'lanes':te.lane_id.nunique(),
                 'mse_constant_total':mse0,'mse_size_total':mse1,'size_model_improves':mse1<mse0,
                 'heldout_mse_reduction_fraction':1-mse1/mse0})
cv=pd.DataFrame(cv_rows)
cv.to_csv(ROOT/'E494_leave_date_out_mse.csv',index=False,float_format='%.15g')

# Hold out each experimental date: probabilities P(top daughter brighter | length class, schema)
# Only a categorical optical log-score, no causal info or physical transfer bits.
score_by_date=[]
for date in sorted(d.date.unique()):
 tr=d[d.date!=date]
 te=d[d.date==date]
 for schema in ('dual','single'):
  t=tr[tr.schema==schema];v=te[te.schema==schema]
  if not len(v):continue
  p0=(int(t.top_brighter.sum())+1)/(len(t)+2)
  counts=t.groupby('size_category').top_brighter.agg(['sum','size'])
  for cat,vv in v.groupby('size_category'):
   cs=int(counts.loc[cat,'sum']) if cat in counts.index else 0
   cn=int(counts.loc[cat,'size']) if cat in counts.index else 0
   p1=(cs+1)/(cn+2)
   o=vv.top_brighter.to_numpy()
   gains=np.where(o==1, np.log2(p1/p0), np.log2((1-p1)/(1-p0)))
   score_by_date.append({'date':str(date),'schema':schema,'size_category':cat,'events':len(vv),
            'sum_gain_bit':float(gains.sum()),'mean_gain_bit_per_event':float(gains.mean())})
score=pd.DataFrame(score_by_date)
score.to_csv(ROOT/'E494_leave_date_out_logscore.csv',index=False,float_format='%.15g')

# Within-lane (camera/microchannel) effect: pooled within lane mean-centered correlation.
a=d.fraction_top_length.to_numpy();b=d.fraction_top_fluorescence.to_numpy()
lane_to_indexes=[]
for _, idx in d.groupby('lane_id',sort=True).indices.items():lane_to_indexes.append(np.asarray(idx,dtype=int))
ac=np.empty_like(a);bc=np.empty_like(b)
for ix in lane_to_indexes:
 ac[ix]=a[ix]-a[ix].mean();bc[ix]=b[ix]-b[ix].mean()
within_corr=corr(ac,bc)
# Permute pF within same lane (including date and imaging schema via lane ID), preserve all distributions.
rng=np.random.default_rng(901494)
perm_correlations=[]
# Scalar dot product faster since permutation preserves norm in each lane.
den=float(np.linalg.norm(ac)*np.linalg.norm(bc))
for _ in range(199):
 shuffled=np.empty_like(bc)
 for ix in lane_to_indexes:
  shuffled[ix]=bc[ix][rng.permutation(len(ix))]
 perm_correlations.append(float(np.dot(ac,shuffled)/den))
# Fixed random seed + direct analytic possible; uncertainty is *conditional lane reassignment*, not wetlab p.
perm_tail=(1+sum(v>=within_corr for v in perm_correlations))/(len(perm_correlations)+1)

# Day-cluster resampling, only descriptive variation within this one laboratory study.
per_date_arrays={day:g[['fraction_top_length','fraction_top_fluorescence']].to_numpy() for day,g in d.groupby('date')}
date_list=np.array(sorted(per_date_arrays));rng2=np.random.default_rng(90149402)
block_corr=[]
for _ in range(1000):
 ids=rng2.choice(date_list,size=len(date_list),replace=True)
 z=np.concatenate([per_date_arrays[k] for k in ids],axis=0)
 block_corr.append(corr(z[:,0],z[:,1]))

by_schema={}
for schema,s in d.groupby('schema'):
 no_tie=s[s.size_category!='LENGTH_EQUAL']
 by_schema[schema]={'n':len(s),'Pearson_r_size_fluo_split':corr(s.fraction_top_length,s.fraction_top_fluorescence),
                    'n_non_tied_lengths':len(no_tie),'n_size_fluo_sign_agree':int(((no_tie.fraction_top_length>.5)==(no_tie.top_brighter==1)).sum()),
                    'fraction_size_fluo_sign_agree':float(((no_tie.fraction_top_length>.5)==(no_tie.top_brighter==1)).mean()),
                    'mean_signed_optical_minus_length':float(s.fluorescence_minus_length_fraction.mean()),
                    'median_abs_optical_minus_length':float(s.fluorescence_minus_length_fraction.abs().median())}

by_date=[]
for date, s in d.groupby('date'):
 r=corr(s.fraction_top_length,s.fraction_top_fluorescence)
 by_date.append({'date':str(date),'n':len(s),'lanes':s.lane_id.nunique(),
   'pearson_r':r,'mean_optical_minus_length':float(s.fluorescence_minus_length_fraction.mean())})
pd.DataFrame(by_date).to_csv(ROOT/'E494_by_date_source_summary.csv',index=False,float_format='%.15g')

result={
 'stage':'E494','source_inputs_sha256':{P.name:EXP_P,Q.name:EXP_Q},
 'source':source_n,'measured_variables':'fluorescence amplitude, estimated pixel length at first daughter frame; no absolute molecule count, no volume',
 'n_non_tied_daughter_lengths':len(size_sign),
 'n_exact_tie_daughter_lengths':int((d.size_category=='LENGTH_EQUAL').sum()),
 'split_pearson_r':corr(a,b),
 'split_pearson_r_squared':corr(a,b)**2,
 'split_size_sign_fluorescence_sign_agreement_non_tie':float(((size_sign.fraction_top_length>.5)==(size_sign.top_brighter==1)).mean()),
 'split_size_sign_fluorescence_sign_agreement_n':int(((size_sign.fraction_top_length>.5)==(size_sign.top_brighter==1)).sum()),
 'mean_abs_fluorescence_minus_length_split':float(np.abs(d.fluorescence_minus_length_fraction).mean()),
 'median_abs_fluorescence_minus_length_split':float(np.median(abs(d.fluorescence_minus_length_fraction))),
 'mean_signed_fluorescence_minus_length_split':float(d.fluorescence_minus_length_fraction.mean()),
 'by_schema':by_schema,
 'daughter_length_category_counts_and_top_bright_fraction':{str(cat):{'events':len(g),'n_top_brighter':int(g.top_brighter.sum()),'fraction_top_brighter':float(g.top_brighter.mean())} for cat,g in d.groupby('size_category')},
 'exclude_author_special_gain_dates_20190605_20190614':{
  'n_events':int(len(d[~d.date.isin((20190605,20190614))])),
  'n_dates':int(d[~d.date.isin((20190605,20190614))].date.nunique()),
  'n_lanes':int(d[~d.date.isin((20190605,20190614))].lane_id.nunique()),
  'pearson_r':corr(d.loc[~d.date.isin((20190605,20190614)),'fraction_top_length'],d.loc[~d.date.isin((20190605,20190614)),'fraction_top_fluorescence'])},
 'out_of_date_linear_prediction':{'source_days':len(cv),'positive_days':int(cv.size_model_improves.sum()),
         'pooled_mse_baseline':float(cv.mse_constant_total.sum()/len(d)),
         'pooled_mse_size':float(cv.mse_size_total.sum()/len(d)),
         'pooled_relative_mse_reduction':float(1-cv.mse_size_total.sum()/cv.mse_constant_total.sum()),
         'min_max_day_reduction':[float(cv.heldout_mse_reduction_fraction.min()),float(cv.heldout_mse_reduction_fraction.max())]},
 'out_of_date_optical_logscore':{'bits_gain_per_event':float(score.sum_gain_bit.sum()/len(d)),'n_events':len(d),
          'positive_days':int((score.groupby('date').sum_gain_bit.sum()>0).sum()),
          'n_days':len(cv),'binomial_Laplace1':'train remaining dates and heldout source date; category based on top vs bottom pixel length incl ties, schema stratification'},
 'within_lane_conditioned':{'r':within_corr,'n_permutations':len(perm_correlations),'permutation_min':min(perm_correlations),
          'permutation_median':float(np.median(perm_correlations)),'permutation_max':max(perm_correlations),'finite_randomization_reference_upper_tail':perm_tail,
          'null':'shuffle fluorescence daughter split *within same lane* while preserving within-lane source distributions'},
 'date_block_resampling_1000_descriptive_r_quantiles':np.quantile(block_corr,[.025,.5,.975]).tolist(),
 'by_date_correlation_all_positive':all(x['pearson_r']>0 for x in by_date),
 'experimental_confound':'Common cell size, optical background, probe conditions and strain/lineage clustering. The processed fluorescence may include size-proportional autofluorescence; no causal partition or molecule count claim.',
 'special_2019_exposure_dates':'sensitivity source exists in E488; cross-source 2018 lambda unverified for this instrument/experimental context',
 'physical_molecule_transfer_bit':'NOT_IDENTIFIED','bit_per_second':'NOT_IDENTIFIED','bit_per_joule':'NOT_IDENTIFIED','publication':'PUBLICATION_NO_GO'
}
# Do not store per-cell label or genealogy in public prose, reproducible derived table is in the zip.
cols=keys+['schema','fraction_top_fluorescence','fraction_top_length','fluorescence_minus_length_fraction','size_category','top_brighter']
d[cols].sort_values(keys).to_csv(ROOT/'E494_frozen_true_divisions_size_fluorescence_split.csv',index=False,float_format='%.15g')
(ROOT/'E494_results.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf8')
print(json.dumps({k:result[k] for k in ('source','split_pearson_r','n_non_tied_daughter_lengths','n_exact_tie_daughter_lengths','split_size_sign_fluorescence_sign_agreement_non_tie','out_of_date_linear_prediction','out_of_date_optical_logscore','within_lane_conditioned','date_block_resampling_1000_descriptive_r_quantiles')},indent=2,ensure_ascii=False))
