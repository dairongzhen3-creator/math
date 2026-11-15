#!/usr/bin/env python3
"""901 E495: new E494-frozen observational optical proxy test + bounded source gate.

Does NOT re-fit 2018 lambda, infer molecule counts, or rerun E487-E494 models.
Every statistic is computed only from 42,658 E494 frozen derived events, preexisting
E494 date-heldout error tables, the original 2018 JATS and Zenodo source probe.
"""
import hashlib
import json
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

P = Path(__file__).resolve().parent
SOURCE_SHA = {
 'E494_frozen_true_divisions_size_fluorescence_split.csv':'57cf0104f01fd93ba1ddde9da085c01b25d96851499828a84eca373a04a47c7a',
 'E494_leave_date_out_mse.csv':'a74c4419a1efbae1a9b90d425f405051ace1e89833c869b5f97a50fd8bbf4a46',
 'E494_results.json':'e1ee9497b3be0492281fede89d6f353d67f9ae5605c673ad8505d7fea0d8dd9d',
 'PMC5768764_2018_fullTextXML.xml':'a3e5574d61551fad47044d61200455bac079847785f0927a24a3da4d3e8869d8',
 'zenodo_746230_ReadMe.md':'4aa76d90746c50dbbdf11812962a38c55586428103bd97619f91e9c1162dbf30',
}
for name, expect in SOURCE_SHA.items():
    actual = hashlib.sha256((P/name).read_bytes()).hexdigest()
    assert actual == expect, (name,actual,expect)

src = pd.read_csv(P/'E494_frozen_true_divisions_size_fluorescence_split.csv',dtype={'date':str,'lane_file':str},low_memory=False)
pre = pd.read_csv(P/'E494_leave_date_out_mse.csv',dtype={'date':str})
old = json.loads((P/'E494_results.json').read_text())
assert len(src)==42658 and src['date'].nunique()==27 and src.lane_file.nunique()==857
assert len(pre)==27 and pre.date.nunique()==27
assert int(pre.events.sum())==42658
assert not src.duplicated(['date','lane_file','parent_id','parent_end_frame']).any()
assert set(src.date)==set(pre.date)
assert set(src.schema)=={'single','dual'}
x=src.fraction_top_length.to_numpy(dtype=float)
y=src.fraction_top_fluorescence.to_numpy(dtype=float)
assert np.isfinite(x).all() and np.isfinite(y).all() and np.all((x>0)&(x<1)) and np.all((y>0)&(y<1))
assert np.max(np.abs(y-x-src.fluorescence_minus_length_fraction.to_numpy(dtype=float)))<1e-13


def slope(x,y):
    cx=x-x.mean();cy=y-y.mean();return float(np.dot(cx,cy)/np.dot(cx,cx))

bydate=[]
for date, g in src.groupby('date',sort=True):
    xx=g.fraction_top_length.to_numpy(dtype=float)
    yy=g.fraction_top_fluorescence.to_numpy(dtype=float)
    beta=slope(xx,yy)
    alpha=float(yy.mean()-beta*xx.mean())
    ideal_sse=float(np.sum((yy-xx)**2))
    bydate.append({'date':str(date),'n_events':len(g),'n_lanes':int(g.lane_file.nunique()),
          'source_observed_optical_slope':beta,'source_observed_optical_intercept':alpha,
          'same_cell_identity_qF_equal_qL_SSE':ideal_sse,
          'same_cell_identity_qF_equal_qL_MSE':ideal_sse/len(g),
          'source_split_pearson_r':float(np.corrcoef(xx,yy)[0,1])})

check=pd.DataFrame(bydate).merge(pre,left_on='date',right_on='date',validate='one_to_one')
assert (check.n_events==check.events).all()
assert (check.n_lanes==check.lanes).all()
assert np.isclose(check.mse_constant_total.sum()/len(src),old['out_of_date_linear_prediction']['pooled_mse_baseline'],atol=1e-14)
assert np.isclose(check.mse_size_total.sum()/len(src),old['out_of_date_linear_prediction']['pooled_mse_size'],atol=1e-14)
assert len(check)==27
check['E494_frozen_fitted_better_than_direct_identity']=check.mse_size_total < check.same_cell_identity_qF_equal_qL_SSE
check['direct_identity_better_than_E494_frozen_constant']=check.same_cell_identity_qF_equal_qL_SSE < check.mse_constant_total
check.to_csv(P/'E495_by_date_direct_identity_gate.csv',index=False,float_format='%.16g')

schema={}
for kind,g in src.groupby('schema',sort=True):
   a=g.fraction_top_length.to_numpy(dtype=float);b=g.fraction_top_fluorescence.to_numpy(dtype=float)
   schema[kind]={'events':len(g),'observed_slope_OLS':slope(a,b),
       'observed_intercept_OLS':float(b.mean()-slope(a,b)*a.mean()),
       'ideal_identity_MSE':float(np.mean((b-a)**2)),
       'source_optical_fluorescence_std':float(b.std()),
       'no_physical_volume_calibration':True}

identity_sse=float(check.same_cell_identity_qF_equal_qL_SSE.sum())
base_sse=float(check.mse_constant_total.sum())
fitted_sse=float(check.mse_size_total.sum())
assert fitted_sse < base_sse < identity_sse

probe=json.loads((P/'E495_partial_official_tar_index.json').read_text())
assert probe['http_status']==206 and probe['content_range']=='bytes 0-50331647/4430365592'
assert probe['compressed_bytes_observed']==48*1024*1024
assert all('/preproc/' not in r['name'] for r in probe['members'])
assert len(probe['members'])==6
assert probe['outcome']=='STREAM_CAPPED_OR_TRUNCATED'
assert probe['members'][3]['name'].startswith('DIMM_MoMA_data/MoMA/lactose/')

root=ET.parse(P/'PMC5768764_2018_fullTextXML.xml').getroot()
article_text=' '.join(''.join(root.itertext()).split())
assert '357 sibling pairs' in article_text and '0.0361' in article_text
assert 'GFP levels at birth were estimated in each daughter cell as average' in article_text
assert 'Using the same criteria as in the decay analysis' in article_text

result={
 'project':901,'stage':'E495','date':'2026-11-15',
 'study':'Kaiser et al. (2018) Nat Commun 9:212 DOI 10.1038/s41467-017-02505-0',
 'frozen_E494_input_sha256':SOURCE_SHA,
 'experimental_observation_scope':{'observed_divisions':len(src),'dates':27,'lanes':857,'schema_counts':{k:v['events'] for k,v in schema.items()},'E494_cross_date_MSE_not_rerun':True},
 'new_identity_gate':{
   'formula':'strict raw processed optical fraction proxy model qF_pred=qL',
   'total_sse_direct_identity':identity_sse,'total_sse_E494_frozen_date_heldout_constant':base_sse,
   'total_sse_E494_frozen_date_heldout_regression':fitted_sse,
   'mse_identity':identity_sse/len(src),'mse_old_E494_baseline':base_sse/len(src),
   'mse_old_E494_linear':fitted_sse/len(src),
   'relative_identity_error_increase_over_old_baseline':identity_sse/base_sse-1,
   'relative_identity_error_increase_over_old_fitted':identity_sse/fitted_sse-1,
   'n_dates_fitted_better_than_identity':int(check.E494_frozen_fitted_better_than_direct_identity.sum()),
   'n_dates_identity_better_than_baseline':int(check.direct_identity_better_than_E494_frozen_constant.sum()),
   'slope_by_schema':schema,
   'daily_observed_slope_min_med_max':list(map(float,check.source_observed_optical_slope.agg(['min','median','max']).to_list())),
   'n_days_daily_slope_positive':int((check.source_observed_optical_slope>0).sum()),
   'n_days_daily_slope_below_unity':int((check.source_observed_optical_slope<1).sum()),
   'new_scientific_reading':'The exact-length-proxy equality qF=qL has worse out-of-date squared error than old source-only regression; not a test of true molecular binomial model, as qL is estimated length not calibrated cytoplasm volume.'
 },
 'bounded_original_archive_access':{'zenodo_record':824793,
   'dataset_file':'DIMM_MoMA_data.tar.gz','full_official_bytes':4430365592,
   'full_official_md5':'299f477f71067e764419381556ac8183',
   'http_206_range_supported':True,
   'compressed_bytes_probed':probe['compressed_bytes_observed'],
   'archive_entries_seen':len(probe['members']),
   'first_three_regular_source_entries':[x['name'] for x in probe['members'] if x['type']=='0'],
   'preproc_or_357_pair_events_accessible_in_probed_prefix':False,
   'full_archive_downloaded':False,'full_archive_checksum_verified':False,
   'caveat':'gzip-compressed tar is not demonstrated seekable at arbitrary member offsets; HTTP range alone does not yield a verified inner preproc file index. Do not claim data absent.'},
 'original_paper_357_pair_test':{'N':357,'author_original_fitted_lambda':0.0361,
   'author_95pct_posterior_interval':[0.026,0.112],
   'source_357_pair_individual_values':'NOT_RECEIVED',
   'author_method':'glucose-phase nonproduction; per-cell GFP extrapolated to birth after separately estimated loss; measured cytoplasm fraction and measurement noise',
   'lambda_independent_replication':'NOT_DONE',
   'cross_source_lambda_transfer':'NOT_QUALIFIED'},
 'physical_molecule_transfer_bit':'NOT_IDENTIFIED',
 'protein_molecule_absolute_count':'NOT_IDENTIFIED',
 'bit_per_second':'NOT_IDENTIFIED','bit_per_joule':'NOT_IDENTIFIED',
 'research_status':'PUBLICATION_NO_GO'
}
(P/'E495_results.json').write_text(json.dumps(result,indent=2,ensure_ascii=False))
print(json.dumps(result['new_identity_gate'],indent=2))
print(json.dumps(result['bounded_original_archive_access'],indent=2))
