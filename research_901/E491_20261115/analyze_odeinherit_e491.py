#!/usr/bin/env python3
"""E491 retrospective same-source-ID-prefix and aligned-time identity diagnostic; never molecule-transfer or causal MI.
Original processed six-protein ODEinherit source is denoised and 5x interpolated from 48 recorded times.
No independent optical technical repeats are available. No claim of point-in-time prediction.
"""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import hashlib, json, csv
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
RAW=ROOT/'cell_lineage_data.rda'
META=ROOT/'metadata_original.csv'
LONG=ROOT/'time_series_original_long.csv'
SOURCE_SHA='4e8c7b89e607f56d9e77c67f6479768f9883513c914f7da65330447f41306f11'
PROTEINS=['Cdc10','Stb3','CLB5','Whi5','Xbp1','Tup1']
assert hashlib.sha256(RAW.read_bytes()).hexdigest()==SOURCE_SHA
meta=pd.read_csv(META,dtype={'cell_id':'string','mother_id':'string'})
raw=pd.read_csv(LONG)
assert list(meta.columns)==['cell_id','mother_id','cell_birth_timepoint']
assert len(meta)==85 and meta.cell_id.nunique()==85 and (meta.mother_id=='root').sum()==25
assert len(raw)==85*6*240 and raw.source_value.notna().sum()==96570
assert list(raw.protein.drop_duplicates())==PROTEINS
# The source was published as exactly 240 interpolated points, not independent measurements.
raw['protein']=pd.Categorical(raw['protein'],categories=PROTEINS,ordered=True)
raw=raw.sort_values(['cell_index','protein','time_index'])
flat=raw.source_value.to_numpy(dtype=float)
arr=flat.reshape((85,6,240))
# Assert ordering block-level exact (avoid silent reshape errors).
assert raw.cell_index.to_numpy()[0]==1 and raw.cell_index.to_numpy()[-1]==85
assert raw.loc[raw.cell_index.eq(1),'protein'].iloc[:240].eq(PROTEINS[0]).all()
assert (meta.cell_id.to_numpy()==raw.groupby('cell_index').cell_id.first().to_numpy()).all()
id2idx={str(v):i for i,v in enumerate(meta.cell_id)}

# Root ancestry is the biological experimental unit for cluster descriptions.
def root_id(cell):
    visited=set(); cur=cell
    while True:
       if cur in visited:raise RuntimeError('lineage cycle')
       visited.add(cur)
       r=meta.iloc[id2idx[cur]]
       if r.mother_id=='root':return cur
       cur=str(r.mother_id)

roots={cell:root_id(str(cell)) for cell in meta.cell_id}
parents=set(str(x) for x in meta.mother_id if x!='root')
assert all(p in id2idx for p in parents)
results=[]
for rid,row in meta.iterrows():
   if row.mother_id=='root':continue
   b=int(row.cell_birth_timepoint)
   if b<=5:continue # source first available at 5; no preceding maternal source before b=5.
   p=str(row.mother_id); child=str(row.cell_id)
   assert b>=10 and b<=240
   field=child.split('-')[0]
   mom=id2idx[p]
   # All alternative cells share the identifier prefix (not verified same movie/batch), and are known mothers,
   # not cells known only as children. This is a competition audit, not randomized counterfactual.
   candidate=[id2idx[z] for z in parents if z.split('-')[0]==field
              and z!=child and int(meta.iloc[id2idx[z]].cell_birth_timepoint)<=b-1
              and np.all(np.isfinite(arr[id2idx[z],:,b-2:b]))]
   candidate=sorted(set(candidate))
   assert mom in candidate
   assert len(candidate)>=2
   target=arr[rid,:,b-1]
   assert np.all(np.isfinite(target))
   for typ,tcol in [('prebirth_adjacent_interpolated',b-2),('coincident_birth_interpolated',b-1)]:
     cand_arr=arr[candidate,:,tcol]
     assert np.isfinite(cand_arr).all()
     absdiff=np.abs(cand_arr-target[None,:])
     pos=candidate.index(mom)
     # percentile = fraction of alternative candidate mothers further away. Fraction between 0 and 1.
     for k,protein in enumerate(PROTEINS):
        dist=float(absdiff[pos,k]); alt=np.delete(absdiff[:,k],pos)
        results.append(dict(child=child,mother=p,root_family=roots[child],field=field,
           birth_timepoint=b,time_alignment=typ,protein=protein,
           matching_n_candidates=len(candidate),n_impostors=len(alt),
           matched_absdiff=dist,matched_closest=int(dist<min(alt)),
           fraction_impostors_farther=float(np.mean(alt>dist)),
           fraction_impostors_at_least_as_far=float(np.mean(alt>=dist)),
           average_impostor_absdiff=float(np.mean(alt)),
           median_impostor_absdiff=float(np.median(alt))))

out=pd.DataFrame(results)
assert len(out)==48*6*2
out.to_csv(ROOT/'E491_parent_identity_per_protein.csv',index=False,float_format='%.15g')
summary={
 'stage':'901_E491','source':'Wu et al ODEinherit, bioRxiv 2025 DOI 10.1101/2024.11.23.624995',
 'source_RData_sha256':SOURCE_SHA,'source_cell_count':85,'source_founder_count':25,
 'source_protein_names':PROTEINS,'source_daughters_total':60,
 'excluded_birth_at_first_available_grid_n':12,
 'source_subsequent_daughters_evaluable':48,
 'source_input_raw_timepoints':48,
 'source_input_interpolated_timepoints':240,
 'source_pretreatment':'functional PCA denoised and local polynomial 5x interpolated',
 'source_identifier_prefix_groups_not_certified_batches':3,
 'observation_kind':'postprocessed reporter levels in six proteins; fluorescence units not absolute molecule counts',
 'analysis_type':'retrospective same-prefix same-time parent identity negative-control; prefix is not independently certified experimental batch, not point-in-time or causal',
 'definitions':{'metric':'per protein fraction of same-source-prefix contemporaneous alternative mothers with larger absolute expression difference',
                'chance_reference':'0.5 under artificial exchangeability of candidate parent labels; not a biological null',
                'comparison':'biological mother observed at b-1 or b vs daughter at b, for b>5; 240 points are interpolated',
                'uncertainty':'descriptive family-bootstrap if reported, not independent experimental replication'},
 'by_alignment':{},'cross_linked_source':False,
 'real_molecular_transfer_bits':'NOT_IDENTIFIED',
 'physical_photon_noise_identified':False,
 'independent_exposure_repeat_verified':False,
 'publication_status':'PUBLICATION_NO_GO',
}
for key,sub in out.groupby('time_alignment'):
 groups={}
 for prot,s in sub.groupby('protein',sort=False):
    groups[prot]={'n_daughters':len(s),'n_root_families':int(s.root_family.nunique()),
      'n_candidate_mother_median':float(s.matching_n_candidates.median()),
      'mean_fraction_impostors_farther':float(s.fraction_impostors_farther.mean()),
      'matched_strict_closest_fraction':float(s.matched_closest.mean()),
      'artificial_label_expected_strict_closest':float((1/s.matching_n_candidates).mean()),
      'mean_absdiff_matched':float(s.matched_absdiff.mean()),
      'mean_absdiff_impostors':float(s.average_impostor_absdiff.mean())}
 summary['by_alignment'][key]=groups
# Founder block resampling only: frozen seed, no biological CI claim.
np_rng=np.random.default_rng(901491)
extra={}
for key,sub in out.groupby('time_alignment'):
  # 25 source founders; some founders may have only one child born at b=5 and thus be absent.
  founder_codes=sorted(sub.root_family.unique());by_protein={}
  for protein,s in sub.groupby('protein'):
    g={n:val.fraction_impostors_farther.to_numpy() for n,val in s.groupby('root_family')}
    means=[]
    for _ in range(2000):
       ids=np_rng.choice(founder_codes,size=len(founder_codes),replace=True)
       concat=[g[k] for k in ids if k in g]
       if concat:means.append(float(np.concatenate(concat).mean()))
    by_protein[protein]=list(np.quantile(means,[.025,.5,.975]))
  extra[key]=by_protein
summary['within_data_root_bootstrap_2000_quantiles_2p5_50_97p5']=extra
(ROOT/'E491_parent_identity_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'audit_events':48,'founder_count':25,'source_dates':3,
 'postdivision':summary['by_alignment']['coincident_birth_interpolated'],
 'prebirth':summary['by_alignment']['prebirth_adjacent_interpolated']},ensure_ascii=False,indent=2))
