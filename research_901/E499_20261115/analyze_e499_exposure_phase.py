#!/usr/bin/env python3
"""901 E499 source-pinned Gervais 2026 exposure-phase identifiability analysis.

Reads immutable E498 source-export and E498 per-cell WLS values. Adds new imaging
cadence, per-scheduled-frame diagnostics and same-cell matched phase identities.
Does not estimate physical molecular counts or independent camera bleaching yield.
"""
from pathlib import Path
import pandas as pd, numpy as np
from collections import Counter
from hashlib import sha256
from math import exp
import json

P=Path(__file__).resolve().parent
SRC=P.parent/'p901_E498' if (P.parent/'p901_E498').exists() else P/'source_E498'
SOURCE_SHA={
 'myframes_bleaching.feather.csv.gz':'0a9b630b95af875ae6bba5be5bac5430acbf96cfcf4bc767576965e59f8251db',
 'myframes_bleaching.feather':'6e1733d6b641f5a47dca2851520fb9a0edfed2aad6faecbaf1ec5eece99917f7',
 'E498_bleaching_per_cell_weighted_slopes.csv':'3e2a7c5a620d767ad9f0f2192f72008d52af90e0d0fd30d3eac188c0f0d0b45a',
 'Gervais_et_al_analysis.Rmd':'7f915a34db6c3326a0be833bc1c04e57127d90684ce281c9fe4b0cd65a48de70',
 'E498_results.json':'644838bc37b9a6a3a519f468470a12d622b3bbe07c435c65ba66dc19e3ac7632'}
for f,h in SOURCE_SHA.items():
 actual=sha256((SRC/f).read_bytes()).hexdigest()
 assert actual==h,(f,actual,h)

b=pd.read_csv(SRC/'myframes_bleaching.feather.csv.gz',dtype={'date':str,'cell':str},low_memory=False)
assert len(b)==91838 and set(b.idx)=={1,2}
assert not b.duplicated(['date','strain','Fluo_sig','exp_series','idx','cell','time_sec']).any()
b.sort_values(['date','strain','Fluo_sig','exp_series','idx','cell','time_sec'],inplace=True)
b['delta_min']=b.groupby(['date','strain','Fluo_sig','exp_series','idx','cell'],sort=False).time_sec.diff()/60.
assert b['delta_min'].dropna().gt(0).all()
cadence=[]
for (date,strain,fluo,series,phase),v in b.groupby(['date','strain','Fluo_sig','exp_series','idx'],sort=True):
 gaps=v.delta_min.dropna()
 mode=float(gaps.mode().iat[0])
 planned=3 if phase==1 else 12
 assert mode==planned,(date,strain,fluo,phase,mode)
 cadence.append({'date':date,'strain':strain,'Fluo_sig':fluo,'exp_series':int(series),'phase_idx':int(phase),
                 'phase':'growth' if phase==1 else 'starvation',
                 'native_planned_minutes':planned,'native_actual_mode_minutes':mode,
                 'source_frame_rows':int(len(v)),'source_unique_cells':int(v.cell.nunique()),
                 'time_gaps_n':int(len(gaps)),'gaps_exact_planned_n':int((gaps==planned).sum()),
                 'fraction_gaps_at_planned_minutes':float((gaps==planned).mean()),
                 'median_actual_interval_min':float(gaps.median()),
                 'gaps_outside_planned_n':int((gaps!=planned).sum())})
cad=pd.DataFrame(cadence).sort_values(['date','Fluo_sig','phase_idx'])
assert len(cad)==14 and cad.date.nunique()==3 and cad.groupby('phase_idx').size().to_dict()=={1:7,2:7}
cad.to_csv(P/'E499_source_exposure_cadence.csv',index=False,float_format='%.16g')

c=pd.read_csv(SRC/'E498_bleaching_per_cell_weighted_slopes.csv',dtype={'date':str,'cell':str},low_memory=False)
assert len(c)==1659 and set(c.phase_idx)=={1,2}
assert not c.duplicated(['date','experiment_series','phase_idx','cell','Fluo_sig']).any()
assert c.source_author_selection_positive.astype(str).str.lower().isin(['true','false']).all()
c['planned_interval_minutes']=np.where(c.phase_idx==1,3,12)
c['nominal_frames_per_hour']=60/c.planned_interval_minutes
c['apparent_log_loss_per_scheduled_frame']=c.raw_degradation_corrected_signal_rate_per_h/c.nominal_frames_per_hour
c['apparent_frac_loss_per_scheduled_frame']=1-np.exp(-c.apparent_log_loss_per_scheduled_frame)
c.to_csv(P/'E499_all_cell_phase_exposure_normalized_slopes.csv',index=False,float_format='%.16g')

p=[]
for (date,fluo),g in c.groupby(['date','Fluo_sig'],sort=True):
 out={'date':date,'Fluo_sig':fluo,'experiment_series':int(g.experiment_series.iloc[0])}
 for code,label in [(1,'growth'),(2,'starvation')]:
  z=g[g.phase_idx==code]
  assert len(z)>=6
  allrates=z.raw_degradation_corrected_signal_rate_per_h.to_numpy()
  zr=allrates[allrates>0]
  out.update({f'{label}_n':len(z),f'{label}_n_positive':len(zr),
              f'{label}_median_all_per_h':float(np.median(allrates)),
              f'{label}_median_positive_per_h':float(np.median(zr)),
              f'{label}_planned_min':3 if code==1 else 12,
              f'{label}_median_positive_per_planned_frame':float(np.median(zr)*(3 if code==1 else 12)/60)})
 out['starve_over_growth_positive_median_hourly_ratio']=out['starvation_median_positive_per_h']/out['growth_median_positive_per_h']
 out['starve_over_growth_positive_median_per_frame_ratio']=4*out['starve_over_growth_positive_median_hourly_ratio']
 out['starvation_minus_growth_quarter_positive_median_per_h']=out['starvation_median_positive_per_h']-out['growth_median_positive_per_h']/4
 out['shared_nonnegative_dark_loss_invariance_diagnostic_PASS']=out['starvation_minus_growth_quarter_positive_median_per_h']>=0
 p.append(out)
phase=pd.DataFrame(p)
assert len(phase)==7 and (phase.starve_over_growth_positive_median_hourly_ratio<.25).all()
phase.to_csv(P/'E499_source_phase_probe_frequency_normalization.csv',index=False,float_format='%.16g')

paired=[]
for (date,fluo,series),g in c.groupby(['date','Fluo_sig','experiment_series'],sort=True):
 a=g[g.phase_idx==1][['cell','raw_degradation_corrected_signal_rate_per_h']].rename(columns={'raw_degradation_corrected_signal_rate_per_h':'growth_k_per_h'})
 d=g[g.phase_idx==2][['cell','raw_degradation_corrected_signal_rate_per_h']].rename(columns={'raw_degradation_corrected_signal_rate_per_h':'starvation_k_per_h'})
 merged=a.merge(d,on='cell',how='inner',validate='one_to_one')
 for row in merged.itertuples(index=False):
  kg,ks=float(row.growth_k_per_h),float(row.starvation_k_per_h)
  paired.append({'date':date,'Fluo_sig':fluo,'exp_series':int(series),'cell':row.cell,
                 'growth_k_per_h':kg,'starvation_k_per_h':ks,
                 'growth_log_loss_per_planned_frame':kg/20,'starvation_log_loss_per_planned_frame':ks/5,
                 'growth_k_greater_than_starve_k':kg>ks,
                 'growth_per_frame_greater_than_starve_per_frame':kg/20>ks/5,
                 'both_phase_positive_rates':kg>0 and ks>0})
matched=pd.DataFrame(paired)
assert len(matched)==16 and set(matched.Fluo_sig)=={'GFP','YFP','CFP'}
matched.to_csv(P/'E499_same_cell_both_phase_overlap.csv',index=False,float_format='%.16g')

# Verify design rank and nonidentifiability without optimizing any new parameter.
X=np.column_stack((np.ones(len(c)), (c.phase_idx==2).astype(int),c.nominal_frames_per_hour))
assert np.linalg.matrix_rank(X)==2
# Source geometry: only 2 phase levels, and f=20-15*I(starvation), exact rank deficiency.
# For every probe-date combination, k_phase = beta*frequency_phase + d_phase.
# Two observed k's for 3 unknowns, even with nonnegative beta,d_phase infinite admissible solutions for >0 k.
# For the stronger extra assumption d_growth=d_starve=d>=0, inferred d=(4*k_s-k_g)/3.
res={'project':901,'stage':'E499','date':'2026-11-15',
 'source_SHA256':SOURCE_SHA,'original_bleaching_rows':len(b),
 'source_cadence':{'n_probe_date_groups':len(phase),'n_phase_groups':len(cad),
 'growth_nominal_min':3,'starvation_nominal_min':12,
 'fraction_native_gap_growth':float(cad.loc[cad.phase_idx==1,'gaps_exact_planned_n'].sum()/cad.loc[cad.phase_idx==1,'time_gaps_n'].sum()),
 'fraction_native_gap_starvation':float(cad.loc[cad.phase_idx==2,'gaps_exact_planned_n'].sum()/cad.loc[cad.phase_idx==2,'time_gaps_n'].sum()),
 'duplicate_cell_time_keys':0,
 'source_actual_recorded_growth_to_starvation_frequency_ratio':4},
 'source_author_existing_bleach_WLS':{'n_cell_phase':len(c),'n_positive':int((c.raw_degradation_corrected_signal_rate_per_h>0).sum()),
 'exact_phase_probe_date_groups':len(phase)},
 'new_positive_only_comparison':{'n_groups':len(phase),'group_hourly_starve_over_growth_ratio_min_median_max':[float(phase.starve_over_growth_positive_median_hourly_ratio.min()),float(phase.starve_over_growth_positive_median_hourly_ratio.median()),float(phase.starve_over_growth_positive_median_hourly_ratio.max())],
 'group_per_frame_starve_over_growth_ratio_min_median_max':[float(phase.starve_over_growth_positive_median_per_frame_ratio.min()),float(phase.starve_over_growth_positive_median_per_frame_ratio.median()),float(phase.starve_over_growth_positive_median_per_frame_ratio.max())],
 'groups_below_constant_per_frame_0p25_hourly':int((phase.starve_over_growth_positive_median_hourly_ratio<.25).sum()),
 'groups_with_nonnegative_shared_dark_loss_under_same_cell_model':int(phase.shared_nonnegative_dark_loss_invariance_diagnostic_PASS.sum())},
 'same_exact_cell_cross_phase':{'n':len(matched),'by_probe':{k:int(v) for k,v in matched.Fluo_sig.value_counts().to_dict().items()},
 'growth_hourly_greater':int(matched.growth_k_greater_than_starve_k.sum()),
 'growth_per_frame_greater':int(matched.growth_per_frame_greater_than_starve_per_frame.sum()),
 'both_positive':int(matched.both_phase_positive_rates.sum()),
 'same_cell_cross_phase_was_observed_in_3_of_7_probe_date_groups':True},
 'structural_identifiability':{'design':'k_{date,probe,phase}=beta_{date,probe}*f_phase + d_{date,probe,phase}, with f_growth=20/h, f_starvation=5/h',
 'linear_dependence':'f_phase = 20 - 15*I(starvation); phase and frequency perfectly collinear within dataset',
 'actual_design_matrix_rank':2,'requested_design_columns':3,
 'physical_beta_vs_phase_loss':'NOT_SEPARATELY_IDENTIFIED',
 'source_published_camera_exposure_time_and_excitation':'article Methods reports nominal settings but no per-cell photon dose/technical repeated same state at both frequencies',
 'strict_extra_assumption_d_growth_equals_d_starvation_nonnegative':'not supported by source; negative residual for positive-only group medians is a sensitivity diagnostic, not physical negative rate'},
 'interpretation':'Source signal is author degradation-corrected; k per h is optical signal decrease, per-planned-frame representation only; neither true photobleaching yield per exposure nor protein degradation corrected in independent wet-lab manner',
 'molecular_count_transfer_bits':'NOT_IDENTIFIED','bit_per_joule':'NOT_IDENTIFIED','publication':'PUBLICATION_NO_GO'}
(P/'E499_results.json').write_text(json.dumps(res,ensure_ascii=False,indent=2))
print(json.dumps({'source_cadence':res['source_cadence'],'new_positive_only':res['new_positive_only_comparison'],'matched':res['same_exact_cell_cross_phase']},ensure_ascii=False,indent=2))
