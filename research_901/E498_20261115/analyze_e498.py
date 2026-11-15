#!/usr/bin/env python3
"""901 E498 primary: source-pinned Gervais 2026 Feather-export re-analysis.
Original Feather bytes are in authenticated E498 GitHub archive; gzipped CSV
conversions preserve rows and float precision. Every assertion is source-specific.
NOT independent wet-lab work. NOT protein copy numbers or molecular MI.
"""
from pathlib import Path
from collections import Counter
import csv, gzip, hashlib, json, math, statistics
import numpy as np
import pandas as pd

P=Path(__file__).resolve().parent
expect={
 'mydegrates.feather.csv.gz':'bcf2a0e4854b17b7c1bbe5dd95f527320874e80c46ccf5c8cd2d64bd4dd7c447',
 'myframes_degradation.feather.csv.gz':'65665284ff303ef6c88f0115d0037fb36465fde42131b200632bd9f3044c5455',
 'myframes_bleaching.feather.csv.gz':'0a9b630b95af875ae6bba5be5bac5430acbf96cfcf4bc767576965e59f8251db',
 'myproteomicsdata_multistrain.feather.csv.gz':'dcd7e982e3774a06c15372e967ff79a2083f344e8172c35615fbd64dcf9ff4ed'}
for name,h in expect.items():assert hashlib.sha256((P/name).read_bytes()).hexdigest()==h,(name,h)
author=pd.read_csv(P/'mydegrates.feather.csv.gz')
assert list(author.columns)==['Fluo1','gamma_dec_Fluo','gamma_init'] and len(author)==4
assert set(author.Fluo1)=={'GFP','YFP','RFP','CFP'}
# Model-only integral of published author fitted rate gamma0*exp(-omega*t).
surv=[]
for _,x in author.iterrows():
    a,w=float(x.gamma_init),float(x.gamma_dec_Fluo)
    assert a>0 and w>0
    rec={'Fluo1':x.Fluo1,'author_gamma_init_per_h':a,'author_omega_per_h':w,
         'decay_rate_halving_time_h':math.log(2)/w,'model_only_infinite_time_extrapolation':math.exp(-a/w)}
    for t in (5,10,30,60): rec[f'model_old_fluorescent_label_survival_{t}h']=math.exp(-(a/w)*(-math.expm1(-w*t)))
    surv.append(rec)
surv=pd.DataFrame(surv).sort_values('Fluo1')
surv.to_csv(P/'E498_author_model_derived_survival.csv',index=False,float_format='%.15g')

# Source degradation: every biological time point has 4 spectral-channel rows.
d=pd.read_csv(P/'myframes_degradation.feather.csv.gz',dtype={'cell':str,'date':str},low_memory=False)
assert len(d)==1116240 and d.date.nunique()==2 and set(d.Fluo)=={'GFP','YFP','RFP','CFP'}
map_f={'pTrc99a_mVenus':'YFP','pTrc99a_mCerulean':'CFP','pTrc99a_mGFP':'GFP','pTrc99a_mCherry2':'RFP'}
assert set(d.strain)==set(map_f)
framekey=['date','cell','time_sec']
frames=d.groupby(framekey).Fluo.nunique()
assert len(frames)==279060 and frames.eq(4).all() and not d.duplicated(framekey+['Fluo']).any()
observed=d[d.Fluo.eq(d.strain.map(map_f))].copy()
assert len(observed)==279060
starved=observed[observed.c_step.eq('SPM.1')].copy()
starved_positive=starved[starved.Ampl.gt(0)&starved.time_sec.gt(starved.m_start+3600)].copy()
starved_positive.sort_values(['date','strain','cell','time_sec'],inplace=True)
G=starved_positive.groupby(['date','strain','cell'],sort=False)
starved_positive['previous_time_sec']=G.time_sec.shift(1)
starved_positive['previous_positive_Ampl']=G.Ampl.shift(1)
v=starved_positive[starved_positive.previous_positive_Ampl.gt(0)&starved_positive.previous_time_sec.notna()].copy()
v['dt_hours']=(v.time_sec-v.previous_time_sec)/3600.
v=v[v.dt_hours.between(2,7)].copy()
v['apparent_log_fluorescence_loss_per_hour']=-np.log(v.Ampl/v.previous_positive_Ampl)/v.dt_hours
v['t_mid_after_starvation_h']=(v.time_sec+v.previous_time_sec-2*v.m_start)/7200.
assert np.isfinite(v[['apparent_log_fluorescence_loss_per_hour','t_mid_after_starvation_h']]).all().all()
assert len(v)==3517 and v.cell.nunique()==439
v['window']=pd.cut(v.t_mid_after_starvation_h,[0,3,12,18,30,36,60,100],labels=['0_3','3_12','12_18','18_30','30_36','36_60','after60'],include_lowest=True,right=False)
columns=['date','strain','Fluo','cell','previous_time_sec','time_sec','m_start','dt_hours','t_mid_after_starvation_h','apparent_log_fluorescence_loss_per_hour','window']
v[columns].to_csv(P/'E498_degradation_verified_time_intervals.csv',index=False,float_format='%.15g')
date_reporter=v.groupby(['date','strain','Fluo'],observed=True).agg(n_intervals=('cell','size'),n_unique_cells=('cell','nunique'),
            median_apparent_loss_per_h=('apparent_log_fluorescence_loss_per_hour','median'),
            mean_apparent_loss_per_h=('apparent_log_fluorescence_loss_per_hour','mean'),
            fraction_apparent_loss_positive=('apparent_log_fluorescence_loss_per_hour',lambda x:float((x>0).mean()))).reset_index()
date_reporter.to_csv(P/'E498_degradation_source_experiments.csv',index=False,float_format='%.15g')
period=v.groupby(['date','strain','Fluo','window'],observed=True).agg(n_intervals=('cell','size'),n_unique_cells=('cell','nunique'),
       median_loss_per_hour=('apparent_log_fluorescence_loss_per_hour','median'),
       mean_loss_per_hour=('apparent_log_fluorescence_loss_per_hour','mean')).reset_index()
period.to_csv(P/'E498_degradation_by_time_window.csv',index=False,float_format='%.15g')

# A paired within-cell comparison (early 3-12h vs late 18-30h), one median per cell per window.
earlylate=v[v.window.isin(['3_12','18_30'])].groupby(['date','strain','Fluo','cell','window'],observed=True).apparent_log_fluorescence_loss_per_hour.median().unstack('window')
paired=earlylate.dropna(subset=['3_12','18_30']).reset_index()
paired['early_minus_late_per_hour']=paired['3_12']-paired['18_30']
paired['early_greater_than_late']=(paired['early_minus_late_per_hour']>0).astype(int)
assert len(paired)==403 and paired.early_greater_than_late.sum()==368
paired.to_csv(P/'E498_degradation_403_within_cell_paired.csv',index=False,float_format='%.15g')
paired_by=paired.groupby(['date','strain','Fluo']).agg(n_cells=('cell','size'),n_early_greater=('early_greater_than_late','sum'),
    median_early_loss_per_h=('3_12','median'),median_late_loss_per_h=('18_30','median'),
    median_within_cell_early_minus_late_per_h=('early_minus_late_per_hour','median')).reset_index()
paired_by.to_csv(P/'E498_degradation_paired_by_source.csv',index=False,float_format='%.15g')

# Low-frequency degradation has 2 acquisition dates, not 1.1 million independent experiments.
# A cell with multiple time intervals is still one tracked physical cell.

# Source bleaching is already degradation-corrected fluorescence signal according to author README.
b=pd.read_csv(P/'myframes_bleaching.feather.csv.gz',dtype={'date':str,'cell':str},low_memory=False)
assert len(b)==91838 and b.date.nunique()==3 and b.signal.gt(0).all()
b['idx']=b.idx.astype(int)
b['t_hours']=b.time_sec/3600.
bleach=[]
for (date,exp_series,idx,cell,Fluo_sig),g in b.groupby(['date','exp_series','idx','cell','Fluo_sig'],sort=False):
    if len(g)<=5:continue
    t=g.t_hours.to_numpy(dtype=float);y=np.log(g.signal.to_numpy(dtype=float));weights=1./g.signal.to_numpy(dtype=float)
    weights=weights/weights.sum();mu_t=float(np.dot(weights,t));mu_y=float(np.dot(weights,y))
    var_t=float(np.dot(weights,(t-mu_t)**2));cov_ft=float(np.dot(weights,(t-mu_t)*(y-mu_y)))
    if var_t<=0:continue
    k=-cov_ft/var_t
    bleach.append({'date':date,'experiment_series':exp_series,'phase_idx':int(idx),'cell':cell,'Fluo_sig':Fluo_sig,
                   'n_source_frames':len(g),'raw_degradation_corrected_signal_rate_per_h':k,
                   'source_author_selection_positive':bool(k>0),'weighted_var_time':var_t})
bleach=pd.DataFrame(bleach)
assert len(bleach)==1659
bleach.sort_values(['date','experiment_series','Fluo_sig','phase_idx','cell']).to_csv(P/'E498_bleaching_per_cell_weighted_slopes.csv',index=False,float_format='%.15g')
bleach_group=bleach.groupby(['date','Fluo_sig','phase_idx']).agg(n_cells=('cell','size'),n_positive=('source_author_selection_positive','sum'),
        median_all_source_rate_per_h=('raw_degradation_corrected_signal_rate_per_h','median'),
        median_positive_only_per_h=('raw_degradation_corrected_signal_rate_per_h',lambda s:float(s[s>0].median())),
        n_source_frames=('n_source_frames','sum')).reset_index()
bleach_group['fraction_positive']=bleach_group.n_positive/bleach_group.n_cells
bleach_group['phase']=np.where(bleach_group.phase_idx==1,'growth','starvation')
bleach_group.to_csv(P/'E498_bleaching_by_date_fluor_phase.csv',index=False,float_format='%.15g')
# YFP in all three actual source dates provides within-study treatment/phase sensitivity.
yfp=bleach_group[bleach_group.Fluo_sig.eq('YFP')].pivot(index='date',columns='phase',values='median_positive_only_per_h')
assert len(yfp)==3 and (yfp.growth>yfp.starvation).all()

# Proteomics data: 10 biological reporter strain backgrounds; no independent technical calibration.
q=pd.read_csv(P/'myproteomicsdata_multistrain.feather.csv.gz',low_memory=False)
assert len(q)==130200 and q.Genes.nunique()==2604 and q.GFP_prom.nunique()==10
assert set(q['sample'])=={1,2,3,4,5} and set(q.time_h)=={-2.5,-1.5,1.4,21.4,45.4}
assert not q.duplicated(['GFP_prom','Genes','sample']).any()
n_na=int(q.Log_Quantity_norm.isna().sum())
base=q[q['sample'].isin([1,2])].groupby(['GFP_prom','Genes']).Log_Quantity_norm.mean().rename('source_baseline_log')
q=q.join(base,on=['GFP_prom','Genes'])
z=q[q['sample'].isin([3,4,5])].copy()
z['delta_ln_normalized_signal']=(z.Log_Quantity_norm-z.source_baseline_log)*math.log(2.) # author Rmd has this conversion
z=z[z.delta_ln_normalized_signal.notna()].copy()
assert z['GFP_prom'].nunique()==10
# Per-gene cross-reporter-strain dispersion is descriptive not same-strain technical repeat variance.
pergene=z.groupby(['Genes','time_h']).delta_ln_normalized_signal.agg(n_backgrounds='count',sd_source_across_backgrounds='std',
    observed_span_ln='ptp',mean_source_delta_ln='mean') if False else None
pergene=z.groupby(['Genes','time_h']).delta_ln_normalized_signal.agg(
    n_backgrounds='count',sd_source_across_backgrounds='std',
    observed_span_ln=lambda s:float(s.max()-s.min()),mean_source_delta_ln='mean').reset_index()
pergene.to_csv(P/'E498_proteomics_per_gene_time_cross_strain_dispersion.csv',index=False,float_format='%.15g')
prot_group=[]
for t,g in z.groupby('time_h'):
    gpg=pergene[(pergene.time_h==t)&(pergene.n_backgrounds>=5)].copy()
    pivot=g.pivot(index='Genes',columns='GFP_prom',values='delta_ln_normalized_signal')
    spearman=pivot.corr(method='spearman',min_periods=500).to_numpy()
    upper=spearman[np.triu_indices(len(spearman),k=1)]
    upper=upper[np.isfinite(upper)]
    prot_group.append({'time_h':float(t),'source_genes_with_at_least_5_strain_backgrounds':len(gpg),
        'source_genes_with_all_10_strain_backgrounds':int((pergene[(pergene.time_h==t)].n_backgrounds==10).sum()),
        'median_across_backgrounds_gene_std_ln':float(gpg.sd_source_across_backgrounds.median()),
        'median_across_backgrounds_gene_span_ln':float(gpg.observed_span_ln.median()),
        'median_pairwise_spearman_between_strain_backgrounds':float(np.median(upper)),
        'pairwise_source_strain_comparisons':len(upper)})
pdg=pd.DataFrame(prot_group)
pdg.to_csv(P/'E498_proteomics_cross_strain_summary.csv',index=False,float_format='%.15g')

result={'project':901,'stage':'E498','date':'2026-11-15',
 'study':'Gervais et al 2026 Molecular Systems Biology DOI 10.1038/s44320-026-00226-5',
 'official_zenodo_record_id':18059846,'official_large_zip_fully_downloaded':False,
 'source_feather_original_verified_sha256':{z['name']:z['feather_sha256'] for z in json.loads((P/'E498_original_to_csv_gzip_receipt.json').read_text())['files']},
 'source_gzip_csv_sha256':expect,
 'degradation':{'all_four_channel_rows':len(d),'unique_cell_time_keys_four_dyes_each':len(frames),
                'source_study_acquisition_dates':sorted(d.date.unique().tolist()),
                'reporter_matched_spectral_rows':len(observed),
                'matched_starvation_source_rows':len(starved),
                'matched_starvation_nonmissing_positive_after_first_hour':len(starved_positive),
                'source_independent_low_exposure_valid_intervals_2_to_7h':len(v),
                'tracked_physical_cells_with_valid_intervals':int(v.cell.nunique()),
                'date_strain_reporter_groups':len(date_reporter),
                'source_paired_cells_3to12h_vs_18to30h':len(paired),
                'source_paired_cells_earlier_rate_greater':int(paired.early_greater_than_late.sum()),
                'within_cell_rate_positive_delta_fraction':float(paired.early_greater_than_late.mean()),
                'group_summary':date_reporter.to_dict('records'),
                'per_group_paired':paired_by.to_dict('records')},
 'author_source_degradation_model':{'form':'gamma(t)=gamma0 exp(-omega*t); S(t)=exp[-gamma0/omega*(1-exp(-omega*t))]',
                                    'rates_unit':'per hour as in source Rmd time_sec/3600',
                                    'all_four_reporters':surv.to_dict('records'),
                                    'extrapolation_warning':'Model-defined survival of pre-existing tagged fluorescent reporter, not directly observed molecule count or full lifetime survival; no absolute copy calibration.'},
 'bleaching':{'source_rows':len(b),'source_dates':sorted(b.date.unique().tolist()),
             'cell_phase_weighted_slopes_n6plus':len(bleach),
             'positive_author_selection_cell_phase_slopes':int(bleach.source_author_selection_positive.sum()),
             'positive_fraction':float(bleach.source_author_selection_positive.mean()),
             'yfp_growth_greater_than_starvation_median_source_phase_dates':len(yfp),
             'model_scope':'Signal already source degradation-corrected; apparent photobleaching slopes include imaging schedule, not per-exposure photon quantum yield.',
             'by_group':bleach_group.to_dict('records')},
 'proteomics':{'source_rows':len(q),'source_genes':int(q.Genes.nunique()),'source_strain_backgrounds':int(q.GFP_prom.nunique()),
               'source_times_h':sorted(float(x) for x in q.time_h.unique()),
               'source_log_signal_missing_rows':n_na,
               'missing_fraction':n_na/len(q),
               'reported_signal_not_absolute_copy':True,
               'per_source_time':pdg.to_dict('records'),
               'source_limit':'Per-gene concentration ratios judged insufficiently precise by original paper; strain differences are not same-strain independent technical repeats.'},
 'raw_mother_daughter_transfer_information_bits':'NOT_IDENTIFIED',
 'biological_independent_experiment_count_not_equal_to_rows':True,
 'no_merge_with_Kscheschinski_2026_lineages':True,
 'publication_status':'PUBLICATION_NO_GO'}
(P/'E498_results.json').write_text(json.dumps(result,indent=2,ensure_ascii=False,default=str),encoding='utf8')
print(json.dumps({k:result[k] for k in ('degradation','bleaching','proteomics')},ensure_ascii=False,indent=2)[:22000])
