#!/usr/bin/env python3
"""901 E500: source-pinned illumination/no-fluorescence-excitation regrowth audit.
Original source: published Gervais2026 Gervais2026_data/Exported_data/myframes_control.feather.
The exact original Feather is already independently extracted with ZIP member CRC/SHA in GitHub Actions.
No claim of randomized exposure or exposure-specific molecular photobleaching, transfer bits.
"""
from pathlib import Path
from collections import Counter
import hashlib,json,math
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import StratifiedTable

P=Path(__file__).resolve().parent
F=P/'myframes_control.feather'
expected_sha='5cb1850cb127a64bf5e3dffa4d4e13e7cc648f73af176f6d03cf4d0e1fa9bb7c'
assert F.stat().st_size==302295912 and hashlib.sha256(F.read_bytes()).hexdigest()==expected_sha
source=json.loads((P/'E500_original_source_receipt.json').read_text())
assert source['source_member_crc_pass'] is True and source['original_member_sha256']==expected_sha
invent=json.loads((P/'E500_original_control_source_inventory.json').read_text())
assert invent['source_rows']==1011012
cases=pd.read_csv(P/'E500_original_illumination_mother_fate_candidates.csv',dtype={'date':str,'cell':str,'gl_id':str,'strain':str,'condition':str})
assert len(cases)==238 and cases.date.nunique()==1 and set(cases.date)=={'20221220'}
assert cases.cell.nunique()==238 and not cases.duplicated(['date','condition','cell']).any()
assert cases.strain.nunique()==4
assert set(cases.fate)=={'Survive','No Growth','Unknown'}
assert (cases.fate.eq('Survive')==cases.end_type.eq('div')).all()
assert (cases.starvation_nominal_h=='60h').all()
assert set(cases.condition)=={'M9_M9zero_60h','M9_M9zero_60h_NO_LIGHT'}
assert (cases.illuminated == ~cases.condition.str.contains('NO_LIGHT')).all()
assert (cases.time_start<cases.m_start).all()
assert cases.gl_id.str.split('.').str[0].eq('20221220').all()
cases['microscope_position']=cases.gl_id.str.split('.').str[1]
assert cases.groupby('microscope_position').illuminated.nunique().max()==1
assert (cases.groupby('microscope_position').strain.nunique()==1).all()
cases['regrew']=cases.fate=='Survive'
cases.to_csv(P/'E500_source_238_mother_cell_fate_with_positions.csv',index=False)

counts=pd.crosstab(cases.illuminated,cases.fate)
assert int(counts.loc[True,'Survive'])==85 and int(counts.loc[False,'Survive'])==43
assert int(counts.loc[True,'Unknown'])==6 and int(counts.loc[False,'Unknown'])==0
lit=cases[cases.illuminated];dark=cases[~cases.illuminated]
ratio_lit=lit.regrew.mean();ratio_dark=dark.regrew.mean();diff=ratio_dark-ratio_lit
p=fisher_exact([[int(lit.regrew.sum()),int((~lit.regrew).sum())],[int(dark.regrew.sum()),int((~dark.regrew).sum())]]).pvalue
# Source biological units are positions in one acquisition date, not independent experimental dates.
strains=[];tabs=[]
for name,g in cases.groupby('strain',sort=True):
 L=g[g.illuminated];D=g[~g.illuminated]
 a,b=int(L.regrew.sum()),len(L)-int(L.regrew.sum())
 c,d=int(D.regrew.sum()),len(D)-int(D.regrew.sum())
 assert min(len(L),len(D))>0
 strains.append({'strain':name,'illuminated_mothers':len(L),'illuminated_regrew':a,'illuminated_unknown':int(L.fate.eq('Unknown').sum()),
  'no_fluorescence_excitation_mothers':len(D),'no_fluorescence_excitation_regrew':c,'absolute_regrowth_gap_no_light_minus_lit':c/(c+d)-a/(a+b)})
 tabs.append(np.array([[c,d],[a,b]]))
pd.DataFrame(strains).to_csv(P/'E500_strain_stratified_illumination_regrowth.csv',index=False)
cmh=StratifiedTable(tabs)
fields=cases.groupby(['condition','illuminated','microscope_position','strain'],sort=True).agg(n_mothers=('cell','nunique'),n_survived=('regrew','sum'),n_unknown=('fate',lambda s:int(s.eq('Unknown').sum()))).reset_index()
fields['regrowth_fraction']=fields.n_survived/fields.n_mothers
fields.to_csv(P/'E500_position_stratified_illumination_regrowth.csv',index=False)

# simple position block bootstrap within illuminated and no-light conditions.
# Descriptive clustering over microscope positions; insufficient independent acquisition dates for formal cross-experiment CI.
rng=np.random.default_rng(901500)
arrays={}
for illum in [True,False]:
 g=fields[fields.illuminated==illum]
 arrays[illum]=g[['n_mothers','n_survived']].to_numpy(dtype=int)
assert len(arrays[True])==18 and len(arrays[False])==6
B=10000;boot=np.empty(B)
for i in range(B):
 rates={}
 for illum,vals in arrays.items():
  draw=vals[rng.integers(0,len(vals),size=len(vals))]
  rates[illum]=draw[:,1].sum()/draw[:,0].sum()
 boot[i]=rates[False]-rates[True]
# censoring sensitivity: if all 6 UNKNOWN in illuminated group survive, still large difference.
unknown=int(lit.fate.eq('Unknown').sum())
ranges={'unknown_as_nonregrow':diff,'unknown_as_regrow':ratio_dark-(lit.regrew.sum()+unknown)/len(lit),'unknown_excluded':ratio_dark-(lit.regrew.sum()/(len(lit)-unknown))}
report={
 'study':'Gervais et al 2026 Molecular Systems Biology DOI 10.1038/s44320-026-00226-5',
 'stage':'E500','official_original_sha256':expected_sha,'official_original_feather_rows':invent['source_rows'],
 'acquisition_dates':sorted(cases.date.unique()),'illuminated':{'n':len(lit),'survived':int(lit.regrew.sum()),'unknown':unknown,'rate':ratio_lit,'microscope_positions':len(arrays[True])},
 'no_fluorescence_excitation':{'n':len(dark),'survived':int(dark.regrew.sum()),'unknown':int(dark.fate.eq('Unknown').sum()),'rate':ratio_dark,'microscope_positions':len(arrays[False])},
 'overall_regrowth_fraction_gap_dark_minus_illuminated':diff,'risk_ratio_dark_to_illuminated':ratio_dark/ratio_lit,
 'source_4_strains':strains,'all_strain_conditional_gaps_positive':all(x['absolute_regrowth_gap_no_light_minus_lit']>0 for x in strains),
 'stratified_mantel_haenszel_odds_ratio_no_light_over_lit':float(cmh.oddsratio_pooled),
 'stratified_mantel_haenszel_test_p_iid_assumption':float(cmh.test_null_odds().pvalue),
 'source_cell_level_fisher_exact_p_iid_assumption':float(p),
 'unknown_fate_treatment_sensitivity':ranges,
 'position_bootstrap_10000_descriptive_gap_quantiles':list(map(float,np.quantile(boot,[0.025,0.5,0.975]))),
 'position_resampling_seed':901500,
 'physical_position_assignment_confounded':True,
 'treatment':'NO_LIGHT=without fluorescence excitation during starvation; source still has phase contrast time-lapse imaging',
 'source_assumption':'Source author classified GM.2 resumed divisions from cells existing before regrowth switch as regrowth; no source randomization of microscope positions demonstrated.',
 'inference':'Strong same-date observational fluorescence illumination/regrowth association across four reporter strains, NOT independent biological replication, NOT causal per-photon phototoxicity, NOT bleaching rate, not molecule count.',
 'molecular_transfer_MI':'NOT_IDENTIFIED','bit_per_second':'NOT_IDENTIFIED','bit_per_joule':'NOT_IDENTIFIED','publication':'PUBLICATION_NO_GO'
}
(P/'E500_results.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
print(json.dumps({k:report[k] for k in ('illuminated','no_fluorescence_excitation','overall_regrowth_fraction_gap_dark_minus_illuminated','unknown_fate_treatment_sensitivity','position_bootstrap_10000_descriptive_gap_quantiles','stratified_mantel_haenszel_odds_ratio_no_light_over_lit')},indent=2))
