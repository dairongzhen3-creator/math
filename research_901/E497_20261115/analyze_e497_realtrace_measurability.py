#!/usr/bin/env python3
"""901 E497: original RealTrace fit parameter identity and observed-vs-posterior sibling.
Reads two SOURCE-SHA-pinned official ZIPs and E496 immutable per-event cells.
Only scalar numeric NPZ arrays via numpy.load(allow_pickle=False), never object arrays.
Does NOT estimate independent photon measurement noise or protein molecular transfer.
"""
from pathlib import Path
from collections import defaultdict
import csv,io,json,hashlib,math,statistics,zipfile
import numpy as np
import pandas as pd

P=Path(__file__).resolve().parent
SRC=P.parent/'p901_E496'
EV=SRC/'E496_2026_original_verified_sibling_events.csv'
RAW={'mother_machine_data_sets.zip':'a5787b437d9d136fd308feab5068fd558790f81929265444b42451287bc7b421',
     'replicates_hi1_rrnB.zip':'e6378d36486757c807f03a644e9653d131c24b00ba279472b2e512a0d43a1b6f'}
PARAM_NAMES=['mean_lambda','gamma_lambda','var_lambda','mean_q','gamma_q','var_q','beta','var_x','var_g','var_dx','var_dg']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert len(RAW)==2
for name,want in RAW.items(): assert sha(SRC/name)==want
D=pd.read_csv(EV,dtype={'date':str},low_memory=False)
assert len(D)==10058 and set(D.cohort)=={'main','replicate'}
assert set(D.date)=={'20190529','20220921','20221111','20230223','20230228','20230302','20230413','20230420','20230428','20230504'} or D.date.nunique()==10
EFILE=hashlib.sha256(EV.read_bytes()).hexdigest()

# Explicitly certify all 31 original ZIP file parameter tables vs corresponding NPZ numeric scalars.
param_rows=[]
for archive,expected in RAW.items():
    with zipfile.ZipFile(SRC/archive) as z:
        assert z.testzip() is None
        preds=sorted(n for n in z.namelist() if n.endswith('_prediction.csv'))
        assert len(preds)==(23 if archive.startswith('mother_machine') else 8)
        for name in preds:
            txt=io.TextIOWrapper(z.open(name),encoding='utf8',newline='')
            lines=csv.reader(txt)
            fields=next(lines)
            assert fields==['no','name','type','init','step','lower_bound','upper_bound','final']
            params=[dict(zip(fields,next(lines))) for _ in range(11)]
            assert next(lines)==[]
            assert next(lines)[0]=='cell_id'
            assert [x['name'] for x in params]==PARAM_NAMES,(name,[q['name'] for q in params])
            partner=name.replace('_autofl_corrected_','_rawdata_').replace('_prediction.csv','_correlations.npz')
            assert partner in z.namelist(),(name,partner)
            zip_bytes=z.read(partner)
            with np.load(io.BytesIO(zip_bytes),allow_pickle=False) as ar:
                scalar_names=set(ar.files)-{'correlations'}
                assert scalar_names==set(PARAM_NAMES),(name,scalar_names)
                for row in params:
                    v=np.asarray(ar[row['name']]);assert v.dtype.kind=='f' and v.ndim==0
                    actual=float(v)
                    claim=float(row['final'])
                    assert abs(claim-actual)<1e-12*max(1,abs(actual)),(name,row['name'],claim,actual)
                # the object array must exist; its bytes are never unpickled.
                assert 'correlations' in ar.files
            p={x['name']:float(x['final']) for x in params}
            typ={x['name']:x['type'] for x in params}
            param_rows.append({'archive':archive,'cohort':'main' if archive.startswith('mother') else 'replicate',
                               'file':Path(name).name,'model_params_npz_sha256':hashlib.sha256(zip_bytes).hexdigest(),
                               'all_11_numeric_parameters_agree_csv_with_npz':True,
                               'var_g':p['var_g'],'var_g_fit_role':typ['var_g'],'var_x':p['var_x'],
                               'var_dg':p['var_dg'],'beta':p['beta'],'mean_q':p['mean_q']})
params=pd.DataFrame(param_rows)
assert len(params)==31 and params['file'].is_unique
assert (params['var_g']>0).all()
params.to_csv(P/'E497_fit_parameters_31files.csv',index=False,float_format='%.15g')

# Novel diagnostic: raw observation vs same-raw-data posterior estimated daughter share.
f=D.fraction_first_raw_fluorescence.to_numpy(dtype=float)
g=D.fraction_first_model_fluorescence.to_numpy(dtype=float)
x=D.fraction_first_length.to_numpy(dtype=float)
assert np.all((f>0)&(f<1)) and np.all((g>0)&(g<1)) and np.all((x>0)&(x<1))
assert np.all((D.child1_model_cov_gg>=0)&(D.child2_model_cov_gg>=0))
D['model_minus_raw_split']=g-f
D['abs_model_minus_raw_split']=abs(g-f)
D['bright_daughter_label_changed']=(g>.5)!=(f>.5)
D['raw_margin_to_half']=abs(f-.5)
D['model_margin_to_half']=abs(g-.5)
D['relative_model_sigma_child1']=np.sqrt(D.child1_model_cov_gg)/D.child1_model_g
D['relative_model_sigma_child2']=np.sqrt(D.child2_model_cov_gg)/D.child2_model_g
assert np.isfinite(D[['relative_model_sigma_child1','relative_model_sigma_child2']]).all().all()

def desc(s):
    a=s.fraction_first_raw_fluorescence.to_numpy();b=s.fraction_first_model_fluorescence.to_numpy();v=s.fraction_first_length.to_numpy()
    return {'n':len(s),'n_dates':int(s.date.nunique()),
            'raw_vs_model_fraction_pearson_r':float(np.corrcoef(a,b)[0,1]),
            'mean_abs_raw_minus_model_fraction':float(np.abs(a-b).mean()),
            'median_abs_raw_minus_model_fraction':float(np.median(abs(a-b))),
            'n_bright_daughter_label_changed':int(s.bright_daughter_label_changed.sum()),
            'fraction_bright_daughter_label_changed':float(s.bright_daughter_label_changed.mean()),
            'length_vs_raw_fraction_r':float(np.corrcoef(v,a)[0,1]),
            'length_vs_posterior_fraction_r':float(np.corrcoef(v,b)[0,1]),
            'median_posterior_relative_sigma_child1':float(s.relative_model_sigma_child1.median()),
            'median_raw_margin_to_half_when_label_changes':float(s.loc[s.bright_daughter_label_changed,'raw_margin_to_half'].median()) if s.bright_daughter_label_changed.any() else None,
            'median_raw_margin_to_half_when_stable':float(s.loc[~s.bright_daughter_label_changed,'raw_margin_to_half'].median()) if (~s.bright_daughter_label_changed).any() else None}

# The 10-date blocked fitting compares two different targets, so DO NOT equate error scale with molecules.
def lodo(target):
    out=[]
    for date in sorted(D.date.unique()):
        tr=D[D.date!=date];te=D[D.date==date]
        xx=tr.fraction_first_length.to_numpy(dtype=float)
        yy=tr[target].to_numpy(dtype=float)
        alpha_slope=np.polyfit(xx,yy,1)
        slope=float(alpha_slope[0]);intercept=float(alpha_slope[1]);baseline=float(yy.mean())
        testx=te.fraction_first_length.to_numpy(dtype=float)
        testy=te[target].to_numpy(dtype=float)
        e0=float(((testy-baseline)**2).sum());e1=float(((testy-slope*testx-intercept)**2).sum())
        out.append({'date':str(date),'target':target,'n':len(te),'mse_base_sse':e0,'mse_size_sse':e1,
                    'improved':bool(e1<e0),'heldout_rel_MSE_reduction':1-e1/e0})
    return out
fitrows=lodo('fraction_first_raw_fluorescence')+lodo('fraction_first_model_fluorescence')
pd.DataFrame(fitrows).to_csv(P/'E497_cross_date_raw_vs_model_predictions.csv',index=False,float_format='%.15g')
byfile=D.groupby(['cohort','file'],sort=True).apply(desc,include_groups=False).reset_index(name='v')
byfile=pd.concat([byfile.drop(columns='v'),pd.json_normalize(byfile.v)],axis=1)
byfile=byfile.merge(params,on=['cohort','file'],validate='one_to_one')
assert len(byfile)==31
byfile.to_csv(P/'E497_by_file_observability.csv',index=False,float_format='%.15g')
# event-level new audit explicitly separate from already frozen E496 data.
COLS=['cohort','file','date','parent_id','child1_id','child2_id','parent_last_time','daughter_birth_time',
      'fraction_first_length','fraction_first_raw_fluorescence','fraction_first_model_fluorescence','model_minus_raw_split',
      'abs_model_minus_raw_split','bright_daughter_label_changed','raw_margin_to_half','model_margin_to_half',
      'relative_model_sigma_child1','relative_model_sigma_child2']
D[COLS].to_csv(P/'E497_raw_and_model_sibling_source_identity.csv',index=False,float_format='%.15g')
raw_fit=[r for r in fitrows if r['target']=='fraction_first_raw_fluorescence']
model_fit=[r for r in fitrows if r['target']=='fraction_first_model_fluorescence']
res={'project':901,'stage':'E497','source_2026_Kscheschinski_Zenodo':'10.5281/zenodo.18129208',
     'official_archives_sha256':RAW,'frozen_event_csv_sha256':EFILE,'n_strict_source_events':len(D),
     'source_parameter_31files_npz_11_numeric_scalars_checked':31*len(PARAM_NAMES),
     'npz_object_correlations_payload_deserialized':False,
     'var_g_estimation':'one optimized RealTrace likelihood parameter per file; NOT independently measured camera noise',
     'RealTrace_code_sources':{
       'setup_py':'https://github.com/nimwegenLab/RealTrace/blob/main/RealTrace_Tutorial/RealTrace_toolbox/setup.py',
       'setup_py_blob_sha':'f2e102e3dbbdde10ee22888643e038f63a606eb8',
       'likelihood_h':'https://github.com/nimwegenLab/RealTrace/blob/main/src/likelihood.h',
       'likelihood_h_blob_sha':'4d6bb32d7df405c247d1878b23eac67195a0380f',
       'predictions_h':'https://github.com/nimwegenLab/RealTrace/blob/main/src/predictions.h',
       'predictions_h_blob_sha':'3afc2aa20bf6916761b0be1d18defad7d09daeba'},
     'model_observation_noise_likelihood':'Depending on noise_model flag, D_gg=var_g or abs[var_g*(state_g+fp_auto)]; fit parameter in same source likelihood, not independent technical repeat',
     'model_posterior_cov_gg':'Kalman/posterior latent state conditional on source observations, same source data, NOT independent optical validation',
     'var_g_31file_summary':{'n':len(params),'median':float(params.var_g.median()),
                             'min':float(params.var_g.min()),'max':float(params.var_g.max()),
                             'models_free_count':int((params.var_g_fit_role=='free').sum()),
                             'models_fixed_count':int((params.var_g_fit_role=='fixed').sum())},
     'pooled':desc(D),'by_cohort':{str(name):desc(s) for name,s in D.groupby('cohort')},
     'by_10_date':{str(name):desc(s) for name,s in D.groupby('date')},
     'lodo_raw':{'total_base_sse':sum(q['mse_base_sse'] for q in raw_fit),'total_model_sse':sum(q['mse_size_sse'] for q in raw_fit),
                 'n_positive_dates':sum(q['improved'] for q in raw_fit)},
     'lodo_model_target':{'total_base_sse':sum(q['mse_base_sse'] for q in model_fit),'total_model_sse':sum(q['mse_size_sse'] for q in model_fit),
                          'n_positive_dates':sum(q['improved'] for q in model_fit)},
     'calibrated_absolute_molecule_count':'NOT_IDENTIFIED','physical_mother_daughters_MI_bits':'NOT_IDENTIFIED',
     'bit_per_second':'NOT_IDENTIFIED','bit_per_joule':'NOT_IDENTIFIED','publication':'PUBLICATION_NO_GO'}
for typ in ['lodo_raw','lodo_model_target']:
    z=res[typ];z['mse_base']=z['total_base_sse']/len(D);z['mse_size']=z['total_model_sse']/len(D);z['reduction_fraction']=1-z['total_model_sse']/z['total_base_sse']
(P/'E497_results.json').write_text(json.dumps(res,ensure_ascii=False,indent=2))
print(json.dumps({'params':res['var_g_31file_summary'],'pooled':res['pooled'],'by_cohort':res['by_cohort'],'raw_lodo':res['lodo_raw'],'posterior_lodo':res['lodo_model_target']},ensure_ascii=False,indent=2))
