#!/usr/bin/env python3
"""E497 separate source-original verification via Python standard library only.
Re-reads ALL 31 source CSVs and original NPZ float scalars; NEVER deserializes object arrays.
Independent population summaries, all 10058 source observation row keys, date-heldout.
"""
from pathlib import Path
from collections import defaultdict,Counter
import csv,io,json,zipfile,hashlib,struct,ast,statistics,math
P=Path(__file__).resolve().parent
SRC=P.parent/'p901_E496'
PARAM=['mean_lambda','gamma_lambda','var_lambda','mean_q','gamma_q','var_q','beta','var_x','var_g','var_dx','var_dg']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
RAW={'mother_machine_data_sets.zip':'a5787b437d9d136fd308feab5068fd558790f81929265444b42451287bc7b421',
     'replicates_hi1_rrnB.zip':'e6378d36486757c807f03a644e9653d131c24b00ba279472b2e512a0d43a1b6f'}
for name,h in RAW.items():assert sha(SRC/name)==h
# Exact E496 source parents, daughters, frozen lineage and observed model states.
rows=[];byfile=defaultdict(list);expect=defaultdict(dict)
with (SRC/'E496_2026_original_verified_sibling_events.csv').open(newline='',encoding='utf8') as f:
 for row in csv.DictReader(f):
  rows.append(row);key=(row['cohort'],row['file']);byfile[key].append(row)
  for name,time,prefix in [(row['parent_id'],row['parent_last_time'],'mother'),(row['child1_id'],row['daughter_birth_time'],'child1'),(row['child2_id'],row['daughter_birth_time'],'child2')]:
   key2=(name,float(time))
   if prefix=='mother':vals=(float(row['mother_raw_fp']),float(row['mother_model_g']),float(row['mother_model_cov_gg']))
   else:vals=(float(row[prefix+'_raw_fp']),float(row[prefix+'_model_g']),float(row[prefix+'_model_cov_gg']))
   v=expect[(row['cohort'],row['file'])].get(key2)
   if v is not None:assert all(abs(a-b)<=1e-11*max(1,abs(a)) for a,b in zip(vals,v)),(key,key2)
   expect[(row['cohort'],row['file'])][key2]=vals
assert len(rows)==10058 and len(byfile)==31

def read_float_npy(b):
 assert b[:6]==b'\x93NUMPY'
 major,minor=b[6],b[7]
 if major==1:
  n=struct.unpack('<H',b[8:10])[0];start=10
 elif major in (2,3):
  n=struct.unpack('<I',b[8:12])[0];start=12
 else:raise AssertionError(('unknown NPY version',major))
 header=ast.literal_eval(b[start:start+n].decode('latin-1'))
 assert header['shape']==() and header['fortran_order'] in (False,True)
 assert header['descr'] in ('<f8','|f8','=f8','float64')
 assert len(b)==start+n+8
 return struct.unpack('<d',b[-8:])[0]

allrowcount=0;original_observed_matches=0;all_params=0;var_g=[];unmatched=[]
for source in RAW:
 cohort='main' if source.startswith('mother_machine') else 'replicate'
 with zipfile.ZipFile(SRC/source) as z:
  assert z.testzip() is None
  preds=sorted(n for n in z.namelist() if n.endswith('_prediction.csv'))
  for path in preds:
   filename=Path(path).name;key=(cohort,filename)
   assert key in byfile
   rr=csv.reader(io.TextIOWrapper(z.open(path),newline='',encoding='utf8'))
   h=next(rr);assert h==['no','name','type','init','step','lower_bound','upper_bound','final']
   params=[dict(zip(h,next(rr))) for _ in range(11)]
   assert [x['name'] for x in params]==PARAM
   assert next(rr)==[]
   src_names=next(rr);assert src_names[0:5]==['cell_id','parent_id','time','log_length','fp']
   partner=path.replace('_autofl_corrected_','_rawdata_').replace('_prediction.csv','_correlations.npz')
   npzbytes=z.read(partner)
   with zipfile.ZipFile(io.BytesIO(npzbytes)) as nz:
    keys=set(x.removesuffix('.npy') for x in nz.namelist())
    assert keys==set(PARAM)|{'correlations'}
    for p in params:
     v=read_float_npy(nz.read(p['name']+'.npy'))
     assert abs(v-float(p['final']))<1e-12*max(1,abs(v))
     all_params+=1
    # DO NOT READ correlations.npy: may contain untrusted object/pickle payloads
    var_g.append(float(params[8]['final']))
   targets=expect[key];matches=set()
   for x in rr:
    if not x:continue
    allrowcount+=1
    pair=(x[0],float(x[2]))
    if pair in targets:
     obs=(float(x[4]),float(x[6]),float(x[13]))
     v=targets[pair]
     assert all(abs(a-b)<=1e-11*max(1,abs(b)) for a,b in zip(obs,v)),(key,pair,obs,v)
     original_observed_matches+=1;matches.add(pair)
   assert matches==set(targets),(filename,len(matches),len(targets))
assert len(var_g)==31 and all_params==341 and allrowcount==682679
# Implement independent score reductions with sufficient-statistic regression, without numpy.
def corr(x,y):
 n=len(x);a=sum(x);b=sum(y);aa=sum(t*t for t in x);bb=sum(t*t for t in y);ab=sum(u*v for u,v in zip(x,y))
 return (ab-a*b/n)/math.sqrt((aa-a*a/n)*(bb-b*b/n))
bydate=defaultdict(list)
f=[];g=[];X=[];n_flip=0
for r in rows:
 x=float(r['fraction_first_length']);z=float(r['fraction_first_raw_fluorescence']);w=float(r['fraction_first_model_fluorescence'])
 a=float(r['child1_raw_fp']);b=float(r['child2_raw_fp']);c=float(r['child1_model_g']);d=float(r['child2_model_g'])
 assert abs(z-a/(a+b))<1e-14 and abs(w-c/(c+d))<1e-14
 X.append(x);f.append(z);g.append(w);bydate[r['date']].append((x,z,w))
 n_flip+=((z>.5)!=(w>.5))
assert len(bydate)==10
res=json.loads((P/'E497_results.json').read_text())
assert abs(corr(f,g)-res['pooled']['raw_vs_model_fraction_pearson_r'])<1e-10
assert n_flip==res['pooled']['n_bright_daughter_label_changed']==923
assert abs(statistics.median([abs(a-b) for a,b in zip(f,g)])-res['pooled']['median_abs_raw_minus_model_fraction'])<1e-12
assert abs(statistics.median(var_g)-res['var_g_31file_summary']['median'])<1e-10

source_check={}
for idx,typ in [(1,'lodo_raw'),(2,'lodo_model_target')]:
 ss0=ss1=0.;good=0
 for date,tests in sorted(bydate.items()):
  train=[row for other,vals in bydate.items() if other!=date for row in vals]
  xs=[v[0] for v in train];ys=[v[idx] for v in train]
  n=len(xs);sx=sum(xs);sy=sum(ys);sxx=sum(x*x for x in xs);sxy=sum(a*b for a,b in zip(xs,ys))
  slope=(sxy-sx*sy/n)/(sxx-sx*sx/n)
  intercept=sy/n-slope*sx/n
  e0=sum((r[idx]-sy/n)**2 for r in tests);e1=sum((r[idx]-intercept-slope*r[0])**2 for r in tests)
  ss0+=e0;ss1+=e1;good+=int(e1<e0)
 want=res[typ]
 assert abs(ss0-want['total_base_sse'])<1e-7,(typ,ss0,want['total_base_sse'])
 assert abs(ss1-want['total_model_sse'])<1e-7,(typ,ss1,want['total_model_sse'])
 assert good==want['n_positive_dates']
 source_check[typ]={'source_date_groups':10,'positive_folds':good,'mse_base':ss0/len(rows),'mse_size':ss1/len(rows)}
qa={'project':901,'stage':'E497','passed':True,
    'author_zip_source_SHA256_verified':True,'original_csv_rows_reparsed':allrowcount,
    'unique_raw_source_cell_frame_keys_exactly_matched':original_observed_matches,
    'experiment_csv_files':len(byfile),'all_31_file_parameter_headings_confirmed':True,
    'source_NPZ_numeric_parameters_scalar_byte_reparsed':all_params,
    'NPZ_object_array_unpickled':False,
    'source_events_verified':len(rows),'independent_raw_vs_posterior_label_flips':n_flip,
    'independent_raw_vs_posterior_pearson':corr(f,g),'source_posterior_vs_raw_median_abs_split_difference':statistics.median([abs(a-b) for a,b in zip(f,g)]),
    'median_source_var_g':statistics.median(var_g),'independent_LODO_both_targets':source_check,
    'source_code_semantics':'Source var_g in RealTrace likelihood, cov_gg model posterior conditional on same source; not independent noise calibration',
    'source_inference':'NO_CALIBRATED_MOLECULAR_TRANSFER_BITS'}
(P/'E497_independent_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2))
print(json.dumps(qa,ensure_ascii=False,indent=2))
