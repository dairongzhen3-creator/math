#!/usr/bin/env python3
"""901 E498 independent all-source numerical audit.
Pure Python standard library; no pandas/numpy, and no imports from the primary.
Validates original Feather source SHA independently and source-row-preserving CSV gz.
This is independent computational implementation, NOT independent new experiment.
"""
from pathlib import Path
from collections import Counter,defaultdict
from decimal import Decimal,localcontext
import csv,hashlib,gzip,json,math,statistics
P=Path(__file__).resolve().parent
ORIG={
 'mydegrates.feather':'2a36a01facd1bf79ab0568e88cf07ff6cfebf099ecf9419ba85a0ead19cf7676',
 'myframes_degradation.feather':'e005351f269c85804c57f5f28d893f1db3f353ed139da4069a87820d491a4131',
 'myframes_bleaching.feather':'6e1733d6b641f5a47dca2851520fb9a0edfed2aad6faecbaf1ec5eece99917f7',
 'myproteomicsdata_multistrain.feather':'d0e83674a86d248f29cf7404a86251faca98b76cb4102e528ab3c4d48cc37386'}
REC=json.loads((P/'E498_original_to_csv_gzip_receipt.json').read_text())
for rec in REC['files']:
 name=rec['name'];src=P/name;csvname=P/(name+'.csv.gz')
 assert hashlib.sha256(src.read_bytes()).hexdigest()==ORIG[name]
 assert hashlib.sha256(csvname.read_bytes()).hexdigest()==rec['gzip_csv_sha256']
source=json.loads((P/'E498_results.json').read_text())
with gzip.open(P/'mydegrates.feather.csv.gz','rt',newline='') as fh:author=list(csv.DictReader(fh))
assert len(author)==4
with localcontext() as ctx:
 ctx.prec=48
 for a in author:
  g=Decimal(a['gamma_init']);w=Decimal(a['gamma_dec_Fluo'])
  ref=next(r for r in source['author_source_degradation_model']['all_four_reporters'] if r['Fluo1']==a['Fluo1'])
  for hour in (5,10,30,60):
   t=Decimal(hour);res=(-(g/w)*(1-(-w*t).exp())).exp()
   assert abs(res-Decimal(str(ref[f'model_old_fluorescent_label_survival_{hour}h'])))<Decimal('1e-13')

# Read all 1,116,240 genuine source spectral rows and construct only the actual matching reporter.
map_f={'pTrc99a_mVenus':'YFP','pTrc99a_mCerulean':'CFP','pTrc99a_mGFP':'GFP','pTrc99a_mCherry2':'RFP'}
frames=Counter();expected=defaultdict(list);nrows=0;nmatched=0;nstarved=0;npositive=0
with gzip.open(P/'myframes_degradation.feather.csv.gz','rt',newline='') as fh:
 for r in csv.DictReader(fh):
  nrows+=1;frames[(r['date'],r['cell'],r['time_sec'])]+=1
  if r['Fluo']!=map_f[r['strain']]:continue
  nmatched+=1
  if r['c_step']!='SPM.1':continue
  nstarved+=1
  if not r['Ampl'] or float(r['Ampl'])<=0 or float(r['time_sec'])<=float(r['m_start'])+3600:continue
  npositive+=1
  expected[(r['date'],r['strain'],r['Fluo'],r['cell'])].append((float(r['time_sec']),float(r['m_start']),float(r['Ampl'])))
assert nrows==1116240 and len(frames)==279060 and all(v==4 for v in frames.values())
assert nmatched==279060 and nstarved==122974 and npositive==3992
n_intervals=0;all_cells=set();percell=defaultdict(lambda:defaultdict(list))
for key,values in expected.items():
 rows=sorted(values)
 for (t0,m0,a0),(t1,m1,a1) in zip(rows[:-1],rows[1:]):
  dt=(t1-t0)/3600.
  if not 2<=dt<=7:continue
  n_intervals+=1;all_cells.add(key[3]);h=(t0+t1-2*m1)/7200.
  rate=-math.log(a1/a0)/dt
  window='3_12' if 3<=h<12 else '18_30' if 18<=h<30 else None
  if window:percell[key][window].append(rate)
assert n_intervals==3517 and len(all_cells)==439
paired=[]
for key,parts in percell.items():
 if '3_12' in parts and '18_30' in parts:
  a=statistics.median(parts['3_12']);b=statistics.median(parts['18_30']);paired.append((*key,a,b,a>b))
assert len(paired)==403 and sum(row[-1] for row in paired)==368
pg=defaultdict(list)
for date,strain,fluor,cell,early,late,gt in paired:pg[(date,strain,fluor)].append((early,late,gt))
prim={(str(z['date']),z['strain'],z['Fluo']):z for z in source['degradation']['per_group_paired']}
assert set(pg)==set(prim)
for k,v in pg.items():
 z=prim[k];assert len(v)==int(z['n_cells'])
 assert sum(bool(q[-1]) for q in v)==int(z['n_early_greater'])
 assert abs(statistics.median(q[0]-q[1] for q in v)-z['median_within_cell_early_minus_late_per_h'])<1e-12

# Independently calculate weighted apparent bleaching slopes per source cell phase.
bleach=defaultdict(list);bleach_rows=0
with gzip.open(P/'myframes_bleaching.feather.csv.gz','rt',newline='') as fh:
 for r in csv.DictReader(fh):
  bleach_rows+=1
  key=(r['date'],r['exp_series'],int(float(r['idx'])),r['cell'],r['Fluo_sig'])
  bleach[key].append((float(r['time_sec'])/3600.,float(r['signal'])))
assert bleach_rows==91838
bleach_per=defaultdict(list);bleach_count=0;bleach_pos=0
for (date,ser,idx,cell,f),records in bleach.items():
 if len(records)<=5:continue
 t=[v[0] for v in records];y=[math.log(v[1]) for v in records];w=[1/v[1] for v in records]
 W=sum(w);mt=sum(a*b for a,b in zip(w,t))/W;my=sum(a*b for a,b in zip(w,y))/W
 v_t=sum(wi*(ti-mt)**2 for wi,ti in zip(w,t))/W
 if v_t<=0:continue
 cov=sum(wi*(ti-mt)*(yi-my) for wi,ti,yi in zip(w,t,y))/W
 alpha=-cov/v_t
 bleach_per[(date,f,idx)].append(alpha)
 bleach_count+=1;bleach_pos+=int(alpha>0)
assert bleach_count==1659 and bleach_pos==1446
bref={(str(q['date']),q['Fluo_sig'],int(q['phase_idx'])):q for q in source['bleaching']['by_group']}
for k,vals in bleach_per.items():
 z=bref[k];assert len(vals)==int(z['n_cells'])
 assert sum(x>0 for x in vals)==int(z['n_positive'])
 assert abs(statistics.median(vals)-z['median_all_source_rate_per_h'])<1e-11
 assert abs(statistics.median(x for x in vals if x>0)-z['median_positive_only_per_h'])<1e-11

# Independent proteomics source dimension, missingness and median cross-background spread.
rawcount=0;missing=0;genes=set();proms=set();by=defaultdict(dict)
with gzip.open(P/'myproteomicsdata_multistrain.feather.csv.gz','rt',newline='') as fh:
 for r in csv.DictReader(fh):
  rawcount+=1;genes.add(r['Genes']);proms.add(r['GFP_prom'])
  key=(r['Genes'],r['GFP_prom'])
  if not r['Log_Quantity_norm'] or r['Log_Quantity_norm'].lower()=='nan':missing+=1;continue
  sample=int(float(r['sample']));by[key][sample]=float(r['Log_Quantity_norm'])
assert rawcount==130200 and missing==16543 and len(genes)==2604 and len(proms)==10
across=defaultdict(list)
for (gene,prom),row in by.items():
 pre=[row[s] for s in (1,2) if s in row]
 if not pre:continue
 start=sum(pre)/len(pre)
 for sample,time in [(3,1.4),(4,21.4),(5,45.4)]:
  if sample in row:across[(gene,time)].append((row[sample]-start)*math.log(2.))
groups={}
for time in (1.4,21.4,45.4):
 values=[max(v)-min(v) for (gene,t),v in across.items() if t==time and len(v)>=5]
 total=sum(len(v)==10 for (gene,t),v in across.items() if t==time)
 groups[time]={'count5':len(values),'count10':total,'median_span':statistics.median(values)}
 rec=next(r for r in source['proteomics']['per_source_time'] if r['time_h']==time)
 assert len(values)==rec['source_genes_with_at_least_5_strain_backgrounds']
 assert total==rec['source_genes_with_all_10_strain_backgrounds']
 assert abs(statistics.median(values)-rec['median_across_backgrounds_gene_span_ln'])<1e-11

qa={'project':901,'stage':'E498','PASS':True,
    'source_original_4_Feather_SHA256_pass':True,'source_export_4_gzip_SHA256_pass':True,
    'source_degradation_rows':nrows,'source_unique_cell_time':len(frames),
    'source_spectral_rows_reporter_matched':nmatched,'source_starvation_rows':nstarved,
    'source_positive_post1h_starvation_rows':npositive,'valid_2to7h_degradation_intervals':n_intervals,
    'source_unique_cells_in_intervals':len(all_cells),'paired_early_late_within_cell':len(paired),
    'paired_early_greater':sum(row[-1] for row in paired),
    'source_bleaching_rows':bleach_rows,'independent_WLS_cells_n_gt5':bleach_count,
    'positive_WLS_bleaching_cells':bleach_pos,
    'source_proteomics_rows':rawcount,'source_proteomics_genes':len(genes),
    'source_proteomics_backgrounds':len(proms),'source_proteomics_missing_logs':missing,
    'source_proteomics_per_time_cross_background_span':groups,
    'author_gamma_survival_model_decimal48_check':'PASS',
    'source_note':'Numerical source audit using independently implemented stdlib streaming from verified original exported rows, not independent wet lab or protein count calibration',
    'physical_molecular_MI':'NOT_IDENTIFIED'}
(P/'E498_independent_QA.json').write_text(json.dumps(qa,indent=2,ensure_ascii=False),encoding='utf8')
print(json.dumps(qa,ensure_ascii=False,indent=2))
