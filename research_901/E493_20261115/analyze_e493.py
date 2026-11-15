#!/usr/bin/env python3
"""901 E493: official TrIPP isotope-label half-life table and three-cell optical contrast.
Scientific scope: measured population heavy/light SiLAC inputs, author-fitted
positive degradation half-lives, and single-lineage raw spectral contrast.
No inference of cross-division absolute molecule counts or transmitted bits.
"""
from pathlib import Path
from artifact_tool import Blob, SpreadsheetFile
import csv, hashlib, json, math, statistics, collections, xml.etree.ElementTree as ET

P=Path(__file__).resolve().parent
SILAC=P/'PMC7889982_TableS2_SILAC_original_mmc2.xlsx'
XML=P/'PMC8165573_fullTextXML.xml'
EXPECTED_SILAC_SHA='2ea57bedff6ccdd71ea982e695eea638cb6a815a1d2c8347b7ae91d8d8c21893'
EXPECTED_SILAC_MD5='4d2104fd2a40a51bf21b8dddfc0f3428'
EXPECTED_XML_SHA='6a59cfb29e1e933577a40b99de88e5ee4ed7cdf1e6b0355e371d8a08fd9b6667'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
assert len(SILAC.read_bytes())==842092 and sha(SILAC)==EXPECTED_SILAC_SHA
assert hashlib.md5(SILAC.read_bytes()).hexdigest()==EXPECTED_SILAC_MD5
assert sha(XML)==EXPECTED_XML_SHA

wb=SpreadsheetFile.import_xlsx(Blob.load(str(SILAC)))
s=wb.worksheets.get_item('TableS2-Half-lives-SilAC')
raw=s.get_range('A4:BC1341').values
assert len(raw)==1338 and all(len(x)==55 for x in raw)
head=s.get_range('A3:BC3').values[0]
assert head[0]=='Protein IDs' and head[51]=='half life (min)' and head[54]=='half life (min)'

def num(x):
 try:
  z=float(x)
  return z if math.isfinite(z) else None
 except (TypeError,ValueError):return None

rows=[]
for r in raw:
 obj={'protein_id':str(r[0]),'gene':str(r[3] or ''),'name':str(r[2] or ''),
      'unique_peptides':num(r[6]),
      'cbf1_slope_per_min':num(r[49]),'cbf1_kd_per_min':num(r[50]),
      'cbf1_half_life_min':num(r[51]),
      'hho1_slope_per_min':num(r[52]),'hho1_kd_per_min':num(r[53]),
      'hho1_half_life_min':num(r[54])}
 for context,base in [('cbf1',7),('hho1',14)]:
  for j,t in enumerate([0,30,60,90,120,150,180]):
   obj[f'{context}_source_HL_t{t}min']=num(r[base+j])
 rows.append(obj)
assert len(set(r['protein_id'] for r in rows))==1338
outcsv=P/'E493_SILAC_source_locked_1338.csv'
with outcsv.open('w',newline='',encoding='utf-8') as fh:
 w=csv.DictWriter(fh,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

statuses={}
for c in ('cbf1','hho1'):
 kd=c+'_kd_per_min';h=c+'_half_life_min';slope=c+'_slope_per_min'
 valid=[r for r in rows if r[h] is not None and r[h]>0 and r[kd] is not None and r[kd]>0]
 invalid=[{'protein_id':r['protein_id'],'gene':r['gene'],h:r[h],kd:r[kd]} for r in rows if r not in valid]
 discrepancy=max(abs(r[h]*r[kd]-math.log(2)) for r in rows if r[kd] not in (None,0) and r[h] is not None)
 offset=[r[slope]-r[kd] for r in rows if r[slope] is not None and r[kd] is not None]
 assert max(offset)-min(offset)<1e-11
 statuses[c]={'valid_positive_count':len(valid),'nonpositive_or_invalid_count':len(invalid),'invalids':invalid,
             'fitted_slope_minus_kd_offset_per_min':statistics.median(offset),
             'offset_range_maxmin':max(offset)-min(offset),'max_half_life_times_kd_ln2_residual':discrepancy,
             'positive_half_life_median_min':statistics.median(r[h] for r in valid)}

paired=[r for r in rows if all(r[c+'_half_life_min'] is not None and r[c+'_half_life_min']>0 and r[c+'_kd_per_min'] is not None and r[c+'_kd_per_min']>0 for c in ('cbf1','hho1'))]
assert len(paired)==1335
x=[r['cbf1_half_life_min'] for r in paired]
y=[r['hho1_half_life_min'] for r in paired]
ratio=[b/a for a,b in zip(x,y)]

def ranks(values):
 indexed=sorted((v,i) for i,v in enumerate(values));out=[None]*len(values);j=0
 while j<len(values):
  k=j+1
  while k<len(values) and indexed[k][0]==indexed[j][0]:k+=1
  rank=((j+1)+k)/2
  for _,ix in indexed[j:k]:out[ix]=rank
  j=k
 return out

def pearson(a,b):
 aa=statistics.mean(a);bb=statistics.mean(b)
 za=[t-aa for t in a];zb=[t-bb for t in b]
 return sum(p*q for p,q in zip(za,zb))/math.sqrt(sum(p*p for p in za)*sum(q*q for q in zb))

quant=lambda v,p: sorted(v)[0] if p==0 else __import__('numpy').quantile(v,p).item()
summary={'stage':'E493','source':'Auboiron et al iScience 2021 DOI 10.1016/j.isci.2021.102075',
 'TableS2_source_sha256':EXPECTED_SILAC_SHA,'TableS2_source_md5':EXPECTED_SILAC_MD5,
 'TableS2_original_bytes':SILAC.stat().st_size,'TableS2_original_protein_rows':len(rows),
 'experimental_contexts':['CBF1 Dendra2 strain','HHO1 Dendra2 strain'],
 'raw_ratios_time_grid_minutes':[0,30,60,90,120,150,180],
 'source_author_fitted_results_not_new_raw_protein_copy_counts':True,
 'source_fit_validity_by_context':statuses,'paired_positive_half_life_count':len(paired),
 'paired_half_life_median_min':{'CBF1':statistics.median(x),'HHO1':statistics.median(y)},
 'paired_half_life_quartiles_min':{
    'CBF1':[quant(x,.25),quant(x,.75)],'HHO1':[quant(y,.25),quant(y,.75)]},
 'paired_source_halflife_spearman':pearson(ranks(x),ranks(y)),
 'paired_source_halflife_pearson_log':pearson([math.log(v) for v in x],[math.log(v) for v in y]),
 'fraction_positive_pairs_hho1_half_life_shorter':sum(b<a for a,b in zip(x,y))/len(paired),
 'paired_hho1_to_cbf1_halflife_median_ratio':statistics.median(ratio),
 'fraction_paired_within_twofold':sum(.5<=r<=2 for r in ratio)/len(ratio),
 'fraction_paired_gt_twofold_divergence':sum(r<.5 or r>2 for r in ratio)/len(ratio),
 'under_fitted_exponential_model_median_old_label_survival_at_120min':{
  'CBF1':statistics.median(math.exp(-math.log(2)*120/a) for a in x),
  'HHO1':statistics.median(math.exp(-math.log(2)*120/b) for b in y)},
 'data_interpretation':'Measured population pulse-SiLAC heavy/light ratios and derived apparent protein degradation half-lives; no per-mother daughter pairing. No assumption of two biological replicates.',
 'molecular_mother_daughter_bit':'NOT_IDENTIFIED'}

# New optical double-ratio avoids assuming cell ROI size or total fluorescence,
# but does NOT remove an additive background, spectral crosstalk, or photobleach.
root=ET.parse(XML).getroot()
tables=[t for t in root.findall('.//table-wrap') if t.get('id')=='tbl1']
assert len(tables)==1
tr=tables[0].findall('.//tr')
obs={}
for row in tr[1:]:
 a=[' '.join(''.join(c.itertext()).split()) for c in row if c.tag in ('td','th')]
 if not a:continue
 cell=int(a[0]);t=float(a[1]);red=float(a[4]);green=float(a[5]);
 assert red>0 and green>0
 obs.setdefault(cell,{})[t]={'red':red,'green':green}
assert list(map(lambda c:len(obs[c]),(1,2,3)))==[61,50,34]
all_t=sorted(set(obs[1])&set(obs[2])&set(obs[3]));assert len(all_t)==34
pair_rows=[]
for t in sorted(set(obs[1])&set(obs[2])):
 if t>=min(obs[3]):continue
 m,d=obs[1][t],obs[2][t]
 q=math.log2((m['red']/m['green'])/(d['red']/d['green']))
 pair_rows.append({'time_min':t,'daughter':2,'same_frame_contrast_log2_RG_mother_to_daughter':q,'red_mean_m':m['red'],'green_mean_m':m['green'],'red_mean_d':d['red'],'green_mean_d':d['green']})
assert len(pair_rows)==16
for t in all_t:
 for c in (2,3):
  m,d=obs[1][t],obs[c][t]
  q=math.log2((m['red']/m['green'])/(d['red']/d['green']))
  pair_rows.append({'time_min':t,'daughter':c,'same_frame_contrast_log2_RG_mother_to_daughter':q,
                    'red_mean_m':m['red'],'green_mean_m':m['green'],'red_mean_d':d['red'],'green_mean_d':d['green']})
with (P/'E493_TrIPP_same_frame_optical_double_ratio.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=list(pair_rows[0]));w.writeheader();w.writerows(pair_rows)
summary['tripp_single_lineage_optical_double_ratio']={
 'source_xml_sha256':EXPECTED_XML_SHA,'physical_lineages':1,'matched_three_cell_frames':len(all_t),
 'matched_total_m_daughter_frame_pairs':len(pair_rows),
 'metric':'log2[(raw_red_m/raw_green_m)/(raw_red_d/raw_green_d)]',
 'interpretation':'A dimensionless raw spectral contrast at matched time; not per-cell old/new molecule or survival ratio',
 'source_window_groups':{}}
for label,subset in [('daughter1_before_daughter2',[r for r in pair_rows if r['daughter']==2 and r['time_min']<229.5]),
 ('daughter1_after_daughter2',[r for r in pair_rows if r['daughter']==2 and r['time_min']>=229.5]),
 ('daughter2_after_birth',[r for r in pair_rows if r['daughter']==3])]:
 z=[r['same_frame_contrast_log2_RG_mother_to_daughter'] for r in subset]
 summary['tripp_single_lineage_optical_double_ratio']['source_window_groups'][label]={'n_matched_frames':len(z),'median_log2_spectral_contrast':statistics.median(z),'positive_count':sum(v>0 for v in z),'median_positive_ratio':2**statistics.median(z)}
summary['tripp_single_lineage_optical_double_ratio']['identification_theorem']='With only raw per-cell mean intensities I_m,I_d>0, even if background=0 and gains equal, total mother share I_m/(I_m+A_d/A_m*I_d) takes every value strictly between 0 and 1 as unknown positive ROI area ratio varies. Thus no nontrivial molecular fraction bound follows.'
summary['tripp_single_lineage_optical_double_ratio']['unresolved']=['ROI area','raw photon background','photoconversion efficiency','per-cell gain','absolute molecule number','technical repeats','n independent lineages']
summary['study_public_data_statement']='All microscopy images and Perl and R scripts available upon request from Lead Contact; 18-protein BudJ raw CSV was not exposed by E492 public materials'
summary['publication']='PUBLICATION_NO_GO'
(P/'E493_results.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({k:v for k,v in summary.items() if k in ('paired_positive_half_life_count','paired_half_life_median_min','paired_source_halflife_spearman','fraction_positive_pairs_hho1_half_life_shorter','tripp_single_lineage_optical_double_ratio')},indent=2))
