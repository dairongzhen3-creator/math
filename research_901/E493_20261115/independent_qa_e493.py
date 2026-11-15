#!/usr/bin/env python3
"""Independent E493 source verification: standard-library OOXML reader (no artifact_tool)
+ SciPy rank-correlation and Decimal direct original-JATS spectral contrast.
Uses original published bytes rather than E493's CSV export.
"""
from pathlib import Path
import zipfile,xml.etree.ElementTree as ET,statistics,hashlib,json,math,decimal
from collections import Counter
from scipy import stats
P=Path(__file__).resolve().parent
p=P/'PMC7889982_TableS2_SILAC_original_mmc2.xlsx'
xml=P/'PMC8165573_fullTextXML.xml'
assert len(p.read_bytes())==842092
assert hashlib.md5(p.read_bytes()).hexdigest()=='4d2104fd2a40a51bf21b8dddfc0f3428'
assert hashlib.sha256(p.read_bytes()).hexdigest()=='2ea57bedff6ccdd71ea982e695eea638cb6a815a1d2c8347b7ae91d8d8c21893'
assert hashlib.sha256(xml.read_bytes()).hexdigest()=='6a59cfb29e1e933577a40b99de88e5ee4ed7cdf1e6b0355e371d8a08fd9b6667'
ns='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
with zipfile.ZipFile(p) as z:
 assert z.testzip() is None
 sh=ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
 strings=ET.fromstring(z.read('xl/sharedStrings.xml'))
 sst=[''.join(e.text or '' for e in si.iter(ns+'t')) for si in strings.findall(ns+'si')]
 def cellval(e):
  v=e.find(ns+'v')
  if v is None:return None
  if e.get('t')=='s':return sst[int(v.text)]
  return float(v.text)
 selected=[]
 for r in sh.iter(ns+'row'):
  rn=int(r.get('r'))
  if not 4<=rn<=1341:continue
  columns={''.join(c for c in cell.get('r') if c.isalpha()):cellval(cell) for cell in r.findall(ns+'c')}
  selected.append((columns.get('A'),columns.get('D'),columns.get('AY'),columns.get('AZ'),columns.get('BB'),columns.get('BC'),columns.get('AX'),columns.get('BA')))
assert len(selected)==1338
assert len(set(v[0] for v in selected))==1338
cbf=[v[3] for v in selected if v[3] is not None and v[3]>0 and v[2] is not None and v[2]>0]
hho=[v[5] for v in selected if v[5] is not None and v[5]>0 and v[4] is not None and v[4]>0]
paired=[v for v in selected if v[3] is not None and v[3]>0 and v[2] is not None and v[2]>0 and v[5] is not None and v[5]>0 and v[4] is not None and v[4]>0]
assert (len(cbf),len(hho),len(paired))==(1337,1336,1335)
negative=[(str(v[1]),float(v[2]),float(v[4])) for v in selected if v[2]<0 or v[4]<0]
assert {v[0] for v in negative}=={'DBP8','LRP1','GYP5'}
for v in selected:
 assert abs(v[2]*v[3]-math.log(2))<1e-10
 assert abs(v[4]*v[5]-math.log(2))<1e-10

X=[v[3] for v in paired];Y=[v[5] for v in paired]
rho=float(stats.spearmanr(X,Y).statistic)
summary=json.loads((P/'E493_results.json').read_text())
assert abs(rho-summary['paired_source_halflife_spearman'])<1e-13
assert abs(statistics.median(X)-summary['paired_half_life_median_min']['CBF1'])<1e-12
assert abs(statistics.median(Y)-summary['paired_half_life_median_min']['HHO1'])<1e-12
assert abs(sum(b<a for a,b in zip(X,Y))/len(X)-summary['fraction_positive_pairs_hho1_half_life_shorter'])<1e-12

with decimal.localcontext() as ctx:
 ctx.prec=48
 root=ET.parse(xml).getroot()
 tbl=[t for t in root.findall('.//table-wrap') if t.get('id')=='tbl1']
 assert len(tbl)==1
 values={}
 for row in tbl[0].findall('.//tr')[1:]:
  a=[' '.join(''.join(c.itertext()).split()) for c in row if c.tag in ('td','th')]
  if not a:continue
  c=int(a[0]);t=decimal.Decimal(a[1]);red=decimal.Decimal(a[4]);green=decimal.Decimal(a[5]);
  values.setdefault(c,{})[t]=(red,green)
 assert [len(values[c]) for c in [1,2,3]]==[61,50,34]
 post=sorted(set(values[1])&set(values[2])&set(values[3]))
 assert len(post)==34
 d1pre=sorted((set(values[1])&set(values[2]))-set(post))
 assert len(d1pre)==16
 def contrast(m,d):
  rM,gM=m;rD,gD=d
  return ((rM*gD)/(gM*rD)).ln()/decimal.Decimal(2).ln()
 def med(xs):return statistics.median(xs)
 original={}
 for label,ts,dc in [('daughter1_before_daughter2',d1pre,2),('daughter1_after_daughter2',post,2),('daughter2_after_birth',post,3)]:
  z=[contrast(values[1][t],values[dc][t]) for t in ts]
  original[label]={'n':len(z),'median':float(med(z)),'count_positive':sum(t>0 for t in z)}
  prev=summary['tripp_single_lineage_optical_double_ratio']['source_window_groups'][label]
  assert original[label]['n']==prev['n_matched_frames']
  assert original[label]['count_positive']==prev['positive_count']
  assert abs(original[label]['median']-prev['median_log2_spectral_contrast'])<1e-13

qa={'stage':'E493','source_independent':True,
 'source_SILAC':'Original exact PMC official mmc2.xlsx byte SHA and MD5; directly re-parsed ZIP+XML without artifact_tool',
 'source_TrIPP':'Source JATS Table1 directly re-read using Decimal48',
 'protein_entries':len(selected),'CBF1_positive':len(cbf),'HHO1_positive':len(hho),'paired_positive':len(paired),
 'nonpositive_fitted_ids':negative,'median_CBF1_min':statistics.median(X),'median_HHO1_min':statistics.median(Y),
 'sciPy_spearman':rho,'optical_source_crosscheck':original,
 'scientific_scope':'source numerical verification, NOT a biological replicate, NOT per-cell molecule count calibration',
 'QA_PASS':True}
(P/'E493_independent_QA.json').write_text(json.dumps(qa,indent=2,ensure_ascii=False))
print(json.dumps(qa,indent=2,ensure_ascii=False))
