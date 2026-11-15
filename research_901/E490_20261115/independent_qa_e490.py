#!/usr/bin/env python3
"""Independent E490 rule verification directly from the original publisher Git LFS
workbook with artifact_tool. No primary CSV-derived calculation is imported.
Does not certify optical calibration, false-negative rate, or physical molecules.
"""
from artifact_tool import Blob, SpreadsheetFile
from pathlib import Path
from collections import Counter
from decimal import Decimal, localcontext
import hashlib, json, csv, math

ROOT=Path(__file__).resolve().parent
name='FLCS_lacInduction.xlsx'
expected='97f93074b72214bd8a20439ff508f6d8d47111ac115b6be6c06538457731212d'
blob=(ROOT/name).read_bytes()
assert len(blob)==33583 and hashlib.sha256(blob).hexdigest()==expected
book=SpreadsheetFile.import_xlsx(Blob.load(str(ROOT/name)))
data=book.worksheets.get_item('MotherMachine').get_range('A2:P75').values
assert len(data)==74
# Author code pin: .33 um probing width; visible plot dotted threshold 2.55,
# fitted number <20; diffusion time in [3,30] ms.
num=lambda x: float(x) if isinstance(x,(int,float)) and math.isfinite(float(x)) else None
median_candidates=[num(v[14]) for v in data]
median_candidates=sorted(x for x in median_candidates if x is not None)
assert len(median_candidates)==66
median=(median_candidates[32]+median_candidates[33])/2
assert median==17.5
observed=Counter();recon=Counter();mismatched=[];all_complete=[];robust=[]
for idx,row in enumerate(data,start=2):
  d,filename=(str(row[0]),str(row[1]))
  label=str(row[2]) if row[2] is not None else 'NA'
  # author workbook column indices K-L fitted N, M-N D, O-P CPM
  n,nse,D,Dse,cpm,cpmse=[num(row[i]) for i in (10,11,12,13,14,15)]
  if label=='-1.0':label='-1'
  if label=='1.0':label='1'
  if label=='0.0':label='0'
  observed[label]+=1
  if any(x is None for x in (n,D,cpm)):
    recon['INCOMPLETE']+=1
    continue
  tau=1000 * Decimal('0.33')**2/(Decimal(4)*Decimal(str(D))) if D>0 else None
  high=(cpm/median)>2.55
  n_small=n<20
  valid_diff=bool(tau is not None and Decimal('3')<=tau<=Decimal('30'))
  flag='1' if high and n_small and valid_diff else '-1' if high else '0'
  recon[flag]+=1
  all_complete.append((d,filename,label,flag))
  if label!=flag:mismatched.append((idx,filename,label,flag))
  if flag=='1':
    # No confidence statements implied: sensitivity of frozen rule to one reported fit SE.
    c_ok=cpmse is not None and cpm-cpmse>2.55*median
    n_ok=nse is not None and n+nse<20
    d_ok=bool(Dse is not None and D-Dse>0 and
              3<=1000*.33**2/(4*(D+Dse)) and 1000*.33**2/(4*(D-Dse))<=30)
    robust.append((d,filename,n_ok,c_ok,d_ok))
assert observed==Counter({'0':45,'1':18,'-1':5,'NA':6}),observed
assert len(all_complete)==66 and not mismatched
assert recon==Counter({'0':43,'1':18,'-1':5,'INCOMPLETE':8}),recon
assert len(robust)==18
assert sum(x[2] for x in robust)==18 and sum(x[3] for x in robust)==18 and sum(x[4] for x in robust)==16
# Check main per-source receipt: QA's author XLSX and primary CSV identities agree.
with (ROOT/'E490_author_rule_per_cell.csv').open(encoding='utf-8',newline='') as f:
   baseline={(x['date'],x['ptu_filename_reference']):x for x in csv.DictReader(f)}
assert len(baseline)==74
assert all(baseline[d,filename]['published_gfp_label']==label for d,filename,label,flag in all_complete)
assert all(baseline[d,filename]['three_rule_reconstruction']==flag for d,filename,label,flag in all_complete)
report={
 'stage':'901_E490','QA_source':'original Git LFS 33583-byte XLSX directly imported through artifact_tool',
 'sha256_xlsx':hashlib.sha256(blob).hexdigest(),'artifact_tool_original_sheets_reparsed':True,
 'source_records':len(data),'source_categories':dict(observed),
 'median_CPM':median,'reconstructed_fully_audited':len(all_complete),
 'agreement':len(all_complete)-len(mismatched),'mismatches':mismatched,
 'source_fit_incomplete':len(data)-len(all_complete),
 'reconstructed_by_class':dict(recon),
 'reported_SE_diagnostic':{'18_positive_particles_margin':sum(x[2] for x in robust),
                           '18_positive_CPM_margin':sum(x[3] for x in robust),
                           '18_positive_diffusion_margin':sum(x[4] for x in robust)},
 'independent_sources_for_false_negative_rate':0,
 'instrumental_molecule_copy_noise_independently_identified':False,
 'method':'Independent artifact_tool XLSX decode + separate Decimal diffusion rule + 74 row comparison; shares same authored workbook but not primary calculation',
 'pass':True
}
(ROOT/'E490_independent_QA.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(report,indent=2,ensure_ascii=False))
