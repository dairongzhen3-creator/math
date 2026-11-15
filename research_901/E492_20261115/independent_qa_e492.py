#!/usr/bin/env python3
"""Separate stdlib XML parser QA. Reconstruct row values from source JATS, not CSV input."""
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
import json, hashlib, xml.etree.ElementTree as ET, csv
P=Path(__file__).resolve().parent
x=ET.parse(P/'PMC8165573_fullTextXML.xml').getroot()
tables=[t for t in x.findall('.//table-wrap') if t.get('id')=='tbl1']
assert len(tables)==1
trs=tables[0].findall('.//tr')
original=[]
for tr in trs[1:]:
 vals=[' '.join(''.join(c.itertext()).split()) for c in tr if c.tag in ('th','td')]
 if vals:original.append(vals)
assert len(original)==145
by=defaultdict(dict)
for a in original:
 assert len(a)==6
 cell=int(a[0]);t=Decimal(a[1]);red=Decimal(a[4]);green=Decimal(a[5]);
 assert t not in by[cell];by[cell][t]=(red,green)
assert list(map(len,[by[1],by[2],by[3]]))==[61,50,34]
from_stat={}
for d in [2,3]:
 vals=[]
 for t in sorted(set(by[1])&set(by[d])):
  if d==2 and t>=min(by[3]):continue
  vals.append((t,by[1][t][0]/sum((by[1][t][0],by[d][t][0])),by[1][t][1]/sum((by[1][t][1],by[d][t][1]))))
 c=lambda a: (a[(len(a)-1)//2]+a[len(a)//2])/Decimal(2)
 from_stat[str(d)]={'n':len(vals),'median_red_fraction':str(c(sorted(q[1] for q in vals))),
 'median_green_fraction':str(c(sorted(q[2] for q in vals)))}
assert from_stat['2']['n']==16 and from_stat['3']['n']==34
res=json.loads((P/'E492_results.json').read_text())
for d in [2,3]:
 got=res['optical_descriptive_summaries'][str(d)]
 assert abs(Decimal(from_stat[str(d)]['median_red_fraction'])-Decimal(str(got['median_raw_red_mean_fraction_mother'])))<Decimal('1e-14')
 assert abs(Decimal(from_stat[str(d)]['median_green_fraction'])-Decimal(str(got['median_raw_green_mean_fraction_mother'])))<Decimal('1e-14')
with (P/'PMC8165573_table_tbl1.csv').open(newline='') as f:
 rows=list(csv.reader(f))
assert len(rows)==146 and rows[1:]==original
qa={'stage':'E492','source':'original verified PMC8165573 JATS XML Table 1 read independently of derived CSV',
'rows_exact_reparsed':len(original),'rows_match_csv':True,
'daughter_optical_pair_check':from_stat,
'claim':'source numerical QA only, not independently acquired biology and not physically calibrated molecule counts',
'pass':True}
(P/'E492_independent_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2))
print(json.dumps(qa,indent=2))
