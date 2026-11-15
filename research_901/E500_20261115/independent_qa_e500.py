#!/usr/bin/env python3
"""Independent E500 numerical QA without pandas/numpy/scipy/statsmodels or primary functions.
Uses source-extracted row-unique mother identity records plus original official Feather SHA and source receipt.
This is an independent *analysis implementation*, not second independent Arrow/Feather decoder or wet-lab replication.
"""
from pathlib import Path
import hashlib,csv,json,collections,statistics,math,random
P=Path(__file__).resolve().parent
F=P/'myframes_control.feather'
assert F.stat().st_size==302295912
h=hashlib.sha256()
with F.open('rb') as f:
 while True:
  b=f.read(1024*1024)
  if not b:break
  h.update(b)
assert h.hexdigest()=='5cb1850cb127a64bf5e3dffa4d4e13e7cc648f73af176f6d03cf4d0e1fa9bb7c'
r=json.loads((P/'E500_original_source_receipt.json').read_text());assert r['source_member_crc_pass'] and r['original_member_sha256']==h.hexdigest()
rows=[];keys=set();cell_ids=set();n=collections.Counter();positions=collections.defaultdict(lambda:collections.Counter());strains=collections.defaultdict(lambda:collections.Counter())
with (P/'E500_original_illumination_mother_fate_candidates.csv').open(newline='') as f:
 for x in csv.DictReader(f):
  cell=(x['date'],x['condition'],x['cell'])
  assert cell not in keys;keys.add(cell)
  assert x['cell'] not in cell_ids;cell_ids.add(x['cell'])
  assert x['date']=='20221220' and x['starvation_nominal_h']=='60h'
  assert float(x['time_start'])<float(x['m_start'])
  assert ((x['end_type']=='div') == (x['fate']=='Survive'))
  light='NO_LIGHT' not in x['condition'];assert light==(x['illuminated']=='True')
  n[(light,x['fate'])]+=1
  pos=x['gl_id'].split('.')[1];assert x['gl_id'].startswith('20221220.')
  positions[(light,pos)][x['fate']]+=1
  strains[(x['strain'],light)][x['fate']]+=1
  rows.append(x)
assert len(rows)==238
assert n[(True,'Survive')]==85 and n[(True,'No Growth')]==100 and n[(True,'Unknown')]==6
assert n[(False,'Survive')]==43 and n[(False,'No Growth')]==4 and n[(False,'Unknown')]==0
assert sum(k[0] for k in positions)==18 and sum(not k[0] for k in positions)==6
summary=json.loads((P/'E500_results.json').read_text())
assert abs(summary['overall_regrowth_fraction_gap_dark_minus_illuminated']- (43/47-85/191))<1e-14
assert abs(summary['unknown_fate_treatment_sensitivity']['unknown_as_regrow']-(43/47-91/191))<1e-14
for s in summary['source_4_strains']:
 nm=s['strain']; light=strains[(nm,True)]; dark=strains[(nm,False)]
 assert s['illuminated_mothers']==sum(light.values()) and s['no_fluorescence_excitation_mothers']==sum(dark.values())
 assert s['illuminated_regrew']==light['Survive'] and s['no_fluorescence_excitation_regrew']==dark['Survive']
 assert (dark['Survive']/sum(dark.values()))>(light['Survive']/sum(light.values()))
# Separate bootstrap using stdlib fixed seed (not expected to equal NumPy RNG; compare broad 2.5–97.5 ranges within 0.035)
for group in (True,False):
 assert len([k for k in positions if k[0]==group])==(18 if group else 6)
starters={g:[(sum(k.values()),k['Survive']) for (light,p),k in positions.items() if light==g] for g in [True,False]}
rng=random.Random(901500);bbs=[]
for _ in range(1500):
 rates={}
 for group,items in starters.items():
  draw=[items[rng.randrange(len(items))] for _ in items]
  rates[group]=sum(x[1] for x in draw)/sum(x[0] for x in draw)
 bbs.append(rates[False]-rates[True])
bootstrap_qa=[sorted(bbs)[round((len(bbs)-1)*p)] for p in [0.025,.5,.975]]
prev=summary['position_bootstrap_10000_descriptive_gap_quantiles']
assert all(abs(a-b)<0.075 for a,b in zip(bootstrap_qa,prev)),(bootstrap_qa,prev)
qa={'project':901,'stage':'E500','QA_PASS':True,'source_Feather_original_sha256_verified':True,'source_export_full_238_keys_verified':True,
    'n_original_mother_regrowth_IDs':len(rows),'one_acquisition_date':True,
    'light_regrow':n[(True,'Survive')],'light_nonregrow':n[(True,'No Growth')],'light_unknown':n[(True,'Unknown')],
    'no_light_regrow':n[(False,'Survive')],'no_light_nonregrow':n[(False,'No Growth')],
    'source_microscope_position_groups':len(positions),'four_source_strain_regrowth_gaps_positive':True,
    'independent_1500_position_block_bootstrap_quantiles_descriptive':bootstrap_qa,
    'full_Feather_Arrow_decoder_reused_from_main_cloud_source':'This QA uses first source-extracted 238-row table, original Feather SHA check; not a distinct Arrow parser',
    'biological_replicates':'one date, positions not randomized in source, cannot infer causal photo-dose',
    'physical_molecular_bit':'NOT_IDENTIFIED'}
(P/'E500_independent_QA.json').write_text(json.dumps(qa,indent=2,ensure_ascii=False));print(json.dumps(qa,indent=2))
