#!/usr/bin/env python3
"""Independent source-Decimal cost audit for E525 vs E524 reported second-best.

Re-reads archived Gurobi action states and original exact cost strings, without
using E525 source parsed floating costs or MILP objective. It cannot by itself
establish complete MPS feasibility; that is checked in the E525 matrix audit.
"""
import csv,zipfile,hashlib,json,decimal,pathlib,argparse
P=argparse.ArgumentParser();P.add_argument('--base',default='/mnt/data');P.add_argument('--e525',default='/mnt/data/901_E525_20261117');P.add_argument('--output',required=True);a=P.parse_args()
BASE=pathlib.Path(a.base);IN=pathlib.Path(a.e525);O=pathlib.Path(a.output);O.mkdir(parents=True,exist_ok=True)
SHA={'Pos5_GL35':'77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6','Pos12_GL22':'73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2'}
SOURCE={'Pos5_GL35':'20221220_rplN_GFP_pl','Pos12_GL22':'20221220_E08_GFP_pl'}
decimal.getcontext().prec=50
result={}
for gl in SHA:
 data=(BASE/f'901_E512_{gl}_full_original_components.zip').read_bytes();assert hashlib.sha256(data).hexdigest()==SHA[gl]
 with zipfile.ZipFile(BASE/f'901_E512_{gl}_full_original_components.zip') as z:
  assert z.testzip() is None
  base=SOURCE[gl]+'/'+gl+'/'
  vals={p[0]:int(p[1]) for line in z.read(base+'cost__gurobi_model.sol').decode().splitlines() if line.strip() and not line.startswith('#') for p in [line.split()]}
  exact={r['id']:decimal.Decimal(r['cost']) for r in csv.DictReader(x for x in z.read(base+'cost__assignment_costs.csv').decode().splitlines() if x.strip() and not x.startswith('#'))}
  found=list(csv.DictReader((IN/('final1_'+gl)/'E525_all_feasible_local_timing_shift_5_action_proposals.csv').open()))
  times=[960] if gl=='Pos5_GL35' else [45,118]
  examples=[]
  for t in times:
   good=[r for r in found if int(r['archived_solver_Div_t'])==t and r['direction']=='early']
   assert len(good)==1,(gl,t,len(good))
   r=good[0];changes={}
   for text in r['five_flipped_binary_actions'].split(';'):
    name,value=text.rsplit('=',1);new=int(value)
    assert name in vals and name in exact and vals[name]!=new
    changes[name]=new
   assert len(changes)==5
   delta=sum((exact[k]*decimal.Decimal(changes[k]-vals[k]) for k in changes),decimal.Decimal(0))
   assert str(delta)==r['source_cost_change_decimal'],(delta,r['source_cost_change_decimal'])
   examples.append({'event_frame_zero_based':t,'new_division_frame_zero_based':t-1,'changed_actions':len(changes),'decimal_delta_exact_source':str(delta),'original_source_zip_sha256':SHA[gl]})
  result[gl]=examples
assert decimal.Decimal(result['Pos12_GL22'][1]['decimal_delta_exact_source'])<decimal.Decimal(result['Pos12_GL22'][0]['decimal_delta_exact_source'])
result['Pos12_E524_secondbest_not_global_source_exact_decimal']={'older_f45':result['Pos12_GL22'][0]['decimal_delta_exact_source'],'new_earlier_f118':result['Pos12_GL22'][1]['decimal_delta_exact_source'],
 'lower_by':str(decimal.Decimal(result['Pos12_GL22'][0]['decimal_delta_exact_source'])-decimal.Decimal(result['Pos12_GL22'][1]['decimal_delta_exact_source']))}
result['limits']='This checks exact original saved COST decimal and variable state toggles; feasibility proven separately only against reconstructed necessary model, not full original MPS.'
(O/'E525_independent_decimal_second_best_correction.json').write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf8')
print('E525_INDEPENDENT_DECIMAL_SOURCE_COST_CORRECTION_PASS',result['Pos12_E524_secondbest_not_global_source_exact_decimal'],flush=True)
