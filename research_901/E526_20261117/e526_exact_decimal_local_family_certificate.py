#!/usr/bin/env python3
"""E526 original-source exact integer-lattice certificate for E525 *local* 5-action moves.

This does not certify a global best-distinct solution of the reconstructed MILP,
and does not establish original MPS or biological division ground truth.
"""
from __future__ import annotations
import collections,csv,decimal,hashlib,io,json,pathlib,zipfile
from decimal import Decimal

BASE=pathlib.Path('/mnt/data')
OUT=BASE/'901_E526_20261117'/'results';OUT.mkdir(exist_ok=True,parents=True)
Z525=zipfile.ZipFile(BASE/'901_E525_20261117'/'901_E525_相邻帧分裂时间竞争与E524第二优值修正_全套复算_20261117.zip')
assert Z525.testzip() is None
FROZEN={
 'Pos5_GL35':('20221220_rplN_GFP_pl','77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6',60,91),
 'Pos12_GL22':('20221220_E08_GFP_pl','73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2',38,60)}
sha=lambda b:hashlib.sha256(b).hexdigest()
rows_all=[];summary={}; receipts=[]
for gl,(src,expect_sha,n_div,nfeasible) in FROZEN.items():
 path=BASE/f'901_E512_{gl}_full_original_components.zip';raw=path.read_bytes()
 assert sha(raw)==expect_sha,('official source ZIP SHA',gl)
 with zipfile.ZipFile(io.BytesIO(raw)) as z:
  assert z.testzip() is None
  prefix=f'{src}/{gl}/';rc=json.loads(z.read(prefix+'E512_original_member_SHA_CRC_critical_source_audit.json'))['complete_source_authentication']
  def orig(kind):
   content=z.read(prefix+'cost__'+kind)
   items=[x for x in rc if x['version']=='cost' and x['kind']==kind]
   assert len(items)==1 and sha(content)==items[0]['original_sha256'] and len(content)==items[0]['uncompressed_bytes'],kind
   return content
  csv_data=orig('assignment_costs.csv').decode()
  costs={r['id']:Decimal(r['cost']) for r in csv.DictReader(ln for ln in csv_data.splitlines() if ln.strip() and not ln.startswith('#'))}
  source_prices_txt={r['id']:r['cost'] for r in csv.DictReader(ln for ln in csv_data.splitlines() if ln.strip() and not ln.startswith('#'))}
  sol={ls.split()[0]:int(ls.split()[1]) for ls in orig('gurobi_model.sol').decode().splitlines() if ls.strip() and not ls.startswith('#')}
  assert len(costs)==(44231 if gl=='Pos5_GL35' else 47879)
  assert all(k in sol for k in costs)
  max_places=max(-v.as_tuple().exponent for v in costs.values())
  S=10**max_places
  assert max_places<=13
  intcost={k:int(v*S) for k,v in costs.items()}
  assert all(Decimal(n)==costs[k]*S for k,n in intcost.items())
  # Numeric audit of prior E525 float->str(float)->Decimal step.
  fmt_mismatch=[k for k,txt in source_prices_txt.items() if Decimal(str(float(txt)))!=Decimal(txt)]
  archived_total=sum(intcost[k] for k in costs if sol[k]==1)
  rows=list(csv.DictReader(io.StringIO(Z525.read(f'results/{gl}/E525_all_feasible_local_timing_shift_5_action_proposals.csv').decode())))
  assert len(rows)==nfeasible
  unique=set();n_changed=0
  for q in rows:
   changes={}
   for a in q['five_flipped_binary_actions'].split(';'):
    k,new=a.rsplit('=',1);assert k in costs and k not in changes and new in ('0','1')
    changes[k]=int(new)
   assert len(changes)==5 and all(sol[k]!=v for k,v in changes.items())
   fp=tuple(sorted(changes.items()));assert fp not in unique;unique.add(fp)
   d=sum((new-sol[k])*intcost[k] for k,new in changes.items())
   assert d>0
   report=Decimal(q['source_cost_change_decimal'])
   assert Decimal(d)/S == report,(gl,q['archived_solver_Div_t'],str(report),str(Decimal(d)/S))
   n_changed+=len(changes)
   rows_all.append({'GL':gl,'source_division_t':int(q['archived_solver_Div_t']),
       'direction':q['direction'],'delta_integer_lattice_units':d,'lattice_scale':S,
       'delta_exact_decimal':format(Decimal(d)/S,'f'),'five_flipped_binary_actions':q['five_flipped_binary_actions'],
       'validated_original_saved_solution':True,'source_S_BIAD2506_zip_sha256':expect_sha,
       'scope':'feasible E524 reconstructed constraints local 5-actions only; global next best not certified'})
  rr=sorted([q for q in rows_all if q['GL']==gl],key=lambda x:(x['delta_integer_lattice_units'],x['source_division_t'],x['direction']))
  best=rr[0];second=rr[1];thresholds=['0.00001','0.00002','0.0001','0.001','0.01']
  summary[gl]={
   'source_zip_sha256':expect_sha,'original_cost_member_sha256':sha(csv_data.encode()),
   'max_original_decimal_places':max_places,'exact_integer_lattice_scale':S,
   'source_costs_that_change_on_float_str_roundtrip':len(fmt_mismatch),
   'archived_alternate_total_exact_decimal':str(Decimal(archived_total)/S),
   'local_feasible_count':len(rr),'local_5bit_swaps_verified':n_changed,'minimum_gap':best['delta_exact_decimal'],
   'lowest_two_local_5flip_candidates':[{'t':q['source_division_t'],'direction':q['direction'],'delta':q['delta_exact_decimal']} for q in rr[:2]],
   'low_gap_counts_by_threshold':{t:sum(Decimal(q['delta_exact_decimal'])<=Decimal(t) for q in rr) for t in thresholds},
   'minimum_per_cost_coefficient_robustness_epsilon_for_saved_vs_local_5flips':str(Decimal(best['delta_exact_decimal'])/Decimal(5)),
   'does_not_certify_unrestricted_second_best':True,'model_limitation':'original full MPS not available',
   'physical_cell_truth':'UNRESOLVED'}
  receipts.append({'GL':gl,'zip_sha256':expect_sha,'zip_crc_pass':True,'original_member_sha_and_size_pass':True,'n_members':len(z.namelist())})
# Compare Pos12 118 to 45: disjoint support -> precise coefficient epsilon threshold Δdifference/10.
a={q['source_division_t']:q for q in rows_all if q['GL']=='Pos12_GL22' and q['direction']=='early' and q['source_division_t'] in (118,45)}
assert 118 in a and 45 in a
A={s.rsplit('=',1)[0] for s in a[118]['five_flipped_binary_actions'].split(';')};B={s.rsplit('=',1)[0] for s in a[45]['five_flipped_binary_actions'].split(';')}
assert len(A)==len(B)==5 and len(A|B)==10 and not A&B
D=Decimal(a[45]['delta_exact_decimal'])-Decimal(a[118]['delta_exact_decimal'])
assert D==Decimal('0.00000022')
summary['Pos12_GL22']['f118_vs_f45_original_cost_difference_exact_decimal']=str(D)
summary['Pos12_GL22']['f118_vs_f45_cost_order_per_coefficient_epsilon_10_disjoint_action_terms']=str(D/Decimal(10))
rows_all.sort(key=lambda q:(q['GL'],q['delta_integer_lattice_units'],q['source_division_t']))
with (OUT/'E526_151_local_moves_exact_integer_ranked.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=list(rows_all[0]));w.writeheader();w.writerows(rows_all)
with (OUT/'E526_source_cost_exact_lattice_certificate.json').open('w',encoding='utf-8') as f:json.dump({'exact_source_cost_integer_lattice':summary,'receipts':receipts,'input_family':'E525 151 feasibility-verified 5-bit swaps, independently replayed original cost prices','does_not_certify_global_second_best':True,'biological_truth':'UNRESOLVED'},f,indent=2,ensure_ascii=False,sort_keys=True)
print('E526_EXACT_LATTICE_CERTIFIED',json.dumps(summary,ensure_ascii=False,sort_keys=True))
