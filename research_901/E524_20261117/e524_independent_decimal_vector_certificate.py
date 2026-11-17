#!/usr/bin/env python3
"""901 E524 independent Decimal source-vector audit, no solver and no biological labels.

Original E512 ZIP source files were previously SHA/CRC and per-member authenticated.
Re-read the exact original Gurobi solution and cost CSV; validate original software
cost, reproduced uncapped optimal action vector, no-good runner-up, and terminal
zero-cost Exit discrepancy caused by an omitted t=1035 continuity boundary.
Not a full archived MPS feasibility check and not a cell-biology truth label.
"""
from pathlib import Path
import hashlib,zipfile,csv,io,json,re
from decimal import Decimal
BASE=Path('/mnt/data')
OUT=BASE/'901_E524_20261117'/'results';OUT.mkdir(parents=True,exist_ok=True)
SRC={'Pos5_GL35':'20221220_rplN_GFP_pl','Pos12_GL22':'20221220_E08_GFP_pl'}
SHA={'Pos5_GL35':'77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6','Pos12_GL22':'73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2'}
EXPECT={'Pos5_GL35':{'N':44231,'hamming':35,'cost':Decimal('-813.431953329'),'next_gap':Decimal('0.00000519'),'second_frame':[959,960],'missed_exit':4},'Pos12_GL22':{'N':47879,'hamming':184,'cost':Decimal('-923.86545359'),'next_gap':Decimal('0.00001587'),'second_frame':[44,45],'missed_exit':5}}
def load_sol(s):
 return {line.split()[0]:int(float(line.split()[1])) for line in s.splitlines() if line.strip() and not line.startswith('#')}
def read_rows(p):
 with p.open(newline='',encoding='utf8') as f:return list(csv.DictReader(f))
def sha(b):return hashlib.sha256(b).hexdigest()
findings={};rows=[]
for gl,stem in SRC.items():
 e=EXPECT[gl];zz=BASE/f'901_E512_{gl}_full_original_components.zip';raw=zz.read_bytes();assert sha(raw)==SHA[gl]
 with zipfile.ZipFile(io.BytesIO(raw)) as z:
  assert z.testzip() is None
  root=f'{stem}/{gl}/';rc=json.loads(z.read(root+'E512_original_member_SHA_CRC_critical_source_audit.json'))['complete_source_authentication']
  def member(v,k):
   b=z.read(root+v+'__'+k);a=[r for r in rc if r['version']==v and r['kind']==k]
   assert len(a)==1 and sha(b)==a[0]['original_sha256'] and len(b)==a[0]['uncompressed_bytes']
   return b
  main=load_sol(member('main','gurobi_model.sol').decode());alt=load_sol(member('cost','gurobi_model.sol').decode())
  lines=member('cost','assignment_costs.csv').decode().splitlines()
  cost={r['id']:Decimal(r['cost']) for r in csv.DictReader(line for line in lines if line.strip() and not line.startswith('#'))}
  names=[k for k in main if k in cost]
  assert len(names)==e['N'] and set(names)==set(cost)
  assert set(main)==set(alt)
  def total(v):return sum((cost[k]*v[k] for k in names),Decimal(0))
  archived_cost=total(alt);assert archived_cost==e['cost'],(gl,archived_cost)
  def reconstruct(p):
   v={k:main[k] for k in names};rr=read_rows(p)
   for d in rr:
    k=d['action'];assert k in v and v[k]==int(d['main']);v[k]=int(d['new'])
    assert int(d['saved_alt'])==alt[k]
   return v,rr
  folder='pos5' if gl.startswith('Pos5') else 'pos12'
  primary,change=reconstruct(BASE/'901_E524_20261117'/('terminal_'+folder)/f'E524_{gl}_cost_uncapped_source_crossing_exit_changes.csv')
  assert len(change)==e['hamming'] and all(primary[k]==alt[k] for k in names)
  runner,change2=reconstruct(BASE/'901_E524_20261117'/('nogood_'+folder)/f'E524_{gl}_cost_best_distinct_from_saved_alt_source_crossing_exit_changes.csv')
  diff=[k for k in names if runner[k]!=alt[k]];assert len(diff)==5
  runner_cost=total(runner);gap=runner_cost-archived_cost;assert gap==e['next_gap'],(gl,gap)
  ticks=sorted(set(int(re.match(r'^(?:Map|Div|Exit|Lys)T(\d+)_',k).group(1)) for k in diff));assert ticks==e['second_frame']
  # Earlier run with missing terminal continuity gave a false zero-cost degeneracy.
  uncorrected,_=reconstruct(BASE/'901_E524_20261117'/('probe_'+folder)/f'E524_{gl}_cost_uncapped_source_crossing_exit_changes.csv')
  omissions=[k for k in names if uncorrected[k]!=alt[k]]
  assert len(omissions)==e['missed_exit'] and all(k.startswith('ExitT1035_') and cost[k]==0 and alt[k]==1 and uncorrected[k]==0 for k in omissions)
  assert total(uncorrected)==archived_cost
  # Both E523 bounded source vectors maintained original terminal Exit, so the
  # E524 correction does not automatically invalidate their successful solutions.
  with zipfile.ZipFile(BASE/'901_E523_20261117'/'901_E523_交叉与退出约束_双生长线优化闭合与复算全套_20261117.zip') as zz523:
   for mode in ['min_hamming','cost_radius']:
    name=f'results/{gl}/E523_{gl}_{mode}_crossing_exit_augmented_changes.csv'
    q=list(csv.DictReader(io.StringIO(zz523.read(name).decode())))
    assert all(not r['action'].startswith('ExitT1035_') for r in q)
  # Numerical HiGHS results: optimization parameters and certified bounds.
  for group,mode,prefix in [('terminal','cost_uncapped','uncapped'),('nogood','cost_best_distinct_from_saved_alt','nogood')]:
   p=BASE/'901_E524_20261117'/(f'{group}_{folder}')/f'E524_{gl}_source_crossing_exit_{prefix}_audit.json'
   d=json.loads(p.read_text())['milp_partial_source_crossing_exit_'+prefix][mode]
   assert d['success'] and d['status']==0 and d['mip_gap']==0
   assert all(d[x]==0 for x in ['crossing_violations','exit_violations','tree_violations','flow_violations','forced_violations'])
   v=primary if group=='terminal' else runner
   assert abs(Decimal(str(d['archived_cost']))-total(v))<Decimal('0.000000001')
  findings[gl]={'sha256_E512_zip':sha(raw),'zip_CRC_all_pass':True,'all_used_original_member_SHA_and_size_pass':True,'actions':len(names),'source_archived_cost_decimal':str(archived_cost),'uncapped_solved_matches_archived_actions':True,'uncapped_hamming_from_standard':e['hamming'],'incorrect_terminal_exits_in_previous_uncapped_probe':len(omissions),'terminal_exits_repaired':True,'best_distinct_source_cost_decimal':str(runner_cost),'best_distinct_increment_decimal':str(gap),'best_distinct_differs_from_archived_actions':len(diff),'best_distinct_frames':ticks,'full_original_MPS_available':False,'biological_division_truth':'UNRESOLVED'}
  for k in diff:
   rows.append({'GL':gl,'variable':k,'event_frame':int(re.match(r'^(?:Map|Div|Exit|Lys)T(\d+)_',k).group(1)), 'archived_alt':alt[k],'best_distinct':runner[k],'source_cost_exact_decimal':str(cost[k])})
(OUT/'E524_independent_source_decimal_vector_certificate.json').write_text(json.dumps(findings,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf8')
with (OUT/'E524_best_distinct_10_action_changes.csv').open('w',encoding='utf8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('E524_INDEPENDENT_DECIMAL_SOURCE_VECTOR_PASS',json.dumps({k:{'cost':v['source_archived_cost_decimal'],'gap':v['best_distinct_increment_decimal'],'action_count':v['actions'],'second_best_changes':v['best_distinct_differs_from_archived_actions']} for k,v in findings.items()},ensure_ascii=False),flush=True)
