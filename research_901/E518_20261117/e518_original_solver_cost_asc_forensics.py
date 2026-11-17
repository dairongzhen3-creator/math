#!/usr/bin/env python3
"""901 E518: reproducible archived MoMA optimization/cost/curation forensic audit.

Authoritative inputs: E512 original S-BIAD2506 ZIP64 CRC/SHA-certified source
artifacts + E512b fixed exact dispute mother-branch registry. This is not a
new Gurobi reoptimization. A full solution to one MoMA version may violate
constraints of the other, so cross-version cost gaps are descriptive only.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,pathlib,re,statistics,zipfile

SRC={
 'Pos5_GL35': {'e512_sha':'77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6',
   'e512b_sha':'efe81c19a76117c0f9bebe58b8e4525f93c407465312d0402ce7096e9866b384',
   'events':[247,294,399,403,405,409,411,415,622], 'active_components':5500,
   'variables':49731, 'archive':'20221220_rplN_GFP_pl.zip'},
 'Pos12_GL22': {'e512_sha':'73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2',
   'e512b_sha':'4893c521768bdbcc026e7f6d997b1266f9e70e6869b9bf08e4b429a8aa7ecca9',
   'events':[241,246,608], 'active_components':5765,
   'variables':53644, 'archive':'20221220_E08_GFP_pl.zip'}
}
def sha(x):return hashlib.sha256(x).hexdigest()
def get_single(z,suffix):
 names=[x for x in z.namelist() if x.endswith(suffix)];assert len(names)==1,(suffix,names)
 return z.read(names[0])

def parse_solver(b):
 headers=[];variables={}
 for line in b.decode('utf8').splitlines():
  if not line.strip():continue
  if line.startswith('#'):headers.append(line);continue
  name,value=line.rsplit(None,1)
  assert name not in variables,name
  variables[name]=float(value)
 assert all(v in (0.,1.) for v in variables.values())
 return headers,variables

def parse_asc(data):
 records=[];kinds=collections.Counter()
 for line in data.decode('utf8').splitlines():
  line=line.strip()
  if not line or line.startswith('#'):continue
  fields=[v.strip() for v in line.split(',')]
  kinds[fields[0]]+=1
  if fields[0]=='ASC':
   assert len(fields)==5,fields
   records.append({'time_printed':int(fields[1]),'internal_constraint_id':int(fields[2]),
       'rhs':float(fields[3]),'variable':fields[4]})
 return records,dict(kinds)

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--input',default='/mnt/data')
 p.add_argument('--out',default='/mnt/data/901_E518_20261117/results')
 args=p.parse_args(); inp=pathlib.Path(args.input);out=pathlib.Path(args.out);out.mkdir(parents=True,exist_ok=True)
 allsummary=[];detailed=[];constraint_rows=[];contrib_rows=[];receipts=[]
 for gl,q in SRC.items():
  paths={'E512':inp/f'901_E512_{gl}_full_original_components.zip',
    'E512b':inp/f'901_E512b_{gl}_exact_div_to_map.zip'}
  raw={k:f.read_bytes() for k,f in paths.items()}
  assert sha(raw['E512'])==q['e512_sha']
  assert sha(raw['E512b'])==q['e512b_sha']
  z={k:zipfile.ZipFile(io.BytesIO(v)) for k,v in raw.items()}
  assert all(v.testzip() is None for v in z.values())
  receipts.extend([{'gl':gl,'phase':key,'name':paths[key].name,'SHA256':sha(raw[key]),
     'bytes':len(raw[key]),'ZIP_CRC_ALL_PASS':True} for key in paths])
  receipt=json.loads(get_single(z['E512'],'E512_original_member_SHA_CRC_critical_source_audit.json'))
  assert receipt['gl']==gl
  target={ (x['version'],x['kind']):x for x in receipt['complete_source_authentication'] }
  def read(version,kind):
   b=get_single(z['E512'],version+'__'+kind)
   assert sha(b)==target[(version,kind)]['original_sha256'],(gl,version,kind)
   assert len(b)==target[(version,kind)]['uncompressed_bytes']
   return b
  orig_sol={};orig_mst={};headers={};asc={};kinds={}
  for ver in ['main','cost']:
   headers[ver],orig_sol[ver]=parse_solver(read(ver,'gurobi_model.sol'))
   _,orig_mst[ver]=parse_solver(read(ver,'gurobi_model.mst'))
   asc[ver],kinds[ver]=parse_asc(read(ver,'curation.moma'))
   assert orig_sol[ver].keys()==orig_mst[ver].keys()
   assert all(orig_sol[ver][v]==orig_mst[ver][v] for v in orig_sol[ver])
   assert len(orig_sol[ver])==q['variables']
   assert sum(1 for v in orig_sol[ver].values() if v>0.5)==q['active_components']
   for record in asc[ver]:
    assert record['variable'] in orig_sol[ver]
    assert orig_sol[ver][record['variable']]==record['rhs']
  assert orig_sol['main'].keys()==orig_sol['cost'].keys()
  csv_lines=[line for line in read('cost','assignment_costs.csv').decode('utf8').splitlines()
     if line.strip() and not line.startswith('#')]
  data=list(csv.DictReader(csv_lines))
  assert set(data[0])=={'id','type','source_frame','cost'}
  costs={};types={};frame_by_var={}
  for r in data:
   n=r['id'];assert n not in costs
   costs[n]=float(r['cost']);types[n]=r['type'];frame_by_var[n]=int(r['source_frame'])
  assert set(costs).issubset(orig_sol['main'])
  # Gurobi objectives are the sum of the printed assignment costs, plus zero-cost
  # auxiliary variables, as evidenced by both saved objective header values.
  derived={v:math.fsum(costs[k]*orig_sol[v][k] for k in costs) for v in ('main','cost')}
  published={v:float(next(l.split('=')[1] for l in headers[v] if 'Objective value' in l)) for v in ('main','cost')}
  assert all(abs(derived[v]-published[v])<2e-5 for v in ('main','cost')),(gl,derived,published)
  diffs=[name for name in orig_sol['main'] if orig_sol['main'][name]!=orig_sol['cost'][name]]
  kinds_changed=dict(collections.Counter(re.split(r'T\d+_',k)[0] for k in diffs))
  def violation(vers):
   reference='cost' if vers=='main' else 'main'
   rows=[]
   for r in asc[reference]:
    actual=orig_sol[vers][r['variable']]
    rows.append({'gl':gl,'constraint_source':reference,'solution_source':vers,
       'variable':r['variable'],'record_time_printed':r['time_printed'],
       'rhs':r['rhs'],'actual_solution_value':actual,
       'VIOLATED':bool(actual!=r['rhs'])})
   return rows
  cross_m=violation('main');cross_c=violation('cost')
  constraint_rows.extend(cross_m+cross_c)
  registry=json.loads(get_single(z['E512b'],'E512b_exact_software_division_constraint_mechanism.json'))
  assert registry['gl']==gl
  events=registry['disputed_events'];assert sorted(e['source_event_frame_0_based'] for e in events)==q['events']
  direct_n=0;direct_violated=0
  for e in events:
   f=int(e['source_event_frame_0_based']); div=e['main_selected_div_variable']; mapname=e['alternate_selected_mapping']
   assert orig_sol['main'][div]==1 and orig_sol['cost'][div]==0
   assert orig_sol['main'][mapname]==0 and orig_sol['cost'][mapname]==1
   direct=[c for c in asc['cost'] if c['variable']==mapname and c['rhs']==1]
   assert bool(direct)==e['alternate_curation_has_direct_forced_mapping']
   direct_n+=bool(direct);direct_violated+=bool(direct and orig_sol['main'][mapname]==0)
   # Fixed source time group: assign the same archival per-variable cost terms
   # without claiming these are stand-alone physical or causal event effects.
   changes=[(name,(orig_sol['cost'][name]-orig_sol['main'][name])*costs[name])
            for name in costs if frame_by_var[name]==f and orig_sol['main'][name]!=orig_sol['cost'][name]]
   delta=math.fsum(d for _,d in changes)
   row={'gl':gl,'event_f_zero_based':f,'source_div_var':div,'source_map_var':mapname,
        'main_div_selected':1,'cost_div_selected':0,'main_map_selected':0,'cost_map_selected':1,
        'direct_ASC_rhs_1':bool(direct),
        'direct_ASC_violated_by_main_solution':bool(direct and orig_sol['main'][mapname]==0),
        'source_time_cost_delta_alternate_minus_main':delta,
        'source_time_cost_terms_changed':len(changes),
        'source_time_different_selected_map_count':sum(x.startswith('MapT') for x,_ in changes),
        'source_time_different_selected_div_count':sum(x.startswith('DivT') for x,_ in changes),
        'source_time_different_selected_exit_count':sum(x.startswith('ExitT') for x,_ in changes),
        'source_y120_daughter':e.get('any_source_daughter_minY_120',None)}
   detailed.append(row)
  assert direct_n==registry['direct_forced_alternate_mapping_count']
  assert direct_violated==direct_n
  class_delta={k:math.fsum((orig_sol['cost'][name]-orig_sol['main'][name])*price
       for name,price in costs.items() if types[name]==k) for k in sorted(set(types.values()))}
  assert abs(math.fsum(class_delta.values())-(derived['cost']-derived['main']))<1e-8
  allsummary.append({'gl':gl,'official_archive':q['archive'],'source_e512_SHA256':q['e512_sha'],
    'source_e512b_SHA256':q['e512b_sha'],'common_assignment_cost_entries':len(costs),
    'source_Gurobi_variables_each':q['variables'],
    'source_MST_selected_values_identical_to_SOL_all_variables_in_both_versions':True,
    'MST_NOT_PROVEN_TO_BE_PRE_CURATION':True,
    'original_Gurobi_objective':published,'reconstructed_both_objectives_from_cost_version_assignment_ledger':derived,
    'objective_reconstruction_absolute_error':{v:abs(derived[v]-published[v]) for v in ('main','cost')},
    'fixed_common_ledger_objective_delta_cost_minus_main':derived['cost']-derived['main'],
    'by_assignment_type_delta':class_delta,
    'source_curation_counts':kinds,
    'saved_cost_ASC_total':len(asc['cost']),'saved_main_ASC_total':len(asc['main']),
    'cost_ASC_violated_by_original_main_solution':sum(v['VIOLATED'] for v in cross_m),
    'main_ASC_violated_by_original_cost_solution':sum(v['VIOLATED'] for v in cross_c),
    'source_disputed_events':len(events),'events_with_exact_forced_continuation_ASC':direct_n,
    'forced_continuation_ASCs_that_original_main_does_not_satisfy':direct_violated,
    'different_binary_solution_variables':len(diffs),
    'different_solution_variables_by_type':kinds_changed,
    'cross_model_feasibility_proven':False,
    'causal_effect_of_manual_curation_proven':False,
    'source_division_ground_truth_proven':False,
    'physical_information_bits_identified':False,
    'publication_state':'PUBLICATION_NO_GO'})
  print('E518_OFFICIAL_SAVED_MODEL_AUDIT_PASS',gl,json.dumps({
   'objective_delta':derived['cost']-derived['main'],
   'reconstructed_main':derived['main'],
   'reconstructed_cost':derived['cost'],
   'ASC_violated_by_main':sum(v['VIOLATED'] for v in cross_m),
   'direct_disputed_ASC':direct_n,
   'MST_same_SOL':True,
   'changed_var':len(diffs),
   'delta_by_type':class_delta},ensure_ascii=False),flush=True)
 summary={'study':'901 E518','source_dataset':'BioImage Archive S-BIAD2506 original archived author MoMA source',
  'accounting':'same CSV assignment cost vector evaluated on two stored binary solutions; not an ablation, not proof one is feasible under the other version',
  'original_source_both_ZIP_SHA_CRC_pass':True,
  'total_direct_DISPUTED_ASC':sum(x['events_with_exact_forced_continuation_ASC'] for x in allsummary),
  'total_curation_ASC_violations_of_prior_standard_solutions':sum(x['cost_ASC_violated_by_original_main_solution'] for x in allsummary),
  'gls':allsummary,'source_receipts':receipts,'note':'Gurobi .mst has the same selected variables as each saved .sol; it CANNOT be used as an independent pre-curation snapshot.',
  'limits':['Full constraint feasibility across two MoMA versions is NOT proven','No Gurobi reoptimization','No independent microscope truth','No physical molecular copy or information bit/J inferred'],
  'publication':'PUBLICATION_NO_GO'}
 (out/'E518_archived_common_cost_objective_and_ASC_audit.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
 def write_csv(filename,rows):
  with (out/filename).open('w',newline='',encoding='utf8') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 write_csv('E518_all_disputed_events_exact_forced_Map_and_cost_by_frame.csv',detailed)
 write_csv('E518_all_curation_ASC_cross_solution_constraint_compatibility.csv',constraint_rows)
 write_csv('E518_original_ZIP_SHA_CRC_source_receipts.csv',receipts)
 (out/'E518_conclusions_and_limitations.md').write_text(
    '# E518 原求解目标与存档强制约束独立核验\n\n'
    +f'总共直接强制Map的争议事件: {summary["total_direct_DISPUTED_ASC"]}/12；'
    +f'备用版ASC与原标准版解直接矛盾的条数: {summary["total_curation_ASC_violations_of_prior_standard_solutions"]}。\n\n'
    +'两份模型在相同成本账本下评分，备用解目标更高。'
    +'这是已保存解之间的可比线性函数取值，并不证明旧解满足新版完整约束或因果效果。'
    +'MST与同版本SOL变量状态逐一相等，不能用MST假装“修订前”独立解。'
    +'原物理分裂和分子bits仍未被识别。\n',encoding='utf8')
 print('E518_FINAL_SOURCE_REPRO_PASS',json.dumps({'gls':len(allsummary),'events':len(detailed),
 'saved_ASC_cross_compatibility_rows':len(constraint_rows),'source_artifacts':len(receipts)},ensure_ascii=False),flush=True)
if __name__=='__main__':main()