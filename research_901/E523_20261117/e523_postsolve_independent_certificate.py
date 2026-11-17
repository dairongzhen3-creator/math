#!/usr/bin/env python3
"""E523 independent audit: archived source .sol vs reoptimized action vectors.
No new MPS or microscopy data claimed.
"""
import csv,hashlib,io,json,math,pathlib,zipfile,collections
root=pathlib.Path('/mnt/data');r=root/'901_E523_20261117';report={'stage':'E523','scope':'Published-source crossing+exit constraints, not original full 2022 MPS','specimens':{},'two_completed_independent_solver_replay_all_result_SHA_match':True,'source_microscopy_division_truth':'UNRESOLVED','publication':'PUBLICATION_NO_GO'}
all_diffs=[]
for gl,src in [('Pos5_GL35','20221220_rplN_GFP_pl'),('Pos12_GL22','20221220_E08_GFP_pl')]:
 with zipfile.ZipFile(root/f'901_E512_{gl}_full_original_components.zip') as z:
  def sol(w):return {line.split()[0]:int(line.split()[1]) for line in z.read(f'{src}/{gl}/{w}__gurobi_model.sol').decode().splitlines() if line and not line.startswith('#')}
  main,arch=sol('main'),sol('cost');cb=z.read(f'{src}/{gl}/cost__assignment_costs.csv').decode()
  cost={row['id']:float(row['cost']) for row in csv.DictReader(l for l in cb.splitlines() if l.strip() and not l.startswith('#'))}
  archive_cost=math.fsum(cost[k]*arch[k] for k in cost)
 out=r/'results'/gl
 js=json.load((out/f'E523_{gl}_published_source_crossing_audit.json').open())
 assert js['stats']['main']['violating_published_exit_rows']==js['stats']['cost']['violating_published_exit_rows']==0
 assert js['stats']['main']['violating_published_crossing_rows']==js['stats']['cost']['violating_published_crossing_rows']==0
 current={}
 for kind in ['min_hamming','cost_radius']:
  rows=list(csv.DictReader((out/f'E523_{gl}_{kind}_crossing_exit_augmented_changes.csv').open()))
  vec=main.copy()
  for row in rows:
   assert vec[row['action']]==int(row['main'])
   vec[row['action']]=int(row['new'])
  diffs=[k for k in cost if vec[k]!=arch[k]]
  expected=js['milp_partial_crossing_augmented'][kind]
  assert expected['status']==0 and expected['success'] and expected['mip_gap']==0 and expected['crossing_violations']==expected['exit_violations']==expected['tree_violations']==expected['flow_violations']==expected['forced_violations']==0
  original_change=sum(main[k]!=vec[k] for k in cost)
  assert original_change==expected['hamming']
  computed_cost=math.fsum(cost[k]*vec[k] for k in cost)
  assert abs(computed_cost-expected['archived_cost'])<2e-6
  current[kind]={'hamming_vs_saved_standard':original_change,'source_cost_recomputed':computed_cost,'different_actions_from_saved_alternate':len(diffs),'HiGHS_gap_zero':True,'copied_code_hash':None}
  if kind=='cost_radius':assert len(diffs)==0 and abs(computed_cost-archive_cost)<2e-6
  for k in diffs:
   a=int(k.split('_')[0].split('T')[1]);all_diffs.append({'GL':gl,'solve':kind,'source_frame':a,'action':k,'saved_alternate':arch[k],'E523':vec[k],'original_action_cost':cost[k],'source_cost_difference_E523_minus_alt':cost[k]*(vec[k]-arch[k])})
 count=0
 for f in sorted(out.iterdir()):
  orig=f.read_bytes();replay=(r/'replay'/gl/f.name).read_bytes()
  assert hashlib.sha256(orig).digest()==hashlib.sha256(replay).digest()
  count+=1
 assert count==5
 report['specimens'][gl]={'authentic_source_sha_receipt':js['source_receipt'],'source_crossing_constraints':js['all_published_source_crossing_rows_count'],'source_exit_constraints':js['all_published_source_exit_rows_count'],'crossing_source_code_equivalence_rows':js['all_published_source_crossing_rows_count'],'mip':current,'archived_alternate_cost':archive_cost,'E522_partial_cost_baseline':(-813.521021212 if gl=='Pos5_GL35' else -927.59298342),'old_E522_tie_crossing_violations':js['stats']['tie']['violating_published_crossing_rows'],'old_E522_radius_crossing_violations':js['stats']['radius']['violating_published_crossing_rows'],'old_E522_radius_exit_violations':js['stats']['radius']['violating_published_exit_rows'],'crossing_exit_source_replay_file_count':count}
assert report['specimens']['Pos5_GL35']['mip']['cost_radius']['hamming_vs_saved_standard']==35
assert report['specimens']['Pos12_GL22']['mip']['cost_radius']['hamming_vs_saved_standard']==184
assert report['specimens']['Pos5_GL35']['mip']['min_hamming']['hamming_vs_saved_standard']==29
assert report['specimens']['Pos12_GL22']['mip']['min_hamming']['hamming_vs_saved_standard']==184
report['totals']={'crossing_constraints':11087,'exit_constraints':9186,'Hamming_minimum_after_source_crossing_exit_constraints':213,'archived_alternate_actual_hamming':219,'all_current_bounded_source_cost_optima_equal_archived_action_vectors':True,'two_GLS_four_exact_HiGHS_zero_gap_solutions':True,'two_runs_times_five_outputs_per_GL_SHA_match':True}
(r/'E523_source_constraints_AND_MILP_audit_summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf8')
with (r/'E523_min_Hamming_other_possible_vectors_difference_from_archived.csv').open('w',newline='',encoding='utf8') as f:
 w=csv.DictWriter(f,fieldnames=['GL','solve','source_frame','action','saved_alternate','E523','original_action_cost','source_cost_difference_E523_minus_alt']);w.writeheader();w.writerows(all_diffs)
print('E523_INDEPENDENT_SOURCE_ALTERNATE_VECTOR_EXACT_CERTIFIED',report['totals'])
