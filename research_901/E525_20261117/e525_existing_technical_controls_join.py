#!/usr/bin/env python3
"""901 E525: Join E510 technical software-agreement controls to true-source
5-action timing-swap feasibility and original phase-only uncalibrated signal.
These are NOT biological ground truth. All 37 correlated windows are technical.
"""
import csv,zipfile,json,pathlib,argparse,hashlib,statistics,decimal
P=argparse.ArgumentParser();P.add_argument('--base',default='/mnt/data');P.add_argument('--e525',default='/mnt/data/901_E525_20261117');P.add_argument('--output',required=True);a=P.parse_args()
B=pathlib.Path(a.base);E=pathlib.Path(a.e525);O=pathlib.Path(a.output);O.mkdir(parents=True,exist_ok=True)
SHA={'Pos5_GL35':'745815a73bc1e2eaa01a2cf99f4edd9c530064d10c22bc5421f6ff5189629d8d','Pos12_GL22':'17f80af916f376dfa380b30e8a8b2e8459358a4b66e13e21e23f637aaa08cfa9'}
phase=list(csv.DictReader((E/'phase_run1/E525_37_true_phase_observability_vs_1frame_timing.csv').open()))
phase={(r['GL'],int(r['software_event_t'])):r for r in phase}
rows=[]
for gl in SHA:
 raw=(B/f'901_E510_{gl}_original_pixels.zip').read_bytes();assert hashlib.sha256(raw).hexdigest()==SHA[gl]
 with zipfile.ZipFile(B/f'901_E510_{gl}_original_pixels.zip') as z:
  assert z.testzip() is None
  key=json.loads(z.read([k for k in z.namelist() if k.endswith('E510_sealed_event_key_do_not_open_for_blind_review.json')][0]))
 registered=[int(r['frame']) for r in key if r['class']=='software_agreement']
 for item in csv.DictReader((E/f'final1_{gl}/E525_each_archived_solver_division_local_timing_availability.csv').open()):
  frame=int(item['frame_zero_based'])
  if frame not in registered:continue
  pr=phase[(gl,frame)]
  assert pr['technical_software_class']=='software_agreement'
  entry={'GL':gl,'frame_zero_based':frame,
       'source_solver_selected_div':item['software_solver_Div'],
       'source_software_agreement_technical_control':True,
       'number_local_source_feasible_5flip_shift_candidates':int(item['number_feasible_in_reconstructed_model']),
       'smallest_source_cost_delta_of_feasible_local_shift':item['min_feasible_positive_source_cost_delta'] if item['min_feasible_positive_source_cost_delta'] else 'NA',
       'direction_of_smallest_valid_local_shift':item['winning_shift_direction'] or 'NA',
       'phase_only_early_over_saved_frame_change':pr['early_over_saved_phase_change'],
       'phase_shift_not_verified_organism_division':True}
  rows.append(entry)
 assert len([x for x in rows if x['GL']==gl])==len(registered) and len(registered)==(4 if gl=='Pos5_GL35' else 9)
rows.sort(key=lambda x:(x['GL'],x['frame_zero_based']))
with (O/'E525_13_agreed_software_controls_timing_cost_phase_join.csv').open('w',newline='',encoding='utf8') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
vals=[decimal.Decimal(r['smallest_source_cost_delta_of_feasible_local_shift']) for r in rows if r['smallest_source_cost_delta_of_feasible_local_shift']!='NA']
summary={'total_technical_agreed_export_events':13,'n_with_at_least_one_valid_source_local_shift':sum(r['number_local_source_feasible_5flip_shift_candidates']>0 for r in rows),
 'n_with_local_best_margin_at_most_0_001':sum(v<=decimal.Decimal('0.001') for v in vals),
 'n_with_local_best_margin_at_most_0_01':sum(v<=decimal.Decimal('0.01') for v in vals),
 'quantities_are_not_cell_division_probabilities_or_independent_biological_replicates':True,
 'stated_phase_roi':'full_GL_x_0_105_y_80_254','phase_before_t_vs_t_to_t1_is_not_a_validated_septum_marker':True,
 'E524_E525_close_960_45_and_118_not_in_E510_phase_cache':True,'physical_truth':'UNRESOLVED'}
(O/'E525_13_agreed_software_controls_scope_summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
print('E525_TECHNICAL_CONTROLS_JOIN_PASS',summary,flush=True)
