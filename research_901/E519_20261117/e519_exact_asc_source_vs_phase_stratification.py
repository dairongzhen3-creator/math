#!/usr/bin/env python3
"""901/E519: reproducible source decision ↔ phase-only technical controls audit.

All software event classifications come solely from the authenticated archived
MoMA analysis and never from independent microscopic biological adjudication.
The phase-only gradient was computed before this cross-layer join in E517; its
ROI was post-hoc after y=120 was identified. No significance tests, biological
false division decisions, absolute protein bits or energy estimates.

Input: E518 / E515 / E517 frozen CSVs; E512 original JSON components and
E518 archived ASC cross-solution list. Source SHA and ZIP CRC are checked.
Outputs are segregated from the blind reviewer packet and do not contain
any private blind case identifiers/keys.
"""
import argparse, hashlib, json, pathlib, zipfile, csv
from collections import Counter
import pandas as pd
import numpy as np

P=argparse.ArgumentParser()
P.add_argument('--base',default='/mnt/data')
P.add_argument('--output',default='/mnt/data/901_E519_20261117/results')
a=P.parse_args();B=pathlib.Path(a.base);OUT=pathlib.Path(a.output);OUT.mkdir(parents=True,exist_ok=True)

files={
 'E518_12_events':B/'901_E518_20261117/results/E518_all_disputed_events_exact_forced_Map_and_cost_by_frame.csv',
 'E518_ASC_cross_solution':B/'901_E518_20261117/results/E518_all_curation_ASC_cross_solution_constraint_compatibility.csv',
 'E515_source_candidates':B/'901_E515_20261117/results/E515_source_event_probabilities_and_preexistence.csv',
 'E517_phase_metric_all_frames':B/'901_E517_20261117/E517_all_203_phc_threshold_free_y120_crossing_measurements.csv',
 'E512_Pos5_source':B/'901_E512_Pos5_GL35_full_original_components.zip',
 'E512_Pos12_source':B/'901_E512_Pos12_GL22_full_original_components.zip'
}
sha_expect={
 'E518_12_events':'0d17b98280a98258147315d7ae3cfe537f7f646c18a84ba72676d01f0a37614f',
 'E518_ASC_cross_solution':'1acb518570f5561b931d03acec65688eb07c0bc6a12764eace470fc147bf5c56',
 'E515_source_candidates':'207020cf2f650451584805ae27b11ed1640aa982e1db84a2a9ee194474280f62',
 'E517_phase_metric_all_frames':'6891433bac9e78a51e9aa59411d286772d851b2fff3e8c7665c471381561a679',
 'E512_Pos5_source':'77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6',
 'E512_Pos12_source':'73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2'
}
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
receipts={}
for name,p in files.items():
 assert p.is_file(),(name,p)
 h=sha(p)
 if name in sha_expect:assert h==sha_expect[name],(name,'SHA_CHANGED',h)
 receipts[name]={'path_relative_to_base':str(p.relative_to(B)),'bytes':p.stat().st_size,'sha256':h}
 if p.suffix=='.zip':
  with zipfile.ZipFile(p) as z:
   assert z.testzip() is None,(name,'CRC_FAILURE')
   receipts[name]['ZIP_CRC_PASS']=True
   receipts[name]['members']=len(z.namelist())

D=pd.read_csv(files['E518_12_events']);A=pd.read_csv(files['E518_ASC_cross_solution'])
C=pd.read_csv(files['E515_source_candidates']);I=pd.read_csv(files['E517_phase_metric_all_frames'])
assert len(D)==12 and len(C)==25 and len(I)==203
assert set(D.gl)==set(C.gl)==set(I.gl)=={'Pos5_GL35','Pos12_GL22'}
assert set(D.event_f_zero_based)=={247,294,399,403,405,409,411,415,622,241,246,608}
assert C[['gl','frame']].duplicated().sum()==0
assert I[['gl','frame_t_zero_based']].duplicated().sum()==0
assert C.software_label.value_counts().to_dict()=={'software_agreement':13,'disputed':12}
D['true_image_target_frame']=D.event_f_zero_based+1
C['true_image_target_frame']=C.frame+1
joined=D.merge(C.drop(columns=['true_image_target_frame']),left_on=['gl','event_f_zero_based'],right_on=['gl','frame'],validate='one_to_one',how='left')
assert len(joined)==12 and (joined.software_label=='disputed').all()
joined=joined.merge(I,left_on=['gl','true_image_target_frame'],right_on=['gl','frame_t_zero_based'],validate='many_to_one',how='left')
assert len(joined)==12 and joined['y119_to_y120_over_local_median_x42_66'].notna().all()

classify=[]
for row in joined.itertuples():
 c=str(row.gl);f=int(row.event_f_zero_based)
 direct=bool(row.direct_ASC_rhs_1)
 kind='same_mother_ASC_required_map' if direct else 'NO_DIRECT_SAME_MOTHER_ASC'
 small=int(row.small_candidate_area_px);g=float(row.y119_to_y120_over_local_median_x42_66)
 rec={'gl':c,'f_0_based':f,'f_plus_1_image_0_based':f+1,'recorded_mechanism_group':kind,
  'direct_same_mother_forced_map':direct,
  'small_software_candidate_pixel_area':small,
  'small_software_candidate_ymin':int(row.small_candidate_ymin),
  'small_software_candidate_unnormalized_UNet_mean':float(row.small_MoMA_probtiff_mean_uncalibrated),
  'software_predecessor_max_IoU_five_frames':float(row.best_past_5_frames_source_component_iou),
  'independent_phase_y119_120_edge_ratio_x42_65':g,
  'independent_phase_y119_120_edge_ratio_x30_74':float(row.y119_to_y120_over_local_median_x30_75),
  'alternate_minus_main_local_assignment_cost':float(row.source_time_cost_delta_alternate_minus_main),
  'n_changed_costed_assignments_same_frame':int(row.source_time_cost_terms_changed),
  'optical_signal_is_NOT_biological_truth':True,
  'event_is_biologically_verified':False}
 classify.append(rec)
classify=sorted(classify,key=lambda x:(x['gl'],x['f_0_based']))
assert sum(x['direct_same_mother_forced_map'] for x in classify)==10
pos=[x for x in classify if x['direct_same_mother_forced_map']]
neg=[x for x in classify if not x['direct_same_mother_forced_map']]
assert sorted(x['f_0_based'] for x in neg)==[247,294]
assert all(x['small_software_candidate_pixel_area']<=61 for x in pos)
assert sorted(x['small_software_candidate_pixel_area'] for x in neg)==[103,195]
assert all(x['small_software_candidate_ymin']==120 for x in pos)
assert all(x['independent_phase_y119_120_edge_ratio_x42_65']<1 for x in pos)
assert all(x['independent_phase_y119_120_edge_ratio_x42_65']>1 for x in neg)

# Preserve software-agreement events as *technical* controls, not biology ground truth.
controls=C[C.software_label=='software_agreement'].merge(I,left_on=['gl','true_image_target_frame'],right_on=['gl','frame_t_zero_based'],how='left',validate='one_to_one')
assert len(controls)==13 and controls.frame_t_zero_based.notna().all()
ccols=['gl','frame','small_candidate_area_px','small_candidate_ymin','best_past_5_frames_source_component_iou',
       'small_MoMA_probtiff_mean_uncalibrated','y119_to_y120_over_local_median_x42_66']
controls=controls[ccols].sort_values(['gl','frame'])

# Compare the one same-time conflicting constraint of P5 f247 to the actual
# dispute mother's original source pixels, component-tree ancestry and root.
with zipfile.ZipFile(files['E512_Pos5_source']) as z:
 path='20221220_rplN_GFP_pl/Pos5_GL35/'
 receipt=json.loads(z.read(path+'E512_original_member_SHA_CRC_critical_source_audit.json'))
 b=z.read(path+'main__component_forests.json')
 expected=next(q['original_sha256'] for q in receipt['complete_source_authentication'] if q['version']=='main' and q['kind']=='component_forests.json')
 assert hashlib.sha256(b).hexdigest()==expected
 comps={q['stringId']:q for q in json.loads(b)['components']}
 assert len(comps)==5500
p247=D[(D.gl=='Pos5_GL35')&(D.event_f_zero_based==247)].iloc[0]
asc=A[(A.gl=='Pos5_GL35')&(A.constraint_source=='cost')&(A.solution_source=='main')&(A.record_time_printed==248)]
assert len(asc)==1 and bool(asc.iloc[0].VIOLATED)
forced=str(asc.iloc[0].variable)
main_mom=str(p247.source_div_var).split('_')[1];forced_mom=forced.split('_')[1]
assert main_mom in comps and forced_mom in comps and main_mom!=forced_mom

def root_chain(k):
 a=[];vis=set()
 while k not in (None, "NA", "null", ""):
  assert k not in vis and k in comps,(k,vis)
  vis.add(k);a.append(k)
  k=comps[k]['parentStringId']
 return a
chain1=root_chain(main_mom);chain2=root_chain(forced_mom)
geom1=set(zip(comps[main_mom]['xCoordinates'],comps[main_mom]['yCoordinates']))
geom2=set(zip(comps[forced_mom]['xCoordinates'],comps[forced_mom]['yCoordinates']))
y1=[y for x,y in geom1];y2=[y for x,y in geom2]
spatial={'GL':'Pos5_GL35','software_event_f':247,'forced_ASC_record_time_printed':248,
 'disputed_source_div_mother_component_id':main_mom,
 'other_ASC_forced_map_mother_component_id':forced_mom,
 'disputed_mother_y_extent':[min(y1),max(y1)],
 'ASC_map_other_mother_y_extent':[min(y2),max(y2)],
 'disputed_mother_component_pixel_count':len(geom1),
 'other_forced_mother_component_pixel_count':len(geom2),
 'pixel_intersection_count':len(geom1&geom2),
 'same_component_tree_ancestor':len(set(chain1)&set(chain2))>0,
 'disputed_component_tree_root_id':chain1[-1],
 'forced_other_component_tree_root_id':chain2[-1],
 'vertical_empty_gap_in_rows':min(y1)-max(y2)-1,
 'same_frame_ASC_does_not_directly_target_disputed_mother':True,
 'not_proof_absence_of_global_optimization_coupling':True}
assert spatial['pixel_intersection_count']==0
assert not spatial['same_component_tree_ancestor']
assert spatial['vertical_empty_gap_in_rows']==23

# E518 ASC record time is original assignment source_time+1; do not confuse.
for d in classify:
 gl=d['gl'];f=d['f_0_based'];t=f+1
 aa=A[(A.gl==gl)&(A.constraint_source=='cost')&(A.solution_source=='main')]
 assert len(aa)==(25 if gl=='Pos5_GL35' else 216)
 v=np.abs(aa.record_time_printed.astype(int).values-t)
 d['closest_any_ASC_record_time_distance_frames']=int(v.min())
 d['same_frame_other_ASC_count']=int(sum(aa.record_time_printed==t))
 d['same_frame_violated_ASC_count']=int(sum((aa.record_time_printed==t)&aa.VIOLATED))
 d['number_any_cost_ASC_within_5_frames']=int(sum(v<=5))
 if d['direct_same_mother_forced_map']:
  assert d['closest_any_ASC_record_time_distance_frames']==0
f294=next(x for x in classify if x['gl']=='Pos5_GL35' and x['f_0_based']==294)
assert f294['closest_any_ASC_record_time_distance_frames']==47
assert f294['number_any_cost_ASC_within_5_frames']==0
assert f294['alternate_minus_main_local_assignment_cost']<0
f247=next(x for x in classify if x['gl']=='Pos5_GL35' and x['f_0_based']==247)
assert f247['closest_any_ASC_record_time_distance_frames']==0
assert f247['same_frame_violated_ASC_count']==1
assert f247['alternate_minus_main_local_assignment_cost']>0

# A simple optical threshold applied post-hoc must NEVER be scored against the
# source software labels as actual biological sensitivity/specificity.
for x in classify:x['optical_edge_lower_than_adjacent_row_median']=bool(x['independent_phase_y119_120_edge_ratio_x42_65']<1)
N=len(controls);nlo=int((controls['y119_to_y120_over_local_median_x42_66']<1).sum())
assert N==13 and nlo==7
aggregate={
 'task':'901 E519 source-mechanism strata versus source phase-only diagnostic, NOT biological truth',
 'source_integrity_sha256':receipts,
 'n_source_disputed':len(classify),'n_direct_same_mother_ASC':len(pos),'n_no_direct_ASC':len(neg),
 'direct_ASC_group':{
  'n':len(pos),'candidate_area_min_max':[min(x['small_software_candidate_pixel_area'] for x in pos),max(x['small_software_candidate_pixel_area'] for x in pos)],
  'source_short_candidates_at_y120':sum(x['small_software_candidate_ymin']==120 for x in pos),
  'independent_phase_edge_ratio_less_than_one':sum(x['independent_phase_y119_120_edge_ratio_x42_65']<1 for x in pos),
  'phase_edge_ratio_median':float(np.median([x['independent_phase_y119_120_edge_ratio_x42_65'] for x in pos]))},
 'no_direct_ASC_group':{
  'n':len(neg),'events':[{'gl':x['gl'],'f':x['f_0_based']} for x in neg],
  'candidate_area_values':sorted(x['small_software_candidate_pixel_area'] for x in neg),
  'independent_phase_edge_ratio_less_than_one':sum(x['independent_phase_y119_120_edge_ratio_x42_65']<1 for x in neg),
  'phase_edge_ratio_median':float(np.median([x['independent_phase_y119_120_edge_ratio_x42_65'] for x in neg]))},
 'software_agreement_only_technical_control':{'n':N,'under_1_phase_edge_ratio':nlo,
    'candidate_area_min':int(controls.small_candidate_area_px.min())},
 'not_a_biological_ground_truth_threshold':True,
 'n_independent_biological_experiments_NOT_12_or_37':True,
 'same_gl_six_consecutive_short_windows_are_one_correlated_block':True,
 'Omnipose_true_pretrained_weights_in_current_container':False,
 'Omnipose_actual_inference_this_stage':False,
 'external_blind_expert_ratings_this_stage':0,
 'physical_division_truth': 'UNRESOLVED_12_OF_12',
 'physical_bit_per_second': 'NOT_IDENTIFIED',
 'physical_bit_per_joule': 'NOT_IDENTIFIED',
 'publication':'PUBLICATION_NO_GO'
}

# Write only outputs that are reproducible from declared data. Avoid timestamps.
with (OUT/'E519_all_12_software_dispute_ASC_and_phase_source_strata.csv').open('w',newline='',encoding='utf8') as f:
 w=csv.DictWriter(f,fieldnames=list(classify[0]));w.writeheader();w.writerows(classify)
controls.to_csv(OUT/'E519_all_13_agreed_software_controls_source_image_metrics.csv',index=False,float_format='%.15g')
for file,obj in [('E519_aggregate_group_boundary_and_negatives.json',aggregate),
 ('E519_Pos5_f247_ASC_other_tree_provenance.json',spatial),
 ('E519_Pos5_f294_ASC_distance_negative_control.json',f294)]:
 (OUT/file).write_text(json.dumps(obj,indent=2,ensure_ascii=False,sort_keys=True),encoding='utf8')
(OUT/'E519_explicit_audit_limitations.txt').write_text(
 'E519 is a retrospective cross-layer audit of two related growth lanes from a single source.\n'
 'The direct ASC=1 map requirement is an observed archived optimizer constraint, not proof of who curated or of biological truth.\n'
 'All size, IoU and U-Net values come from MoMA source; original phase-edge gradient was computed without source masks in E517, but ROI selected after seeing y=120.\n'
 'The six f399-415 intervals are temporally dependent, not six independent wet experiments.\n'
 'Original biological true divisions 12/12 unresolved; external pretrained model inference and expert blind ratings still missing.\n'
 'No hypothesis testing or trained classification, no live physical molecular bits or energy. PUBLICATION_NO_GO.\n',encoding='utf8')
print('E519_REPRO_PASS',json.dumps({'n':len(classify),'direct_ASC':len(pos),'no_direct_ASC':len(neg),
 'direct_ASC_area':aggregate['direct_ASC_group']['candidate_area_min_max'],
 'non_ASC_area':aggregate['no_direct_ASC_group']['candidate_area_values'],
 'phase_edge_less_than_one_direct':aggregate['direct_ASC_group']['independent_phase_edge_ratio_less_than_one'],
 'phase_edge_less_than_one_nonASC':aggregate['no_direct_ASC_group']['independent_phase_edge_ratio_less_than_one'],
 'software_agree_controls_below_one':nlo,
 'f247_competing_ASC_mother_same_tree':spatial['same_component_tree_ancestor'],
 'f294_nearest_ASC_frame_gap':f294['closest_any_ASC_record_time_distance_frames']}),flush=True)