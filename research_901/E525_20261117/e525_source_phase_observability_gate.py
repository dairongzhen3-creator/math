#!/usr/bin/env python3
"""901 E525: zero-based source phase temporal contrast at 37 frozen cases.

This is a PHASE-ONLY low-information observational control, NOT a division
classifier, independent human annotation, or probability of cell fission.
E524 f960 and f45 are NOT in the frozen E510 set; no image is fabricated.
"""
import argparse,csv,hashlib,io,json,pathlib,statistics,zipfile,collections
import numpy as np
P=argparse.ArgumentParser();P.add_argument('--base',default='/mnt/data');P.add_argument('--output',required=True);args=P.parse_args()
B=pathlib.Path(args.base);O=pathlib.Path(args.output);O.mkdir(parents=True,exist_ok=True)
SHA={'Pos5_GL35':'745815a73bc1e2eaa01a2cf99f4edd9c530064d10c22bc5421f6ff5189629d8d','Pos12_GL22':'17f80af916f376dfa380b30e8a8b2e8459358a4b66e13e21e23f637aaa08cfa9'}
NEAR={'Pos5_GL35':[960],'Pos12_GL22':[45,118]}

def sha(b):return hashlib.sha256(b).hexdigest()
rows=[];summary={}
for gl in SHA:
 raw=(B/f'901_E510_{gl}_original_pixels.zip').read_bytes();assert sha(raw)==SHA[gl]
 with zipfile.ZipFile(io.BytesIO(raw)) as z:
  assert z.testzip() is None
  names=[n for n in z.namelist() if n.endswith('E510_sealed_event_key_do_not_open_for_blind_review.json')]
  assert len(names)==1
  key=json.loads(z.read(names[0]));assert len(key)==(22 if gl=='Pos5_GL35' else 15)
  frames={}
  for case in key:
   cid=case['blind_case_id'];f=int(case['frame'])
   found=[n for n in z.namelist() if n.endswith(cid+'_exact_float32_phase_frames.npz')];assert len(found)==1
   with np.load(io.BytesIO(z.read(found[0])),allow_pickle=False) as v:
    px=v['phase_pixels'];labels=v['frames'];offsets=v['offsets'].tolist()
    assert px.shape==(9,300,106) and px.dtype==np.float32
    assert offsets==[-3,-2,-1,0,1,2,3,5,8]
    assert labels.tolist()==[f+t for t in offsets]
    for l,img in zip(labels,px):
     t=int(l)-1  # E513 source-confirmed TIFF zero-based time mapping
     if t in frames:assert np.array_equal(frames[t],img)
     else:frames[t]=img.copy()
  assert len(frames)==(119 if gl=='Pos5_GL35' else 84)
  assert all(t not in frames for focus in NEAR[gl] for t in (focus-1,focus,focus+1)), ('new E524/E525 timing images unexpectedly present',gl)
  for case in key:
   f=int(case['frame']);a=frames[f-1][80:255,:].astype(np.float64);b=frames[f][80:255,:].astype(np.float64);c=frames[f+1][80:255,:].astype(np.float64)
   d_early=float(np.mean(np.abs(b-a)));d_saved=float(np.mean(np.abs(c-b)))
   assert d_saved>0
   rows.append({'GL':gl,'software_event_t':f,'technical_software_class':case['class'],
         'relative_early_frame_zero_based':f-1,'relative_saved_frame_zero_based':f,
         'mean_abs_phase_delta_before_t':format(d_early,'.11f'),'mean_abs_phase_delta_t_to_t1':format(d_saved,'.11f'),
         'early_over_saved_phase_change':format(d_early/d_saved,'.9f'),
         'view': 'full_GL_x_0_105_y_80_254','SOURCE_PHASE_ONLY_NOT_BIOLOGICAL_DIVISION_TRUTH':True})
 counts=collections.Counter(k['class'] for k in key)
 assert counts==({'software_agreement':4,'disputed':9,'software_nondivision':9} if gl=='Pos5_GL35' else {'software_agreement':9,'disputed':3,'software_nondivision':3})
 summary[gl]={'zip_sha256':sha(raw),'original_zipped_member_crc_all_pass':True,'distinct_true_zero_based_source_frames':len(frames),'all_technical_cases':len(key),'technical_classes':dict(counts),'software_near_tie_original_div_t_absent_from_E510':NEAR[gl],
              'these_near_tie_frames_present_in_E510':False,
              'class_descriptive_median_early_over_saved_phase_change':{cls:statistics.median(float(x['early_over_saved_phase_change']) for x in rows if x['GL']==gl and x['technical_software_class']==cls) for cls in counts}}
rows.sort(key=lambda r:(r['GL'],r['software_event_t']))
with (O/'E525_37_true_phase_observability_vs_1frame_timing.csv').open('w',newline='',encoding='utf8') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary.update({'phase_controls_not_validated_for_ground_truth':True,'all_37_not_independent_biological_replicates':True,'E524_near_ties_missing_raw_observation':True,'future_physical_truth_UNRESOLVED':True})
(O/'E525_phase_only_data_coverage_and_control_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf8')
print('E525_SOURCE_PHASE_ONLY_COVERAGE_PASS',json.dumps(summary,ensure_ascii=False,sort_keys=True),flush=True)
