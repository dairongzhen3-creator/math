#!/usr/bin/env python3
"""901 E513: audit the time-index offset against source phase pixels and CellStats.

Inputs: immutable E508 source CSV, E510 original float32 phase frame NPZ, E512
original component coordinates and SHA receipts. No new microscopic pixels or
human ground truth are invented. Source E510 frame LABEL was one too high:
the saved array from label L was source TIFF time INDEX L-1.
All existing E510 files remain unmodified; corrected metrics are new artifacts.
"""
import argparse, collections, csv, hashlib, json, pathlib, statistics, numpy as np

A=argparse.ArgumentParser()
A.add_argument('--gl',required=True,choices=['Pos5_GL35','Pos12_GL22'])
A.add_argument('--e508',default='input_E508')
A.add_argument('--e510',default='input_E510')
A.add_argument('--e512',default='input_E512')
a=A.parse_args()
CONFIG={
'Pos5_GL35':'20221220_rplN_GFP_pl',
'Pos12_GL22':'20221220_E08_GFP_pl'
}
arc=CONFIG[a.gl]
source=pathlib.Path(a.e512)/arc/a.gl
phasefolder=pathlib.Path(a.e510)/a.gl
statfolder=pathlib.Path(a.e508)/'e507_alternate_QA'/arc
sha_receipts=pathlib.Path(a.e508)/'e508_source_conflicts'/arc/'E508_original_source_csv_sha256.json'
OUT=pathlib.Path('e513_zero_based_time_audit')/a.gl
OUT.mkdir(parents=True,exist_ok=True)

receipt=json.loads((source/'E512_original_member_SHA_CRC_critical_source_audit.json').read_text())
forest_path=source/'main__component_forests.json'
source_sha=next(x['original_sha256'] for x in receipt['complete_source_authentication']
    if x['version']=='main' and x['kind']=='component_forests.json')
assert hashlib.sha256(forest_path.read_bytes()).hexdigest()==source_sha
fobj=json.loads(forest_path.read_text())
assert 'components' in fobj and fobj['version']=='0.1.0'
comps=fobj['components']
by_geom=collections.defaultdict(list)
for c in comps:
    xx,yy=c['xCoordinates'],c['yCoordinates']
    assert len(xx)==len(yy)>0 and min(xx)>=0 and max(xx)<106
    assert min(yy)>=0 and max(yy)<300
    by_geom[(int(c['frameNumber']),len(xx),min(yy),max(yy))].append(c)
csv_path=statfolder/(a.gl+'_main_CellStats.csv')
meta=json.loads(sha_receipts.read_text())
met=[x for x in meta if x['GL']==a.gl and x['version']=='main']
assert len(met)==1 and hashlib.sha256(csv_path.read_bytes()).hexdigest()==met[0]['sha256']
rows=list(csv.DictReader(x for x in csv_path.read_text(encoding='utf-8-sig').splitlines() if x.strip() and not x.startswith('#')))

# Independently verify the old pipeline's saved numerical slices' provenance,
# using exact original labels, with no new/source-model-generated phase pixels.
key=json.loads((phasefolder/'E510_sealed_event_key_do_not_open_for_blind_review.json').read_text())
key_by_id={x['blind_case_id']:x for x in key}
assert len(key_by_id)==len(key)
phase_by_t={}
frame_sha={}
records=[]
OFFSET_OLD=[-3,-2,-1,0,1,2,3,5,8]
OFFSET_TRUE=[x-1 for x in OFFSET_OLD]
for case in key:
    cid=case['blind_case_id']; t=int(case['frame'])
    npz=phasefolder/(cid+'_exact_float32_phase_frames.npz')
    with np.load(npz,allow_pickle=False) as d:
        labels=d['frames'].astype('int32').tolist()
        saved=d['phase_pixels']
        assert labels==[t+i for i in OFFSET_OLD] and d['offsets'].tolist()==OFFSET_OLD
        assert saved.shape==(9,300,106) and saved.dtype==np.float32
        for label,image in zip(labels,saved):
            time_index=int(label)-1
            if time_index in phase_by_t:
                assert np.array_equal(phase_by_t[time_index],image),('INCONSISTENT_SAME_TIME_PHASE',a.gl,time_index)
            else:
                phase_by_t[time_index]=image.copy()
                frame_sha[time_index]=hashlib.sha256(image.tobytes()).hexdigest()
        # In the immutable E510 cache, x[3] was TIFF f-1, x[4] was TIFF f.
        # Correct physical candidate division f->f+1 is x[4] -> x[5].
        # Use five other adjacent transitions as within-card controls.
        def calc(region):
            sel=saved[:,region[0],region[1]].astype('float64')
            steps=[float(np.mean(np.abs(sel[j+1]-sel[j]))) for j in range(6)]
            assert all(np.isfinite(s) for s in steps)
            old=steps[3]; corrected=steps[4]
            old_base=np.median([v for i,v in enumerate(steps) if i!=3])
            cor_base=np.median([v for i,v in enumerate(steps) if i!=4])
            return {'old_candidate_mean_abs_delta':old,
                    'old_ratio':float(old/old_base) if old_base>0 else None,
                    'corrected_candidate_mean_abs_delta':corrected,
                    'corrected_ratio':float(corrected/cor_base) if cor_base>0 else None,
                    'old_is_top2_of_six':int(sum(x>old for x in steps)<2),
                    'corrected_is_top2_of_six':int(sum(x>corrected for x in steps)<2),
                    'six_adjacent_deltas':steps}
        r=calc((slice(80,255),slice(None)))
        central=calc((slice(80,255),slice(20,70)))
        top=calc((slice(110,150),slice(None)))
        records.append({'GL':a.gl,'blind_case_id':cid,'software_class':case['class'],
                        'source_cellstats_frame_index':t,
                        'corrected_original_tiff_zero_based_event_frame':t,
                        'E510_saved_label_center':t,
                        'E510_originally_displayed_TIFF_index':t-1,
                        'E513_true_candidate_TIFF_index':t,
                        'E510_stored_labels':labels,
                        'corrected_offsets_relative_to_source_event':OFFSET_TRUE,
                        'area_full':r,'area_center':central,'area_y110_149':top})

# Official source CellStats uses bounding-box pixel sum (MoMA Imglib2Utils
# getIntensityTotal accepts LabelRegion as Interval and Views.interval makes
# the bounding rectangle). Reproduce the actual CSV float to 5 decimals.
# No use of software division or control labels in selecting intensity checks.
sum_checks=[]
nonunique=0; unmatched=0
for row in rows:
    f=int(row['frame'])
    if f not in phase_by_t or f-1 not in phase_by_t:continue
    k=(f,int(row['area_px']),int(row['bbox_top_px']),int(row['bbox_bottom_px']))
    possible=by_geom.get(k,[])
    if len(possible)==0:unmatched+=1;continue
    if len(possible)>1:nonunique+=1;continue
    comp=possible[0]
    xx=comp['xCoordinates']; yy=comp['yCoordinates']
    sy=slice(min(yy),max(yy)+1);sx=slice(min(xx),max(xx)+1)
    correct_sum=float(np.sum(phase_by_t[f][sy,sx],dtype=np.float64))
    prev_sum=float(np.sum(phase_by_t[f-1][sy,sx],dtype=np.float64))
    original=float(row['phc_total_intensity_au'])
    sum_checks.append({'gl':a.gl,'cell_ID':int(row['cell_ID']),'source_frame':f,
       'component_ID':comp['stringId'],'original_CSV_phc_intensity':original,
       'correct_TIFF_f_bbox_sum':correct_sum,'incorrect_TIFF_f_minus_1_bbox_sum':prev_sum,
       'correct_absolute_error':abs(correct_sum-original),
       'incorrect_absolute_error':abs(prev_sum-original),
       'correct_matches_CSV_5_decimals':round(correct_sum,5)==round(original,5),
       'source_component_bbox_size':[max(xx)-min(xx)+1,max(yy)-min(yy)+1],
       'cellstats_area':int(row['area_px'])})
assert len(sum_checks)>=25,(a.gl,len(sum_checks),unmatched,nonunique)
correct_err=np.asarray([x['correct_absolute_error'] for x in sum_checks])
wrong_err=np.asarray([x['incorrect_absolute_error'] for x in sum_checks])
assert float(np.max(correct_err))<0.00025,(a.gl,float(np.max(correct_err)))
assert float(np.median(wrong_err))>0.02,(a.gl,float(np.median(wrong_err)))
assert (correct_err < wrong_err).mean()>=.90,(a.gl,(correct_err < wrong_err).mean())

with (OUT/'E513_original_cell_phc_export_roundtrip.csv').open('w',newline='',encoding='utf8') as h:
    w=csv.DictWriter(h,fieldnames=list(sum_checks[0].keys()));w.writeheader();w.writerows(sum_checks)
with (OUT/'E513_corrected_software_window_metrics.csv').open('w',newline='',encoding='utf8') as h:
    cols=['GL','blind_case_id','software_class','source_cellstats_frame_index',
          'old_ratio','corrected_ratio','old_is_top2_of_six','corrected_is_top2_of_six',
          'old_center_ratio','corrected_center_ratio','old_boundary_ratio','corrected_boundary_ratio']
    w=csv.DictWriter(h,fieldnames=cols);w.writeheader()
    for z in records:
        w.writerow({'GL':a.gl,'blind_case_id':z['blind_case_id'],'software_class':z['software_class'],
          'source_cellstats_frame_index':z['source_cellstats_frame_index'],
          'old_ratio':z['area_full']['old_ratio'],'corrected_ratio':z['area_full']['corrected_ratio'],
          'old_is_top2_of_six':z['area_full']['old_is_top2_of_six'],
          'corrected_is_top2_of_six':z['area_full']['corrected_is_top2_of_six'],
          'old_center_ratio':z['area_center']['old_ratio'],'corrected_center_ratio':z['area_center']['corrected_ratio'],
          'old_boundary_ratio':z['area_y110_149']['old_ratio'],'corrected_boundary_ratio':z['area_y110_149']['corrected_ratio']})
by_class=collections.defaultdict(list)
for z in records:by_class[z['software_class']].append(z)
group={}
for name,arr in sorted(by_class.items()):
    group[name]={'n':len(arr),
        'old_pixel_ratio_median':statistics.median(x['area_full']['old_ratio'] for x in arr),
        'corrected_pixel_ratio_median':statistics.median(x['area_full']['corrected_ratio'] for x in arr),
        'corrected_pixel_center_ratio_median':statistics.median(x['area_center']['corrected_ratio'] for x in arr),
        'corrected_pixel_boundary_ratio_median':statistics.median(x['area_y110_149']['corrected_ratio'] for x in arr),
        'old_top2':sum(x['area_full']['old_is_top2_of_six'] for x in arr),
        'corrected_top2':sum(x['area_full']['corrected_is_top2_of_six'] for x in arr)}
audit={'gl':a.gl,'source_archive':arc+'.zip',
       'source_e508_csv_sha':met[0]['sha256'],
       'source_e512_components_sha':source_sha,
       'source_e512_frozen_probability_tiff_sha':receipt['map_properties']['main']['source_file_byte_sha256'],
       'original_e510_frame_mapping':{'E510_label_L_is_TIFF_index': 'L-1',
           'mistaken_event_tiff_index':'f-1','correct_event_tiff_index':'f',
           'old_offsets':OFFSET_OLD,'correct_offsets':OFFSET_TRUE,
           'old_npz_immutable':True},
       'complete_original_source_phase_tiff_already_authenticated_in_E510':True,
       'same_TIFF_frame_overlapping_cards_bitwise_equal':True,
       'unique_source_phase_frames_in_E510_cache':len(phase_by_t),
       'original_source_export_roundtrip_n':len(sum_checks),
       'unmatched_ambiguous_skipped':nonunique,'unmatched_missing_component':unmatched,
       'max_correct_bbox_phc_abs_error':float(np.max(correct_err)),
       'median_correct_bbox_phc_abs_error':float(np.median(correct_err)),
       'median_f_minus_1_bbox_phc_abs_error':float(np.median(wrong_err)),
       'number_with_exact_5dp_roundtrip':sum(x['correct_matches_CSV_5_decimals'] for x in sum_checks),
       'correct_error_smaller_than_wrong_count':int(sum(correct_err < wrong_err)),
       'by_software_class':group,
       'cases':records,
       'source_support':'public original MoMA CellStatsExporter.java frame=timeStep and ImageProvider.getColorChannelAtTime(0,timeStep), empirically checked pixel bbox sums',
       'physical_division_ground_truth':'UNRESOLVED',
       'molecular_information_bit_per_joule':'NOT_IDENTIFIED',
       'publication':'PUBLICATION_NO_GO'}
(OUT/'E513_roundtrip_and_corrected_event_metrics.json').write_text(
    json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
(OUT/'E513_unique_TIFF_frame_sha256.json').write_text(json.dumps({
    'gl':a.gl,'mapping':'true TIFF zero-based','sha':frame_sha
},indent=2),encoding='utf8')
print('E513_FRAME_ZERO_BASED_SOURCE_ROUNDTRIP_PASS',json.dumps({
    'gl':a.gl,'cells':len(sum_checks),'max_correct_error':float(np.max(correct_err)),
    'median_wrong_error':float(np.median(wrong_err)),
    'correct_better':int(sum(correct_err<wrong_err)),
    'exact_5dp':audit['number_with_exact_5dp_roundtrip'],
    'unique_pixels':len(phase_by_t),'class_medians':group},ensure_ascii=False),flush=True)
