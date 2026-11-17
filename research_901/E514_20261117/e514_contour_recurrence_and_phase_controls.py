#!/usr/bin/env python3
"""901 E514 scientific audit: repeat source candidate masks and non-MoMA phase-only controls.

Inputs: prior officially CRC/SHA-authenticated E510 + E512 + E512b immutable
source artifacts, downloaded and re-verified from user's Google Drive.
There is NO true-divisions classifier trained by this script. The source
component overlap and optical heuristics are two DISTINCT evidence tiers.

Generated optical PNGs use only the original TIFF float32 pixels and are
indicated as 'source mask overlays'; overlays are NOT independent labels.
"""
import argparse,collections,csv,hashlib,io,json,math,pathlib,statistics,zipfile
import numpy as np
from scipy.ndimage import gaussian_filter,gaussian_filter1d,label,binary_erosion
from scipy.signal import find_peaks
from skimage.filters import threshold_otsu
from PIL import Image,ImageDraw,ImageFont

P=argparse.ArgumentParser()
P.add_argument('--input',default='/mnt/data',help='Dir holding 6 frozen source ZIPs')
P.add_argument('--output',default='/mnt/data/901_E514_20261117/results')
a=P.parse_args()
INPUT=pathlib.Path(a.input);OUTPUT=pathlib.Path(a.output);OUTPUT.mkdir(parents=True,exist_ok=True)
assert len(str(OUTPUT))<300
known_sha={
'Pos5_GL35':{
 'E510':'745815a73bc1e2eaa01a2cf99f4edd9c530064d10c22bc5421f6ff5189629d8d',
 'E512':'77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6',
 'E512b':'efe81c19a76117c0f9bebe58b8e4525f93c407465312d0402ce7096e9866b384'},
'Pos12_GL22':{
 'E510':'17f80af916f376dfa380b30e8a8b2e8459358a4b66e13e21e23f637aaa08cfa9',
 'E512':'73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2',
 'E512b':'4893c521768bdbcc026e7f6d997b1266f9e70e6869b9bf08e4b429a8aa7ecca9'
}}
SOURCE={'Pos5_GL35':'20221220_rplN_GFP_pl','Pos12_GL22':'20221220_E08_GFP_pl'}
SIX=[399,403,405,409,411,415]
SOURCE_ROI_Y=120
CAMERA_OFFSETS=[-4,-3,-2,-1,0,1,2,4,7]  # E513 verified: E510 stored label is +1
LEVEL_NOTE='SOURCE_SOFTWARE_CANDIDATES_NOT_INDEPENDENT_BIOLOGICAL_GROUND_TRUTH'
ALL={};source_entries=[];frame_data={};mask_by_event={};original_source_receipts=[];keys_by_gl={}

def sha(data):return hashlib.sha256(data).hexdigest()
def zipread(z,suffix):
 found=[n for n in z.namelist() if n.endswith(suffix)]
 assert len(found)==1,(suffix,found)
 return z.read(found[0])

def source_mask(q):
 return set(zip(q['xCoordinates'],q['yCoordinates']))

def rectangular_signal(frames,idx,region):
 if idx not in frames or idx+1 not in frames:return (None,None)
 ys,xs=region
 v,w=frames[idx][ys,xs].astype('float64'),frames[idx+1][ys,xs].astype('float64')
 mean_abs=float(np.mean(np.abs(w-v)))
 corr=float(np.corrcoef(v.ravel(),w.ravel())[0,1]) if np.std(v)>1e-8 and np.std(w)>1e-8 else None
 return mean_abs,corr

for gl,src in SOURCE.items():
 z={}
 for key in ['E510','E512','E512b']:
  path=INPUT/(f'901_{key}_{gl}_'+ ('original_pixels.zip' if key=='E510' else 'full_original_components.zip' if key=='E512' else 'exact_div_to_map.zip'))
  raw=path.read_bytes();assert sha(raw)==known_sha[gl][key],(gl,key,'ZIP_SHA_CHANGED')
  z[key]=zipfile.ZipFile(io.BytesIO(raw));assert z[key].testzip() is None,(gl,key,'ZIP_CRC')
  original_source_receipts.append({'GL':gl,'source_era':key,'archive_name':path.name,'sha256':sha(raw),'bytes':len(raw),'zip_crc_pass':True})
 key=json.loads(zipread(z['E510'],'E510_sealed_event_key_do_not_open_for_blind_review.json'))
 assert len(key)==(22 if gl=='Pos5_GL35' else 15)
 keys_by_gl[gl]=key
 frames={}
 for case in key:
  cid=case['blind_case_id']
  raw=zipread(z['E510'],f'{cid}_exact_float32_phase_frames.npz')
  with np.load(io.BytesIO(raw),allow_pickle=False) as data:
   assert data['offsets'].tolist()==[-3,-2,-1,0,1,2,3,5,8]
   labels=data['frames'];pix=data['phase_pixels']
   assert pix.shape==(9,300,106) and pix.dtype==np.float32
   assert labels.tolist()==[int(case['frame'])+i for i in [-3,-2,-1,0,1,2,3,5,8]]
   for lab,p in zip(labels,pix):
    real_t=int(lab)-1
    if real_t in frames:
     assert np.array_equal(frames[real_t],p),('overlapping source frame mismatch',gl,real_t)
    else:frames[real_t]=p.copy()
 assert len(frames)==(119 if gl=='Pos5_GL35' else 84),(gl,len(frames))
 frame_data[gl]=frames
 root=SOURCE[gl]+'/'+gl+'/'
 receipt=json.loads(z['E512'].read(root+'E512_original_member_SHA_CRC_critical_source_audit.json'))
 forest_bytes=z['E512'].read(root+'main__component_forests.json')
 forest_expected=next(x['original_sha256'] for x in receipt['complete_source_authentication'] if x['version']=='main' and x['kind']=='component_forests.json')
 assert sha(forest_bytes)==forest_expected
 comps=json.loads(forest_bytes)['components'];byid={c['stringId']:c for c in comps};assert len(byid)==len(comps)
 costforests=json.loads(z['E512'].read(root+'cost__component_forests.json'))['components']
 assert len(comps)==len(costforests) and set(byid)=={x['stringId'] for x in costforests}
 assert receipt['map_properties']['main']['source_file_byte_sha256']==receipt['map_properties']['cost']['source_file_byte_sha256']
 def solve(which):
  ss=z['E512'].read(root+which+'__gurobi_model.sol').decode()
  out={}
  for line in ss.splitlines():
   if line.startswith('#') or not line.strip():continue
   k,v=line.split();assert k not in out;out[k]=float(v)
  return out
 main,cost=solve('main'),solve('cost')
 assert main.keys()==cost.keys()
 earlier=json.loads(zipread(z['E512b'],'E512b_exact_software_division_constraint_mechanism.json'))
 assert earlier['gl']==gl and earlier['component_count_each_version']==len(comps)
 for case in sorted(key,key=lambda x:x['frame']):
  if case['class']=='software_nondivision':continue
  t=int(case['frame'])
  cand=[k for k,v in main.items() if k.startswith(f'DivT{t}_') and v>0.5]
  assert len(cand)==1,(gl,t,'ambiguous component')
  div=cand[0];parts=div.split('_');mother_id=parts[1];daughter_ids=parts[2:]
  assert len(daughter_ids)==2 and all(q in byid for q in daughter_ids)
  girls=[byid[q] for q in daughter_ids]
  assert all(int(q['frameNumber'])==t+1 for q in girls)
  small=min(girls,key=lambda q:len(q['xCoordinates']));large=max(girls,key=lambda q:len(q['xCoordinates']))
  sm=source_mask(small);lg=source_mask(large)
  assert not sm&lg,('same-frame hypothesis overlap',gl,t)
  key2=(gl,t)
  mask_by_event[key2]=sm
  att={'GL':gl,'software_event_frame_0based':t,'software_class':case['class'],
       'selected_div_variable':div,'cost_pipeline_also_selects_div':cost[div]>0.5,
       'software_source_component_small':small['stringId'],
       'software_source_component_large':large['stringId'],
       'small_candidate_area_px':len(sm),'large_candidate_area_px':len(lg),
       'small_candidate_ymin':min(y for x,y in sm),'small_candidate_ymax':max(y for x,y in sm),
       'small_candidate_xmin':min(x for x,y in sm),'small_candidate_xmax':max(x for x,y in sm),
       'large_candidate_ymin':min(y for x,y in lg),'large_candidate_ymax':max(y for x,y in lg),
       'source_same_frame_daughter_masks_disjoint':True,
       'small_touches_ROI_y120':min(y for x,y in sm)==SOURCE_ROI_Y,
       'small_candidate_fraction_pixels_phase_above_global_image_otsu':None,
       'source_evidence_level':LEVEL_NOTE}
  source_entries.append(att)
  # Source MoMA masks are used ONLY to sample values from source pixels and
  # never as phase-only independent segmentation labels.
  xx=np.asarray([x for x,y in sm],dtype='int32');yy=np.asarray([y for x,y in sm],dtype='int32')
  if t in frames and t+1 in frames:
   before,after=frames[t][yy,xx],frames[t+1][yy,xx]
   att['phc_at_software_small_mask_t']=float(np.mean(before))
   att['phc_at_software_small_mask_t_plus1']=float(np.mean(after))
   att['phc_abs_mean_change_in_mask']=float(np.mean(np.abs(after-before)))
 # Verify Gurobi frame241/242 case at exact variable level, not one-frame shift inference.
 if gl=='Pos12_GL22':
  common=[k for k,v in main.items() if k.startswith('DivT242_') and v>0.5 and cost[k]>0.5]
  assert len(common)==1
  assert any(x['software_event_frame_0based']==241 and not x['cost_pipeline_also_selects_div'] for x in source_entries if x['GL']==gl)
  (OUTPUT/'E514_pos12_f241_f242_solver_source_disambiguation.json').write_text(json.dumps({
   'f241_main_division_only':True,'f242_identical_active_division_both_versions':common[0],
   'interpretation':'An additional top clipped branch at f241 in standard software; common f242 division is NOT the same shifted event',
   'source':'original Gurobi active variable names and values; not biological ground truth'
  },indent=2),encoding='utf8')
 # All original phase frames used without candidate label spatial restriction
 ALL[gl]={'source_GL':gl,'timepoints_unique':len(frames),'source_division_candidates':len([e for e in key if e['class']!='software_nondivision']),
          'phc_original_zero_based':True,'expected_original_probability_tiff_sha':receipt['map_properties']['main']['source_file_byte_sha256'],
          'component_forest_original_sha':forest_expected}

# Source candidate mask repeated occupancy from f399–415; the numerical mask
# recurrence is source-layer mechanism evidence, NOT independent living cell proof.
iou_rows=[]
for ta,tb in zip(SIX,SIX[1:]):
 a,b=mask_by_event['Pos5_GL35',ta],mask_by_event['Pos5_GL35',tb]
 iou=len(a&b)/len(a|b)
 iou_rows.append({'GL':'Pos5_GL35','event_A':ta,'event_B':tb,'IoU_small_software_candidate_masks':iou,
                  'pixels_intersection':len(a&b),'pixels_union':len(a|b),
                  'basis':'two source-MoMA component masks at distinct correct t+1 frames'})
a1,b1=mask_by_event['Pos5_GL35',SIX[0]],mask_by_event['Pos5_GL35',SIX[-1]]
first_last=len(a1&b1)/len(a1|b1)
assert min(r['IoU_small_software_candidate_masks'] for r in iou_rows)>.75

# Pre-declared image-only optical checks. Not a cell segmentation reference;
# no contours or source model are used to construct these features.
# 1) fixed strip x[46:62], y[120:185], profiles smooth sigma1.3, valleys >=5px;
# 2) 2D bright-component Otsu per GL from all available chronological images,
# with fixed x[42:66], y[108:197], Gaussian sigma0.9 and min-area10.
# Optical feature is considered persistent only if absent at t-2,t-1,t and
# present t+1,t+2,t+4 at positions within 4px (1D) or 6/8px (2D).
profile_prom=[.04,.07,.10,.15]
seg_factors=[.8,1.0,1.2]
phase_records=[];phase_summaries={};lane_events={}
for gl,frames in frame_data.items():
 glkeys=keys_by_gl[gl]
 image_cache={}; otsu_roi=(slice(108,197),slice(42,66))
 pixel_input=np.concatenate([arr[otsu_roi].ravel() for arr in frames.values()])
 threshold=float(threshold_otsu(pixel_input))
 phase_summaries[gl]={'otsu_global_phase_2d_threshold':threshold}
 for t,image in frames.items():
  yy=image[120:185,46:62].astype('float64')
  profile=gaussian_filter1d(np.mean(yy,axis=1),sigma=1.3)
  image_cache[t]={'profile':profile}
  for fac in seg_factors:
   m=gaussian_filter(image[otsu_roi].astype('float64'),sigma=0.9)>(threshold*fac)
   labs,n=label(m,structure=np.ones((3,3),dtype=bool))
   comp2=[]
   for k in range(1,n+1):
    coords=np.nonzero(labs==k)
    if len(coords[0])<10:continue
    # true coords original preprocessed phase TIFF index [0..299,0..105]
    comp2.append({'centroid_y':float(np.mean(coords[0]))+108,'centroid_x':float(np.mean(coords[1]))+42,'area':len(coords[0])})
   image_cache[t][fac]=comp2
 phase_summaries[gl]['image_only_nominal_threshold_source']={'Y':[108,197],'X':[42,66],'n_frames':len(frames),
   'selected_without_using_moma_component_mask':True}
 def profile_vals(t,prom):
  p=image_cache[t]['profile']
  locs,_=find_peaks(-p,prominence=prom,distance=5)
  return [int(j+120) for j in locs]
 def persistent_valley(t,prom):
  vals={k:profile_vals(k,prom) for k in [t-2,t-1,t,t+1,t+2,t+4]}
  cand=[]
  for y in vals[t+1]:
   if any(abs(y-z)<=3 for q in [t-2,t-1,t] for z in vals[q]):continue
   if not all(any(abs(y-z)<=4 for z in vals[q]) for q in [t+2,t+4]):continue
   cand.append(y)
  return cand
 def persistent_2d(t,fac):
  c=image_cache[t+1][fac];before=image_cache[t][fac]
  f2=image_cache[t+2][fac];f4=image_cache[t+4][fac]
  def dist(a,b):return math.hypot(a['centroid_y']-b['centroid_y'],a['centroid_x']-b['centroid_x'])
  cand=[]
  for q in c:
   if not 115<=q['centroid_y']<=185:continue
   if any(dist(q,k)<6 for k in before):continue
   if not any(dist(q,k)<6 for k in f2):continue
   if not any(dist(q,k)<8 for k in f4):continue
   cand.append(q)
  return cand
 for c in glkeys:
  t=int(c['frame'])
  needed={t-2,t-1,t,t+1,t+2,t+4}
  assert needed.issubset(frames),('missing fixed image-time control frames',gl,t,sorted(needed-set(frames)))
  record={'GL':gl,'source_event_frame':t,'software_category_posthoc':c['class'],
          'real_TIFF_zero_based_times_analyzed':sorted(needed)}
  for prm in profile_prom:
   record[f'image_only_persistent_dark_gap_prom_{prm}']=persistent_valley(t,prm)
  for fac in seg_factors:
   record[f'image_only_2d_connected_new_bright_component_fac_{fac}']=persistent_2d(t,fac)
  topdelta,corr=rectangular_signal(frames,t,(slice(120,128),slice(45,64)))
  lowdelta,lowcorr=rectangular_signal(frames,t,(slice(128,154),slice(45,64)))
  record.update({'fixed_roi_top_frame_delta':topdelta,'fixed_roi_top_frame_corr':corr,
    'fixed_adjacent_lower_frame_delta':lowdelta,'fixed_adjacent_lower_frame_corr':lowcorr})
  phase_records.append(record)
 # Summaries must respect correlated windows; do NOT calculate false biological p-values
 for family,thresholds in [('persistent_dark_gap_prom',profile_prom),('2d_connected_new_bright_component_fac',seg_factors)]:
  for p in thresholds:
   field='image_only_'+family+'_'+str(p)
   phase_summaries[gl][field]={s:{'hits':sum(bool(x[field]) for x in phase_records if x['GL']==gl and x['software_category_posthoc']==s),
                   'total':sum(x['GL']==gl and x['software_category_posthoc']==s for x in phase_records)}
        for s in ['software_agreement','disputed','software_nondivision']}

# Generate actual-source corrected-TIFF scientific image diagnostics. Source
# overlays shown only in turquoise at the t+1 images, clearly declared as
# software hypotheses and excluded from optical features.

def visualize(gl,event_set,path):
 frames=frame_data[gl]
 out=Image.new('RGB',(5*248,len(event_set)*228),(240,240,240))
 draw=ImageDraw.Draw(out)
 for row,t in enumerate(event_set):
  for col,tt in enumerate([t-2,t-1,t,t+1,t+2]):
   # All displayed images are unaltered numerical source phase arrays before
   # display scaling; x[36:76],y[108:165] fixed per GL and event category.
   src=frames[tt][108:165,36:76]
   pix=np.uint8(np.clip((src+0.11)/(1.8+0.11),0,1)*255)
   image=Image.fromarray(pix,mode='L').resize((200,171),Image.Resampling.NEAREST).convert('RGB')
   if col==3 and (gl,t) in mask_by_event:
    # Annotate the SOURCE segmentation hypothesis; not independent pixel evidence.
    mask=np.zeros((57,40),dtype=bool)
    for x,y in mask_by_event[(gl,t)]:
     if 108<=y<165 and 36<=x<76:mask[y-108,x-36]=True
    bd=mask & (~binary_erosion(mask))
    m=np.repeat(np.repeat(bd,3,axis=0),5,axis=1)
    aimg=np.asarray(image).copy();aimg[m]=[0,255,240];image=Image.fromarray(aimg)
   x=col*248+6;y=row*228+24
   out.paste(image,(x,y))
   draw.text((x,y-18),f'source t={tt}'+ (' [MoMA mask]' if col==3 else ''),fill=(20,20,20))
   # yellow input cell detection ROI y=120 as pre-registered original config
   draw.line(((x,y+(120-108)*3),(x+200,y+(120-108)*3)),fill=(255,220,0),width=1)
  draw.text((1230,row*228+60),str(t),fill=(0,0,0))
 out.save(path)

visualize('Pos5_GL35',SIX,OUTPUT/'E514_Pos5_six_divisions_correct_time_original_phase_with_MoMA_mask_annotated.png')
visualize('Pos12_GL22',[241,242,246],OUTPUT/'E514_Pos12_f241_f242_f246_original_phase_with_MoMA_mask_annotated.png')

# Write machine-readable event and control audit
with (OUTPUT/'E514_source_MoMA_candidate_daughter_geometry.csv').open('w',newline='',encoding='utf-8') as f:
 cols=sorted({key for e in source_entries for key in e});w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(source_entries)
with (OUTPUT/'E514_six_repeat_top_small_daughter_iou.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=list(iou_rows[0]));w.writeheader();w.writerows(iou_rows)
with (OUTPUT/'E514_image_only_37_case_controls.json').open('w',encoding='utf-8') as f:
 json.dump({'protocol':{'x_phase_profile':[46,62],'y_phase_profile':[120,185],
   'profile_smoothing_sigma':1.3,'profile_peak_min_distance':5,'prominence_sensitivity':profile_prom,
   'phase_only_2d_global_otsu_x':[42,66],'phase_only_2d_global_otsu_y':[108,197],
   'phase_only_2d_sigma':0.9,'min_bright_component_area':10,
   'phase_2d_otsu_sensitivity_factors':seg_factors,'fixed_raw_source_phase_TIFF_time_zero_based':True,
   'independent_cell_segmentation_ground_truth':False,
   'avoid_threshold_selection_based_on_software_event_class':True},
   'gl_summary':phase_summaries,'cases':phase_records},f,ensure_ascii=False,indent=2)

summary={
 'E514_type':'original phase-image independent signal controls PLUS separate existing MoMA hypothesis reanalysis',
 'original_archive_integrity':original_source_receipts,
 'E513_corrected_time_axis_reused':True,
 'registered_focus_Pos5_frames':SIX,
 'registered_focus_Pos12_frames':[241,242,246],
 'fully_separate_software_source_vs_phase_image_evidence':True,
 'Pos5_six_source_moma_small_candidate_areas':[next(x for x in source_entries if x['GL']=='Pos5_GL35' and x['software_event_frame_0based']==t)['small_candidate_area_px'] for t in SIX],
 'Pos5_six_source_moma_small_candidate_minY':[next(x for x in source_entries if x['GL']=='Pos5_GL35' and x['software_event_frame_0based']==t)['small_candidate_ymin'] for t in SIX],
 'Pos5_six_pairwise_neighbor_MoMA_mask_iou':[x['IoU_small_software_candidate_masks'] for x in iou_rows],
 'Pos5_six_first_last_MoMA_mask_iou':first_last,
 'Pos5_six_source_IoU_no_biological_phenotype_inferred':True,
 'Pos12_f242_division_is_both_software_pipelines_same_exact_solver_variable':True,
 'image_only_37_case_controls':phase_summaries,
 'Pos5_six_true_source_phase_top_patch_mean_absolute_changes':[
    next(x for x in phase_records if x['GL']=='Pos5_GL35' and x['source_event_frame']==t)['fixed_roi_top_frame_delta'] for t in SIX],
 'Pos5_six_true_source_phase_top_patch_frame_correlations':[
    next(x for x in phase_records if x['GL']=='Pos5_GL35' and x['source_event_frame']==t)['fixed_roi_top_frame_corr'] for t in SIX],
 'interpretation':'Repetitive y=120 source-MoMA candidate masks are software tracking/ROI signature. Image-only brightness island and valley methods have inadequate sensitivity on software-agreement controls. Neither determines physical truth.',
 'software_division_disputed_n':12,
 'physical_division_ground_truth_unresolved_n':12,
 'publication':'PUBLICATION_NO_GO',
 'physical_biological_molecular_bit':'NOT_IDENTIFIED',
 'source_code':'e514_contour_recurrence_and_phase_controls.py',
 'no_prior_results_overwritten':True,
}
(OUTPUT/'E514_frozen_integrity_and_findings_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf8')
print('E514_ANALYSIS_DONE',json.dumps({
 'E510_frame_reuse_unique_P5':len(frame_data['Pos5_GL35']),
 'E510_frame_reuse_unique_P12':len(frame_data['Pos12_GL22']),
 'six_small_areas':summary['Pos5_six_source_moma_small_candidate_areas'],
 'six_iou':summary['Pos5_six_pairwise_neighbor_MoMA_mask_iou'],
 'six_first_last_iou':first_last,
 'phase_only_2d':{gl:phase_summaries[gl]['image_only_2d_connected_new_bright_component_fac_1.0'] for gl in phase_summaries},
 'physical_ground_truth':'UNRESOLVED'},ensure_ascii=False),flush=True)