#!/usr/bin/env python3
"""901 E515: immutable author-source temporal-candidate audit + software-blind phase microscopy packet.

Reuses ONLY official BioImage Archive source artifacts E510/E512/E512b with their frozen SHA.
Analysis has two incompatible evidence tiers:
A) MoMA candidates and probability map: source/software hypotheses only, NOT biological ground truth.
B) Raw preprocessed phase-contrast images: no MoMA geometry/labels in reviewer media.
All event classes remain hidden from blind packet; separate sealed-key file never bundled in its ZIP.
No model trained on MoMA output; no biological error-rate or molecular bits inferred.
"""
import argparse,collections,csv,hashlib,io,itertools,json,math,pathlib,statistics,zipfile
import numpy as np
from PIL import Image,ImageDraw,ImageFont

P=argparse.ArgumentParser();P.add_argument('--input',default='/mnt/data');P.add_argument('--output',default='/mnt/data/901_E515_20261117/replay')
args=P.parse_args();SRC=pathlib.Path(args.input);OUT=pathlib.Path(args.output);OUT.mkdir(parents=True,exist_ok=True)
PUBLIC=OUT/'blind_media';PUBLIC.mkdir(exist_ok=True);PRIVATE=OUT/'private_key';PRIVATE.mkdir(exist_ok=True)
known_sha={
 'Pos5_GL35':{'E510':'745815a73bc1e2eaa01a2cf99f4edd9c530064d10c22bc5421f6ff5189629d8d','E512':'77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6','E512b':'efe81c19a76117c0f9bebe58b8e4525f93c407465312d0402ce7096e9866b384'},
 'Pos12_GL22':{'E510':'17f80af916f376dfa380b30e8a8b2e8459358a4b66e13e21e23f637aaa08cfa9','E512':'73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2','E512b':'4893c521768bdbcc026e7f6d997b1266f9e70e6869b9bf08e4b429a8aa7ecca9'}
}
ARCH={'Pos5_GL35':'20221220_rplN_GFP_pl','Pos12_GL22':'20221220_E08_GFP_pl'}
expected={'Pos5_GL35':{'disputed':9,'software_agreement':4,'software_nondivision':9},'Pos12_GL22':{'disputed':3,'software_agreement':9,'software_nondivision':3}}
SOURCE_GL=[];phase_by_gl={};cases=[];source_events=[];integrity=[];components_by_gl={}
OFFSETS=[-2,-1,0,1,2,4]

def check_zip(gl,kind):
  suffix={'E510':'original_pixels.zip','E512':'full_original_components.zip','E512b':'exact_div_to_map.zip'}[kind]
  path=SRC/(f'901_{kind}_{gl}_'+suffix)
  assert path.is_file(),path
  raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
  assert digest==known_sha[gl][kind],('frozen_source_zip_modified',gl,kind)
  z=zipfile.ZipFile(io.BytesIO(raw));assert z.testzip() is None
  integrity.append({'GL':gl,'source_phase':kind,'filename':path.name,'sha256':digest,'size':len(raw),'all_zip_crc_pass':True})
  return z

def member(z,suffix):
  names=[n for n in z.namelist() if n.endswith(suffix)]
  assert len(names)==1,(suffix,names)
  return z.read(names[0])

def mask(c):return frozenset(zip(c['xCoordinates'],c['yCoordinates']))

def valsort(v):
  return int(v['frame'])

for gl in ARCH:
  Z={kind:check_zip(gl,kind) for kind in ['E510','E512','E512b']}
  cr=ARCH[gl]+'/'+gl+'/'
  r=json.loads(member(Z['E512'],'E512_original_member_SHA_CRC_critical_source_audit.json'))
  assert r['gl']==gl and r['archive']==ARCH[gl]+'.zip'
  assert r['map_properties']['main']['source_file_byte_sha256']==r['map_properties']['cost']['source_file_byte_sha256']
  forestbytes=Z['E512'].read(cr+'main__component_forests.json')
  forestsha=next(x['original_sha256'] for x in r['complete_source_authentication'] if x['version']=='main' and x['kind']=='component_forests.json')
  assert hashlib.sha256(forestbytes).hexdigest()==forestsha
  comp=json.loads(forestbytes)['components'];compmap={c['stringId']:c for c in comp};assert len(compmap)==len(comp)
  byframe=collections.defaultdict(list)
  for c in comp:
    # Use exact original source component 2D pixel list.
    c['_m']=mask(c);byframe[int(c['frameNumber'])].append(c)
  components_by_gl[gl]=byframe
  s={}
  for line in Z['E512'].read(cr+'main__gurobi_model.sol').decode().splitlines():
    if not line.strip() or line.startswith('#'):continue
    k,v=line.split();assert k not in s;s[k]=float(v)
  # Under E513 correction the E512 saved frame number LABEL L is time index L-1.
  src=Z['E512'].read(cr+'main__frozen_probability_frames.npz')
  with np.load(io.BytesIO(src),allow_pickle=False) as n:
    pframes={int(f)-1:a.copy() for f,a in zip(n['original_frame_numbers'],n['probability_values'])}
  assert len(pframes)==(81 if gl=='Pos5_GL35' else 66)
  key=json.loads(member(Z['E510'],'E510_sealed_event_key_do_not_open_for_blind_review.json'))
  assert collections.Counter(c['class'] for c in key)==expected[gl]
  original={}
  for case in key:
    prefix=case['blind_case_id'];npz=member(Z['E510'],f'{prefix}_exact_float32_phase_frames.npz')
    with np.load(io.BytesIO(npz),allow_pickle=False) as dat:
      ims=dat['phase_pixels'];t=int(case['frame']);oldf=dat['frames']
      assert ims.shape==(9,300,106) and ims.dtype==np.float32
      assert oldf.tolist()==[t+i for i in [-3,-2,-1,0,1,2,3,5,8]]
      for label,pix in zip(oldf,ims):
        actual=int(label)-1
        if actual in original:assert np.array_equal(original[actual],pix),('same_time_mismatch',gl,actual)
        else:original[actual]=pix.copy()
  assert len(original)==(119 if gl=='Pos5_GL35' else 84)
  phase_by_gl[gl]=original
  # Fixed GL-wide grayscale from all E510 unique, genuine phase pixels, not MoMA labels.
  sample=np.concatenate([im[::3,::2].ravel() for im in original.values()])
  graylo,grayhi=map(float,np.percentile(sample,[1,99]));assert grayhi>graylo
  SOURCE_GL.append({'gl':gl,'archive':ARCH[gl]+'.zip','source_PhC_time_zero_based':True,
    'source_float_2D_frames':len(original),'source_probability_timepoints':len(pframes),
    'GL_image_normalization_p01':graylo,'GL_image_normalization_p99':grayhi,
    'source_forest_sha':forestsha,'original_prob_tiff_sha':r['map_properties']['main']['source_file_byte_sha256']})
  for case in key:
    f=int(case['frame']);classif=case['class']
    public_id=hashlib.sha256(('901E515_fixed_order_phc_only_blind_20261117_'+gl+'_'+str(f)+'_'+case['blind_case_id']).encode()).hexdigest()[:10].upper()
    assert all(f+d in original for d in OFFSETS),(gl,f)
    cases.append({'public_id':'R'+public_id,'gl':gl,'f':f,'class':classif,'old_case_id':case['blind_case_id'],
       'scale':[graylo,grayhi]})
    if classif=='software_nondivision':continue
    selected=[k for k,v in s.items() if k.startswith(f'DivT{f}_') and v>.5]
    assert len(selected)==1,(gl,f,selected)
    ds=selected[0].split('_')[2:];assert len(ds)==2 and all(k in compmap for k in ds)
    kids=[compmap[k] for k in ds]
    small=min(kids,key=lambda c:len(c['_m']));large=max(kids,key=lambda c:len(c['_m']))
    assert int(small['frameNumber'])==f+1 and int(large['frameNumber'])==f+1
    n=small['_m'];xx=np.array(small['xCoordinates'],dtype=np.int32);yy=np.array(small['yCoordinates'],dtype=np.int32)
    assert f+1 in pframes,(gl,f,'missing_correct_probability_frame')
    confidence=float(pframes[f+1][yy,xx].mean())
    prior=[]
    for dt in range(-5,0):
      t=f+1+dt
      cand=[(len(n&m['_m'])/len(n|m['_m']),m['stringId']) for m in byframe[t] if len(n)*.5<=len(m['_m'])<=len(n)*2]
      best=max(cand,key=lambda x:x[0],default=(0,None))
      prior.append({'relative_to_small_candidate_frame':dt,'candidate_frame':t,'best_source_mask_iou':best[0],'closest_prior_component_id':best[1]})
    topmatch=max(prior,key=lambda v:v['best_source_mask_iou'])
    source_events.append({'gl':gl,'frame':f,'software_label':classif,'original_neutral_case_id':case['blind_case_id'],
      'software_division_variable':selected[0],'source_small_candidate_id':small['stringId'],
      'small_candidate_area_px':len(n),'small_candidate_ymin':int(min(yy)),'small_candidate_ymax':int(max(yy)),
      'small_candidate_xmin':int(min(xx)),'small_candidate_xmax':int(max(xx)),
      'small_MoMA_probtiff_mean_uncalibrated':confidence,
      'small_MoMA_probtiff_median_uncalibrated':float(np.median(pframes[f+1][yy,xx])),
      'small_MoMA_probtiff_gt_0_5_fraction':float(np.mean(pframes[f+1][yy,xx]>.5)),
      'best_past_5_frames_source_component_iou':topmatch['best_source_mask_iou'],
      'best_past_source_frame':topmatch['candidate_frame'],
      'best_past_source_component_id':topmatch['closest_prior_component_id'],
      'all_five_prior_component_matches':prior,
      'uses_source_MoMA_hypotheses':True,
      'is_independent_biological_truth':False})

# No software features, candidate source masks or event categories enter image rendering.
for c in cases:
  gl=c['gl'];f=c['f'];src=phase_by_gl[gl];lo,hi=c['scale'];caseid=c['public_id']
  def image_at(t,subset,scale):
    data=src[t][subset]
    a=np.uint8(np.clip((data-lo)/(hi-lo),0,1)*255)
    h,w=a.shape
    return Image.fromarray(a,mode='L').resize((w*scale,h*scale),Image.Resampling.NEAREST).convert('RGB')
  # Full 300x106 field, fixed contrast for both software classes. Same size for all reviewers.
  for slot,sub,scale in [('full',np.s_[:,:],2),('focus',np.s_[100:185,33:78],4)]:
    tile=[image_at(f+o,sub,scale) for o in OFFSETS]
    width=tile[0].width;height=tile[0].height
    sheet=Image.new('RGB',(width*len(tile),height+34),(242,242,242));d=ImageDraw.Draw(sheet)
    for j,png in enumerate(tile):
      sheet.paste(png,(j*width,34))
      d.text((j*width+6,8),f'{caseid}  relative frame {OFFSETS[j]:+d}',fill=(0,0,0))
    path=PUBLIC/f'{caseid}_{slot}_original_PhC.png';sheet.save(path)

# Source-free pixel optical time-series control: fixed strips selected because y120 was
# already a source-derived hypothesis; they are NOT independent of ROI choice.
strip_summary={}
GL='Pos5_GL35';fr=phase_by_gl[GL];focus_frames=[399+1,403+1,405+1,409+1,411+1,415+1]
assert all(t in fr for t in focus_frames)
for key,y0,y1,x0,x1 in [('edge_120',120,128,45,64),('near_128',128,136,45,64),('near_144',144,152,45,64)]:
  cc=[]
  for t1,t2 in itertools.combinations(focus_frames,2):
    u=fr[t1][y0:y1,x0:x1].ravel().astype('float64');v=fr[t2][y0:y1,x0:x1].ravel().astype('float64')
    cc.append(float(np.corrcoef(u,v)[0,1]))
  strip_summary[key]={'n_correlated_pairs':len(cc),'median_pearson_r':statistics.median(cc),'range':[min(cc),max(cc)],
    'geometry':{'y': [y0,y1], 'x':[x0,x1]},'source_phase_only':True,'fixed_rectangles_but_selected_after_y120_source_hypothesis':True}

stat={}
for gl in ARCH:
  stat[gl]={}
  for cl in ['disputed','software_agreement']:
    v=[e for e in source_events if e['gl']==gl and e['software_label']==cl]
    assert v
    stat[gl][cl]={'n':len(v),'small_candidate_probability_mean_median':statistics.median(x['small_MoMA_probtiff_mean_uncalibrated'] for x in v),
      'small_candidate_probability_mean_minmax':[min(x['small_MoMA_probtiff_mean_uncalibrated'] for x in v),max(x['small_MoMA_probtiff_mean_uncalibrated'] for x in v)],
      'preexistent_component_max_IoU_median':statistics.median(x['best_past_5_frames_source_component_iou'] for x in v),
      'preexistent_component_IoU_ge_0_8_count':sum(x['best_past_5_frames_source_component_iou']>=.8 for x in v),
      'software_event_samples_not_independent':True}
assert len(source_events)==25 and len(cases)==37
six=[next(e for e in source_events if e['gl']=='Pos5_GL35' and e['frame']==f) for f in [399,403,405,409,411,415]]
assert all(x['best_past_5_frames_source_component_iou']>.85 for x in six)
results={'study':'901 E515 source-mask cross-time counterfactual + uncalibrated original UNet probability','all_source_zip_receipts':integrity,
  'source_GL':SOURCE_GL,'events':source_events,'grouped_descriptive':stat,'Pos5_six_preexisting_IoU':[x['best_past_5_frames_source_component_iou'] for x in six],
  'Pos5_six_preexisting_means_all_gt_0_85':True,'Pos5_source_phase_fixed_strip_correlations':strip_summary,
  'strength':'source-software candidate preexistence and probability response; does not independently identify biological mothers/daughters',
  'independent_human_annotations_completed':0,'independent_externally_validated_model_completed':False,
  'physical_divisions_truth':'UNRESOLVED','physical_molecular_information_bit':'NOT_IDENTIFIED','publication':'PUBLICATION_NO_GO'}
(OUT/'E515_source_time_probabilities_and_preexistence.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf8')
# flattened table
with (OUT/'E515_source_event_probabilities_and_preexistence.csv').open('w',newline='',encoding='utf8') as h:
  fields=['gl','frame','software_label','small_candidate_area_px','small_candidate_ymin','small_candidate_ymax',
   'small_candidate_xmin','small_candidate_xmax','small_MoMA_probtiff_mean_uncalibrated',
   'small_MoMA_probtiff_median_uncalibrated','best_past_5_frames_source_component_iou','best_past_source_frame',
   'source_small_candidate_id','best_past_source_component_id']
  w=csv.DictWriter(h,fieldnames=fields);w.writeheader();w.writerows({k:v[k] for k in fields} for v in sorted(source_events,key=lambda v:(v['gl'],v['frame'])))

# Hash-ordered, public mode, no event labels or chronological absolute times.
order=sorted(cases,key=lambda c:c['public_id'])
public_manifest={'study':'901 E515 blinded phase-only external annotations','media':'unchanged original 2D preprocessed phase pixels; globally fixed grayscale percentile per GL',
 'not_original_camera_full_fov':True,'no_source_software_mask_or_division_label_in_media':True,
 'offsets_real_TIFF_zero_based_relative_to_unknown_event':OFFSETS,'number_of_cases':37,
 'source_time_stamp_hidden':True,'ordered_cases':[{'case_id':c['public_id'],'images':[f"{c['public_id']}_full_original_PhC.png",f"{c['public_id']}_focus_original_PhC.png"]} for c in order],
 'avoid_treating_consensus_or_plain_time_window_as_verified_ground_truth':True}
(PUBLIC/'E515_public_blind_case_manifest.json').write_text(json.dumps(public_manifest,ensure_ascii=False,indent=2),encoding='utf8')
with (PUBLIC/'E515_blank_independent_reviewer_ratings.csv').open('w',newline='',encoding='utf8') as h:
  fs=['blind_case_id','reviewer_anonymous_id','classification','image_quality_good_fair_poor',
      'distinct_mother_seen_before_t','two_persistent_daughters_visible_by_t_plus_2_or_4','evidence_relative_frame_offsets',
      'evidence_Y_or_full_ROI','uncertainty_or_clipping_reasons','free_text_observation']
  w=csv.DictWriter(h,fieldnames=fs);w.writeheader()
  w.writerows({'blind_case_id':c['public_id']} for c in order)
with (PRIVATE/'E515_private_unblind_key_DO_NOT_SEND_TO_REVIEWER.csv').open('w',newline='',encoding='utf8') as h:
  fs=['blind_case_id','original_GL','source_zero_based_frame','E510_old_case_id','original_software_class','user_recorded_ground_truth']
  w=csv.DictWriter(h,fieldnames=fs);w.writeheader()
  w.writerows({'blind_case_id':c['public_id'],'original_GL':c['gl'],'source_zero_based_frame':c['f'],
      'E510_old_case_id':c['old_case_id'],'original_software_class':c['class'],
      'user_recorded_ground_truth':'UNRESOLVED'} for c in order)

protocol='''# 901 E515 — Phase-contrast independent reviewer protocol (blind media only)

## Scope and limitations

This packet contains **37 pseudonymized image cases** from two growth lanes of one published experiment. The TIFF pixels are original **preprocessed** phase-contrast data, not unprocessed camera frames, and have been checked against source hashes and MoMA CellStats intensities in prior certified rounds. There is **no software lineage, segmentation mask, candidate component, event class, original frame number, or solver metadata** shown in these cards. The full-lane and focus views are exactly the same source float32 pixels with one fixed GL-wise display intensity transformation; the zoom is nearest-neighbor only. The focus crop is identical for every case (y=100..184, x=33..77, zero-based). No borders or colored prediction annotations are placed on the images. Use the full field as primary evidence and the focus crop only to clarify visibility.

The six columns represent original TIFF zero-based times relative to a **concealed central time**, at offsets -2, -1, 0, +1, +2, +4. The frame-0 offset is an observation window anchor, not a certified biological division timestamp. Temporal windows can overlap within the same growth lane. Treat cases as correlated image contexts, **not 37 independent biological samples**.

## Required per-case outcome

Choose one of exactly these four outcome codes (write it in the blank rating CSV):

- `VISIBLE_DIVISION`: A single identifiable physical mother before the anchor separates into two spatially distinct daughter bodies after it, with both credible offspring continuing for the next available frames. Record boundary and time offsets. A bright line alone is not enough.
- `VISIBLE_CONTINUITY_NO_DIVISION`: A clearly seen physical body stays continuous through the full follow-up and its visibility and available frame coverage exclude the proposed split. Do not choose this simply because a gap is hard to see.
- `TIMING_SHIFT`: The same region provides a plausible real physical split outside the candidate anchor; record the relative offset and evidence for continuity.
- `UNRESOLVED`: Ambiguous morphology, overlapping bodies, poor depth of focus, clipping at the channel top, insufficient resolution, or uncertainty about a mother-to-two-daughter lineage. This is a valid and expected result.

Record image quality, whether the mother is unambiguously identifiable before time 0, whether two offspring can be followed to +2 or +4, approximate image y region and relevant offsets, and free text evidence. Do not infer an absolute 16S copy count or Shannon bits from a phase image.

## Independence rules

Annotators should be familiar with phase-contrast microscopy of bacteria and **must not be shown** E508–E515 source/software event reports, the private sealed key or raw MoMA labels before scoring. They must record ratings before class identities are disclosed. Disagreements require a second independent read and an explicit adjudication trail. A single algorithm or ChatGPT interpretation is not independent expert two-person ground truth. This study has not obtained independent human labels yet; the CSV intentionally contains blanks.

The observations from two GLs and repeated times are correlated, and all originate from a single 2022-12-20 source. Software-version consensus is *not* positive biological truth, and an ordinary software window is *not* negative biological truth. No sensitivity, specificity, binary accuracy, p-values or calibrated molecule-level information rates should be reported before an independent reference exists.

## Provenance

Reconstruction script: E515 `e515_sources_and_blind_media.py`. E513 verified source `CellStats.frame=f` ↔ TIFF zero-based index `f`. This packet is a source-verified external-review instrument, **not** an external review already conducted.
'''
(PUBLIC/'E515_blind_review_protocol.md').write_text(protocol,encoding='utf8')

# Do NOT zip the key with the reviewer media! Append checksum manifest separately.
manifest={}
for p in sorted(OUT.rglob('*')):
  if p.is_file() and p.name!='E515_output_sha256_manifest.json':
    manifest[str(p.relative_to(OUT))]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(OUT/'E515_output_sha256_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf8')
print('E515_ANALYSIS_DONE',json.dumps({'samples':len(cases),'source_events':len(source_events),
  'group':stat,'six_preexisting_IoU':[x['best_past_5_frames_source_component_iou'] for x in six],
  'blinded_pngs':len(list(PUBLIC.glob('*.png'))),'raw_phase_unique':[x['source_float_2D_frames'] for x in SOURCE_GL],
  'independent_microscopy_truth':'UNRESOLVED'},ensure_ascii=False),flush=True)