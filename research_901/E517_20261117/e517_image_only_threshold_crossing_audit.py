#!/usr/bin/env python3
"""Image-only crossing features; no pretrained model; not biological ground truth.
Independent of source MoMA masks. Otsu base threshold from all original TIFF source phase pixels.
The ROI geometry (y=120) was selected after previous MoMA findings and is not blind.
"""
import zipfile,io,json,hashlib,statistics,collections,argparse
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter,label
from skimage.filters import threshold_otsu
P=argparse.ArgumentParser();P.add_argument('--input-root',default='/mnt/data');P.add_argument('--output',default='/mnt/data/901_E517_20261117')
A=P.parse_args();R=Path(A.input_root);O=Path(A.output);O.mkdir(parents=True,exist_ok=True)
sha={'Pos5_GL35':'745815a73bc1e2eaa01a2cf99f4edd9c530064d10c22bc5421f6ff5189629d8d','Pos12_GL22':'17f80af916f376dfa380b30e8a8b2e8459358a4b66e13e21e23f637aaa08cfa9'}
result={}
for gl,expected in sha.items():
 raw=(R/f'901_E510_{gl}_original_pixels.zip').read_bytes();assert hashlib.sha256(raw).hexdigest()==expected
 with zipfile.ZipFile(io.BytesIO(raw)) as z:
  assert z.testzip() is None
  k=json.loads(next(z.read(n) for n in z.namelist() if n.endswith('E510_sealed_event_key_do_not_open_for_blind_review.json')))
  frames={}
  for name in z.namelist():
   if not name.endswith('.npz'):continue
   with np.load(io.BytesIO(z.read(name)),allow_pickle=False) as q:
    for s,p in zip(q['frames'],q['phase_pixels']):
     t=int(s)-1
     if t in frames:assert np.array_equal(frames[t],p)
     else:frames[t]=p.copy()
 thr=threshold_otsu(np.concatenate([f[95:155,38:74].ravel() for f in frames.values()]))
 raw_rows=[]
 for t,frame in frames.items():
  data=frame[100:143,38:74]
  for factor in (.8,1.,1.2):
   mm=gaussian_filter(data,sigma=.8)>thr*factor
   labs,n=label(mm)
   crossing=[]
   # top cutoff relative crop 100 => y120 is rel20
   for cc in range(1,n+1):
    yy,xx=np.where(labs==cc)
    if len(yy)<12:continue
    if (yy<=18).any() and (yy>=22).any():
     crossing.append({'area_px':len(yy),'miny':int(min(yy)+100),'maxy':int(max(yy)+100),'width':int(max(xx)-min(xx)+1)})
   raw_rows.append({'t':t,'factor':factor,'n_all_components':n,'n_crossing':len(crossing),'crossing':crossing})
 # save image-only results before class merge
 (O/f'E517_{gl}_image_only_roi_crossing_before_label_join.json').write_text(json.dumps(raw_rows),encoding='utf-8')
 cats=collections.defaultdict(dict)
 for c in k:cats[int(c['frame'])+1][c['class']]=True
 by={}
 for c in ['disputed','software_agreement','software_nondivision']:
  by[c]={}
  for factor in (.8,1.0,1.2):
   rr=[r for r in raw_rows if r['factor']==factor and c in cats[r['t']]]
   by[c][str(factor)]={'n':len(rr),'crossing_any':sum(r['n_crossing']>0 for r in rr),'median_count':statistics.median(r['n_crossing'] for r in rr)}
 all_gl_cross={str(factor):{'n':len([r for r in raw_rows if r['factor']==factor]),'n_with_bright_connection_crossing_y120':sum(r['n_crossing']>0 for r in raw_rows if r['factor']==factor)} for factor in (.8,1.0,1.2)}
 result[gl]={'global_OTSU_source_input':float(thr),'all_unique_frames':len(frames),'by_software_posthoc':by,
  'all_GL_frames_cross_by_factor':all_gl_cross,
  'gl_cropped_source_tiff_ROI_y100_143_x38_74':'posthoc engineering check',
  'strong_external_segmentation_truth':'NOT_PROVIDED','model_download_or_omnipose_inference':'NOT_RUN'}
 print('E517_CROSSING',gl,json.dumps(by))
(O/'E517_image_only_bright_connected_crossing_summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')