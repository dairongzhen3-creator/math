#!/usr/bin/env python3
"""901 E517 independent threshold-free phase boundary optical diagnostic.

No MoMA masks or software labels are used to derive the optical measurements.
Numerical original float32 microscopy is the only input to score computation.
Software categories are joined after each frame has been scored, for
calibration as technical controls rather than biological ground truth.
Post-hoc ROI x[42:66] and known MoMA ROI edge y=120 are hypothesis-driven,
not a preregistered prospective test.
"""
import io,json,hashlib,zipfile,collections,statistics,csv,argparse
from pathlib import Path
import numpy as np
P=argparse.ArgumentParser();P.add_argument('--input-root',default='/mnt/data');P.add_argument('--output',default='/mnt/data/901_E517_20261117')
A=P.parse_args();R=Path(A.input_root);O=Path(A.output);O.mkdir(parents=True,exist_ok=True)
sha={
'Pos5_GL35':'745815a73bc1e2eaa01a2cf99f4edd9c530064d10c22bc5421f6ff5189629d8d',
'Pos12_GL22':'17f80af916f376dfa380b30e8a8b2e8459358a4b66e13e21e23f637aaa08cfa9'
}
rows=[];summary={}
for gl in sha:
 raw=(R/f'901_E510_{gl}_original_pixels.zip').read_bytes()
 assert hashlib.sha256(raw).hexdigest()==sha[gl]
 f={};cases=[]
 with zipfile.ZipFile(io.BytesIO(raw)) as z:
  assert z.testzip() is None
  cases=json.loads(next(z.read(n) for n in z.namelist() if n.endswith('E510_sealed_event_key_do_not_open_for_blind_review.json')))
  for name in z.namelist():
   if not name.endswith('_exact_float32_phase_frames.npz'):continue
   with np.load(io.BytesIO(z.read(name)),allow_pickle=False) as npz:
    for label,frame in zip(npz['frames'],npz['phase_pixels']):
     t=int(label)-1
     if t in f:assert np.array_equal(f[t],frame)
     else:f[t]=frame.copy()
 assert len(f)==(119 if gl=='Pos5_GL35' else 84)
 events={int(c['frame'])+1:[] for c in cases}
 for c in cases:events[int(c['frame'])+1].append(c['class'])
 blind=[]
 for t,im in sorted(f.items()):
  p={'gl':gl,'frame_t_zero_based':t,'phase_frame_sha256':hashlib.sha256(im.tobytes()).hexdigest()}
  for lo,hi in [(42,66),(30,75)]:
   arr=im[:,lo:hi].astype(np.float64)
   g=np.mean(np.abs(np.diff(arr,axis=0)),axis=1)
   base=float(np.median(g[109:149]));edge=float(g[119])
   p[f'y119_to_y120_edge_mean_abs_x{lo}_{hi}']=edge
   p[f'y119_to_y120_over_local_median_x{lo}_{hi}']=float(edge/base)
   # compare to nearest +/- 5 rows; no ROI model boundary predictions used.
   p[f'y119_to_y120_over_neighbor_median_x{lo}_{hi}']=float(edge/np.median(np.r_[g[114:119],g[120:125]]))
  blind.append(p)
 # Save all case-free optical measurements before mapping any software class.
 (O/f'E517_{gl}_blind_all_frames_threshold_free_optical_metrics.json').write_text(json.dumps(blind,indent=2),encoding='utf8')
 for rec in blind:
  rec['source_software_event_categories_posthoc']=events.get(rec['frame_t_zero_based'],[])
  rows.append(rec)
 bycat={}
 for cls in ['disputed','software_agreement','software_nondivision']:
  subset=[a for a in blind if cls in a['source_software_event_categories_posthoc']]
  bycat[cls]={
   'n_windows':len(subset),
   'median_edge_gradient_local_ratio_x42_66':float(statistics.median(x['y119_to_y120_over_local_median_x42_66'] for x in subset)),
   'median_edge_gradient_neighbor_ratio_x42_66':float(statistics.median(x['y119_to_y120_over_neighbor_median_x42_66'] for x in subset)),
   'median_edge_gradient_local_ratio_x30_75':float(statistics.median(x['y119_to_y120_over_local_median_x30_75'] for x in subset)),
  }
 allstats={k:{'median':float(statistics.median(x[k] for x in blind)),'p10':float(np.quantile([x[k] for x in blind],.10)),'p90':float(np.quantile([x[k] for x in blind],.90))}
     for k in ['y119_to_y120_over_local_median_x42_66','y119_to_y120_over_local_median_x30_75']}
 target=[x for x in blind if gl=='Pos5_GL35' and x['frame_t_zero_based'] in (400,404,406,410,412,416)]
 summary[gl]={'source_sha256':sha[gl],'n_unique_frames':len(blind),'n_fixed_windows':len(cases),'by_technical_software_class':bycat,
              'whole_GL_frame_distribution':allstats,
              'six_pos5_focus':[{k:v for k,v in x.items() if 'x42_66' in k or k=='frame_t_zero_based'} for x in target],
              'technical_interpretation':'phase optical row gradients at y120, not independent physical cell division truth',
              'unblinding_after_metric_calculation':True,'x_subregions_chosen_posthoc':True}
 print('E517_OPTICAL_ROW_RESULT',gl,json.dumps(summary[gl],ensure_ascii=False))
with (O/'E517_all_203_phc_threshold_free_y120_crossing_measurements.csv').open('w',newline='',encoding='utf8') as f:
 fieldnames=list(rows[0]);w=csv.DictWriter(f,fieldnames=fieldnames);w.writeheader();w.writerows(rows)
(O/'E517_optical_boundary_technical_group_descriptives.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf8')
print('E517_OPTICAL_COMPLETED',len(rows))