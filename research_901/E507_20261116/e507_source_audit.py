#!/usr/bin/env python3
"""901 E507: 9 pre-frozen challenge GLs + 3 E506 source controls.
Original BioImage S-BIAD2506 ZIP64 central-directory SHA and exact member CRC.
Pixel-only division surprise index uses no parent box/ID; independent source metadata
is compared only AFTER scoring. Null time points are correlated technical controls,
NOT a biological p-value or a camera-independent image segmentation.
"""
import argparse, pathlib, requests, struct, hashlib, time, json, zlib, csv, io, collections
import numpy as np
import tifffile

A=argparse.ArgumentParser();A.add_argument('--archive',required=True);ARG=A.parse_args()
ARCHIVES={
 '20221220_E08_GFP.zip': ('352c0f14d8cc41f64942c3b850cde4309971c9654e814a722440f0b7c2ce3dbc',78458),
 '20221220_E08_GFP_pl.zip': ('5bd93ce604256d00233d4c5bf91865205b44b0d7fc2f52fefe0f68ef0e23e94c',141104),
 '20221220_rplN_GFP.zip': ('6d91b5e59329712b6920aa1bb78fd06c8f838197be38cf546dbddaace35b1402',65624),
 '20221220_rplN_GFP_pl.zip': ('7f55aaa61d7025260773c33aa662e13ae2b2266dea300384f9e09a1667dc5ac6',193307),
 '20221220_NO_LIGHT_E08_GFP.zip': ('d9a00f4ba085b3445448eec94a92c290e6edc66f9b2a1905255879460c4348ee',38345),
 '20221220_NO_LIGHT_E08_GFP_pl.zip': ('bae990f0288f2898ab3c9525181cbceb3637cf2a08f1a120ec5292fad597f7fb',16967),
 '20221220_NO_LIGHT_rplN_GFP.zip': ('79fbbe50faef4d417813e81162370848e631d9d4975ab1f696daeacc5e65aeeb',27476),
}
# Immutable prospective selection: gl, old parent ID, last parent frame, first daughter frame,
# reason; selected from E505 authenticated ledger prior to image-only score calculation.
CASES={
 '20221220_E08_GFP.zip':[
 ('Pos27_GL8',75,239,240,'new_early_light_survive'),
 ('Pos28_GL37',64,269,270,'new_illum_unknown'),
 ],
 '20221220_E08_GFP_pl.zip':[
 ('Pos10_GL30',44,239,240,'new_early_light_no_growth'),
 ('Pos10_GL16',75,242,243,'E506_control_light')],
 '20221220_rplN_GFP.zip':[
 ('Pos21_GL33',56,378,379,'new_very_late_unknown'),
 ],
 '20221220_rplN_GFP_pl.zip':[
 ('Pos4_GL13',21,544,545,'new_latest_70h_two_daughters'),
 ('Pos4_GL21',60,344,345,'new_late_30h'),
 ],
 '20221220_NO_LIGHT_E08_GFP.zip':[
 ('Pos25_GL36',49,252,253,'E506_control_discordant_sisters'),
 ],
 '20221220_NO_LIGHT_E08_GFP_pl.zip':[
 ('Pos16_GL30',29,239,240,'new_early_dark_survive'),
 ('Pos16_GL16',92,274,275,'new_late_dark_survive'),
 ('Pos16_GL14',89,258,259,'E506_control_dark')],
 '20221220_NO_LIGHT_rplN_GFP.zip':[
 ('Pos24_GL14',90,270,271,'new_mid_dark_survive'),
 ]
}
assert ARG.archive in ARCHIVES and ARG.archive in CASES
OUT=pathlib.Path('e507_source_challenges')/ARG.archive[:-4]
OUT.mkdir(parents=True,exist_ok=True)
expected_cd, expected_count=ARCHIVES[ARG.archive]
BASE='https://www.ebi.ac.uk/biostudies/files/S-BIAD2506/'
s=requests.Session();s.headers['User-Agent']='901-E507-fixed-source-kymo-negative-control-Audit-20261116'
r=s.get(BASE+'Statphase_Transition_GLs/deepMoMA.json',timeout=(15,50));r.raise_for_status()
assert hashlib.sha256(r.content).hexdigest()=='f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f'
m={pathlib.PurePosixPath(x['path']).name:x for x in json.loads(r.content)}
item=m[ARG.archive];N=int(item['size']);url=BASE+item['path']

def one_range(lo,hi):
 assert 0<=lo<=hi<N and hi-lo+1<=4194304
 for it in range(3):
  try:
   rr=s.get(url,headers={'Range':f'bytes={lo}-{hi}'},stream=True,timeout=(15,80))
   if rr.status_code!=206 or rr.headers.get('Content-Range')!=f'bytes {lo}-{hi}/{N}':raise ValueError('not strict HTTP206')
   b=bytearray()
   for chunk in rr.iter_content(131072):
    if chunk:b.extend(chunk)
    if len(b)>hi-lo+1:raise ValueError('range overflow')
   if len(b)!=hi-lo+1:raise ValueError('range length not met')
   return bytes(b)
  except Exception as e:
   err=e; time.sleep(1+it)
 raise err

def range_read(lo,hi):
 return b''.join(one_range(x,min(x+4194303,hi)) for x in range(lo,hi+1,4194304))

def zip64(ex,raw,com,off):
 if all(x!=0xffffffff for x in (raw,com,off)):return raw,com,off
 i=0
 while i<len(ex):
  tag,n=struct.unpack_from('<HH',ex,i);data=ex[i+4:i+4+n];i+=4+n
  if tag==1:
   j=0
   if raw==0xffffffff:raw=struct.unpack_from('<Q',data,j)[0];j+=8
   if com==0xffffffff:com=struct.unpack_from('<Q',data,j)[0];j+=8
   if off==0xffffffff:off=struct.unpack_from('<Q',data,j)[0];j+=8
   return raw,com,off
 raise ValueError('ZIP64 extra unavailable')

tail=one_range(N-131072,N-1);i=tail.rfind(b'PK\x05\x06');assert i>=0
x=struct.unpack_from('<4s4H2IH',tail,i);count,length,offset=x[4],x[5],x[6]
if count==65535 or length==0xffffffff or offset==0xffffffff:
 loc=tail.rfind(b'PK\x06\x07',0,i);assert loc>=0
 z=one_range(struct.unpack_from('<Q',tail,loc+8)[0],struct.unpack_from('<Q',tail,loc+8)[0]+55)
 assert z[:4]==b'PK\x06\x06'
 count=struct.unpack_from('<Q',z,32)[0];length=struct.unpack_from('<Q',z,40)[0];offset=struct.unpack_from('<Q',z,48)[0]
assert count==expected_count and 0<length<67108864
cd=range_read(offset,offset+length-1);assert hashlib.sha256(cd).hexdigest()==expected_cd
found={q[0]:{} for q in CASES[ARG.archive]}; i=0;n=0
while i<len(cd):
 h=struct.unpack_from('<4s6H3I5H2I',cd,i);assert h[0]==b'PK\x01\x02'
 nl,el,cl=h[10:13];name=cd[i+46:i+46+nl].decode('utf8' if h[3]&2048 else 'cp437','replace')
 ex=cd[i+46+nl:i+46+nl+el]
 if '__BKP' not in name:
  for gl,parent,f0,f1,group in CASES[ARG.archive]:
   if f'/{gl}/' not in name:continue
   base=pathlib.PurePosixPath(name).name
   if base.endswith('_kymo.tif'):kind='kymo'
   elif base.startswith('CellStats__') and base.endswith('.csv'):kind='stats'
   elif base.startswith('glf_') and base.endswith('.png') and base[4:-4].isdigit() and f0-2<=int(base[4:-4])<=f1+2:kind=base
   else:continue
   raw,com,off=zip64(ex,h[9],h[8],h[16])
   found[gl].setdefault(kind,[]).append({'name':name,'raw':raw,'compressed':com,'off':off,'crc':h[7],'method':h[4]})
 i+=46+nl+el+cl;n+=1
assert n==count and i==length

def materialize(q):
 assert q['raw']>0 and q['raw']<15000000 and q['compressed']<15000000
 off=q['off'];hd=struct.unpack('<4s5H3I2H',one_range(off,off+29));assert hd[0]==b'PK\x03\x04'
 nl,el=hd[-2:];hdr=one_range(off+30,off+30+nl+el-1)
 assert hdr[:nl].decode('utf-8' if hd[2]&2048 else 'cp437','replace')==q['name']
 data=range_read(off+30+nl+el,off+30+nl+el+q['compressed']-1)
 raw=zlib.decompress(data,-15) if q['method']==8 else data if q['method']==0 else None
 assert raw is not None and len(raw)==q['raw'] and zlib.crc32(raw)==q['crc']
 return raw

results=[];raw_receipts=[]
for gl,parent,f0,f1,group in CASES[ARG.archive]:
 dicts=found[gl];kinds=['stats','kymo']+[f'glf_{t}.png' for t in range(f0-2,f1+3)]
 assert all(len(dicts.get(k,[]))==1 for k in kinds),(gl,[(k,len(dicts.get(k,[]))) for k in kinds])
 input_bytes={}
 for name in kinds:
  q=dicts[name][0];b=materialize(q);input_bytes[name]=b
  raw_receipts.append({'archive':ARG.archive,'GL':gl,'kind':name,'zip_offset':q['off'], 'src_path':q['name'],'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'crc32':f'{q["crc"]:08x}'})
  if name=='stats' or name=='kymo':(OUT/f'{gl}__{name}.csv' if name=='stats' else OUT/f'{gl}__kymo.tif').write_bytes(b)
  # Preserve six original PNGs too, in a nested directory to avoid ambiguous names
  elif name.endswith('.png'):
   d=OUT/gl;d.mkdir(exist_ok=True);(d/name).write_bytes(b)
 # Verify frozen source mother-daughter identity; no labels used to derive the pixel-only score.
 text=input_bytes['stats'].decode('utf-8-sig')
 rows=[t for t in text.splitlines() if t.strip() and not t.startswith('#')]
 parsed=list(csv.DictReader(rows));father=[x for x in parsed if int(x['cell_ID'])==parent and int(x['frame'])==f0]
 children=[x for x in parsed if int(x['frame'])==f1 and int(x['parent_ID'])==parent]
 assert len(father)==1 and father[0]['type_of_end']=='div' and len(children)==2
 fatherY=(int(father[0]['bbox_top_px']),int(father[0]['bbox_bottom_px']))
 childY=sorted([(int(c['bbox_top_px']),int(c['bbox_bottom_px'])) for c in children])
 seam=(childY[0][1]+childY[1][0])/2
 # Image-blind event score: only full GL phase-contrast kymo, no cell IDs, ROIs or masks.
 # Kymo column index is glf frame-1, established and frozen in E506.
 img=tifffile.imread(io.BytesIO(input_bytes['kymo']))
 assert img.ndim==3 and img.shape[0]==3 and img.shape[-1]==1036
 phase=img[0].astype('float64');assert np.isfinite(phase).all()
 phase=phase-np.median(phase,axis=0,keepdims=True)
 diff=np.abs(np.diff(phase,axis=1))
 global_score=np.mean(diff,axis=0)
 width=20
 # maximum 20px spatial window (independent of source bbox location)
 kernel=np.ones(width)/width
 local_score=np.array([np.convolve(diff[:,j],kernel,mode='valid').max() for j in range(diff.shape[1])])
 event_index=f0-1
 candidates=[t for t in range(f0-21,f0+20) if 0<=t<len(global_score) and abs(t-event_index)>2]
 assert len(candidates)==36
 def get_stat(scores):
  e=float(scores[event_index]);null=scores[candidates]
  quantile=(1+int(np.count_nonzero(null<=e)))/(1+len(null))
  return {'event_value':e,'null_median':float(np.median(null)),'within_GL_percentile_rank':float(quantile),'strong_signal_at_95pct':bool(quantile>=.95), 'null_n':len(null)}
 local_at_event=np.convolve(diff[:,event_index],kernel,mode='valid')
 maxloc=float(np.argmax(local_at_event)+(width-1)/2)
 rec={'archive':ARG.archive,'GL':gl,'selection_reason':group,'pre_parent_ID':parent,'pre_parent_frame':f0,
      'first_child_frame':f1,'source_two_children':sorted([int(c['cell_ID']) for c in children]),
      'author_bbox_parent_y':list(fatherY),'author_bbox_daughters_y':childY,'source_seam_y':seam,
      'pixel_only_spatial_window_max_center_y':maxloc,'pixel_hotspot_dist_to_MoMA_seam_y':abs(maxloc-seam),
      'pixel_hotspot_inside_author_parent_bbox':bool(fatherY[0]<=maxloc<=fatherY[1]),
      'phase_kymo_shape':list(img.shape),'full_line_absdiff_stat':get_stat(global_score),
      'best_20px_window_absdiff_stat':get_stat(local_score),'score_details_blind_to_source_lineage':True,
      'source_exposure_type':'NO_LIGHT' if '_NO_LIGHT_' in ARG.archive else 'illuminated',
      'source_is_preprocessed':True,'source_is_raw_camera_FOV':False}
 results.append(rec)
 print('E507_GL_SCORE',json.dumps({'GL':gl,'group':group,'event_frame':f0,'p_global':rec['full_line_absdiff_stat']['within_GL_percentile_rank'],'p_20px':rec['best_20px_window_absdiff_stat']['within_GL_percentile_rank'],'hotspot_error_y':rec['pixel_hotspot_dist_to_MoMA_seam_y']}),flush=True)
with (OUT/'E507_frozen_GL_scores.json').open('w') as f:json.dump(results,f,indent=2,ensure_ascii=False)
with (OUT/'E507_source_authentic_member_receipts.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(raw_receipts[0]));w.writeheader();w.writerows(raw_receipts)
(OUT/'E507_acquisition_identity.json').write_text(json.dumps({'round':'E507','archive':ARG.archive,'source_manifest_sha256':'f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f', 'central_dir_sha256':expected_cd,'archive_member_count':n,'gl_cases':len(results), 'verified_files':len(raw_receipts),'no_full_image_download':True,'time_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())},indent=2))
print('E507_ARCHIVE_COMPLETE',ARG.archive,len(results),len(raw_receipts),flush=True)
