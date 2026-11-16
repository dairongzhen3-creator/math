#!/usr/bin/env python3
"""901 E507 secondary source experiment: compare official main versus cost_function_analysis_1
CellStats by geometry and frame, never same numeric cell_ID across pipelines.
Archives and candidate GLs frozen from E502 authenticated ZIP central directories.
"""
import argparse, pathlib,requests,struct,hashlib,zlib,io,csv,json,time,collections
P=argparse.ArgumentParser();P.add_argument('--archive',required=True);ARG=P.parse_args()
ARCH={
'20221220_E08_GFP.zip':('352c0f14d8cc41f64942c3b850cde4309971c9654e814a722440f0b7c2ce3dbc',78458,['Pos28_GL37']),
'20221220_E08_GFP_pl.zip':('5bd93ce604256d00233d4c5bf91865205b44b0d7fc2f52fefe0f68ef0e23e94c',141104,['Pos12_GL22']),
'20221220_rplN_GFP_pl.zip':('7f55aaa61d7025260773c33aa662e13ae2b2266dea300384f9e09a1667dc5ac6',193307,['Pos4_GL4','Pos5_GL35']),
}
assert ARG.archive in ARCH
exp_sha,exp_n,gls=ARCH[ARG.archive];OUT=pathlib.Path('e512_verified_original_sources')/ARG.archive[:-4];OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session();S.headers['User-Agent']='901-E507-alt-cost-source-method-sensitivity-20261116'
BASE='https://www.ebi.ac.uk/biostudies/files/S-BIAD2506/'
a=S.get(BASE+'Statphase_Transition_GLs/deepMoMA.json',timeout=(15,55));a.raise_for_status()
assert hashlib.sha256(a.content).hexdigest()=='f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f'
m={pathlib.PurePosixPath(i['path']).name:i for i in json.loads(a.content)};N=int(m[ARG.archive]['size']);URL=BASE+m[ARG.archive]['path']
def get(a,b):
 assert 0<=a<=b<N
 out=bytearray()
 for begin in range(a,b+1,4194304):
  end=min(b,begin+4194303)
  for attempt in range(3):
   try:
    r=S.get(URL,headers={'Range':f'bytes={begin}-{end}'},timeout=(15,80),stream=True)
    if r.status_code!=206 or r.headers.get('Content-Range')!=f'bytes {begin}-{end}/{N}':raise RuntimeError('HTTP strict 206 failed')
    x=bytearray()
    for c in r.iter_content(131072):
     if c:x.extend(c)
     if len(x)>end-begin+1:raise ValueError('overflow')
    if len(x)!=end-begin+1:raise ValueError('truncated range')
    out.extend(x);break
   except Exception as e:
    fail=e;time.sleep(1+attempt)
  else:raise fail
 return bytes(out)
def wide(ex,raw,com,off):
 if all(i!=0xffffffff for i in (raw,com,off)):return raw,com,off
 i=0
 while i<len(ex):
  tag,n=struct.unpack_from('<HH',ex,i);d=ex[i+4:i+4+n];i+=4+n
  if tag==1:
   j=0
   if raw==0xffffffff:raw=struct.unpack_from('<Q',d,j)[0];j+=8
   if com==0xffffffff:com=struct.unpack_from('<Q',d,j)[0];j+=8
   if off==0xffffffff:off=struct.unpack_from('<Q',d,j)[0];j+=8
   return raw,com,off
 raise ValueError('ZIP64 missing')
tail=get(N-131072,N-1);i=tail.rfind(b'PK\x05\x06');assert i>=0
x=struct.unpack_from('<4s4H2IH',tail,i);count,length,offset=x[4],x[5],x[6]
if count==65535 or length==0xffffffff or offset==0xffffffff:
 loc=tail.rfind(b'PK\x06\x07',0,i);assert loc>=0
 off=struct.unpack_from('<Q',tail,loc+8)[0];s=get(off,off+55);assert s[:4]==b'PK\x06\x06'
 count=struct.unpack_from('<Q',s,32)[0];length=struct.unpack_from('<Q',s,40)[0];offset=struct.unpack_from('<Q',s,48)[0]
assert count==exp_n and length<67108864
cd=get(offset,offset+length-1);assert hashlib.sha256(cd).hexdigest()==exp_sha

# Full byte authentication of all original small component/optimizer state sources,
# and complete original probability TIFF. The full TIFF is not duplicated into Drive.
import numpy as np
import tifffile
GLS={"20221220_E08_GFP_pl.zip":["Pos12_GL22"],"20221220_rplN_GFP_pl.zip":["Pos5_GL35"]}
FREEZE={
 'Pos5_GL35':{'disputed':[247,294,399,403,405,409,411,415,622],
             'software_agreement':[203,225,232,244]},
 'Pos12_GL22':{'disputed':[241,246,608],
             'software_agreement':[201,203,208,223,226,227,242,247,260]},
}
assert ARG.archive in GLS
needed=['probability_maps.tif','component_forests.json','curation.moma','file_format.json','gurobi_model.sol','gurobi_model.mst']
cost_extra=['hypothesis_states.csv','assignment_states.csv','assignment_costs.csv','assignment_filter_intensities.csv']
sources=collections.defaultdict(list);j=0;n=0
while j<len(cd):
 h=struct.unpack_from('<4s6H3I5H2I',cd,j)
 assert h[0]==bytes.fromhex('504b0102')
 nl,el,cl=h[10:13];name=cd[j+46:j+46+nl].decode('utf-8' if h[3]&2048 else 'cp437','replace')
 for gl in GLS[ARG.archive]:
  if f'/{gl}/' not in name or '__BKP' in name:continue
  bname=pathlib.PurePosixPath(name).name
  if bname not in needed+cost_extra or '/track_data__' not in name:continue
  if '/20221220_analysis/track_data__20221220_analysis/' in name:v='main'
  elif '/cost_function_analysis_1/track_data__cost_function_analysis_1/' in name:v='cost'
  else:continue
  raw,compressed,loc=wide(cd[j+46+nl:j+46+nl+el],h[9],h[8],h[16])
  q={'gl':gl,'version':v,'kind':bname,'path':name,'raw':raw,'compressed':compressed,'loc':loc,'crc':h[7],'method':h[4]}
  sources[(gl,v,bname)].append(q)
 j+=46+nl+el+cl;n+=1
assert j==length and n==count
def materialize(e):
 assert 0<e['raw']<200000000 and 0<e['compressed']<60000000
 loc=e['loc'];lh=struct.unpack('<4s5H3I2H',get(loc,loc+29))
 assert lh[0]==bytes.fromhex('504b0304')
 nl,el=lh[-2:]
 fn=get(loc+30,loc+30+nl+el-1)[:nl].decode('utf-8' if lh[2]&2048 else 'cp437','replace')
 assert fn==e['path']
 off=loc+30+nl+el
 payload=get(off,off+e['compressed']-1)
 assert len(payload)==e['compressed']
 original=zlib.decompress(payload,-15) if e['method']==8 else payload if e['method']==0 else None
 assert original is not None and len(original)==e['raw']
 assert zlib.crc32(original)==e['crc'],('CRC_MISMATCH',e['path'])
 return original

for gl in GLS[ARG.archive]:
 gd=OUT/gl;gd.mkdir(exist_ok=True)
 source_rows=[];maps={};meta={};forests={}
 selected=FREEZE[gl]
 frozen_events=sorted(set(selected['disputed']+selected['software_agreement']))
 frozen_frames=sorted(set([t+off for t in frozen_events for off in [-3,-2,-1,0,1,2,3,5,8] if 1<=t+off<=1036]))
 for ver in ['main','cost']:
  for kind in needed+ (cost_extra if ver=='cost' else []):
   rr=sources[(gl,ver,kind)]
   if kind in cost_extra and not rr:
    continue
   assert len(rr)==1,(gl,ver,kind,len(rr))
   e=rr[0]; b=materialize(e); digest=hashlib.sha256(b).hexdigest()
   rec={'gl':gl,'version':ver,'source_path':e['path'],'kind':kind,
        'source_zip_member_CRC32':f'{e["crc"]:08x}',
        'source_zip_local_offset':e['loc'],'uncompressed_bytes':len(b),
        'original_sha256':digest,'CRC_AND_SHA_COMPUTED':True,
        'compressed_bytes':e['compressed']}
   if kind=='probability_maps.tif':
    with tifffile.TiffFile(io.BytesIO(b)) as tt:
     ser=tt.series[0]
     info={'series_axes':ser.axes,'series_shape':list(ser.shape),'series_dtype':str(ser.dtype),
           'tiff_page_count':len(tt.pages),'source_file_byte_sha256':digest}
     arr=ser.asarray()
     assert arr.dtype.kind in 'fiu'
     if arr.ndim==3 and arr.shape[0]==1036 and arr.shape[-1]==106:
      pass
     elif arr.ndim==4 and arr.shape[0]==1036 and arr.shape[-1]==106 and 1 in arr.shape:
      arr=arr.squeeze()
     else:raise AssertionError(('MAP_TIME_DIM_OR_SHAPE_UNKNOWN',gl,ver,arr.shape,ser.axes))
     assert np.isfinite(arr).all()
     indices=np.asarray([t-1 for t in frozen_frames],dtype='int32')
     selected_array=np.ascontiguousarray(arr[indices])
     np.savez_compressed(gd/f'{ver}__frozen_probability_frames.npz',
       original_frame_numbers=np.asarray(frozen_frames,dtype='int32'),
       probability_values=selected_array)
     info.update({'frozen_unique_frames':len(frozen_frames),
        'selected_numeric_dtype':str(selected_array.dtype),
        'selected_shape':list(selected_array.shape),
        'array_min':float(arr.min()),'array_max':float(arr.max()),
        'array_mean':float(arr.mean()),
        'selected_array_sha256':hashlib.sha256(selected_array.tobytes()).hexdigest(),
        'selection_is_frozen_prior_to_pixel_read':True})
     meta[ver]=info
     print('E512_AUTHENTIC_PROBABILITY_TIFF',gl,ver,len(b),digest,ser.axes,arr.shape,len(frozen_frames),flush=True)
    del arr,selected_array,b
   else:
    save=gd/f'{ver}__{kind}'
    save.write_bytes(b)
    if kind=='component_forests.json':
     try:
      obj=json.loads(b)
      if isinstance(obj,dict):
       desc={'container':'dict','keys':list(obj)[:25],'key_count':len(obj),
         'member_types':{str(k):type(v).__name__ for k,v in list(obj.items())[:15]}}
      elif isinstance(obj,list):
       desc={'container':'list','length':len(obj),'first_type':type(obj[0]).__name__ if obj else None,
         'first_keys':list(obj[0])[:25] if obj and isinstance(obj[0],dict) else []}
      else:desc={'container':type(obj).__name__}
      forests[ver]=desc
      print('E512_AUTHENTIC_COMPONENT_FOREST',gl,ver,len(b),digest,json.dumps(desc,ensure_ascii=False),flush=True)
     except (ValueError,UnicodeError) as e:
      forests[ver]={'json_parse_failed':str(e)}
    elif kind=='curation.moma':
     rec['utf8_first_500']=b[:500].decode('utf8','replace')
    elif kind.startswith('assignment') or kind=='hypothesis_states.csv':
     rec['csv_header']=b.decode('utf8','replace').splitlines()[0][:1000] if b else ''
    elif kind in ('gurobi_model.sol','gurobi_model.mst'):
     rec['first_10_lines']=b.decode('utf8','replace').splitlines()[:10]
   source_rows.append(rec)
  assert 'main' in meta and 'cost' in meta
 print('E512_SOURCE_PAIR_COMPARISON',json.dumps({
   'gl':gl,'probability_map_exact_sha_match':meta['main']['source_file_byte_sha256']==meta['cost']['source_file_byte_sha256'],
   'probability_map_exact_selected_array_sha_match':meta['main']['selected_array_sha256']==meta['cost']['selected_array_sha256'],
   'component_forests_sha_match':next(x['original_sha256'] for x in source_rows if x['version']=='main' and x['kind']=='component_forests.json')==next(x['original_sha256'] for x in source_rows if x['version']=='cost' and x['kind']=='component_forests.json'),
   'curation_sha_match':next(x['original_sha256'] for x in source_rows if x['version']=='main' and x['kind']=='curation.moma')==next(x['original_sha256'] for x in source_rows if x['version']=='cost' and x['kind']=='curation.moma')
 },ensure_ascii=False),flush=True)
 receipt={'gl':gl,'archive':ARG.archive,
  'official_source_manifest_sha256':'f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f',
  'central_directory_sha256':exp_sha,'central_directory_entries':n,
  'frozen_event_source':'E509 authenticated true MoMA software source; not independently validated physical divisions',
  'frozen_events':selected,'frozen_frame_numbers':frozen_frames,'map_properties':meta,
  'component_forest_structure':forests,'complete_source_authentication':source_rows,
  'observational_inference_scope':'source pipeline data only; source model segmentation and software decisions are not independent physical split ground truth',
  'independent_molecular_bits_identified':False}
 (gd/'E512_original_member_SHA_CRC_critical_source_audit.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False),encoding='utf8')
print('E512_SOURCE_AUTH_COMPLETE',ARG.archive,len(sources),flush=True)
