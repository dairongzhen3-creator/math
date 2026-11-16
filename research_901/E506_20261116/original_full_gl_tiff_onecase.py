#!/usr/bin/env python3
"""E506: bounded single originally archived preprocessed GL phase-contrast TIFF.
Strict HTTP 206 ranges, original ZIP64 central-directory and member CRC,
source authenticity before reading selected frame pixels. No full archive.
"""
import pathlib,requests,struct,hashlib,zlib,json,tempfile,os,time,io
import tifffile,numpy as np
from PIL import Image
O=pathlib.Path('e506_one_full_preprocessed_GL');O.mkdir(exist_ok=True)
BASE='https://www.ebi.ac.uk/biostudies/files/S-BIAD2506/'
ARCH='20221220_NO_LIGHT_E08_GFP.zip';GL='Pos25_GL36';FR=[250,251,252,253,254,255]
EXPECTED_CD='d9a00f4ba085b3445448eec94a92c290e6edc66f9b2a1905255879460c4348ee'
S=requests.Session();S.headers['User-Agent']='901-E506-one-GL-targeted-camera-derivative-source-audit-20261116'
man=S.get(BASE+'Statphase_Transition_GLs/deepMoMA.json',timeout=(15,55));man.raise_for_status()
assert hashlib.sha256(man.content).hexdigest()=='f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f'
entries={pathlib.PurePosixPath(x['path']).name:x for x in json.loads(man.content)};source=entries[ARCH]
N=int(source['size']);url=BASE+source['path']
def read_part(start,end):
 assert 0<=start<=end<N and end-start+1<=4194304
 for attempt in range(4):
  try:
   r=S.get(url,headers={'Range':f'bytes={start}-{end}'},timeout=(20,95),stream=True)
   if r.status_code!=206 or r.headers.get('Content-Range')!=f'bytes {start}-{end}/{N}':raise ValueError('BAD_206_OR_CONTENT_RANGE')
   bb=bytearray()
   for c in r.iter_content(131072):
    if c:bb.extend(c)
    if len(bb)>end-start+1:raise ValueError('RANGE_OVERFLOW')
   assert len(bb)==end-start+1
   return bytes(bb)
  except Exception as e:
   ex=e;time.sleep(0.5*(attempt+1))
 raise ex
def zip64_extra(ex,raw,comp,offset):
 if all(v!=0xffffffff for v in (raw,comp,offset)):return raw,comp,offset
 j=0
 while j<len(ex):
  tag,n=struct.unpack_from('<HH',ex,j);data=ex[j+4:j+4+n];j+=4+n
  if tag!=1:continue
  k=0
  if raw==0xffffffff:raw=struct.unpack_from('<Q',data,k)[0];k+=8
  if comp==0xffffffff:comp=struct.unpack_from('<Q',data,k)[0];k+=8
  if offset==0xffffffff:offset=struct.unpack_from('<Q',data,k)[0];k+=8
  return raw,comp,offset
 raise ValueError('ZIP64_EXTRA_NOT_FOUND')
tail=read_part(N-131072,N-1);i=tail.rfind(b'PK\x05\x06');assert i>=0
end=struct.unpack_from('<4s4H2IH',tail,i);cnt,sz,at=end[4],end[5],end[6]
if cnt==65535 or sz==0xffffffff or at==0xffffffff:
 loc=tail.rfind(b'PK\x06\x07',0,i);assert loc>=0
 off=struct.unpack_from('<Q',tail,loc+8)[0];z=read_part(off,off+55);assert z[:4]==b'PK\x06\x06'
 cnt=struct.unpack_from('<Q',z,32)[0];sz=struct.unpack_from('<Q',z,40)[0];at=struct.unpack_from('<Q',z,48)[0]
assert cnt==38345 and 0<sz<50000000
cd=bytearray()
for start in range(at,at+sz,4194304):cd.extend(read_part(start,min(at+sz-1,start+4194303)))
assert hashlib.sha256(cd).hexdigest()==EXPECTED_CD
wanted=f'20221220_NO_LIGHT_E08_GFP/{GL}/20221220_glu_spcm_1_MMStack_{GL}.tif'
found=[];i=0;n=0
while i<len(cd):
 h=struct.unpack_from('<4s6H3I5H2I',cd,i);assert h[0]==b'PK\x01\x02'
 nl,el,cl=h[10:13];name=cd[i+46:i+46+nl].decode('utf-8' if h[3]&2048 else 'cp437','replace')
 if name==wanted:
  raw,com,off=zip64_extra(cd[i+46+nl:i+46+nl+el],h[9],h[8],h[16]);found.append({'raw':raw,'compressed':com,'off':off,'method':h[4],'crc':h[7]})
 i+=46+nl+el+cl;n+=1
assert n==cnt and i==sz and len(found)==1
f=found[0];assert f['method']==8 and f['raw']==395733178 and f['compressed']==321736153
at=f['off'];hdr=struct.unpack('<4s5H3I2H',read_part(at,at+29));assert hdr[0]==b'PK\x03\x04'
nl,el=hdr[-2:];name=read_part(at+30,at+30+nl+el-1)[:nl].decode('utf-8' if hdr[2]&2048 else 'cp437','replace');assert name==wanted
start=at+30+nl+el;stop=start+f['compressed']-1
rawlen=0;sha=hashlib.sha256();crc=0;compress_sha=hashlib.sha256();unpacker=zlib.decompressobj(-15)
outpath=pathlib.Path('/tmp/901_e506_true_single_source_gl.tif')
with outpath.open('wb') as dest:
 for partstart in range(start,stop+1,4194304):
  b=read_part(partstart,min(stop,partstart+4194303));compress_sha.update(b)
  data=unpacker.decompress(b)
  if data:dest.write(data);rawlen+=len(data);sha.update(data);crc=zlib.crc32(data,crc)
 data=unpacker.flush()
 if data:dest.write(data);rawlen+=len(data);sha.update(data);crc=zlib.crc32(data,crc)
assert unpacker.eof and rawlen==f['raw'] and crc==f['crc']
print('E506_TRUE_PREPROCESSED_TIFF_CRC_PASS',rawlen,sha.hexdigest(),flush=True)
with tifffile.TiffFile(outpath) as tt:
 meta={'pages':len(tt.pages),'series':[{'axes':s.axes,'shape':s.shape,'dtype':str(s.dtype)} for s in tt.series],
       'page0_shape':tt.pages[0].shape,'page0_dtype':str(tt.pages[0].dtype),'tiff_is_imagej':bool(tt.is_imagej),
       'tiff_imagej_metadata':str(tt.imagej_metadata)[:1500] if tt.is_imagej else None}
 print('E506_TIFF_HEADER',json.dumps(meta,ensure_ascii=False),flush=True)
 axes=tt.series[0].axes;shape=tt.series[0].shape
 assert any(x==1036 for x in shape),('SOURCE_NOT_1036_FRAMES',axes,shape)
 arr=tt.series[0].asarray() # <= 450 MiB documented original source memory
 assert arr.nbytes<600000000
 t_axis=axes.find('T') if 'T' in axes else next(i for i,d in enumerate(shape) if d==1036)
 arr=np.moveaxis(arr,t_axis,0)
 channel_axis=axes.find('C')
 if channel_axis>=0:
  moved_c=channel_axis+(1 if channel_axis<t_axis else 0) # index after T move to axis 0
  # safer source axes list mutation
  q=list(axes);v=q.pop(t_axis);q.insert(0,v)
  new_c=q.index('C');arr=np.moveaxis(arr,new_c,1)
 else:
  if arr.ndim==3:arr=arr[:,None,:,:]
 assert arr.shape[0]==1036 and arr.shape[-2]>=250 and arr.shape[-1]>=20,(arr.shape,axes)
 print('E506_TIFF_PAGES_AXES',arr.shape,flush=True)
 outviews=[]
 for fr in FR:
  p=arr[fr-1,0]
  assert p.ndim==2
  # save exact numeric pixel data as lossless compressed .npz, only selected frames
  outviews.append(p.copy())
  q=p.astype('float64')
  lo,hi=np.quantile(q,[.01,.99]); hi=max(hi,lo+1e-12)
  png=np.clip((q-lo)*255/(hi-lo),0,255).astype('uint8')
  Image.fromarray(png).save(O/f'Pos25_GL36_GM_preprocessed_PhC_frame_{fr:04d}.png')
 np.savez_compressed(O/'E506_Pos25_GL36_authentic_selected_PhC_numeric_frames.npz',frames=np.array(FR),pixel_values=np.stack(outviews))
 meta['source_selection_shape']=list(arr.shape);meta['selected_frame_numbers']=FR
meta.update({'source_archive':ARCH,'gl':GL,'tiff_zip_member':wanted,'compression_method':8,
 'compressed_member_bytes':f['compressed'],'original_raw_bytes':f['raw'],'original_source_zip_crc32':format(f['crc'],'08x'),
 'full_member_original_crc_pass':True,'original_member_sha256':sha.hexdigest(),'compressed_member_sha256':compress_sha.hexdigest(),
 'source':'BIA S-BIAD2506, preprocessed GL TIFF, not original camera FOV or unprocessed hardware recording',
 'independent_raw_camera_validation':False,'intended_evidence':'Manually inspect real phase contrast pixels, compare with MoMA segmentation and lineage on same preprocessed source',
 'timestamp_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())})
(O/'E506_original_preprocessed_TIFF_receipt.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2))
print('E506_EXTRACTED_SIX_PHASE_FRAMES',json.dumps({'shape':meta['source_selection_shape'],'files':len(FR),'tiff_sha':meta['original_member_sha256']}),flush=True)
try:outpath.unlink()
except:pass
