import argparse, pathlib, requests, hashlib, struct, zlib, json, csv, io, time
p=argparse.ArgumentParser();p.add_argument('--archive',required=True);a=p.parse_args()
CASES={
 '20221220_E08_GFP_pl.zip':('Pos10_GL16',242,243,'5bd93ce604256d00233d4c5bf91865205b44b0d7fc2f52fefe0f68ef0e23e94c',141104),
 '20221220_NO_LIGHT_E08_GFP_pl.zip':('Pos16_GL14',258,259,'bae990f0288f2898ab3c9525181cbceb3637cf2a08f1a120ec5292fad597f7fb',16967),
 '20221220_NO_LIGHT_E08_GFP.zip':('Pos25_GL36',252,253,'d9a00f4ba085b3445448eec94a92c290e6edc66f9b2a1905255879460c4348ee',38345)
}
assert a.archive in CASES
gl,f0,f1,want_cd,want_count=CASES[a.archive]
O=pathlib.Path('e506_bounded')/a.archive[:-4];O.mkdir(parents=True,exist_ok=True)
S=requests.Session();S.headers['User-Agent']='901-E506-original-image-source-bounded-range-20261116'
base='https://www.ebi.ac.uk/biostudies/files/S-BIAD2506/'
man=S.get(base+'Statphase_Transition_GLs/deepMoMA.json',timeout=(15,45));man.raise_for_status()
assert hashlib.sha256(man.content).hexdigest()=='f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f'
source={pathlib.PurePosixPath(v['path']).name:v for v in json.loads(man.content)}
t=source[a.archive];N=int(t['size']);url=base+t['path']
def get(start,end):
 assert 0<=start<=end<N and end-start+1<=4194304
 for it in range(3):
  try:
   r=S.get(url,headers={'Range':f'bytes={start}-{end}'},timeout=(15,65),stream=True)
   if r.status_code!=206 or r.headers.get('Content-Range')!=f'bytes {start}-{end}/{N}':raise RuntimeError('HTTP206_OR_RANGE_FAILED')
   b=bytearray()
   for c in r.iter_content(65536):
    if c:b+=c
    if len(b)>end-start+1:raise ValueError('OVERSIZE')
   assert len(b)==end-start+1
   return bytes(b)
  except Exception as exc:
   last=exc;time.sleep(0.7*(it+1))
 raise last
def get_big(start,end):
 return b''.join(get(i,min(i+4194303,end)) for i in range(start,end+1,4194304))
def z64(extra,raw,compressed,offset):
 if all(x!=0xffffffff for x in (raw,compressed,offset)):return raw,compressed,offset
 i=0
 while i<len(extra):
  kind,n=struct.unpack_from('<HH',extra,i);data=extra[i+4:i+4+n];i+=4+n
  if kind!=1:continue
  j=0
  if raw==0xffffffff:raw=struct.unpack_from('<Q',data,j)[0];j+=8
  if compressed==0xffffffff:compressed=struct.unpack_from('<Q',data,j)[0];j+=8
  if offset==0xffffffff:offset=struct.unpack_from('<Q',data,j)[0];j+=8
  return raw,compressed,offset
 raise ValueError('ZIP64_MISSING')
tail=get(N-131072,N-1);i=tail.rfind(b'PK\x05\x06');assert i>=0
h=struct.unpack_from('<4s4H2IH',tail,i);cnt,sz,at=h[4],h[5],h[6]
if cnt==65535 or sz==0xffffffff or at==0xffffffff:
 loc=tail.rfind(b'PK\x06\x07',0,i);assert loc>=0
 off=struct.unpack_from('<Q',tail,loc+8)[0];zz=get(off,off+55)
 assert zz[:4]==b'PK\x06\x06'
 cnt=struct.unpack_from('<Q',zz,32)[0];sz=struct.unpack_from('<Q',zz,40)[0];at=struct.unpack_from('<Q',zz,48)[0]
assert cnt==want_count and 0<sz<=67108864
cd=get_big(at,at+sz-1)
assert len(cd)==sz and hashlib.sha256(cd).hexdigest()==want_cd
targetframes=set(range(f0-2,f1+3));entries=[];i=0;n=0
while i<len(cd):
 h=struct.unpack_from('<4s6H3I5H2I',cd,i);assert h[0]==b'PK\x01\x02'
 nl,el,cl=h[10:13];name=cd[i+46:i+46+nl].decode('utf8' if h[3]&2048 else 'cp437','replace')
 ex=cd[i+46+nl:i+46+nl+el];s=pathlib.PurePosixPath(name).name
 gl_hit=f'/{gl}/' in name and '__BKP' not in name
 selected=gl_hit and ((s.startswith('glf_') and s.endswith('.png') and s[4:-4].isdigit() and int(s[4:-4]) in targetframes) or s.startswith('CellStats__') and s.endswith('.csv') or s.endswith('_kymo.tif'))
 tif=gl_hit and (s.endswith('.tif') and ('MMStack' in s or 'CellMasks' in s))
 if selected or tif:
  raw,comp,off=z64(ex,h[9],h[8],h[16])
  entries.append({'path':name,'raw_bytes':raw,'compressed_bytes':comp,'local_offset':off,'method':h[4],'crc':h[7],'selected':bool(selected),'tiff':bool(tif)})
 i+=46+nl+el+cl;n+=1
assert n==cnt and i==sz
png=[e for e in entries if e['selected'] and e['path'].endswith('.png')]
stat=[e for e in entries if e['selected'] and e['path'].endswith('.csv')]
kymos=[e for e in entries if e['selected'] and e['path'].endswith('_kymo.tif')]
assert len(png)==len(targetframes) and len(stat)==1 and len(kymos)==1,(len(png),len(stat),len(kymos))
rec={'project':901,'round':'E506','GL':gl,'archive':a.archive,'archive_bytes':N,'source_cdir_sha':want_cd,'source_entries':n,'expected_parent_last_frame':f0,'expected_daughter_first_frame':f1,'source_tiff_metadata':[e for e in entries if e['tiff']],'raw_png_sources':[],'raw_cellstats':None,'original_gl_kymograph':None,'note':'glf PNGs are MoMA processed diagnostic frames; raw TIFF provenance assessed separately; not independent original camera data'}
for e in sorted(png+stat+kymos,key=lambda x:x['path']):
 assert e['raw_bytes']<15000000 and e['compressed_bytes']<4194304
 off=e['local_offset'];lh=struct.unpack('<4s5H3I2H',get(off,off+29));assert lh[0]==b'PK\x03\x04'
 nl,el=lh[-2:];head=get(off+30,off+30+nl+el-1)
 assert head[:nl].decode('utf8' if lh[2]&2048 else 'cp437','replace')==e['path']
 start=off+30+nl+el;compressed=get(start,start+e['compressed_bytes']-1)
 if e['method']==8:raw=zlib.decompress(compressed,-15)
 elif e['method']==0:raw=compressed
 else:raise ValueError('UNSUPPORTED_METHOD')
 assert len(raw)==e['raw_bytes'] and zlib.crc32(raw)==e['crc']
 name=pathlib.PurePosixPath(e['path']).name
 (O/name).write_bytes(raw)
 item={'path':e['path'],'sha256':hashlib.sha256(raw).hexdigest(),'crc32':format(e['crc'],'08x'),'bytes':len(raw),'zip_member_offset':off}
 if name.endswith('.png'):
  assert raw[:8]==b'\x89PNG\r\n\x1a\n';rec['raw_png_sources'].append(item)
 elif name.endswith('_kymo.tif'):
  assert raw[:2] in (b'II',b'MM');rec['original_gl_kymograph']=item
 else:
  lines=[l for l in raw.decode('utf-8-sig').splitlines() if l.strip() and not l.startswith('#')]
  item['columns']=list(csv.DictReader(lines).fieldnames)
  rec['raw_cellstats']=item
(O/'E506_source_kymograph_receipt.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2))
print('E506_PREFLIGHT_PASS',json.dumps({'GL':gl,'archive':a.archive,'pngs':len(png),'tiff_members':[(x['path'].split('/')[-1],x['raw_bytes'],x['compressed_bytes'],x['method']) for x in rec['source_tiff_metadata']]},ensure_ascii=False),flush=True)
