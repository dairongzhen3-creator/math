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
exp_sha,exp_n,gls=ARCH[ARG.archive];OUT=pathlib.Path('e511_source_sidecars')/ARG.archive[:-4];OUT.mkdir(parents=True,exist_ok=True)
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

# Read only original small per-GL settings / diagnostics from both official
# analysis paths. No central-source duplicate __BKP, no image downloads.
import re
GLS={"20221220_E08_GFP_pl.zip":["Pos12_GL22"],"20221220_rplN_GFP_pl.zip":["Pos5_GL35"]}
assert ARG.archive in GLS
gls=GLS[ARG.archive]
candidates=collections.defaultdict(list);i=0;n=0
all_gl_names=collections.defaultdict(list)
while i<len(cd):
 h=struct.unpack_from('<4s6H3I5H2I',cd,i);assert h[0]==bytes.fromhex('504b0102')
 nl,el,cl=h[10:13]
 path=cd[i+46:i+46+nl].decode('utf-8' if h[3]&2048 else 'cp437','replace')
 for gl in gls:
  if f'/{gl}/' not in path or '__BKP' in path:continue
  all_gl_names[gl].append(path)
  base=pathlib.PurePosixPath(path).name.lower()
  if base not in ('mm.properties','analysis_metadata.json','moma.log') and not (
     base.endswith('.properties') or (base.endswith('.json') and any(t in base for t in ('config','moma','analysis','meta','param')))):
   continue
  if '/20221220_analysis/' in path:version='main'
  elif '/cost_function_analysis_1/' in path:version='cost_function'
  else:version='other'
  raw,compressed,loc=wide(cd[i+46+nl:i+46+nl+el],h[9],h[8],h[16])
  if not (raw>0 and raw<7500000 and compressed>0 and compressed<7500000):
   continue
  candidates[gl].append({'path':path,'version':version,'raw':raw,'compressed':compressed,'off':loc,'method':h[4],'crc':h[7]})
 i+=46+nl+el+cl;n+=1
assert i==length and n==count

def decode(e):
 off=e['off'];h=struct.unpack('<4s5H3I2H',get(off,off+29));assert h[0]==bytes.fromhex('504b0304')
 nl,el=h[-2:]
 name=get(off+30,off+30+nl+el-1)[:nl]
 assert name.decode('utf-8' if h[2]&2048 else 'cp437','replace')==e['path']
 start=off+30+nl+el
 payload=get(start,start+e['compressed']-1)
 raw=zlib.decompress(payload,-15) if e['method']==8 else payload if e['method']==0 else None
 assert raw is not None and len(raw)==e['raw'] and zlib.crc32(raw)==e['crc']
 return raw

records=[]
for gl in gls:
 entries=candidates[gl]
 assert entries and len(entries)<70,(gl,len(entries))
 (OUT/gl).mkdir(exist_ok=True)
 listall=OUT/(gl+'_entire_GL_official_members_listing.json')
 listall.write_text(json.dumps({'gl':gl,'original_zip_member_count':len(all_gl_names[gl]),
   'filename_list':all_gl_names[gl]},indent=2,ensure_ascii=False),encoding='utf8')
 for q in entries:
  b=decode(q)
  basename=pathlib.PurePosixPath(q['path']).name
  version=q['version']
  dest=OUT/gl/(version+'__'+str(len(records))+'__'+basename)
  dest.write_bytes(b)
  content=b.decode('utf-8-sig','replace')
  lines=content.splitlines()
  whole=re.compile(r'(?i)(roi|crop|segment|lane|cell|top|bottom|boundary|offset|gl_width|time_range|120|model|optimizer)')
  exact120=re.compile(r'(?<![0-9])120(?![0-9])')
  hits=[{'line':j+1,'text':line[:280]} for j,line in enumerate(lines) if whole.search(line)]
  coords=[{'line':j+1,'text':line[:280]} for j,line in enumerate(lines) if exact120.search(line)]
  matches={'gl':gl,'archive':ARG.archive,'version':version,'path':q['path'],
    'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'original_crc32':f'{q["crc"]:08x}',
    'original_zip_local_offset':q['off'],
    'lines':len(lines),'matched_line_count':len(hits),'matched_lines_head':hits[:80],
    'exact_y120_line_count':len(coords),'exact_y120_lines_head':coords[:60],
    'contains_120_text':bool(coords),
    'source_is_original_authentic_per_GL_sidecar':True,
    'source_is_original_camera_acquisition_settings':False}
  records.append(matches)
  print('E511_VERIFIED_SIDECAR',json.dumps({k:matches[k] for k in ('gl','version','path','bytes','sha256','exact_y120_line_count','matched_line_count')}),flush=True)
(OUT/'E511_all_original_sidecar_source_findings.json').write_text(json.dumps({
  'archive':ARG.archive,
  'official_study':'BIA S-BIAD2506',
  'central_zip_directory_sha256':exp_sha,
  'original_manifest_sha256':'f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f',
  'gls':gls,'sources':records,
  'source_level_only_not_independent_image_segmentation':True,
  'method':'literal exact source per-GL settings and logs, case-insensitive token strings, no inference of y=120 unless explicit documented source fields'
 },indent=2,ensure_ascii=False),encoding='utf8')
print('E511_ORIGINAL_SOURCE_PARAMETER_AUDIT_PASS',ARG.archive,len(records),flush=True)
