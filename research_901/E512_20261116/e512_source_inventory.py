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
exp_sha,exp_n,gls=ARCH[ARG.archive];OUT=pathlib.Path('e512_verified_inventory')/ARG.archive[:-4];OUT.mkdir(parents=True,exist_ok=True)
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

# E512 targets only frozen official per-GL binary/tracking source members.
GLS={"20221220_E08_GFP_pl.zip":["Pos12_GL22"],"20221220_rplN_GFP_pl.zip":["Pos5_GL35"]}
assert ARG.archive in GLS
wanted_names={
  'probability_maps.tif', 'component_forests.json', 'curation.moma',
  'hypothesis_states.csv','assignment_states.csv','assignment_costs.csv',
  'assignment_filter_intensities.csv','file_format.json','gurobi_model.sol','gurobi_model.mst'
}
inventory=[];by_key=collections.defaultdict(list);j=0;n=0
while j<len(cd):
 h=struct.unpack_from('<4s6H3I5H2I',cd,j);assert h[0]==bytes.fromhex('504b0102')
 nl,el,cl=h[10:13];path=cd[j+46:j+46+nl].decode('utf8' if h[3]&2048 else 'cp437','replace')
 for gl in GLS[ARG.archive]:
  if f'/{gl}/' not in path or '__BKP' in path:continue
  basename=pathlib.PurePosixPath(path).name
  if basename not in wanted_names or '/track_data__' not in path:continue
  if '/20221220_analysis/track_data__20221220_analysis/' in path:version='main'
  elif '/cost_function_analysis_1/track_data__cost_function_analysis_1/' in path:version='cost'
  else:continue
  raw,comp,loc=wide(cd[j+46+nl:j+46+nl+el],h[9],h[8],h[16])
  q={'archive':ARG.archive,'gl':gl,'version':version,'kind':basename,'source_path':path,'raw_bytes':raw,'compressed_bytes':comp,'original_crc32':f'{h[7]:08x}','zip_method':h[4],'zip_offset':loc}
  inventory.append(q);by_key[(gl,version,basename)].append(q)
 j+=46+nl+el+cl;n+=1
assert j==length and n==count
for gl in GLS[ARG.archive]:
 for v in ['main','cost']:
  for kind in ['probability_maps.tif','component_forests.json','curation.moma','file_format.json']:
   assert len(by_key[(gl,v,kind)])==1,(gl,v,kind)
  print('E512_OFFICIAL_ORIGINAL_MEMBER_INVENTORY',json.dumps({
    'gl':gl,'version':v,'entries':[{'kind':t['kind'],'raw_bytes':t['raw_bytes'],'compressed_bytes':t['compressed_bytes'],'crc':t['original_crc32']} for t in inventory if t['gl']==gl and t['version']==v]},ensure_ascii=False),flush=True)
receipt={'archive':ARG.archive,'official_deepMOMA_sha256':'f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f','central_directory_sha256':exp_sha,
 'central_dir_verified_member_count':n,'complete_central_dir_sha_verified':True,
 'source_member_header_not_yet_read':True,
 'raw_member_crc_not_yet_recomputed':True,'items':inventory,
 'scope':'source inventory, binary payload must still be authenticated before scientific use'}
(OUT/'E512_verified_original_member_inventory.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False),encoding='utf-8')
print('E512_INVENTORY_PASS',ARG.archive,len(inventory),flush=True)
