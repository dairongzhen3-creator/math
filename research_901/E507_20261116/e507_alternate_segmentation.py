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
exp_sha,exp_n,gls=ARCH[ARG.archive];OUT=pathlib.Path('e507_alternate_QA')/ARG.archive[:-4];OUT.mkdir(parents=True,exist_ok=True)
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
candidates=collections.defaultdict(dict);i=0;n=0
while i<len(cd):
 h=struct.unpack_from('<4s6H3I5H2I',cd,i);assert h[0]==b'PK\x01\x02'
 nl,el,cl=h[10:13];path=cd[i+46:i+46+nl].decode('utf-8' if h[3]&2048 else 'cp437','replace')
 for gl in gls:
  if f'/{gl}/' not in path or '__BKP' in path:continue
  if path.endswith('.csv') and '/CellStats__' in path:
   for version,segment in [('main','/20221220_analysis/export_data__20221220_analysis/'),('cost','/cost_function_analysis_1/export_data__cost_function_analysis_1/')]:
    if segment in path:
     raw,com,off=wide(cd[i+46+nl:i+46+nl+el],h[9],h[8],h[16]);candidates[gl].setdefault(version,[]).append({'path':path,'raw':raw,'compressed':com,'off':off,'method':h[4],'crc':h[7]})
 i+=46+nl+el+cl;n+=1
assert i==length and n==count

def decode(entry):
 assert entry['compressed']<12000000 and entry['raw']<22000000
 off=entry['off'];h=struct.unpack('<4s5H3I2H',get(off,off+29));assert h[0]==b'PK\x03\x04'
 nl,el=h[-2:];d=get(off+30,off+30+nl+el-1)
 assert d[:nl].decode('utf-8' if h[2]&2048 else 'cp437','replace')==entry['path']
 compressed=get(off+30+nl+el,off+30+nl+el+entry['compressed']-1)
 original=zlib.decompress(compressed,-15) if entry['method']==8 else compressed if entry['method']==0 else None
 assert original is not None and len(original)==entry['raw'] and zlib.crc32(original)==entry['crc']
 return original

def source_events(rows):
 bycell=collections.defaultdict(list)
 for rec in rows:bycell[int(rec['cell_ID'])].append(rec)
 divs=[]
 for pid,pr in bycell.items():
  pr.sort(key=lambda r:int(r['frame']))
  for r in pr:
   if r['type_of_end']!='div':continue
   fr=int(r['frame']);dau=[k for k,v in bycell.items() if any(int(x['parent_ID'])==pid and int(x['frame'])==fr+1 for x in v)]
   if len(dau)!=2:continue
   lo=int(r['bbox_top_px']);hi=int(r['bbox_bottom_px']);center=.5*(lo+hi)
   divs.append({'parent':pid,'frame':fr,'bbox_top':lo,'bbox_bottom':hi,'center_y':center,'length':float(r['length_px']),'daughters':sorted(dau)})
 return divs

outputs=[];receipts=[]
for gl in gls:
 c=candidates[gl];assert len(c.get('main',[]))==1 and len(c.get('cost',[]))==1,(gl,{k:len(v) for k,v in c.items()})
 per={}
 for version in ('main','cost'):
  entry=c[version][0];b=decode(entry)
  f=OUT/f'{gl}_{version}_CellStats.csv';f.write_bytes(b)
  lines=[t for t in b.decode('utf8-sig').splitlines() if t.strip() and not t.startswith('#')]
  data=list(csv.DictReader(lines));assert len(data)>100
  events=source_events(data)
  per[version]={'records':len(data),'unique_cells':len({d['cell_ID'] for d in data}), 'unique_frames':len({d['frame'] for d in data}),'division_events':events}
  receipts.append({'archive':ARG.archive,'gl':gl,'version':version,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'original_crc32':f'{entry["crc"]:08x}','source_path':entry['path'],'original_offset':entry['off']})
 # Use no numeric ID matching across analyses; only equal time and 1D physical coordinate tolerance (5 px center/8 px length).
 # Focus only f in [200, 700] to match full relevant starvation interval and avoid single-frame left censoring.
 main=[x for x in per['main']['division_events'] if 200<=x['frame']<=700]
 cost=[x for x in per['cost']['division_events'] if 200<=x['frame']<=700]
 pairs=[];unused=set(range(len(cost)))
 for ev in main:
  acceptable=[(abs(ev['center_y']-cost[i]['center_y'])+abs(ev['length']-cost[i]['length']),i) for i in unused if ev['frame']==cost[i]['frame'] and abs(ev['center_y']-cost[i]['center_y'])<=5 and abs(ev['length']-cost[i]['length'])<=8]
  if acceptable:
   _,i=min(acceptable);unused.remove(i);pairs.append({'frame':ev['frame'],'main_parent':ev['parent'],'cost_parent':cost[i]['parent'],'center_delta_px':ev['center_y']-cost[i]['center_y'],'length_delta_px':ev['length']-cost[i]['length']})
 out={'GL':gl,'archive':ARG.archive,'main_rows':per['main']['records'],'cost_rows':per['cost']['records'],'main_unique_cell_IDs':per['main']['unique_cells'],'cost_unique_cell_IDs':per['cost']['unique_cells'], 'main_total_divisions_all_frames':len(per['main']['division_events']),'cost_total_divisions_all_frames':len(per['cost']['division_events']),'main_divisions_frames200to700':len(main),'cost_divisions_frames200to700':len(cost),'matched_geometry_divisions':len(pairs),'unmatched_main':len(main)-len(pairs),'unmatched_cost':len(cost)-len(pairs), 'pairs':pairs,'comparison':'Exact frame and matched physical bbox center <=5 pixels and length<=8 pixels; same ID across algorithms is never assumed; matching threshold pre-frozen diagnostic, not truth.'}
 outputs.append(out)
 print('E507_ALTERNATE',json.dumps({k:v for k,v in out.items() if k!='pairs'}),flush=True)
with (OUT/'E507_alternate_algorithm_frame_geometry_results.json').open('w') as f:json.dump(outputs,f,indent=2,ensure_ascii=False)
with (OUT/'E507_alt_original_member_source_sha.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(receipts[0]));w.writeheader();w.writerows(receipts)
print('E507_ALT_ARCHIVE_PASS',ARG.archive,'GLs',len(outputs))
