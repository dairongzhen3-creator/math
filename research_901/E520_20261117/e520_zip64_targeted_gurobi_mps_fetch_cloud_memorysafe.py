#!/usr/bin/env python3
"""901 E520: strictly source-authenticated remote ZIP64 extraction of archived native Gurobi .mps/.lp.

Never claim acquisition from synthetic/offline tests. The only scientifically valid
payload is from the original BioImage Archive original large ZIP64 whose central
directory hash, manifest hash and each member local header/CRC32/SHA are checked.

Run on networked cloud machine:
  python e520_zip64_targeted_gurobi_mps_fetch.py --archive 20221220_rplN_GFP_pl.zip --output results
Or deterministic TEST ONLY, never as source data:
  python e520_zip64_targeted_gurobi_mps_fetch.py --archive 20221220_rplN_GFP_pl.zip --local-fixture /tmp/fixture.zip --output /tmp/test
"""
import argparse,hashlib,json,pathlib,struct,zlib,time,zipfile,io,sys

A=argparse.ArgumentParser()
A.add_argument('--archive',required=True,choices=['20221220_rplN_GFP_pl.zip','20221220_E08_GFP_pl.zip'])
A.add_argument('--output',required=True)
A.add_argument('--local-fixture',default=None,help='synthetic/offline TEST ONLY; cannot establish original source identity')
a=A.parse_args()
CFG={
 '20221220_rplN_GFP_pl.zip':{'gl':'Pos5_GL35','count':193307,'cdir_sha':'7f55aaa61d7025260773c33aa662e13ae2b2266dea300384f9e09a1667dc5ac6'},
 '20221220_E08_GFP_pl.zip':{'gl':'Pos12_GL22','count':141104,'cdir_sha':'5bd93ce604256d00233d4c5bf91865205b44b0d7fc2f52fefe0f68ef0e23e94c'}
}
C=CFG[a.archive];GL=C['gl'];OUT=pathlib.Path(a.output);OUT.mkdir(parents=True,exist_ok=True)
MAN_SHA='f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f'
BASE='https://www.ebi.ac.uk/biostudies/files/S-BIAD2506/'
S=None
if a.local_fixture:
    LOCAL=pathlib.Path(a.local_fixture).open('rb');N=LOCAL.seek(0,2);URL='LOCAL_SYNTHETIC_TEST_ONLY'
else:
    import requests
    S=requests.Session(); S.headers['User-Agent']='901-E520-original-moma-gurobi-MPS-LP-targeted-zip64-source-20261117'
    resp=S.get(BASE+'Statphase_Transition_GLs/deepMoMA.json',timeout=(12,65));resp.raise_for_status()
    assert hashlib.sha256(resp.content).hexdigest()==MAN_SHA
    m={pathlib.PurePosixPath(v['path']).name:v for v in json.loads(resp.content)}
    assert a.archive in m and m[a.archive]['path'].endswith(a.archive)
    N=int(m[a.archive]['size']);URL=BASE+m[a.archive]['path']

def chunks(start,stop):
    assert 0<=start<stop<=N
    for first in range(start,stop,4194304):
        last=min(first+4194304,stop)
        if a.local_fixture:
            LOCAL.seek(first);b=LOCAL.read(last-first)
            if len(b)!=last-first:raise IOError('local read truncation')
            yield b;continue
        err=None
        for attempt in range(3):
            try:
                r=S.get(URL,headers={'Range':f'bytes={first}-{last-1}'},stream=True,timeout=(20,75))
                if r.status_code!=206 or r.headers.get('Content-Range')!=f'bytes {first}-{last-1}/{N}':
                    raise RuntimeError(('strict HTTP206 range failed',r.status_code,r.headers.get('Content-Range')))
                b=bytearray()
                for part in r.iter_content(131072):
                    if part:b.extend(part)
                    if len(b)>last-first:raise RuntimeError('oversized HTTP206 response')
                if len(b)!=last-first:raise RuntimeError('truncated HTTP206 response')
                yield bytes(b);break
            except Exception as e:
                err=e;time.sleep(1+attempt)
        else:raise err

def read(start,stop):return b''.join(chunks(start,stop))

def wide(extra,raw,comp,off):
    if all(v!=0xffffffff for v in (raw,comp,off)):return raw,comp,off
    i=0
    while i<len(extra):
        tag,n=struct.unpack_from('<HH',extra,i);d=extra[i+4:i+4+n];i+=4+n
        if tag==1:
            j=0
            if raw==0xffffffff:raw=struct.unpack_from('<Q',d,j)[0];j+=8
            if comp==0xffffffff:comp=struct.unpack_from('<Q',d,j)[0];j+=8
            if off==0xffffffff:off=struct.unpack_from('<Q',d,j)[0];j+=8
            return raw,comp,off
    raise ValueError('Missing ZIP64 extra field')

tailstart=max(0,N-131072);tail=read(tailstart,N)
i=tail.rfind(b'PK\x05\x06');assert i>=0,'No EOCD'
x=struct.unpack_from('<4s4H2IH',tail,i);count,length,offset=x[4],x[5],x[6]
if count==65535 or length==0xffffffff or offset==0xffffffff:
    loc=tail.rfind(b'PK\x06\x07',0,i);assert loc>=0,'No ZIP64 locator'
    z64off=struct.unpack_from('<Q',tail,loc+8)[0];h=read(z64off,z64off+56);assert h[:4]==b'PK\x06\x06'
    count=struct.unpack_from('<Q',h,32)[0]
    length=struct.unpack_from('<Q',h,40)[0]
    offset=struct.unpack_from('<Q',h,48)[0]
assert 0<count<1000000 and 0<length<100000000 and 0<=offset<N and offset+length<=N
cd=read(offset,offset+length);dhash=hashlib.sha256(cd).hexdigest()
if not a.local_fixture:
    assert count==C['count'] and dhash==C['cdir_sha'],(count,dhash)
entries={};j=0;n=0
while j<len(cd):
    h=struct.unpack_from('<4s6H3I5H2I',cd,j);assert h[0]==b'PK\x01\x02'
    nl,el,cl=h[10:13]
    name=cd[j+46:j+46+nl].decode('utf8' if h[3]&2048 else 'cp437','replace')
    if name.endswith(('gurobi_model.lp','gurobi_model.mps')) and f'/{GL}/' in name and '__BKP' not in name:
        version='main' if '/20221220_analysis/track_data__20221220_analysis/' in name else 'cost' if '/cost_function_analysis_1/track_data__cost_function_analysis_1/' in name else None
        if version is not None:
            raw,compressed,local_offset=wide(cd[j+46+nl:j+46+nl+el],h[9],h[8],h[16])
            kind=name.rsplit('.',1)[-1]; key=(version,kind)
            assert key not in entries
            entries[key]={'kind':kind,'version':version,'path':name,'raw_bytes':raw,
                  'compressed_bytes':compressed,'crc32':h[7],'local_header_offset':local_offset,'method':h[4],
                  'name_flags':h[3]}
    j+=46+nl+el+cl;n+=1
assert j==length and n==count and len(entries)==4,(n,len(entries),entries)
# Require exact previously inventoried paths, not basename-only matches.
for (v,kind),e in entries.items():
    expected=a.archive[:-4]+'/'+GL+'/'+('20221220_analysis/track_data__20221220_analysis' if v=='main' else 'cost_function_analysis_1/track_data__cost_function_analysis_1')+'/gurobi_model.'+kind
    assert e['path']==expected,(e['path'],expected)
    assert e['method'] in [0,8] and e['raw_bytes']<2_000_000_000 and e['compressed_bytes']<1_000_000_000

receipt={'source_study':'S-BIAD2506','source_archive':a.archive,'source_url':URL,
    'biostudies_manifest_sha256':MAN_SHA if not a.local_fixture else 'SKIPPED_SYNTHETIC_TEST',
    'verified_original_central_directory_sha256':dhash,'central_directory_entry_count':n,
    'identity':'UNAUTHENTICATED_SYNTHETIC_TEST' if a.local_fixture else 'OFFICIAL_ARCHIVE_SHA_VERIFIED',
    'members':[]}
for (ver,kind),e in sorted(entries.items()):
    lo=e['local_header_offset'];h=struct.unpack('<4s5H3I2H',read(lo,lo+30));assert h[0]==b'PK\x03\x04'
    name_len,extra_len=h[-2:];raw_name=read(lo+30,lo+30+name_len)
    name=raw_name.decode('utf8' if h[2]&2048 else 'cp437','replace');assert name==e['path']
    first=lo+30+name_len+extra_len
    dst=OUT/(GL+'__'+ver+'__gurobi_model.'+kind)
    hasher=hashlib.sha256();crc=0;nraw=0
    dec=zlib.decompressobj(-15) if e['method']==8 else None
    with dst.open('wb') as f:
        for block in chunks(first,first+e['compressed_bytes']):
            un=dec.decompress(block) if dec else block
            f.write(un);nraw+=len(un);crc=zlib.crc32(un,crc);hasher.update(un)
            assert nraw<=e['raw_bytes'],'decompressed beyond declared original'
        if dec:
            un=dec.flush();f.write(un);nraw+=len(un);crc=zlib.crc32(un,crc);hasher.update(un)
            assert dec.eof, 'compressed member not ended'
    assert nraw==e['raw_bytes'] and crc==e['crc32'],('SRC_CRC_OR_LENGTH_MISMATCH',e)
    if kind=='mps':
        with dst.open('rb') as f: sample=f.read(8192)
        assert b'NAME' in sample and b'ROWS' in sample,('MPS_FORMAT_HEADER',str(dst))
    if kind=='lp':
        with dst.open('rb') as f: sample=f.read(1024)
        assert any(x in sample.lower() for x in [b'minimize',b'maximize',b'subject to']),('LP_FORMAT_HEADER',str(dst))
    rec=dict(e);rec.update({'original_zip_member_crc32_verified':True,'original_zip_member_crc32_hex':f'{crc:08x}',
        'full_original_SHA256':hasher.hexdigest(),'decompressed_size_reverified':nraw,'saved_filename':dst.name})
    receipt['members'].append(rec)
    print('E520_NATIVE_MPS_LP_SOURCE_PASS',json.dumps({'gl':GL,'version':ver,'format':kind,'bytes':nraw,'sha256':hasher.hexdigest(),'original_crc':f'{crc:08x}','status':receipt['identity']}),flush=True)
(OUT/(GL+'__original_gurobi_model_MPS_LP_ZIP64_source_receipt.json')).write_text(json.dumps(receipt,indent=2,ensure_ascii=False),encoding='utf8')
if a.local_fixture:print('E520_MPS_LP_SYNTHETIC_TEST_ONLY_SUCCESS',GL,flush=True)
else:print('E520_MPS_LP_REAL_OFFICIAL_SOURCE_AUTHENTICATED',GL,flush=True)
