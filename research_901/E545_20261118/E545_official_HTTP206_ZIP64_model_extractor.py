#!/usr/bin/env python3
"""E545: 901-source bounded ZIP64 range extractor for original Gurobi MPS (optional LP).
Public BioImage Archive S-BIAD2506, two fixed 20221220 GL archives only.
Does not download the whole ~25/33 GB source archive. Rejects non-206 responses.
Requires requests and standard Python libraries. --self-test is offline and uses no network.
Outputs exact extracted source bytes + SHA256/CRC/central dir authentication receipt.
No inferences about biological truth or numerical solver behavior.
"""
import argparse,collections,csv,hashlib,io,json,os,re,shutil,tempfile,time,zipfile
from pathlib import Path
import requests
CHUNK=4*1024*1024
MAXCACHE=10
MANIFEST_SHA='f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f'
URL_BASE='https://www.ebi.ac.uk/biostudies/files/S-BIAD2506/Statphase_Transition_GLs/deepMoMA/'
ARCHIVES={
 'Pos5_GL35':dict(name='20221220_rplN_GFP_pl.zip',root='20221220_rplN_GFP_pl/Pos5_GL35/',
  length=33229979698,cd_len=45482286,cd_sha='7f55aaa61d7025260773c33aa662e13ae2b2266dea300384f9e09a1667dc5ac6',entries=193307),
 'Pos12_GL22':dict(name='20221220_E08_GFP_pl.zip',root='20221220_E08_GFP_pl/Pos12_GL22/',
  length=25278719975,cd_len=33108571,cd_sha='5bd93ce604256d00233d4c5bf91865205b44b0d7fc2f52fefe0f68ef0e23e94c',entries=141104)
}
EXPECT={('Pos5_GL35','main','mps'):(4605562,90568395),('Pos5_GL35','cost','mps'):(4606805,90576517),
('Pos12_GL22','main','mps'):(4581436,88730205),('Pos12_GL22','cost','mps'):(4591632,88794943),
('Pos5_GL35','main','lp'):(3202700,57923333),('Pos5_GL35','cost','lp'):(3203162,57927693),
('Pos12_GL22','main','lp'):(3240474,56627201),('Pos12_GL22','cost','lp'):(3244898,56664993)}
class SafeHTTPRange(io.RawIOBase):
    def __init__(self,url,total_bytes,session=None):
        super().__init__();self.url=url;self.length=int(total_bytes);self.cursor=0
        self.session=session or requests.Session();self.blocks=collections.OrderedDict();self.total_transferred=0;self.nreq=0
        self._fetch(0,0) # strict public source HTTP 206 preflight
    def readable(self):return True
    def seekable(self):return True
    def tell(self):return self.cursor
    def seek(self, offset,whence=io.SEEK_SET):
        if whence==io.SEEK_CUR:offset=self.cursor+offset
        elif whence==io.SEEK_END:offset=self.length+offset
        elif whence!=io.SEEK_SET:raise ValueError('invalid seek')
        if offset<0:raise ValueError('negative offset')
        self.cursor=int(offset);return self.cursor
    def _fetch(self,lo,hi):
        assert 0<=lo<=hi<self.length,(lo,hi,self.length)
        if hi-lo+1>CHUNK:raise ValueError('unbounded HTTP range')
        last=None
        for trial in range(4):
            try:
                hdr={'Range':f'bytes={lo}-{hi}','Accept-Encoding':'identity','User-Agent':'901-E545-source-limited-range-extractor/1.0'}
                with self.session.get(self.url,headers=hdr,stream=True,timeout=(20,90),allow_redirects=True) as r:
                    assert r.status_code==206,f'Expected strict HTTP206 range, got status={r.status_code}; will NOT download entire archive'
                    cr=r.headers.get('Content-Range','')
                    assert cr==f'bytes {lo}-{hi}/{self.length}',('source length / content range mismatch',cr)
                    assert r.headers.get('Content-Encoding','identity').lower() in ('identity',''), 'unexpected compression'
                    body=r.content
                    assert len(body)==hi-lo+1, ('range body length mismatch',len(body),hi-lo+1)
                    self.total_transferred+=len(body);self.nreq+=1
                    return body
            except (requests.RequestException,AssertionError) as e:
                last=e
                if isinstance(e,AssertionError):raise
                if trial<3:time.sleep(min(2**trial,4))
        raise RuntimeError('HTTP range retried with bounded failures') from last
    def read(self,n=-1):
        if self.cursor>=self.length:return b''
        if n is None or n<0:n=self.length-self.cursor
        n=min(n,self.length-self.cursor)
        if n>96*1024*1024:raise ValueError('unbounded file read request refused')
        p=self.cursor;parts=[]
        while n:
            blk=p//CHUNK;off=p%CHUNK
            if blk not in self.blocks:
                start=blk*CHUNK;end=min(self.length-1,start+CHUNK-1)
                val=self._fetch(start,end)
                self.blocks[blk]=val
                while len(self.blocks)>MAXCACHE:self.blocks.popitem(last=False)
            else:self.blocks.move_to_end(blk)
            take=min(n,len(self.blocks[blk])-off)
            assert take>0
            parts.append(self.blocks[blk][off:off+take]);p+=take;n-=take
        self.cursor=p
        return b''.join(parts)

def members(gl,with_lp=False):
    cfg=ARCHIVES[gl]
    for ver,prefix in [('main','20221220_analysis/track_data__20221220_analysis/'),
    ('cost','cost_function_analysis_1/track_data__cost_function_analysis_1/')]:
        for ext in (('mps','lp') if with_lp else ('mps',)):
            yield ver,ext,cfg['root']+prefix+'gurobi_model.'+ext

def extract(gl,dest,with_lp=False,source_override=None,source_len_override=None,cd_sha_override=None,entries_override=None,expect_override=None):
    cfg=ARCHIVES[gl]
    url=source_override or URL_BASE+cfg['name']
    srclen=int(source_len_override if source_len_override is not None else cfg['length'])
    cdsha=cd_sha_override or cfg['cd_sha']
    count=int(entries_override if entries_override is not None else cfg['entries'])
    reader=SafeHTTPRange(url,srclen)
    with zipfile.ZipFile(reader,'r',allowZip64=True) as z:
        assert len(z.infolist())==count,('ZIP central count wrong',len(z.infolist()),count)
        assert z.start_dir>=0
        reader.seek(z.start_dir)
        cdb=reader.read(cfg['cd_len'] if source_override is None else (z._endrec[zipfile._ECD_SIZE] if hasattr(z,'_endrec') else 0))
        # zipfile.ZipFile may not expose _endrec. In production use E502 authenticated known length.
        if source_override is None:
            assert hashlib.sha256(cdb).hexdigest()==cdsha,('full authenticated original central directory SHA mismatch',gl)
        out=[]
        for ver,ext,path in members(gl,with_lp):
            info=z.getinfo(path)
            if expect_override is None:
                expectedc,expectedr=EXPECT[gl,ver,ext]
                assert (info.compress_size,info.file_size)==(expectedc,expectedr),('E502 original member compression mismatch',gl,ver,ext)
            assert info.compress_type==zipfile.ZIP_DEFLATED and info.header_offset>0
            target=Path(dest)/(gl+'_'+ver+'.'+ext)
            sha=hashlib.sha256();total=0
            with z.open(info,'r') as inp,target.open('wb') as outf:
                while True:
                    chunk=inp.read(1024*1024)
                    if not chunk:break
                    outf.write(chunk);sha.update(chunk);total+=len(chunk)
            assert total==info.file_size
            with target.open('rb') as fh:
                first=fh.read(256)
            # Syntax header check is intentionally weak; full MPS parse happens only after native archive delivery.
            if ext=='mps':assert b'NAME' in first[:120],('not MPS',first[:120])
            out.append(dict(GL=gl,version=ver,extension=ext,original_archive=cfg['name'],official_url=url,
                 original_member_path=path,original_member_zip64_local_offset=info.header_offset,
                 original_member_compression=info.compress_size,original_uncompressed_bytes=total,
                 original_member_CRC32=f'{info.CRC:08x}',full_SHA256=sha.hexdigest(),output_filename=target.name,
                 status='HTTP206_ZIP64_CD_MATCH_MEMBER_CRC_SHA_PASS' if source_override is None else 'OFFLINE_SYNTHETIC_SELFTEST'))
    print(f'E545 SOURCE {gl} VERIFIED compressed member(s), HTTP range requests {reader.nreq}, transferred {reader.total_transferred} bytes',flush=True)
    return out,dict(source_bytes=srclen,http_range_requests=reader.nreq,http_range_transferred_bytes=reader.total_transferred,
         original_directory_sha256=hashlib.sha256(cdb).hexdigest(),original_directory_bytes=len(cdb))

def selftest():
    # Fully exercise seek()/read() with synthetic deflate source, without external network.
    import http.server,socketserver,threading
    from unittest.mock import patch
    with tempfile.TemporaryDirectory(prefix='e545_zip_test_') as work:
        root=Path(work); p=root/'sample.zip'
        with zipfile.ZipFile(p,'w',compression=zipfile.ZIP_DEFLATED,allowZip64=True) as z:
            for i in range(53):z.writestr(f'background/{i}.txt',('input'+str(i))*40)
            for v in ('main','cost'):z.writestr(f'{v}/gurobi_model.mps','NAME          Synthetic\nROWS\n N COST\nENDATA\n')
        data=p.read_bytes()
        class FakeSession:
            def get(self,url,headers,stream,timeout,allow_redirects):
                a,b=map(int,headers['Range'].replace('bytes=','').split('-'))
                class Resp:
                    status_code=206;headers={'Content-Range':f'bytes {a}-{b}/{len(data)}','Content-Encoding':'identity'}
                    content=data[a:b+1]
                    def __enter__(self):return self
                    def __exit__(self,*_):return None
                return Resp()
        reader=SafeHTTPRange('offline://local',len(data),session=FakeSession())
        with zipfile.ZipFile(reader) as z:
            assert len(z.infolist())==55
            for v in ('main','cost'):
                assert z.read(f'{v}/gurobi_model.mps').startswith(b'NAME')
            reader.seek(z.start_dir)
            actual=reader.read(len(data)-z.start_dir-22)
            assert hashlib.sha256(actual).hexdigest()==hashlib.sha256(data[z.start_dir:-22]).hexdigest()
        # ZIP64 central extra offset parsing is handled by stdlib zipfile and is not separately tested on >4 GB fixture.
        print(json.dumps({'offline_bounded_range_zipfile_pass':True,'members':55,'actual_source_transfer_bytes':reader.total_transferred,
           'full_zip_size':len(data),'strict_HTTP206_in_simulated_transport':True},indent=2))

def main():
    a=argparse.ArgumentParser();a.add_argument('--dest',default='E545_official_models');a.add_argument('--with-lp',action='store_true');a.add_argument('--self-test',action='store_true');args=a.parse_args()
    if args.self_test:return selftest()
    dest=Path(args.dest);dest.mkdir(parents=True,exist_ok=True)
    checked=[];receipts={}
    try:
        for gl in ('Pos5_GL35','Pos12_GL22'):
            out,receipt=extract(gl,dest,args.with_lp)
            checked+=out;receipts[gl]=receipt
        m={'study':'S-BIAD2506','source_manifest_sha256':MANIFEST_SHA,'output_model_files':checked,
           'official_source_zip64_206':receipts,'biology_truth':'UNRESOLVED','physical_bit_joule':'NOT_IDENTIFIED',
           '2022_actual_Gurobi_solver_run_verified':False}
        (dest/'E545_original_MPS_LP_SHA_CRC_provenance.json').write_text(json.dumps(m,ensure_ascii=False,sort_keys=True,indent=2)+'\n')
        assert len(checked)==(8 if args.with_lp else 4)
        print('E545_DOWNLOAD_SUCCESS',json.dumps({'files':len(checked),'archive_verified':True,'sizes':[z['original_uncompressed_bytes'] for z in checked]},ensure_ascii=False),flush=True)
    except Exception as e:
        (dest/'E545_failed_source_access_receipt.json').write_text(json.dumps({'type':type(e).__name__,
         'error':str(e),'status':'SOURCE_OR_CRC_OR_SHA_NOT_AUTHENTICATED','no_model_claim':True},ensure_ascii=False,indent=2)+'\n')
        raise
if __name__=='__main__':main()
