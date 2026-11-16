#!/usr/bin/env python3
"""901 E505: constrained original MoMA CellTracks / CellStats source inspection.
Only 4 source GLs. Strict 206 Content-Range, ZIP64 offset, CRC32 and SHA.
Source-identity questions are separate from purported biological pedigree.
"""
import requests,struct,hashlib,json,zlib,re,pathlib,time,csv,io,collections
O=pathlib.Path("e505_pilot_originals");O.mkdir(exist_ok=True)
S=requests.Session();S.headers["User-Agent"]="901-E505-authentic-MoMA-track-csv-bounded-20261116"
BASE="https://www.ebi.ac.uk/biostudies/files/S-BIAD2506/"
MAN_URL=BASE+"Statphase_Transition_GLs/deepMoMA.json"
m=S.get(MAN_URL,timeout=(15,45));m.raise_for_status()
assert len(m.content)==37313 and hashlib.sha256(m.content).hexdigest()=="f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f"
manifest=json.loads(m.content)
targets={
"20221220_E08_GFP_pl.zip":{"expected_cd":"5bd93ce604256d00233d4c5bf91865205b44b0d7fc2f52fefe0f68ef0e23e94c","n":141104,"GLs":["Pos10_GL16","Pos10_GL11"]},
"20221220_NO_LIGHT_E08_GFP_pl.zip":{"expected_cd":"bae990f0288f2898ab3c9525181cbceb3637cf2a08f1a120ec5292fad597f7fb","n":16967,"GLs":["Pos16_GL14","Pos16_GL27"]}}
T={pathlib.PurePosixPath(v["path"]).name:v for v in manifest if pathlib.PurePosixPath(v["path"]).name in targets}
assert set(T)==set(targets)
rec={"study":"S-BIAD2506","date":"20221220","stage":"E505_pilot","source_manifest_sha":hashlib.sha256(m.content).hexdigest(),"datetime_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),"files":[],"archives":[],"max_single_range_bytes":4194304,"download_whole_images":False}
def GET(url,a,b,total):
    assert 0<=a<=b<total and b-a+1<=4194304
    rr=S.get(url,headers={"Range":f"bytes={a}-{b}"},stream=True,timeout=(15,65))
    if rr.status_code!=206: raise ValueError(f"HTTP_NON206_{rr.status_code}")
    if rr.headers.get("Content-Range")!=f"bytes {a}-{b}/{total}":raise ValueError("CONTENT_RANGE_BAD")
    buf=bytearray()
    for part in rr.iter_content(131072):
        if part:buf.extend(part)
        if len(buf)>b-a+1:raise ValueError("RANGE_OVER")
    assert len(buf)==b-a+1
    return bytes(buf)
def resolve64(extra,raw,compressed,offset):
    if not any(v==0xffffffff for v in (raw,compressed,offset)):return raw,compressed,offset
    i=0
    while i<len(extra):
        typ,n=struct.unpack_from("<HH",extra,i);v=extra[i+4:i+4+n];i+=n+4
        if typ!=1:continue
        j=0
        if raw==0xffffffff:raw=struct.unpack_from("<Q",v,j)[0];j+=8
        if compressed==0xffffffff:compressed=struct.unpack_from("<Q",v,j)[0];j+=8
        if offset==0xffffffff:offset=struct.unpack_from("<Q",v,j)[0];j+=8
        return raw,compressed,offset
    raise ValueError("ZIP64_EXTRA_NOT_FOUND")
for an,expect in targets.items():
    meta=T[an];n=int(meta["size"]);url=BASE+meta["path"]
    tail=GET(url,n-131072,n-1,n)
    pos=tail.rfind(b"PK\x05\x06");assert pos>=0
    end=struct.unpack_from("<4s4H2IH",tail,pos)
    count,dlen,dstart=end[4],end[5],end[6]
    if count==65535 or dlen==0xffffffff or dstart==0xffffffff:
        loc=tail.rfind(b"PK\x06\x07",0,pos);assert loc>=0
        zoff=struct.unpack_from("<Q",tail,loc+8)[0]
        z64=GET(url,zoff,zoff+55,n)
        assert z64[:4]==b"PK\x06\x06"
        count=struct.unpack_from("<Q",z64,32)[0]
        dlen=struct.unpack_from("<Q",z64,40)[0]
        dstart=struct.unpack_from("<Q",z64,48)[0]
    assert 0<dlen<=67108864 and count==expect["n"]
    cd=bytearray()
    for start in range(dstart,dstart+dlen,4194304):
        cd.extend(GET(url,start,min(dstart+dlen-1,start+4194303),n))
    assert len(cd)==dlen and hashlib.sha256(cd).hexdigest()==expect["expected_cd"]
    entries=[];idx=0
    while idx<len(cd):
        head=struct.unpack_from("<4s6H3I5H2I",cd,idx);assert head[0]==b"PK\x01\x02"
        nl,el,cl=head[10:13]
        name=cd[idx+46:idx+46+nl].decode("utf-8" if head[3]&2048 else "cp437","replace")
        extra=cd[idx+46+nl:idx+46+nl+el]
        if (
          any(f"/{gl}/" in name for gl in expect["GLs"])
          and "/20221220_analysis/export_data__20221220_analysis/" in name
          and re.search(r"/Cell(?:Tracks|Stats)__[^/]+\.csv$",name)
        ):
            raw,com,off=resolve64(extra,head[9],head[8],head[16])
            entries.append({"name":name,"raw":raw,"compressed":com,"offset":off,"method":head[4],"crc":head[7]})
        idx+=46+nl+el+cl
    assert idx==dlen and count==expect["n"]
    assert len(entries)==len(expect["GLs"])*2,(an,len(entries),[x["name"] for x in entries])
    rec["archives"].append({"archive":an,"source_size":n,"zip64_cdir_sha256":expect["expected_cd"],"members_checked":count,"selected_files":len(entries)})
    for entry in entries:
        assert entry["compressed"]<4194304 and entry["raw"]<12000000
        off=entry["offset"];h=struct.unpack("<4s5H3I2H",GET(url,off,off+29,n))
        assert h[0]==b"PK\x03\x04"
        nl,el=h[-2:]
        hh=GET(url,off+30,off+30+nl+el-1,n)
        assert hh[:nl].decode("utf8" if h[2]&2048 else "cp437","replace")==entry["name"]
        start=off+30+nl+el
        payload=GET(url,start,start+entry["compressed"]-1,n)
        if entry["method"]==8:data=zlib.decompress(payload,-15)
        elif entry["method"]==0:data=payload
        else:raise ValueError("ZIP_UNSUPPORTED")
        assert len(data)==entry["raw"] and zlib.crc32(data)==entry["crc"]
        parts=entry["name"].split("/")
        gl=next(x for x in expect["GLs"] if f"/{x}/" in entry["name"])
        typ="CellTracks" if "/CellTracks__" in entry["name"] else "CellStats"
        simple=f"{an[:-4]}__{gl}__{typ}.csv"
        (O/simple).write_bytes(data)
        # evidence: present source CSV column names and first three data records, no guessed schema
        samples=list(csv.reader(io.StringIO(data.decode("utf8-sig","replace"))))[:4]
        rr={"source_archive":an,"gl":gl,"kind":typ,"original_path":entry["name"],
            "original_csv_bytes":len(data),"sha256":hashlib.sha256(data).hexdigest(),"crc32":format(zlib.crc32(data),"08x"),
            "csv_header":samples[0] if samples else [],"first_3_rows":samples[1:4],"pilot_only":True}
        rec["files"].append(rr)
        print("E505_EXACT_PILOT",json.dumps({k:rr[k] for k in ["source_archive","gl","kind","original_csv_bytes","sha256","csv_header"]},ensure_ascii=False),flush=True)
(O/"E505_pilot_raw_origin_and_schema.json").write_text(json.dumps(rec,ensure_ascii=False,indent=2))
print("E505_PILOT_FINISHED",len(rec["files"]),"CSV_FULL_CRC",flush=True)
