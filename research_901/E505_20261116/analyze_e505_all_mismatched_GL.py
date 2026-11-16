#!/usr/bin/env python3
"""901 E505 source-verified original MoMA parent_ID lineage for E504's 118 GL ID changes.
One original BIA archive per job; original CSV bytes SHA/CRC & ZIP64 central-directory pass.
No source image pixels, no unverified string genealogy, no causal or molecular-bit claims.
"""
import argparse,collections,csv,hashlib,io,json,pathlib,re,struct,time,zlib
import requests
p=argparse.ArgumentParser();p.add_argument("--archive",required=True);args=p.parse_args()
AR=args.archive
EXPECTED={
"20221220_E08_GFP.zip":("352c0f14d8cc41f64942c3b850cde4309971c9654e814a722440f0b7c2ce3dbc",78458),
"20221220_E08_GFP_pl.zip":("5bd93ce604256d00233d4c5bf91865205b44b0d7fc2f52fefe0f68ef0e23e94c",141104),
"20221220_rplN_GFP.zip":("6d91b5e59329712b6920aa1bb78fd06c8f838197be38cf546dbddaace35b1402",65624),
"20221220_rplN_GFP_pl.zip":("7f55aaa61d7025260773c33aa662e13ae2b2266dea300384f9e09a1667dc5ac6",193307),
"20221220_NO_LIGHT_E08_GFP.zip":("d9a00f4ba085b3445448eec94a92c290e6edc66f9b2a1905255879460c4348ee",38345),
"20221220_NO_LIGHT_E08_GFP_pl.zip":("bae990f0288f2898ab3c9525181cbceb3637cf2a08f1a120ec5292fad597f7fb",16967),
"20221220_NO_LIGHT_rplN_GFP.zip":("79fbbe50faef4d417813e81162370848e631d9d4975ab1f696daeacc5e65aeeb",27476),
"20221220_NO_LIGHT_rplN_GFP_pl.zip":("dbb84066aa1754450a98e772b73b6d597683e13d438e51e578db99b50d3da4d6",33593)}
assert AR in EXPECTED
O=pathlib.Path("e505_full_verified")/AR[:-4];O.mkdir(parents=True,exist_ok=True)
E=pathlib.Path("e504_frozen")
def pinned(name,sha):
    f=E/name; b=f.read_bytes();assert hashlib.sha256(b).hexdigest()==sha,(name,len(b));return list(csv.DictReader(io.StringIO(b.decode("utf-8-sig"))))
r0=pinned("E504_pre_exposure_last_frame_rank0_exact_IDs.csv","5f104a4cc5b1c005a8d874887a83c23a4fc92ad96317757078a751ef03d60510")
sel=pinned("E504_238_selected_exact_phase_identity_bridge.csv","f7d36caff67a9b14b655e39f8efd0e5822ea27ef4efb926e7a81bbf05914b342")
assert len(r0)==231 and len(sel)==238
K=lambda r:(r["condition"],r["gl_id"])
old={K(r):r for r in r0}
new=collections.defaultdict(list)
for r in sel:new[K(r)].append(r)
assert set(old)==set(new) and len(new)==231
gl_archive={}
for k,ls in new.items():
    assert len({z["strain"] for z in ls})==1
    strain=ls[0]["strain"]
    gl_archive[k]="20221220_"+("NO_LIGHT_" if k[0].endswith("NO_LIGHT") else "")+strain+".zip"
assert set(gl_archive.values())==set(EXPECTED)
mismatch=[k for k in old if old[k]["GM2_author_selected_same_exact_cell"]=="False"]
assert len(mismatch)==118
matched=[k for k in old if k not in mismatch]
assert len(matched)==113
pilot_controls={"20221220.10.11","20221220.16.27"}
one_per_archive={}
for k in sorted(matched):
    one_per_archive.setdefault(gl_archive[k],k)
targets=[k for k in old if gl_archive[k]==AR and (k in mismatch or k in one_per_archive.values() or k[1] in pilot_controls)]
assert len(set(targets))==len(targets)
targetnames={f"Pos{int(k[1].split('.')[1])}_GL{int(k[1].split('.')[2])}":k for k in targets}
assert len(targetnames)==len(targets)
S=requests.Session();S.headers["User-Agent"]="901-E505-original-S-BIAD2506-MoMA-GM1-GM2-parent-20261116"
BASE="https://www.ebi.ac.uk/biostudies/files/S-BIAD2506/"
a=S.get(BASE+"Statphase_Transition_GLs/deepMoMA.json",timeout=(15,45));a.raise_for_status()
assert hashlib.sha256(a.content).hexdigest()=="f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f"
meta={pathlib.PurePosixPath(v["path"]).name:v for v in json.loads(a.content)}
assert AR in meta
m=meta[AR];N=int(m["size"]);url=BASE+m["path"]
def GET(start,end):
    assert 0<=start<=end<N
    out=bytearray()
    # each request <= 4 MiB, hard 4 MiB reads
    for s in range(start,end+1,4194304):
        e=min(end,s+4194304-1)
        last=None
        for attempt in range(3):
            try:
                r=S.get(url,headers={"Range":f"bytes={s}-{e}"},stream=True,timeout=(15,70))
                if r.status_code!=206:raise RuntimeError(f"HTTP {r.status_code} expected206")
                if r.headers.get("Content-Range")!=f"bytes {s}-{e}/{N}":raise RuntimeError("Content-Range mismatch")
                b=bytearray()
                for chunk in r.iter_content(131072):
                    if chunk:b+=chunk
                    if len(b)>e-s+1:raise RuntimeError("oversize")
                assert len(b)==e-s+1
                out.extend(b);last=None;break
            except Exception as exc:
                last=exc;time.sleep(1+attempt)
        if last:raise last
    assert len(out)==end-start+1
    return bytes(out)
tail=GET(N-131072,N-1)
pos=tail.rfind(b"PK\x05\x06");assert pos>=0
end=struct.unpack_from("<4s4H2IH",tail,pos)
cnt,dlen,dstart=end[4],end[5],end[6]
if cnt==65535 or dlen==0xffffffff or dstart==0xffffffff:
    loc=tail.rfind(b"PK\x06\x07",0,pos);assert loc>=0
    off=struct.unpack_from("<Q",tail,loc+8)[0]
    z=GET(off,off+55);assert z[:4]==b"PK\x06\x06"
    cnt=struct.unpack_from("<Q",z,32)[0]
    dlen=struct.unpack_from("<Q",z,40)[0]
    dstart=struct.unpack_from("<Q",z,48)[0]
assert cnt==EXPECTED[AR][1] and 0<dlen<67108864
data=GET(dstart,dstart+dlen-1)
assert hashlib.sha256(data).hexdigest()==EXPECTED[AR][0]
def z64(extra,raw,compressed,offset):
    if all(a!=0xffffffff for a in [raw,compressed,offset]):return raw,compressed,offset
    i=0
    while i<len(extra):
        tag,l=struct.unpack_from("<HH",extra,i);v=extra[i+4:i+4+l];i+=4+l
        if tag!=1:continue
        j=0
        if raw==0xffffffff:raw=struct.unpack_from("<Q",v,j)[0];j+=8
        if compressed==0xffffffff:compressed=struct.unpack_from("<Q",v,j)[0];j+=8
        if offset==0xffffffff:offset=struct.unpack_from("<Q",v,j)[0];j+=8
        return raw,compressed,offset
    raise ValueError("ZIP64_EXTRA_NOT_FOUND")
entries={};idx=0;total=0
while idx<len(data):
    hd=struct.unpack_from("<4s6H3I5H2I",data,idx)
    if hd[0]!=b"PK\x01\x02":raise ValueError(f"CD signature index{idx}")
    nl,el,cl=hd[10:13]
    nm=data[idx+46:idx+46+nl].decode("utf-8" if hd[3]&2048 else "cp437","replace")
    ex=data[idx+46+nl:idx+46+nl+el]
    mt=re.search(r"/(Pos\d+_GL\d+)/20221220_analysis/export_data__20221220_analysis/CellStats__[^/]+\.csv$",nm)
    if mt and mt.group(1) in targetnames:
        raw,com,loc=z64(ex,hd[9],hd[8],hd[16])
        entries.setdefault(mt.group(1),[]).append({"path":nm,"raw":raw,"compressed":com,"offset":loc,"method":hd[4],"crc":hd[7]})
    idx+=46+nl+el+cl;total+=1
assert idx==dlen and total==cnt
assert all(len(entries.get(g,[]))==1 for g in targetnames),[(g,len(entries.get(g,[]))) for g in targetnames if len(entries.get(g,[]))!=1]
def cellnum(row):
    parts=row["cell"].split(".")
    assert len(parts)==4 and parts[0]=="20221220"
    assert re.fullmatch(r"\d+",parts[-1]) and row["gl_id"]==".".join(parts[:3])
    return int(parts[-1])
def age_chain(target,parent,stop):
    visited=set();chain=[target]
    while True:
        if chain[-1]==stop:return {"reaches_old":True,"ids":chain,"hops":len(chain)-1}
        cur=chain[-1]
        if cur in visited or cur not in parent:return {"reaches_old":False,"ids":chain,"hops":None,"reason":"CYCLE_OR_MISSING"}
        visited.add(cur)
        anc=parent[cur]
        if anc<0:return {"reaches_old":False,"ids":chain,"hops":None,"reason":"NO_PARENT"}
        chain.append(anc)
        if len(chain)>5000:return {"reaches_old":False,"ids":chain,"hops":None,"reason":"CHAIN_DEPTH_LIMIT"}
allres=[];src_receipts=[];total_uncompressed=0
for glname,k in sorted(targetnames.items()):
    ent=entries[glname][0]
    assert 0<ent["compressed"]<12582912 and 0<ent["raw"]<25000000
    off=ent["offset"]
    lh=struct.unpack("<4s5H3I2H",GET(off,off+29));assert lh[0]==b"PK\x03\x04"
    nl,el=lh[-2:]
    head=GET(off+30,off+30+nl+el-1)
    assert head[:nl].decode("utf-8" if lh[2]&2048 else "cp437","replace")==ent["path"]
    body=GET(off+30+nl+el,off+30+nl+el+ent["compressed"]-1)
    if ent["method"]==8:raw=zlib.decompress(body,-15)
    elif ent["method"]==0:raw=body
    else:raise RuntimeError("ZIP method unsupported")
    assert len(raw)==ent["raw"] and zlib.crc32(raw)==ent["crc"]
    total_uncompressed+=len(raw)
    assert total_uncompressed<250000000
    text=raw.decode("utf-8-sig")
    txt="\n".join(l for l in text.splitlines() if not l.startswith("#") and l.strip())
    reader=csv.DictReader(io.StringIO(txt))
    req={"lane_ID","cell_ID","frame","cell_rank","genealogy","type_of_end","parent_ID"}
    assert req.issubset(set(reader.fieldnames or [])),reader.fieldnames
    record=collections.defaultdict(lambda:{"parent_ids":set(),"genealogies":set(),"end_types":set(),"frames":[],"rank_at239":[],"rank_after555":0})
    parsed=0
    for row in reader:
        parsed+=1
        cid=int(row["cell_ID"]);pa=int(row["parent_ID"]) if row["parent_ID"].strip() else -1
        frame=int(row["frame"]);rank=int(row["cell_rank"])
        v=record[cid]
        v["parent_ids"].add(pa);v["genealogies"].add(row["genealogy"])
        if row["type_of_end"].strip():v["end_types"].add(row["type_of_end"].strip())
        v["frames"].append(frame)
        if frame==239:v["rank_at239"].append(rank)
        if frame>=555 and rank==0:v["rank_after555"]+=1
        assert row["lane_ID"].lower().replace("_","")==glname.lower().replace("_",""),row["lane_ID"]
    invalid_parent={cid:sorted(x["parent_ids"]) for cid,x in record.items() if len(x["parent_ids"])!=1}
    invalid_genealogy={cid:sorted(x["genealogies"]) for cid,x in record.items() if len(x["genealogies"])!=1}
    parent={cid:next(iter(v["parent_ids"])) for cid,v in record.items() if len(v["parent_ids"])==1}
    gene={cid:next(iter(v["genealogies"])) for cid,v in record.items() if len(v["genealogies"])==1}
    oldid=cellnum(old[k]);newids=[cellnum(x) for x in new[k]]
    oldin=oldid in record
    o239=(0 in record[oldid]["rank_at239"]) if oldin else False
    # trace explicit MoMA parent IDs to original pre-starvation mother ID
    mapped=[]
    for r in new[k]:
        cid=cellnum(r);link=age_chain(cid,parent,oldid)
        link["mother_cell_ID"]=cid
        link["fate"]=r["fate"]
        link["time_start_s"]=float(r["time_start"])
        link["post_rank0"]=record.get(cid,{}).get("rank_after555",0)>0
        link["genealogical_string_source"]=gene.get(cid)
        link["old_genealogical_string_source"]=gene.get(oldid)
        edgeChecks=[]
        if link["reaches_old"] and len(link["ids"])>1:
            for ch,par in zip(link["ids"][:-1],link["ids"][1:]):
                if ch not in record or par not in record:
                    edgeChecks.append({"child":ch,"parent":par,"valid":False,"reason":"NOT_IN_SOURCE"});continue
                pv=record[par];cv=record[ch]
                after=(min(cv["frames"])==max(pv["frames"])+1)
                parent_div="div" in pv["end_types"]
                prefix=bool(gene.get(par)) and bool(gene.get(ch)) and gene[ch].startswith(gene[par])
                edgeChecks.append({"child":ch,"parent":par,"valid":bool(after and parent_div and prefix),
                                   "frame_contiguous":after,"parent_end_div":parent_div,"genealogy_prefix":prefix,
                                   "parent_last_frame":max(pv["frames"]),"child_first_frame":min(cv["frames"])})
        link["link_edges"]=edgeChecks
        link["all_edges_pass"]=bool(link["reaches_old"] and all(x["valid"] for x in edgeChecks))
        mapped.append(link)
    rec={
      "archive":AR,"condition":k[0],"gl_id":k[1],"archive_cellstats_path":ent["path"],
      "source_raw_sha256":hashlib.sha256(raw).hexdigest(),"source_raw_crc32":format(ent["crc"],"08x"),
      "source_raw_bytes":len(raw),"source_compressed_bytes":ent["compressed"],"source_zip_offset":ent["offset"],
      "number_of_cell_instances":parsed,"number_of_cell_IDs":len(record),
      "source_parent_conflicting_IDs":len(invalid_parent),"source_genealogy_conflicting_IDs":len(invalid_genealogy),
      "pre_mother_ID":oldid,"pre_mother_in_original_stats":oldin,"pre_mother_rank0_at_exact_frame239":o239,
      "GM2_selected_n":len(mapped),"GM2_selected_IDs":newids,"pre_mother_identical_selected":oldid in newids,
      "GM2_selected":[{"id":v["mother_cell_ID"],"fate":v["fate"],"descendant_of_pre_mother":v["reaches_old"],
         "hops":v["hops"],"edges_valid":v["all_edges_pass"],"post_rank0":v["post_rank0"]} for v in mapped],
      "any_parent_chain_to_pre_mother":any(v["reaches_old"] for v in mapped),
      "any_strict_direct_parent":any(v["reaches_old"] and v["hops"]==1 and v["all_edges_pass"] for v in mapped),
      "all_parent_chains_from_pre_mother":all(v["reaches_old"] for v in mapped),
      "every_linked_chain_valid":all(v["all_edges_pass"] for v in mapped if v["reaches_old"]),
      "parent_source_type":"MoMA algorithmal lineage; not independent physical image revalidation",
      "details":mapped
    }
    allres.append(rec)
    src_receipts.append({a:rec[a] for a in ("archive","condition","gl_id","source_raw_sha256","source_raw_crc32","source_raw_bytes","source_compressed_bytes","source_zip_offset","number_of_cell_instances","number_of_cell_IDs")})
    print("E505_GL_TRUE_SOURCE",json.dumps({a:rec[a] for a in ["gl_id","pre_mother_ID","GM2_selected_IDs","any_parent_chain_to_pre_mother","any_strict_direct_parent","pre_mother_rank0_at_exact_frame239"]},ensure_ascii=False),flush=True)
# known true-source pilot positive controls
for entry in allres:
  if entry["gl_id"]=="20221220.10.16":
    assert entry["pre_mother_ID"]==75 and 83 in entry["GM2_selected_IDs"] and entry["any_strict_direct_parent"]
  if entry["gl_id"]=="20221220.16.14":
    assert entry["pre_mother_ID"]==89 and 99 in entry["GM2_selected_IDs"] and entry["any_strict_direct_parent"]
  if entry["gl_id"]=="20221220.10.11":assert entry["pre_mother_ID"] in entry["GM2_selected_IDs"]
  if entry["gl_id"]=="20221220.16.27":assert entry["pre_mother_ID"] in entry["GM2_selected_IDs"]
assert len(allres)==len(targets)
summ={"study":"S-BIAD2506 2022-12-20 source-authenticated MoMA data","stage":"901 E505",
      "archive":AR,"manifest_sha":"f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f",
      "zip_cdir_sha":EXPECTED[AR][0],"n_entries":EXPECTED[AR][1],"n_target_GL":len(targets),
      "n_total_mismatched_GL":sum(k in mismatch for k in targets),
      "n_total_matched_QC_GL":sum(k in matched for k in targets),
      "all_pre_ID_original_stats_present":all(x["pre_mother_in_original_stats"] for x in allres),
      "all_pre_ID_rank0_frame239":all(x["pre_mother_rank0_at_exact_frame239"] for x in allres),
      "all_selected_ID_in_original_stats":all(v["post_rank0"] for x in allres for v in x["GM2_selected"]),
      "any_descendant_GL":sum(x["any_parent_chain_to_pre_mother"] for x in allres if not x["pre_mother_identical_selected"]),
      "direct_valid_GL":sum(x["any_strict_direct_parent"] for x in allres if not x["pre_mother_identical_selected"]),
      "missing_genealogy_or_parent_conflicts":sum(x["source_parent_conflicting_IDs"]+x["source_genealogy_conflicting_IDs"] for x in allres),
      "no_full_image_download":True,
      "timestamp_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime())}
(O/"E505_full_original_CSV_parent_chain_per_GL.json").write_text(json.dumps(allres,ensure_ascii=False,indent=2))
(O/"E505_source_authentication_receipts.csv").write_text("\n".join([",".join(src_receipts[0].keys())]+[",".join(str(x) for x in v.values()) for v in src_receipts])+"\n" if src_receipts else "no_targets\n")
(O/"E505_per_archive_summary.json").write_text(json.dumps(summ,ensure_ascii=False,indent=2))
print("E505_ARCHIVE_SUMMARY",json.dumps(summ,ensure_ascii=False),flush=True)
