#!/usr/bin/env python3
"""901 E510: authenticated original preprocessed 2D phase-contrast TIFF pixel intake.

At most one complete archived TIFF is downloaded per job from official BIA
S-BIAD2506 via strict bounded HTTP206 ZIP64 member offsets. All frames are
selected from frozen E509 events *before* pixels are examined. Blinded mosaic
PNGs and exact float32 phase pixel NPZs are preserved for independent review.
MoMA labels are NOT independent pixel ground truth. No calibration to photons.
"""
import argparse, pathlib, hashlib, json, struct, zlib, time, io, collections
import requests, numpy as np, tifffile
from PIL import Image, ImageDraw

P=argparse.ArgumentParser()
P.add_argument("--gl",required=True,choices=["Pos5_GL35","Pos12_GL22"])
args=P.parse_args()
GL=args.gl
ARCHIVES={
 "Pos5_GL35":{
    "archive":"20221220_rplN_GFP_pl.zip",
    "cd_sha256":"7f55aaa61d7025260773c33aa662e13ae2b2266dea300384f9e09a1667dc5ac6",
    "count":193307,"expected_disputed":9,"expected_matched":4,
 },
 "Pos12_GL22":{
    "archive":"20221220_E08_GFP_pl.zip",
    "cd_sha256":"5bd93ce604256d00233d4c5bf91865205b44b0d7fc2f52fefe0f68ef0e23e94c",
    "count":141104,"expected_disputed":3,"expected_matched":9,
 }
}
CFG=ARCHIVES[GL];ARCH=CFG["archive"]
OUT=pathlib.Path("e510_true_2d_phase_pixels")/GL;OUT.mkdir(parents=True,exist_ok=True)
BLIND=OUT/"blinded_cards";BLIND.mkdir(exist_ok=True)

# Freeze event list from E507/E508 source algorithm outputs (frozen methods).
base=pathlib.Path("e507_alternate_QA")/ARCH[:-4]
group=json.loads((base/"E507_alternate_algorithm_frame_geometry_results.json").read_text())
group=[v for v in group if v["GL"]==GL]
assert len(group)==1
group=group[0]
assert group["matched_geometry_divisions"]==CFG["expected_matched"]
p=pathlib.Path("e508_source_conflicts")/ARCH[:-4]
ev=json.loads((p/"E508_event_level_version_conflicts.json").read_text())
ev=[r for r in ev if r["GL"]==GL]
assert len(ev)==CFG["expected_disputed"]
matched=[{"frame":int(r["frame"]),"parent":int(r["main_parent"])}
         for r in group["pairs"]]
assert len({(x["frame"],x["parent"]) for x in matched+ev})==len(matched)+len(ev)
controls=[]
all_divisions=[int(x["frame"]) for x in matched+ev]
used=set(all_divisions)
# Frozen software-non-dividing controls: the earliest eligible frame in the
# [event+10,event+60] window, at least 6 frames from ANY source DIV. This
# group is not certified negative truth; it tests visibly ordinary time windows.
for r in sorted(ev,key=lambda v:(v["frame"],v["parent"])):
    f=int(r["frame"]);choice=None
    for c in range(f+10,min(f+61,700)):
        if c in used or min(abs(c-t) for t in all_divisions)<=6:continue
        choice=c;break
    assert choice is not None
    controls.append({"frame":choice,"parent":None,"anchor_disputed_frame":f})
    used.add(choice)

records=([{"frame":int(x["frame"]),"parent":int(x["parent"]),"class":"disputed"} for x in ev]
    +[{"frame":int(x["frame"]),"parent":int(x["parent"]),"class":"software_agreement"} for x in matched]
    +[{"frame":int(x["frame"]),"parent":None,"class":"software_nondivision"} for x in controls])
# Neutral order, selected without viewing image data. Sealed key saved separately.
for r in records:
    r["sort_digest"]=hashlib.sha256(("901E510-frozen-"+GL+"-"+str(r["frame"])+":"+str(r["parent"])).encode()).hexdigest()
records.sort(key=lambda r:r["sort_digest"])
for i,r in enumerate(records,1):
    r["blind_case_id"]=f"B{i:03d}"
assert len(records)==2*CFG["expected_disputed"]+CFG["expected_matched"]
(OUT/"E510_sealed_event_key_do_not_open_for_blind_review.json").write_text(
    json.dumps(records,ensure_ascii=False,indent=2),encoding="utf8"
)
print("E510_FROZEN_CASE_COUNTS",GL,json.dumps(dict(collections.Counter(r["class"] for r in records))),flush=True)
print("E510_BLIND_RULE_PRE_FROZEN","time offsets [-3,-2,-1,0,1,2,3,5,8]; no software bbox overlays; same scale across GL",flush=True)

# Authenticate original file listing.
BASE="https://www.ebi.ac.uk/biostudies/files/S-BIAD2506/"
S=requests.Session();S.headers["User-Agent"]="901-E510-original-preprocessed-2D-phase-source-CRC-SHA-20261116"
man=S.get(BASE+"Statphase_Transition_GLs/deepMoMA.json",timeout=(15,55))
man.raise_for_status()
MAN_SHA=hashlib.sha256(man.content).hexdigest()
assert MAN_SHA=="f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f"
index={pathlib.PurePosixPath(i["path"]).name:i for i in json.loads(man.content)}
entry=index[ARCH];SIZE=int(entry["size"]);URL=BASE+entry["path"]
def get_part(lo,hi):
    assert 0<=lo<=hi<SIZE and hi-lo+1<=4194304
    error=None
    for attempt in range(4):
        try:
            rr=S.get(URL,headers={"Range":f"bytes={lo}-{hi}"},timeout=(18,110),stream=True)
            if rr.status_code!=206 or rr.headers.get("Content-Range")!=f"bytes {lo}-{hi}/{SIZE}":
                raise RuntimeError("source not strict HTTP206")
            buf=bytearray()
            for chunk in rr.iter_content(131072):
                if chunk:buf.extend(chunk)
                if len(buf)>hi-lo+1:raise RuntimeError("source bytes overflow")
            if len(buf)!=hi-lo+1:raise RuntimeError("source truncated")
            return bytes(buf)
        except Exception as e:
            error=e; time.sleep(1+attempt)
    raise error
def get_range(lo,hi):
    return b"".join(get_part(q,min(q+4194303,hi)) for q in range(lo,hi+1,4194304))
def decode_zip64(ex,raw,compressed,offset):
    if raw!=0xffffffff and compressed!=0xffffffff and offset!=0xffffffff:
        return raw,compressed,offset
    i=0
    while i<len(ex):
        tag,size=struct.unpack_from("<HH",ex,i)
        d=ex[i+4:i+4+size];i+=4+size
        if tag!=1:continue
        j=0
        if raw==0xffffffff:raw=struct.unpack_from("<Q",d,j)[0];j+=8
        if compressed==0xffffffff:compressed=struct.unpack_from("<Q",d,j)[0];j+=8
        if offset==0xffffffff:offset=struct.unpack_from("<Q",d,j)[0];j+=8
        return raw,compressed,offset
    raise RuntimeError("missing ZIP64 extra")
tail=get_part(SIZE-131072,SIZE-1)
pos=tail.rfind(bytes.fromhex("504b0506"))
assert pos>=0
h=struct.unpack_from("<4s4H2IH",tail,pos)
count,cd_size,cd_offset=h[4],h[5],h[6]
if count==65535 or cd_size==0xffffffff or cd_offset==0xffffffff:
    pos64=tail.rfind(bytes.fromhex("504b0607"),0,pos)
    assert pos64>=0
    off=struct.unpack_from("<Q",tail,pos64+8)[0]
    v=get_part(off,off+55)
    assert v[:4]==bytes.fromhex("504b0606")
    count=struct.unpack_from("<Q",v,32)[0]
    cd_size=struct.unpack_from("<Q",v,40)[0]
    cd_offset=struct.unpack_from("<Q",v,48)[0]
assert count==CFG["count"] and 0<cd_size<67108864
cd=get_range(cd_offset,cd_offset+cd_size-1)
assert hashlib.sha256(cd).hexdigest()==CFG["cd_sha256"]
wanted=f"{ARCH[:-4]}/{GL}/20221220_glu_spcm_1_MMStack_{GL}.tif"
found=[];off=0;n=0
while off<len(cd):
    h=struct.unpack_from("<4s6H3I5H2I",cd,off)
    assert h[0]==bytes.fromhex("504b0102")
    nl,el,cl=h[10:13]
    name=cd[off+46:off+46+nl].decode("utf8" if h[3]&2048 else "cp437","replace")
    if name==wanted:
        raw,compressed,loc=decode_zip64(cd[off+46+nl:off+46+nl+el],h[9],h[8],h[16])
        found.append({"path":name,"raw":raw,"compressed":compressed,"offset":loc,"crc":h[7],"method":h[4]})
    off+=46+nl+el+cl;n+=1
assert n==count and off==cd_size and len(found)==1,(n,count,len(found),wanted)
z=found[0]
assert z["method"]==8 and 0<z["raw"]<650_000_000 and 0<z["compressed"]<550_000_000,z
loc=z["offset"]
lhead=get_part(loc,loc+29)
lh=struct.unpack_from("<4s5H3I2H",lhead)
assert lh[0]==bytes.fromhex("504b0304")
nl,el=lh[-2:]
assert get_part(loc+30,loc+30+nl+el-1)[:nl].decode("utf8" if lh[2]&2048 else "cp437","replace")==wanted
start=loc+30+nl+el;end=start+z["compressed"]-1
# Stream compressed member into an ephemeral runner file; never upload full TIFF.
temp=pathlib.Path("/tmp/901_e510_"+GL+"_official_preprocessed_GL.tif")
sha=hashlib.sha256();crc=0;raw_bytes=0;unpack=zlib.decompressobj(-15)
with temp.open("wb") as f:
    for left in range(start,end+1,4194304):
        block=get_part(left,min(end,left+4194303))
        b=unpack.decompress(block)
        if b:f.write(b);sha.update(b);crc=zlib.crc32(b,crc);raw_bytes+=len(b)
    b=unpack.flush()
    if b:f.write(b);sha.update(b);crc=zlib.crc32(b,crc);raw_bytes+=len(b)
assert unpack.eof and raw_bytes==z["raw"] and crc==z["crc"],("CRC_FAIL",raw_bytes,z)
print("E510_ORIGINAL_COMPLETE_TIFF_SOURCE_AUTH_PASS",GL,raw_bytes,sha.hexdigest(),f"{crc:08x}",flush=True)
with tifffile.TiffFile(temp) as tif:
    series=tif.series[0]
    axes=list(series.axes)
    arr=series.asarray()
    # source original E506 reference: 1036 time steps, 3 channels, Y x 106
    if "T" in axes:tix=axes.index("T")
    else:
        tix=next(k for k,size in enumerate(arr.shape) if size==1036)
        axes[tix]="T"
    arr=np.moveaxis(arr,tix,0);axes.insert(0,axes.pop(tix))
    if "C" in axes:
        cix=axes.index("C");arr=np.moveaxis(arr,cix,1);axes.insert(1,axes.pop(cix))
    else:
        assert arr.ndim==3;arr=arr[:,None];axes.insert(1,"C")
    assert arr.ndim==4 and arr.shape[0]==1036 and arr.shape[1]>=3 and arr.shape[2]>=200 and arr.shape[3]==106,(axes,arr.shape)
    phase=arr[:,0,:,:]
    assert phase.dtype.kind=="f" and np.isfinite(phase[::17]).all(),(phase.dtype,axes)
    # Frame number in source MoMA is 1-based: source TIFF arr[frame - 1].
    # Same mapping is confirmed independently by E506 original TIFF/kymo source audit.
    # A frame selected by source E509 should be visible directly without bbox transformations.
    low,high=np.quantile(phase[::23,::4,::3].astype("float64"),[.01,.99])
    assert high>low

    offsets=[-3,-2,-1,0,1,2,3,5,8]
    summary=[]
    for r in records:
        center=int(r["frame"])
        frames=[center+d for d in offsets]
        assert min(frames)>0 and max(frames)<=1036
        source=np.stack([phase[f-1].copy() for f in frames])   # exact float source data
        assert source.dtype==np.float32 and source.shape[1:]==tuple(phase.shape[1:])
        cardid=r["blind_case_id"]
        np.savez_compressed(OUT/(cardid+"_exact_float32_phase_frames.npz"),
                            frames=np.asarray(frames,dtype="int32"),
                            offsets=np.asarray(offsets,dtype="int32"),
                            phase_pixels=source)
        # Blind visualization: same source-level 1/99 percentile grayscale scale
        # across every candidate GL frame. No overlays, labels or ROIs from software.
        ylo,yhi=80,min(source.shape[1],255)
        tiles=[]
        for i in range(len(offsets)):
            img=np.uint8(np.clip((source[i,ylo:yhi]-low)/(high-low),0,1)*255)
            # upscale with nearest-neighbor for reproducible pixel inspection
            pil=Image.fromarray(img,mode="L").resize((318,(yhi-ylo)*3),Image.Resampling.NEAREST).convert("RGB")
            tiles.append(pil)
        canvas=Image.new("RGB",(318*len(tiles),(yhi-ylo)*3+45),(235,235,235))
        draw=ImageDraw.Draw(canvas)
        for j,tile in enumerate(tiles):
            canvas.paste(tile,(j*318,45))
            draw.text((j*318+8,8),f"{cardid}  offset={offsets[j]:+d}",fill=(20,20,20))
        canvas.save(BLIND/(cardid+"_blind_preprocessed_PhC.png"))
        # Summary independent of software bbox, no event classification conclusions:
        raw_t=np.asarray(source[3],dtype="float64")
        raw_tnext=np.asarray(source[4],dtype="float64")
        frame_change=float(np.mean(np.abs(raw_tnext-raw_t)))
        summary.append({"blind_case_id":cardid,"pixel_frames_saved":len(frames),
                        "phase_shape_yx":list(raw_t.shape),"finite":bool(np.isfinite(source).all()),
                        "mean_abs_delta_t_to_tplus1":frame_change,
                        "pixel_float32_sha256":hashlib.sha256(source.tobytes()).hexdigest(),
                        "visualization_low_global":float(low),"visualization_high_global":float(high),
                        "visualization_crop_rows":[ylo,yhi],"visualization_nearest_neighbor_3x":True})
        print("E510_AUTHENTIC_PIXEL_CASE",cardid,"original_frame",center,"shape",source.shape,flush=True)
    (OUT/"E510_blind_per_case_pixel_checks.json").write_text(json.dumps(summary,indent=2),encoding="utf8")
    receipt={"evidence_stage":"SOURCE_2D_TIFF_PIXEL_INTAKE_ONLY",
             "source":"BioImage Archive S-BIAD2506",
             "archive":ARCH,"GL":GL,"verified_study_manifest_sha256":MAN_SHA,
             "verified_zip_cdir_sha256":CFG["cd_sha256"],
             "verified_zip_member_count":n,"member":z,
             "original_uncompressed_tiff_sha256":sha.hexdigest(),
             "original_uncompressed_tiff_bytes":raw_bytes,
             "original_member_crc32":f"{crc:08x}",
             "original_member_crc_pass":True,
             "axes_after_reordering":axes,"full_image_shape":list(arr.shape),
             "original_pixel_dtype":str(source.dtype),"same_source_preprocessed_not_raw_camera":True,
             "no_moma_bbox_used_to_select_image_roi":True,
             "not_independent_pixel_segmentation_or_biological_ground_truth":True,
             "case_count":len(records),
             "no_full_image_upload":True,
             "case_input_source":"E509/E508 event register; same fixed event source ledger",
             "time_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime())}
    (OUT/"E510_original_full_preprocessed_2D_TIFF_source_receipt.json").write_text(
        json.dumps(receipt,ensure_ascii=False,indent=2),encoding="utf8")
    print("E510_ALL_2D_PIXEL_CASES_PRESERVED",GL,len(records),"CRC_SHA_PASS",flush=True)
temp.unlink()
