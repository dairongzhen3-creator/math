#!/usr/bin/env python3
"""E529: read only two original MoMA probability pixels from official S-BIAD2506 ZIP64 members.
No original MPS, microscopy truth, physical information rate or molecular measurement is inferred.
"""
import datetime, hashlib, io, json, math, pathlib, struct, sys, time, traceback, zlib
import requests
import numpy as np
import tifffile

BASE = "https://www.ebi.ac.uk/biostudies/files/S-BIAD2506/"
MANIFEST_SHA = "f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f"
OUT = pathlib.Path("E529_outputs")
OUT.mkdir(exist_ok=True)
CFG = [
    dict(gl="Pos5_GL35", archive="20221220_rplN_GFP_pl.zip",
         path="20221220_rplN_GFP_pl/Pos5_GL35/cost_function_analysis_1/track_data__cost_function_analysis_1/probability_maps.tif",
         loc=16898524769, packed=23663755, unpacked=131932501,
         crc="6a9ae969", sha="2d06799a61c6b37c1494360d745024168ff951e4bfab0f8a0b6bf9d1db7ba9c8",
         frame=960, x=55, y=160, parent_component_delta="0.00000522",
         original_5_action_delta="0.00000519", e527_conditional_p=0.518261873145),
    dict(gl="Pos12_GL22", archive="20221220_E08_GFP_pl.zip",
         path="20221220_E08_GFP_pl/Pos12_GL22/cost_function_analysis_1/track_data__cost_function_analysis_1/probability_maps.tif",
         loc=15926971746, packed=23461517, unpacked=131932501,
         crc="10595976", sha="bea1570256e7043a8af2a0430865b64771ca940d379612af7efdd5c6e15aa7bd",
         frame=118, x=50, y=197, parent_component_delta="0.00001563",
         original_5_action_delta="0.00001565", e527_conditional_p=0.554487757202)
]
S = requests.Session()
S.headers.update({"User-Agent":"901-E529-source-only-reproducibility/20261117",
                  "Accept":"application/octet-stream","Accept-Encoding":"identity"})
SHA=lambda b:hashlib.sha256(b).hexdigest()
def jwrite(name, obj):
    (OUT / name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")

def exact_range(url, start, stop, total):
    assert 0 <= start <= stop < total
    data=bytearray()
    for a in range(start,stop+1,4*1024*1024):
        b=min(stop,a+4*1024*1024-1)
        failure=None
        for attempt in range(3):
            try:
                with S.get(url,headers={"Range":f"bytes={a}-{b}"},
                           timeout=(15,95),stream=True) as r:
                    if r.status_code!=206:
                        raise RuntimeError(f"REQUIRED_HTTP_206, received {r.status_code} for {a}-{b}")
                    wanted=f"bytes {a}-{b}/{total}"
                    if r.headers.get("Content-Range")!=wanted:
                        raise RuntimeError(f"WRONG_CONTENT_RANGE {r.headers.get('Content-Range')} expected {wanted}")
                    chunk=bytearray()
                    for part in r.iter_content(chunk_size=131072):
                        if part: chunk.extend(part)
                        if len(chunk)>b-a+1:raise ValueError("range overflow")
                    if len(chunk)!=b-a+1:raise ValueError(f"range truncated: {len(chunk)}")
                data.extend(chunk)
                break
            except Exception as ex:
                failure=ex
                time.sleep(attempt+1)
        else:raise failure
    assert len(data)==stop-start+1
    return bytes(data)

def read_tiff_member(entry, manifest):
    item=manifest[entry["archive"]]
    N=int(item["size"])
    assert N > entry["loc"]+entry["packed"] and N < 50_000_000_000
    url=BASE+item["path"]
    loc=entry["loc"]
    header=exact_range(url,loc,loc+29,N)
    h=struct.unpack("<4s5H3I2H",header)
    assert h[0]==b"PK\x03\x04", "ZIP local header invalid"
    method=h[3]
    assert method==8, ("compression method not DEFLATE",method)
    filename_bytes=exact_range(url,loc+30,loc+30+h[-2]+h[-1]-1,N)
    filename=filename_bytes[:h[-2]].decode("utf-8" if (h[2]&2048) else "cp437")
    assert filename==entry["path"], ("wrong source member path",filename)
    payload_offset=loc+30+h[-2]+h[-1]
    compressed=exact_range(url,payload_offset,payload_offset+entry["packed"]-1,N)
    source=zlib.decompress(compressed,-15)
    assert len(source)==entry["unpacked"], ("TIFF byte count",len(source))
    assert format(zlib.crc32(source)&0xffffffff,"08x")==entry["crc"], "original ZIP central CRC mismatch"
    assert SHA(source)==entry["sha"], "original E512 full TIFF SHA mismatch"
    return source,dict(archive_url=url,archive_size=N,original_member_path=filename,
                       source_local_offset=loc,compressed_byte_range=[payload_offset,payload_offset+entry["packed"]-1],
                       compressed_bytes=len(compressed),compressed_sha256=SHA(compressed),
                       full_tiff_bytes=len(source),full_tiff_sha256=SHA(source),
                       source_zip_member_crc32=entry["crc"],
                       strict_partial_content=True)

def main():
    req=S.get(BASE+"Statphase_Transition_GLs/deepMoMA.json",timeout=(15,75))
    req.raise_for_status()
    assert SHA(req.content)==MANIFEST_SHA, "official manifest source identity changed"
    manifest={pathlib.PurePosixPath(row["path"]).name:row for row in json.loads(req.content)}
    results=[]
    for c in CFG:
        data,receipt=read_tiff_member(c,manifest)
        with tifffile.TiffFile(io.BytesIO(data)) as tif:
            assert len(tif.pages)==1036
            shape=list(tif.series[0].shape)
            assert tuple(shape)==(1036,300,106),shape
            pp={}
            for offset in [-2,-1,0,1,2]:
                t=c["frame"]+offset
                a=tif.pages[t].asarray()
                assert a.shape==(300,106) and a.dtype==np.dtype("float32")
                val=float(a[c["y"],c["x"]])
                assert math.isfinite(val) and 0 <= val <= 1
                pp[str(offset)]={"original_tiff_array_index_zero_based":t,
                                 "source_frame_number_one_based":t+1,
                                 "probability_float32":val,
                                 "raw_float32_hex":a[c["y"],c["x"]].tobytes().hex()}
        p=pp["0"]["probability_float32"]
        pc=min(max(p,1e-7),0.9999999)
        kappa=0.2/2800
        actual=kappa*math.log(pc/(1-pc))
        node_delta=float(c["parent_component_delta"])
        predicted_p=1/(1+math.exp(-node_delta/kappa))
        error=actual-node_delta
        result=dict(gl=c["gl"],component_frame_zero_based=c["frame"],
          requested_tiff_array_index_zero_based=c["frame"],source_pixel_xy=[c["x"],c["y"]],
          original_probability=p, model_component_cost_from_original_pixel=actual,
          original_component_price_difference=node_delta,
          original_five_action_cost_difference=c["original_5_action_delta"],
          component_cost_residual=error, absolute_component_cost_residual=abs(error),
          conditional_inverted_probability_E527=c["e527_conditional_p"],
          source_pixel_minus_E527_conditional_p=p-c["e527_conditional_p"],
          inverted_probability_from_rounded_component_cost=predicted_p,
          source_5_frame_controls=pp, provenance=receipt,
          e528_formula="kappa * log(p/(1-p)); kappa=0.2/2800; cost source is separate archive",
          source_agrees_to_e528_tolerance_1e_minus_7=bool(abs(error)<1e-7))
        results.append(result)
        jwrite(f"E529_{c['gl']}_official_one_pixel.json",result)
        print("E529_SOURCE_PIXEL",json.dumps({k:result[k] for k in
            ("gl","source_pixel_xy","original_probability","original_component_price_difference",
             "model_component_cost_from_original_pixel","absolute_component_cost_residual",
             "source_agrees_to_e528_tolerance_1e_minus_7")}),flush=True)
    summary=dict(project=901,stage="E529",timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source="EMBL-EBI BioStudies public S-BIAD2506",
        official_manifest_sha256=MANIFEST_SHA,exact_tiff_full_sha_crc_for_both=True,
        cases=results, passed_e528_component_cost_tolerance=sum(x["source_agrees_to_e528_tolerance_1e_minus_7"] for x in results),
        n_cases=2, scope="original software model probability and component price only",
        2022_binary_feature_flag_independently_verified=False,
        independent_biological_division_ground_truth="UNRESOLVED",
        physical_bits_per_second_per_joule="NOT_IDENTIFIED",
        original_full_gurobi_mps="UNAVAILABLE")
    jwrite("E529_original_two_pixel_official_verification.json",summary)
    print("E529_FINAL",json.dumps({"n":summary["n_cases"],"pass":summary["passed_e528_component_cost_tolerance"]}),flush=True)
if __name__=="__main__":
    try:main()
    except Exception as exc:
        jwrite("E529_source_failure_receipt.json",
          {"error_type":type(exc).__name__,"error":str(exc),"timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
           "trace":traceback.format_exc()[-2500:],"stage":"E529","verified":False})
        raise
