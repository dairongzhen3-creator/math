#!/usr/bin/env python3
"""901 E554: narrow, source-authenticated original MoMA curation and microscopy intake.
Public BioImage Archive S-BIAD2506, original 20221220 ZIP64 only.
This does not infer biological truth or historical human operator identity.
"""
import argparse, hashlib, importlib.util, json, os, re, shutil, sys, tempfile, traceback, zipfile
from pathlib import Path
import numpy as np
import tifffile
from PIL import Image

E545_SOURCE = Path(__file__).resolve().parents[1] / "E545_20261118" / "E545_official_HTTP206_ZIP64_model_extractor.py"
E545_GIT_BLOB_SHA = "74a157ca970935abec108c957e2f281e8c33e5eb"
EXPECTED_CURATIONS = {
    ("Pos5_GL35","main"):"7b8a49d1b99996e3074f38cff21cad284c471b2879e3fe7419f7e7b5257c1033",
    ("Pos5_GL35","cost"):"00d590779b499d5e4fa1b5053f140b339531dbc84400fd84997927eb7eb65253",
    ("Pos12_GL22","main"):"ba09876bdf2951b1d462725a987969af5d1768d147966fe39e415608dad66732",
    ("Pos12_GL22","cost"):"bda41445be2a901354bea2442180612bd3d073d90e90629093aaa80f66e4cdc9",
}
ROIS = {
 "Pos5_GL35":{"focus":960, "anchors":[200,294,630],"phase_compressed":348447146,"prob_compressed":23663755},
 "Pos12_GL22":{"focus":118, "anchors":[198,240,616],"phase_compressed":336844734,"prob_compressed":23461517}
}
def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4*1024*1024),b""):h.update(b)
    return h.hexdigest()
def save_json(path,x):
    Path(path).write_text(json.dumps(x,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
def normalize_preview(arr):
    a=np.asarray(arr,dtype=np.float64)
    good=a[np.isfinite(a)]
    if not len(good):return Image.fromarray(np.zeros(a.shape[-2:],dtype=np.uint8))
    lo,hi=np.percentile(good,[1,99])
    if hi<=lo:hi=lo+1
    return Image.fromarray(np.clip((a-lo)/(hi-lo)*255,0,255).astype(np.uint8))
def extract_file(z,info,out,source,small_limit=None):
    assert info.file_size <= (small_limit or 450_000_000),("size cap",info.filename,info.file_size)
    h=hashlib.sha256();total=0
    with z.open(info,"r") as inp,Path(out).open("wb") as fd:
        while True:
            b=inp.read(1024*1024)
            if not b:break
            fd.write(b);h.update(b);total+=len(b)
            assert total<=info.file_size
    assert total==info.file_size
    rec={"member_path":info.filename,"source_archive":source,"compressed_bytes":info.compress_size,
         "uncompressed_bytes":info.file_size,"CRC32":f"{info.CRC:08x}",
         "extracted_sha256":h.hexdigest()}
    return rec
def frames_tif(path,gl,role,output):
    config=ROIS[gl]
    with tifffile.TiffFile(str(path)) as tif:
        pages=tif.pages
        series=tif.series[0]
        focus=config["focus"]
        frames=set(range(focus-5,focus+6))
        for a in config["anchors"]:frames.update([a-1,a,a+1])
        if role=="phase_contrast":
            # Verified by E510 official TIFF receipt: 1036 times, 3 channels, 300 x 106.
            # MoMA preprocessing uses channel index 0 for phase contrast.
            assert len(series.shape)==4 and tuple(series.shape[:2])==(1036,3),(gl,role,series.shape)
            assert series.shape[-2:]==(300,106),(gl,series.shape)
            assert len(pages)==1036*3,(gl,role,len(pages),series.shape)
            n=1036
            indices=sorted(i for i in frames if 0<=i<n)
            arrays=[np.asarray(pages[3*i].asarray()) for i in indices]
            page_map={str(i):int(3*i) for i in indices}
            channel=0
        else:
            n=len(pages)
            assert n==1036,(gl,role,n,series.shape)
            indices=sorted(i for i in frames if 0<=i<n)
            arrays=[np.asarray(pages[i].asarray()) for i in indices]
            page_map={str(i):int(i) for i in indices}
            channel=None
        assert n>=focus+6
        shape=arrays[0].shape
        assert all(a.shape==shape and a.ndim==2 for a in arrays),(gl,role,shape)
        assert shape==(300,106),(gl,role,shape)
        data=np.stack(arrays)
        npz=output/(gl+"_"+role+"_selected_source_frames.npz")
        np.savez_compressed(npz,original_frame_numbers=np.asarray(indices,dtype=np.int32),images=data)
        for k in [focus-1,focus,focus+1]:
            if k in indices:
                normalize_preview(data[indices.index(k)]).save(output/(gl+"_"+role+"_frame"+str(k)+".png"))
        return {"role":role,"frame_count":n,"source_page_count":len(pages),
                "source_page_shape":list(shape),"source_series_shape":list(series.shape),
                "source_series_axes":str(series.axes),"source_channel_zero_based":channel,
                "original_frame_to_tiff_page":page_map,
                "source_dtype":str(data.dtype),"selected_frame_indices":indices,
                "source_frames_npz_sha256":sha256(npz),"npz_bytes":npz.stat().st_size}

def main(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    source_bytes=E545_SOURCE.read_bytes()
    git_sha=hashlib.sha1(b"blob "+str(len(source_bytes)).encode()+b"\x00"+source_bytes).hexdigest()
    assert git_sha==E545_GIT_BLOB_SHA,("E545 source hash mismatch",git_sha)
    spec=importlib.util.spec_from_file_location("p901_e545",str(E545_SOURCE))
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    record={"study":"901","stage":"E554","status":"SOURCE_INTAKE_RUNNING",
            "source_extractor_git_blob":git_sha,"library":"S-BIAD2506",
            "biological_division_truth":"UNRESOLVED","operator_causality":"UNRESOLVED",
            "groups":[]}
    with tempfile.TemporaryDirectory(prefix="901_e554_") as temp:
        temp=Path(temp)
        for gl in ["Pos5_GL35","Pos12_GL22"]:
            cfg=mod.ARCHIVES[gl]
            reader=mod.SafeHTTPRange(mod.URL_BASE+cfg["name"],cfg["length"])
            result={"GL":gl,"original_archive":cfg["name"],"members":[],"image_samples":[]}
            with zipfile.ZipFile(reader,"r",allowZip64=True) as z:
                assert len(z.infolist())==cfg["entries"]
                reader.seek(z.start_dir)
                cdb=reader.read(cfg["cd_len"])
                assert hashlib.sha256(cdb).hexdigest()==cfg["cd_sha"]
                result["zip_central_directory_sha256"]=hashlib.sha256(cdb).hexdigest()
                root=cfg["root"]
                prefix_main=root+"20221220_analysis/track_data__20221220_analysis/"
                prefix_cost=root+"cost_function_analysis_1/track_data__cost_function_analysis_1/"
                selected=[]
                for info in z.infolist():
                    path=info.filename
                    if not path.startswith(root) or info.is_dir():continue
                    leaf=path.rsplit("/",1)[-1]
                    path_in_model=path.startswith(prefix_main) or path.startswith(prefix_cost) or (
                        "/track_data__" in path and "/cost_function_analysis_1/" in path)
                    if path_in_model and (leaf in {"curation.moma","moma.log","mm.properties"} or leaf.startswith("curation.moma__BKP")):
                        if info.file_size<=650000:selected.append(info)
                    if path==root+"20221220_glu_spcm_1_MMStack_"+gl+".csv":
                        selected.append(info)
                assert any(x.filename==prefix_main+"curation.moma" for x in selected)
                assert any(x.filename==prefix_cost+"curation.moma" for x in selected)
                for info in selected:
                    relative=info.filename[len(root):].replace("/","__")
                    smallpath=output/(gl+"__"+relative)
                    rec=extract_file(z,info,smallpath,cfg["name"],small_limit=650000)
                    if info.filename in (prefix_main+"curation.moma",prefix_cost+"curation.moma"):
                        tag="main" if info.filename.startswith(prefix_main) else "cost"
                        assert rec["extracted_sha256"]==EXPECTED_CURATIONS[(gl,tag)]
                        lines=smallpath.read_text(encoding="utf-8").splitlines()
                        rec["ASC_count"]=sum(line.strip().startswith("ASC,") for line in lines)
                        rec["PR_count"]=sum(line.strip().startswith("PR,") for line in lines)
                        rec["source_version_header"]=next(line for line in lines if "moma_version:" in line)
                    result["members"].append(rec)
                paths=[("main_probability",prefix_main+"probability_maps.tif"),
                       ("cost_probability",prefix_cost+"probability_maps.tif"),
                       ("phase_contrast",root+"20221220_glu_spcm_1_MMStack_"+gl+".tif")]
                for role,p in paths:
                    info=z.getinfo(p)
                    if role=="phase_contrast":assert info.compress_size==ROIS[gl]["phase_compressed"]
                    else:assert info.compress_size==ROIS[gl]["prob_compressed"]
                    tmp=temp/(gl+"__"+role+".tif")
                    rec=extract_file(z,info,tmp,cfg["name"])
                    obs=frames_tif(tmp,gl,role,output)
                    rec.update(obs);result["image_samples"].append(rec)
                    tmp.unlink()
                assert reader.total_transferred<1_500_000_000,"unbounded source transfer cap"
                result["source_HTTP206_requests"]=reader.nreq
                result["source_HTTP206_bytes"]=reader.total_transferred
            record["groups"].append(result)
            save_json(output/"E554_progress_receipt.json",record)
    record["status"]="OFFICIAL_ZIP64_CENTRAL_SHA_CRC_CURATIONS_AND_SOURCE_IMAGE_FRAMES_PASS"
    save_json(output/"E554_full_source_manifest.json",record)
    print("E554_OFFICIAL_SOURCE_PASS",json.dumps({"groups":[{"GL":x["GL"],"sidecars":len(x["members"]),"images":len(x["image_samples"]),"download_bytes":x["source_HTTP206_bytes"]} for x in record["groups"]]},ensure_ascii=False),flush=True)
if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--output",required=True);args=ap.parse_args()
    try:main(Path(args.output))
    except Exception as exc:
        out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
        save_json(out/"E554_negative_failure_receipt.json",{"status":"FAIL_CLOSED_NOT_SOURCE_VALIDATED","error_type":type(exc).__name__,"error":str(exc),"traceback":traceback.format_exc()})
        raise
