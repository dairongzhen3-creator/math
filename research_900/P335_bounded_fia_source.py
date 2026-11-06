#!/usr/bin/env python3
"""900 P335: bounded read-only pinned public FIA raw RDS source and tabular schema."""
import hashlib, json, urllib.request, urllib.parse
from pathlib import Path
import pandas as pd
import pyreadr
ROOT = Path("P335_FIA_RAW_SOURCE")
ROOT.mkdir(exist_ok=True)
SOURCE = "https://raw.githubusercontent.com/PMI-Basel/Joller_et_al_Glass_Clay/b9844d740f33f8420a78d68b7b08272ee4b0bba5/B_statistical_analysis/03_Timeseries_Glass/Input/fia_peaks.RDS"
SIZE = 2825204
GIT_BLOB = "476e87993300c65d93a8c43afd0101610341ee0e"
RECEIPT = {"project":900,"stage":"P335C","source_url":SOURCE,"expected_bytes":SIZE,"expected_git_blob_sha1":GIT_BLOB,"status":"UNVERIFIED"}
def finish():
 (ROOT/"P335_FIA_RAW_SOURCE_RECEIPT.json").write_text(json.dumps(RECEIPT,ensure_ascii=False,sort_keys=True,indent=2,default=str)+"\n",encoding="utf8")
 print(json.dumps({"status":RECEIPT["status"],"rows":RECEIPT.get("rows"),"columns":RECEIPT.get("columns"),"error":RECEIPT.get("error")},sort_keys=True),flush=True)
class RedirectGuard(urllib.request.HTTPRedirectHandler):
 def redirect_request(self,req,fp,code,msg,headers,url):
  p=urllib.parse.urlsplit(url)
  if p.scheme!="https" or p.hostname!="raw.githubusercontent.com":
   raise RuntimeError("SOURCE_REDIRECT_NOT_ALLOWED")
  return super().redirect_request(req,fp,code,msg,headers,url)
try:
 req=urllib.request.Request(SOURCE,headers={"User-Agent":"900-P335-pinned-public-source-qa/1.0"})
 with urllib.request.build_opener(RedirectGuard()).open(req,timeout=35) as h:
  if h.status!=200 or urllib.parse.urlsplit(h.geturl()).hostname!="raw.githubusercontent.com":raise RuntimeError("BAD_HTTP_STATUS_OR_HOST")
  if h.headers.get("Content-Length") and int(h.headers["Content-Length"])>SIZE:raise RuntimeError("BAD_SOURCE_CONTENT_LENGTH")
  raw=h.read(SIZE+1)
 if len(raw)!=SIZE:raise RuntimeError("SOURCE_BYTE_COUNT_MISMATCH")
 blob=hashlib.sha1(b"blob "+str(len(raw)).encode("ascii")+bytes([0])+raw).hexdigest()
 if blob!=GIT_BLOB:raise RuntimeError("PINNED_GIT_BLOB_MISMATCH")
 target=ROOT/"fia_peaks.RDS"
 target.write_bytes(raw)
 RECEIPT.update(status="RAW_SHA_LOCKED",actual_bytes=len(raw),git_blob_sha1=blob,sha256=hashlib.sha256(raw).hexdigest())
 try:
  objects=pyreadr.read_r(str(target))
  RECEIPT["object_keys"]=[str(k) for k in objects]
  frames=[d for d in objects.values() if isinstance(d,pd.DataFrame)]
  if len(frames)!=1:raise RuntimeError("RDS_NOT_EXACTLY_ONE_DATAFRAME")
  df=frames[0]
  if len(df)>2000000 or len(df.columns)>300:raise RuntimeError("RDS_EXPORT_DIMENSION_LIMIT")
  RECEIPT["rows"]=len(df);RECEIPT["columns"]=list(map(str,df.columns))
  RECEIPT["column_dtypes"]={str(c):str(df[c].dtype) for c in df.columns}
  out=ROOT/"fia_peaks_source_export.csv.gz"
  df.to_csv(out,index=False,na_rep="",compression={"method":"gzip","compresslevel":6,"mtime":0})
  RECEIPT["csv_gz_bytes"]=out.stat().st_size
  RECEIPT["csv_gz_sha256"]=hashlib.sha256(out.read_bytes()).hexdigest()
  cat_cols=[c for c in df.columns if str(c).lower() in ("plate","sample","jslname","timepoint","solvent","type","ionidx","well","injection_order","injectionorder","runorder","run_order","sampleid")]
  RECEIPT["observed_matching_columns"]=list(map(str,cat_cols))
  for c in cat_cols:
   RECEIPT["distinct_"+str(c)]=int(df[c].nunique(dropna=True))
  if "Plate" in df.columns and "sample" in df.columns:
   pairs=df[["Plate","sample"]].drop_duplicates()
   RECEIPT["plate_sample_unique_pairs"]=len(pairs)
   RECEIPT["plate_levels"]=sorted(str(x) for x in df["Plate"].dropna().unique())[:30]
   tab=pairs.groupby("Plate",dropna=False).size().reset_index(name="unique_samples")
   tab.to_csv(ROOT/"P335_PLATE_SAMPLE_COUNTS.csv",index=False)
  RECEIPT["status"]="RAW_AND_TABULAR_EXPORT_VERIFIED"
 except Exception as e:
  RECEIPT["status"]="RAW_SHA_LOCKED_BUT_RDS_CONVERSION_FAILED"
  RECEIPT["conversion_error_type"]=type(e).__name__
  RECEIPT["error"]=str(e)[:300]
except Exception as e:
 RECEIPT["status"]="RAW_SOURCE_REJECTED"
 RECEIPT["error_type"]=type(e).__name__
 RECEIPT["error"]=str(e)[:300]
finish()
