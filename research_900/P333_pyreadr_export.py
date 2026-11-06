#!/usr/bin/env python3
"""Convert accepted Joller public RDS data to inspectable CSV, no remote code execution."""
import csv,hashlib,json,platform
from pathlib import Path
import pyreadr
SRC=Path("P333_JOLLER_PINNED_SOURCES")
receipt=json.loads((SRC/"P333_SOURCE_RECEIPT.json").read_text())
if receipt.get("status")!="ALL_10_PINNED_SOURCES_RECEIVED" or receipt.get("accepted")!=10:raise RuntimeError("UPSTREAM_VERIFICATION_NOT_ACCEPTED")
OUT=SRC/"derived";OUT.mkdir(exist_ok=True)
records=[]
for rel in ("Input/syn_melt.RDS","Input/fia_ions.RDS","fia_targetCompounds.RDS","Input/CFURoot.RDS","Input/DESIGN_T18.RDS"):
 res={"file":rel,"status":"UNAVAILABLE"}
 try:
  row=next(x for x in receipt["files"] if x["relative_path"].endswith(rel))
  raw=(SRC/rel).read_bytes()
  sha=hashlib.sha256(raw).hexdigest()
  if row["status"]!="ACCEPTED_GIT_BLOB_PINNED" or row["sha256"]!=sha:raise RuntimeError("HASH_MISMATCH")
  data=pyreadr.read_r(str(SRC/rel))
  res["objects"]=[str(x) for x in data.keys()]
  dfs=[x for x in data.values() if hasattr(x,"to_csv") and hasattr(x,"columns")]
  if len(dfs)!=1:raise RuntimeError("NOT_ONE_TABULAR_OBJECT")
  df=dfs[0]
  if len(df)>200000 or len(df.columns)>500:raise RuntimeError("CONVERSION_TOO_LARGE")
  fname=Path(rel).stem+"_SOURCE_EXPORT.csv"
  out=OUT/fname
  df.to_csv(out,index=False,na_rep="NA",float_format="%.17g")
  res.update(status="TABLE_EXPORTED",rows=len(df),cols=len(df.columns),columns=list(map(str,df.columns)),output=fname,output_bytes=out.stat().st_size,output_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),source_sha256=sha)
 except Exception as e:
  res.update(status="CONVERSION_FAILED_OR_NON_TABULAR",error_type=type(e).__name__,error=str(e)[:220])
 records.append(res)
 print(json.dumps(res,sort_keys=True,default=str),flush=True)
(OUT/"P333_RDS_PYREADR_EXPORT_RECEIPT.json").write_text(json.dumps({"tool":"pyreadr","version":getattr(pyreadr,"__version__","unavailable"),"python":platform.python_version(),"files":records},ensure_ascii=False,sort_keys=True,indent=2)+"\n")
