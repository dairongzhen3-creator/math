#!/usr/bin/env python3
"""P333: tightly bounded read-only collection of pinned public Joller 2026 source bytes."""
import hashlib,json,urllib.request,urllib.parse
from pathlib import Path
OWNER="PMI-Basel"; REPO="Joller_et_al_Glass_Clay"
COMMIT="b9844d740f33f8420a78d68b7b08272ee4b0bba5"
BASE="B_statistical_analysis/03_Timeseries_Glass"
FILES={
"Input/Re_PQN_normalized.csv":("2331c8ce2201c01c222e7b54c6bec2efdce80197",166051),
"Input/Re_normalized.csv":("04b6c2ed5a509684b0e61ac962c6abe651bb4f21",183628),
"Input/DESIGN_T18.xlsx":("67b1bc1c7310403ac164fbf48cb9b6630146ce37",14581),
"Input/DESIGN_T18.RDS":("16002a09a3e1087e9de4e9a69967244199467759",1163),
"Input/syn_melt.RDS":("63fd82722d7d129c3f02f5cec61ede1419808aab",28704),
"Input/fia_ions.RDS":("b3b3a7857c14b2d0398fea3807c32e117b88a6d4",65468),
"fia_targetCompounds.RDS":("01afe6ed8b5fa9c8b8d4c923dcafb825fd4d3f8a",14949),
"Input/CFURoot.RDS":("e5dabf19ac358802646bdef0ad97c8135fde8e21",3468),
"Input/PlateLayouts.xlsx":("ca126fc6098ef769cfce5d9e71716b49c523f739",16363),
"Input/Mer_PQN_normalized.csv":("b3e36046b317db59908421f37341fd409d591f47",1326191)
}
ROOT=Path("P333_JOLLER_PINNED_SOURCES")
ROOT.mkdir(exist_ok=True)
report={"stage":"P333C","upstream":"https://github.com/"+OWNER+"/"+REPO,"commit":COMMIT,"status":"SOURCE_NOT_RECEIVED","files":[],"aggregate_byte_cap":3_000_000}
class Guard(urllib.request.HTTPRedirectHandler):
 def redirect_request(self,req,fp,code,msg,headers,url):
  p=urllib.parse.urlsplit(url)
  if p.scheme!="https" or p.hostname!="raw.githubusercontent.com": raise ValueError("DISALLOWED_REDIRECT")
  return super().redirect_request(req,fp,code,msg,headers,url)
opener=urllib.request.build_opener(Guard())
aggregate=0
for rel,(expected,size) in FILES.items():
 path=BASE+"/"+rel
 url=f"https://raw.githubusercontent.com/{OWNER}/{REPO}/{COMMIT}/{path}"
 record={"relative_path":path,"expected_git_blob":expected,"expected_bytes":size,"status":"NOT_RECEIVED"}
 try:
  if size>2_000_000 or aggregate+size>3_000_000:raise ValueError("BOUNDS_EXCEEDED")
  req=urllib.request.Request(url,headers={"User-Agent":"900-P333-bounded-source/1.0"})
  with opener.open(req,timeout=25) as h:
   if h.status!=200 or urllib.parse.urlsplit(h.geturl()).hostname!="raw.githubusercontent.com":raise ValueError("HOST_STATUS")
   if h.headers.get("Content-Length") and int(h.headers.get("Content-Length"))>size:raise ValueError("HEADER_SIZE_EXCEEDED")
   raw=h.read(size+1)
  if len(raw)!=size:raise ValueError(f"BYTE_COUNT_MISMATCH:{len(raw)}")
  gitsha=hashlib.sha1(b"blob "+str(len(raw)).encode("ascii")+bytes([0])+raw).hexdigest()
  if gitsha!=expected:raise ValueError("PINNED_GIT_BLOB_SHA_MISMATCH")
  dest=ROOT/rel;dest.parent.mkdir(exist_ok=True,parents=True);dest.write_bytes(raw)
  aggregate+=len(raw)
  record.update(status="ACCEPTED_GIT_BLOB_PINNED",bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),git_blob_sha1=gitsha)
 except Exception as e:
  record.update(status="REJECTED",error_type=type(e).__name__,error=str(e)[:180])
 report["files"].append(record)
 print(json.dumps(record,sort_keys=True),flush=True)
report["total_bytes"]=aggregate
report["accepted"]=sum(x["status"]=="ACCEPTED_GIT_BLOB_PINNED" for x in report["files"])
report["status"]="ALL_10_PINNED_SOURCES_RECEIVED" if report["accepted"]==len(FILES) else "PARTIAL_OR_REJECTED"
(ROOT/"P333_SOURCE_RECEIPT.json").write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf8")
print(json.dumps({"status":report["status"],"accepted":report["accepted"],"total_bytes":aggregate},sort_keys=True),flush=True)
