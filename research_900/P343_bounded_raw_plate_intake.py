#!/usr/bin/env python3
"""P343 bounded read-only acquisition; no author script execution, no source overwrites."""
from pathlib import Path
from urllib import request, parse
from zipfile import ZipFile
import hashlib, json, io, sys, time
ROOT=Path("P343_AUTHOR_RAW");ROOT.mkdir(exist_ok=True)
COMMIT="ee3a2079a4a7e47c090a7ad9945b053c72723490"
REPO="amichaibk/nutrient-dependent-auxotrophy_and_CF"
FILES=[
 ("Data/Crossfeeding Experiments/Bs_data_Ex1.xlsx",60054,"6dba3a358cbe449105daad027b20084b3273c81a"),
 ("Data/Crossfeeding Experiments/Bs_data_Ex2.xlsx",67370,"b9fcd786d6ec94801a356eba07a328d1c75ff964"),
 ("Data/Crossfeeding Experiments/Bs_data_Ex3.xlsx",68005,"b1078c055ca19323997c06b4608a8e52ddf4cd64"),
 ("Data/Crossfeeding Experiments/Ec_data_Ex1.xlsx",67457,"5ceb226dd5d9838b669cc9697ef2568c053a45c6"),
 ("Data/Crossfeeding Experiments/Ec_data_Ex2.xlsx",65764,"a2366f713165539048b908bfe63156dd12e45fad"),
 ("Data/Crossfeeding Experiments/Ec_data_Ex3.xlsx",67856,"212dfd2a15a98ad618bfd904e0773ff4d895d478"),
 ("Data/Sample Names.xlsx",11012,"cf6122517da99c5f03e9e9d91512ab8b52b02b4e"),
 ("Generated Tables/big_norm_df.csv",940382,"3621051892373a51ff7cddc1d4ab86fa5fb21edc"),
 ("Generated Tables/filtered_df.csv",740494,"9fd952175ebed9c4aa3bdb0fe600a8ce0a5e3e7f")
]
class ExactHostRedirect(request.HTTPRedirectHandler):
 def redirect_request(self,req,fp,code,msg,headers,url):
  p=parse.urlsplit(url)
  if p.scheme!="https" or p.hostname!="raw.githubusercontent.com":raise ValueError("REDIRECT_HOST_GATE")
  return super().redirect_request(req,fp,code,msg,headers,url)
opener=request.build_opener(ExactHostRedirect())
def sha1_git(blob):
 return hashlib.sha1(b"blob "+str(len(blob)).encode("ascii")+bytes([0])+blob).hexdigest()
receipt={"stage":"900_P343C","repository":REPO,"commit":COMMIT,"bytes_cap":3_500_000,
 "source_status":"NOT_CHECKED","file_count_expected":len(FILES),"results":[],"only_exact_public_files":True}
used=0
for rel,size,expected in FILES:
 record={"relative_path":rel,"expected_bytes":size,"expected_git_blob_sha1":expected,"status":"UNAVAILABLE"}
 try:
  if not size or size>1_200_000 or used+size>receipt["bytes_cap"]:raise ValueError("BYTE_CAP")
  safe_path="/".join(parse.quote(v,safe="") for v in rel.split("/"))
  url=f"https://raw.githubusercontent.com/{REPO}/{COMMIT}/{safe_path}"
  last_error=None
  for attempt in range(2):
   try:
    req=request.Request(url,headers={"User-Agent":"900-P343-source-qualification/1.0"})
    with opener.open(req,timeout=28) as resp:
     if resp.status!=200 or parse.urlsplit(resp.geturl()).hostname!="raw.githubusercontent.com":raise ValueError("HTTP_OR_HOST_GATE")
     cl=resp.headers.get("Content-Length")
     if cl and int(cl)>size:raise ValueError("RESPONSE_TOO_LARGE")
     blob=resp.read(size+1)
    if len(blob)!=size:raise ValueError(f"BYTE_COUNT_GATE:{len(blob)}")
    if sha1_git(blob)!=expected:raise ValueError("PINNED_GIT_BLOB_MISMATCH")
    break
   except Exception as e:
    last_error=e
    if attempt==0:time.sleep(1)
  else:raise last_error
  if rel.endswith(".xlsx"):
   if not blob.startswith(b"PK"):raise ValueError("XLSX_MAGIC_INVALID")
   with ZipFile(io.BytesIO(blob)) as z:
    if z.testzip() is not None:raise ValueError("XLSX_CRC_ERROR")
    if not {"[Content_Types].xml","xl/workbook.xml"}.issubset(set(z.namelist())):raise ValueError("XLSX_MISSING_WORKBOOK")
    record["ooxml_members"]=len(z.namelist())
  target=ROOT/rel
  target.parent.mkdir(exist_ok=True,parents=True)
  if target.exists():raise ValueError("REFUSE_REPLACE_EXISTING")
  target.write_bytes(blob)
  used+=size
  record.update(status="SHA1_GIT_BLOB_AND_SIZE_ACCEPTED",sha256=hashlib.sha256(blob).hexdigest(),actual_bytes=size)
 except Exception as e:
  record.update(status="REJECTED_OR_UNAVAILABLE",error_type=type(e).__name__,error=str(e)[:180])
 receipt["results"].append(record)
 print(json.dumps(record,ensure_ascii=False,sort_keys=True),flush=True)
receipt["accepted"]=sum(x["status"]=="SHA1_GIT_BLOB_AND_SIZE_ACCEPTED" for x in receipt["results"])
receipt["bytes_received_verified"]=used
receipt["source_status"]="ALL_9_PINNED" if receipt["accepted"]==len(FILES) else "PARTIAL_OR_NO_SOURCE"
(ROOT/"P343_SOURCE_RECEIPT.json").write_text(json.dumps(receipt,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf8")
print(json.dumps({"source_status":receipt["source_status"],"accepted":receipt["accepted"],"bytes":used},sort_keys=True),flush=True)
