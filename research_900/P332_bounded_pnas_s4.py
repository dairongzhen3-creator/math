#!/usr/bin/env python3
"""900 P332C: one bounded read-only official PNAS Dataset S4 receipt."""
from pathlib import Path
from urllib import request,parse
from zipfile import ZipFile
import io,json,hashlib
OUT=Path("P332_OFFICIAL_S4");OUT.mkdir(exist_ok=True)
URL="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13143037/supplementaryFiles"
S3_SHA="ee8badd7d32ca8da477ef9140e9319d229a01b1e6c7d1d0af2f347b261f29e30"
TAR="pnas.2532423123.sd04.xlsx"; S3="pnas.2532423123.sd03.xlsx";CAP=10_000_000
def sha(b):return hashlib.sha256(b).hexdigest()
class RedirectGuard(request.HTTPRedirectHandler):
 def redirect_request(self,req,fp,code,msg,headers,url):
  p=parse.urlsplit(url)
  if p.scheme!="https" or p.hostname not in ("www.ebi.ac.uk","ebi.ac.uk"):raise ValueError("REDIRECT_NOT_APPROVED")
  return super().redirect_request(req,fp,code,msg,headers,url)
receipt={"stage":"P332C","source":URL,"target":TAR,"expected_signed_s3_sha256":S3_SHA,"status":"NO_SOURCE","total_byte_cap":CAP}
try:
 req=request.Request(URL,headers={"User-Agent":"900-P332-bounded-scientific-dataset-verification/1.0"})
 with request.build_opener(RedirectGuard()).open(req,timeout=35) as r:
  if r.status!=200:raise ValueError("HTTP_NOT_200")
  p=parse.urlsplit(r.geturl())
  if p.scheme!="https" or p.hostname not in ("www.ebi.ac.uk","ebi.ac.uk"):raise ValueError("BAD_FINAL_HOST")
  if r.headers.get("Content-Length") and int(r.headers.get("Content-Length"))>CAP:raise ValueError("HEADER_TOO_LARGE")
  body=r.read(CAP+1)
 if len(body)>CAP:raise ValueError("BUNDLE_TOO_LARGE")
 with ZipFile(io.BytesIO(body)) as z:
  if z.testzip() is not None:raise ValueError("BUNDLE_CRC_ERROR")
  s3=[n for n in z.namelist() if n.lower().endswith(S3)]
  s4=[n for n in z.namelist() if n.lower().endswith(TAR)]
  if len(s3)!=1 or len(s4)!=1:raise ValueError("MISSING_OR_DUPLICATE_S3_S4")
  if z.getinfo(s3[0]).file_size>500_000:raise ValueError("S3_TOO_LARGE")
  if sha(z.read(s3[0]))!=S3_SHA:raise ValueError("SIGNED_SOURCE_CONTINUITY_FAILED")
  if z.getinfo(s4[0]).file_size>2_000_000:raise ValueError("S4_TOO_LARGE")
  bb=z.read(s4[0])
  if not bb.startswith(b"PK"):raise ValueError("S4_NOT_XLSX")
  with ZipFile(io.BytesIO(bb)) as f:
   if f.testzip() is not None:raise ValueError("S4_CRC_ERROR")
   if not {"[Content_Types].xml","xl/workbook.xml"}.issubset(f.namelist()):raise ValueError("S4_OOXML_MISSING")
   inner_members=len(f.namelist())
  (OUT/TAR).write_bytes(bb)
  receipt.update(status="S4_ACCEPTED_WITH_S3_SHA_CONTINUITY",bundle_bytes=len(body),bundle_sha256=sha(body),zip_members=len(z.namelist()),s4_bytes=len(bb),s4_sha256=sha(bb),s4_ooxml_members=inner_members,s4_member=s4[0])
except Exception as e:
 receipt.update(status="S4_REJECTED_OR_UNAVAILABLE",error_type=type(e).__name__,error=str(e)[:200])
(OUT/"P332_OFFICIAL_S4_RECEIPT.json").write_text(json.dumps(receipt,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf8")
print(json.dumps(receipt,ensure_ascii=False,sort_keys=True),flush=True)
