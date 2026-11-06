#!/usr/bin/env python3
"""P330: bounded official PNAS S2 retrieval using S3 member SHA continuity."""
from pathlib import Path
from urllib import request, parse
from zipfile import ZipFile
import io, json, hashlib

OUT=Path("P330_PNAS_S2_OFFICIAL"); OUT.mkdir(exist_ok=True)
URL="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13143037/supplementaryFiles"
S3_SHA="ee8badd7d32ca8da477ef9140e9319d229a01b1e6c7d1d0af2f347b261f29e30"
S3="pnas.2532423123.sd03.xlsx"
S2="pnas.2532423123.sd02.xlsx"
CAP=10_000_000
def sha(b):return hashlib.sha256(b).hexdigest()
class SafeRedirect(request.HTTPRedirectHandler):
 def redirect_request(self,req,fp,code,msg,headers,url):
  p=parse.urlsplit(url)
  if p.scheme!="https" or p.hostname not in ("www.ebi.ac.uk","ebi.ac.uk"):raise ValueError("UNAPPROVED_REDIRECT")
  return super().redirect_request(req,fp,code,msg,headers,url)
report={"stage":"P330C","source":URL,"target":S2,"signed_S3_expected_sha256":S3_SHA,"status":"SOURCE_NOT_RECEIVED","max_bytes":CAP}
try:
 req=request.Request(URL,headers={"User-Agent":"900-P330-official-S2/1.0","Accept":"application/zip,application/octet-stream"})
 with request.build_opener(SafeRedirect()).open(req,timeout=35) as r:
  if r.status!=200:raise ValueError("HTTP_NOT_200")
  p=parse.urlsplit(r.geturl())
  if p.scheme!="https" or p.hostname not in ("www.ebi.ac.uk","ebi.ac.uk"):raise ValueError("BAD_FINAL_HOST")
  if r.headers.get("Content-Length") and int(r.headers["Content-Length"])>CAP:raise ValueError("HEADER_TOO_LARGE")
  bundle=r.read(CAP+1)
 if len(bundle)>CAP or not bundle.startswith(b"PK"):raise ValueError("BUNDLE_BAD_OR_OVERSIZED")
 with ZipFile(io.BytesIO(bundle)) as z:
  if z.testzip() is not None:raise ValueError("BUNDLE_CRC_FAIL")
  ss=[n for n in z.namelist() if n.lower().endswith(S3)]
  ts=[n for n in z.namelist() if n.lower().endswith(S2)]
  if len(ss)!=1 or len(ts)!=1:raise ValueError("MISSING_OR_DUPLICATE_MEMBER")
  if z.getinfo(ss[0]).file_size>500000 or sha(z.read(ss[0]))!=S3_SHA:raise ValueError("S3_SIGNATURE_FAIL")
  if z.getinfo(ts[0]).file_size>100000:raise ValueError("S2_TOO_LARGE")
  s2=z.read(ts[0])
  with ZipFile(io.BytesIO(s2)) as x:
   if x.testzip() is not None:raise ValueError("S2_CRC_FAIL")
   names=set(x.namelist())
   if not {"[Content_Types].xml","xl/workbook.xml"}<=names:raise ValueError("S2_NOT_XLSX")
   n=len(names)
  (OUT/S2).write_bytes(s2)
  report.update(status="S2_VERIFIED_S3_CONTINUITY",bundle_bytes=len(bundle),bundle_sha256=sha(bundle),s2_bytes=len(s2),s2_sha256=sha(s2),s2_member=ts[0],s3_member=ss[0],s2_ooxml_members=n,zip_members=len(z.namelist()))
except Exception as e:
 report.update(status="SOURCE_REJECTED_OR_UNAVAILABLE",error_type=type(e).__name__,error=str(e)[:200])
(OUT/"P330_S2_OFFICIAL_SOURCE_RECEIPT.json").write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf8")
print(json.dumps(report,ensure_ascii=False,sort_keys=True),flush=True)
