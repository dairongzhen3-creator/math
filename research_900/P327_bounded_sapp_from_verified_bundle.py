#!/usr/bin/env python3
"""Read the public PNAS supplementary bundle; accept only the sha-locked copy, extract SI Appendix PDF."""
import hashlib,io,json,urllib.request,urllib.parse,zipfile
from pathlib import Path
OUT=Path("P327_OFFICIAL_SI_SOURCE"); OUT.mkdir(exist_ok=True)
URL="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13143037/supplementaryFiles"
EXPECTED_S3_SHA="ee8badd7d32ca8da477ef9140e9319d229a01b1e6c7d1d0af2f347b261f29e30"
MAX=10_000_000
report={"source":URL,"status":"SOURCE_NOT_RECEIVED","target":"pnas.2532423123.sapp.pdf","expected_signed_S3_sha256":EXPECTED_S3_SHA}
try:
 req=urllib.request.Request(URL,headers={"User-Agent":"900-P327-read-only-source-validation/1.0"})
 with urllib.request.urlopen(req,timeout=35) as response:
  if response.status!=200:raise ValueError("BAD_STATUS")
  host=urllib.parse.urlsplit(response.geturl()).hostname
  if host not in ("www.ebi.ac.uk","ebi.ac.uk"):raise ValueError("BAD_FINAL_HOST")
  if response.headers.get("Content-Length") and int(response.headers["Content-Length"])>MAX:raise ValueError("TOO_LARGE_HEADER")
  payload=response.read(MAX+1)
 if len(payload)>MAX:raise ValueError("TOO_LARGE")
 actual=hashlib.sha256(payload).hexdigest()
 with zipfile.ZipFile(io.BytesIO(payload)) as z:
  s3=[n for n in z.namelist() if n.lower().endswith("pnas.2532423123.sd03.xlsx")]
  if len(s3)!=1 or hashlib.sha256(z.read(s3[0])).hexdigest()!=EXPECTED_S3_SHA: raise ValueError("SIGNED_S3_MEMBER_SHA256_MISMATCH")
  if z.testzip():raise ValueError("ZIP_BAD_CRC")
  names=[n for n in z.namelist() if n.lower().endswith("pnas.2532423123.sapp.pdf")]
  if len(names)!=1:raise ValueError("PDF_AMBIGUOUS")
  item=z.getinfo(names[0])
  if item.file_size>3000000:raise ValueError("PDF_LIMIT")
  pdf=z.read(names[0])
  if not pdf.startswith(b"%PDF-"):raise ValueError("PDF_NOT_PDF")
  (OUT/"pnas.2532423123.sapp.pdf").write_bytes(pdf)
  report.update(status="SI_APPENDIX_VERIFIED",bundle_bytes=len(payload),bundle_sha256=actual,zip_members=len(z.namelist()),pdf_bytes=len(pdf),pdf_sha256=hashlib.sha256(pdf).hexdigest(),archive_member_name=names[0])
except Exception as exc:
 report.update(status="INTAKE_STOP",error_type=type(exc).__name__,error=str(exc)[:200])
(OUT/"P327_SI_OFFICIAL_RECEIPT.json").write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf8")
print(json.dumps(report,ensure_ascii=False,sort_keys=True),flush=True)
