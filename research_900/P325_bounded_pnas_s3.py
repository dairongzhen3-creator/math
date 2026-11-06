#!/usr/bin/env python3
"""900 P325C single bounded publisher-side evidence fetch. No third-party code runs."""
import hashlib, json, urllib.request, urllib.parse, zipfile
from pathlib import Path
OUT=Path("P325_PUBLISHER_S3_RECEIPT"); OUT.mkdir(exist_ok=True)
URL="https://pmc.ncbi.nlm.nih.gov/articles/instance/13143037/bin/pnas.2532423123.sd03.xlsx"
MAX=2500000
class StrictRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        z=urllib.parse.urlparse(newurl)
        if z.scheme!="https" or not z.hostname.endswith(".ncbi.nlm.nih.gov"):
            raise ValueError("redirect rejected")
        return super().redirect_request(req, fp, code, msg, headers, newurl)
record={"project":"900","phase":"P325C","source":URL,"status":"SOURCE_UNAVAILABLE","max_bytes":MAX}
try:
    opener=urllib.request.build_opener(StrictRedirect())
    req=urllib.request.Request(URL,headers={"User-Agent":"900-reproducibility-public-access/1.0","Accept":"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"})
    with opener.open(req,timeout=30) as response:
        final=urllib.parse.urlparse(response.geturl())
        if final.scheme!="https" or not final.hostname.endswith(".ncbi.nlm.nih.gov"): raise ValueError("final host")
        if response.status!=200: raise ValueError("bad HTTP status")
        length=response.headers.get("Content-Length")
        if length and int(length)>MAX: raise ValueError("content length exceeds bound")
        data=response.read(MAX+1)
    if len(data)>MAX: raise ValueError("download exceeds bound")
    if not data.startswith(b"PK"):raise ValueError("not OOXML zip")
    with zipfile.ZipFile(__import__("io").BytesIO(data)) as z:
        bad=z.testzip()
        if bad: raise ValueError("zip CRC failure")
        names=z.namelist()
        if "[Content_Types].xml" not in names or not any(x.startswith("xl/workbook") for x in names): raise ValueError("not Excel")
    (OUT/"pnas.2532423123.sd03.xlsx").write_bytes(data)
    record.update(status="SOURCE_BYTES_VERIFIED",bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),members=len(names),host=final.hostname)
except Exception as exc:
    record.update(error_class=type(exc).__name__,error=str(exc)[:300])
(OUT/"P325_S3_SOURCE_RECEIPT.json").write_text(json.dumps(record,sort_keys=True,indent=2)+"\n")
print(json.dumps(record,sort_keys=True),flush=True)
