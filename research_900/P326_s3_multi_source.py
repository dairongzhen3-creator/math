#!/usr/bin/env python3
"""Bounded, read-only public source validation for PNAS S3; no executing downloaded content."""
import io, json, hashlib, urllib.request, urllib.parse, zipfile
from pathlib import Path
OUT=Path('P326_PNAS_S3_SOURCE'); OUT.mkdir(exist_ok=True)
URLS=[
 'https://www.pnas.org/doi/suppl/10.1073/pnas.2532423123/suppl_file/pnas.2532423123.sd03.xlsx',
 'https://pmc.ncbi.nlm.nih.gov/articles/PMC13143037/bin/pnas.2532423123.sd03.xlsx',
 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC13143037/bin/pnas.2532423123.sd03.xlsx',
 'https://europepmc.org/articles/PMC13143037/bin/pnas.2532423123.sd03.xlsx',
 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13143037/supplementaryFiles',
]
MAX=10_000_000
ALLOW=('pnas.org','ncbi.nlm.nih.gov','europepmc.org','ebi.ac.uk')
def host_ok(h): return h and any(h==s or h.endswith('.'+s) for s in ALLOW)
class SafeRedirect(urllib.request.HTTPRedirectHandler):
 def redirect_request(self, req, fp, code, msg, headers, url):
  p=urllib.parse.urlparse(url)
  if p.scheme!='https' or not host_ok(p.hostname):raise RuntimeError('UNSAFE_REDIRECT')
  return super().redirect_request(req,fp,code,msg,headers,url)

def grab(url):
 req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 (public-academic-source-verification/1.0)','Accept':'application/octet-stream,*/*'})
 with urllib.request.build_opener(SafeRedirect()).open(req,timeout=28) as r:
  u=urllib.parse.urlparse(r.geturl());assert u.scheme=='https' and host_ok(u.hostname)
  if r.status!=200:raise RuntimeError('HTTP_NOT_200')
  if r.headers.get('Content-Length') and int(r.headers['Content-Length'])>MAX:raise RuntimeError('TOO_LARGE_BY_HEADER')
  b=r.read(MAX+1)
  if len(b)>MAX: raise RuntimeError('TOO_LARGE')
  return b,{'response_url':r.geturl(),'mime':r.headers.get('Content-Type'), 'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}

def readbook(b):
 if not b.startswith(b'PK'):return False,0
 try:
  with zipfile.ZipFile(io.BytesIO(b)) as z:
   if z.testzip() is not None:return False,0
   names=z.namelist()
   return ('[Content_Types].xml' in names and 'xl/workbook.xml' in names),len(names)
 except (zipfile.BadZipFile,ValueError):return False,0

result={'project':'900','stage':'P326C','max_per_request':MAX,'records':[],'found_valid_S3':False}
for i,url in enumerate(URLS):
 row={'candidate':i,'url':url,'status':'UNAVAILABLE'}
 try:
  b,metadata=grab(url);row.update(metadata)
  ok,n=readbook(b)
  if ok:
   row.update(status='S3_XLSX_VERIFIED',members=n)
   (OUT/'pnas.2532423123.sd03.xlsx').write_bytes(b)
   result['found_valid_S3']=True
  elif b.startswith(b'PK'):
   with zipfile.ZipFile(io.BytesIO(b)) as z:
    names=[name for name in z.namelist() if name.lower().endswith('pnas.2532423123.sd03.xlsx')]
    row.update(status='ZIP_WITHOUT_TARGET',zip_members=len(z.namelist()))
    if names:
     sub=z.read(names[0]);ok,n=readbook(sub)
     if ok:
      row.update(status='S3_XLSX_VERIFIED_FROM_ZIP',file_inside=names[0],internal_sha256=hashlib.sha256(sub).hexdigest(),internal_bytes=len(sub),members=n)
      (OUT/'pnas.2532423123.sd03.xlsx').write_bytes(sub);result['found_valid_S3']=True
  else:
   row.update(status='NOT_XLSX_OR_ZIP',first_32_bytes_hex=b[:32].hex())
 except Exception as e:
  row.update(status='SOURCE_ERROR',exception=type(e).__name__,error=str(e)[:140])
 result['records'].append(row)
 print(json.dumps(row,sort_keys=True),flush=True)
 if result['found_valid_S3']:break
(OUT/'P326_S3_RETRIEVAL_RECEIPT.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print('SOURCE_VERIFIED',result['found_valid_S3'],flush=True)
