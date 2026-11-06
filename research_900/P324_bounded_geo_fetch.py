#!/usr/bin/env python3
"""P324C one-shot public NCBI GEO intake, bounded; no DE and no secret tokens."""
import csv, gzip, hashlib, io, json, sys, time, urllib.request, urllib.parse
from pathlib import Path

SOURCES = [
  ('GSE208376', 'GSE208376_all_sample_plant_cnt.tsv.gz', 'GSE208nnn', 9),
  ('GSE270548', 'GSE270548_all_compare.txt.gz', 'GSE270nnn', 6),
]
MAX_TRANSFER = 4_194_304
MAX_DECODED = 64_000_000
ACCEPTED_HOSTS = {'ftp.ncbi.nlm.nih.gov', 'www.ncbi.nlm.nih.gov'}
EXPECTED = {
 'GSE208376': ['mock_1','mock_2','mock_3','flg22_1','flg22_2','flg22_3','2_AA_1','2_AA_2','2_AA_3'],
 'GSE270548': ['Mock1','Mock2','Mock3','AA1','AA2','AA3'],
}

class SafeRedirect(urllib.request.HTTPRedirectHandler):
 def redirect_request(self, request, fp, code, msg, headers, newurl):
  target=urllib.parse.urlparse(newurl)
  if target.scheme!='https' or target.hostname not in ACCEPTED_HOSTS:
   raise ValueError('UNEXPECTED_REDIRECT')
  return super().redirect_request(request, fp, code, msg, headers, newurl)

def fetch_bounded(url):
 opener=urllib.request.build_opener(SafeRedirect())
 request=urllib.request.Request(url, headers={'User-Agent':'900-P324C-reproducibility-pilot/1.0', 'Accept':'application/gzip'})
 with opener.open(request,timeout=45) as response:
  final=urllib.parse.urlparse(response.geturl())
  if final.scheme!='https' or final.hostname not in ACCEPTED_HOSTS:
   raise ValueError('UNEXPECTED_HOST')
  if response.status!=200: raise ValueError('HTTP_NOT_200')
  content_length=response.headers.get('Content-Length')
  if content_length and int(content_length)>MAX_TRANSFER: raise ValueError('CONTENT_LENGTH_TOO_LARGE')
  data=bytearray()
  while True:
   chunk=response.read(min(65536, MAX_TRANSFER+1-len(data)))
   if not chunk: break
   data+=chunk
   if len(data)>MAX_TRANSFER: raise ValueError('TRANSFER_LIMIT')
 if len(data)<2 or data[:2]!=b'\x1f\x8b': raise ValueError('NOT_GZIP')
 return bytes(data)

def inspect(data, series):
 with gzip.GzipFile(fileobj=io.BytesIO(data)) as gz:
  decoded=gz.read(MAX_DECODED+1)
  if len(decoded)>MAX_DECODED: raise ValueError('DECODED_LIMIT')
 charset='utf-16' if decoded.startswith((b'\xff\xfe',b'\xfe\xff')) else 'utf-8-sig'
 text=decoded.decode(charset)
 rows=csv.reader(io.StringIO(text,newline=''),delimiter='\t')
 header=next(rows)
 if len(header)!=len(set(x for x in header if x)) + sum(not x for x in header):
  raise ValueError('DUPLICATE_NONBLANK_HEADER')
 if series=='GSE208376':
  required=['mock_1','mock_2','mock_3','flg22_1','flg22_2','flg22_3',"2'AA_1","2'AA_2","2'AA_3"]
  if [x for x in header if x]!=required or len(header)!=15: raise ValueError('PLANT_HEADER_MISMATCH')
  cols=required;use_recovered_header=True;expected_width=10
 else:
  required=['Mock1_count','Mock2_count','Mock3_count','AA1_count','AA2_count','AA3_count']
  if not all(x in header for x in required): raise ValueError('BACT_COUNT_COLS_MISSING')
  cols=required;use_recovered_header=False;expected_width=len(header)
 indexes=[(i+1) if use_recovered_header else header.index(x) for i,x in enumerate(cols)]
 totals=[0]*len(cols);zeros=[0]*len(cols);ids=set();n=0
 for row in rows:
  if len(row)!=expected_width: raise ValueError('RAGGED_ROWS')
  if not row[0] or row[0] in ids: raise ValueError('MISSING_OR_DUPLICATED_GENE_ID')
  ids.add(row[0]);n+=1
  for i,c in enumerate(indexes):
   val=row[c].strip()
   if not val.isdecimal(): raise ValueError('NON_INTEGER_OR_MISSING_COUNTS')
   k=int(val);totals[i]+=k;zeros[i]+=(k==0)
 if n<100: raise ValueError('TOO_FEW_GENES')
 return {'status':'SOURCE_SCHEMA_PASS','data_class':'AUTHOR_PROVIDED_PER_SAMPLE_INTEGER_GENE_COUNTS',
         'source_char_encoding':charset,'header_width_original':len(header),
         'source_header_has_extra_blank_slots':use_recovered_header,
         'observed_raw_count_columns':cols,'gene_rows':n,
         'sample_library_totals':dict(zip(cols,totals)),
         'sample_zero_gene_counts':dict(zip(cols,zeros)),
         'decoded_bytes':len(decoded),'no_biological_de_claim':True}

def run(out):
 out=Path(out);out.mkdir(parents=True,exist_ok=True)
 report={'project':'900','phase':'P324C','public_sources_only':True,'utc_run':'source-runner','no_rnaseq_de':True,
         'source_count':len(SOURCES),'results':[]}
 for acc,basename,group,n in SOURCES:
  url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{group}/{acc}/suppl/{basename}'
  row={'accession':acc,'filename':basename,'official_url':url,'expected_samples':n,'status':'SOURCE_NOT_RECEIVED'}
  try:
   data=fetch_bounded(url)
   row.update(sha256=hashlib.sha256(data).hexdigest(),bytes=len(data),source_received=True)
   (out/basename).write_bytes(data)
   row.update(inspect(data,acc))
  except Exception as e:
   row.update(status='INTAKE_STOP',failure_class=type(e).__name__,failure_message=str(e)[:180],source_received=bool(row.get('source_received',False)))
  report['results'].append(row)
  print(json.dumps(row,ensure_ascii=False),flush=True)
 (out/'P324_GEO_RECEIPT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf8')
 return report

if __name__=='__main__':
 run(sys.argv[1] if len(sys.argv)>1 else 'P324_GEO_PUBLIC_INTAKE')
