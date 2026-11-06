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
 text=decoded.decode('utf-8-sig')
 rows=csv.reader(io.StringIO(text,newline=''),delimiter='\t')
 header=next(rows)
 if len(header)!=len(set(header)) or len(header)<2: raise ValueError('BAD_HEADER')
 names=EXPECTED[series]
 # Only exact published sample labels or exact GSM accessions are unambiguous.
 gm={'GSE208376':[f'GSM{i}' for i in range(6342032,6342041)],'GSE270548':[f'GSM{i}' for i in range(8346505,8346511)]}[series]
 found=names if set(names).issubset(header) else gm if set(gm).issubset(header) else []
 sample_index=[header.index(x) for x in found]
 gene_count=0
 valid_int=True; all_nonnegative=True
 totals=[0]*len(found); zeros=[0]*len(found)
 preview=[]
 for row in rows:
  if len(row)!=len(header): raise ValueError('RAGGED_TABLE')
  if not any(v.strip() for v in row): continue
  gene_count+=1
  if gene_count<=3: preview.append(row[0][:60])
  if found:
   for j,col in enumerate(sample_index):
    val=row[col].strip()
    if not val.isdecimal(): valid_int=False; continue
    iv=int(val); totals[j]+=iv; zeros[j]+=(iv==0)
 summary={'header_preview':header[:24], 'n_columns':len(header),'gene_or_table_rows':gene_count,'id_preview':preview,
          'published_sample_labels_exact':bool(found),'matched_column_names':found, 'integer_and_nonnegative':valid_int and all_nonnegative,
          'totals_if_sample_counts':dict(zip(found,totals)) if found and valid_int else {},
          'zeros_if_sample_counts':dict(zip(found,zeros)) if found and valid_int else {},
          'decoded_bytes':len(decoded)}
 if found and valid_int and gene_count>=100 and len(found)==len(names):
  summary['data_class']='REPORTED_SAMPLE_COUNT_MATRIX_PROVISIONAL';summary['status']='SCHEMA_PASS_CULTURE_AND_BATCH_PENDING'
 else:
  summary['data_class']='UNVERIFIED_DERIVED_OR_INCOMPATIBLE';summary['status']='NO_RAW_COUNT_CERTIFICATE'
 return summary

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
