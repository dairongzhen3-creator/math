#!/usr/bin/env python3
"""Recombine source artifacts from Google Drive E498 multipart originals.
Preserves source bytes; never trusts file names without matching SHA256.
Usage: python reassemble_verify_e498_source_artifacts.py --dir /path/to/source_parts --out /path/to/recovered
"""
from pathlib import Path
import zipfile,hashlib,json,argparse
p=argparse.ArgumentParser();p.add_argument('--dir',default='source_parts');p.add_argument('--out',default='recovered_originals');ns=p.parse_args()
source=Path(ns.dir);out=Path(ns.out);out.mkdir(parents=True,exist_ok=True)
manifest=json.loads((source/'E498_source_artifacts_split_manifest.json').read_text())
sha=lambda b:hashlib.sha256(b).hexdigest()
for record in manifest['artifacts']:
 chunks=[]
 for x in record['parts']:
  b=(source/x['filename']).read_bytes()
  assert len(b)==x['bytes'] and sha(b)==x['sha256'],x['filename']
  chunks.append(b)
 all_=b''.join(chunks)
 assert len(all_)==record['full_bytes'] and sha(all_)==record['full_sha256']
 p=out/record['name'];p.write_bytes(all_)
 with zipfile.ZipFile(p) as z: assert z.testzip() is None
 print('PASS',p,len(all_),sha(all_))
