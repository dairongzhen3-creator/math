#!/usr/bin/env python3
"""901 E117 Publisher official source data: direct HTTPS bounded read, read-only tabulation, no fake samples."""
import urllib.request,zipfile,io,hashlib,xml.etree.ElementTree as ET,json,collections,re
from pathlib import Path
PREFIX="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41928-026-01642-z/MediaObjects/41928_2026_1642_MOESM"
targets=[(9,"Published Fig.4 Mimosa electrical response/source data"),
         (16,"Extended Data Fig.8 synchronized sensing, signalling and movement response"),
         (18,"Extended Data Fig.10 choy sum Modes I/II")]
ns="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
rid_ns="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
rows=[]
for num,role in targets:
    url=PREFIX+str(num)+"_ESM.xlsx"
    r={"publication":"10.1038/s41928-026-01642-z","file_number":num,"source_role":role,"source_link":url}
    try:
        req=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0 (901 scientific source audit)","Accept":"application/octet-stream"})
        with urllib.request.urlopen(req,timeout=45) as res:
            r["status"]=res.status
            r["content_type"]=res.headers.get("Content-Type")
            b=res.read(4000001)
        if len(b)>4000000:raise ValueError("BOUND_EXCEEDED_4MB")
        if b[:2]!=b"PK": raise ValueError("NOT_XLSX_ZIP")
        r["sha256"]=hashlib.sha256(b).hexdigest()
        r["raw_bytes"]=len(b)
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            bad=z.testzip()
            if bad:raise ValueError("CRC_MISMATCH "+str(bad))
            ss=[]
            if "xl/sharedStrings.xml" in z.namelist():
                root=ET.fromstring(z.read("xl/sharedStrings.xml"))
                ss=["".join(s.text or "" for s in item.iter() if s.tag==ns+"t") for item in root.findall(ns+"si")]
            wb=ET.fromstring(z.read("xl/workbook.xml"))
            rel=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
            ids={m.get("Id"):m.get("Target") for m in rel}
            pages=[]
            for elem in wb.iter(ns+"sheet"):
                name=elem.attrib.get("name")
                target=ids.get(elem.get(rid_ns+"id"),"")
                path=target.lstrip("/") if target.startswith("/") else "xl/"+target
                if path not in z.namelist(): pages.append({"name":name,"unresolved_target":path});continue
                s=ET.fromstring(z.read(path))
                cells={}
                formulas=0;nonblank=0
                for c in s.iter(ns+"c"):
                    v=c.find(ns+"v"); inline=c.find(ns+"is");f=c.find(ns+"f")
                    if f is not None:formulas+=1
                    if v is not None and v.text is not None:
                        x=v.text
                        if c.get("t")=="s":x=ss[int(x)] if int(x)<len(ss) else "<missing shared>"
                    elif inline is not None:x="".join(t.text or "" for t in inline.iter() if t.tag==ns+"t")
                    else:x=""
                    if x!="":cells[c.get("r","")]=x;nonblank+=1
                # visible head and complete compact rows if publication source is short
                rkeys=sorted({int(re.findall(r"\d+",k)[0]) for k in cells},key=int)
                head=[]
                for i in rkeys[:32]:
                    matches=sorted((k,v) for k,v in cells.items() if re.fullmatch(r"[A-Z]+"+str(i),k))
                    head.append({"excel_row":i,"cells":matches[:40]})
                keys=list(cells)
                labels=[{"cell":k,"value":v[:200]} for k,v in cells.items() if re.search("time|sensing|signal|response|touch|hot|cold|movement|petiole|stimulus|plant|average|mean|SD|angle|second|speed|amplitude|sample",str(v),re.I) and len(v)<220][:65]
                pages.append({"sheet":name,"cells":len(cells),"formulas":formulas,"row_count_with_values":len(rkeys),"max_row":max(rkeys) if rkeys else 0,
                              "rows_first_32":head,"matched_labels":labels})
            r["pages"]=pages
            r["xlsx_CRC_all_pass"]=True
    except Exception as e:
        r["error"]=type(e).__name__
        r["reason"]=str(e)[:300]
    rows.append(r)
    print("E117_PUBLISHER_SOURCE",json.dumps({k:v for k,v in r.items() if k!="pages"},ensure_ascii=False),flush=True)
    for p in r.get("pages",[]):
        print("E117_TAB",json.dumps({"source_number":num,"sheet":p.get("sheet"),"cell_count":p.get("cells"),"row_count":p.get("row_count_with_values"),"preview":p.get("rows_first_32",[])[:12]},ensure_ascii=False)[:15000],flush=True)
p=Path("901_E117_publisher_source_xlsx_gate.json")
p.write_text(json.dumps({"stage":"E117","publication":"Li et al Nature Electronics 2026","sources":rows,"no_author_matlab_executed":True,"new_biological_experiments":0},ensure_ascii=False,indent=2),encoding="utf-8")
print("E117_PUBLISHER_JSON_SHA",hashlib.sha256(p.read_bytes()).hexdigest(),flush=True)
