#!/usr/bin/env python3
"""901 E117: exact source-level audit of independent publisher-author raw multielectrode files.
No third-party MATLAB or source code is executed. Does not establish plant/experiment independence.
"""
import csv, hashlib, io, json, math, re, urllib.request, zipfile, xml.etree.ElementTree as ET
from decimal import Decimal, getcontext
from pathlib import Path
getcontext().prec=40
REPO="Wenlong0-0/NE2026_plant_mapping"
DIR="Fig4 Mimosa touch bottom signal speed and amp calculation"
API="https://api.github.com/repos/"+REPO+"/contents/"+urllib.parse.quote(DIR)
USER_AGENT="901-source-integrity-AUDIT/1.0"
def read(url,maxb=800000):
    req=urllib.request.Request(url,headers={"User-Agent":USER_AGENT,"Accept":"application/vnd.github+json" if "/api.github.com/" in url else "application/octet-stream"})
    with urllib.request.urlopen(req,timeout=40) as r:
        s=r.read(maxb+1); status=r.status;final=r.geturl()
    if len(s)>maxb:raise ValueError("FILE_EXCEEDS_BOUND "+str(len(s)))
    return s,status,final
source,_,_=read(API,100000)
meta=json.loads(source)
files=sorted([m for m in meta if m.get("name","").endswith(".txt") and "replicate" in m.get("name","").lower()],key=lambda a:int(re.search(r'replicate\s*(\d+)',a["name"]).group(1)))
assert len(files)==8,{"actual":len(files),"names":[f["name"] for f in files]}
result={"project":"901","stage":"E117","paper_doi":"10.1038/s41928-026-01642-z","source_repository":REPO,
        "source_commit_ref":"main (content SHA per original file recorded)","author_matlab_execution":False,
        "audit_method":"read published numeric raw text; follow published MATLAB global-max peak rule; independent Decimal index QC",
        "channels_per_record_expected":11,"fs_Hz_from_author_script":100,"channel_distance_mm_from_author_script":2,
        "files":[],"author_workbook":{}}
allrows=[]
for m in files:
    content,code,final=read(m["download_url"],600000)
    assert code==200
    sha_blob=hashlib.sha1(b"blob "+str(len(content)).encode()+b"\0"+content).hexdigest()
    assert sha_blob==m["sha"],("GIT_BLOB_MISMATCH",m["name"])
    assert len(content)==m["size"]
    # header absent, each row has electrical potentials of 11 simultaneous channels in author's order
    text=content.decode("utf-8-sig")
    vals=list(csv.reader(io.StringIO(text)))
    assert vals and all(len(r)==11 for r in vals)
    numbers=[[float(t) for t in r] for r in vals]
    assert all(math.isfinite(x) for row in numbers for x in row)
    n=len(numbers)
    vmax=[max(row[j] for row in numbers) for j in range(11)]
    idx=[next(i for i,row in enumerate(numbers) if row[j]==vmax[j]) for j in range(11)]
    idx_dec=[]
    for j in range(11):
        col=[Decimal(row[j]) for row in vals]
        best=max(col);idx_dec.append(col.index(best))
    assert idx==idx_dec
    deltas=[idx[j+1]-idx[j] for j in range(10)]
    speed=[200/d if d else None for d in deltas]  # 2mm * 100Hz / index lag
    negative=sum(x is not None and x<0 for x in speed)
    zero=sum(x is None for x in speed)
    record={"name":m["name"],"github_blob_sha1":sha_blob,"raw_sha256":hashlib.sha256(content).hexdigest(),
            "bytes":len(content),"rows":n,"columns":11,"recording_seconds":n/100,
            "first_line_all_zero":all(float(s)==0 for s in vals[0]),
            "max_time_samples":idx,"author_rule_speed_mm_per_s":speed,"max_amplitude_mV_assuming_original_author_plot_units":vmax,
            "negative_velocity_edges_count":negative,"zero_lag_edges_count":zero}
    result["files"].append(record)
    for j,(d,v) in enumerate(zip(deltas,speed),1):
        allrows.append({"replicate":m["name"],"edge":j,"max_index_source":idx[j-1],"max_index_dest":idx[j],
                        "delta_index":d,"velocity_mm_s_author_rule":v})
# Fetch published author's compiled XLSX workbook; map each original row/column to its display cell
m=[x for x in meta if x.get("name")=="WaveSpeed_Amplitude_Matrices.xlsx"]
assert len(m)==1
file=m[0]; raw,code,_=read(file["download_url"],200000)
assert hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()==file["sha"]
result["author_workbook"]={"sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw),"expected_sheets":["Speed_mm_per_s","Amplitude_mV"]}
N="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
R="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
with zipfile.ZipFile(io.BytesIO(raw)) as z:
    assert z.testzip() is None
    shared=[]
    if "xl/sharedStrings.xml" in z.namelist():
        sr=ET.fromstring(z.read("xl/sharedStrings.xml"))
        shared=["".join(t.text or "" for t in s.iter() if t.tag==N+"t") for s in sr.findall(N+"si")]
    wb=ET.fromstring(z.read("xl/workbook.xml"))
    rels=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    lookup={x.attrib.get("Id"):x.attrib.get("Target") for x in rels}
    sheets={}
    for elem in wb.iter(N+"sheet"):
        name=elem.attrib.get("name");rid=elem.attrib.get(R+"id")
        target=lookup.get(rid,"").lstrip("/")
        path=target if target.startswith("xl/") else "xl/"+target
        if path not in z.namelist(): continue
        root=ET.fromstring(z.read(path))
        cval={}
        for c in root.iter(N+"c"):
            rawv=c.find(N+"v")
            if rawv is not None and rawv.text is not None:
                v=rawv.text
                if c.attrib.get("t")=="s":v=shared[int(v)]
            else:
                v="".join(t.text or "" for t in c.iter() if t.tag==N+"t")
            cval[c.attrib["r"]]=v
        sheets[name]=cval
    result["author_workbook"]["sheets"]=list(sheets)
    assert all(x in sheets for x in result["author_workbook"]["expected_sheets"])
    alphabet="ABCDEFGH"
    mismatches=[]
    successes=0; missing=0;author_non_numeric=0
    for h,src in enumerate(result["files"]):
        col=alphabet[h]
        name=sheets["Speed_mm_per_s"].get(col+"1","")
        if ("replicate "+str(h+1)) not in name.lower():
            mismatches.append({"header_column":col,"header":name,"expected":src["name"]})
        for j,v in enumerate(src["author_rule_speed_mm_per_s"],2):
            got=sheets["Speed_mm_per_s"].get(col+str(j),"")
            if v is None:
                if got=="" or got.lower() in ("nan","#num!","na"):successes+=1
                else:mismatches.append({"cell":col+str(j),"orig":"undefined","workbook":got})
            else:
                try:
                    if math.isclose(float(got),v,rel_tol=1e-8,abs_tol=1e-8):successes+=1
                    else:mismatches.append({"cell":col+str(j),"computed":v,"workbook":got})
                except ValueError:
                    mismatches.append({"cell":col+str(j),"computed":v,"workbook":got})
    result["author_workbook"]["compared_speed_cells"]=80
    result["author_workbook"]["match_count"]=successes
    result["author_workbook"]["mismatches"]=mismatches[:35]
result["source_scope"]="8 author-labelled recordings, NOT independently certified 8 biological plants; 100Hz, 11 voltage channels; no independent petiole angle/movement Y in this directory"
result["aggregate"]={
    "recordings":len(result["files"]),
    "raw_rows":sum(x["rows"] for x in result["files"]),
    "channels":11,
    "wave_edges":len(allrows),
    "defined_wave_edges":sum(x["velocity_mm_s_author_rule"] is not None for x in allrows),
    "negative_signed_edges":sum(x["velocity_mm_s_author_rule"] is not None and x["velocity_mm_s_author_rule"]<0 for x in allrows),
    "zero_delta_edges":sum(x["velocity_mm_s_author_rule"] is None for x in allrows),
    "independent_parent_n":"not documented in checked source"
}
out=Path("901_E117_fig4_source_audit.json")
out.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
with open("901_E117_fig4_speed_by_edge.csv","w",newline="",encoding="utf-8") as f:
    wr=csv.DictWriter(f,fieldnames=allrows[0]);wr.writeheader();wr.writerows(allrows)
print("E117_FIG4_SUMMARY",json.dumps({"aggregate":result["aggregate"],"author_workbook":result["author_workbook"],"by_file":[{k:r[k] for k in ("name","rows","recording_seconds","negative_velocity_edges_count","zero_lag_edges_count")} for r in result["files"]]},ensure_ascii=False),flush=True)
print("E117_FIG4_RESULTS_SHA",hashlib.sha256(out.read_bytes()).hexdigest(),flush=True)
