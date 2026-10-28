#!/usr/bin/env python3
"""E117: independently check three-stage per-sample latencies from publisher source figure ED8c.
The observations are published author measurements, not new experiments or independent predictive validation.
"""
import csv,urllib.request,zipfile,hashlib,io,math,statistics,json,decimal,xml.etree.ElementTree as ET
from pathlib import Path
url="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41928-026-01642-z/MediaObjects/41928_2026_1642_MOESM16_ESM.xlsx"
sha_expected="8a62524f8137fe438c799929f6cf376a1d5139a190418f0769f4cc5b3d9725d4"
req=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0","Accept":"application/octet-stream"})
with urllib.request.urlopen(req,timeout=45) as res:b=res.read(200001)
assert len(b)<=200000 and b[:2]==b"PK"
assert hashlib.sha256(b).hexdigest()==sha_expected
N="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
with zipfile.ZipFile(io.BytesIO(b)) as z:
    assert z.testzip() is None
    path="xl/worksheets/sheet1.xml"
    root=ET.fromstring(z.read(path))
    vals={}
    for c in root.iter(N+"c"):
        v=c.find(N+"v")
        if v is not None and v.text is not None:vals[c.get("r")]=v.text
        else:
            ins=c.find(N+"is")
            if ins is not None:vals[c.get("r")]="".join(t.text or "" for t in ins.iter() if t.tag==N+"t")
groups=[("Touch","A","B","C",12),("Cold","E","F","G",11),("Hot","I","J","K",14)]
headers=[("Touch Sensing","Touch Signalling","Touch Response"),("Cold Sensing","Cold Signalling","Cold Response"),("Hot Sensing","Hot Signalling","Hot Response")]
obs=[]
validation=[]
for (group,a,b,c,expected),head in zip(groups,headers):
    assert tuple(vals.get(col+"1") for col in (a,b,c))==head
    data=[]
    for row in range(2,16):
        raw=[vals.get(col+str(row)) for col in (a,b,c)]
        if all(x is None for x in raw):continue
        assert all(x is not None and x!="" for x in raw),(group,row,raw)
        vec=[float(x) for x in raw]
        assert all(math.isfinite(x) and x>0 for x in vec)
        data.append((row,raw,vec))
    assert len(data)==expected,(group,len(data),expected)
    for n,(row,raw,vec) in enumerate(data,1):
        sens,signal,response=vec
        total=sum(vec)
        obs.append({"stimulus":group,"source_excel_row":row,"within_stimulus_observation_no":n,
            "sensing_s":sens,"signalling_s":signal,"response_s":response,
            "total_s":total,"fraction_signalling":signal/total,"signalling_minus_response_s":signal-response})
    # independent exact-decimal parsing and sums cross-check primary float:
    decimals=[[decimal.Decimal(s) for s in raw] for _,raw,_ in data]
    exact=[sum(v[i] for v in decimals)/len(decimals) for i in range(3)]
    main=[statistics.mean(x[i] for _,raw,x in data) for i in range(3)]
    assert max(abs(float(exact[i])-main[i]) for i in range(3))<1e-12
    responses=[v[2] for _,raw,v in data]
    totals=[sum(v) for _,raw,v in data]
    xs=[v for _,raw,v in data]
    validation.append({"group":group,"author_reported_n":expected,"observed_paired_triplets":len(data),"stage_names":head,
        "n_compared_cells":3*expected,"means_s":{"sensing":main[0],"signalling":main[1],"response":main[2]},
        "sample_sd_s":{"sensing":statistics.stdev(v[0] for v in xs),
                       "signalling":statistics.stdev(v[1] for v in xs),
                       "response":statistics.stdev(v[2] for v in xs)},
        "total_time_mean_s":statistics.mean(totals),
        "signalling_fraction_mean":statistics.mean(v[1]/sum(v) for v in xs),
        "response_at_or_under_0_5_s_count":sum(v<=.5 for v in responses),
        "response_over_0_5_s_count":sum(v>.5 for v in responses),
        "all_signal_stage_longer_than_response":all(v[1]>v[2] for v in xs),
        "within_row_data_sha256":hashlib.sha256(("\n".join(",".join(raw) for _,raw,_ in data)).encode()).hexdigest()})
assert len(obs)==37
r={"stage":"901_E117","origin":"publisher-original Nature Electronics Extended Data Fig.8 Source Data",
"publisher_url":url,"original_xlsx_sha256":sha_expected,"original_file_bytes":len(b),
"source_validation":"OOXML literal per-row raw cells + independent Decimal recomputation (111 values)",
"original_biological_n_certified":"author says independent samples 12+11+14, no individual plant_id/vessel/date keys in this source",
"new_wet_laboratory_observations":0,"new_independent_external_HA3_predictions":0,
"real_measured_physiological_outcome":"same-study petiole drooping/angle-based response onset, shown as per-trial response-stage delay",
"future_drought_24_48h":"NOT_TESTED","direct_energy_bit_per_J":"NOT_IDENTIFIED",
"conditions":validation,
"overall":{"measured_trial_groups":37,"stage_triplet_cells":111,
 "response_below_or_equal_0_5_s":sum(d["response_at_or_under_0_5_s_count"] for d in validation),
 "response_above_0_5_s":sum(d["response_over_0_5_s_count"] for d in validation),
 "all_response_positive":all(v["response_s"]>0 for v in obs),
 "all_signalling_longer_than_response":all(v["signalling_minus_response_s"]>0 for v in obs)},
"limitations":["source table is latency intervals not raw voltage or movie-linked plantID-by-sample time series",
               "published individual experimental trials are not authenticated as independent plants",
               "unrandomized experiments cannot isolate causal added predictive information in electricity relative to stimulus-only baseline",
               "no external test split for forecast loss, and no 24h/48h physiological target"]}
out=Path("901_E117_ED8c_latency_audit.json")
out.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding="utf8")
with open("901_E117_ED8c_37_observations.csv","w",newline="",encoding="utf8") as f:
    w=csv.DictWriter(f,fieldnames=list(obs[0]));w.writeheader();w.writerows(obs)
print("E117_ED8C_SUMMARY",json.dumps({"per_condition":validation,"overall":r["overall"],"source_sha":sha_expected},ensure_ascii=False),flush=True)
print("E117_ED8C_SHA",hashlib.sha256(out.read_bytes()).hexdigest(),flush=True)
