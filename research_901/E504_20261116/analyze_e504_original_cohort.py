#!/usr/bin/env python3
"""901 E504: genuine pre-starvation risk-set and cross-phase identity audit.
Original Gervais et al. 2026; no causality, dose, or physical-molecule claims.
Only exact source IDs; never infer genealogy from dotted-string prefixes.
"""
import os,csv,json,hashlib,collections,pathlib,math,time
import numpy as np
import pandas as pd
import pyarrow.feather as feather

OUT=pathlib.Path("e504_audit");OUT.mkdir(exist_ok=True)
F=pathlib.Path("original/myframes_control.feather")
HASH="5cb1850cb127a64bf5e3dffa4d4e13e7cc648f73af176f6d03cf4d0e1fa9bb7c"
def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(2**20),b""):h.update(b)
    return h.hexdigest()
assert F.stat().st_size==302295912 and sha256(F)==HASH
A=pathlib.Path("original/E500_original_illumination_mother_fate_candidates.csv")
assert A.exists()
cols=["date","condition","strain","gl_id","cell","c_step","b_rank","end_type","time_start","time_end","m_start","m_end","time_sec"]
tab=feather.read_table(F,columns=cols)
assert tab.num_rows==1011012
df=tab.to_pandas()
for c in ["date","condition","strain","gl_id","cell","c_step"]:
    df[c]=df[c].astype(str)
K=["date","condition","gl_id","cell"]
L=["date","condition","gl_id"]
CUT=43200.0
assert set(df.c_step)=={"GM.1","SPM.1","GM.2"}
assert set(df.date)=={"20221220"}
assert set(df.condition)=={"M9_M9zero_60h","M9_M9zero_60h_NO_LIGHT"}
def keyset(q):
    return set(q[K].itertuples(index=False,name=None))
def group_count(q):
    if not len(q):return {}
    return {str(k):int(v) for k,v in q.groupby("condition",dropna=False).size().items()}
def status(q, who):
    z=q.assign(label=who).copy();return z
phase_members={p:keyset(df.loc[df.c_step==p,K].drop_duplicates()) for p in ("GM.1","SPM.1","GM.2")}
cross={p:len(phase_members[p]) for p in phase_members}
cross["GM1_SPM1_exact"]=len(phase_members["GM.1"] & phase_members["SPM.1"])
cross["GM1_GM2_exact"]=len(phase_members["GM.1"] & phase_members["GM.2"])
cross["SPM1_GM2_exact"]=len(phase_members["SPM.1"] & phase_members["GM.2"])
cross["all_three_exact"]=len(set.intersection(*phase_members.values()))
cross["union_all_exact"]=len(set.union(*phase_members.values()))
# only contemporaneous observed rows define eligibility: no future time_end selections
pre=df[(df.c_step=="GM.1")&(df.time_sec<CUT)].copy()
assert len(pre)>0
last_time=pre.groupby(L,sort=False,dropna=False).time_sec.transform("max")
pre_last=pre.loc[pre.time_sec.eq(last_time)].copy()
t0_at=pre_last[L+["time_sec"]].drop_duplicates()
assert t0_at.groupby(L).size().max()==1
t0_delta=CUT-t0_at.time_sec
r0_all=pre_last[K].drop_duplicates()
r0_mother=pre_last.loc[pre_last.b_rank.eq(0),K].drop_duplicates()
pre_any_mother=pre.loc[pre.b_rank.eq(0),K].drop_duplicates()
pre_any_all=pre[K].drop_duplicates()
sp=df[df.c_step=="SPM.1"]
sp_early=sp.loc[sp.time_sec.eq(43200),K].drop_duplicates()
gm2=df[df.c_step=="GM.2"]
gm2_before=gm2.loc[(gm2.time_start<gm2.m_start)&(gm2.b_rank==0)].copy()
gm2_before["fate"]=np.select(
  [gm2_before.end_type.eq("div"),gm2_before.end_type.isin(["eod","lysis"])],
  ["Survive","No Growth"],default="Unknown")
gm2_before.loc[(gm2_before.fate=="Unknown")&(gm2_before.time_end>gm2_before.m_start+72000),"fate"]="No Growth"
sel=gm2_before[K+["strain","end_type","fate","time_start","time_end","m_start"]].drop_duplicates()
assert len(sel)==238 and len(keyset(sel))==238
old=pd.read_csv(A,dtype={"date":str,"cell":str,"gl_id":str})
old["date"]=old.date.astype(str)
assert len(old)==238
selected_cols=["date","condition","gl_id","cell","strain","end_type","fate","time_start","time_end","m_start"]
m=sel.merge(old[selected_cols],on=K,suffixes=("_new","_old"),validate="one_to_one",indicator=True)
assert len(m)==238 and m._merge.eq("both").all()
for c in ["strain","end_type","fate","time_start","time_end","m_start"]:
    new=m[c+"_new"];old_col=m[c+"_old"]
    assert (np.isclose(new,old_col,atol=1e-9).all() if c in ("time_start","time_end","m_start") else new.eq(old_col).all()),c
R0=keyset(r0_mother); R_all=keyset(r0_all); P_pre=keyset(pre_any_mother)
gm2_any=phase_members["GM.2"]; SP=phase_members["SPM.1"]; GM1=phase_members["GM.1"]
selected=keyset(sel)
assert selected.issubset(gm2_any)
members={"GM1_any":GM1,"SPM1_any":SP,"GM2_any":gm2_any,"GM1_pre_allrank":keyset(pre_any_all),
         "GM1_pre_mother_any":P_pre,"GM1_last_frame_anyrank":R_all,
         "GM1_last_frame_mother_rank0":R0,"SPM1_first_frame_43200":keyset(sp_early),
         "GM2_author_selected":selected}
overlap={}
for nm,ss in members.items():
    overlap[nm]={"unique_exact_cell_keys":len(ss),
        "intersect_GM2_author_238":len(ss & selected),
        "intersect_GM2_any":len(ss & gm2_any),
        "conditions":dict(collections.Counter(k[1] for k in ss))}
selected_out=sel.copy()
for nm in ("GM1_pre_allrank","GM1_pre_mother_any","GM1_last_frame_anyrank","GM1_last_frame_mother_rank0","SPM1_any","SPM1_first_frame_43200"):
    ss=members[nm];selected_out[nm]=[tuple(x) in ss for x in selected_out[K].itertuples(index=False,name=None)]
selected_out["tracked_start_before_exposure"]=selected_out.time_start<CUT
selected_out["tracked_start_after_or_at_exposure"]=selected_out.time_start>=CUT
selected_out.sort_values(K).to_csv(OUT/"E504_238_selected_exact_phase_identity_bridge.csv",index=False)
# Exposure-start cohort is "last observed pre-exposure frame" per gl_id, not biological exposure-randomized.
r0=r0_mother.copy()
r0["GM2_observed_same_exact_cell"]=[tuple(x) in gm2_any for x in r0[K].itertuples(index=False,name=None)]
r0["SPM1_observed_same_exact_cell"]=[tuple(x) in SP for x in r0[K].itertuples(index=False,name=None)]
r0["GM2_author_selected_same_exact_cell"]=[tuple(x) in selected for x in r0[K].itertuples(index=False,name=None)]
r0["last_pre_exposure_time_sec"]=r0.set_index(L).index.map(t0_at.set_index(L).time_sec.to_dict())
r0.sort_values(K).to_csv(OUT/"E504_pre_exposure_last_frame_rank0_exact_IDs.csv",index=False)
summ=[]
for (date,cond,gl),v in pre_last.groupby(L,sort=True):
    this=(date,cond,gl)
    r0q=v.loc[v.b_rank.eq(0),K].drop_duplicates()
    g0key=keyset(r0q)
    rows0=v[K].drop_duplicates()
    ss=sel[(sel.date==date)&(sel.condition==cond)&(sel.gl_id==gl)]
    summ.append({"date":date,"condition":cond,"gl_id":gl,
        "last_pre_exposure_sec":float(v.time_sec.max()),"seconds_to_starvation":CUT-float(v.time_sec.max()),
        "last_frame_allrank_unique_cells":len(rows0),
        "last_frame_rank0_unique_cells":len(g0key),
        "last_frame_rank0_exact_ID_present_in_SPM1":len(g0key & SP),
        "last_frame_rank0_exact_ID_present_in_GM2":len(g0key & gm2_any),
        "last_frame_rank0_exact_ID_author_selected":len(g0key & selected),
        "GM2_author_selected_cells_at_position":len(ss),
        "GM2_author_selected_pretracked":int((ss.time_start<CUT).sum()),
        "GM2_author_selected_posttracked":int((ss.time_start>=CUT).sum())})
by_gl=pd.DataFrame(summ).sort_values(L)
by_gl.to_csv(OUT/"E504_pre_exposure_riskset_by_GL.csv",index=False)
assert by_gl.seconds_to_starvation.min()>=0
# Timing-aware breakdown: descriptions, not an unconditional causal denominator.
all_res={}
for p in ["M9_M9zero_60h","M9_M9zero_60h_NO_LIGHT"]:
    x=selected_out[selected_out.condition==p]
    prex=x[x.tracked_start_before_exposure]
    late=x[~x.tracked_start_before_exposure]
    all_res[p]={
      "author_selected":len(x),"selected_with_pretracked_time":len(prex),"selected_after_or_at_t0":len(late),
      "selected_pretrained_exact_GM1_key":int(x.GM1_pre_allrank.sum()),
      "selected_exact_GM1_lastframe_rank0":int(x.GM1_last_frame_mother_rank0.sum()),
      "selected_exact_SPM1_key":int(x.SPM1_any.sum()),
      "selected_fate":x.fate.value_counts().to_dict(),
      "selected_pretracked_fate":prex.fate.value_counts().to_dict(),
      "selected_latertracked_fate":late.fate.value_counts().to_dict()}
# independent native Arrow batch traversal (not pandas groupby or merges): IDs at last frames and selection
last_by_gl={(r.date,r.condition,r.gl_id):float(r.last_pre_exposure_sec) for r in by_gl.itertuples()}
alt_phases=collections.defaultdict(set)
alt_pre_mother=set();alt_last_mother=set();alt_selected={}
rowcount=0
for batch in tab.to_batches(max_chunksize=16384):
    for rec in batch.to_pylist():
        rowcount+=1
        k=(str(rec["date"]),str(rec["condition"]),str(rec["gl_id"]),str(rec["cell"]))
        phase=rec["c_step"];alt_phases[phase].add(k)
        if phase=="GM.1" and rec["time_sec"]<CUT and rec["b_rank"]==0:
            alt_pre_mother.add(k)
            glk=k[:3]
            if float(rec["time_sec"])==last_by_gl.get(glk):alt_last_mother.add(k)
        if phase=="GM.2" and rec["time_start"]<rec["m_start"] and rec["b_rank"]==0:
            fate=("Survive" if rec["end_type"]=="div" else
                  "No Growth" if rec["end_type"] in ("eod","lysis") or (rec["end_type"] not in ("div","eod","lysis") and rec["time_end"]>rec["m_start"]+72000) else "Unknown")
            if k in alt_selected:assert alt_selected[k]==fate
            alt_selected[k]=fate
assert rowcount==1011012
assert all(alt_phases[p]==phase_members[p] for p in alt_phases)
assert alt_pre_mother==P_pre and alt_last_mother==R0
assert set(alt_selected)==selected
assert all(alt_selected[tuple(r[c] for c in K)]==r["fate"] for _,r in sel.iterrows())
rec={
 "study":"Gervais et al 2026; original Zenodo 18059846","stage":"901 E504",
 "source_Feather_sha256":HASH,"source_zip_archive_sha256":"43f88ea9986b25a6810712b732fb6ff5435ee9ced55d30d901aa02e1de3f464d",
 "rows_authenticated":rowcount,"date":"20221220","stages":["GM.1","SPM.1","GM.2"],
 "starvation_transition_sec":CUT,"counterfactual_or_exposure_randomization":"NOT_IDENTIFIED",
 "phase_id_overlap":cross,
 "cross_eligibility_overlap":overlap,
 "by_condition":all_res,
 "last_pre_exposure_frame_Gl_count":len(by_gl),
 "last_pre_exposure_frames_sec_to_transition_range":[float(by_gl.seconds_to_starvation.min()),float(by_gl.seconds_to_starvation.max())],
 "last_frame_rank0_count":len(R0),
 "last_frame_rank0_to_author_238_exact_match":len(R0 & selected),
 "independent_arrow_recordbatch_qa_pass":True,
 "unknown_later_outcome_if_absent_from_GM2":"NOT_DEAD",
 "source_semantics":"Cell identity exact composite date,condition,gl_id,cell only; b_rank=0 means author-defined mother rank; no inferred genealogy from names.",
 "science_limits":["pre-exposure last-frame cohort is conditional on source trackability","same location never contains both illumination states","any GM2 loss is tracking nonobservation not biological death","no independent absolute photon-dose and protein copy-number calibration","no bit/J or physical molecular information"],"QA_PASS":True,
 "timestamp_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime())}
(OUT/"E504_authenticated_pre_exposure_cohort_results.json").write_text(json.dumps(rec,indent=2,ensure_ascii=False,default=str))
(OUT/"README_E504_REPRODUCE.md").write_text("""Original source pinned to source SHA256 and historic E500 official archive. This audit reads and verifies all 1,011,012 authentic Feather records, then computes per exact cell ID and physical GL which cell is seen in the last growth-phase frame before t=43200s. Rank0 refers to AUTHOR mother rank and is not a random cohort. Independent Arrow RecordBatch iteration revalidates exact source IDs, all stages, pre-risk and author E500 238 IDs. No mother-daughter genealogy inferred from IDs; no causal/energy/information-rate claims. Compare JSON and per-GL/per-ID CSVs, retain previous negative results. No full raw image retrieval or biological replication performed.\n""")
print("E504_FINAL",json.dumps({"QA":rec["QA_PASS"],"phase_id_overlap":cross,"by_condition":all_res,"exact_pre_r0":len(R0),"exact_pre_to_238":len(R0 & selected),"near_t0_range_sec":rec["last_pre_exposure_frames_sec_to_transition_range"]},ensure_ascii=False),flush=True)
