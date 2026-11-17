#!/usr/bin/env python3
"""901 E520 exact original-source audit: selected MoMA Div variables vs CellStats after PR.

Source: officially ZIP64/CRC/SHA validated E508 and E512 original source artifacts.
This identifies *software* objective/export differences, not actual cell-division truth.
"""
import argparse,csv,hashlib,io,json,zipfile,pathlib,collections,re
P=argparse.ArgumentParser();P.add_argument('--base',default='/mnt/data');P.add_argument('--output',required=True)
a=P.parse_args();B=pathlib.Path(a.base);OUT=pathlib.Path(a.output);OUT.mkdir(exist_ok=True,parents=True)
C={
 'Pos12_GL22':{'base':'20221220_E08_GFP_pl','sha508':'9207589c8be47c613415bb196b3be6251fbec59b6f2219ae5b65c196368f9832','sha512':'73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2'},
 'Pos5_GL35':{'base':'20221220_rplN_GFP_pl','sha508':'4818288be9f78536ac6ef71d718f3f529a784c70c0979161fd8b6ad64e2b9375','sha512':'77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6'},
}
def sha(b):return hashlib.sha256(b).hexdigest()
def open_checked(path,ref):
 buf=path.read_bytes();assert sha(buf)==ref,(path,sha(buf));z=zipfile.ZipFile(io.BytesIO(buf));assert z.testzip() is None
 return z

def find(z,end):
 all=[n for n in z.namelist() if n.endswith(end)]
 assert len(all)==1,(end,all)
 return all[0],z.read(all[0])
def sel_sol(text):
 m={}
 for s in text.decode().splitlines():
  if not s.strip() or s.startswith('#'):continue
  kv=s.split();assert len(kv)==2
  k,val=kv
  assert k not in m
  m[k]=float(val)
 return m

def active_div(s):return {k for k,v in s.items() if k.startswith('DivT') and v>.5}
def tvar(k):
 m=re.match(r'^(?:Map|Div|Lys|Exit)T(\d+)_',k)
 assert m,k
 return int(m[1])
def original(z,rec,v,name):
 n,data=find(z,v+'__'+name)
 matches=[r for r in rec if r['kind']==name and r['version']==v]
 assert len(matches)==1
 assert sha(data)==matches[0]['original_sha256'] and len(data)==matches[0]['uncompressed_bytes']
 return data

def cells(z,gl,v,cfg):
 filename=f'{gl}_{v}_CellStats.csv'
 n,b=find(z,filename)
 record_name,rr=find(z,'E508_original_source_csv_sha256.json')
 receipt=json.loads(rr)
 row=[r for r in receipt if r['GL']==gl and r['version']==v]
 assert len(row)==1 and sha(b)==row[0]['sha256'] and len(b)==row[0]['size']
 rows=list(csv.DictReader(x for x in b.decode('utf-8-sig').splitlines() if x.strip() and not x.startswith('#')))
 assert len(rows)>100
 # E508 frozen source registration: mother at end=div and exactly two distinct daughter
 # cell IDs with parent_ID==mother cell ID one source frame later.
 track=collections.defaultdict(list)
 for r in rows:track[int(r['cell_ID'])].append(r)
 out=[]
 for cellid,cr in track.items():
  for r in cr:
   if r['type_of_end']!='div':continue
   f=int(r['frame']);daughters=sorted({int(x['cell_ID']) for x in rows if int(x['parent_ID'])==cellid and int(x['frame'])==f+1})
   if len(daughters)!=2:continue
   if not 200<=f<=700:continue
   out.append({'f':f,'cell_id':cellid,'parent_ID':int(r['parent_ID']),'ymin':int(r['bbox_top_px']),'ymax':int(r['bbox_bottom_px']),'daughters':daughters})
 return rows,out

all_reports={};source_checks=[]
for gl,cfg in C.items():
 z8=open_checked(B/'e508_gdrive_readback'/('901_E508_'+cfg['base'].replace('20221220_','')+'_source_conflicts.zip'),cfg['sha508'])
 z12=open_checked(B/('901_E512_'+gl+'_full_original_components.zip'),cfg['sha512'])
 receipt=json.loads(z12.read(cfg['base']+'/'+gl+'/E512_original_member_SHA_CRC_critical_source_audit.json'))
 assert receipt['archive']==cfg['base']+'.zip'
 rec=receipt['complete_source_authentication']
 sms={};stat={};events={};curations={};selected={}
 for v in ['main','cost']:
  sms[v]=sel_sol(original(z12,rec,v,'gurobi_model.sol'))
  stat[v],events[v]=cells(z8,gl,v,cfg)
  curations[v]=original(z12,rec,v,'curation.moma').decode()
  selected[v]=active_div(sms[v])
 assert sms['main'].keys()==sms['cost'].keys()
 good={v:{q['f'] for q in events[v]} for v in ['main','cost']}
 summary={
    'source_GL':gl,'original_E508_archive_sha':cfg['sha508'],'original_E512_archive_sha':cfg['sha512'],
    'selected_solver_Div_all_times':{v:len(selected[v]) for v in ['main','cost']},
    'selected_solver_Div_times_200_700':{v:sum(200<=tvar(k)<=700 for k in selected[v]) for v in ['main','cost']},
    'selected_solver_Div_main_only':[tvar(x) for x in sorted(selected['main']-selected['cost'],key=tvar)],
    'CellStats_qualifying_export_div_times_200_700':{v:len(events[v]) for v in ['main','cost']},
    'CellStats_qualifying_main_only':[t for t in sorted(good['main']-good['cost'])],
    'CellStats_definition':'MoMA mother type_of_end == div and exactly two distinct daughter cell_ID with parent_ID=mother ID at frame+1, restricted frames 200..700; this is NOT original physical independent truth',
 }
 all_reports[gl]=summary
 source_checks.append({'gl':gl,'original_E508_archive_sha256':cfg['sha508'],'original_E512_archive_sha256':cfg['sha512'],'zip_member_CRC_pass':True,'original_source_sol_curation_CRC_SHA_pass':True,'source_CellStats_Sha_pass':True})
 if gl=='Pos5_GL35':
  assert len(selected['main']-selected['cost'])==9
  assert len(events['main'])==13 and len(events['cost'])==4
  assert good['main']-good['cost']=={247,294,399,403,405,409,411,415,622}
  z8.close();z12.close();continue
 assert len(selected['main']-selected['cost'])==17
 assert len(events['main'])==12 and len(events['cost'])==9
 assert good['main']-good['cost']=={241,246,608}
 # Canonical source pruning root: source curations were exported at frame+1.
 pr=[line.strip() for line in curations['main'].splitlines() if line.strip().startswith('PR,')]
 assert len(pr)==1,pr
 pp=[t.strip() for t in pr[0].split(',',3)];assert int(pp[1])==262 and pp[3]=='HypT261T127B150L51R60H1554741365'
 prroot=pp[3]
 rootrows=[r for r in stat['main'] if r['type_of_end']=='pruned' and int(r['frame'])==261 and int(r['bbox_top_px'])==127 and int(r['bbox_bottom_px'])==150]
 assert len(rootrows)==1 and int(rootrows[0]['cell_ID'])==79
 # All original solver selected Map/Div edges across both branches.
 graphs={};outgoing={}
 for v in ['main','cost']:
  graph=collections.defaultdict(list); out=collections.defaultdict(list)
  for var,value in sms[v].items():
   if value<.5:continue
   parts=var.split('_');kind=parts[0]
   if kind.startswith('MapT') and len(parts)==3:dest=parts[2:]
   elif kind.startswith('DivT') and len(parts)==4:dest=parts[2:]
   else:continue
   for d in dest:
    graph[parts[1]].append(d);out[parts[1]].append((var,d))
  graphs[v]=graph;outgoing[v]=out
 seen={prroot};q=[prroot]
 while q:
  node=q.pop()
  for child in graphs['main'].get(node,[]):
   if child not in seen:seen.add(child);q.append(child)
 div_main_only=sorted(selected['main']-selected['cost'],key=tvar)
 subset=[k for k in div_main_only if k.split('_')[1] in seen]
 assert len(subset)==14 and {tvar(x) for x in subset}=={261,274,369,373,377,379,381,386,391,395,398,401,403,405}
 remaining=[k for k in div_main_only if k not in subset]
 assert {tvar(k) for k in remaining}=={241,246,608}
 # Check exact absence from exported mother-two-daughter list (not just filename).
 assert all(tvar(k) not in good['main'] for k in subset)
 assert all(tvar(k) in good['main'] for k in remaining)
 # Traverse one complete uninterrupted original cost-branch selected mapping line.
 mom='HypT246T126B149L48R60H1589543408'
 chain=[]
 for f in range(246,607):
  opts=outgoing['cost'].get(mom,[])
  assert len(opts)==1,(f,mom,len(opts))
  var,child=opts[0]
  assert var.startswith(f'MapT{f}_') and var.split('_')[1]==mom and var.split('_')[2]==child
  chain.append({'gl':gl,'source_frame':f,'source_mother_component':mom,'cost_selected_map_variable':var,'next_component':child,
   'main_selected_same_mapping':bool(sms['main'][var]>.5),'main_selected_div_same_mother':next((k for k in selected['main'] if tvar(k)==f and k.split('_')[1]==mom),'')})
  mom=child
 assert len(chain)==361 and re.match(r'^HypT607',mom)
 # Source all cost-version ASC equality rows, and check if each forced exactly one of this continuous chain.
 allforced=[]
 for line in curations['cost'].splitlines():
  line=line.strip()
  if not line.startswith('ASC,'):continue
  fields=[s.strip() for s in line.split(',',4)];assert len(fields)==5 and float(fields[3])==1
  allforced.append((int(fields[1])-1,fields[4]))
 assert len(allforced)==216
 forced={var for _,var in allforced}
 forcedchain=[r for r in chain if r['cost_selected_map_variable'] in forced]
 assert len(forcedchain)==160 and [r['source_frame'] for r in forcedchain]==list(range(246,406))
 main_different=[r for r in forcedchain if not r['main_selected_same_mapping']]
 assert len(main_different)==15
 assert {r['source_frame'] for r in main_different}=={246,261,274,369,373,377,379,381,386,391,395,398,401,403,405}
 assert all(r['main_selected_div_same_mother'] for r in main_different)
 # Authenticate the real CellStats exported durations, not merely the MIP source variable graph.
 mainpruned=[r for r in stat['main'] if r['cell_ID']=='79']
 costexport=[r for r in stat['cost'] if r['cell_ID']=='76']
 assert len(mainpruned)==15 and sorted(int(r['frame']) for r in mainpruned)==list(range(247,262))
 assert mainpruned[-1]['type_of_end']=='pruned'
 assert len(costexport)==365 and sorted(int(r['frame']) for r in costexport)==list(range(243,608))
 assert costexport[-1]['type_of_end']=='exit'
 # These 14 extra source-main divisions are downstream of a flagged pruned root
 # and absent in the source CellStats list of actual-exported division events.
 lineage={r['source_frame']:r for r in chain}
 for k in subset:
  f=tvar(k)
  assert lineage[f]['main_selected_div_same_mother']==k
 summary.update({
 'main_version_pruning_record_original':pr[0],
 'main_CellStats_source_pruned_parent_id':79,
 'main_CellStats_pruned_track_frames_247_261_count':len(mainpruned),
 'cost_CellStats_unbroken_exported_track_cell_id':76,
 'cost_CellStats_unbroken_exported_track_frames_243_607_count':len(costexport),
 'cost_CellStats_unbroken_exported_track_end_label':'exit',
 'pruned_main_solver_division_selected_only_n':14,
 'pruned_main_solver_division_event_frames':[tvar(x) for x in subset],
 'non_pruned_main_solver_divisions_that_are_CellStats_extra':[tvar(x) for x in remaining],
 'cost_selected_unbroken_Map_start_frame':246,'cost_selected_unbroken_Map_end_frame':607,
 'cost_selected_unbroken_Map_transition_count':len(chain),
 'cost_unbroken_Mapping_forced_ASC_start_frame':246,'cost_unbroken_Mapping_forced_ASC_end_frame':405,
 'cost_unbroken_Mapping_consecutive_ASC_count':len(forcedchain),
 'cost_unbroken_Mapping_ASC_not_accepted_by_main':len(main_different),
 'cost_unbroken_Mapping_ASC_not_accepted_by_main_frames':[r['source_frame'] for r in main_different],
 '14_other_P12_ASC_conflicts_on_same_cost_mapping_chain':True,
 '14_source_main_activated_division_variables_below_PR_root_but_unexported':True,
 'biological_truth':'NOT_IDENTIFIED_SOFTWARE_GRAPH_ONLY'})
 # Detailed original-main Div by original CellStats status and pruned graph ancestry.
 with (OUT/'E520_Pos12_17_original_extra_solver_DIV_vs_CellStats_export.csv').open('w',newline='',encoding='utf8') as fp:
  cols=['gl','frame','original_selected_div_variable','source_main_selected_1','source_cost_selected_0','source_main_PR_root_descendant','source_CellStats_qualified_export','source_mother_component']
  writer=csv.DictWriter(fp,fieldnames=cols);writer.writeheader()
  for k in div_main_only:
   writer.writerow({'gl':gl,'frame':tvar(k),'original_selected_div_variable':k,'source_main_selected_1':True,'source_cost_selected_0':True,
          'source_main_PR_root_descendant':k.split('_')[1] in seen,
          'source_CellStats_qualified_export':tvar(k) in good['main'],
          'source_mother_component':k.split('_')[1]})
 with (OUT/'E520_Pos12_361_cost_selected_unbroken_map_chain.csv').open('w',newline='',encoding='utf8') as fp:
  for r in chain:r.update({'cost_ASC_forced':r['cost_selected_map_variable'] in forced,'source_main_pruning_f261_long_tail':r['source_frame']>=261})
  writer=csv.DictWriter(fp,fieldnames=list(chain[0]));writer.writeheader();writer.writerows(chain)
 with (OUT/'E520_Pos12_160_consecutive_cost_forced_mappings.csv').open('w',newline='',encoding='utf8') as fp:
  writer=csv.DictWriter(fp,fieldnames=list(chain[0]));writer.writeheader();writer.writerows(forcedchain)
 z8.close();z12.close()

(OUT/'E520_original_source_solver_vs_CellStats_reconciliation.json').write_text(json.dumps({'sources':source_checks,'comparisons':all_reports,
 'evidence_note':'MoMA saved optimizer binary assignments and its own CellStats export, NOT independent biological truth',
 'E508_overall_42_vs_30_event_ledger_unchanged':True,'the_source_Div_variable_counts_are_a_distinct_measurement_tier':True,
 'MPS_LP_full_matrix_still_unextracted':True,'independent_human_ground_truth':'NOT_AVAILABLE','publication':'PUBLICATION_NO_GO'},ensure_ascii=False,indent=2),encoding='utf8')
print('E520_PRUNING_EXPORT_SOURCE_RECONCILIATION_PASS',json.dumps({x:{'solver_main':v['selected_solver_Div_times_200_700']['main'],'solver_cost':v['selected_solver_Div_times_200_700']['cost'],'valid_main':v['CellStats_qualifying_export_div_times_200_700']['main'],'valid_cost':v['CellStats_qualifying_export_div_times_200_700']['cost']} for x,v in all_reports.items()},ensure_ascii=False),flush=True)
