#!/usr/bin/env python3
"""901 E521: archived source Gurobi-assignment graph, true continuity/degree subset, minimal Hamming repair.

Scientific tier: exact mixed-integer optimization on an EXPLICIT SUBSET of the
original MoMA MILP constraints. Original MPS/LP is NOT available; other original
constraints (in particular tree path-blocking, crossings, etc.) are omitted.
Therefore the minimizer and objective are NOT the original MoMA optimum, and
this does not classify true biological division.

Frozen source: E512 source-authenticated original full Gurobi .sol,
curation.moma, component_forests.json. Source ZIP SHA, all member ZIP CRC,
and bytes of each used original member verified against source-original hashes.

Run: python e521_source_graph_lowerbound_and_flow.py --base /mnt/data --output /mnt/data/901_E521_20261117/results
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,pathlib,re,zipfile
import numpy as np
from scipy.sparse import coo_matrix,vstack
from scipy.optimize import milp,Bounds,LinearConstraint

P=argparse.ArgumentParser();P.add_argument('--base',default='/mnt/data');P.add_argument('--output',required=True);P.add_argument('--time_limit',type=float,default=20.0);args=P.parse_args()
ROOT=pathlib.Path(args.base);OUT=pathlib.Path(args.output);OUT.mkdir(parents=True,exist_ok=True)
GLS={
 'Pos5_GL35':{'src':'20221220_rplN_GFP_pl','sha':'77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6','forest_sha':''},
 'Pos12_GL22':{'src':'20221220_E08_GFP_pl','sha':'73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2','forest_sha':''}}
def sha(b):return hashlib.sha256(b).hexdigest()
def write_json(name,x):
 (OUT/name).write_text(json.dumps(x,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf-8')
def write_csv(name,rows):
 assert rows
 with (OUT/name).open('w',newline='',encoding='utf-8') as f:
  writer=csv.DictWriter(f,fieldnames=sorted(rows[0]));writer.writeheader();writer.writerows(rows)

def parse_action(k):
 ss=k.split('_');typ=re.match(r'^(Map|Div|Exit|Lys)T(\d+)$',ss[0]);assert typ
 n={'Map':1,'Div':2,'Exit':0,'Lys':0}[typ.group(1)];assert len(ss)==n+2
 return typ.group(1),int(typ.group(2)),ss[1],ss[2:]
def t_of_node(k):return int(re.match(r'^HypT(\d+)T',k).group(1))

def load(gl):
 cfg=GLS[gl];src=cfg['src'];zpath=ROOT/f'901_E512_{gl}_full_original_components.zip';bb=zpath.read_bytes()
 assert sha(bb)==cfg['sha'],('source outer SHA',gl)
 with zipfile.ZipFile(io.BytesIO(bb)) as z:
  assert z.testzip() is None,('ZIP CRC',gl)
  rec=json.loads(z.read(f'{src}/{gl}/E512_original_member_SHA_CRC_critical_source_audit.json'))
  rr=rec['complete_source_authentication']
  def original(v,kind):
   b=z.read(f'{src}/{gl}/{v}__{kind}')
   q=[x for x in rr if x['version']==v and x['kind']==kind]
   assert len(q)==1 and sha(b)==q[0]['original_sha256'] and len(b)==q[0]['uncompressed_bytes'],(gl,v,kind)
   return b
  mains=original('main','gurobi_model.sol').decode();alts=original('cost','gurobi_model.sol').decode()
  cur=original('cost','curation.moma').decode();fore=original('main','component_forests.json')
  f=json.loads(fore)['components'];assert len(f)==(5500 if gl=='Pos5_GL35' else 5765)
  forests={x['stringId']:x for x in f}
  def readsol(s):
   x={}
   for line in s.splitlines():
    if not line.strip() or line.startswith('#'):continue
    k,v=line.split();assert k not in x and v in ('0','1'),(k,v)
    x[k]=int(v)
   return x
  main,alt=readsol(mains),readsol(alts)
  assert set(main)==set(alt)
  forced=[]
  for line in cur.splitlines():
   if not line.strip().startswith('ASC,'):continue
   p=[x.strip() for x in line.strip().split(',',4)]
   assert len(p)==5 and p[3]=='1.0' and p[4] in main and alt[p[4]]==1
   assert int(p[1])==parse_action(p[4])[1]+1
   forced.append(p[4])
  assert len(forced)==(25 if gl=='Pos5_GL35' else 216) and len(set(forced))==len(forced)
  return main,alt,forced,forests,{
   'GL':gl,'original_E512_outer_SHA256':sha(bb),'original_component_forest_SHA256':sha(fore),
   'all_source_member_ZIP_CRC_pass':True,'original_used_member_SHA256_pass':True,
   'original_variable_count':len(main),'hypotheses':len(f)}

def build(gl):
 main,alt,forced,forest,receipt=load(gl)
 names=[k for k in main if k.startswith(('MapT','DivT','ExitT','LysT'))]
 indexed={k:i for i,k in enumerate(names)};N=len(names)
 original=np.asarray([main[k] for k in names],dtype=np.int8)
 archived_alt=np.asarray([alt[k] for k in names],dtype=np.int8)
 edge=[parse_action(k) for k in names]
 bysource=collections.defaultdict(list);bydest=collections.defaultdict(list)
 for i,(typ,t,src,dests) in enumerate(edge):
  assert t_of_node(src)==t and src in forest
  for d in dests:assert t_of_node(d)==t+1 and d in forest
  bysource[src].append(i)
  for d in dests:bydest[d].append(i)
 inn_nodes=set(bydest);out_nodes=set(bysource);allnodes=inn_nodes|out_nodes
 internal={node:i for i,node in enumerate(sorted(q for q in allnodes if 0<t_of_node(q)<1035))}
 outr={node:i+len(internal) for i,node in enumerate(sorted(out_nodes))}
 inr={node:i+len(internal)+len(outr) for i,node in enumerate(sorted(inn_nodes))}
 row=[];col=[];data=[]
 for j,(typ,t,src,dests) in enumerate(edge):
  if src in internal:row.append(internal[src]);col.append(j);data.append(-1.)
  if src in outr:row.append(outr[src]);col.append(j);data.append(1.)
  for d in dests:
   if d in internal:row.append(internal[d]);col.append(j);data.append(1.)
   if d in inr:row.append(inr[d]);col.append(j);data.append(1.)
 R=len(internal)+len(outr)+len(inr)
 A=coo_matrix((np.asarray(data),(row,col)),shape=(R,N),dtype=np.float64).tocsr()
 lower=np.zeros(R);upper=np.r_[np.zeros(len(internal)),np.ones(len(outr)+len(inr))]
 def valid(x):
  vals=A@x
  return int(np.count_nonzero((vals<lower-1e-7)|(vals>upper+1e-7)))
 assert valid(original)==valid(archived_alt)==0
 assert len(names)==(44231 if gl=='Pos5_GL35' else 47879)
 def flow(x):
  incoming=collections.defaultdict(list);outgoing=collections.defaultdict(list)
  for j in np.flatnonzero(x):
   typ,t,src,dests=edge[j]
   outgoing[src].append(names[j])
   for d in dests:incoming[d].append(names[j])
  return incoming,outgoing
 in_main,out_main=flow(original)
 contrad=[k for k in forced if not main[k]]
 assert len(contrad)==(8 if gl=='Pos5_GL35' else 17)
 naïve=original.copy();local=[]
 for var in contrad:
  j=indexed[var];typ,t,src,dest=edge[j];assert typ=='Map' and len(dest)==1
  old=out_main[src];assert len(old)==1
  k=old[0];typ_old,_,old_src,old_tar=parse_action(k)
  assert typ_old in ('Div','Exit') and old_src==src and k!=var
  assert naïve[j]==0 and naïve[indexed[k]]==1
  naïve[j]=1;naïve[indexed[k]]=0
  extra=[x for x in old_tar if x!=dest[0]]
  row={'GL':gl,'f_zero_based':t,'forced_ASC_map':var,'old_selected_outgoing':k,
       'old_selected_outgoing_type':typ_old,'mother_source_id':src,
       'forced_map_target_id':dest[0],'old_target_incoming_count':len(in_main[dest[0]]),
       'other_daughter_id':extra[0] if extra else '',
       'other_daughter_minY':int(min(forest[extra[0]]['yCoordinates'])) if extra else None,
       'other_daughter_pixel_area':len(forest[extra[0]]['xCoordinates']) if extra else None}
  local.append(row)
 inc_n,out_n=flow(naïve)
 imbalance=[]
 for node in sorted(allnodes):
  t=t_of_node(node)
  if t==0 or t==1035:continue
  before=len(in_main[node])-len(out_main[node])
  assert before==0,(gl,'original main mismatch',node,before)
  now=len(inc_n[node])-len(out_n[node])
  if now:
   type_='ORPHAN_OUTGOING_WITHOUT_INCOMING' if now==-1 else 'DOUBLE_INCOMING_ONE_OUTGOING' if now==1 else 'OTHER'
   imbalance.append({'GL':gl,'t_zero_based':t,'candidate_id':node,'in_degree':len(inc_n[node]),
       'out_degree':len(out_n[node]),'imbalance_type':type_,
       'candidate_minY':min(forest[node]['yCoordinates']),
       'candidate_area_pixels':len(forest[node]['xCoordinates']),
       'source_main_incoming':';'.join(in_main[node]),'source_main_outgoing':';'.join(out_main[node]),
       'naive_modified_incoming':';'.join(inc_n[node]),'naive_modified_outgoing':';'.join(out_n[node])})
 # source-specific complete flow and at-most-one constraints validated separately
 assert len(imbalance)==len(contrad) and valid(naïve)==len(contrad)+(gl=="Pos5_GL35")
 typ_counts=collections.Counter(q['imbalance_type'] for q in imbalance)
 assert typ_counts['ORPHAN_OUTGOING_WITHOUT_INCOMING']==(7 if gl=='Pos5_GL35' else 17)
 assert typ_counts['DOUBLE_INCOMING_ONE_OUTGOING']==(1 if gl=='Pos5_GL35' else 0)
 assert all(q['candidate_minY']==120 for q in imbalance if q['imbalance_type']=='ORPHAN_OUTGOING_WITHOUT_INCOMING')

 # Tight, original-file-grounded but RELAXED MIP: all archived possible
 # Map/Div/Exit/Lys actions and original ASC=1; equal in/out at all internal
 # hypothesis times; <=1 incoming and <=1 outgoing per hypothesis.
 # Missing original source constraints (e.g. component-path blocking) keep
 # this as a *lower bound* on changes needed in the complete original ILP.
 L=np.zeros(N);U=np.ones(N);L[[indexed[k] for k in forced]]=1
 c=np.asarray(1-2*original,dtype=np.float64)
 res=milp(c=c,integrality=np.ones(N),bounds=Bounds(L,U),constraints=LinearConstraint(A,lower,upper),
          options={'time_limit':float(args.time_limit),'mip_rel_gap':0.0,'presolve':True,'disp':False})
 assert res.success and res.status==0 and res.mip_gap==0 and res.x is not None,(gl,res.message)
 y=np.rint(res.x).astype(np.int8)
 assert np.max(np.abs(res.x-y))<1e-5 and valid(y)==0
 assert all(y[indexed[k]]==1 for k in forced)
 bound=int(np.abs(y-original).sum())
 assert abs((int(original.sum())+res.fun)-bound)<1e-6
 assert abs((int(original.sum())+res.mip_dual_bound)-bound)<1e-6
 assert bound==(29 if gl=='Pos5_GL35' else 72)
 actual=int(np.abs(archived_alt-original).sum())
 assert actual==(35 if gl=='Pos5_GL35' else 184)
 assert bound<=actual
 flips=[]
 for j in np.flatnonzero(y!=original):
  typ,t,src,dests=edge[j]
  flips.append({'GL':gl,'frame_t':t,'action':names[j], 'type':typ,'main_value':int(original[j]),'relaxed_value':int(y[j]),
      'cost_version_value':int(archived_alt[j]),'forced_ASC':names[j] in forced})

 # A separate ORIGINAL-MoMA source component-forest path-block audit.
 # Active hypotheses in every time frame cannot be both ancestor and
 # descendant in the same original component tree. This check is kept
 # outside the relaxation and used as an explicit missing-constraint witness.
 parents={sid:rec['parentStringId'] for sid,rec in forest.items()}
 def path_conflicts(x):
  active={edge[j][2] for j in np.flatnonzero(x)}
  bad=set()
  for node in active:
   p=parents[node];seen=set()
   while p!='NA':
    assert p in parents and p not in seen,(gl,p)
    seen.add(p)
    if p in active:
     bad.add((p,node))
    p=parents[p]
  return sorted(bad)
 path0=path_conflicts(original);path1=path_conflicts(archived_alt);path2=path_conflicts(y)
 assert not path0 and not path1
 tree_rows=[{'GL':gl,'ancestor_id':a,'descendant_id':b,'t':t_of_node(a)} for a,b in path2]

 # Second relaxation adds original source component-forest tree exclusions.
 # If two hypotheses are ancestor/descendant in one original source forest,
 # they may not both have active outgoing assignments at the same source time.
 # We impose pairwise exclusions for ALL ancestor pairs, not just parent links.
 # These remain a relaxation: full original MPS crossings/cell counts and
 # official objective are NOT included.
 pair_rows=[]
 for node in sorted(bysource):
  anc=parents[node];seen=set()
  while anc!='NA':
   assert anc in parents and anc not in seen
   seen.add(anc)
   if anc in bysource:pair_rows.append((anc,node))
   anc=parents[anc]
 p_r=[];p_c=[]
 for k,(anc,node) in enumerate(pair_rows):
  for idx in bysource[anc]:p_r.append(k);p_c.append(idx)
  for idx in bysource[node]:p_r.append(k);p_c.append(idx)
 B=coo_matrix((np.ones(len(p_r)),(p_r,p_c)),shape=(len(pair_rows),N)).tocsr()
 assert not np.any(B@original>1+1e-9) and not np.any(B@archived_alt>1+1e-9)
 A_tree=vstack([A,B],format='csr')
 low_tree=np.r_[lower,np.zeros(len(pair_rows))]
 upper_tree=np.r_[upper,np.ones(len(pair_rows))]
 res_tree=milp(c=c,integrality=np.ones(N),bounds=Bounds(L,U),
      constraints=LinearConstraint(A_tree,low_tree,upper_tree),
      options={'time_limit':float(args.time_limit),'mip_rel_gap':0.0,'presolve':True,'disp':False})
 assert res_tree.success and res_tree.status==0 and res_tree.mip_gap==0 and res_tree.x is not None,(gl,'tree repair incomplete',res_tree.message)
 y_tree=np.rint(res_tree.x).astype(np.int8)
 assert np.max(np.abs(res_tree.x-y_tree))<1e-6
 assert valid(y_tree)==0 and not np.any(B@y_tree>1+1e-8)
 assert not path_conflicts(y_tree)
 assert all(y_tree[indexed[k]]==1 for k in forced)
 tree_bound=int(np.abs(y_tree-original).sum());assert tree_bound>=bound and tree_bound<=actual
 assert abs(int(original.sum())+res_tree.fun-tree_bound)<1e-6
 assert abs(int(original.sum())+res_tree.mip_dual_bound-tree_bound)<1e-6
 tree_flips=[]
 for j in np.flatnonzero(y_tree!=original):
  typ,t,src,dests=edge[j]
  tree_flips.append({'GL':gl,'frame_t':t,'action':names[j], 'type':typ,'main_value':int(original[j]),'tree_relaxed_value':int(y_tree[j]),
     'cost_version_value':int(archived_alt[j]),'forced_ASC':names[j] in forced})
 tree_vector_sha=sha('\n'.join(f'{names[i]} {int(y_tree[i])}' for i in range(N)).encode())

 # Input receipts, proof-carrying vector checksum for deterministic rerun.
 sha_vec=sha('\n'.join(f'{names[i]} {int(y[i])}' for i in range(N)).encode())
 summary={'GL':gl,'source':receipt,'original_action_variables':N,'original_all_gurobi_variables':len(main),
   'original_ASC_count':len(forced),'original_ASC_conflicts':len(contrad),'archived_main_graph_violations':0,
   'archived_cost_graph_violations':0,'forced_mapping_plus_old_outgoing_cancel_binary_toggles':len(contrad)*2,
   'naive_flow_violations':len(imbalance),'naive_orphan_top_y120':typ_counts['ORPHAN_OUTGOING_WITHOUT_INCOMING'],
   'naive_double_incoming':typ_counts['DOUBLE_INCOMING_ONE_OUTGOING'],
   'relaxed_model':'ASC=1 + original available Map/Div/Exit/Lys actions + per-hypothesis in-out continuity + each in/out degree <=1; excludes original tree path blocking, crossings, hypothesis budget etc',
   'min_hamming_lower_bound_for_original_FULL_MIP':bound,
   'relaxed_mip_optimal':True,'relaxed_mip_gap':float(res.mip_gap),'relaxed_mip_dual_bound':float(res.mip_dual_bound+original.sum()),
   'archived_cost_vs_main_assignment_flips':actual,
   'relaxed_optimum_flip_0to1':int(((original==0)&(y==1)).sum()),
   'relaxed_optimum_flip_1to0':int(((original==1)&(y==0)).sum()),
   'relaxed_optimum_source_pathblocking_violation_pairs':len(tree_rows),
   'archived_main_source_pathblocking_violation_pairs':0,'archived_cost_source_pathblocking_violation_pairs':0,
   'relaxed_solution_SHA256_over_sorted_action_lines':sha_vec,
   'tree_pairwise_original_source_ancestor_exclusions':len(pair_rows),
   'min_tree_constrained_hamming_lower_bound':tree_bound,
   'tree_relaxed_mip_gap':float(res_tree.mip_gap),
   'tree_relaxed_mip_dual_bound':float(int(original.sum())+res_tree.mip_dual_bound),
   'tree_relaxed_solution_SHA256_over_sorted_action_lines':tree_vector_sha,
   'true_original_mps_lp_bytes_available':False,'original_Gurobi_full_model_reoptimized':False,
   'independent_biological_division_ground_truth':'UNRESOLVED','molecular_bit_per_sec_per_joule':'NOT_IDENTIFIED',
   'publication':'PUBLICATION_NO_GO'}
 return summary,local,imbalance,flips,tree_rows,tree_flips

summary=[];local=[];imbalance=[];flips=[];tree=[];tree_flips_all=[]
for gl in GLS:
 a,b,c,d,e,f=build(gl)
 summary.append(a);local.extend(b);imbalance.extend(c);flips.extend(d);tree.extend(e);tree_flips_all.extend(f)
 print('E521_GLU_PASS',json.dumps({k:a[k] for k in ['GL','original_action_variables','original_ASC_conflicts','naive_flow_violations','naive_orphan_top_y120','naive_double_incoming','min_hamming_lower_bound_for_original_FULL_MIP','archived_cost_vs_main_assignment_flips','relaxed_optimum_source_pathblocking_violation_pairs','min_tree_constrained_hamming_lower_bound']},ensure_ascii=False),flush=True)
assert len(local)==25 and len(imbalance)==25 and len(flips)==101
write_csv('E521_25_original_ASC_forced_map_substitution_ledger.csv',sorted(local,key=lambda q:(q['GL'],q['f_zero_based'])))
write_csv('E521_25_naive_graph_continuity_violations.csv',sorted(imbalance,key=lambda q:(q['GL'],q['t_zero_based'])))
write_csv('E521_101_relaxed_minimum_hamming_assignment_changes.csv',sorted(flips,key=lambda q:(q['GL'],q['frame_t'],q['action'])))
write_csv('E521_tree_constrained_minimum_hamming_assignment_changes.csv',sorted(tree_flips_all,key=lambda q:(q['GL'],q['frame_t'],q['action'])))
if tree:write_csv('E521_missing_original_component_tree_pathblocking_witness.csv',sorted(tree,key=lambda q:(q['GL'],q['t'],q['ancestor_id'])))
else:(OUT/'E521_missing_original_component_tree_pathblocking_witness.csv').write_text('No nested active source hypotheses found in graph-only optimal repairs; other original MPS constraints still omitted.\n',encoding='utf8')
write_json('E521_source_inference_and_lower_bound_certificate.json',{'stage':'E521','results':summary,'totals':{
 'original_ASCs_forced_violated':25,'naive_selected_action_binary_flips':50,
 'naive_remaining_flow_violations':25,'naive_orphan_top_daughters':24,'naive_double_incoming_targets':1,
 'sum_relaxed_MIP_Hamming_lower_bounds':sum(v['min_hamming_lower_bound_for_original_FULL_MIP'] for v in summary),
 'sum_archived_cost_vs_main_action_changes':sum(v['archived_cost_vs_main_assignment_flips'] for v in summary),
 'sum_tree_constrained_MIP_Hamming_lower_bounds':sum(v['min_tree_constrained_hamming_lower_bound'] for v in summary)},
 'warning':'A lower bound on original full MIP assignment-variable changes, NOT original MPS reoptimization nor cellular ground truth.'})
(OUT/'E521_scientific_limits.txt').write_text('E521 is a strictly source-authenticated graph-relaxation MIP, not the real complete MPS: it omits MoMA component-tree path blocking, crossing constraints, selection criteria, and source cost objective. Objective is minimum count of assignment bits changed from the SAVED STANDARD solution, not source assignment cost or biophysical energy. All 12 actual microbial divisions are UNRESOLVED.\n',encoding='utf8')
print('E521_GRAPH_LOWERBOUND_CERTIFIED',json.dumps({'sum_lower_bound':sum(v['min_hamming_lower_bound_for_original_FULL_MIP'] for v in summary),'main_source_violations':0,'cost_source_violations':0,'tree_pathblocking_relaxed_witnesses':len(tree),'tree_constrained_lower_bound':sum(v['min_tree_constrained_hamming_lower_bound'] for v in summary)},ensure_ascii=False),flush=True)
