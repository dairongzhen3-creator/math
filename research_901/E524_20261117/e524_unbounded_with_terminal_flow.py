#!/usr/bin/env python3
"""901 E524 TERMINAL-FIXED: uncapped original-cost minimization in source-code-derived MoMA subset of constraints.

EVIDENCE: archived original MoMA component forests and official saved Gurobi assignments,
backed by E512 ZIP SHA, ZIP CRC, original member SHA/size; the published MoMA source
GrowthlaneTrackingILP.addCrossingConstraints and AdvancedComponent.getOrdinalValue.
This reconstructs the crossing *part* of MoMA's constraints, NOT the complete original MPS.
The original MPS and independent microscopy ground truth remain unavailable.
"""
import argparse,collections,csv,hashlib,io,json,math,pathlib,re,time,zipfile
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import coo_matrix,vstack,csr_matrix
p=argparse.ArgumentParser();p.add_argument('--base',default='/mnt/data');p.add_argument('--output',default='/mnt/data/901_E523_probe/results');p.add_argument('--run-milp',action='store_true');p.add_argument('--time-limit',type=int,default=55);p.add_argument('--gl',default='Pos12_GL22',choices=['Pos5_GL35','Pos12_GL22']);args=p.parse_args()
ROOT=pathlib.Path(args.base);OUT=pathlib.Path(args.output);OUT.mkdir(parents=True,exist_ok=True)
SHA={'Pos5_GL35':'77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6','Pos12_GL22':'73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2'}
SRC={'Pos5_GL35':'20221220_rplN_GFP_pl','Pos12_GL22':'20221220_E08_GFP_pl'}
def sha(b):return hashlib.sha256(b).hexdigest()
def parse(k):
 s=k.split('_');m=re.fullmatch(r'(Map|Div|Exit|Lys)T(\d+)',s[0]);assert m
 assert len(s)=={'Map':3,'Div':4,'Exit':2,'Lys':2}[m.group(1)]
 return m.group(1),int(m.group(2)),s[1],s[2:]
def node_time(k):return int(re.match(r'^HypT(\d+)T',k).group(1))
def write_json(name,obj):
 (OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+'\n',encoding='utf8')
def write_csv(name,rows):
 if not rows:return (OUT/name).write_text('no_rows\n',encoding='utf8')
 with (OUT/name).open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in rows))));w.writeheader();w.writerows(rows)
def load(gl):
 src=SRC[gl];path=ROOT/f'901_E512_{gl}_full_original_components.zip';bb=path.read_bytes();assert sha(bb)==SHA[gl],gl
 # Verify source mm.properties crossing feature flags for both original analyses.
 p511=ROOT/f'901_E511_{gl}_original_moma_sidecars.zip';raw511=p511.read_bytes()
 assert sha(raw511)=={'Pos5_GL35':'1dad2aee7728757330fef6f16127faa5b5ab1648b78080db78a24621b2ffa596',
                    'Pos12_GL22':'be8695b97cb1a2c946a8f667b201a40574e064e8b1b203124a40cb4bca5adcba'}[gl]
 settings={}
 with zipfile.ZipFile(io.BytesIO(raw511)) as zz:
  assert zz.testzip() is None
  audit=json.loads(zz.read(f'{src}/E511_all_original_sidecar_source_findings.json'))
  for ver,short in [('main','main__2__mm.properties'),('cost_function','cost_function__5__mm.properties')]:
   props=zz.read(f'{src}/{gl}/{short}')
   receipt=[x for x in audit['sources'] if x['version']==ver and x['path'].endswith('/mm.properties')]
   assert len(receipt)==1 and sha(props)==receipt[0]['sha256'] and len(props)==receipt[0]['bytes']
   assert 'USE_FEATURE_CROSSING_CONSTRAINTS=1' in props.decode('utf8').splitlines(),(gl,ver)
   settings[ver]={'original_mm_properties_sha256':sha(props),'USE_FEATURE_CROSSING_CONSTRAINTS':1}
 with zipfile.ZipFile(io.BytesIO(bb)) as z:
  assert z.testzip() is None
  prefix=src+'/'+gl+'/'
  rc=json.loads(z.read(prefix+'E512_original_member_SHA_CRC_critical_source_audit.json'))['complete_source_authentication']
  def original(version,kind):
   b=z.read(prefix+version+'__'+kind)
   x=[a for a in rc if a['version']==version and a['kind']==kind]
   assert len(x)==1 and sha(b)==x[0]['original_sha256'] and len(b)==x[0]['uncompressed_bytes']
   return b
  f={x['stringId']:x for x in json.loads(original('main','component_forests.json'))['components']}
  def sol(v):
   s=original(v,'gurobi_model.sol').decode()
   return {ls.split()[0]:int(ls.split()[1]) for ls in s.splitlines() if ls.strip() and not ls.startswith('#')}
  main,alt=sol('main'),sol('cost');assert set(main)==set(alt)
  c={r['id']:float(r['cost']) for r in csv.DictReader(ln for ln in original('cost','assignment_costs.csv').decode().splitlines() if ln.strip() and not ln.startswith('#'))}
  cur=original('cost','curation.moma').decode()
  forced=[]
  for line in cur.splitlines():
   if line.lstrip().startswith('ASC,'):
    q=[a.strip() for a in line.split(',',4)];assert len(q)==5 and q[3]=='1.0' and alt[q[4]]==1
    forced.append(q[4])
  assert len(forced)==(25 if gl=='Pos5_GL35' else 216)
 return f,main,alt,c,forced,{'gl':gl,'original_E512_SHA256':sha(bb),'original_E511_SHA256':sha(raw511),'verified_original_crossing_feature_flags':settings,'original_members_and_crc_pass':True,'original_member_sha256_pass':True}
def ordinal_rank(f):
 # Source AdvancedComponentForest sorts root and sibling components by first moment y;
 # getRankRelativeToLeafComponent counts all leaves spatially below; ordinal=2^rank.
 byframe=collections.defaultdict(list)
 for k,v in f.items():byframe[v['frameNumber']].append(k)
 ords,span,nleaf={}, {},{};order_ok=[]
 for t,ids in sorted(byframe.items()):
  def mean_y(k):return sum(f[k]['yCoordinates'])/len(f[k]['yCoordinates'])
  roots=sorted([k for k in ids if f[k]['parentStringId']=='NA'],key=mean_y)
  seen=[]
  def dfs(k):
   children=f[k]['childrenStringIds']
   if children:
    assert children==sorted(children,key=mean_y),('source children not center sorted',t,k)
    for cc in children:dfs(cc)
   else:seen.append(k)
  for root in roots:dfs(root)
  assert set(seen)=={k for k in ids if not f[k]['childrenStringIds']}
  assert seen==sorted(seen,key=mean_y),('source leaf centroid global order mismatch',t)
  assert seen==sorted(seen,key=lambda k:min(f[k]['yCoordinates'])),('source leaf bounding order mismatch',t)
  n=len(seen);nleaf[t]=n;pos={k:j for j,k in enumerate(seen)}
  def rec(k):
   if k in ords:return ords[k],span[k]
   children=f[k]['childrenStringIds']
   if not children:
    i=pos[k];o=2**(n-1-i);seg=(i,i)
   else:
    subs=[rec(cc) for cc in children];o=sum(q[0] for q in subs);seg=(min(q[1][0] for q in subs),max(q[1][1] for q in subs))
    assert sorted(range(len(subs)),key=lambda x:subs[x][1][0])==list(range(len(subs)))
   ords[k]=o;span[k]=seg;return o,seg
  for k in ids:rec(k)
  assert set(ids)<=ords.keys()
  order_ok.append({'frame':t,'leaf_count':n,'bigM_for_constraints_into_this_frame':2**n})
 return ords,span,nleaf,order_ok

def load_vectors(gl,main):
 result={}
 for tag in ('tie','radius'):
  path=ROOT/f'901_E522_20261117/results/E522_{gl}_{tag}_changed_actions.csv'
  assert path.exists(),path
  v=main.copy()
  with path.open() as f:
   for q in csv.DictReader(f):
    k=q['action'];assert v[k]==int(q['main']);v[k]=int(q['new_partial'])
  result[tag]=v
 return result

def analyze(gl):
 forest,main,alt,costs,forced,receipt=load(gl)
 ords,span,nleaf,leafcheck=ordinal_rank(forest)
 names=[k for k in main if k.startswith(('MapT','DivT','ExitT','LysT'))]
 assert len(names)==(44231 if gl=='Pos5_GL35' else 47879) and set(costs)==set(names)
 idx={k:i for i,k in enumerate(names)};N=len(names)
 edge=[parse(k) for k in names]
 # Complete candidate outgoing crossing terms; only Map/Div contribute to crossing.
 outs=collections.defaultdict(list);mt=collections.defaultdict(list);all_by_t=collections.defaultdict(list)
 for j,(typ,t,s,d) in enumerate(edge):
  outs[s].append(j)
  if typ in ('Map','Div'):
   mt[s].append((j,sum(ords[x] for x in d)))
   all_by_t[t].append(s)
 for t in all_by_t:all_by_t[t]=sorted(set(all_by_t[t]))
 # The source rule, rewritten in <= form, for each candidate H at time t:
 # sum_{Map/Div from H} (M(t+1)-ordinal(dest))*x +
 # sum_{Map/Div from sources below H} ordinal(dest)*x <= M(t+1).
 # "below" exactly matches lower-neighbor-subtree rule in published AdvancedComponent.
 cr,cc,cv,ub,records=[],[],[],[],[]
 for t,ids in sorted(all_by_t.items()):
  M=2**nleaf[t+1]
  for h in ids:
   below=[b for b in ids if span[b][0]>span[h][1]]
   if not mt[h]:continue
   rowid=len(ub); ub.append(M);records.append((t,h,M,below))
   for j,o in mt[h]:cr.append(rowid);cc.append(j);cv.append(M-o)
   for b in below:
    for j,o in mt[b]:cr.append(rowid);cc.append(j);cv.append(o)
 X=coo_matrix((np.asarray(cv,float),(cr,cc)),shape=(len(ub),N)).tocsr()
 bigM=np.asarray(ub,float)
 # Exact exit boundary constraints from published MoMA ExitAssignment.addConstraintsToILP():
 # |Hup(h)|*Exit(h) + sum_{upper u < h, all a != Exit(u)} a <= |Hup(h)|.
 # LpUtils.getHup is defined by strict first-moment Y position, not by tree rank;
 # the static ConfigurationManager.DISABLE_EXIT_CONSTRAINTS is FALSE in this official source.
 out_all=collections.defaultdict(list)
 for j,(typ,t,source,dest) in enumerate(edge):out_all[source].append(j)
 bysrcframe=collections.defaultdict(list)
 for sid,rec in forest.items():bysrcframe[rec['frameNumber']].append(sid)
 means={sid:sum(rec['yCoordinates'])/len(rec['yCoordinates']) for sid,rec in forest.items()}
 er,ec,ev,eupper,exitmetadata=[],[],[],[],[]
 for j,(typ,t,src,dest) in enumerate(edge):
  if typ!='Exit':continue
  hup=[sid for sid in bysrcframe[t] if means[sid]<means[src]]
  if not hup:continue
  m=len(hup)
  rr=len(eupper);nonexit=[]
  for sid in hup:
   nonexit.extend(k for k in out_all[sid] if edge[k][0]!='Exit')
  if not nonexit:continue
  er.append(rr);ec.append(j);ev.append(float(m))
  for k in nonexit:er.append(rr);ec.append(k);ev.append(1.0)
  eupper.append(float(m));exitmetadata.append((t,src,names[j],m))
 E=coo_matrix((np.asarray(ev,float),(er,ec)),shape=(len(eupper),N)).tocsr()
 ex_bound=np.asarray(eupper,float)
 assert X.min()>=0
 def vzero(sol):return np.asarray([sol[k] for k in names],dtype=np.int8)
 vectors={'main':vzero(main),'cost':vzero(alt)}
 for tag,vec in load_vectors(gl,main).items():vectors[tag]=vzero(vec)
 stats={};source_details=[];exit_details=[]
 for tag,v in vectors.items():
  val=X@v;bad=np.where(val>bigM+1e-7)[0]
  evval=E@v;exitbad=np.where(evval>ex_bound+1e-7)[0]
  stats[tag]={'violating_published_crossing_rows':len(bad),'evaluated_source_code_equivalent_rows':len(ub),'max_violation':float(max([0]+list(val[bad]-bigM[bad]))),
   'violating_published_exit_rows':len(exitbad),'evaluated_published_exit_rows':len(eupper),'max_exit_violation':float(max([0]+list(evval[exitbad]-ex_bound[exitbad])))}
  for eidx in exitbad:
   t,src,action,m=exitmetadata[int(eidx)]
   exit_details.append({'GL':gl,'source_vector':tag,'frame_zero_based':t,'source_hypothesis':src,'exit_action':action,
                        'exit_source_meanY':means[src],'Hup_count':m,'exit_constraint_lhs':float(evval[eidx]),'RHS':float(ex_bound[eidx]),
                        'excess':float(evval[eidx]-ex_bound[eidx])})
  for rr in bad:
   t,h,M,below=records[int(rr)];selected=names
   focus=[(names[j],o) for j,o in mt[h] if v[j]]
   selected_below=[(names[j],o) for b in below for j,o in mt[b] if v[j]]
   assert len(focus)<=1
   source_details.append({'GL':gl,'source_vector':tag,'frame_zero_based':t,'hypothesis':h,'source_Ymin':min(forest[h]['yCoordinates']), 'bigM':M,
     'crossing_row_lhs':float(val[rr]),'violation_excess':float(val[rr]-M),
     'active_interest':focus[0][0] if focus else None,'active_interest_target_ordinal':focus[0][1] if focus else None,
     'below_selected_actions':';'.join(k for k,_ in selected_below),'below_destination_ordinal_sum':sum(o for _,o in selected_below)})
 assert stats['main']['violating_published_crossing_rows']==stats['cost']['violating_published_crossing_rows']==0
 assert stats['main']['violating_published_exit_rows']==stats['cost']['violating_published_exit_rows']==0
 write_csv('E523_crossing_violation_witnesses.csv',source_details)
 write_csv('E523_exit_violation_witnesses.csv',exit_details)
 out={'scope':'Crossing constraints reconstructed from published MoMA GrowthlaneTrackingILP + AdvancedComponent source and original E512 component forest / solver vectors; not full MPS',
      'GL':gl,'source_receipt':receipt,'model_flags_source_mm_properties': 'Verified for both standard and alternate original mm.properties via E511 original-member SHA receipts',
      'all_published_source_crossing_rows_count':len(ub),'all_published_source_exit_rows_count':len(eupper),'total_source_binary_actions':N,'crossing_row_nnz':int(X.nnz),'exit_row_nnz':int(E.nnz),
      'original_leaves_centroid_sorted_bbox_sorted_and_DFS_sorted_in_all_1036_frames':True,
      'stats':stats,'independent_biological_truth':'UNRESOLVED'}
 if args.run_milp:
  # partial original constraints: ASC + all available graph edges flow + source forest path blocking + NEW crossing equations.
  # True original model may include other constraints (e.g. auxiliary/cellcount/segmentation). No full MPS.
  orig=vectors['main'];arch=vectors['cost'];
  bysource=collections.defaultdict(list);bydest=collections.defaultdict(list)
  for j,(typ,t,s,dests) in enumerate(edge):
   bysource[s].append(j)
   for d in dests:bydest[d].append(j)
  internal=sorted(x for x in set(bysource)|set(bydest) if 0<node_time(x)<=1035)
  R=len(internal)+len(bysource)+len(bydest)
  rowix={x:i for i,x in enumerate(internal)}
  rowout={x:i+len(internal) for i,x in enumerate(bysource)}
  rowin={x:i+len(internal)+len(bysource) for i,x in enumerate(bydest)}
  rr=[];cc=[];vv=[]
  for j,(typ,t,s,dests) in enumerate(edge):
   if s in rowix:rr.append(rowix[s]);cc.append(j);vv.append(-1)
   rr.append(rowout[s]);cc.append(j);vv.append(1)
   for d in dests:
    if d in rowix:rr.append(rowix[d]);cc.append(j);vv.append(1)
    rr.append(rowin[d]);cc.append(j);vv.append(1)
  A=coo_matrix((np.asarray(vv,float),(rr,cc)),shape=(R,N)).tocsr()
  upper=np.r_[np.zeros(len(internal)),np.ones(len(bysource)+len(bydest))]
  lower=np.zeros(R)
  pairs=[]
  for n in bysource:
   par=forest[n]['parentStringId'];seen=set()
   while par!='NA':
    assert par not in seen and par in forest
    seen.add(par)
    if par in bysource:pairs.append((par,n))
    par=forest[par]['parentStringId']
  rr=[];cc=[]
  for row,(an,ch) in enumerate(pairs):
   for j in bysource[an]:rr.append(row);cc.append(j)
   for j in bysource[ch]:rr.append(row);cc.append(j)
  B=coo_matrix((np.ones(len(rr),float),(rr,cc)),shape=(len(pairs),N)).tocsr()
  tree_ub=np.ones(len(pairs))
  assert not np.any(B@orig>tree_ub) and not np.any(B@arch>tree_ub)
  assert not np.any(A@orig>upper) and not np.any(A@arch>upper)
  csource=np.array([costs[k] for k in names],float)
  forced_idx=[idx[k] for k in forced];L=np.zeros(N);L[forced_idx]=1;U=np.ones(N)
  hamming_sign=np.array(1-2*orig,float);hconst=int(orig.sum());
  M=vstack([A,B,X,E],format='csr')
  lo=np.r_[lower,np.zeros(len(pairs)),np.full(len(bigM),-np.inf),np.full(len(ex_bound),-np.inf)]
  hi=np.r_[upper,tree_ub,bigM,ex_bound]
  assert len(lo)==M.shape[0]
  assert M.shape[0]>10000
  print('E524_UNCAPPED_MILP_PREP',json.dumps({'GL':gl,'N':N,'rows':M.shape[0],'nnz':M.nnz,'crossing_rows':X.shape[0],'exit_rows':E.shape[0]}),flush=True)
  solver={}
  # Optimize minimal changes over all crossing rows (global validated bound).
  # The same model with source costs can test whether E522 cost optimum survives.
  for typ in ('cost_uncapped',):
   cvec=csource
   res=milp(c=cvec,integrality=np.ones(N),bounds=Bounds(L,U),constraints=LinearConstraint(M,lo,hi),
           options={'time_limit':float(args.time_limit),'mip_rel_gap':0,'presolve':True,'disp':False})
   d={'status':res.status,'message':res.message,'success':res.success,'mip_gap':None if res.mip_gap is None else float(res.mip_gap),
      'mip_dual_bound':None if res.mip_dual_bound is None else float(res.mip_dual_bound),'walltime_limit':args.time_limit}
   if res.x is not None:
    vec=np.rint(res.x).astype('int8');d.update({'hamming':int(np.abs(vec-orig).sum()),'archived_cost':float(csource@vec),
          'primal_fun':float(res.fun),'crossing_violations':int(np.count_nonzero(X@vec>bigM+1e-7)),
          'exit_violations':int(np.count_nonzero(E@vec>ex_bound+1e-7)),'tree_violations':int(np.count_nonzero(B@vec>1+1e-7)),'flow_violations':int(np.count_nonzero(A@vec>upper+1e-7)),
          'forced_violations':int(sum(vec[j]!=1 for j in forced_idx))})
    assert d['crossing_violations']==d['exit_violations']==d['tree_violations']==d['flow_violations']==d['forced_violations']==0
    if res.success:
     with (OUT/f'E524_{gl}_{typ}_source_crossing_exit_changes.csv').open('w',newline='',encoding='utf8') as f:
      w=csv.DictWriter(f,fieldnames=['action','source_time','main','saved_alt','new','source_cost']);w.writeheader()
      for j in np.flatnonzero(vec!=orig):w.writerow({'action':names[j],'source_time':edge[j][1],'main':int(orig[j]),'saved_alt':int(arch[j]),'new':int(vec[j]),'source_cost':float(csource[j])})
   solver[typ]=d;print('E524_UNCAPPED_MILP_RESULT',json.dumps({'GL':gl,'type':typ,**d}),flush=True)
  out['milp_partial_source_crossing_exit_uncapped']=solver
 write_json('E524_'+gl+'_source_crossing_exit_uncapped_audit.json',out)
 print('E524_UNCAPPED_SOURCE_AUDIT_DONE',json.dumps({'GL':gl,'baseline_crossing':{k:v['violating_published_crossing_rows'] for k,v in stats.items()},'baseline_exit':{k:v['violating_published_exit_rows'] for k,v in stats.items()},'number_rows':len(ub),'has_MILP':args.run_milp}),flush=True)

analyze(args.gl)
