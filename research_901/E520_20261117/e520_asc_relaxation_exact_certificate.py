#!/usr/bin/env python3
"""901 E520 source-authenticated certificate of minimum ASC release for fixed standard MoMA solution.

E520 is NOT a reoptimized Gurobi/HiGHS counterfactual and does NOT provide biological truth.
Reads original official-ZIP-authenticated E511/512/512b material; every source ZIP and
used original member checked before conclusions. Never accesses the private blind key.
"""
import argparse,csv,hashlib,io,json,pathlib,re,zipfile,collections,statistics

P=argparse.ArgumentParser()
P.add_argument('--base',default='/mnt/data')
P.add_argument('--output',required=True)
a=P.parse_args()
B=pathlib.Path(a.base); OUT=pathlib.Path(a.output);OUT.mkdir(parents=True,exist_ok=True)
CFG={
 'Pos5_GL35': {'id':'20221220_rplN_GFP_pl','central':'7f55aaa61d7025260773c33aa662e13ae2b2266dea300384f9e09a1667dc5ac6','count':193307,
     'sha': {'E511':'1dad2aee7728757330fef6f16127faa5b5ab1648b78080db78a24621b2ffa596','E512':'77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6','E512b':'efe81c19a76117c0f9bebe58b8e4525f93c407465312d0402ce7096e9866b384'}},
 'Pos12_GL22':{'id':'20221220_E08_GFP_pl','central':'5bd93ce604256d00233d4c5bf91865205b44b0d7fc2f52fefe0f68ef0e23e94c','count':141104,
     'sha': {'E511':'be8695b97cb1a2c946a8f667b201a40574e064e8b1b203124a40cb4bca5adcba','E512':'73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2','E512b':'4893c521768bdbcc026e7f6d997b1266f9e70e6869b9bf08e4b429a8aa7ecca9'}},
}
MAN_SHA='f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f'

def sha(b):return hashlib.sha256(b).hexdigest()
def archive(path,expected):
 raw=path.read_bytes();assert sha(raw)==expected,(path,'OUTER_SHA');z=zipfile.ZipFile(io.BytesIO(raw))
 assert z.testzip() is None,(path,'ZIP_CRC');return z

def member(z,suffix):
 matched=[x for x in z.namelist() if x.endswith(suffix)]
 assert len(matched)==1,(suffix,matched)
 return matched[0],z.read(matched[0])

def original(z,receipts,version,name):
 path,data=member(z,version+'__'+name)
 rr=[r for r in receipts if r['kind']==name and r['version']==version]
 assert len(rr)==1,(name,version)
 r=rr[0];assert sha(data)==r['original_sha256'] and len(data)==r['uncompressed_bytes']
 return data

def solve(s):
 x={}
 for line in s.decode().splitlines():
  if not line.strip() or line.startswith('#'):continue
  p=line.split()
  assert len(p)==2,p
  name,val=p;assert name not in x
  x[name]=float(val)
 assert x and all(-1e-9<=v<=1+1e-9 for v in x.values())
 return x

violations=[];groups={};all_asc=[];extraction=[];receipts=[]
for gl,cfg in CFG.items():
 source=cfg['id'];url_root=source+'/'+gl+'/'
 z1=archive(B/'e511_gdrive_readback'/('901_E511_'+gl+'_original_moma_sidecars.zip'),cfg['sha']['E511'])
 z2=archive(B/('901_E512_'+gl+'_full_original_components.zip'),cfg['sha']['E512'])
 z3=archive(B/('901_E512b_'+gl+'_exact_div_to_map.zip'),cfg['sha']['E512b'])
 receipts.extend([{'gl':gl,'stage':v,'outer_sha256':cfg['sha'][v],'all_crc_pass':True} for v in ['E511','E512','E512b']])
 lstname,lstbytes=member(z1,'_entire_GL_official_members_listing.json')
 lst=json.loads(lstbytes);assert lst['gl']==gl
 manifest=[]
 for ver,location in [('main','20221220_analysis/track_data__20221220_analysis'),('cost','cost_function_analysis_1/track_data__cost_function_analysis_1')]:
  for extension in ('mps','lp'):
   source_path=url_root+location+'/gurobi_model.'+extension
   assert lst['filename_list'].count(source_path)==1,(gl,ver,extension)
   manifest.append({'gl':gl,'source_archive':source+'.zip','analysis':ver,'format':extension,'official_zip_member_path':source_path,
                    'source_zip_central_directory_SHA256':cfg['central'],'original_central_member_count':cfg['count'],
                    'original_member_crc_sha_NOT_YET_CHECKED':True,'model_matrix_available_locally':False})
 extraction.extend(manifest)
 root=url_root+'E512_original_member_SHA_CRC_critical_source_audit.json'
 receipt=json.loads(z2.read(root));assert receipt['archive']==source+'.zip' and receipt['central_directory_sha256']==cfg['central']
 info=receipt['complete_source_authentication']
 standard=solve(original(z2,info,'main','gurobi_model.sol'))
 alternate=solve(original(z2,info,'cost','gurobi_model.sol'))
 assert standard.keys()==alternate.keys()
 cur=original(z2,info,'cost','curation.moma').decode('utf8')
 data=z3.read(gl+'/E512b_exact_software_division_constraint_mechanism.json')
 ev=json.loads(data);assert ev['gl']==gl and ev['all_optimizer_variable_names_identical']
 bydirect={x['alternate_selected_mapping']:x['source_event_frame_0_based'] for x in ev['disputed_events'] if x['alternate_curation_has_direct_forced_mapping']}
 # exact raw curation ASC lines and assignments from archived original software solution, not E518 CSV
 entries=[];seen=set()
 for line in cur.splitlines():
  s=line.strip()
  if not s.startswith('ASC,'):continue
  fields=[v.strip() for v in s.split(',',4)]
  assert len(fields)==5,fields
  word,t,idnumber,rhs,var=fields
  assert float(rhs)==1 and var in standard and var in alternate and var not in seen,(gl,var)
  seen.add(var)
  t=int(t);assert t==int(re.search(r'^(?:Map|Div|Exit|Lys)T(\d+)_',var).group(1))+1
  prior=float(standard[var]);later=float(alternate[var])
  assert abs(later-1)<1e-9,(gl,var,later)
  violated=abs(prior-1)>1e-9
  direct=var in bydirect
  bbox=re.search(r'_HypT\d+T(-?\d+)B(-?\d+)',var)
  rec={'gl':gl,'ASC_time_printed':t,'ASC_source_frame_zero_based':t-1,'assignment_variable':var,'force_rhs':1.0,
       'standard_solution_value':prior,'alternate_solution_value':later,'standard_violates_cost_ASC':violated,
       'directly_same_mother_dispute':direct,'disputed_event_frame':bydirect.get(var,''),
       'mother_component_y_min':int(bbox.group(1)) if bbox else '', 'mother_component_y_max':int(bbox.group(2)) if bbox else '',
       'classification':'DIRECT_DISPUTE_FORCED_MAP' if direct else 'OTHER_ASC_VIOLATED' if violated else 'ASC_ALREADY_COMPATIBLE'}
  entries.append(rec)
  if violated:violations.append(rec)
  all_asc.append(rec)
 n_violate=sum(q['standard_violates_cost_ASC'] for q in entries)
 n_direct=sum(q['directly_same_mother_dispute'] for q in entries)
 assert n_violate==(8 if gl=='Pos5_GL35' else 17)
 assert n_direct==(7 if gl=='Pos5_GL35' else 3)
 assert all(q['standard_violates_cost_ASC'] for q in entries if q['directly_same_mother_dispute'])
 remain=[q for q in entries if q['standard_violates_cost_ASC'] and not q['directly_same_mother_dispute']]
 assert len(remain)==(1 if gl=='Pos5_GL35' else 14)
 groups[gl]={'saved_variables':len(standard),'source_cost_ASC':len(entries),
             'cost_ASC_violated_by_standard_solution':n_violate,
             'violated_direct_dispute_ASC':n_direct,'violated_non_direct_ASC':len(remain),
             'fixed_original_standard_solution_requires_at_least_this_many_ASC_changes':n_violate,
             'minimum_to_make_original_standard_vector_ASC_compatible_with_direct_dispute_ASC_deleted':len(remain),
             'non_direct_violated_source_timepoints':[q['ASC_source_frame_zero_based'] for q in remain],
             'main_solution_full_counterfactual_feasibility':'NOT_TESTED_NO_ORIGINAL_MPS_BYTES',
             'retained_original_source_MPS_LP_members':len(manifest)}
 assert len(entries)==(25 if gl=='Pos5_GL35' else 216)
 z1.close();z2.close();z3.close()
assert len(all_asc)==241 and len(violations)==25 and len(extraction)==8
assert sum(x['directly_same_mother_dispute'] for x in violations)==10
assert sum(not x['directly_same_mother_dispute'] for x in violations)==15
assert sum(x['ASC_source_frame_zero_based'] in [369,373,377,379,381,386,391,395,398,401,403,405] for x in violations if x['gl']=='Pos12_GL22')==12

# A 25-element violation set is an EXACT ASC-only incompatibility certificate:
# Fixing x to the archived standard solution, a distinct equality x_a=1 is violated
# at each of these 25 unique binary assignment variables. Changing fewer than 25
# of the original ASC equalities while fixing x cannot satisfy this ASC subset.
# It does NOT certify feasibility under the remaining thousands of MPS rows.
with (OUT/'E520_full_241_ASC_source_cross_solution_audit.csv').open('w',newline='',encoding='utf8') as fp:
 fields=list(all_asc[0]);w=csv.DictWriter(fp,fieldnames=fields);w.writeheader();w.writerows(sorted(all_asc,key=lambda z:(z['gl'],z['ASC_time_printed'],z['assignment_variable'])))
with (OUT/'E520_minimum_25_required_ASC_changes_for_fixed_standard_solution.csv').open('w',newline='',encoding='utf8') as fp:
 fields=list(violations[0]);w=csv.DictWriter(fp,fieldnames=fields);w.writeheader();w.writerows(sorted(violations,key=lambda z:(z['gl'],z['ASC_time_printed'],z['assignment_variable'])))
with (OUT/'E520_remaining_15_violations_after_removing_direct_ten.csv').open('w',newline='',encoding='utf8') as fp:
 x=[r for r in violations if not r['directly_same_mother_dispute']]
 fields=list(x[0]);w=csv.DictWriter(fp,fieldnames=fields);w.writeheader();w.writerows(sorted(x,key=lambda z:(z['gl'],z['ASC_time_printed'])))
(OUT/'E520_8_original_official_MPS_LP_exact_zip_paths.json').write_text(json.dumps({'manifest_source':'E511 official original ZIP64 fully verified central directories and GL member listings','global_manifest_SHA256':MAN_SHA,'files':extraction,'note':'MPS and LP existence proven in archived source filename listings; payload CRC/SHA not established until separately extracted'},indent=2,ensure_ascii=False),encoding='utf8')
cert={'stage':'901 E520','science_tier':'exact archived source-ASC compatibility certificate only',
      'full_original_MPS_extracted':False,'external_omnipose_weight_bytes_obtained':False,'independent_human_labels_obtained':False,
      'groups':groups,'total_ASCs_checked':len(all_asc),'total_source_cost_ASC_violations_for_standard_vector':len(violations),
      '10_disputed_direct_ASC_mandated_changes':10,'additional_ASC_mandated_changes_for_original_standard_vector':15,
      'minimum_ASC_equalities_to_relax_to_restore_compatibility_of_fixed_main_vector':25,
      'ASC_only_compatibility_after_25_relaxations':True,'full_model_feasibility_after_25_relaxations':'NOT_TESTED',
      'cross_version_full_model_feasible_region_same':'NOT_ESTABLISHED',
      '12_physical_division_ground_truth':'ALL_UNRESOLVED','physical_molecular_bits':'NOT_IDENTIFIED','publication':'PUBLICATION_NO_GO',
      'official_zip_integrity_receipts':receipts,'MPS_LP_original_target_entries':len(extraction),
      'proof':'For each of the 25 distinct original ASC equality rows x_a=1, archived main x_a=0, so any fixed-x ASC-compatible representation must relax/change each such row; exactly 25 adjustments suffice for the ASC-only subset, but no non-ASC constraints assessed.'}
(OUT/'E520_ASC_only_exact_counterfactual_necessity_certificate.json').write_text(json.dumps(cert,ensure_ascii=False,indent=2),encoding='utf8')
(OUT/'E520_limitations_and_precise_next_step.md').write_text('''# 901 E520 original-source findings and exact limitations\n\n- Original E511 authenticated ZIP64 GL listings show 8 genuine `gurobi_model.mps`/`.lp` entries across two GL × two versions. The MPS bytes have **not** been fetched or SHA-certified yet.\n- 241 archived cost-version ASC equalities are `x_assignment=1`; the fixed main-version archived solution violates exactly 25, including 10 direct disputed mother mappings. Removing only these 10 leaves **15 more** direct contradictions to the same fixed original main assignment vector.\n- Of the 15 residual violations, one is Pos5 f247 on a different mother candidate; 14 are Pos12 f261,274,369,373,377,379,381,386,391,395,398,401,403,405; 12 form a dense serial run in f369–405. These are not 15 distinct physical cells.\n- This is an exact combinatorial compatibility certificate for the **ASC subset only**. No known complete feasible solution to the cost-version MPS after ASC deletion is claimed. The true source model constraints, cross constraints, objective coefficients, variable types, and exact MPS matrix must be checked before any solver counterfactual.\n- Original MPS is strongly preferred to LP for coefficient fidelity. Source reconstructed here from original archivist inventory, not guessed filename paths. No Omnipose model output and no independent biologist review.\n- Twelve actual biological fates remain UNRESOLVED. No physical molecular information/bit/s/bit/J is derived. Publication NO_GO.\n''',encoding='utf8')
print('E520_EXACT_SOURCE_ASC_CERTIFICATE_PASS',json.dumps({**groups,'total_ASCs':len(all_asc),'main_violations':len(violations),'remaining_after_direct_ten':15,'original_MPS_LP_targets':len(extraction)},ensure_ascii=False),flush=True)
