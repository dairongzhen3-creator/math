#!/usr/bin/env python3
"""901 E512b: source-authenticated exact original MoMA candidate and ILP decision audit.

Uses E512 official original ZIP64-CRC-SHA-certified Google/GitHub cloud
artifacts, fixed from GitHub Actions run 43827391164. Does not re-download or
change any original source byte. An active ILP division variable is a software
decision, not biological evidence that a true division occurred.
"""
import argparse,collections,csv,hashlib,json,pathlib,re

p=argparse.ArgumentParser()
p.add_argument('--gl',required=True,choices=['Pos5_GL35','Pos12_GL22'])
p.add_argument('--source-root',default='verified_e512')
a=p.parse_args()
config={
 'Pos5_GL35':{'archive':'20221220_rplN_GFP_pl.zip',
  'disputed':[247,294,399,403,405,409,411,415,622],
  'n_components':5500,'n_variables':49731},
 'Pos12_GL22':{'archive':'20221220_E08_GFP_pl.zip',
  'disputed':[241,246,608],
  'n_components':5765,'n_variables':53644}
}
x=config[a.gl];root=pathlib.Path(a.source_root)
folders=list(root.rglob(a.gl))
assert len(folders)==1,(root,a.gl,folders)
source=folders[0]
out=pathlib.Path('e512_optimization_audit')/a.gl
out.mkdir(parents=True,exist_ok=True)
receipt=json.loads((source/'E512_original_member_SHA_CRC_critical_source_audit.json').read_text(encoding='utf8'))
assert receipt['gl']==a.gl and receipt['archive']==x['archive']
assert receipt['map_properties']['main']['source_file_byte_sha256']==receipt['map_properties']['cost']['source_file_byte_sha256']
assert receipt['map_properties']['main']['selected_array_sha256']==receipt['map_properties']['cost']['selected_array_sha256']
authenticated=receipt['complete_source_authentication']
for r in authenticated:
 if r['kind']=='probability_maps.tif':continue  # Full TIFF was authenticated before ephemeral deletion.
 candidate=source/(r['version']+'__'+r['kind'])
 assert candidate.is_file(),candidate
 assert hashlib.sha256(candidate.read_bytes()).hexdigest()==r['original_sha256'],candidate
print('E512B_CLOUD_ARTIFACT_SOURCE_SHA_PASS',a.gl,len(authenticated),flush=True)
def read(kind,version):
 return (source/(version+'__'+kind)).read_bytes()
m=json.loads(read('component_forests.json','main'))['components']
c=json.loads(read('component_forests.json','cost'))['components']
assert len(m)==len(c)==x['n_components']
def idmap(arr):
 d={}
 for row in arr:
  k=row['stringId']
  assert k not in d,k
  d[k]=row
 return d
md,cd=idmap(m),idmap(c)
assert set(md)==set(cd)
compare_fields=collections.Counter()
for k in md:
 for field in set(md[k])|set(cd[k]):
  if md[k].get(field)!=cd[k].get(field):compare_fields[field]+=1
structural=['stringId','frameNumber','label','parentStringId','childrenStringIds','value','xCoordinates','yCoordinates']
assert all(compare_fields[f]==0 for f in structural),compare_fields
null_to_present={f:sum(md[k].get(f) is None and cd[k].get(f) is not None for k in md) for f in compare_fields}
print('E512B_ALL_CANDIDATE_GEOMETRY_IDENTICAL',a.gl,len(md),json.dumps(dict(compare_fields),sort_keys=True),flush=True)

def solutions(v):
 out={}
 for line in read('gurobi_model.sol',v).decode('utf8','replace').splitlines():
  if not line or line.startswith('#'):continue
  pair=line.strip().rsplit(None,1)
  assert len(pair)==2,pair
  name,value=pair
  assert name not in out
  out[name]=float(value)
 return out
main,cost=solutions('main'),solutions('cost')
assert len(main)==len(cost)==x['n_variables'] and main.keys()==cost.keys()
changed=[k for k in main if main[k]!=cost[k]]
def typ(name):return re.split(r'T\d+_',name,maxsplit=1)[0]
changed_types=dict(collections.Counter(typ(k) for k in changed))
# Compare actual source-pipeline assignments and encoded selection constraints,
# but do not infer who wrote the .moma curation file.
def curation(v):
 raw=read('curation.moma',v).decode('utf8','replace')
 groups=collections.Counter()
 constraints=collections.defaultdict(list)
 for line in raw.splitlines():
  stripped=line.strip()
  if not stripped or stripped.startswith('#'):continue
  parts=[s.strip() for s in stripped.split(',')]
  groups[parts[0]]+=1
  if parts[0]=='ASC':
   assert len(parts)>=5
   vname=parts[4]
   constraints[vname].append({'time':int(parts[1]),'record':stripped})
 return {'counts':dict(groups),'constraints':constraints,'sha256':hashlib.sha256(raw.encode('utf8')).hexdigest()}
mc,cc=curation('main'),curation('cost')
event_details=[]
for f in x['disputed']:
 def selected(mode,kind):
  pref=kind+'T'+str(f)+'_'
  return sorted(k for k,v in mode.items() if k.startswith(pref) and v>0.5)
 mdiv=selected(main,'Div');cdiv=selected(cost,'Div')
 added_map=[k for k in selected(cost,'Map') if main[k]<0.5]
 lost_div=[k for k in mdiv if cost[k]<0.5]
 assert len(lost_div)==1 and len(added_map)>=1,(a.gl,f,lost_div,added_map)
 # Match the exact same source mother candidate, not other simultaneous cells.
 mother_component=lost_div[0].split('_')[1]
 maps_for_this_mother=[k for k in added_map if k.split('_')[1]==mother_component]
 assert len(maps_for_this_mother)==1,(a.gl,f,mother_component,maps_for_this_mother)
 mapvar=maps_for_this_mother[0]
 # The source daughter components are encoded directly in the selected Div name.
 daughters=lost_div[0].split('_')[2:]
 assert len(daughters)==2 and all(d in md for d in daughters),(a.gl,f,daughters)
 desc=[{'id':d,'frameNumber':md[d]['frameNumber'],
        'minY':min(md[d]['yCoordinates']),'maxY':max(md[d]['yCoordinates']),
        'numPixels':len(md[d]['yCoordinates']),
        'minX':min(md[d]['xCoordinates']),'maxX':max(md[d]['xCoordinates'])} for d in daughters]
 direct_asc=(mapvar in cc['constraints'])
 if direct_asc:
  assert not mapvar in mc['constraints']
  assert all(c['time']==f+1 for c in cc['constraints'][mapvar]),(f,mapvar,cc['constraints'][mapvar])
 rec={'gl':a.gl,'source_event_frame_0_based':f,
      'mother_component_id':mother_component,
      'main_selected_div_variable':lost_div[0],
      'main_div_variable_value':main[lost_div[0]],'cost_div_variable_value':cost[lost_div[0]],
      'alternate_selected_mapping':mapvar,
      'main_mapping_value':main[mapvar],'cost_mapping_value':cost[mapvar],
      'alternate_curation_has_direct_forced_mapping':direct_asc,
      'alternate_curation_direct_map_records':cc['constraints'].get(mapvar,[]),
      'main_curation_direct_map_records':mc['constraints'].get(mapvar,[]),
      'source_daughter_components_in_main_div':desc,
      'any_source_daughter_minY_120':any(d['minY']==120 for d in desc),
      'evidence_tier':'same biological record, original software source candidates and optimization results; not independent physical division ground truth'}
 event_details.append(rec)
 print('E512B_DISPUTED_DIV_TO_MAP',json.dumps({'gl':a.gl,'frame':f,'ASC':direct_asc,'daughter_minY':[q['minY'] for q in desc]},ensure_ascii=False),flush=True)
assert len(event_details)==len(x['disputed'])
out_json={'gl':a.gl,'source_archive':x['archive'],
 'probability_map_SHA256_main_equals_cost':True,
 'original_probability_sha256':receipt['map_properties']['main']['source_file_byte_sha256'],
 'component_count_each_version':len(md),
 'original_component_IDS_exactly_equal':True,
 'component_field_changes':dict(compare_fields),
 'component_geometry_parent_child_exactly_equal':True,
 'original_variable_set_count_each_version':len(main),
 'all_optimizer_variable_names_identical':True,
 'variables_different_selected_values':len(changed),
 'changed_variable_types':changed_types,
 'original_curation_counts_main':mc['counts'],
 'original_curation_counts_cost':cc['counts'],
 'curation_constraint_origin_human_vs_auto':'NOT_IDENTIFIED',
 'disputed_events':event_details,
 'direct_forced_alternate_mapping_count':sum(x['alternate_curation_has_direct_forced_mapping'] for x in event_details),
 'non_matched_alternate_mapping_framenumbers':[x['source_event_frame_0_based'] for x in event_details if not x['alternate_curation_has_direct_forced_mapping']],
 'source_receipt_sha':hashlib.sha256((source/'E512_original_member_SHA_CRC_critical_source_audit.json').read_bytes()).hexdigest(),
 'science_limit':'No ablation of optimizer/curation and no independent physical daughter pixel identification; do not equate source constraints with cause or biology',
 'publication_state':'PUBLICATION_NO_GO'}
(out/'E512b_exact_software_division_constraint_mechanism.json').write_text(json.dumps(out_json,ensure_ascii=False,indent=2),encoding='utf8')
print('E512B_SOURCE_DECISIONS_PASS',json.dumps({'gl':a.gl,'changed_var':len(changed),'events':len(event_details),'direct_ASC':out_json['direct_forced_alternate_mapping_count'],'main_ASC':mc['counts'].get('ASC',0),'cost_ASC':cc['counts'].get('ASC',0)},ensure_ascii=False),flush=True)
