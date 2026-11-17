#!/usr/bin/env python3
"""901 E527: frozen real-source component-vs-geometry cost audit of 151 one-frame MoMA swaps.

Scientific boundaries: original saved cost coefficients and component candidates are
source-software inputs, NOT independent image noise or biological division truth.
Public MoMA source (pinned df7a3dd92...) gives the compatibility-cost construction;
this script verifies it against archived 2022 cost CSV before interpretation.

Usage: python e527_original_cost_mechanism_exact_decomposition.py \
   --base /mnt/data --output /mnt/data/901_E527_20261117/results
"""
import argparse,csv,hashlib,io,json,math,pathlib,statistics,zipfile
from collections import Counter,defaultdict
from decimal import Decimal,getcontext
import numpy as np
getcontext().prec=46

P=argparse.ArgumentParser();P.add_argument('--base',default='/mnt/data');P.add_argument('--output',default='/mnt/data/901_E527_20261117/results');p=P.parse_args()
BASE=pathlib.Path(p.base);OUT=pathlib.Path(p.output);OUT.mkdir(parents=True,exist_ok=True)
ZS={"Pos5_GL35":"77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6",
    "Pos12_GL22":"73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2"}
RANKED_SHA="b04584801314d8d7a5b9b848ae8bec038a18aa8570e777348111af22464928bb"
SIDE_SHA={"Pos5_GL35":"1dad2aee7728757330fef6f16127faa5b5ab1648b78080db78a24621b2ffa596",
          "Pos12_GL22":"be8695b97cb1a2c946a8f667b201a40574e064e8b1b203124a40cb4bca5adcba"}
RANKED=BASE/'901_E526_20261117/results/E526_151_local_moves_exact_integer_ranked.csv'
D0=Decimal('0');D05=Decimal('0.5');D1=Decimal('1')
def sha(b):return hashlib.sha256(b).hexdigest()
def require_members(z,receipt,version,kind):
    item=[x for x in receipt['complete_source_authentication'] if x['version']==version and x['kind']==kind]
    assert len(item)==1,(version,kind)
    r=item[0];pref=z.namelist()[0].split('/')[0]+'/'+receipt['gl']+'/'
    key=pref+version+'__'+kind
    assert key in z.namelist(),key
    content=z.read(key);assert sha(content)==r['original_sha256'] and len(content)==r['uncompressed_bytes'],key
    return content,r

def calc_growth(source,target,touches):
    v=(target-source)/340.
    if v>0:
       v=max(0.,v-0.05);expo=4
    else:
       if touches:return 0.
       expo=40
    v=abs(v)
    return v*(1+v)**expo

def calc_migration(a,b):
    # Exact branch semantics of source CostFactory.getMigrationCost (downward
    # Math.max() return is discarded, so no subtract-0.01 effect).
    v=(a-b)/340.
    if v>0:
       v=max(0.,v-0.05);exp=3
    else:
       exp=6
    v=abs(v)
    return v*(1+v)**exp

def cost_parts(action,components):
    arr=action.split('_');t=arr[0];typ=t.split('T')[0]
    if typ not in ('Map','Div','Exit','Lys'):raise ValueError(typ)
    cc=[components[k] for k in arr[1:]]
    if typ=='Exit':return D0,0.,0.
    if typ=='Lys':return D0,0.,10.
    assert len(cc)==(2 if typ=='Map' else 3)
    ys=[(min(x['yCoordinates']),max(x['yCoordinates'])) for x in cc]
    sz=[b-a for a,b in ys]
    comp=D05*(cc[0]['cost']+sum((x['cost'] for x in cc[1:]),D0))
    touches=ys[1][0]<=120
    growth=calc_growth(sz[0],sum(sz[1:]),touches)
    migration=0.5*calc_migration(ys[0][0],ys[1][0])+0.5*calc_migration(ys[0][1],ys[-1][1])
    return comp,growth,migration

assert sha(RANKED.read_bytes())==RANKED_SHA
ranked=list(csv.DictReader(RANKED.open(encoding='utf8',newline='')))
assert len(ranked)==151 and Counter(x['GL'] for x in ranked)=={'Pos5_GL35':91,'Pos12_GL22':60}
source_receipts=[];source_details={};output=[];deltas={};by_component=[]
for gl in ZS:
    side=BASE/f'901_E511_{gl}_original_moma_sidecars.zip'
    assert sha(side.read_bytes())==SIDE_SHA[gl]
    with zipfile.ZipFile(side) as zz:
        assert zz.testzip() is None
        props=[]
        for which in ('main__2__mm.properties','cost_function__5__mm.properties'):
            names=[n for n in zz.namelist() if n.endswith(which)]
            assert len(names)==1
            rawprop=zz.read(names[0])
            table={k.strip():v.strip() for l in rawprop.decode(errors='replace').splitlines() if '=' in l and not l.lstrip().startswith('#') for k,v in [l.split('=',1)]}
            assert table['USE_FEATURE_MIGRATION_COSTS']=='1'
            assert int(table['CELL_DETECTION_ROI_OFFSET_TOP'])==120
            assert table['ASSIGNMENT_COST_CUTOFF']=='2.0'
            props.append({'filename':names[0],'sha256':sha(rawprop),'migration_flag':table['USE_FEATURE_MIGRATION_COSTS'],'roi_y_top':table['CELL_DETECTION_ROI_OFFSET_TOP']})
    source_details[gl]={'official_E511_mm_properties_zip_sha256':SIDE_SHA[gl],
                         'mm_properties_both_versions_immutable':props}
    file=BASE/f'901_E512_{gl}_full_original_components.zip'
    raw=file.read_bytes();assert sha(raw)==ZS[gl]
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
       assert z.testzip() is None
       indexname=next(n for n in z.namelist() if n.endswith('E512_original_member_SHA_CRC_critical_source_audit.json'))
       rec=json.loads(z.read(indexname))
       assert rec['gl']==gl
       component_bytes,rc=require_members(z,rec,'cost','component_forests.json')
       cost_bytes,ra=require_members(z,rec,'cost','assignment_costs.csv')
       sol_bytes,rs=require_members(z,rec,'cost','gurobi_model.sol')
       comps=json.loads(component_bytes,parse_float=Decimal)['components']
       components={c['stringId']:c for c in comps}
       assert len(components)==len(comps)
       costs={row['id']:Decimal(row['cost']) for row in csv.DictReader(cost_bytes.decode().splitlines()[1:])}
       sol={p[0]:float(p[1]) for line in sol_bytes.decode().splitlines() if line and not line.startswith('#') for p in [line.split()]}
       assert len(costs) == (44231 if gl=='Pos5_GL35' else 47879)
       source_receipts.append({'GL':gl,'archive':file.name,'SHA256':sha(raw),'original_cost_member_SHA256':ra['original_sha256'],
         'original_forest_member_SHA256':rc['original_sha256'],'original_solution_member_SHA256':rs['original_sha256'],
         'all_zip_members_CRC':True,'candidate_components':len(comps),'original_assignment_cost_entries':len(costs)})
       group=[row for row in ranked if row['GL']==gl]
       for case in group:
         frame=int(case['source_division_t']);direction=case['direction']; assert direction in ('early','late'),direction
         tokens=[x.strip() for x in case['five_flipped_binary_actions'].split(';')]
         assert len(tokens)==5,(gl,frame,tokens)
         five=[];total=D0;component=D0;residual=D0;growth=0.;migration=0.;max_abs_recon=0.;sum_abs_cost=D0
         for token in tokens:
            action,snew=token.rsplit('=',1);new=int(snew);assert new in (0,1)
            assert action in sol and action in costs and int(sol[action])==1-new,(gl,frame,action)
            sign=1 if new==1 else -1
            cp,gr,mg=cost_parts(action,components)
            original=costs[action]
            if action.startswith('LysT') or action.startswith('ExitT'):
              compat=original-cp
            else:
              compat=original-cp
            c_err=abs(float(compat)-(gr+mg))
            max_abs_recon=max(max_abs_recon,c_err)
            total+=sign*original;component+=sign*cp;residual+=sign*compat
            growth+=sign*gr;migration+=sign*mg
            sum_abs_cost+=abs(original)
            parts=action.split('_')
            recitem={'GL':gl,'frame':frame,'direction':direction,'action':action,'new_bit':new,'archive_original_cost':str(original),
              'source_component_cost_contribution':str(cp), 'observed_non_component_residual':str(compat),
              'reconstructed_geometry_growth':format(gr,'.14g'),'reconstructed_geometry_migration':format(mg,'.14g'),
              'component_cost_count':len(parts)-1,'geometry_residual_error_approx':format(c_err,'.14g')}
            by_component.append(recitem)
            five.append((action,new,sign,parts))
         assert total==Decimal(case['delta_exact_decimal']),(gl,frame,total,case['delta_exact_decimal'])
         assert total==component+residual
         t0=min(int(it[3][0].split('T')[1]) for it in five)
         assert t0==(frame-1 if direction=='early' else frame),(gl,frame,direction,t0)
         changed0=[it for it in five if int(it[3][0].split('T')[1])==t0]
         changed1=[it for it in five if int(it[3][0].split('T')[1])==t0+1]
         assert len(changed0)==2 and len(changed1)==3,(gl,frame,direction,[(x[0].split('_')[0],x[1]) for x in five])
         if direction=='early':
           mp=next(it for it in changed0 if it[0].startswith('Map') and it[1]==0)
           dv=next(it for it in changed0 if it[0].startswith('Div') and it[1]==1)
         else:
           mp=next(it for it in changed0 if it[0].startswith('Map') and it[1]==1)
           dv=next(it for it in changed0 if it[0].startswith('Div') and it[1]==0)
         assert mp[3][1]==dv[3][1],(gl,frame,'same_source')
         unsplit=mp[3][2]
         children=dv[3][2:]
         assert len(children)==2 and all(x in components for x in children)
         assert all(components[x]['frameNumber']==t0+1 for x in [unsplit]+children)
         sgn=1 if direction=='early' else -1
         split_difference=(components[children[0]]['cost']+components[children[1]]['cost']-components[unsplit]['cost'])*sgn
         parent_mask=set(zip(components[unsplit]['xCoordinates'],components[unsplit]['yCoordinates']))
         ch1_mask=set(zip(components[children[0]]['xCoordinates'],components[children[0]]['yCoordinates']))
         ch2_mask=set(zip(components[children[1]]['xCoordinates'],components[children[1]]['yCoordinates']))
         missing=parent_mask-(ch1_mask|ch2_mask)
         extra=(ch1_mask|ch2_mask)-parent_mask
         overlap=ch1_mask&ch2_mask
         if len(missing)==1 and not extra and not overlap:
             logit=float((components[children[0]]['cost']+components[children[1]]['cost']-components[unsplit]['cost'])/Decimal('0.2')*Decimal(2800))
             implied_p=1.0/(1.0+math.exp(-logit))
         else:implied_p=None
         assert split_difference==component,(gl,frame,'component telescoping failed',split_difference,component)
         for cx in changed1:
           if direction=='early':
             if cx[0].startswith('Div'):assert cx[1]==0 and cx[3][1]==unsplit
             else:assert cx[0].startswith('Map') and cx[1]==1 and cx[3][1] in children
           else:
             if cx[0].startswith('Div'):assert cx[1]==1 and cx[3][1]==unsplit
             else:assert cx[0].startswith('Map') and cx[1]==0 and cx[3][1] in children
         is_geom_zero=all(abs(float(recitem['reconstructed_geometry_growth']))<1e-15 and abs(float(recitem['reconstructed_geometry_migration']))<1e-15 for recitem in by_component[-5:])
         out={'GL':gl,'source_division_t':frame,'direction':direction,'delta_exact_source_cost':str(total),
           'component_telescoping_delta':str(component),'noncomponent_exact_archive_residual':str(residual),
           'reconstructed_growth_delta_approx':format(growth,'.14g'), 'reconstructed_migration_delta_approx':format(migration,'.14g'),
           'full_geometry_predicted_delta_approx':format(growth+migration,'.14g'),
           'max_abs_error_source_cost_vs_reconstructed_parts_per_action':format(max_abs_recon,'.14g'),
           'all_five_geometry_costs_zero':is_geom_zero,'component_vs_unsplit_same_frame_t_plus1':t0+1,
           'unsplit_component_id':unsplit,'split_upper_component_id':children[0],
           'split_lower_component_id':children[1],
           'unsplit_pixel_count':len(parent_mask),
           'daughter_union_pixel_count':len(ch1_mask|ch2_mask),
           'daughter_mask_overlap_px':len(overlap),
           'parent_missing_in_daughters_px':len(missing),
           'daughter_pixels_outside_parent_px':len(extra),
           'missing_pixel_coordinates_0based':str(sorted(missing)) if len(missing)<=8 else 'MANY_NOT_LISTED',
           'implied_single_gap_probability_CONDITIONAL_loglikelihood_method':format(implied_p,'.12g') if implied_p is not None else '',
           'unverified_historical_component_cost_mode':True,
           'unsplit_cost_archived':str(components[unsplit]['cost']),
           'split_costs_sum_archived':str(components[children[0]]['cost']+components[children[1]]['cost']),
           'abs_cost_sum_of_five_flips':str(sum_abs_cost),
           'ratio_margin_over_sum_abs_cost':format(float(total/sum_abs_cost),'.14g'),
           'NOT_BIOLOGICAL_TRUTH':True}
         output.append(out)
    print(gl,'source_sha_verified',len(group))

assert len(output)==151 and len(by_component)==755
by_gl={}
for gl in ZS:
    q=[x for x in output if x['GL']==gl]
    q=sorted(q,key=lambda v:Decimal(v['delta_exact_source_cost']))
    errs=[float(v['max_abs_error_source_cost_vs_reconstructed_parts_per_action']) for v in q]
    zeros=[x for x in q if x['all_five_geometry_costs_zero']]
    zero_comp_only=[x for x in q if abs(Decimal(x['noncomponent_exact_archive_residual']))<=Decimal('0.0000002')]
    near=[x for x in q if Decimal(x['delta_exact_source_cost'])<=Decimal('0.001')]
    by_gl[gl]={
      'candidate_count':len(q),'near_delta_le_0_001_count':len(near),
      'fully_geometry_penalty_free_actions_count':len(zeros),
      'residual_delta_at_most_2e_7_count':len(zero_comp_only),
      'max_action_model_cost_reconstruction_residual':max(errs),
      'median_action_model_cost_reconstruction_residual':statistics.median(errs),
      'lowest_margin':q[0],
      'second_margin':q[1],
      'geom_penalty_free_near_count':sum(x['all_five_geometry_costs_zero'] for x in near),
      'near_component_only_residual_le_2e_7_count':sum(abs(Decimal(x['noncomponent_exact_archive_residual']))<=Decimal('0.0000002') for x in near),
      'all_151_one_pixel_gap_count':sum(int(x['parent_missing_in_daughters_px'])==1 and int(x['daughter_pixels_outside_parent_px'])==0 and int(x['daughter_mask_overlap_px'])==0 for x in q),
      'near_one_pixel_gap_count':sum(int(x['parent_missing_in_daughters_px'])==1 and int(x['daughter_pixels_outside_parent_px'])==0 and int(x['daughter_mask_overlap_px'])==0 for x in near),
      'missing_pixels_distribution':dict(Counter(int(x['parent_missing_in_daughters_px']) for x in q)),
      'extra_pixels_distribution':dict(Counter(int(x['daughter_pixels_outside_parent_px']) for x in q))
    }
rows=sorted(output,key=lambda v:(v['GL'],int(v['source_division_t']),v['direction']))
with (OUT/'E527_151_exact_source_cost_component_geometry_decomposition.csv').open('w',newline='',encoding='utf8') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
with (OUT/'E527_755_original_actions_cost_layer_reconstruction.csv').open('w',newline='',encoding='utf8') as f:
 w=csv.DictWriter(f,fieldnames=list(by_component[0]));w.writeheader();w.writerows(by_component)
report={'E527_SCOPED':'SOFTWARE_COMPONENT_COST_AND_GEOMETRIC_COMPATIBILITY_NOT_BIOLOGICAL_TRUTH',
  'source_artifacts':source_receipts,'original_mm_properties':source_details,'ranked_151_sha256':sha(RANKED.read_bytes()),
  'pinned_public_source_commit':'df7a3dd92cd15344c82b19e83f96a1f132e52ab7',
  'component_cost_loglikelihood_mode_in_pinned_public_FeatureFlags_not_2022_runtime_proven':True,
  'source_formula':'Map=.5*Csource+.5*Ctarget+growth+mean_boundary_migration; Div=.5*Csource+.5*(Cup+Clow)+growth(source,sum daughter lengths)+mean_boundary_migration; Exit=0; Lys=config cost',
  'cost_factory_normalizer_px_340':340,'growth_free_up_to_17_px':True,'upward_migration_free_up_to_17_px':True,
  'source_exported_component_cost_is_an_algorithmic_value_not_probability_or_physical_energy':True,
  'all_151_local_five_flips_reproduced_exact_from_original_cost_csv':True,
  'source_component_cost_telescoping_verified_all_151':True,
  'summary_by_GL':by_gl,
  'true_float32_image_component_uncertainty_calibration':'NOT_IDENTIFIED',
  'full_original_MPS':'NOT_ACQUIRED',
  'third_party_segmentation_ground_truth':'NOT_ACQUIRED',
  'physical_division_disputes_unresolved':12,'publication':'NO_GO'}
(OUT/'E527_source_mechanism_and_residuals.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,default=str),encoding='utf8')
print('E527_SOURCE_COMPONENT_DECOMPOSITION_VERIFIED',json.dumps({gl:{k:v for k,v in d.items() if k not in ('lowest_margin','second_margin')} for gl,d in by_gl.items()},ensure_ascii=False),flush=True)
print('E527_TWO_BEST',json.dumps({gl:{'first':by_gl[gl]['lowest_margin'],'second':by_gl[gl]['second_margin']} for gl in by_gl},ensure_ascii=False,default=str),flush=True)
