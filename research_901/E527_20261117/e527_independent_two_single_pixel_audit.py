#!/usr/bin/env python3
"""Independent of the 151-family script: authenticate and rederive two 1-pixel split margins.

Uses 2022 official extracted archive bytes, not source TIFF, not model predictions;
conditional log-likelihood p-inference is NOT a measured source pixel probability.
"""
import csv,hashlib,json,math,pathlib,zipfile
from decimal import Decimal,getcontext
getcontext().prec=45
BASE=pathlib.Path('/mnt/data')
OUT=BASE/'901_E527_20261117/results'
OUT.mkdir(parents=True,exist_ok=True)
archsha={
'Pos5_GL35':'77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6',
'Pos12_GL22':'73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2',
}
cases={
'Pos5_GL35':{'f':960,'expected_delta':Decimal('0.00000519'),'px':(55,160),'source_output':'centered_pos5_scale1'},
'Pos12_GL22':{'f':118,'expected_delta':Decimal('0.00001565'),'px':(50,197),'source_output':'centered_pos12_scale1'},
}
res=[]
for gl,case in cases.items():
 zpath=BASE/('901_E512_'+gl+'_full_original_components.zip');raw=zpath.read_bytes();assert hashlib.sha256(raw).hexdigest()==archsha[gl]
 with zipfile.ZipFile(zpath) as z:
  assert z.testzip() is None
  n=next(n for n in z.namelist() if n.endswith('cost__component_forests.json'))
  compdata=z.read(n);comp={x['stringId']:x for x in json.loads(compdata,parse_float=Decimal)['components']}
  cpath=next(n for n in z.namelist() if n.endswith('cost__assignment_costs.csv'))
  price={r['id']:Decimal(r['cost']) for r in csv.DictReader(z.read(cpath).decode().splitlines()[1:])}
  f=case['f']
  changed=BASE/('901_E526_20261117/results/'+case['source_output']+'/E526_'+gl+'_centered_cost_best_distinct_changed_actions.csv')
  actions=list(csv.DictReader(changed.open()))
  assert len(actions)==5
  assert all(int(r['alt'])==1-int(r['new']) and price[r['action']]==Decimal(r['original_cost_decimal']) for r in actions)
  delta=sum(((Decimal(r['new'])-Decimal(r['alt']))*price[r['action']] for r in actions),Decimal('0'))
  assert delta==case['expected_delta']
  early_map=next(r for r in actions if r['action'].startswith('MapT'+str(f-1)+'_'))
  early_div=next(r for r in actions if r['action'].startswith('DivT'+str(f-1)+'_'))
  p0=early_map['action'].split('_')
  p1=early_div['action'].split('_')
  assert p0[1]==p1[1]
  whole=p0[2];twins=p1[2:4]
  original=[set(zip(comp[key]['xCoordinates'],comp[key]['yCoordinates'])) for key in [whole]+twins]
  assert not(original[1]&original[2])
  union=original[1]|original[2]
  assert union<=original[0]
  gap=sorted(original[0]-union)
  assert gap==[case['px']]
  sepCost=sum((comp[t]['cost'] for t in twins),Decimal('0'))-comp[whole]['cost']
  residual=delta-sepCost
  assert abs(residual)<=Decimal('0.00000004')
  k=Decimal('0.2')/Decimal(2800)
  logit=sepCost/k
  p_conditional=1/(1+math.exp(-float(logit)))
  # Forward expression: if costs were directly computed as pinned UsingLogLikelihoodCost
  # by additive pixel logit, the missing singleton's p would take this value.
  res.append({'GL':gl,'software_original_division_f':f,'source_whole_component':whole,'source_daughter_components':twins,
    'unsplit_component_area_px':len(original[0]),'twins_area_px':[len(original[1]),len(original[2])],
    'twins_overlap_px':len(original[1]&original[2]),'missing_gap_pixel_0based_xy':gap,
    'whole_component_cost_archived':str(comp[whole]['cost']),
    'twins_sum_component_cost_archived':str(comp[twins[0]]['cost']+comp[twins[1]]['cost']),
    'difference_component_source_decimal':str(sepCost),
    'five_action_total_difference_source_decimal':str(delta),
    'archive_cost_vs_saved_component_rounding_compatible_residual_decimal':str(residual),
    'conditional_p_under_pinned_public_UsingLogLikelihoodCost':format(p_conditional,'.12g'),
    'conditional_p_NOT_source_measured':True,
    'source_zip_sha256':archsha[gl],
    'comp_source_member_sha256':hashlib.sha256(compdata).hexdigest(),
    'model_and_semantics':'costs+geometry software audit; no original probability pixel; historical source feature flag not proven'} )
(OUT/'E527_independent_two_single_gap_exact_certificate.json').write_text(json.dumps({'cases':res,'scope':'source cost and candidate mask only; conditional probability not measured'},ensure_ascii=False,indent=2),encoding='utf8')
print('E527_SINGLE_PIXEL_SOURCE_EXACT_PASS',json.dumps([{k:v for k,v in x.items() if k in ('GL','software_original_division_f','missing_gap_pixel_0based_xy','difference_component_source_decimal','five_action_total_difference_source_decimal','archive_cost_vs_saved_component_rounding_compatible_residual_decimal','conditional_p_under_pinned_public_UsingLogLikelihoodCost')} for x in res],ensure_ascii=False),flush=True)
