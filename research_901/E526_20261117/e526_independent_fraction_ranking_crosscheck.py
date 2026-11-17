#!/usr/bin/env python3
"""Independent exact Fraction verification of E526's rational ranking and centered-MILP witnesses.

Only local 5-flip family and reconstructed source-code MILP; no full MPS or biological truth.
"""
from __future__ import annotations
import csv,io,json,zipfile,hashlib,pathlib,collections
from fractions import Fraction
ROOT=pathlib.Path('/mnt/data')
OUT=ROOT/'901_E526_20261117'/'results';OUT.mkdir(exist_ok=True,parents=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
expected={
'Pos5_GL35':('20221220_rplN_GFP_pl','77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6',91),
'Pos12_GL22':('20221220_E08_GFP_pl','73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2',60)}
local=list(csv.DictReader((OUT/'E526_151_local_moves_exact_integer_ranked.csv').open()))
assert len(local)==151
summary={}
for gl,(src,digest,n) in expected.items():
 path=ROOT/f'901_E512_{gl}_full_original_components.zip';assert sha(path)==digest
 with zipfile.ZipFile(path) as z:
  assert z.testzip() is None
  prefix=f'{src}/{gl}/';receipts=json.loads(z.read(prefix+'E512_original_member_SHA_CRC_critical_source_audit.json'))['complete_source_authentication']
  def read(member):
   b=z.read(prefix+'cost__'+member)
   rc=[x for x in receipts if x['version']=='cost' and x['kind']==member]
   assert len(rc)==1 and hashlib.sha256(b).hexdigest()==rc[0]['original_sha256'] and len(b)==rc[0]['uncompressed_bytes']
   return b.decode()
  prices={r['id']:Fraction(r['cost']) for r in csv.DictReader(ln for ln in read('assignment_costs.csv').splitlines() if ln.strip() and not ln.startswith('#'))}
  sol={a.split()[0]:int(a.split()[1]) for a in read('gurobi_model.sol').splitlines() if a.strip() and not a.startswith('#')}
  assert len(prices)==len([k for k in sol if k.startswith(('MapT','DivT','ExitT','LysT'))])
  local_sub=[q for q in local if q['GL']==gl];assert len(local_sub)==n
  ratios=[]
  for row in local_sub:
   flips={x.rsplit('=',1)[0]:int(x.rsplit('=',1)[1]) for x in row['five_flipped_binary_actions'].split(';')}
   assert len(flips)==5 and all(sol[k]==1-v for k,v in flips.items())
   d=sum((v-sol[k])*prices[k] for k,v in flips.items())
   assert d>0 and d==Fraction(row['delta_exact_decimal'])
   ratios.append((d,int(row['source_division_t']),row['direction']))
  ratios.sort()
  summary[gl]={
    'n_local_five_bit_feasible':n,
    'top_three_source_exact_fractions':[{'t':t,'direction':direction,'numerator':d.numerator,'denominator':d.denominator,'as_decimal':str(float(d))} for d,t,direction in ratios[:3]],
    'all_candidates_recomputed_with_Python_Fraction':True,
    'original_source_member_sha_crc_pass':True,
  }
  for factor in ('scale1','scale1e5'):
   p=OUT/f'centered_{gl.replace("_GL","_").lower()}_{factor}'
   # actual path intentionally selected explicitly below, not guessed
  tagged='pos5' if gl=='Pos5_GL35' else 'pos12'
  centered=[]
  for scale in ['scale1','scale1e5']:
   p=OUT/f'centered_{tagged}_{scale}'/f'E526_{gl}_centered_cost_best_distinct_changed_actions.csv'
   if not p.exists():continue
   rr=list(csv.DictReader(p.open()))
   assert rr and len(rr)==5
   moved={x['action']:int(x['new']) for x in rr}
   d=sum((v-sol[k])*prices[k] for k,v in moved.items())
   assert d==ratios[0][0]
   match=[x for x in local_sub if {a for a in x['five_flipped_binary_actions'].split(';')}=={k+'='+str(v) for k,v in moved.items()}]
   assert len(match)==1 and int(match[0]['source_division_t'])==ratios[0][1]
   centered.append({'scale':scale,'exact_delta_as_decimal':match[0]['delta_exact_decimal'],'match_local_optimum':True,'moved_bits':len(rr)})
  summary[gl]['independent_centered_optimizer_witnesses']=centered
(OUT/'E526_independent_fraction_ranking_crosscheck.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,sort_keys=True)+'\n')
print('E526_INDEPENDENT_FRACTION_PASS',json.dumps(summary,ensure_ascii=False,sort_keys=True))
