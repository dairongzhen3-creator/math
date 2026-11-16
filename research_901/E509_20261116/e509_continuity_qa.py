#!/usr/bin/env python3
"""901 E509: independent event-level source-geometry and one-frame continuity audit.

Only authenticated 2022-12-20 BioImage Archive S-BIAD2506 CellStats exported by
frozen E507/E508 source scripts. No image-pixel ground truth and no inference of
phototoxicity causality or molecular information transfer.
"""
import argparse, collections, csv, hashlib, json, pathlib, statistics, sys

A=argparse.ArgumentParser()
A.add_argument('--archive',required=True,choices=[
 '20221220_E08_GFP.zip','20221220_E08_GFP_pl.zip','20221220_rplN_GFP_pl.zip'])
a=A.parse_args()
ROOT=pathlib.Path('e507_alternate_QA')/a.archive[:-4]
CON=pathlib.Path('e508_source_conflicts')/a.archive[:-4]
OUT=pathlib.Path('e509_continuity_qa')/a.archive[:-4];OUT.mkdir(parents=True,exist_ok=True)
FILES=['E507_alternate_algorithm_frame_geometry_results.json']
summ=json.loads((ROOT/FILES[0]).read_text())
conflicts=json.loads((CON/'E508_event_level_version_conflicts.json').read_text())
by_id={(x['GL'],x['frame'],x['parent']):x for x in conflicts}
assert len(by_id)==len(conflicts)
out=[]
source_receipts=[]
for group in summ:
 gl=group['GL']
 def load(v):
  p=ROOT/(gl+'_'+v+'_CellStats.csv')
  lines=[x for x in p.read_text(encoding='utf-8-sig').splitlines() if x.strip() and not x.startswith('#')]
  source_receipts.append({'GL':gl,'version':v,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
  return list(csv.DictReader(lines))
 mr,cr=load('main'),load('cost')
 all_main=collections.defaultdict(list)
 for r in mr:all_main[(int(r['cell_ID']),int(r['frame']))].append(r)
 main_children=collections.defaultdict(list)
 for r in mr:main_children[(int(r['parent_ID']),int(r['frame']))].append(r)
 cost_at={(int(r['cell_ID']),int(r['frame'])):r for r in cr}
 cost_frame=collections.defaultdict(list)
 for r in cr:cost_frame[int(r['frame'])].append(r)
 matched={(int(e['frame']),int(e['main_parent'])) for e in group['pairs']}
 valid_divs=[]
 for r in mr:
  if r['type_of_end']!='div':continue
  t=int(r['frame']);pid=int(r['cell_ID'])
  if not 200<=t<=700:continue
  girls=main_children[(pid,t+1)]
  if len({int(x['cell_ID']) for x in girls})!=2:continue
  if len(girls)!=2:raise AssertionError(('duplicate daughter rows',gl,t,pid))
  large=max(girls,key=lambda x:float(x['length_px']))
  small=min(girls,key=lambda x:float(x['length_px']))
  mother_len=float(r['length_px'])
  m=(t,pid) in matched
  rec={
   'archive':a.archive,'GL':gl,'parent_ID_main':pid,'division_frame_main':t,'cross_algorithm_matched':m,
   'mother_length_px':mother_len,
   'large_daughter_ID_main':int(large['cell_ID']),
   'large_daughter_length_px':float(large['length_px']),
   'large_daughter_to_mother_length_ratio':float(large['length_px'])/mother_len,
   'small_daughter_ID_main':int(small['cell_ID']),
   'small_daughter_length_px':float(small['length_px']),
   'small_daughter_bbox_top_px':int(small['bbox_top_px']),
   'small_daughter_bbox_bottom_px':int(small['bbox_bottom_px']),
   'small_daughter_bbox_top_equals_120':int(small['bbox_top_px'])==120,
  }
  if not m:
   key=(gl,t,pid)
   assert key in by_id,(gl,t,pid)
   nearby=by_id[key]['nearest_three_cost_cells']
   assert nearby,(gl,t,pid)
   nearest=nearby[0]
   cid=int(nearest['cost_cell_id'])
   r0=cost_at.get((cid,t));r1=cost_at.get((cid,t+1))
   assert r0 is not None
   rec.update({
    'nearest_cost_cell_ID':cid,
    'cost_mother_at_t_end_label':r0['type_of_end'],
    'cost_same_cell_ID_observed_at_t_plus_1':r1 is not None,
    'cost_cell_at_t_plus_1_end_label':r1['type_of_end'] if r1 is not None else None,
    'cost_at_t_plus_1_length_px':float(r1['length_px']) if r1 is not None else None,
    'cost_at_t_plus_1_bbox_top_px':int(r1['bbox_top_px']) if r1 is not None else None,
    'cost_at_t_plus_1_bbox_bottom_px':int(r1['bbox_bottom_px']) if r1 is not None else None,
    'cost_same_cell_ID_has_div_label_at_t_plus_1':bool(r1 is not None and r1['type_of_end']=='div'),
    'cost_mother_at_t_length_px':float(r0['length_px']),
    'cost_mother_bbox_identical_to_main_at_t':(
        int(r0['bbox_top_px'])==int(r['bbox_top_px']) and
        int(r0['bbox_bottom_px'])==int(r['bbox_bottom_px'])),
    'cost_mother_length_identical_to_main_at_t':float(r0['length_px'])==mother_len,
    'e508_original_event_id':by_id[key]['event_id']
   })
  out.append(rec);valid_divs.append(rec)
 assert len(valid_divs)==group['main_divisions_frames200to700'],(gl,len(valid_divs),group)
 assert sum(x['cross_algorithm_matched'] for x in valid_divs)==group['matched_geometry_divisions']
 assert sum(not x['cross_algorithm_matched'] for x in valid_divs)==group['unmatched_main']
 assert len({(x['division_frame_main'],x['parent_ID_main']) for x in valid_divs})==len(valid_divs)
mainonly=[x for x in out if not x['cross_algorithm_matched']]
matched=[x for x in out if x['cross_algorithm_matched']]
assert len(mainonly)==len(conflicts)
def describe(xs):
 if not xs:return {'n':0}
 ratios=[x['large_daughter_to_mother_length_ratio'] for x in xs]
 return {'n':len(xs),'dominant_daughter_ratio_min':min(ratios),'median':statistics.median(ratios),
         'max':max(ratios),'short_daughter_top120_n':sum(x['small_daughter_bbox_top_equals_120'] for x in xs)}
report={
 'archive':a.archive,'only_source_level_evidence':True,'not_pixel_ground_truth':True,
 'matched_standard_divisions':describe(matched),'standard_only_divisions':describe(mainonly),
 'all_standard_divisions_200to700':len(out),
 'maximum_matched_dominant_daughter_ratio':max([x['large_daughter_to_mother_length_ratio'] for x in matched],default=None),
 'minimum_disputed_dominant_daughter_ratio':min([x['large_daughter_to_mother_length_ratio'] for x in mainonly],default=None),
 'disputed_alternate_same_ID_continues_tplus1':sum(x['cost_same_cell_ID_observed_at_t_plus_1'] for x in mainonly),
 'disputed_alternate_div_at_tplus1':sum(x['cost_same_cell_ID_has_div_label_at_t_plus_1'] for x in mainonly),
 'notes':'Source main version labels two daughters even where the larger daughter retains near-full mother length; pixel validation required. No biological p-values.',
 'sources':source_receipts
}
(OUT/'E509_division_geometry_and_continuity.json').write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')
(OUT/'E509_source_summary.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print('E509_SOURCE_PER_EVENT',json.dumps([{
 'gl':x['GL'],'frame':x['division_frame_main'],'main_parent':x['parent_ID_main'],
 'daughter_ratio':round(x['large_daughter_to_mother_length_ratio'],5),
 'small_daughter_top':x['small_daughter_bbox_top_px'],
 'cost_same_ID_at_next_frame':x['cost_same_cell_ID_observed_at_t_plus_1'],
 'cost_next_end_label':x['cost_cell_at_t_plus_1_end_label']} for x in mainonly],ensure_ascii=False),flush=True)
print('E509_SOURCE_SUMMARY_PASS',json.dumps(report,sort_keys=True,ensure_ascii=False),flush=True)
