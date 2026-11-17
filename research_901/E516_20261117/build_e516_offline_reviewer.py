#!/usr/bin/env python3
"""Build E516 offline, self-contained, source-faithful 37-case independent reader.
No external dependencies, tracking requests, analytics, MoMA masks, sealed labels,
source frame numbers, or attached source solver metadata in this package.
"""
import zipfile, hashlib, json, pathlib, html, csv, io

ROOT=pathlib.Path('/mnt/data/901_E516_20261117')
SOURCE=pathlib.Path('/mnt/data/901_E515_20261117/901_E515_37组独立相衬盲审包_20261117.zip')
OUTPUT=ROOT/'901_E516_37例离线独立显微读片工作台_不含解盲表_20261117.zip'
EXPECTED=''
def sha(b):return hashlib.sha256(b).hexdigest()
with zipfile.ZipFile(SOURCE) as z:
 assert z.testzip() is None
 source_hash=sha(SOURCE.read_bytes())
 orig=json.loads(z.read('blind_media/E515_public_blind_case_manifest.json'))
 cases=orig['ordered_cases']
 assert len(cases)==37 and len({c['case_id'] for c in cases})==37
 ids=[c['case_id'] for c in cases]
 images={}
 for c in cases:
  assert len(c['images'])==2
  for p in c['images']:
   assert p.startswith(c['case_id']) and p.endswith('.png') and ('full_original_PhC' in p or 'focus_original_PhC' in p)
   images[p]=z.read('blind_media/'+p)
 assert len(images)==74
 protocol=z.read('blind_media/E515_blind_review_protocol.md')
 blank=z.read('blind_media/E515_blank_independent_reviewer_ratings.csv')
 assert blank.decode().count('\n')==38
 assert orig['offsets_real_TIFF_zero_based_relative_to_unknown_event']==[-2,-1,0,1,2,4]
 manifest={'packet_id':'901-E516-independent-human-reader-offline','origin':'E515 preprocessed source PhC, time-corrected E513',
           'source_archive_sha256':source_hash,'case_count':37,'image_count':74,
           'image_sha256':{k:sha(v) for k,v in sorted(images.items())},
           'no_expert_annotations_yet':True,'no_source_event_classifiers_or_MoMA_masks':True,
           'no_private_unblind_key_in_package':True,'no_external_network_requests':True,
           'time_offsets':[-2,-1,0,1,2,4],
           'no_claim_of_biological_truth_or_independence_without_external_expert':True}
 with (ROOT/'E516_blind_media_source_manifest.json').open('w',encoding='utf8') as out:json.dump(manifest,out,ensure_ascii=False,indent=2)

js_cases=json.dumps([{'id':c['case_id'], 'full':c['images'][0] if 'full' in c['images'][0] else c['images'][1],
                  'focus':c['images'][0] if 'focus' in c['images'][0] else c['images'][1]} for c in cases],ensure_ascii=False)
css='''
:root{--bg:#10151f;--pane:#1a2230;--txt:#edf4fc;--sub:#9baebf;--line:#334154;--accent:#8ce2e0}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--txt);font:15px system-ui,-apple-system,Segoe UI,sans-serif}
header{padding:16px 24px;border-bottom:1px solid var(--line);background:#171f2b;display:flex;flex-wrap:wrap;gap:16px;align-items:center;justify-content:space-between}
h1{font-size:21px;font-weight:650;margin:0}.tag{font-size:12px;color:#081419;background:var(--accent);padding:4px 7px;border-radius:6px}
main{padding:16px 24px;max-width:1650px;margin:auto}.row{display:flex;gap:12px;flex-wrap:wrap;align-items:center}.stack{display:grid;gap:12px}
.muted{color:var(--sub)}.panel{background:var(--pane);border:1px solid var(--line);border-radius:10px;padding:15px}
button,select,input,textarea{font:inherit}button{background:#2c3b50;color:var(--txt);border:1px solid #48627d;border-radius:7px;padding:9px 14px;cursor:pointer}button:hover{background:#344a62}
button.primary{background:#256e71;color:#fff}.button-disabled{opacity:.5}label{display:block}.small{font-size:12px}
input[type=text],input[type=number],textarea,select{color:var(--txt);background:#101b28;border:1px solid #3d536c;padding:9px;border-radius:6px}
textarea{min-height:82px;width:100%}#imgwrap{overflow:auto;background:#0d121a;padding:6px;border:1px solid #26384c;border-radius:9px;max-height:620px}
#mainImage{max-width:none;display:block;transform-origin:top left;image-rendering:pixelated;}
fieldset{border:1px solid var(--line);border-radius:8px;padding:11px}legend{padding:0 6px;color:var(--sub)}
.radios{display:flex;gap:14px;flex-wrap:wrap}.radios label{cursor:pointer;display:flex;gap:6px;align-items:center;padding:6px;border-radius:6px}.radios label:hover{background:#27354a}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:12px}
#counter{font-size:16px;font-weight:600}.hint{font-size:12px;color:var(--sub)}
#notice{color:#8ce2e0}a{color:#a8ece9}.code{font-family:ui-monospace,SFMono-Regular,monospace}
@media(max-width:680px){main{padding:10px}header{padding:10px}.radios{display:block}}
'''
html_doc='''<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>901 E516 | Independent phase-contrast review</title><style>__STYLE__</style></head><body>
<header><div><h1>901 · Independent phase-contrast reviewer</h1><div class="muted small">37 blinded source-preprocessed bacterial image sequences · six relative timepoints per case · no segmentation overlays</div></div><span class="tag">OFFLINE · NO DATA UPLOAD</span></header>
<main class="stack">
<section class="panel row"><label>Anonymous reviewer ID <input id="reviewer" type="text" placeholder="e.g., A01" maxlength="80"/></label><label>Case <select id="caseSelector" aria-label="Select blinded case"></select></label><button id="prev">Previous ←</button><button id="next">Next →</button><span id="counter"></span><span id="notice" role="status"></span></section>
<section class="panel"><div class="row"><strong>Preprocessed phase-contrast source pixels</strong><label><input type="radio" name="view" value="full" checked/> Full lane</label><label><input type="radio" name="view" value="focus"/> Fixed focus view</label><label>Zoom <input id="zoom" type="range" min="60" max="160" value="90" step="10"/> <output id="zoomVal">90%</output></label></div>
<p class="hint">Panels are real TIFF indices at relative −2, −1, 0, +1, +2, +4. Time 0 is a concealed observation anchor, NOT verified biological division. Scale, ROI, and nearest-neighbor enlargement were frozen before review. No MoMA masks or software answers are shown.</p>
<div id="imgwrap"><img id="mainImage" alt="Original preprocessed phase microscopy sequence"/></div>
</section>
<section class="panel stack">
<fieldset><legend>Independent reviewer outcome — required, no default</legend><div id="classGroup" class="radios">
<label><input type="radio" name="classification" value="VISIBLE_DIVISION"/> Visible physical mother → two persistent daughters</label>
<label><input type="radio" name="classification" value="VISIBLE_CONTINUITY_NO_DIVISION"/> Visible continuity without division</label>
<label><input type="radio" name="classification" value="TIMING_SHIFT"/> Visible division at shifted time</label>
<label><input type="radio" name="classification" value="UNRESOLVED"/> Unresolved — uncertain/obscured</label>
</div></fieldset>
<div class="grid"><label>Image quality<br/><select id="quality"><option value="">Choose…</option><option value="good">Good</option><option value="fair">Fair</option><option value="poor">Poor</option></select></label>
<label>Identifiable mother before anchor<br/><select id="mother"><option value="">Choose…</option><option value="yes">Yes</option><option value="no">No</option><option value="uncertain">Uncertain</option></select></label>
<label>Two daughters persistent to +2 or +4<br/><select id="daughters"><option value="">Choose…</option><option value="yes">Yes</option><option value="no">No</option><option value="uncertain">Uncertain</option></select></label>
<label>Evidence relative offsets<br/><input type="text" id="offsets" placeholder="e.g. −1, +1, +2, +4"/></label>
<label>Evidence y region / field<br/><input type="text" id="region" placeholder="e.g. y≈130–170; full lane"/></label></div>
<label>Clipping, overlap, uncertain lineage, limitations<br/><textarea id="uncertainty" placeholder="Describe why the outcome is supported or unresolved. Be specific."></textarea></label>
<label>Direct observations (not software predictions)<br/><textarea id="notes" placeholder="Record image evidence and two-daughter persistence; do not infer molecular information."></textarea></label>
<div class="row"><button class="primary" id="saveNext">Save case and next</button><button id="exportCsv">Export CSV</button><button id="exportJson">Export JSON</button><label class="small">Restore exported JSON <input id="restore" type="file" accept=".json,application/json"/></label><span id="progress" class="muted"></span></div>
<p class="hint">Data remain in this browser on your computer, with local recovery storage if permitted. Export CSV/JSON regularly; do not rely only on browser storage. Two reviewers should work independently and keep their files separate until both have submitted ratings. Do not open the private unblinding key.</p>
</section>
<section class="panel"><strong>Scientific scoring boundary</strong><p class="muted">A visible bright line alone is insufficient. A division requires one physically identifiable mother and two distinct, persistently trackable daughter bodies. A missed boundary is not proof of no division. Edge clipping, sub-resolution features, uncertain ancestry, and overlap should be marked <span class="code">UNRESOLVED</span>. The 37 cases overlap in time and come from two lanes of a single experiment; they are not independent biological replicates.</p></section>
</main>
<script>
'use strict';
const CASES=__CASES__;
const STORAGE='901-E516-blind-human-review-v1';
const FIELDS=['classification','quality','mother','daughters','offsets','region','uncertainty','notes'];
const CSV_FIELDS=['blind_case_id','reviewer_anonymous_id','classification','image_quality_good_fair_poor','distinct_mother_seen_before_t','two_persistent_daughters_visible_by_t_plus_2_or_4','evidence_relative_frame_offsets','evidence_Y_or_full_ROI','uncertainty_or_clipping_reasons','free_text_observation'];
const ids=new Set(CASES.map(x=>x.id));
const memo={};
try {const old=JSON.parse(localStorage.getItem(STORAGE)||'{}');
if(old && typeof old==='object') {for(const [k,v] of Object.entries(old.ratings||{}))if(ids.has(k)&&v&&typeof v==='object')memo[k]=v;
if(typeof old.reviewer==='string')document.getElementById('reviewer').value=old.reviewer;}}catch(_){}
const $=id=>document.getElementById(id);
const selector=$('caseSelector');CASES.forEach((c,i)=>{const opt=document.createElement('option');opt.value=String(i);opt.textContent=`${i+1}. ${c.id}`;selector.appendChild(opt);});
let index=0;
function notice(msg){$('notice').textContent=msg;}
function recordNow(){const c=CASES[index];const r=memo[c.id]||{};for(const field of FIELDS){let val='';if(field==='classification')val=document.querySelector('input[name=classification]:checked')?.value||'';else val=$(field).value||'';r[field]=val.trim();}memo[c.id]=r;persist();updateProgress();}
function persist(){try{localStorage.setItem(STORAGE,JSON.stringify({reviewer:$('reviewer').value,ratings:memo}));}catch(_) {}}
function loadAt(i){index=Math.min(Math.max(0,i),CASES.length-1);selector.value=String(index);const c=CASES[index];const rec=memo[c.id]||{};
for(const radio of document.querySelectorAll('input[name=classification]'))radio.checked=(radio.value===rec.classification);
for(const f of FIELDS.filter(x=>x!=='classification'))$(f).value=rec[f]||'';
renderImage();$('counter').textContent=`Case ${index+1}/${CASES.length}`;updateProgress();notice('');}
function renderImage(){const c=CASES[index];const type=document.querySelector('input[name=view]:checked').value;
$('mainImage').src='media/'+(type==='focus'?c.focus:c.full);$('mainImage').style.width=$('zoom').value+'%';$('zoomVal').textContent=$('zoom').value+'%';$('imgwrap').scrollTo(0,0);}
function updateProgress(){let filled=CASES.filter(c=>memo[c.id]?.classification).length;$('progress').textContent=`${filled}/37 cases have a classification`}
function go(i){recordNow();loadAt(i);}
function validateReview(){recordNow();const rid=$('reviewer').value.trim();if(!rid){notice('Enter anonymous reviewer ID before export.');$('reviewer').focus();return false;}return true;}
function rows(){const rid=$('reviewer').value.trim();return CASES.map(c=>{const v=memo[c.id]||{};return {blind_case_id:c.id,reviewer_anonymous_id:rid,
classification:v.classification||'',image_quality_good_fair_poor:v.quality||'',distinct_mother_seen_before_t:v.mother||'',
two_persistent_daughters_visible_by_t_plus_2_or_4:v.daughters||'',evidence_relative_frame_offsets:v.offsets||'',
evidence_Y_or_full_ROI:v.region||'',uncertainty_or_clipping_reasons:v.uncertainty||'',free_text_observation:v.notes||''};});}
function download(txt,name,type){const blob=new Blob([txt],{type});const u=URL.createObjectURL(blob);const el=document.createElement('a');el.href=u;el.download=name;document.body.appendChild(el);el.click();el.remove();URL.revokeObjectURL(u);}
function quote(s){const v=String(s??'');return /[",\r\n]/.test(v)?'"'+v.replaceAll('"','""')+'"':v;}
$('exportCsv').onclick=()=>{if(!validateReview())return;const out=[CSV_FIELDS.join(','),...rows().map(r=>CSV_FIELDS.map(f=>quote(r[f])).join(','))].join('\r\n');download('\ufeff'+out,`E516_review_${$('reviewer').value.trim()}.csv`,'text/csv;charset=utf-8');notice('CSV exported with all 37 rows; unfilled cases remain blank.');};
$('exportJson').onclick=()=>{if(!validateReview())return;download(JSON.stringify({packet:'901-E516-blind-phase-reader',date_local:new Date().toISOString(),reviewer:$('reviewer').value.trim(),ratings:memo},null,2),`E516_review_${$('reviewer').value.trim()}.json`,'application/json');notice('JSON review exported.');};
$('restore').onchange=async function(){const f=this.files?.[0];if(!f)return;try{const v=JSON.parse(await f.text());if(v.packet!=='901-E516-blind-phase-reader')throw Error('Not an E516 review export');for(const [id,r] of Object.entries(v.ratings||{}))if(ids.has(id)&&r&&typeof r==='object')memo[id]=r;if(v.reviewer)$('reviewer').value=String(v.reviewer);persist();loadAt(index);notice('Review restored.');}catch(e){notice('Cannot restore: '+e.message);}};
selector.onchange=()=>go(Number(selector.value));$('prev').onclick=()=>go(index-1);$('next').onclick=()=>go(index+1);$('saveNext').onclick=()=>go(index+1);
$('zoom').oninput=renderImage;document.querySelectorAll('input[name=view]').forEach(r=>r.onchange=renderImage);
$('reviewer').onchange=persist;
document.querySelectorAll('input[name=classification]').forEach(x=>x.onchange=recordNow);
for(const x of FIELDS.filter(y=>y!=='classification'))$(x).oninput=recordNow;
document.addEventListener('keydown',e=>{if(['INPUT','TEXTAREA','SELECT'].includes(document.activeElement?.tagName))return;if(e.key==='ArrowRight')go(index+1);if(e.key==='ArrowLeft')go(index-1);if(['1','2','3','4'].includes(e.key)){const radio=document.querySelectorAll('input[name=classification]')[Number(e.key)-1];if(radio){radio.checked=true;recordNow();}}});
loadAt(0);
</script></body></html>'''.replace('__STYLE__',css).replace('__CASES__',js_cases)
# JS syntax requires literal escape characters, not Python's decoded CR/LF,
# inside JavaScript regexes and CSV-joining string literals.
html_doc=html_doc.replace('/[",\r\n]/',r'/[",\r\n]/')
html_doc=html_doc.replace("join('\r\n')",r"join('\r\n')")
assert 'https://' not in html_doc and 'fetch(' not in html_doc and 'XMLHttpRequest' not in html_doc
readme='''901 E516 | Offline blind phase-contrast microscopy review workstation

1. Extract the whole ZIP to a folder, preserving the media/ directory.
2. Open E516_offline_review_workbench.html in a modern browser. No internet connection or account is needed; no data leaves the device automatically.
3. Enter an anonymous reviewer ID, navigate the 37 pseudonymized cases and record all observations. Use UNRESOLVED when physical mother-daughter ancestry cannot be established.
4. Export both CSV and JSON; save to a safe location. Separate expert reviewers must not view each other's responses before both submit.
5. The PRIVATE unblinding key is stored separately by the study owner and must not be given to reviewers.

This packet contains original preprocessed phase images, not raw camera frames or MoMA masks. It contains no biological truths or software consensus labels. An external microscope reviewer is required for biological adjudication. No software inference has been automatically carried out merely by opening this HTML tool.
'''

fixed=(1980,1,1,0,0,0)
with zipfile.ZipFile(OUTPUT,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=7) as w:
    items={'E516_offline_review_workbench.html':html_doc.encode(),
           'README_INDEPENDENT_REVIEW.txt':readme.encode(),
           'E515_blind_review_protocol.md':protocol,
           'E515_blank_independent_reviewer_ratings.csv':blank,
           'E516_public_media_source_hash_manifest.json':json.dumps(manifest,ensure_ascii=False,indent=2).encode()}
    items.update({'media/'+n:data for n,data in images.items()})
    for name,dat in sorted(items.items()):
        info=zipfile.ZipInfo(name,date_time=fixed);info.compress_type=zipfile.ZIP_DEFLATED
        info.create_system=3
        w.writestr(info,dat,compress_type=zipfile.ZIP_DEFLATED,compresslevel=7)
with zipfile.ZipFile(OUTPUT) as z:
 assert z.testzip() is None
 assert len(z.namelist())==79, len(z.namelist())
 assert all(sha(z.read('media/'+n))==sha(data) for n,data in images.items())
 assert 'PRIVATE' not in '\n'.join(z.namelist())
 print('E516_OFFLINE_REVIEW_WORKBENCH_PASS',json.dumps({'cases':len(cases),'images':len(images),'zip_bytes':OUTPUT.stat().st_size,'zip_sha256':sha(OUTPUT.read_bytes()),'HTML_bytes':len(html_doc),'source_zip_sha256':source_hash}))