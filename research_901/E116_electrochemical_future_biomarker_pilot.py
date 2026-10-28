# 901 E116 independent-source electrochemical time-series pilot.
# Source original checksum pinned; no author-supplied scripts executed.
import os,io,math,hashlib,zipfile,urllib.request,xml.etree.ElementTree as ET,json,random,statistics
from pathlib import Path
from collections import Counter, defaultdict
URL='https://zenodo.org/records/19386720/files/Supplementary%20data%20and%20code.zip?download=1'
MD5='9ceba7a9ad2c868eaa9ef422e50f091e'
NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
import sys
if len(sys.argv)>2 and sys.argv[1]=='--input':
    archive=Path(sys.argv[2]).read_bytes()
else:
    req=urllib.request.Request(URL,headers={'User-Agent':'901-E116-read-only-plant-source'})
    with urllib.request.urlopen(req,timeout=60) as r: archive=r.read(14_000_001)
assert len(archive)<=14_000_000 and hashlib.md5(archive).hexdigest()==MD5
with zipfile.ZipFile(io.BytesIO(archive)) as z:
    xnames=[n for n in z.namelist() if n.endswith('/plant_data.xlsx')]
    assert len(xnames)==1
    data=z.read(xnames[0]); input_sha=hashlib.sha256(data).hexdigest()
    paths=z.namelist(); src_script=z.read(next(n for n in paths if n.endswith('/t0_process_features.py'))).decode('utf-8','replace')
plants=defaultdict(dict); groups={}; nrows=0; duplicate=0; skipped=0
with zipfile.ZipFile(io.BytesIO(data)) as z, z.open('xl/worksheets/sheet1.xml') as f:
    shared_root=ET.fromstring(z.read('xl/sharedStrings.xml'))
    shared_strings=[''.join(t.text or '' for t in si.iter() if t.tag==NS+'t') for si in shared_root.findall(NS+'si')]
    for event,row in ET.iterparse(f,events=('end',)):
        if row.tag!=NS+'row': continue
        rowid=int(row.attrib.get('r','0'))
        if rowid<=1: row.clear(); continue
        v={}
        for c in row.findall(NS+'c'):
            col=''.join(q for q in c.get('r','') if q.isalpha())
            e=c.find(NS+'v')
            if e is None or e.text is None: continue
            v[col]=e.text
        row.clear()
        try:
            t=int(round(float(v['A'])));raw_p=int(float(v['B'])); raw_g=v['C']; group=shared_strings[int(raw_g)] if raw_g.isdigit() else raw_g; p=f'{group}:{raw_p}'
            k=float(v['D']);ph=float(v['E']);h=float(v['F'])
            if not all(math.isfinite(q) for q in (k,ph,h)): raise ValueError('nonfinite')
            if t in plants[p]: duplicate+=1
            plants[p][t]=(k,ph,h)
            if p in groups and groups[p]!=group: raise ValueError('group drift')
            groups[p]=group;nrows+=1
        except (ValueError,KeyError,TypeError): skipped+=1
nplant=len(plants)
assert nplant==60 and duplicate==0 and skipped==0,(nplant,duplicate,skipped)
bygroup=defaultdict(list)
for p,g in groups.items():bygroup[g].append(p)
def pct(x,q):
    a=sorted(x)
    return a[min(len(a)-1,max(0,int(round(q*(len(a)-1)))))] if a else None
monos=[];second_d=[];tensec=[];missing=[]
for p,trace in plants.items():
    times=sorted(trace)
    missing.append(1-sum(times[i+1]-times[i]==10 for i in range(len(times)-1))/max(1,len(times)-1))
    dh=[trace[times[i+1]][2]-trace[times[i]][2] for i in range(len(times)-1) if times[i+1]-times[i]==10]
    monos.append(sum(x<0 for x in dh)/len(dh) if dh else None)
    second_d.extend(abs(dh[i+1]-dh[i]) for i in range(len(dh)-1))
    tensec.extend(dh)
allgroups=sorted(bygroup)
foldmap={p:i%5 for g in allgroups for i,p in enumerate(sorted(bygroup[g]))}
start_t=600;horizon=3600;stride=600
records=[]
for p,trace in plants.items():
    for t in range(start_t, max(trace)-horizon+1, stride):
        if t not in trace or t-600 not in trace or t+horizon not in trace:continue
        k,ph,h=trace[t]; kp,php,hp=trace[t-600]; future=trace[t+horizon][2]
        if not all(math.isfinite(v) for v in (k,ph,h,kp,php,hp,future)):continue
        extras=(k,ph,k-kp,ph-php)
        group=groups[p]
        base=(h,h-hp,t/28800)+tuple(1.0 if group==x else 0.0 for x in allgroups[1:])
        records.append({'p':p,'g':group,'t':t,'y':future,'persist':h,'base':base,'extra':extras})
def solve(a,b):
    n=len(b);a=[list(aa)+[bb] for aa,bb in zip(a,b)]
    for i in range(n):
        j=max(range(i,n),key=lambda k:abs(a[k][i]))
        a[i],a[j]=a[j],a[i];v=a[i][i]
        if abs(v)<1e-10: raise ArithmeticError('ill-conditioned normal equation')
        a[i]=[u/v for u in a[i]]
        for j in range(i+1,n):
            fac=a[j][i]
            if fac:a[j]=[u-fac*w for u,w in zip(a[j],a[i])]
    x=[0.]*n
    for i in range(n-1,-1,-1):x[i]=a[i][-1]-sum(a[i][j]*x[j] for j in range(i+1,n))
    return x
def fit(train,fn,penalty=3.0):
    raw=[fn(x) for x in train];n=len(raw);k=len(raw[0]);mean=[];sd=[]
    for j in range(k):
        a=[r[j] for r in raw]; m=sum(a)/n; s=(sum((v-m)**2 for v in a)/n)**.5
        mean.append(m);sd.append(s if s>1e-9 else 1)
    xmat=[[1]+[(v-mean[j])/sd[j] for j,v in enumerate(row)] for row in raw]
    L=k+1
    mat=[[0.]*L for _ in range(L)]; vec=[0.]*L
    for x,rec in zip(xmat,train):
        for i in range(L):
            vec[i]+=x[i]*rec['y']
            for j in range(i+1):mat[i][j]+=x[i]*x[j]
    for i in range(L):
        for j in range(i):mat[j][i]=mat[i][j]
    for i in range(1,L):mat[i][i]+=penalty
    coef=solve(mat,vec)
    return (mean,sd,coef)
def predict(model,x):
    m,sd,c=model
    return c[0]+sum(c[j+1]*(a-m[j])/sd[j] for j,a in enumerate(x))
placebo={}
for g,pids in bygroup.items():
    a=sorted(pids)
    for i,p in enumerate(a):placebo[p]=a[(i+3)%len(a)] if len(a)>1 else p
by_p={};fold_results=[]
for fold in range(5):
    train=[a for a in records if foldmap[a['p']]!=fold]
    test=[a for a in records if foldmap[a['p']]==fold]
    if not train or not test:continue
    f0=fit(train,lambda x:x['base'])
    f1=fit(train,lambda x:x['base']+x['extra'])
    def sham(rec):
        q=placebo[rec['p']];tr=plants[q];t=rec['t']
        if t not in tr or t-600 not in tr:return rec['extra']
        kk,pp,_=tr[t];kp,ph_p,_=tr[t-600]
        return (kk,pp,kk-kp,pp-ph_p)
    fsh=fit(train,lambda x:x['base']+sham(x))
    errors={'persist':[],'M0':[],'M1':[],'placebo':[]}
    for rec in test:
        p=rec['p'];y=rec['y'];vals={'persist':rec['persist'], 'M0':predict(f0,rec['base']), 'M1':predict(f1,rec['base']+rec['extra']),'placebo':predict(fsh,rec['base']+sham(rec))}
        by_p.setdefault(p,{'g':rec['g'],'e':{name:[] for name in vals}})
        for name,val in vals.items():
            err=abs(val-y);errors[name].append(err);by_p[p]['e'][name].append(err)
    fold_results.append({'fold':fold,'train_plants':len({r['p'] for r in train}),'test_plants':len({r['p'] for r in test}),'train_windows':len(train),'test_windows':len(test),'mae_by_model':{name:sum(err)/len(err) for name,err in errors.items()}})
perplant=[]
for p,rec in sorted(by_p.items()):
    m={name:sum(q)/len(q) for name,q in rec['e'].items() if q}
    perplant.append({'plant_id':p,'group':rec['g'],'n_eval':len(rec['e']['M0']),'mean_abs_err':m,'gain_M0_minus_M1':m['M0']-m['M1']})
rng=random.Random(901116)
diff=[r['gain_M0_minus_M1'] for r in perplant]
obs=sum(diff)/len(diff) if diff else None
boot=[]
for _ in range(2000):boot.append(sum(diff[rng.randrange(len(diff))] for j in range(len(diff)))/len(diff))
boot.sort();ci=(boot[49],boot[1949]) if boot else None
means={name:statistics.mean(r['mean_abs_err'][name] for r in perplant) for name in ('persist','M0','M1','placebo')}
report={'stage':'E116','study':'Zhou lettuce 60 plant K/pH/H2O2 source','origin':'https://zenodo.org/records/19386720','original_archive_md5':MD5,'xlsx_sha256':input_sha,'source_scripts_executed':False,'columns':['Time (s)','Plant_ID','Group','K_Ratio','pH_Value','H2O2_Conc'],
 'observed_rows':nrows,'duplicates_plant_time':duplicate,'unparsable_rows':skipped,'unique_plants':nplant,'physical_parent_key':'Group + Plant_ID (Plant_ID reused across 6 groups)',
 'plants_by_group':{g:len(ids) for g,ids in bygroup.items()},'trace_sizes':{str(p):len(v) for p,v in plants.items()},
 'regular_10s_gap_missing_fraction_median':pct(missing,.5),'H2O2_first_diff_fraction_negative_median_by_plant':pct([x for x in monos if x is not None],.5),
 'H2O2_abs_second_difference_q50':pct(second_d,.5),'H2O2_abs_second_difference_q95':pct(second_d,.95),
 'H2O2_first_difference_q05_q50_q95':[pct(tensec,q) for q in (.05,.5,.95)],
 'prediction':'future H2O2_Conc at t+3600 seconds','not_Bombadil_membrane_potential':True,'not_drought_ground_truth':True,
 'group_held_out':'5fold plant-level split stratified by published group','model_M0':'current H2O2, prior 10min H2O2 change, clock, known group',
 'model_M1':'M0 plus concurrent K_Ratio, pH and each preceding 10min change',
 'placebo':'M0 plus same-treatment other-plant K/pH features','anchor_stride_seconds':stride,
 'pred_horizon_seconds':horizon,'n_eligible_plant_windows':len(records), 'folds':fold_results,
 'mean_MAE_by_plant':means,'mean_plant_gain_M0_minus_M1':obs,'gain_bootstrap_95pct_plant_level':ci,
 'perplant':perplant,
 'limitations':['assays are distinct electrochemical biomarkers, not plant electrical membrane-potential recordings','future H2O2 is a measured downstream biomarker, not a visually confirmed drought/wilt endpoint', 'time series may be transformed or unusually smooth; true raw-instrument independent provenance not reverified','one species, known treatment and controlled conditions; results not cross-site or cross-species generalization']}
Path('901_e116_pilot_results.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print('E116_RESULTS',json.dumps({k:report[k] for k in report if k not in ('perplant','trace_sizes','folds')},ensure_ascii=False),flush=True)
print('E116_FOLDS',json.dumps(fold_results),flush=True)
print('E116_SHA',hashlib.sha256(Path('901_e116_pilot_results.json').read_bytes()).hexdigest(),flush=True)
