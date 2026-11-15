#!/usr/bin/env python3
"""901 E492 authentic TrIPP protocol Table 1 optical source audit.
No molecular absolute number, parent->daughter MI, or technical repeat inferred.
Inputs: PMC8165573 XML and published Table 1 in original raw form.
"""
from pathlib import Path
from collections import Counter
import csv, hashlib, json, statistics, zipfile, io, xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
SOURCE_SHA = '6a59cfb29e1e933577a40b99de88e5ee4ed7cdf1e6b0355e371d8a08fd9b6667'
FULL_BUNDLE_SHA='2c03b1356cd696dd7d06f0a820463f388942e08c702fe8aa5f32e02848e4c9b8'
INNER_DATA_S3_SHA='7dc3bce68688852585c11e2dff07351c874cc3cf118c7505ef5358fc674ebe47'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    xml=ROOT/'PMC8165573_fullTextXML.xml'
    source_csv=ROOT/'PMC8165573_table_tbl1.csv'
    bundle=ROOT/'PMC8165573_DataS3_mmc3_original.zip'
    assert digest(xml)==SOURCE_SHA
    assert digest(bundle)==FULL_BUNDLE_SHA
    root=ET.parse(xml).getroot()
    tables=[t for t in root.findall('.//table-wrap') if t.get('id')=='tbl1']
    assert len(tables)==1
    caption=' '.join(''.join(tables[0].find('caption').itertext()).split())
    assert 'mother cell and its two daughters' in caption
    with zipfile.ZipFile(bundle) as z:
        assert z.testzip() is None
        inner=z.read('mmc3.zip')
        assert hashlib.sha256(inner).hexdigest()==INNER_DATA_S3_SHA
        with zipfile.ZipFile(io.BytesIO(inner)) as innerzip:
            assert innerzip.testzip() is None
            names=innerzip.namelist()
            assert len(names)==13 and 'Readme_before_running_TrIPPed.txt' in names
            readme=innerzip.read('Readme_before_running_TrIPPed.txt').decode('utf8','replace')
            assert 'csv output files from BudJ' in readme
            (ROOT/'PMC8165573_DataS3_13scripts_original.zip').write_bytes(inner)
    with source_csv.open(newline='',encoding='utf8') as f: raw=list(csv.DictReader(f))
    assert len(raw)==145
    data={}
    for r in raw:
        c=int(r['Cell (1:Mother: 2:Daughter 1; 3:Daughter 2)'])
        t=float(r['Time (min)'])
        obs={'red_mean':float(r['FL1 (red) FL mean']), 'green_mean':float(r['FL2 (green) FL mean']),
             'x':float(r['x (coordinate for the position of the mother cell in the image)']),
             'y':float(r['y (coordinate of the position of the mother cell in the image)'])}
        data.setdefault(c,{})[t]=obs
    assert {c:len(x) for c,x in data.items()}=={1:61,2:50,3:34}
    assert len(raw)==sum(map(len,data.values()))
    births={2:min(data[2]),3:min(data[3])}
    assert births=={2:93.5,3:229.5}
    for c,v in data.items():
        ts=sorted(v)
        assert all(abs(ts[k]-ts[k-1]-8.5)<1e-9 for k in range(1,len(ts)))
    rows=[]
    for daughter in (2,3):
        first=births[daughter]
        upper=births[3] if daughter==2 else float('inf')
        for t in sorted(data[1].keys() & data[daughter].keys()):
            if not first<=t<upper:continue
            m=data[1][t];d=data[daughter][t]
            rows.append({'daughter':daughter,'time_min':t,'red_mother_mean_raw_AU':m['red_mean'],
                         'red_daughter_mean_raw_AU':d['red_mean'], 'green_mother_mean_raw_AU':m['green_mean'],
                         'green_daughter_mean_raw_AU':d['green_mean'],
                         'raw_red_mean_fraction_mother':m['red_mean']/(m['red_mean']+d['red_mean']),
                         'raw_green_mean_fraction_mother':m['green_mean']/(m['green_mean']+d['green_mean'])})
    assert len(rows)==50
    assert sum(x['daughter']==2 for x in rows)==16
    assert sum(x['daughter']==3 for x in rows)==34
    summary={}
    for d in (2,3):
        sub=[x for x in rows if x['daughter']==d]
        r=[x['raw_red_mean_fraction_mother'] for x in sub]
        g=[x['raw_green_mean_fraction_mother'] for x in sub]
        summary[str(d)]={'n_synchronous_frames':len(sub),'first_min':sub[0]['time_min'],
                         'last_min':sub[-1]['time_min'],
                         'median_raw_red_mean_fraction_mother':statistics.median(r),
                         'median_raw_green_mean_fraction_mother':statistics.median(g),
                         'raw_red_mean_fraction_min_max':[min(r),max(r)],
                         'raw_green_mean_fraction_min_max':[min(g),max(g)]}
    assert abs(summary['2']['median_raw_red_mean_fraction_mother']-0.5141911186524042)<1e-12
    assert abs(summary['3']['median_raw_red_mean_fraction_mother']-0.4882195030483316)<1e-12
    with (ROOT/'E492_raw_synchronous_optical_means.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    result={
      'stage':'E492','source':'Auboiron et al., Star Protocols 2021 TrIPP protocol DOI 10.1016/j.xpro.2021.100557',
      'protocol_pmcid':'PMC8165573','study_pmcid':'PMC7889982',
      'protocol_jats_source_sha256':SOURCE_SHA,'full_official_supplement_bundle_sha256':FULL_BUNDLE_SHA,
      'dataS3_inner_zip_sha256':INNER_DATA_S3_SHA,'dataS3_source_scripts_n':len(names),
      'source_table_title':caption,'n_raw_rows':145,'cell_ids':{'mother':'1','daughter_1':'2','daughter_2':'3'},
      'n_rows_by_cell':{'mother':61,'daughter1':50,'daughter2':34},
      'first_visible_daughter_birth_times_min':{'daughter1':93.5,'daughter2':229.5},
      'time_grid_step_min':8.5,'synchronous_observation_count_for_optical_comparison':len(rows),
      'optical_descriptive_summaries':summary,
      'source_protein':'H4-Dendra2','channel_units':'raw per-cell mean fluorescence, uncorrected, arbitrary optical units',
      'background_subtraction_independent_reconstruction':'NOT_AVAILABLE',
      'per_cell_total_fluorescence':'NOT_AVAILABLE',
      'per_cell_absolute_protein_copy':'NOT_IDENTIFIED',
      'mother_to_daughter_molecular_transfer_MI_bit':'NOT_IDENTIFIED',
      'N_independent_mother_lineages':1,
      'no_population_mi_inference':True,
      'warning':'Raw M/(M+D) ratios of per-cell mean intensities are NOT fractional conserved molecules, not after-background measured red protein transfer, not independent replicates, and not comparable to 136-mother Figure 2 distribution.',
      'source_qualification':'ONE_AUTHENTIC_LINEAGE_WITH_RAW_OPTICAL_MEAN_ONLY',
      'publication_status':'PUBLICATION_NO_GO'
    }
    (ROOT/'E492_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({k:result[k] for k in ('n_raw_rows','n_rows_by_cell','time_grid_step_min','synchronous_observation_count_for_optical_comparison','optical_descriptive_summaries','source_qualification')},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
