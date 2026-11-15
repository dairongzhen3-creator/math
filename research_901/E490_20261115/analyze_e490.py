#!/usr/bin/env python3
"""E490 published-rule FLCS observability audit, no refitting or hypothesis selection.
Inputs are frozen official XLSX-derived CSV copied from E489; source SHA pinned.
This counts threshold agreement, operational uncertainty, and gives two latent-state
counterexamples. Does NOT estimate true molecule-copy transmission.
"""
from pathlib import Path
from collections import Counter
import csv, math, statistics, hashlib, json

ROOT = Path(__file__).resolve().parent
SRC = ROOT / 'FLCS_lacInduction_MotherMachine_source_cells.csv'
JOIN = ROOT / 'E489_exact_source_join_all.csv'
RMD_SOURCE_URL = 'https://github.com/julou/MoM_lacInduction/blob/plos/src/MoM_lacInduction_FLIM.Rmd'
PINNED_CODE_BLOB_SHA = 'bfc059418ebcc88f326c9ba722abb1a6d5ef52fe'
AUTHOR_FIGURE_NCPM_THRESHOLD = 2.55
AUTHOR_N_MAX = 20
AUTHOR_TAU_MIN_MS = 3
AUTHOR_TAU_MAX_MS = 30
CONFOCAL_W_UM = 0.33

def num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None

def entropy(p):
    return -sum(x * math.log2(x) for x in p if x > 0)

def mi_rows(rows):
    n = len(rows)
    joint = Counter(rows)
    nx = Counter(x for x, y in rows)
    ny = Counter(y for x, y in rows)
    return sum((v/n)*math.log2((v*n)/(nx[x]*ny[y])) for (x,y),v in joint.items())

def main():
    with SRC.open(encoding='utf-8-sig', newline='') as fh:
        raw = list(csv.reader(fh))
    assert len(raw) == 75
    cpm_good = [num(row[14]) for row in raw[1:]]
    cpm_good = [x for x in cpm_good if x is not None]
    assert len(cpm_good) == 66
    cpm_med = statistics.median(cpm_good)
    assert cpm_med == 17.5
    cpm_cut = cpm_med * AUTHOR_FIGURE_NCPM_THRESHOLD
    outputs = []
    for row in raw[1:]:
        date, filename, gfp = str(row[0]), str(row[1]), str(row[2])
        n, nse, D, Dse, cpm, cpmse = [num(row[i]) for i in (10,11,12,13,14,15)]
        fit_complete = all(x is not None for x in (n,D,cpm))
        tau = (CONFOCAL_W_UM ** 2 * 1000 /(4 * D)) if (fit_complete and D > 0) else None
        pass_cpm = bool(fit_complete and cpm / cpm_med > AUTHOR_FIGURE_NCPM_THRESHOLD)
        pass_n = bool(fit_complete and n < AUTHOR_N_MAX)
        pass_tau = bool(fit_complete and tau is not None and AUTHOR_TAU_MIN_MS <= tau <= AUTHOR_TAU_MAX_MS)
        pred = ('1' if pass_cpm and pass_n and pass_tau else '-1' if pass_cpm else '0') if fit_complete else 'INCOMPLETE'
        se_bounded_n = se_bounded_cpm = se_bounded_tau = None
        if gfp == '1' and fit_complete:
            se_bounded_n = bool(nse is not None and n + nse < AUTHOR_N_MAX)
            se_bounded_cpm = bool(cpmse is not None and cpm - cpmse > cpm_cut)
            se_bounded_tau = bool(
                Dse is not None and D-Dse > 0 and
                AUTHOR_TAU_MIN_MS <= CONFOCAL_W_UM**2*1000/(4*(D+Dse)) and
                CONFOCAL_W_UM**2*1000/(4*(D-Dse)) <= AUTHOR_TAU_MAX_MS)
        outputs.append({
           'date':date,'ptu_filename_reference':filename,'published_gfp_label':gfp,
           'FLCS_particle_fit':n,'FLCS_particle_fit_se':nse,
           'FLCS_diffusion_um2_per_s':D,'FLCS_diffusion_se':Dse,
           'FLCS_CPM':cpm,'FLCS_CPM_se':cpmse,'derived_tau_ms':tau,
           'complete_for_rule':fit_complete,'passes_figure_CPM_cutoff':pass_cpm if fit_complete else None,
           'passes_published_particle_cutoff':pass_n if fit_complete else None,
           'passes_published_diffusion_time_range':pass_tau if fit_complete else None,
           'three_rule_reconstruction':pred,'one_reported_SE_CPM_margin':se_bounded_cpm,
           'one_reported_SE_particle_margin':se_bounded_n,'one_reported_SE_diffusion_margin':se_bounded_tau,
           'published_label_rule_agree':pred==gfp if fit_complete else None,
        })
    assert Counter(x['published_gfp_label'] for x in outputs)==Counter({'1':18,'0':45,'-1':5,'NA':6})
    complete=[x for x in outputs if x['complete_for_rule']]
    assert len(complete)==66
    assert all(x['published_label_rule_agree'] for x in complete)
    byflag={k:[x for x in outputs if x['published_gfp_label']==k] for k in ('1','0','-1','NA')}
    positives=byflag['1']
    assert sum(bool(x['one_reported_SE_CPM_margin']) for x in positives)==18
    assert sum(bool(x['one_reported_SE_particle_margin']) for x in positives)==18
    assert sum(bool(x['one_reported_SE_diffusion_margin']) for x in positives)==16
    assert all(x['FLCS_particle_fit']<12 for x in positives)
    assert all(not x['passes_figure_CPM_cutoff'] for x in byflag['0'] if x['complete_for_rule'])
    assert all(x['passes_figure_CPM_cutoff'] and (not x['passes_published_particle_cutoff'] or not x['passes_published_diffusion_time_range']) for x in byflag['-1'])
    # Published E489 56-cell exact joint sample; no new subgroup selection.
    with JOIN.open(encoding='utf-8-sig',newline='') as f:
        matched = list(csv.DictReader(f))
    pairs=[(int(x['gfp']),x['response']) for x in matched if x['gfp'] in {'0','1'} and x['response'] in {'S','L'}]
    assert len(pairs)==56
    counts=Counter(pairs)
    assert counts==Counter({(0,'S'):8,(0,'L'):34,(1,'S'):13,(1,'L'):1})
    observed_mi=mi_rows(pairs)
    assert abs(observed_mi-0.3347760077969534)<1e-14
    # Both constructed latent model states meet the one-sided zero-false-positive
    # condition D=1 => Z=1 while preserving the exact P(D,Y).
    candidate1=[(d,y) for d,y in pairs]  # perfect detection Z=D
    candidate2=[(1,y) for _,y in pairs]  # all molecules present Z=1 but variable sensitivity
    mi1,mi2=mi_rows(candidate1),mi_rows(candidate2)
    assert abs(mi1-observed_mi)<1e-14 and abs(mi2)<1e-14
    fields=list(outputs[0])
    with (ROOT/'E490_author_rule_per_cell.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(outputs)
    result={
      'project':'901','stage':'E490','date':'2026-11-15',
      'source':'Julou et al (2020) PLOS Biology DOI 10.1371/journal.pbio.3000952',
      'author_RMD_blob_sha':PINNED_CODE_BLOB_SHA,'author_RMD_url':RMD_SOURCE_URL,
      'source_workbook_sha256':'97f93074b72214bd8a20439ff508f6d8d47111ac115b6be6c06538457731212d',
      'source_csv_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),
      'source_rows':len(outputs),'source_labels':dict(Counter(x['published_gfp_label'] for x in outputs)),
      'fit_complete':len(complete),'fit_incomplete':len(outputs)-len(complete),
      'fit_complete_by_label':dict(Counter(x['published_gfp_label'] for x in complete)),
      'figure_normalized_CPM_boundary':AUTHOR_FIGURE_NCPM_THRESHOLD,
      'complete_case_CPM_median':cpm_med,'figure_implied_CPM_boundary':cpm_cut,
      'particle_limit_strict':AUTHOR_N_MAX,
      'diffusion_time_ms_bounds':[AUTHOR_TAU_MIN_MS,AUTHOR_TAU_MAX_MS],
      'diffusion_time_conversion':'(.33 um)^2/(4 * D_um2_s) *1000',
      'exact_agreement_count':sum(x['published_label_rule_agree'] is True for x in complete),
      'rule_condition_counts':{k:sum(x['passes_figure_CPM_cutoff'] is True for x in complete) for k in ('any',)},
      'cpm_pass_count':sum(x['passes_figure_CPM_cutoff'] is True for x in complete),
      'cpm_and_n_pass_count':sum(x['passes_figure_CPM_cutoff'] is True and x['passes_published_particle_cutoff'] is True for x in complete),
      'all_three_pass_count':sum(x['three_rule_reconstruction']=='1' for x in complete),
      'dubious_reason_rows':[{k:x[k] for k in ('date','ptu_filename_reference','FLCS_particle_fit','FLCS_diffusion_um2_per_s','derived_tau_ms','passes_published_particle_cutoff','passes_published_diffusion_time_range')} for x in byflag['-1']],
      'one_reported_standard_error_diagnostic_positive_cells':{
          'particle_margin_pass':sum(bool(x['one_reported_SE_particle_margin']) for x in positives),
          'CPM_margin_pass':sum(bool(x['one_reported_SE_CPM_margin']) for x in positives),
          'diffusion_margin_pass':sum(bool(x['one_reported_SE_diffusion_margin']) for x in positives),
          'all_three_margin_pass':sum(all(x[k] is True for k in ('one_reported_SE_particle_margin','one_reported_SE_CPM_margin','one_reported_SE_diffusion_margin')) for x in positives),
          'failed_diffusion_margin_files':[x['ptu_filename_reference'] for x in positives if x['one_reported_SE_diffusion_margin'] is not True],
          'meaning':'author-reported fit SE sensitivity only; NOT independent noise repeat or confidence statement'},
      'fitted_particle_number_positive_range':[min(x['FLCS_particle_fit'] for x in positives),max(x['FLCS_particle_fit'] for x in positives)],
      'E489_published_joint_counts':{f'{k[0]}_{k[1]}':v for k,v in sorted(counts.items())},
      'E489_empirical_observed_MI_bit_per_cell':observed_mi,
      'latent_state_nonidentifiability_even_if_no_false_positive':{
          'observed_sample_n':len(pairs),'model_1_Z_equal_detected_compatible_bit':mi1,
          'model_2_Z_constant_present_compatible_bit':mi2,
          'bound_status':'NOT_IDENTIFIED; counterexamples require no arbitrary detection parameter',
          'assumption':'Only if D=1 then Z=1; sensitivity and dependence on imaging state uncalibrated; illustrative constructions, not biological claims'},
      'physical_molecule_inheritance_MI':'NOT_IDENTIFIED','energy_efficiency_bit_per_J':'NOT_IDENTIFIED',
      'main_paper':'PUBLICATION_NO_GO'
    }
    (ROOT/'E490_results.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('source_rows','fit_complete','exact_agreement_count','cpm_pass_count','cpm_and_n_pass_count','all_three_pass_count','one_reported_standard_error_diagnostic_positive_cells','E489_empirical_observed_MI_bit_per_cell')},indent=2))

if __name__=='__main__': main()
