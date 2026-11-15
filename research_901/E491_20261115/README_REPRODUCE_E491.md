# 901 E491 — Reproduce and provenance

## Authoritative inputs (frozen)

1. Official Wu et al. 2025 bioRxiv source DOI **10.1101/2024.11.23.624995**; repository `WenbinWu2001/ODEinherit/data/cell_lineage_data.rda`, Git blob SHA-1 `7476698785cf75fde518ec0fd29dc59be38ae799`; raw byte SHA-256 `4e8c7b89e607f56d9e77c67f6479768f9883513c914f7da65330447f41306f11` (747,617 bytes). `metadata_original.csv` and `time_series_original_long.csv` were exported by an official-data-loading Rscript, original .rda retained. 85 cells, 25 roots, 60 daughters, 6 proteins. **48 raw acquisition points** transformed to **240 smoothed/interpolated grid points** (see upstream `man/cell_lineage_data.Rd`). No point-in-time claim.
2. Julou E487/E488 user Drive validated input: `E488_same_event_size_and_background.csv`, SHA-256 `7e4dadde0b754c78e21fa22ab03f9fdcca911a608cecefe7ab7b23ce0a527072`; `E488_cv_predictions.csv`, SHA-256 `b06ada09bbeb9b536305d165b08550d62c67a9a39fc7360d0409be2d4438c8a1`. Original 111,442,020 B Julou source already preserved in 8 parts in 901 Drive; E491 does not redownload/re-run E487/E488. Use 42,658 exact locked-event identities from E488.
3. Rosenfeld et al., *Biophys J.* 2006 DOI `10.1529/biophysj.105.073098`: **not** a raw data input. Its binomial partition/fluorescence calibration requires a genuinely pulse-chased non-synthesizing reporter and validated optical error. Public Swain Lab software example `datatree.dat` is **synthetic**, not the reported primary experiment.

## Technical correction to source receipt

The first cloud downloader's `original_source_receipt.json` contained an **incorrect `published_study_doi`**, `10.1080/01621459.2025.2525238`. Preserve that original audit record but **use `original_source_receipt_corrected.json`**, which confirms the bioRxiv DOI is `10.1101/2024.11.23.624995` (PubMed PMID `40631107`). No source bytes were changed.

## Run locally or in an authorized cloud execution environment

Requires Python 3 with `numpy` and `pandas`. Keep the official RData and exported source CSVs in the same folder.

```bash
python analyze_odeinherit_e491.py
python independent_qa_e491.py
```

For `analyze_fluorescence_conservation_e491.py` only, supply frozen E488 source files at `/mnt/data/p901_E488/` in the original controlled environment, or set the script's two source paths to local immutable copies **after matching the specified SHA256**. The canonical 42,658-event derived input `E491_frozen_division_sum_audit.csv` and `E491_division_sum_by_date.csv` are also archived and may be checked directly by independent QA without recomputing E488.

Repeatable tests:
- original RData SHA and Git blob are preserved; source numeric arrays are export products, not additional wet-lab observations;
- lineage identities and six protein labels are authenticated; 60 daughters, 12 excluded at earliest available birth-grid time, 48 have usable preceding grid point;
- 'wrong mother' comparisons use the **same leading cell-ID prefix only**, not an independently authenticated microscope field;
- the six expression scores are **per-marker fractions of alternative mother-candidate distances that exceed the actual mother distance**. Values are not Shannon MI;
- bootstrap draws entire root lineages of the source **within the same study**; quantiles do not convey independent laboratory coverage;
- source fluorescence ratio is not calibrated physical protein count; no common model is fit across E487 and Wu datasets.

## Current scientific status

`SOURCE_BYTES_PASS / SOURCE_IDENTITY_PASS / OBSERVED_FLUORESCENCE_SUM_AUDIT_PASS / INDEPENDENT_NUMERICAL_QA_PASS / TRUE_MOLECULE_COPY_TRANSFER_NOT_IDENTIFIED / FULL_LIFETIME_NOT_IDENTIFIED / PUBLICATION_NO_GO`.
