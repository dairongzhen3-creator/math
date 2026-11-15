# 901 E500 — genuine original-source phototoxicity/regrowth reproducibility

## Exact source

Published Gervais et al. 2026 *Molecular Systems Biology* DOI `10.1038/s44320-026-00226-5`, original Zenodo version `18059846`, exact ZIP-member `Gervais2026_data/Exported_data/myframes_control.feather` (302,295,912 bytes, CRC32 130ceec7, SHA256 `5cb1850cb127a64bf5e3dffa4d4e13e7cc648f73af176f6d03cf4d0e1fa9bb7c`). The complete provider 842MB ZIP was not downloaded. Original Feather was extracted and source-hashed in the frozen GitHub Actions source workflow and is included in this FULL archive.

Author methods `Gervais_et_al_analysis.Rmd` frozen SHA256 `7f915a34db6c3326a0be833bc1c04e57127d90684ce281c9fe4b0cd65a48de70` (source E498); use `phototoxicity` and `lag_phototoxicity` Rmd source blocks to validate model and source selection.

## Computational reproduction

Extract the full ZIP and run from its directory:

```
python analyze_e500_light_mother_fates.py
python independent_qa_e500.py
python independent_qa_e500_regrowth_lag.py
```

Main requires pandas numpy scipy statsmodels; both alternative QA scripts use pure Python standard library. Main and alternative numerical QA consume the 238-source-row CSV exported by the first authenticated Arrow reader; they both verify original Feather byte SHA. **Additional more independent raw decoding** was executed in separate GitHub Action: `E500_independent_Arrow_source_QA.json` certifies full 1,011,012 original Feather source observations and exact 238 identities with second Arrow RecordBatch code (not reused pandas analysis).

The source-vetted 71,755-row mother regrowth time series `.csv.gz` was extracted directly from the original Feather (no data censoring), and the 238-row lag regression CSV frozen from the original author's regression recipe. `independent_qa_e500_regrowth_lag.py` recomputes all 238 original fits from the source-exported data via independent Python standard-library OLS and verifies all 96 eligible lags and exclusion reasons. This remains the same experiment, not a new biological replicate.

`E500_results.json` has original 60h 2×2 and four reporter strain breakdown, 24 microscopy positions, 10,000 position-cluster bootstrap. `E500_mother_regrowth_lag_original_receipt.json` has the full cohort fit quality gate. Source study's published result precedes ours.

The official BioImage Archive `S-BIAD2506` metadata and `deepMoMA.json` (SHA256 `f564c31f5497a70d898dc4988582f086237d6eaeefe1fbd6272743618267826f`) and `legacyMoMA.json` (SHA256 `c3ca3b6903b6837e0860d891a85ba5e867b1224856ba2486d4112b2ba26b58ef`) are preserved; eight matching 20221220 image archives appear in their source registry. Raw 112.96GB images were not downloaded.

Only one acquisition day and noncrossed microscope positions. The NO_LIGHT treatment denotes no *fluorescence excitation* during starvation, not all light absent. Risk difference is observational regrowth endpoint, NOT experimentally calibrated photobleaching per frame. E499 remains rank deficient for chemical/optical-rate separation. No real molecular MI/bit/s/bit/J, no manuscript submission.
