# Reproduce 901 E490

## Source pins
- Original Julou (2020) PLOS Biology DOI 10.1371/journal.pbio.3000952.
- `FLCS_lacInduction.xlsx` original pinned publisher-author Git LFS sha256: `97f93074b72214bd8a20439ff508f6d8d47111ac115b6be6c06538457731212d` (33583 bytes).
- Author publication tag `plos`, source code blob `bfc059418ebcc88f326c9ba722abb1a6d5ef52fe` at `src/MoM_lacInduction_FLIM.Rmd`.
- `FLCS_lags.xlsx` original SHA256 `c927315fd7ce755e3f493a5d3ec152c3ed61a1ffbf09f722391523bcf559c7cf` (12317 bytes); historical E489 observed joint data held fixed.
- Original GitHub IDR0091 path index git blob `4799c7ded57d01c7f2446e465dbe3a446981890b` from `https://github.com/IDR/idr0091-julou-lacinduction/blob/master/experimentA/idr0091-experimentA-filePaths.tsv`. 1252 source paths, 36 historical experiment dates, no 2020 FLCS `.ptu` in this official inventory.

## Run
Use Python standard library for the primary CSV-derived re-analysis, `artifact_tool` for independent original XLSX checks.

```
python analyze_e490.py
python independent_qa_e490.py
```

Inputs include the two E489 source-frozen XLSX files and E489's original source-export CSV and exact join. Outputs: `E490_results.json`, `E490_author_rule_per_cell.csv`, `E490_independent_QA.json`.

No raw `.ptu` was acquired. 66/66 concordance is published-label algorithm reconstruction, not an independent optical measurement or biological experiment. One-SE tests are purely parametric sensitivity diagnostics; no false-positive/false-negative calibration is inferred. Refuse to infer physical molecule-transfer bit, bit/s, bit/J or full-lifetime total from these observations. No manuscript submission authorized.
