# 901 E493 reproducibility — source-frozen 2026-11-15

Source 1: Auboiron et al. (2021) *iScience* DOI `10.1016/j.isci.2021.102075`, published supplementary Table S2 `mmc2.xlsx` via official Europe PMC archive, original `PMC7889982_TableS2_SILAC_original_mmc2.xlsx` (842092 bytes), JATS-published MD5 `4d2104fd2a40a51bf21b8dddfc0f3428`, SHA256 `2ea57bedff6ccdd71ea982e695eea638cb6a815a1d2c8347b7ae91d8d8c21893`. Cloud acquisition and exact source receipt https://github.com/dairongzhen3-creator/math/actions/runs/43712125191.

Source 2: Auboiron et al. *STAR Protocols* (2021) DOI `10.1016/j.xpro.2021.100557`, exact original `PMC8165573_fullTextXML.xml`, SHA256 `6a59cfb29e1e933577a40b99de88e5ee4ed7cdf1e6b0355e371d8a08fd9b6667`. Historical E492 source frozen; this bundle includes the original research input for independent verification.

Run from the folder containing these two original bytes and scripts:

```bash
python analyze_e493.py
python independent_qa_e493.py
```

Main uses `artifact_tool` workbook reader + standard library XML. Independent QA uses **standard-library OOXML ZIP/XML reader instead of artifact_tool**, SciPy rank correlation, and Decimal48 independent original JATS parsing. The source decoder and written CSV are not shared between the two analyses. QA checks all 1338 raw half-life identities and both optical source windows. Python environment uses free libraries only; no paid compute.

Outputs `E493_results.json`, `E493_SILAC_source_locked_1338.csv` (original peptide IDs, raw seven-time H/L ratios and fitted rate columns), `E493_TrIPP_same_frame_optical_double_ratio.csv`, `E493_independent_QA.json`, alongside formal Chinese scientific report and continuation.

Dataset distinctions: fitted half-life is NOT measured mother-daughter molecule copying. Three negative fitted rates remain preserved. Source TrIPP Table1 provides only one physical lineage and uncorrected per-cell mean fluorescence. Single-lineage 84 matched frame pairs do NOT provide 84 independent divisions. `I(molecular mother; daughter)` and bit/s and bit/J remain NOT_IDENTIFIED. No paper submission.
