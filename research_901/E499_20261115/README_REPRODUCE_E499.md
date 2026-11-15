# Reproduce E499 exactly from its own archive

1. Extract this ZIP into a directory, preserving `source_E498/`.
2. `python analyze_e499_exposure_phase.py` (Python 3.11+, pandas and numpy)
3. `python independent_qa_e499.py` (Python 3.11+ standard library only)

Both scripts assert the official E498 source SHA for the original validated Gervais 2026 Feather, exact source-export gzip CSV, original author Rmd, E498 source frozen per-cell WLS and E498 old result JSON. If the source `p901_E498` is in the parent directory, it is used; otherwise both use this package's byte-identical `source_E498/` directory. The original `.feather` is included untouched; numeric calculations read the E498 previously row-preserving exported `.csv.gz`. This original-vs-export equivalence was audited in E498; E499 independent QA validates numerical results with an alternative streaming calculation, not an independent Arrow/Feather decoder.

Inputs: 91,838 real source cell-stage fluorescence observations, 1,659 previously source-verified per-cell WLS rates. Outputs 14 true grouped sampling schedules, 7 positive-selected ratios and 16 true same-cell pairs. Run identity and results were verified in independent directory; no older analysis is retrained.

Interpretation: model parameters and exposure schedules are observational variables; dose/quantum yield per flash are not source-estimated. No physical molecular copy number, mother-to-daughter information in bit, bit/s, bit/J. Main manuscript no-go pending physical calibration and full original paper plus PaperQA requirements. Sources:
- Gervais et al. (2026), https://doi.org/10.1038/s44320-026-00226-5
- https://zenodo.org/records/18059846 (official original member source; full 842MB zip NOT downloaded)
- frozen E498 cloud source, report and pass receipt in 901 Drive E498 folder.
