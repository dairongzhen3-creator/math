# E492 reproducibility / source guard

This project records **only a single published H4-Dendra2 mother/two-daughter source**. Never claim 145 raw time rows are 145 independent biological replicates.

Official published inputs:
- `PMC8165573_fullTextXML.xml` — Europe PMC source, protocol DOI `10.1016/j.xpro.2021.100557`, SHA256 `6a59cfb29e1e933577a40b99de88e5ee4ed7cdf1e6b0355e371d8a08fd9b6667`.
- `PMC7889982_fullTextXML.xml` — main iScience original JATS, DOI `10.1016/j.isci.2021.102075`, SHA256 `42d58f8db16d91fd09ed02ed3303c109d9bb8b84246c54e4d708bb21dab19903`.
- `PMC8165573_DataS3_mmc3_original.zip` — official Europe PMC returned **whole 26-member supplementary** archive, SHA256 `2c03b1356cd696dd7d06f0a820463f388942e08c702fe8aa5f32e02848e4c9b8`, 17,227,220 bytes. This bundle contains nested `mmc3.zip` (17,161 bytes), SHA256 `7dc3bce68688852585c11e2dff07351c874cc3cf118c7505ef5358fc674ebe47`.
- `PMC8165573_table_tbl1.csv` — Table1 JATS extraction, 145 authentic optical mean readouts: 61 mother, 50 daughter1, 34 daughter2. JATS is authoritative; exported CSV is derivative.

When reproducing in authorized cloud VM, supply the official XML, full 17MB ZIP, and CSV in same directory, then run:

```bash
python analyze_e492.py
python independent_qa_e492.py
```

Output: E492_results.json, E492_raw_synchronous_optical_means.csv, E492_independent_QA.json and extracted nested DataS3 script archive. Both primary and QA use only the Python standard library; no background calibration, statistical model, batch exchangeability, or independent microscopy repeat is assumed.

Time-relative comparisons: daughter1 from t=93.5 (inclusive) to t<229.5, daughter2 from t=229.5 to t=510. `F_M/(F_M+F_D)` is **raw per-cell mean fluorescence comparison**, not conserved molecules/old-protein segregation. Source DataS3 documents position-dependent photo-conversion delay (~40min first position). Missing published multi-mother BudJ CSV files prevent population molecular MI, bitrate, energy efficiency, and lifetime total.

Official workflow: `https://github.com/dairongzhen3-creator/math/actions/runs/43710397389`, DataS3 workflow `https://github.com/dairongzhen3-creator/math/actions/runs/43710547733`.
