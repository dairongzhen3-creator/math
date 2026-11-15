# 901 E494 full reproducibility / original source and methodological boundary

1. Frozen exact E488 original source inputs must be together with two scripts:
  - `E488_same_event_size_and_background.csv`, SHA256 `7e4dadde0b754c78e21fa22ab03f9fdcca911a608cecefe7ab7b23ce0a527072`
  - `E488_cv_predictions.csv`, SHA256 `b06ada09bbeb9b536305d165b08550d62c67a9a39fc7360d0409be2d4438c8a1`
  - **Do not rerun E488's original model**. E494 identifies the fixed exact 42658 event keys by inner join and creates a *new* daughter-size-vs-fluorescence audit only.
2. `python analyze_e494_size_partition.py`: needs numpy/pandas; fixed seeds. Outputs per-event observations, per-date summaries, 27-date heldout optical MSE and log scores, 199 within-lane shuffles, date-block descriptive resampling.
3. `python independent_qa_e494.py`: different standard library CSV direct parser and independent sufficient-statistics equations; independently validates 42,658 original keys and all numeric inference; no pandas, numpy or functions from main analysis. Both programs read same source, so this is **computational QA**, not new biological experiment.
4. Original 2018 Nature Communications paper JATS `PMC5768764_2018_fullTextXML.xml` SHA256 `a3e5574d61551fad47044d61200455bac079847785f0927a24a3da4d3e8869d8`. Zenodo official `zenodo_746230_ReadMe.md` SHA256 `4aa76d90746c50dbbdf11812962a38c55586428103bd97619f91e9c1162dbf30`, MD5 `48dd9429111f1bb4ff7cc28b38f86759`. Zenodo concept record resolves to version `824793`; published processed dataset **4.43 GB not downloaded**, source images **34.49 GB not downloaded**. The **357** sibling calibration factor is from published paper, *not independently re-estimated here*.
5. Secondary 2015 localizome official JATS `PMC4309519_localizome_fullTextXML.xml` SHA256 `42b5f5b5d4ca1ef1c7028d5ea8e3ac2dff8790a0ac95427280cf0046af2b9afc`; no original 773-protein per-division full optical matrix downloaded. Its article and 2018 article are external context, not dataset merges.
6. Model scope: pL is **estimated daughter pixel length fraction**, NOT true cytoplasmic volume, and pF is **processed fluorescent amplitude fraction**, NOT molecules. Heldout optical score and size r cannot be called causal protein transmission, physical Shannon channel capacity, bit/s, or bit/J. No paper submission.

Source files and E494 new script outputs are in E494 ZIP. Independent report lists source identity, missing data and exact claims.
