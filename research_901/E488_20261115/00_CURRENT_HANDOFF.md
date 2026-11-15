# 901 E488 — 2026-11-15 current continuation

Research boundary: quantifying real biological internal/external information. This stage continues verified E487 Julou 2020 Zenodo source and does not replace prior evidence.

## New independently audited evidence

- Published original source: Julou et al. 2020 PLOS Biology DOI 10.1371/journal.pbio.3000952; Zenodo 10.5281/zenodo.3894719. Original 111,442,020 B / SHA256 `4b63dfee10ad9839060186671f888178e10d3962a80481c752175fa37d23ccd6`; retained in E487 Google Drive independent 901 data warehouse.
- E488 new field extraction: 42,658 physical mother/two-daughter observations, parent t-1 and t reporter values, daughters t+1, across 27 QC-included dates and 857 microscopy lanes; date-specific, held-out whole lane fivefold model. Same event and fluorescence bytes cross-checked with E487.
- Frozen 4-bin per-date target learned ONLY from each fold's training lanes; covariates standardized on training lanes only. Source baseline includes maternal lengths, previous length, imaging background, lane cell count, date and frame index.
- Heldout log loss in bits/observed division: BASE 1.336575136; PREV 0.663479006; LAST 0.625563535; BOTH 0.602925082. Additional recent maternal signal beyond PREV: +0.060553924 bit/heldout division (27/27 dates, 5/5 folds positive).
- Excluding source dates 20190605 and 20190614 with documented unusual illumination changes and **refitting models**: 25 dates, 40,628 divisions, 794 lanes, additional recent-state gain +0.055730733 bit/observation; all 25 dates and five folds positive.
- Independent standard-library audit verified all identity keys, lane fold boundaries, per-row paired score differences, and totals; software QA PASS.
- Original author's tagged `plos` R source uses `length_um=(vertical_bottom-vertical_top)*0.065`, an autofluorescence baseline `422.8*length_um`, and a conditional proxy calibration. Published paper also includes separate single-molecule FLCS experiments, **not directly joined with these 42,658 division IDs**.
- **Physical molecular transfer `I(Smother;Sdaughter)`, bit/s, bit/J and lifetime total remain NOT_IDENTIFIED.** Reporter log-score gain is neither causal physical transfer nor Shannon channel MI. Main paper PUBLICATION_NO_GO; citation full-text + PaperQA remains unfinished.

## Immutable cloud artifacts

- [Full E488 research report](https://drive.google.com/file/d/1BAJsyQEgB6a3wrYZqChyc1AJSyuHCRba/view)
- [Current E488 continuation](https://drive.google.com/file/d/1uHUEyAVh-l9xJB32jPj497aBd4TLrU-x/view)
- [Source inputs, code, model scores, independent QA, execution manifest and complete reproducible archive](https://drive.google.com/file/d/1Oz6AqFsuhL6OnsvuhdH302RSrvPpdoFg/view)
- [Complete E488 process data folder](https://drive.google.com/drive/folders/1wgoSS1zVduE2_edstSAtf3m_d2uQc9KK)
- [E487 verified immutable 111MB official source chunks](https://drive.google.com/drive/folders/1TmXJZJzxgLHTrSv2uLPqSoapCp-O2jyH)

Next E489: test whether author's FLCS/calibration original identifiers link to mother/daughter same physical division and same-time repeat exposures. Do not invent a noise model, causal claim, absolute molecular transfer, energy efficiency or submit the paper.
