# 901 E499 — latest authoritative cloud handoff | 2026-11-15

**Frozen prior stages**: E1–E498 real original sources, negative controls, original observations and calculations must not be overwritten or silently rerun. 901 main manuscript remains **PUBLICATION_NO_GO** pending physical calibration + mandatory complete full-text citation review/PaperQA.

## Read exact source artifacts first
- [E499 formal Chinese scientific report](https://drive.google.com/file/d/1Ek-2ReG9ZBnB0QoEfY3ZtpiGxZvNHlLD/view)
- [E499 current continuation card](https://drive.google.com/file/d/1m44i_Am8Gb2bzQcDng5ptcYR6tbABxQP/view)
- [E499 full **self-contained original-source and numerical replication** archive, 21 members](https://drive.google.com/file/d/1YjZoNzn_JTGf8wT5tJMlVD2Hp39-cEoT/view)
- [E499 independent process data and separate scripts](https://drive.google.com/drive/folders/162zcjjurDVTC4bfRIhStSUQzP9lzpDrr)
- [E499 verified 21-member ZIP Google Drive re-download SHA-256/CRC receipt](https://drive.google.com/file/d/1W6c3Dqu_lldPfNVwmyCDN2BGHmOH_myb/view)
- [E498 earlier certified original source archive and exact Feather](https://drive.google.com/drive/folders/1Zi7rfZow1EWrPYT5mHyF0givPvq1rVxV)

## Scientific result, source-verified

Gervais et al. 2026 Mol Sys Biol DOI 10.1038/s44320-026-00226-5, Zenodo version 18059846. Direct verified source Gervais `myframes_bleaching.feather` (9,171,952B, SHA `6e1733d6b641f5a47dca2851520fb9a0edfed2aad6faecbaf1ec5eece99917f7`), row-identity certified source `myframes_bleaching.feather.csv.gz` SHA `0a9b630b95af875ae6bba5be5bac5430acbf96cfcf4bc767576965e59f8251db`. All **91,838** real cell-phase source rows read; all unique ID+time. Author analysis Rmd SHA `7f915a34db6c3326a0be833bc1c04e57127d90684ce281c9fe4b0cd65a48de70`, source E498 per-cell WLS 1659 rates SHA `3e2a7c5a620d767ad9f0f2192f72008d52af90e0d0fd30d3eac188c0f0d0b45a`. Independent standard-library calculation separately reconstructed **all 1,659** source WLS values from original row-preserving gz and matches.

3 real source acquisition dates, **7 date×probe**, 14 group×phase. Original source native 3 minute growth: **18,205/18,331=99.31%** exact gaps. Starvation 12 minute: **70,705/71,436=98.98%** exact gaps. Positive-WLS source median 7 group starvation:growth hourly rate ratio median **0.128745** (min0.0829275 max0.183345); re-express per *planned* frame using source 4× interval ratio: median **0.5149803** (range0.331710–0.733381). This is an algebraic optical re-expression, **NOT intrinsic photobleaching per photon**.

Only **16** same physical source cell identities had at least six valid frames in both physiological phases (YFP9, GFP5, CFP2); 14/16 growth source hourly decline greater, 13/16 after scheduled-frame normalization; both phases positive in only 9. Biologically different phases, not independent culture repeats.

Actual design `[1, I_starvation, exposure_frequency]` rank **2 of 3** since the frequency is `20 - 15*I_starvation`. No `same physiological phase × two frame frequencies` factorial contrast. The source data cannot identify photo-bleaching coefficient separately from phase-specific dark/non-photochemical effects. Even with 1,659 measurements, design rank remains deficient. Published article says same confounded imaging schedule 3min growth/12min starve, source signal already degradation-corrected. No independent camera dose, actual fluorophores or protein molecule counts. **Physical transfer MI, bit/s, bit/J NOT_IDENTIFIED.**

No author email, no paid service, no submission; full provider 842MB original not downloaded; only exact source original members retrieved by verified HTTP 206/ZIP CRC and archived E498.

## E500 independent next step
1. Search official Gervais data/BioImage Archive S-BIAD2506 for **within-phase** multiple imaging frequencies or independent same-state photobleach technical exposure calibration; source-verify exact ID, true wavelength/intensity and noise. Do not assume it exists.
2. If no suitable source, formalize the minimal phase×exposure-rate cross-design with actual dose measurements and biological repeated units, or switch to other public data with same-state irradiation.
3. For molecular *transmission*, also need genuine same-division old/new molecule counts and optical calibration; never merge Gervais starvation times with independent Kscheschinski2026 10058-lineage sample or old E487 data to create synthetic bit values.
