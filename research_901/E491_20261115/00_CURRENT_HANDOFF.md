# 901 E491 — authoritative cloud handoff | 2026-11-15

**Read first:** [Project continuation](https://drive.google.com/file/d/1KRh_eGwdBNjt0V9Hu5MC1Pnm3mv7rOcg/view); [formal research report](https://drive.google.com/file/d/1DqJw7UFrVNhC7cwk_7C4pje5muTaXvPQ/view); [full independent-research replication archive](https://drive.google.com/file/d/1IN2wXsjyrLtVeQGDZuNo2gxQKRXI48dT/view); [original six-protein RData and receipts](https://drive.google.com/drive/folders/1yAUn63ZpfBX0T7o6H0doFyGXPAcRRD4J).

This stage **does not** overwrite E1–E490. No paper submission, no paid infrastructure, no author email. 901 research status is `PUBLICATION_NO_GO`.

## New grounded actions

- Retrieved Wu et al.'s publicly released `cell_lineage_data.rda` (747,617 original bytes; Git blob SHA-1 `7476698785cf75fde518ec0fd29dc59be38ae799`; source SHA-256 `4e8c7b89e607f56d9e77c67f6479768f9883513c914f7da65330447f41306f11`). Verified 25 original mothers, 60 daughters, six protein reporter channels. **Only 48 biological observation time points; the 240 processed points are smoothed/interpolated.**
- Among 48 descendants with a source timepoint preceding birth, 24 root-family groups, actual-parent-versus-other-candidate relative distances at the same source ID prefix and birth-time grid: Cdc10 `0.804315476`, CLB5 `0.425223214`; the source prefix is NOT independently verified as biological batch.
- From the already-validated E488 42,658 real bacterial mother/two-daughter events, the daughter-total-to-mother-preceding-frame **observed fluorescent signal** median is `1.0237207884`. `53.1600%` are in [0.95, 1.05]; `65.6313%` above 1. Cannot equate this with conserved protein molecule counts.
- Independent standard-library numerical checks passed for six protein identity metrics and all 42,658 observed fluorescence-sum ratios. Full scientific report discusses all source qualifications.
- Rosenfeld et al 2006 pulse-chase calibration DOI `10.1529/biophysj.105.073098` requires no synthesis/degradation in labeling window, validated response and noise. Its online software demo is synthetic, not the original wet-lab matrix.
- Wu et al paper is **bioRxiv** `10.1101/2024.11.23.624995`, not the incorrect DOI in the initial downloader receipt. Preserve the initial receipt and use `original_source_receipt_corrected.json` in the archived package.

## Reproduction

Source-pinned scripts in this folder:

- `analyze_odeinherit_e491.py` for authentic mother-vs-alternative candidate distances.
- `analyze_fluorescence_conservation_e491.py` for E488 **frozen** event fluorescence ratios, requiring SHA-pinned E488 features.
- `independent_qa_e491.py` for separate stdlib numerical QA.
- `README_REPRODUCE_E491.md` for all prerequisites and the exact source hashes.

The 23-member ZIP holds every small source, code, output, independent QA, published RData byte copy, and original and corrected receipts. **Do not** regard this newly accessed yeast dataset as an independent optical calibration, a measured molecular transfer, or a cross-species rank.

## E492 shortest next step

Find truly measured old-protein/new-protein pulse-separated mother/daughter data with same-division IDs, independent optical noise or technical repeats and absolute photon/molecule calibration. First verify the actual TrIPP Auboiron et al 2021 raw quantitative mother/daughter data (iScience DOI `10.1016/j.isci.2021.102075`); then Rosenfeld2006 **real** experimental lineage numbers rather than its simulated `datatree.dat`. If unavailable, state the precise missing payload; do not repeat E487/E490 calculations or invent physical bit/J.
