# 901 E519 — Archived constraint mechanism × original phase-image audit

**Source of truth:** [E519 latest handoff](https://drive.google.com/file/d/16Qx6qhkBE1tS0CK-KQOvy05wIdBSCE-r/view) and [formal stage report](https://drive.google.com/file/d/1efW8QuWRWwgILhnwIMhM2qQ8nxCkLs_n/view).

## Reproduce

Use the immutable, SHA-validated sources preserved in the earlier project cloud records:

- E518 two archived Gurobi/cost/ASC result CSVs in `901_E518_20261117/results/`
- E515 source-component/probability/5-frame-IoU CSV in `901_E515_20261117/results/`
- E517 original 203-frame **phase-only** y119→y120 gradient CSV in `901_E517_20261117/`
- E512 exact original component-forest ZIPs for Pos5_GL35 and Pos12_GL22 in the working root.

Run: `python research_901/E519_20261117/e519_exact_asc_source_vs_phase_stratification.py --base /mnt/data --output /mnt/data/901_E519_20261117/results`

The frozen script asserts original input SHA and original ZIP CRC, joins software frame `f` to true TIFF image frame `f+1`, checks 12 source events and 13 software-agreement technical controls, and authenticates the exception at Pos5 f247 using original source component-tree ancestry. It confirms f294 has no cost-ASC in the nearest 46 frames (closest 47). The script was actually run twice with six derived output SHA values identical in the ChatGPT isolated cloud runtime; there is no claimed GitHub Actions run here.

**Verified Drive deliverables:** [E519 archive, output tables, code and checksums](https://drive.google.com/drive/folders/1m3pXTiiwHWwmHirGBmtum_tYhJKM5oTR). Archive SHA256: `6c816d107a79329b5bbed42de12a4691712db8e32534b9952498c43e1da8e289`. The complete Drive readback hash and all eleven ZIP member CRC values passed.

**Scientific limitation:** All 12 biological division truths remain UNRESOLVED. The 10 direct ASC events are source-software mechanisms, not verified false divisions; optical ROI y=120 was selected retrospectively. True independent Omnipose inference and external expert ratings are not available. Physical molecular bits and bit/J not identified. Publication NO_GO. Other projects and the old E464/16S branch unchanged.
