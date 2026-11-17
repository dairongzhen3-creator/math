# 901 E520 — Original solver pruning / exported-lineage audit, and native MPS acquisition gate

## Live scientific handoff
- [E520 formal stage report](https://drive.google.com/file/d/1QLCcTTuz0PmIu8m4GRg47DPH09rwkDV3/view)
- [E520 latest handoff card](https://drive.google.com/file/d/1msgzxktiMYWOM6-m3rpejbHQTmvVLamk/view)
- [E520 authenticated source/process archive](https://drive.google.com/drive/folders/1MHW7kpVu3DDjEaVdnfVSaUb0_Xrv24jf)

## What was really run
Two distinct source-audit scripts were executed twice each in the ChatGPT isolated cloud environment, with 11 data outputs matching SHA256 across runs. Inputs were prior frozen E508/E511/E512/E512b official-source-authenticated cloud bundles. Drive upload **and remote raw-file readback SHA256 / ZIP21 members CRC** verified the finished archive with SHA256 `21108c179aca7269c44cfacbae43bbfce79d53d3a29e3187a5a44c1b9e669df4`.

1. `e520_asc_relaxation_exact_certificate.py` proves that of 241 saved alternate ASC=1, the original standard solution violates 25; removing only the ten directly tied to disputed mothers leaves 15 additional violations. This is an **ASC-subsystem** compatibility test, not a complete ILP feasibility/optimality result.
2. `e520_pruned_solver_vs_exported_lineage.py` separates selected optimization `Div` from actual qualified mother+2 daughter CellStats export. Pos12_GL22 f200–700 standard selected Div=26 but exported divisions=12; 14 extra selected Div are in the subtree descended from f261 `PR` (pruning root). Alternate has 9 selected Div, 9 qualified export divisions and a run of 361 selected Map f246–606, including 160 consecutive ASC-fixed Map f246–405. Main CellStats corresponding track cell_ID79 f247–261 ends pruned; alternate track cell_ID76 f243–607 ends exit.
3. `e520_zip64_targeted_gurobi_mps_fetch.py` is strict-HTTP206/source manifest/ZIP64 member local header/CRC/SHA MPS+LP byte retrieval, and `e520_zip64_targeted_gurobi_mps_fetch_cloud_memorysafe.py` limits model-file magic-header checks to the first 8192 bytes. Both synthetic fixture tests pass; **the original 8 official MPS/LP members have NOT been fetched**, owing to EBI network/DNS failure in this runtime. No real MPS SHA, no independent HiGHS solve, no ASC-removal full counterfactual. Synthetic fixture data is excluded from the uploaded official evidence archive.
4. Submitted workflow `.github/workflows/901-e520-fetch-official-mps.yml` offers future authenticated extraction in GitHub Actions. There is **no confirmed E520 GitHub Actions run**; do not claim execution from workflow submission.

## Evidentiary scope
The observations are **original-source software and optimization export behavior**, not biological cell-division ground truth. All 12 original biological disputes UNRESOLVED. Absolute molecule copies, bits, bit/s and bit/J NOT_IDENTIFIED; publication NO_GO. E464/16S old branch, all other paper IDs, user local storage, email and paid services unchanged. Historical E1–E477 individual member review and citation original-paper/PaperQA review remain open.
