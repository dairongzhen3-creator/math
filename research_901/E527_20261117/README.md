# E527 reproducibility — original MoMA source-cost decomposition and boundary pixels

Scope: **source-software cost accounting**, not actual biological division-time labels or any molecular information/energy measurement.

## Inputs

From user-owned previously verified 901 cloud archive, automatically mounted under `/mnt/data`:

- `901_E512_Pos5_GL35_full_original_components.zip` (SHA256 `77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6`)
- `901_E512_Pos12_GL22_full_original_components.zip` (SHA256 `73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2`)
- `901_E511_Pos5_GL35_original_moma_sidecars.zip` (SHA256 `1dad2aee7728757330fef6f16127faa5b5ab1648b78080db78a24621b2ffa596`)
- `901_E511_Pos12_GL22_original_moma_sidecars.zip` (SHA256 `be8695b97cb1a2c946a8f667b201a40574e064e8b1b203124a40cb4bca5adcba`)
- `901_E526_20261117/results/E526_151_local_moves_exact_integer_ranked.csv` (SHA256 `b04584801314d8d7a5b9b848ae8bec038a18aa8570e777348111af22464928bb`)
- `901_E526_20261117/results/centered_pos5_scale1/...changed_actions.csv` and Pos12 counterpart for separate 2-single-pixel source recheck.

## Commands

```bash
python e527_original_cost_mechanism_exact_decomposition.py --base /mnt/data --output /mnt/data/901_E527_20261117/results
python e527_independent_two_single_pixel_audit.py
python e527_plot_boundary_gap_vs_margin.py
```

The main script authenticates original ZIP outer SHA, every member CRC and required member SHA+size from E512 receipts, validates E511 actual `mm.properties`, checks exactly 755 changed source actions and solves a **symbolic source component telescoping identity** for all 151 existing E525 feasible one-frame moves. It does **not** train any model or run a new biological or original full Gurobi MPS optimizer. All source price strings and component costs are parsed as Python Decimal. Two independent scripts are used for source accounting and the two singleton pixel-gap cases.

Public-source code pinned to `df7a3dd92cd15344c82b19e83f96a1f132e52ab7`:

- [CostFactory.java](https://github.com/nimwegenLab/moma/blob/df7a3dd92cd15344c82b19e83f96a1f132e52ab7/src/main/java/com/jug/lp/costs/CostFactory.java): 340 pixel normalization, migration/growth functions, `getLogLikelihoodComponentCost` conditional formula.
- [GrowthlaneTrackingILP.java](https://github.com/nimwegenLab/moma/blob/df7a3dd92cd15344c82b19e83f96a1f132e52ab7/src/main/java/com/jug/lp/GrowthlaneTrackingILP.java): Map and Div component cost weights each 0.5, geometry compatibility.
- [FeatureFlags.java](https://github.com/nimwegenLab/moma/blob/df7a3dd92cd15344c82b19e83f96a1f132e52ab7/src/main/java/com/jug/development/featureflags/FeatureFlags.java): current pinned `UsingLogLikelihoodCost` mode, **unproven as historical 2022 compiled runtime flag**.
- [ComponentTreeUtils.java](https://github.com/nimwegenLab/moma/blob/df7a3dd92cd15344c82b19e83f96a1f132e52ab7/src/main/java/com/jug/util/ComponentTreeUtils.java): vertical size is maxY-minY.

## Limitations / falsification paths

The original probability-map TIFF slices for Pos5 frame960 and Pos12 frame118 have **not** been obtained; inferred single-gap pixel values `p≈0.51826` / `p≈0.55449` are **mathematical consequences IF the pinned log-likelihood component-cost mode was the historical algorithm**, not measured source probabilities. The actual source components and source assignment-cost differences ARE authenticated. No imaging perturbation or biological noise calibration was done. The max all-151 per-action public formula vs archived coefficient discrepancy is ~4.18e-6; for the near-optimal 55 it is ≤3.5e-8. Keep archived source values authoritative.

E527 first run had an assertion failure from incorrectly assigning early candidates' first source frame to f instead of f−1. Fixed, archived in `E527_preflight_frame_alignment_failure_kept.md`. The 151-case run then succeeded twice with identical output SHA. Independent two-singleton audit run twice and produced identical files.

E464/16S unrelated; original full MPS and real third-party Omnipose inference unresolved. 12 biological division disputes remain UNRESOLVED and molecule bit/bits-s/bit-J NOT_IDENTIFIED; publication NO_GO. Full E1–E477 old archive and citation+PaperQA audit remain open.
