# 901 E523 Reproduction — original software crossing and exit constraints

This is a **partially reconstructed source-code-defined MILP**, not the unavailable original Gurobi `.mps/.lp` model. Independent cellular truth is not established.

Required frozen already-authenticated inputs in `--base`:
- `901_E511_Pos5_GL35_original_moma_sidecars.zip` and Pos12 counterpart; the verified 2022 standard/alternate `mm.properties`, each with `USE_FEATURE_CROSSING_CONSTRAINTS=1`.
- `901_E512_Pos5_GL35_full_original_components.zip` and Pos12 counterpart; original `.sol`, `.json` forest, `.moma` ASC and cost `assignment_costs.csv`.
- E522 results containing `E522_<GL>_{tie,radius}_changed_actions.csv`, used to **audit** the candidate witnesses; E523's new optimization does not take those candidate choices as training labels.

Command for each `GL`:

```bash
python e523_source_crossing_exit_constraints.py --base /mnt/data --gl Pos5_GL35 --output results/Pos5_GL35 --run-milp --time-limit 90
python e523_source_crossing_exit_constraints.py --base /mnt/data --gl Pos12_GL22 --output results/Pos12_GL22 --run-milp --time-limit 90
```

Dependencies: Python, `numpy`, `scipy` (tested SciPy 1.17.0). Replays executed twice in the isolated cloud runtime; ten machine-readable per-GL results have identical SHA256. HiGHS results certified with `status=0`, `mip_gap=0` and consistency with original archived vectors.

Original-source mathematical rules:
- [GrowthlaneTrackingILP crossing terms](https://github.com/nimwegenLab/moma/blob/df7a3dd92cd15344c82b19e83f96a1f132e52ab7/src/main/java/com/jug/lp/GrowthlaneTrackingILP.java).
- [AdvancedComponent ordinal and below-leaf rank](https://github.com/nimwegenLab/moma/blob/df7a3dd92cd15344c82b19e83f96a1f132e52ab7/src/main/java/com/jug/util/componenttree/AdvancedComponent.java).
- [ExitAssignment exit ordering](https://github.com/nimwegenLab/moma/blob/df7a3dd92cd15344c82b19e83f96a1f132e52ab7/src/main/java/com/jug/lp/ExitAssignment.java), [LpUtils `getHup`](https://github.com/nimwegenLab/moma/blob/df7a3dd92cd15344c82b19e83f96a1f132e52ab7/src/main/java/com/jug/lp/LpUtils.java).

**Mathematical scope:** The exact original saved candidate actions and original saved objective coefficients are used. The reconstructed model additionally contains ASC=1, flow, in/out degree <=1, original forest ancestor exclusion, published crossing and exit inequalities. Its bounded source-cost optimum under Hamming caps equal to the **archived alternate model's action Hamming** (Pos5 35, Pos12 184) reproduces the archived alternate action vector. There are still original MPS/auxiliary constraints not available. Biological divisions 12/12 unresolved; no physical molecule bit/J estimate; manuscript publication NO_GO.
