# 901 / E521 — Source-authenticated ASC assignment graph and exact MILP lower bounds

**Scientific scope:** Frozen original MoMA source candidate assignments and archived solver vectors. Not biological cell division ground truth, not full original Gurobi MPS optimization.

- [E521 formal report](https://drive.google.com/file/d/1dyXkONapPzxou-a8qqhIeudgDKSZhw2S/view)
- [E521 latest continuation card](https://drive.google.com/file/d/1kpT0XXSoKYSN3faleskUcvCWcQ7fkTNX/view)
- [E521 source/outputs/checksums folder](https://drive.google.com/drive/folders/1SGKfho8Sr7b_jZ1GWvpHbgVB0344IyKr)

## Reproducibility

Retrieve the already authenticated E512 full-original-component ZIPs for Pos5_GL35 and Pos12_GL22 into the same base directory, retaining their verified original filenames. Run

```
python research_901/E521_20261117/e521_source_graph_lowerbound_and_flow.py \
    --base /mnt/data --output /mnt/data/901_E521_20261117/results
```

Python 3 + SciPy 1.17.0 + NumPy. Program verifies two frozen E512 ZIP outer SHA256s, every ZIP-member CRC and the original source individual SOL/curation/forest member SHA and size receipts. It validates both archived main/cost solutions in the reconstructed subset of constraints. Both two-GL MILPs (with and without source-forest ancestry exclusion) achieve **HiGHS optimal gap 0**. Two full replays of seven derived files have matching SHA256.

### Verified results

- Original action-variable universes: Pos5 44,231, Pos12 47,879. Archived alternate ASC=1: 25/216; of these, main archived vector violates 8/17.
- Forcing those 25 maps into main and cancelling each old same-source outgoing action creates **24 lost-incoming-but-active-outgoing top-boundary daughter hypotheses** and **one double-incoming target**, not a valid lineage.
- Exact **minimum assignment-bit Hamming changes** from main, given ASC + per-hypothesis continuity/uniqueness: Pos5 29; Pos12 72. Add original source component-tree ancestor/descendant exclusivity: Pos5 29; Pos12 74. **103 lower-bound total** vs 219 actual assignment-bit differences in archived alternate solutions.
- The first Pos12 relaxation had two tree-ancestry conflicts, at f241 and f246. The stronger relaxation eliminates them.
- Original complete Gurobi MPS/LP bytes remain **unobtained** (external EBI DNS access blocked), so original full constraint matrices and source cost objective were **not** solved. Hamming is a custom methodological objective, not biological or physical energy. This partial formulation omits other original crossing, segmentation and auxiliary constraints.

**Scientific status:** The original twelve disputed divisions have no independent microscopic biological truth (UNRESOLVED), molecule bits/bits-s/bits-J unidentifiable and full manuscript PUBLICATION_NO_GO. Original E508–E520 evidence preserved. E464 16S branch unrelated. Synthetic MPS tests excluded. Full historical E1–E477 ZIP-member and source-paper/PaperQA review still pending.
