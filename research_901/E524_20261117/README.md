# 901 E524 — Uncapped source-cost minimization and best-distinct optimum

This directory preserves **actual original-source-verified partial MoMA optimization experiments**, not the unavailable original complete `.mps/.lp` model and not physical biological cell-division labels.

Required existing frozen inputs under `/mnt/data`:

- `901_E511_Pos5_GL35_original_moma_sidecars.zip`, `901_E511_Pos12_GL22_original_moma_sidecars.zip`
- `901_E512_Pos5_GL35_full_original_components.zip`, `901_E512_Pos12_GL22_full_original_components.zip`
- The four E522 source-candidate CSV files under `901_E522_20261117/results/` for independently checking previous partial witnesses.
- `901_E523_20261117/901_E523_交叉与退出约束_双生长线优化闭合与复算全套_20261117.zip` only for validating earlier terminal Exit states in the Decimal independent certificate.

## Commands

```
python e524_unbounded_with_terminal_flow.py --base /mnt/data --gl Pos5_GL35 --output results/terminal_pos5 --run-milp --time-limit 90
python e524_unbounded_with_terminal_flow.py --base /mnt/data --gl Pos12_GL22 --output results/terminal_pos12 --run-milp --time-limit 90
python e524_best_distinct_solution_probe.py --base /mnt/data --gl Pos5_GL35 --output results/nogood_pos5 --run-milp --time-limit 75
python e524_best_distinct_solution_probe.py --base /mnt/data --gl Pos12_GL22 --output results/nogood_pos12 --run-milp --time-limit 75
python e524_independent_decimal_vector_certificate.py
```

SciPy 1.17.0 with its HiGHS backend was used for actual MILP solutions. In resource-constrained cloud execution, run jobs **sequentially** rather than all four at once, and increase `--time-limit` for retries without changing source/model/optimizer coefficients. Each main experiment has two **successful** optimal solves, and the source assignment CSV SHA256s match on the two successful executions. The extra 4-way parallel stress test had three timeout statuses, retained in `logs/failed_parallel_replay_*.log`, not counted as successful.

**Scope:** Source action variables and per-action cost as recorded in original source ZIP; archived ASC=1, source-forest ancestry exclusion, in/out degree and continuity including true final t=1035, crossing and exit constraints derived from official MoMA sources, not full original historical solver MPS. No Hamming cap in the E524 main objective. `no-good` experiments additionally exclude the single archived alternate action vector, which is a *diagnostic*, not a claim that the author used this row in the original optimization.

The uncorrected `probe_*` results are stored as transparent history: E523 inherited an omission of the terminal continuity constraints. The corrected script forces final-frame continuity. The original source files are neither edited nor replaced; all earlier E508–E523 outputs frozen.

- Source code: https://github.com/nimwegenLab/moma/blob/df7a3dd92cd15344c82b19e83f96a1f132e52ab7/src/main/java/com/jug/lp/GrowthlaneTrackingILP.java
- Detailed scientific report: `901_E524_正式科研阶段报告_20261117.md`
- Study outcome: software-mathematical mechanism and small one-frame second-best timing margins, **not** actual organism division timing; 12 biological cases UNRESOLVED, bits/s and bits/J NOT_IDENTIFIED, main publication NO_GO.
