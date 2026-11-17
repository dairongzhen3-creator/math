# E522 source-authenticated partial MILP — replay guide

From the existing E521 source-authenticated 2 GL bundles, run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python e522_asc_graph_source_cost_tiebreak_and_radius.py \
  --base /mnt/data --output /mnt/data/901_E522_20261117/results --time_limit 40
```

Required inputs (existing cloud files, never synthetic):

- `/mnt/data/901_E512_Pos5_GL35_full_original_components.zip`
- `/mnt/data/901_E512_Pos12_GL22_full_original_components.zip`

This script asserts the ZIP SHA-256 identity, ZIP CRC for all members, original internal SHA/size, full source assignment cost ledger coverage, original saved objective reconstruction, archived main/cost solution feasibility for the **partial** constraints, and binary primal and HiGHS dual optimality for two cost minimizations per GL.

The objective uses the **same archived `cost__assignment_costs.csv`** for both saved Gurobi vectors and the E522 vector. E522 constrains all alt archived ASC=1 plus selected original source flow/degrees, source component-tree ancestry exclusions, and (A) Hamming exactly minimal or (B) Hamming at most saved alternate's changes. It DOES NOT reproduce official complete MPS constraints or original full solver; no source MPS/LP bytes were available.

E521's seven output files are also emitted by the E522 extension as frozen ancestor-integrity diagnostics. The five E522-specific result files (4 changed-action CSVs and one JSON cost certificate) together make twelve files. Two complete runs yielded exactly the same SHA-256 for all twelve; see `E522_replay_source_manifest.json` and original terminal logs. Runtime stopwatch readings are intentionally excluded from hashed science outputs.

All 12 biological divisions remain **UNRESOLVED**, actual molecular info bits/bit/s/bit/J **NOT_IDENTIFIED**, manuscript **PUBLICATION_NO_GO**. E464/16S branch, all other projects and user local disks unchanged.
