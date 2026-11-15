# 901 E497 — reproducibility, scope and source pins

Source authors' original archives stored unmodified as 8 Google Drive chunks in E496 source folder, original and reassembled SHA256:

- `mother_machine_data_sets.zip`: `a5787b437d9d136fd308feab5068fd558790f81929265444b42451287bc7b421` (52,912,016 B).
- `replicates_hi1_rrnB.zip`: `e6378d36486757c807f03a644e9653d131c24b00ba279472b2e512a0d43a1b6f` (16,405,453 B).
- `E496_2026_original_verified_sibling_events.csv`: SHA256 `c772970c39ea34bd677f813ac6d0c2d1d636813346777ab738ab62281e26168e` (10,058 canonical source IDs).

The following two executable scripts require the source directories `p901_E497` adjacent to `p901_E496` exactly as in this project. To reproduce elsewhere adjust `SRC` in both scripts to the path of immutable verified sources only. Both need original author ZIPs, no `npz` object/pickle deserialization.

```bash
python analyze_e497_realtrace_measurability.py
python independent_qa_e497.py
```

The main script uses numpy/pandas to reverify all 31×11 source parameter values against NPY scalar fields and separately derive raw input versus posterior inferred sibling fractions and LODO. Independent QA does not use numpy/pandas, parses raw author ZIP and float64 NPY files with stdlib `csv`,`zipfile`,`struct`,`ast.literal_eval`, independently checks 682,679 rows and 30,174 exact source keys, then recomputes all estimands via sufficient statistics.

Authoritative exact model source blobs: `src/likelihood.h` 4d6bb32d7df405c247d1878b23eac67195a0380f; `src/predictions.h` 3afc2aa20bf6916761b0be1d18defad7d09daeba; `RealTrace_Tutorial/RealTrace_toolbox/setup.py` f2e102e3dbbdde10ee22888643e038f63a606eb8 from `nimwegenLab/RealTrace`.

`fp` in official file is a processed autofluorescence-corrected optical input and `mean_g`, `cov_gg` are latent-state posteriors inferred using the same data. `var_g` has two likelihood noise-model options; a source file's free fit scalar does not uniquely fix source camera photon variance. The posterior is not independent optical technical repeat. No molecular copy transfer bit, bit/s or bit/J is inferred. No paper submission.
