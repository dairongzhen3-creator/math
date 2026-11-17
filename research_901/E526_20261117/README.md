# 901 E526 reproducibility — exact archived-cost lattice and centered no-good MILP

This continues E525 and corrects overstrong E524 Pos12 “strict second-best” wording. **Not** the unavailable original full MPS nor an independent biological cell-division truth assessment.

## Required original authenticated materials

Already mounted in `/mnt/data`:

- `901_E512_Pos5_GL35_full_original_components.zip` (SHA256 `77fd2189e68c2fe07578bb81634d88c626e4e199895232ce194e3569c1af13a6`)
- `901_E512_Pos12_GL22_full_original_components.zip` (SHA256 `73d73c7a74f9df8cf4e3a83d51b171e14f01dbcbc7494a2cd3aaab3d457b60b2`)
- original E511 sidecar ZIPs (source `mm.properties` original SHA as E524)
- E524 full evidence ZIP, E525 full evidence ZIP (both zip CRC passed)
- E522 previously stored `results/E522_<GL>_{tie,radius}_changed_actions.csv`, unpacked from its verified original evidence ZIP, as needed for the inherited E524 source audit script (these are context diagnostics, **not** labels or optimization constraints)

## Commands

```bash
python e526_exact_decimal_local_family_certificate.py
python e526_independent_fraction_ranking_crosscheck.py
python e526_centered_flips_best_distinct_probe.py --gl Pos5_GL35  --base /mnt/data --output results/centered_pos5_scale1 --run-milp --time-limit 105 --cost-scale 1
python e526_centered_flips_best_distinct_probe.py --gl Pos5_GL35  --base /mnt/data --output results/centered_pos5_scale1e5 --run-milp --time-limit 105 --cost-scale 100000
python e526_centered_flips_best_distinct_probe.py --gl Pos12_GL22 --base /mnt/data --output results/centered_pos12_scale1 --run-milp --time-limit 105 --cost-scale 1
python e526_centered_flips_best_distinct_probe.py --gl Pos12_GL22 --base /mnt/data --output results/centered_pos12_scale1e5 --run-milp --time-limit 105 --cost-scale 100000
python e526_independent_fraction_ranking_crosscheck.py
```

Use sequential single-threaded CPU for scarce cloud resources. Original E524 code is retained read-only, and the new `e526_centered_flips_best_distinct_probe.py` contains the complete inherited actual published-source matrix construction **plus a fully algebraically equivalent flipped-variable representation**, with no need to import an inaccessible runtime. Source archive/member hash verifications are still active.

- **Four success runs** (`status=0`, `mip_gap=0`, each complete original-row feasibility replay, no guessed results). Each growthlane's two output 5-action CSV files have identical SHA256.
- **151 local feasible 5-action swaps** independently audited by exact Decimal integer-lattice and Python Fraction against raw original cost strings.
- **Numerical interpretation:** Pos12 centered no-good solves both pick 118→117, exact cost +0.00001565, whereas old uncentered E524 no-good had 45→44, exact +0.00001587. That discrepancy is evidence for scale/conditioning sensitivity; neither decimal arithmetic on the candidate alone nor floating MIP gap=0 proves a rational global MILP optimum.
- **Evidence scope:** Only source-code-reconstructed MoMA ASC, candidate-tree, continuity (incl. final f1035), crossing and exit conditions. Full historical MPS/auxiliary constraints and qualified microscopy labels missing. E464/16S unchanged; 12 biological disputes remain UNRESOLVED, physical bit/J NOT_IDENTIFIED, publication NO_GO.
