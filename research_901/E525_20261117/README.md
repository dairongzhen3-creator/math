# 901 E525 reproduction — source ILP local time-shift envelope and E524 correction

**Not original full MPS, not independent biological truth.** Inputs are already archived and SHA/CRC verified E511, E512, E510 ZIPs under `/mnt/data`; no external model weights or network access are used.

## Commands

```bash
# 1. Each GL two identical source-model replay runs:
python e525_five_action_division_time_shift_scan.py --base /mnt/data --gl Pos5_GL35 --run-scan --output results/Pos5_GL35
python e525_five_action_division_time_shift_scan.py --base /mnt/data --gl Pos12_GL22 --run-scan --output results/Pos12_GL22

# 2. True E513 zero-based phase signal only on existing 37 technical cases:
python e525_source_phase_observability_gate.py --base /mnt/data --output results/phase

# 3. Source CSV independent exact decimal counterexample audit:
python e525_independent_source_decimal_counterexample.py --base /mnt/data --e525 /mnt/data/901_E525_20261117 --output results/decimal

# 4. Previously existing E510 software-agreement technical controls:
python e525_existing_technical_controls_join.py --base /mnt/data --e525 /mnt/data/901_E525_20261117 --output results/technical_controls
```

The source scan constructs only five-bit local time swaps of original selected Div actions and verifies every feasible proposal against the E524 source-derived ASC/flow/terminal/source tree/Crossing/Exit sparse constraint matrix (both sparse row delta and full matrix product). This is a **feasibility survey**, not a fresh HiGHS optimization or full-MPS result.

Main results: 98 source optimizer divisions, 497 syntactically extant five-bit shifts, 151 partial-model feasible changes; 96/98 source Div events have at least one. All 151 cost increments positive, 55 increments ≤0.001 in original software cost units. **E524 numerical correction**: original Decimal CSV Pos12 f118→117 cost delta `0.00001565` is smaller than previously reported f45→44 `0.00001587`; the earlier `mip_gap=0` certificate cannot be claimed as exact global second-best rank at 2.2e-7 precision.

Scientific limitations: selected solver Div includes pruned/non-exported hypotheses. E510 cache covers only 203 unique observed zero-based frames and does not include image f44–45, f117–118 or f959–960. All biological division truth remains UNRESOLVED; absolute physical molecular bits not identified; manuscript publication NO_GO. Early historical ZIP-member and original reference-paper/PaperQA review remains incomplete.

Official-source MoMA code fixed at GitHub commit `df7a3dd92cd15344c82b19e83f96a1f132e52ab7`. The full original `.mps/.lp` remains unacquired. Frozen E508–E524 records and prior failures must not be overwritten.
