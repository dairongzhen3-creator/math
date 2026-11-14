#!/usr/bin/env python3
"""E486: source-frozen algebraic nonidentifiability, not molecular transfer estimation.
Reads only E485's audited 2x2 table. No optimized or arbitrary measurement noise rate.
"""
import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "observed_E485_2x2.csv"


def binary_entropy(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def main():
    with SOURCE.open(newline="", encoding="utf-8") as f:
        data = {r["mother_pole_class"]: r for r in csv.DictReader(f)}
    assert set(data) == {"new", "old"}
    counts = {k: (int(v["new_bright_yes"]), int(v["new_bright_no"]))
              for k, v in data.items()}
    assert counts == {"new": (119, 58), "old": (147, 31)}
    ns = {k: sum(v) for k, v in counts.items()}
    assert ns == {"new": 177, "old": 178}
    n = sum(ns.values())
    q = {k: counts[k][0] / ns[k] for k in ns}
    w = {k: ns[k] / n for k in ns}
    assert n == 355

    def info(epsilon):
        qp = {k: (q[k] - epsilon) / (1 - 2 * epsilon) for k in q}
        assert all(-1e-12 <= v <= 1 + 1e-12 for v in qp.values())
        avg = sum(w[k] * qp[k] for k in qp)
        return binary_entropy(avg) - sum(
            w[k] * binary_entropy(qp[k]) for k in qp
        ), qp

    # Binary symmetric independent observation error is ONLY a conditional,
    # unverified sensitivity model. Feasibility follows latent q in [0,1].
    epsilon_max = min(*(min(val, 1 - val) for val in q.values()))
    observed_mi, q_at_zero = info(0.0)
    max_model_mi, q_at_boundary = info(epsilon_max)
    assert abs(observed_mi - 0.02290337097781335) < 1e-12
    assert abs(epsilon_max - 31 / 178) < 1e-12
    assert max_model_mi >= observed_mi
    # Both distinct models produce precisely the same observed probabilities.
    for epsilon, qp in ((0.0, q_at_zero),
                        (epsilon_max, q_at_boundary)):
        for cls in ("new", "old"):
            observed = epsilon + (1 - 2 * epsilon) * qp[cls]
            assert abs(observed - q[cls]) < 1e-12

    result = {
        "stage": "E486", "date": "2026-11-14",
        "input_type": "E485_FROZEN_DERIVED_2x2",
        "source_355row_csv_sha256": "3fd9763b9b50135b5b3b2a97f5ffcf0936c87fabb2860d033998f0fb36c3a366",
        "source_official_supp5_xlsx_sha256": "d9b5e2f74007d6ffe298bc94c39e2f554b17655b48565a5de96a1f3eb8595baf",
        "n": n,
        "counts": counts,
        "empirical_observed_mi_bit_per_sibling_pair": observed_mi,
        "model": "common_age_independent_symmetric_binary_observation_flip_UNVERIFIED",
        "epsilon_feasible_max": epsilon_max,
        "model_latent_mi_lower_at_epsilon0": observed_mi,
        "model_latent_mi_upper_at_epsilon_max": max_model_mi,
        "upper_model_conditional_new_bright_by_age": q_at_boundary,
        "identifiability": "NOT_POINT_IDENTIFIED",
        "real_molecule_transfer_bit": "NOT_IDENTIFIED",
        "physical_causal_conclusion": False,
        "independent_biological_replicates_certified": False,
        "no_reanalysis_E485": True,
    }
    print(json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
