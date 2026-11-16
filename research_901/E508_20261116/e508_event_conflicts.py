#!/usr/bin/env python3
"""901 E508: event-level audit of E507 main versus alternate MoMA source CSVs.

This source-level audit does NOT certify independent phase-contrast pixels, identify
a correct segmentation version, or infer molecular information. It reconstructs
the 12 standard-only events under a maximum-cardinality bipartite matching,
rather than assuming cross-version numeric cell IDs are comparable.
"""
import argparse
import collections
import csv
import hashlib
import json
import pathlib
import subprocess
import sys

ARCH_GLS = {
    "20221220_E08_GFP.zip": ["Pos28_GL37"],
    "20221220_E08_GFP_pl.zip": ["Pos12_GL22"],
    "20221220_rplN_GFP_pl.zip": ["Pos4_GL4", "Pos5_GL35"],
}
EXPECTED = {
    "Pos28_GL37": (9, 9, 9),
    "Pos12_GL22": (12, 9, 9),
    "Pos4_GL4": (8, 8, 8),
    "Pos5_GL35": (13, 4, 4),
}
ROOT = pathlib.Path("e508_source_conflicts")

def csv_rows(path):
    with open(path, encoding="utf-8-sig", newline="") as src:
        rows = [line for line in src if line.strip() and not line.startswith("#")]
    return list(csv.DictReader(rows))

def event_list(rows):
    by_cell = collections.defaultdict(list)
    for row in rows:
        by_cell[int(row["cell_ID"])].append(row)
    out = []
    for pid, cellrows in by_cell.items():
        for r in cellrows:
            if r["type_of_end"] != "div":
                continue
            frame = int(r["frame"])
            daughters = sorted({
                int(x["cell_ID"])
                for x in rows
                if int(x["parent_ID"]) == pid and int(x["frame"]) == frame + 1
            })
            if len(daughters) != 2:
                continue
            top, bottom = int(r["bbox_top_px"]), int(r["bbox_bottom_px"])
            out.append({
                "parent": pid, "frame": frame, "top": top, "bottom": bottom,
                "center_y": (top + bottom) / 2,
                "length_px": float(r["length_px"]), "daughters": daughters
            })
    return sorted((e for e in out if 200 <= e["frame"] <= 700),
                  key=lambda x: (x["frame"], x["center_y"], x["parent"]))

def compatible(a, b):
    return (a["frame"] == b["frame"]
            and abs(a["center_y"] - b["center_y"]) <= 5
            and abs(a["length_px"] - b["length_px"]) <= 8)

def maximum_matching(main, cost):
    """Maximum cardinality one-to-one matching: independent of greedy order."""
    choices = [
        sorted((j for j, c in enumerate(cost) if compatible(m, c)),
               key=lambda j: (
                   abs(m["center_y"] - cost[j]["center_y"])
                   + abs(m["length_px"] - cost[j]["length_px"]),
                   j))
        for m in main
    ]
    assigned = {}
    def augment(i, used):
        for j in choices[i]:
            if j in used:
                continue
            used.add(j)
            if j not in assigned or augment(assigned[j], used):
                assigned[j] = i
                return True
        return False
    for i in range(len(main)):
        augment(i, set())
    pairs = sorted((i, j) for j, i in assigned.items())
    assert len({i for i, j in pairs}) == len(pairs)
    assert len({j for i, j in pairs}) == len(pairs)
    return pairs

def approximate_cost_mother(main_ev, cost_rows):
    rows = [r for r in cost_rows if int(r["frame"]) == main_ev["frame"]]
    pool = []
    for r in rows:
        top, bottom = int(r["bbox_top_px"]), int(r["bbox_bottom_px"])
        center = (top + bottom) / 2
        length = float(r["length_px"])
        pid = int(r["cell_ID"])
        daughters = sorted({
            int(c["cell_ID"]) for c in cost_rows
            if int(c["frame"]) == main_ev["frame"] + 1
            and int(c["parent_ID"]) == pid
        })
        pool.append({
            "cost_cell_id": pid,
            "cost_end_label": r["type_of_end"],
            "cost_bbox_top": top, "cost_bbox_bottom": bottom,
            "cost_center_y": center,
            "cost_length_px": length,
            "cost_daughter_count_at_next_frame": len(daughters),
            "cost_daughter_ids": daughters,
            "center_delta_px": round(center - main_ev["center_y"], 3),
            "length_delta_px": round(length - main_ev["length_px"], 6),
        })
    pool.sort(key=lambda x: (
        abs(x["center_delta_px"]) + abs(x["length_delta_px"]),
        abs(x["center_delta_px"]), x["cost_cell_id"]
    ))
    nearby = [x for x in pool if abs(x["center_delta_px"]) <= 5
              and abs(x["length_delta_px"]) <= 8]
    if not nearby:
        reason = "NO_COST_MOTHER_WITHIN_FROZEN_GEOMETRY"
    elif nearby[0]["cost_end_label"] != "div":
        reason = "COST_NEARBY_MOTHER_NOT_LABELED_DIV"
    elif nearby[0]["cost_daughter_count_at_next_frame"] != 2:
        reason = "COST_DIV_MOTHER_WITHOUT_TWO_LINKED_DAUGHTERS"
    else:
        reason = "COST_HAS_NEARBY_DIV_SOURCE_NEEDS_PAIRING_AUDIT"
    return reason, nearby[:3], pool[:3]

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--archive", required=True, choices=sorted(ARCH_GLS))
    p.add_argument("--recompute-official-csv", action="store_true")
    args = p.parse_args()
    folder = pathlib.Path("e507_alternate_QA") / args.archive[:-4]
    if args.recompute_official_csv:
        subprocess.run([
            sys.executable,
            "research_901/E507_20261116/e507_alternate_segmentation.py",
            "--archive", args.archive
        ], check=True)
    outdir = ROOT / args.archive[:-4]
    outdir.mkdir(parents=True, exist_ok=True)
    records = []
    per_gl = []
    sources = []
    for gl in ARCH_GLS[args.archive]:
        mpath = folder / (gl + "_main_CellStats.csv")
        cpath = folder / (gl + "_cost_CellStats.csv")
        assert mpath.is_file() and cpath.is_file(), (mpath, cpath)
        mr, cr = csv_rows(mpath), csv_rows(cpath)
        main_e, cost_e = event_list(mr), event_list(cr)
        pairs = maximum_matching(main_e, cost_e)
        assert (len(main_e), len(cost_e), len(pairs)) == EXPECTED[gl], (gl, len(main_e), len(cost_e), len(pairs))
        covered = {i for i, j in pairs}
        cost_covered = {j for i, j in pairs}
        unmatched = [main_e[i] for i in range(len(main_e)) if i not in covered]
        assert len(cost_covered) == len(cost_e), (gl, len(cost_covered), len(cost_e))
        per_gl.append({
            "GL": gl, "main_divisions": len(main_e), "cost_divisions": len(cost_e),
            "maximum_cardinality_matches": len(pairs),
            "main_only": len(unmatched), "cost_only": len(cost_e) - len(cost_covered),
            "source_main_csv_sha256": hashlib.sha256(mpath.read_bytes()).hexdigest(),
            "source_cost_csv_sha256": hashlib.sha256(cpath.read_bytes()).hexdigest(),
        })
        sources.extend([
            {"GL": gl, "version": tag, "filename": str(path),
             "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
             "size": path.stat().st_size}
            for tag, path in [("main", mpath), ("cost", cpath)]
        ])
        for index, ev in enumerate(unmatched, 1):
            reason, exact_cost, closest = approximate_cost_mother(ev, cr)
            record = dict(ev)
            record.update({
                "archive": args.archive, "GL": gl,
                "source_version": "standard",
                "version_disagreement_not_image_ground_truth": True,
                "cost_source_reason": reason,
                "nearby_cost_cells_fixed_geometry": exact_cost,
                "nearest_three_cost_cells": closest,
                "event_id": f"{gl}_f{ev['frame']}_p{ev['parent']}",
            })
            records.append(record)
            print("E508_CONFLICT_EVENT " + json.dumps({
                "gl": gl, "frame": ev["frame"], "parent": ev["parent"],
                "bbox": [ev["top"], ev["bottom"]],
                "reason": reason,
                "nearest_cost_id": closest[0]["cost_cell_id"] if closest else None
            }, sort_keys=True), flush=True)
        print("E508_SOURCE_GL_PASS " + json.dumps(per_gl[-1], sort_keys=True), flush=True)
    for name, obj in [
        ("E508_event_level_version_conflicts.json", records),
        ("E508_gl_maximum_matching_summary.json", per_gl),
        ("E508_original_source_csv_sha256.json", sources)
    ]:
        (outdir / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")
    with (outdir / "E508_event_level_version_conflicts.csv").open("w", encoding="utf-8", newline="") as f:
        cols = ["event_id", "archive", "GL", "frame", "parent", "top", "bottom", "center_y",
                "length_px", "daughters", "cost_source_reason"]
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in records:
            w.writerow({k: json.dumps(r[k]) if k == "daughters" else r[k] for k in cols})
    print("E508_ARCHIVE_SOURCE_QA_PASS " + json.dumps({
        "archive": args.archive, "GL": len(per_gl),
        "unmatched_events": len(records),
        "provenance": "official original BIA ZIP64 source CSVs, verified by E507 unchanged code",
        "interpretation": "algorithmic version disagreement; no ground truth yet"
    }, sort_keys=True), flush=True)

if __name__ == "__main__":
    main()
