# 901 E510 | Frozen blinded 2D phase-contrast review protocol (pre-image interpretation)
Date: 2026-11-16. Scope: 901 only; existing E478–E509 results remain immutable.

## Source and frozen inputs
- Official BioImage Archive S-BIAD2506, full **preprocessed** GL TIFF, not raw camera FOV.
- GLs: Pos5_GL35 from 20221220_rplN_GFP_pl.zip; Pos12_GL22 from 20221220_E08_GFP_pl.zip.
- Frozen E509 source diagnostic defines 12 disputed standard-only division events, with 13 additional same-GL dual-version matched division events serving as *software-consensus positive controls*, not biological truth.
- Add exactly 12 *software-nondivision* descriptive time controls, chosen deterministically as the first available frame in event+10...event+60 whose distance from each E509 reported division is >6 source frames and whose frame has not been selected already; such windows may contain unrecognized divisions and are NOT certified true negatives.
- Source frame conversion: original TIFF position = source MoMA frame number - 1 (established E506 for analogous official source); source confirmation per GL is required if independent temporal provenance conflicts.
- Each case: original 2D phase channel frames at offsets [-3,-2,-1,0,+1,+2,+3,+5,+8], full channel width 106 pixels, original complete vertical dimension in exact float32 NPZ. Visualization always shows the fixed y=80..254 band, nearest-neighbor 3x enlargement, one display scale determined by full-stack temporal/spatial sample. No MoMA bbox, lineage ID, division class, or local optimization over cases appears on blinded cards.
- Preserve source TIFF SHA256, ZIP64 central SHA, original member CRC and exact length; do not replace the raw numeric NPZ with PNG.

## Before unblinding
Inspect cards by neutral case ID. Document each rating without viewing the sealed mapping:
- VISIBLE_DIVISION: clear transition from one *identified physical mother body* to two separate offspring, independent phase-contrast boundary/septation and continued plausible offspring shapes in later frames. Note visible location and frame offset.
- VISIBLE_CONTINUITY_NO_DIVISION: one cell body remains continuous through observation, no new valid cleavage or pair of persistent daughter objects; require sufficient exposure and clear depth of focus.
- TIMING_SHIFT: visible split in same physical region but plausible at another frame offset; classify separately from false/true at software timestamp.
- UNRESOLVED: indistinct cell boundaries, overlap, lighting instability, edge clipping, not enough image evidence, or no unequivocal ancestry.
- Every rating records visibility/quality and free text evidence; *never* force 12 cases into a binary count.

## Controls, unblinding, and interpretation
- Calculate raw pixel checks independently from MoMA labels first. Only after saving blinded ratings, open the sealed key and compare disputed versus matched and ordinary controls.
- Software version agreement is not truth; matching controls contain errors too. Reviewers should not tune thresholds, crop positions, segmentation kernels, or per-event normalization until outcomes look favorable.
- Possible artifact: software daughter at y=120 may be boundary/fragment; this is an E509 hypothesis, NOT a prevalidated visual conclusion. In particular, y=120 evidence cannot itself certify or refute real cell division.
- Only TIFF pixel credibility is tested; even confirmed divisions do not quantify true protein copies, information transmission bit, bit/s, photon budget, or bit/J.
- Preprocessed TIFF is not the raw sensor FOV; camera hardware/time/power calibration remains unverified.
- If full 2D source cannot be certified, mark SOURCE_BLOCKED and keep prior event labels unresolved.
- Human independent annotation is a separate level: an AI/operator examining these images is not an external, independently blinded biologist. Report which review level was actually achieved.

## Frozen outcome records
An event-level CSV should include blind_case_id, classification, evidence_frame_offsets, evidence_region, reviewer_type, image_quality, uncertainty, free_text_evidence, and later the sealed E509 reference after unlock. If no trustworthy direct judgement is possible, all such events remain UNRESOLVED. No p-values or sample-independent CI from 37 technique windows.

## Storage and reproducibility
Save new report/handoff in project 901; large content in 901 separate Drive E510 folder, never write large artifacts to synchronized thesis folder or replace historical records. Verify size, SHA256, ZIP CRC and read-back. Source code and successful/failed Actions runs retain URLs and commit SHA. Not a submission-ready result by itself; PUBLICATION_NO_GO.
