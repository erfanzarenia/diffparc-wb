#!/usr/bin/env python3
"""Merge a native-space cortical parcellation with a warped subcortical atlas into
one hybrid target dseg for ``tck2connectome`` (see rules/target_synthseg_tian.smk).

Cortex comes from SynthSeg ``--parc`` (subject-native, already resliced onto the
node grid); subcortex comes from a template atlas (e.g. Tian S4) warped onto the
same node grid. The two are relabelled into a single contiguous ``1..N`` image
that satisfies the node-builder/​slicer contract (build_voxel_nodes.py /
slice_connectome_block.py): atlas label ``k`` == column ``k`` of the connectivity
matrix, in the order of the target's ``labels:`` list.

Layout (matches the config ``labels:`` order):
    cortex    -> ids 1 .. n_cortex        (fixed by --cortex-lut, NOT by the subject)
    subcortex -> ids n_cortex+1 .. n_cortex+n_sub

Design choices that keep matrices comparable across subjects:
  * Cortical columns are defined by the FIXED lookup, never by ``unique()`` of this
    subject. A parcel absent in a subject is simply an all-zero column; a cortical
    label the segmentation emits that is NOT in the lookup is dropped (and counted),
    so the column set is identical for every subject and every atlas version.
  * At the cortex/subcortex seam the SUBCORTEX wins (it is a precise deep-grey
    atlas), so the cortical ribbon can never bleed into subcortical territory.
"""

import argparse
import csv
import json
import os
import sys

import nibabel as nib
import numpy as np


def load_cortex_lut(path):
    """Return {fs_label:int -> hybrid_idx:int} from a `fs_label  name  idx` TSV."""
    mapping = {}
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            mapping[int(row["fs_label"])] = int(row["idx"])
    if not mapping:
        sys.exit(f"empty cortex lut: {path}")
    return mapping


def same_grid(a, b, atol=1e-3):
    return a.shape == b.shape and np.allclose(
        a.affine, b.affine, atol=atol
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cortex-dseg", required=True,
                    help="Cortical parcellation (SynthSeg --parc) ON THE NODE GRID.")
    ap.add_argument("--cortex-lut", required=True,
                    help="TSV: fs_label<TAB>name<TAB>idx (canonical cortical columns).")
    ap.add_argument("--subcortex-dseg", required=True,
                    help="Subcortical atlas ON THE NODE GRID, labels 1..n_sub.")
    ap.add_argument("--n-cortex", type=int, required=True,
                    help="Number of cortical columns (subcortex is offset by this, "
                         "so columns stay aligned even if a cortical parcel is absent).")
    ap.add_argument("--out-dseg", required=True)
    ap.add_argument("--out-qc")
    args = ap.parse_args()

    cort_img = nib.load(args.cortex_dseg)
    sub_img = nib.load(args.subcortex_dseg)
    if not same_grid(cort_img, sub_img):
        sys.exit(
            "cortex and subcortex are not on the same voxel grid; reslice/warp both "
            "onto the node grid first "
            f"(cortex {cort_img.shape}, subcortex {sub_img.shape})."
        )

    cort = np.rint(np.asarray(cort_img.dataobj)).astype(np.int64)
    sub = np.rint(np.asarray(sub_img.dataobj)).astype(np.int64)
    lut = load_cortex_lut(args.cortex_lut)

    max_lut_idx = max(lut.values())
    if max_lut_idx > args.n_cortex:
        sys.exit(f"cortex lut has idx {max_lut_idx} > --n-cortex {args.n_cortex}")

    out = np.zeros(cort.shape, dtype=np.int32)

    # --- cortex: map fs_label -> hybrid idx via a dense lookup array (fast) ---
    cort_present = np.unique(cort)
    cort_present = cort_present[cort_present > 0]
    max_fs = int(cort_present.max()) if cort_present.size else 0
    remap = np.zeros(max_fs + 1, dtype=np.int32)
    for fs, idx in lut.items():
        if fs <= max_fs:
            remap[fs] = idx
    cort_mapped = remap[np.clip(cort, 0, max_fs)]
    out[cort_mapped > 0] = cort_mapped[cort_mapped > 0]

    # count cortical labels the segmentation produced that are NOT in the lookup
    known = set(lut.keys())
    unmapped = sorted(int(v) for v in cort_present if int(v) not in known)
    n_vox_unmapped = int(np.isin(cort, unmapped).sum()) if unmapped else 0

    # --- subcortex: offset above cortex; SUBCORTEX WINS at any overlap ---
    sub_mask = sub > 0
    overlap = int((sub_mask & (out > 0)).sum())
    out[sub_mask] = sub[sub_mask].astype(np.int32) + args.n_cortex

    # --- write hybrid dseg (fresh header from affine, integer labels) ---
    out_img = nib.Nifti1Image(out, cort_img.affine)
    out_img.set_data_dtype(np.int32)
    os.makedirs(os.path.dirname(args.out_dseg) or ".", exist_ok=True)
    out_img.to_filename(args.out_dseg)

    present = np.unique(out); present = present[present > 0]
    n_sub = int(sub_mask.any() and sub[sub_mask].max())
    stats = {
        "n_cortex_columns": int(args.n_cortex),
        "n_subcortex_labels_seen": int(np.unique(sub[sub_mask]).size),
        "max_subcortex_label": n_sub,
        "n_total_columns": int(args.n_cortex + n_sub),
        "labels_present_in_subject": int(present.size),
        "max_label_id": int(present.max()) if present.size else 0,
        "cortex_voxels": int((out > 0).sum() - sub_mask.sum() + overlap),
        "subcortex_voxels": int(sub_mask.sum()),
        "cortex_subcortex_overlap_voxels": overlap,
        "unmapped_cortical_labels_dropped": unmapped,
        "unmapped_cortical_voxels_dropped": n_vox_unmapped,
    }
    if args.out_qc:
        os.makedirs(os.path.dirname(args.out_qc) or ".", exist_ok=True)
        with open(args.out_qc, "w") as f:
            json.dump(stats, f, indent=2)

    print(
        f"build_hybrid_target: cortex_cols={args.n_cortex} subcortex={n_sub} "
        f"-> {args.n_cortex + n_sub} columns; present_in_subject={present.size}; "
        f"overlap(subcortex_wins)={overlap}vox; "
        f"unmapped_cortical={len(unmapped)} labels/{n_vox_unmapped}vox",
        file=sys.stderr,
    )
    if unmapped:
        print(
            "  NOTE unmapped cortical labels (dropped, not in --cortex-lut): "
            + ", ".join(map(str, unmapped)),
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
