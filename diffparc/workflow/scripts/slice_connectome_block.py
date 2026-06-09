#!/usr/bin/env python3
"""Extract the seed-voxel x atlas-target block from a ``tck2connectome`` matrix.

The connectome is produced by ``tck2connectome`` on the combined node image from
``build_voxel_nodes.py`` (atlas ids 1..n_targets, seed-voxel ids
n_targets+1..n_targets+n_seed). MRtrix sizes the matrix by the maximum node id
and indexes it by node *value*, so:

    rows / cols 0 .. n_targets-1               -> atlas targets (label order)
    rows / cols n_targets .. n_targets+n_seed-1 -> seed voxels (C-order)

We slice rows = seed voxels, cols = atlas targets, and write a voxel x target
CSV with a single header row of target labels -- the same format the previous
custom aggregator emitted, so downstream consumers are unchanged.

This file contains no connectivity math: ``tck2connectome`` already did the
streamline-to-node assignment and SIFT2-weighted accumulation. Here we only
slice and label the resulting matrix.
"""

import argparse
import json
import os
import sys

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--connectome", required=True,
                    help="Full NxN matrix written by tck2connectome (no header).")
    ap.add_argument("--voxel-index", required=True,
                    help="seed_voxel_index.csv (its data-row count gives n_seed).")
    ap.add_argument("--header", required=True,
                    help="Comma-separated target labels (defines n_targets/order).")
    ap.add_argument("--out-matrix", required=True)
    ap.add_argument("--out-qc", required=False)
    ap.add_argument("--empty", action="store_true",
                    help="Degenerate empty tractogram (0 streamlines): write a zero "
                         "(n_seed x n_targets) matrix without reading --connectome.")
    ap.add_argument("--mu-file",
                    help="Optional SIFT2 proportionality coefficient (mu) file (single "
                         "float). When given, --out-matrix is written mu-SCALED "
                         "(fibre-density-calibrated, comparable across subjects/"
                         "conditions); the unscaled block is preserved at "
                         "--out-matrix-raw. mu is a single global scalar, so any later "
                         "row-normalization is unaffected (mu cancels).")
    ap.add_argument("--out-matrix-raw",
                    help="Optional path for the UNSCALED (raw weight-sum) block, kept "
                         "for provenance/reversibility when --mu-file is used.")
    args = ap.parse_args()

    labels = [s.strip() for s in args.header.split(",") if s.strip()]
    n_targets = len(labels)
    if n_targets == 0:
        sys.exit("empty --header; expected comma-separated target labels.")

    with open(args.voxel_index) as f:
        n_seed = sum(1 for _ in f) - 1  # subtract the header line
    if n_seed <= 0:
        sys.exit(f"voxel index has no data rows: {args.voxel_index}")

    n_total = n_targets + n_seed
    if args.empty:
        # 0-streamline tractogram: every seed-voxel x target edge is 0.
        block = np.zeros((n_seed, n_targets), dtype=float)
    else:
        mat = np.atleast_2d(np.loadtxt(args.connectome))
        if mat.shape != (n_total, n_total):
            sys.exit(
                f"connectome is {mat.shape}, expected ({n_total}, {n_total}) "
                f"= n_targets({n_targets}) + n_seed({n_seed}). "
                "A size mismatch usually means a seed voxel id collided with the "
                "atlas range; check that atlas labels are within 1..n_targets."
            )
        # rows = seed voxels, cols = atlas targets
        block = mat[n_targets:n_total, 0:n_targets]

    header = ",".join(labels)

    # Optional SIFT2 mu scaling. mu is a SINGLE GLOBAL scalar from tcksift2, so the
    # mu-scaled matrix is fibre-density-calibrated (comparable across subjects and
    # conditions), while a wrong/ill-defined mu never corrupts the preserved raw
    # block. Because mu multiplies every edge equally, later row-normalization is
    # unchanged (mu cancels exactly).
    mu = None
    if args.mu_file:
        with open(args.mu_file) as f:
            mu = float(f.read().strip().split()[0])
    primary = block if mu is None else (mu * block)

    os.makedirs(os.path.dirname(args.out_matrix) or ".", exist_ok=True)
    np.savetxt(args.out_matrix, primary, delimiter=",", header=header, comments="")

    if args.out_matrix_raw:
        os.makedirs(os.path.dirname(args.out_matrix_raw) or ".", exist_ok=True)
        np.savetxt(args.out_matrix_raw, block, delimiter=",", header=header,
                   comments="")

    if args.out_qc:
        os.makedirs(os.path.dirname(args.out_qc) or ".", exist_ok=True)
        with open(args.out_qc, "w") as f:
            json.dump(
                {
                    "n_seed_voxels": int(n_seed),
                    "n_targets": int(n_targets),
                    "connectome_shape": [n_total, n_total],
                    "block_shape": [int(x) for x in block.shape],
                    "seed_voxels_with_any_connection": int(
                        (block.sum(axis=1) > 0).sum()
                    ),
                    "total_block_weight": float(block.sum()),  # raw weight-sum
                    "sift2_mu": (None if mu is None else float(mu)),
                    "mu_scaled": bool(mu is not None),
                    "total_block_weight_scaled": (
                        None if mu is None else float(primary.sum())
                    ),
                },
                f,
                indent=2,
            )

    scaled_note = "" if mu is None else f" (mu-scaled, mu={mu:.6g})"
    print(f"slice_connectome_block: wrote {block.shape[0]}x{block.shape[1]} "
          f"(seed x target) matrix{scaled_note}", file=sys.stderr)


if __name__ == "__main__":
    main()
