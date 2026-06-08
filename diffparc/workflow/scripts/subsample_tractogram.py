#!/usr/bin/env python3
"""Reproducible random subsampling of an MRtrix .tck tractogram.

Used by the enrichment_sweep module: the maximal ROI tractogram is generated
once, and every lower enrichment level is obtained by random subsampling ONLY
(no re-seeding / no re-running tckgen).

The number of streamlines kept is ``round(n_total * fraction)`` where
``fraction = level / max_level`` (clamped to <= n_total).

Subsampling is *nested*: a single, level-INDEPENDENT random permutation of the
streamline indices is drawn (seeded by ``seed`` only), and each level keeps the
first ``n_keep`` of that fixed permutation. Because the permutation is the same
for every level, the kept set at a lower level is a strict subset of every
higher level (e.g. 1000 ⊂ 3000 ⊂ 5000 — each higher level is the lower one plus
extra streamlines). Indices are sorted before extraction so the output preserves
the original streamline order.

The selected indices (and ``n_total``) are also written to an optional meta JSON
so a companion step can propagate per-streamline SIFT2 weights from the maximal
tractogram to this level without re-fitting SIFT2.

This performs no connectivity math and does not touch SIFT2 — it only writes a
smaller tractogram (and, optionally, the selection metadata).
"""

import argparse
import json
import os
import sys

import numpy as np
import nibabel as nib


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-tck", required=True, help="Input (maximal) ROI tractogram.")
    ap.add_argument("--out-tck", required=True, help="Output subsampled tractogram.")
    ap.add_argument("--level", type=int, required=True, help="This enrichment level.")
    ap.add_argument("--max-level", type=int, required=True, help="Maximal enrichment level.")
    ap.add_argument("--seed", type=int, default=42, help="Base RNG seed (level-independent).")
    ap.add_argument("--out-count", required=False, help="Write kept streamline count here.")
    ap.add_argument("--out-meta", required=False,
                    help="Write selection metadata (n_total, indices, ...) as JSON here.")
    args = ap.parse_args()

    if args.max_level <= 0:
        sys.exit("--max-level must be > 0")
    if args.level <= 0:
        sys.exit("--level must be > 0")
    if args.level > args.max_level:
        sys.exit(f"--level ({args.level}) cannot exceed --max-level ({args.max_level})")

    tf = nib.streamlines.load(args.in_tck)
    streamlines = tf.streamlines
    n_total = len(streamlines)

    fraction = args.level / args.max_level
    n_keep = int(round(n_total * fraction))
    n_keep = max(0, min(n_keep, n_total))

    os.makedirs(os.path.dirname(args.out_tck) or ".", exist_ok=True)

    # Nested selection: one fixed permutation (seeded by `seed` ONLY, not by
    # level), sliced to a prefix of size n_keep. Identical permutation across
    # levels => prefixes are nested => lower levels are strict subsets of higher
    # ones. Sorting keeps the output in original streamline order. At the maximal
    # level n_keep == n_total, so perm[:n_keep] sorted == arange(n_total) and
    # everything is kept (the superset of all lower levels).
    rng = np.random.default_rng(args.seed)
    perm = rng.permutation(n_total)
    idx = np.sort(perm[:n_keep]).astype(np.int64)

    new_streamlines = [streamlines[int(i)] for i in idx]

    out_tractogram = nib.streamlines.Tractogram(
        streamlines=new_streamlines,
        affine_to_rasmm=tf.tractogram.affine_to_rasmm,
    )
    nib.streamlines.save(out_tractogram, args.out_tck, header=tf.header)

    if args.out_count:
        os.makedirs(os.path.dirname(args.out_count) or ".", exist_ok=True)
        with open(args.out_count, "w") as f:
            f.write(str(len(idx)) + "\n")

    if args.out_meta:
        os.makedirs(os.path.dirname(args.out_meta) or ".", exist_ok=True)
        with open(args.out_meta, "w") as f:
            json.dump(
                {
                    "n_total": int(n_total),
                    "n_keep": int(len(idx)),
                    "level": int(args.level),
                    "max_level": int(args.max_level),
                    "seed": int(args.seed),
                    "fraction": float(fraction),
                    "nested": True,
                    "indices": [int(i) for i in idx],
                },
                f,
            )

    print(
        f"subsample_tractogram: level={args.level} max={args.max_level} "
        f"fraction={fraction:.4f} n_total={n_total} n_keep={len(idx)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
