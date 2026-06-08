#!/usr/bin/env python3
"""Reproducible random subsampling of an MRtrix .tck tractogram.

Used by the enrichment_sweep module: the maximal ROI tractogram is generated
once, and every lower enrichment level is obtained by random subsampling ONLY
(no re-seeding / no re-running tckgen).

The number of streamlines kept is ``round(n_total * fraction)`` where
``fraction = level / max_level`` (clamped to <= n_total). The selection is a
uniform random draw without replacement, made reproducible by seeding NumPy's
RNG with ``(seed, level)``. Selected indices are kept in their original order so
the output is a faithful random subset of the input.

This performs no connectivity math and does not touch SIFT2 — it only writes a
smaller tractogram.
"""

import argparse
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
    ap.add_argument("--seed", type=int, default=42, help="Base RNG seed.")
    ap.add_argument("--out-count", required=False, help="Write kept streamline count here.")
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

    if n_keep >= n_total:
        # Maximal level (fraction 1.0): keep everything, no random draw needed.
        idx = np.arange(n_total, dtype=np.int64)
    else:
        rng = np.random.default_rng([args.seed, args.level])
        idx = rng.choice(n_total, size=n_keep, replace=False)
        idx.sort()  # preserve original order

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

    print(
        f"subsample_tractogram: level={args.level} max={args.max_level} "
        f"fraction={fraction:.4f} n_total={n_total} n_keep={len(idx)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
