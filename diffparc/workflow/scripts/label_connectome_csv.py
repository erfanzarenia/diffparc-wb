#!/usr/bin/env python3
"""Attach a parcel-label header to a raw tck2connectome matrix and write CSV.

``tck2connectome`` writes an NxN matrix (whitespace-separated, no header) indexed
by node value 1..N. This wraps it into a comma-separated CSV whose first row is
the atlas parcel labels -- matching the header style of the pipeline's other
connectome CSVs (slice_connectome_block.py). No connectivity math here.
"""

import argparse
import os
import sys

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True,
                    help="Raw NxN matrix written by tck2connectome.")
    ap.add_argument("--labels", required=True,
                    help="Comma-separated atlas parcel labels (defines N / order).")
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()

    labels = [s.strip() for s in args.labels.split(",") if s.strip()]
    n = len(labels)
    if n == 0:
        sys.exit("empty --labels; expected comma-separated parcel labels.")

    mat = np.atleast_2d(np.loadtxt(args.matrix))
    if mat.shape != (n, n):
        sys.exit(
            f"connectome is {mat.shape}, expected ({n}, {n}) to match "
            f"{n} atlas labels. A mismatch usually means the atlas dseg has node "
            "values outside 1..N or the labels list is out of sync with the atlas."
        )

    os.makedirs(os.path.dirname(args.out_csv) or ".", exist_ok=True)
    np.savetxt(args.out_csv, mat, delimiter=",",
               header=",".join(labels), comments="")

    print(f"label_connectome_csv: wrote {n}x{n} matrix -> {args.out_csv}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
