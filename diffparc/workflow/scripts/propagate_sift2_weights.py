#!/usr/bin/env python3
"""Propagate SIFT2 weights from the maximal combined tractogram to one level.

Branch B of the enrichment sweep: SIFT2 is fit ONCE on the maximal, fully
enriched combined tractogram ``WB ++ ROI_max`` (rule ``enrich_sift2`` at the
maximal level). Lower levels are NOT re-fit; instead we reuse those weights.

The maximal-fit weight vector is ordered exactly like the tractogram it was fit
on, i.e. the whole-brain block first, then the ROI_max block:

    weights_max = [ w_wb(0 .. n_wb-1) | w_roimax(0 .. n_total-1) ]

A lower level's combined tractogram is ``WB ++ ROI_subset`` where ROI_subset is
the nested random subset chosen by ``subsample_tractogram.py`` (ascending
original-index order). The propagated weights are therefore::

    out = [ w_wb ]                       (WB baseline condition), or
    out = [ w_wb | w_roimax[indices] ]   (an enrichment level),

which align 1:1 with that level's combined tractogram. ``n_wb`` is recovered as
``len(weights_max) - n_total`` (n_total = ROI_max streamline count, read from the
maximal-level subsample meta), so no separate streamline count has to be passed.

No connectivity math, no SIFT2 re-fit, no tractogram modification — only a
reindexing of an existing weight vector.
"""

import argparse
import json
import sys

import numpy as np


def load_meta(path):
    with open(path) as f:
        return json.load(f)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--weights-max", required=True,
                    help="SIFT2 weights of the maximal combined tractogram (WB ++ ROI_max).")
    ap.add_argument("--roi-max-meta", required=True,
                    help="Subsample meta JSON at the MAXIMAL level; supplies n_total (ROI_max count).")
    ap.add_argument("--indices", default="",
                    help="Subsample meta JSON for THIS level (nested ROI indices). "
                         "Omit/empty for the WB baseline (no ROI streamlines).")
    ap.add_argument("--out-weights", required=True)
    args = ap.parse_args()

    w = np.loadtxt(args.weights_max, dtype=float, ndmin=1)
    n_total = int(load_meta(args.roi_max_meta)["n_total"])
    n_wb = len(w) - n_total
    if n_wb < 0:
        sys.exit(
            f"weights-max length ({len(w)}) < n_total ({n_total}); "
            "the maximal SIFT2 weights and ROI_max meta are inconsistent."
        )
    wb_block = w[:n_wb]
    roi_block = w[n_wb:]

    if args.indices:
        idx = np.asarray(load_meta(args.indices)["indices"], dtype=np.int64)
        if idx.size and (idx.min() < 0 or idx.max() >= n_total):
            sys.exit(
                f"ROI indices out of range [0, {n_total}) for the ROI_max weight block."
            )
        out = np.concatenate([wb_block, roi_block[idx]]) if idx.size else wb_block
        n_roi_kept = int(idx.size)
    else:
        out = wb_block  # WB baseline: no ROI streamlines
        n_roi_kept = 0

    # Single-line, space-separated -- matches MRtrix's own weight-file style and
    # is consumed unchanged by tckedit/tck2connectome -tck_weights_in.
    np.savetxt(args.out_weights, out.reshape(1, -1), fmt="%.8e")

    print(
        f"propagate_sift2_weights: n_wb={n_wb} n_total={n_total} "
        f"n_roi_kept={n_roi_kept} n_out={out.size}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
