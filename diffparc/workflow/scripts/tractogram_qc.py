#!/usr/bin/env python3
"""Generate visual-QC artifacts for one tractogram into the /qc tree.

Produces (all via MRtrix3 except the histogram PNG):
  * TDI            -- track-density image                   (tckmap)
  * endpoint map   -- endpoint density                      (tckmap -ends_only)
  * DEC map        -- direction-encoded-colour TDI, SIFT2-weighted when weights
                      are supplied                           (tckmap -dec)
  * 200k subset    -- reproducible RANDOM subsample for mrview inspection
                      (reservoir sampling; memory-safe via lazy streamline load)
  * length CSV     -- full-tractogram streamline-length histogram
                      (tckstats -histogram)
  * length PNG     -- the same distribution from the subsample, plotted

The track-weighted DEC here and the FOD-based DEC (fod2dec, a separate rule) are
complementary: overlay the tracks / track subset on the FOD DEC for an
anatomical background. No connectivity math; QC only.
"""

import argparse
import os
import subprocess

import numpy as np
import nibabel as nib


def run(cmd):
    """Run a command, echoing it; inherits stdout/stderr (rule redirects to log)."""
    print("+ " + " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True)


def reservoir_subsample(in_tck, out_tck, n_keep, seed):
    """Uniform random subsample of <= n_keep streamlines in a single streaming pass.

    Reservoir sampling keeps memory bounded to n_keep streamlines, so this is safe
    on multi-million-streamline whole-brain tractograms. Reproducible via `seed`.
    """
    tf = nib.streamlines.load(in_tck, lazy_load=True)
    rng = np.random.default_rng(seed)
    reservoir = []
    n_seen = 0
    for s in tf.streamlines:
        if n_seen < n_keep:
            reservoir.append(np.asarray(s).copy())
        else:
            j = int(rng.integers(0, n_seen + 1))
            if j < n_keep:
                reservoir[j] = np.asarray(s).copy()
        n_seen += 1

    out = nib.streamlines.Tractogram(
        streamlines=reservoir,
        affine_to_rasmm=tf.tractogram.affine_to_rasmm,
    )
    nib.streamlines.save(out, out_tck, header=tf.header)
    return reservoir, n_seen


def length_mm(streamline):
    pts = np.asarray(streamline)
    if pts.shape[0] < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())


def plot_length_hist(streamlines, out_png, n_total):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lengths = [length_mm(s) for s in streamlines]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    if lengths:
        ax.hist(lengths, bins=60, color="#4060c0", edgecolor="none")
        ax.set_xlabel("streamline length (mm)")
        ax.set_ylabel("count")
    else:
        ax.text(0.5, 0.5, "no streamlines", ha="center", va="center")
        ax.set_axis_off()
    ax.set_title(
        f"streamline length distribution "
        f"(n={len(lengths)} sampled of {n_total} total)"
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tck", required=True)
    ap.add_argument("--template", required=True, help="Grid image for tckmap (DWI mask).")
    ap.add_argument("--weights", default="",
                    help="SIFT2 weights for the DEC map; '' = unweighted.")
    ap.add_argument("--out-tdi", required=True)
    ap.add_argument("--out-endpoints", required=True)
    ap.add_argument("--out-dec", required=True)
    ap.add_argument("--out-subset", required=True)
    ap.add_argument("--out-lengths-csv", required=True)
    ap.add_argument("--out-lengths-png", required=True)
    ap.add_argument("--n-subsample", type=int, default=200000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--nthreads", type=int, default=1)
    args = ap.parse_args()

    for p in (args.out_tdi, args.out_endpoints, args.out_dec, args.out_subset,
              args.out_lengths_csv, args.out_lengths_png):
        os.makedirs(os.path.dirname(p) or ".", exist_ok=True)

    nt = ["-nthreads", args.nthreads]

    # 1. TDI (unweighted track-count density)
    run(["tckmap", "-force", *nt, "-template", args.template,
         args.tck, args.out_tdi])

    # 2. Endpoint density map
    run(["tckmap", "-force", *nt, "-template", args.template, "-ends_only",
         args.tck, args.out_endpoints])

    # 3. Track-weighted DEC map (SIFT2-weighted when weights are given)
    dec = ["tckmap", "-force", *nt, "-template", args.template, "-dec"]
    if args.weights:
        dec += ["-tck_weights_in", args.weights]
    dec += [args.tck, args.out_dec]
    run(dec)

    # 4. Full-tractogram streamline-length histogram -> CSV
    run(["tckstats", "-force", "-histogram", args.out_lengths_csv, args.tck])

    # 5. Reproducible random 200k subset (+ representative length histogram PNG)
    reservoir, n_total = reservoir_subsample(
        args.tck, args.out_subset, args.n_subsample, args.seed)
    plot_length_hist(reservoir, args.out_lengths_png, n_total)

    print(f"tractogram_qc: n_total={n_total} subset={len(reservoir)}", flush=True)


if __name__ == "__main__":
    main()
