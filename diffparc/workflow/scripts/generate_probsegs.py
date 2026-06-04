#!/usr/bin/env python3
"""
generate_probsegs.py  --  ONE-TIME preprocessing utility.

This script is NOT invoked by any Snakemake rule. Run it once to generate the
hemisphere-split template probability seeds that then sit in
resources/.../CIT168_prob/"""

import argparse
import os
import sys

import nibabel as nib
import numpy as np

# CIT168 label name -> 4D volume index (validated against labels.txt below).
CIT168_INDEX = {"SNc": 6, "SNr": 8, "PBP": 9, "VTA": 10}

# Output seed name -> filename stem written into the prob dir.
OUTPUT_SEEDS = ("snc", "vtapbp", "vtasncpbp")


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def validate_labels(prob_dir, index):
    """Best-effort check that labels.txt agrees with the hard-coded indices."""
    labels_txt = os.path.join(prob_dir, "labels.txt")
    if not os.path.exists(labels_txt):
        log("[warn] %s not found; skipping label validation." % labels_txt)
        return
    idx_to_name = {}
    with open(labels_txt) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                idx_to_name[int(parts[0])] = parts[1]
    for name, idx in index.items():
        found = idx_to_name.get(idx)
        if found != name:
            raise SystemExit(
                "[error] labels.txt mismatch: expected index %d == %s, but found "
                "%r. Update CIT168_INDEX for this atlas version." % (idx, name, found)
            )
    log("[ok] labels.txt matches expected SNc/SNr/PBP/VTA indices.")


def load(path):
    if not os.path.exists(path):
        raise SystemExit("[error] missing input: %s" % path)
    img = nib.load(path)
    if len(img.shape) != 3:
        raise SystemExit("[error] expected 3D image, got shape %r: %s" % (img.shape, path))
    return img


def determine_left_is_negative_x(det_dir, left_name, right_name):
    """
    Decide whether anatomical LEFT corresponds to world-x < 0, by comparing the
    centroid world-x of the atlas's own SNc_left vs SNc_right deterministic
    masks. Centroids are computed via each file's own affine, so this is
    independent of grid/orientation. Falls back to standard RAS if unavailable.
    """
    lpath = os.path.join(det_dir, left_name)
    rpath = os.path.join(det_dir, right_name)
    if not (os.path.exists(lpath) and os.path.exists(rpath)):
        log("[warn] %s / %s not found in %s; assuming standard RAS "
            "(left = world-x < 0)." % (left_name, right_name, det_dir))
        return True

    def centroid_x(p):
        img = nib.load(p)
        d = np.asarray(img.dataobj)
        vox = np.argwhere(d > 0)
        if vox.size == 0:
            return None
        world = nib.affines.apply_affine(img.affine, vox.astype(np.float64))
        return float(world[:, 0].mean())

    lx, rx = centroid_x(lpath), centroid_x(rpath)
    if lx is None or rx is None:
        log("[warn] empty SNc_left/right mask; assuming standard RAS "
            "(left = world-x < 0).")
        return True

    left_is_negative = lx < rx
    log("[ok] L/R sign verified from atlas: SNc_left mean-x=%.2f, "
        "SNc_right mean-x=%.2f -> anatomical LEFT is world-x %s 0"
        % (lx, rx, "<" if left_is_negative else ">"))
    return left_is_negative


def world_x_volume(shape, affine):
    """world-x coordinate of every voxel center, as an (I,J,K) array."""
    I, J, K = shape
    ii = np.arange(I).reshape(I, 1, 1)
    jj = np.arange(J).reshape(1, J, 1)
    kk = np.arange(K).reshape(1, 1, K)
    a = affine
    return (a[0, 0] * ii + a[0, 1] * jj + a[0, 2] * kk + a[0, 3]).astype(np.float64)


def main():
    tpl_dir = os.path.join("resources", "tpl-MNI152NLin2009cAsym")
    default_prob = os.path.join(tpl_dir, "CIT168_prob")
    default_det = os.path.join(tpl_dir, "CIT168")

    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--prob-dir", default=default_prob,
                    help="CIT168_prob dir (per-label prob volumes; default: %(default)s)")
    ap.add_argument("--det-dir", default=default_det,
                    help="CIT168 deterministic dir, for SNc_left/right L/R sign "
                         "verification (default: %(default)s)")
    ap.add_argument("--prob-pattern", default="CIT168toMNI152-2009c_prob_{idx:04d}.nii.gz",
                    help="per-label prob filename pattern (default: %(default)s)")
    ap.add_argument("--snc-left-name", default="SNc_left.nii.gz",
                    help="atlas SNc left mask for sign check (default: %(default)s)")
    ap.add_argument("--snc-right-name", default="SNc_right.nii.gz",
                    help="atlas SNc right mask for sign check (default: %(default)s)")
    ap.add_argument("--out-pattern", default="{seed}_hemi-{hemi}_probseg.nii.gz",
                    help="output filename pattern (default: %(default)s)")
    ap.add_argument("--hemis", nargs="+", default=["L", "R"])
    ap.add_argument("--overwrite", action="store_true",
                    help="replace existing outputs instead of erroring")
    args = ap.parse_args()

    prob_dir = args.prob_dir
    validate_labels(prob_dir, CIT168_INDEX)

    def prob_path(idx):
        return os.path.join(prob_dir, args.prob_pattern.format(idx=idx))

    snc_img = load(prob_path(CIT168_INDEX["SNc"]))
    vta_img = load(prob_path(CIT168_INDEX["VTA"]))
    pbp_img = load(prob_path(CIT168_INDEX["PBP"]))
    snr_img = load(prob_path(CIT168_INDEX["SNr"]))   # for SNr-exclusion guard

    # The per-label volumes are slices of one 4D atlas, so they must share a grid;
    # all voxel math (and the output) is on this shared component grid.
    for name, img in (("VTA", vta_img), ("PBP", pbp_img), ("SNr", snr_img)):
        if img.shape != snc_img.shape or not np.allclose(img.affine, snc_img.affine):
            raise SystemExit(
                "[error] %s grid differs from SNc; per-label prob volumes must "
                "share a grid." % name
            )

    affine = snc_img.affine
    snc = snc_img.get_fdata(dtype=np.float64)
    vta = vta_img.get_fdata(dtype=np.float64)
    pbp = pbp_img.get_fdata(dtype=np.float64)
    snr = snr_img.get_fdata(dtype=np.float64)

    # Probability merges.
    # snc / vtapbp are RAW component probs here (no WTA). They overlap at the
    # boundary on purpose: the winner-take-all is done once in subject space
    # (sweep_binarize_seed) AFTER warping, so we never warp a sharp WTA edge with
    # linear interpolation (which re-smears it). Masking to the union happens in
    # the write loop below.
    a10 = np.clip(vta + pbp, 0.0, 1.0)               # mesolimbic complex
    da_all = np.clip(snc + vta + pbp, 0.0, 1.0)      # full DA midbrain
    snc_ps = snc                                      # raw P_SNc
    vtapbp_ps = a10                                   # raw P_VTA + P_PBP

    seed_full = {"snc": snc_ps, "vtapbp": vtapbp_ps, "vtasncpbp": da_all}

    # Hemisphere/extent selector: full SNc U VTA U PBP union, with an SNr-
    # exclusion guard. SNr is GABAergic and must be excluded entirely, but the
    # broad union (needed to capture all of PBP) otherwise sweeps in a rim of
    # voxels where SNr is the dominant label. Requiring the DA complex to win
    # over SNr (da_all >= P_SNr) drops that rim while keeping ~94% of PBP.
    union = ((snc > 0) | (vta > 0) | (pbp > 0)) & (da_all >= snr)
    n_snr_excluded = int((((snc > 0) | (vta > 0) | (pbp > 0)) & (da_all < snr)).sum())
    log("[info] SNr-exclusion guard removed %d SNr-dominant voxels from the "
        "selector." % n_snr_excluded)
    left_is_negative = determine_left_is_negative_x(
        args.det_dir, args.snc_left_name, args.snc_right_name
    )
    wx = world_x_volume(snc.shape, affine)
    left_side = (wx < 0.0) if left_is_negative else (wx > 0.0)
    hemi_side = {"L": left_side, "R": ~left_side}

    n_union = int(union.sum())
    log("[info] union(SNc,VTA,PBP) voxels = %d  (L=%d, R=%d)"
        % (n_union, int((union & hemi_side["L"]).sum()),
           int((union & hemi_side["R"]).sum())))

    for hemi in args.hemis:
        sel = union & hemi_side[hemi]
        for seed in OUTPUT_SEEDS:
            data = (seed_full[seed] * sel).astype(np.float32)
            out_path = os.path.join(
                prob_dir, args.out_pattern.format(seed=seed, hemi=hemi)
            )
            if os.path.exists(out_path) and not args.overwrite:
                raise SystemExit(
                    "[error] %s exists; pass --overwrite to replace." % out_path
                )
            out_img = nib.Nifti1Image(data, affine)
            out_img.set_data_dtype(np.float32)
            nib.save(out_img, out_path)
            log("[ok] %s  nonzero=%d  max=%.3f"
                % (out_path, int((data > 0).sum()),
                   float(data.max() if data.size else 0.0)))

    # Sanity: snc U vtapbp == vtasncpbp (voxel sets), per hemi. snc and vtapbp
    # are raw components here so they OVERLAP at the boundary (expected); the WTA
    # that makes them disjoint runs later in subject space (sweep_binarize_seed).
    for hemi in args.hemis:
        sel = union & hemi_side[hemi]
        s = (seed_full["snc"] * sel) > 0
        v = (seed_full["vtapbp"] * sel) > 0
        c = (seed_full["vtasncpbp"] * sel) > 0
        if not np.array_equal(s | v, c):
            log("[warn] hemi-%s: (snc U vtapbp) != vtasncpbp voxel set "
                "(unexpected; check merge logic)." % hemi)
        else:
            log("[ok] hemi-%s: snc U vtapbp == vtasncpbp (%d voxels); "
                "snc/vtapbp overlap = %d (expected, resolved by subject-space WTA)."
                % (hemi, int(c.sum()), int((s & v).sum())))

    log("[done] Generated %s hemi probsegs in %s"
        % ("/".join(OUTPUT_SEEDS), prob_dir))


if __name__ == "__main__":
    main()
