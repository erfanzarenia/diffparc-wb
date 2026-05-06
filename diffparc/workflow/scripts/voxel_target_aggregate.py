#!/usr/bin/env python3
import argparse
import sys
import numpy as np
import nibabel as nib
import os


def read_tck_endpoints(tck_path: str):
    """
    Minimal .tck reader returning (first_point, last_point) per streamline.

    MRtrix .tck uses:
      - NaN NaN NaN as streamline delimiter
      - Inf Inf Inf as end-of-file marker

    We must stop on the Inf marker, or we risk counting an extra streamline.
    """
    with open(tck_path, "rb") as f:
        header = b""
        while True:
            line = f.readline()
            if not line:
                raise RuntimeError("Unexpected EOF while reading .tck header")
            header += line
            if line.strip() == b"END":
                break

        header_text = header.decode("utf-8", errors="ignore")
        offset = None
        for ln in header_text.splitlines():
            if ln.lower().startswith("file:"):
                parts = ln.split()
                if len(parts) >= 3:
                    try:
                        offset = int(parts[2])
                    except ValueError:
                        pass
        if offset is None:
            offset = f.tell()

        f.seek(offset, os.SEEK_SET)

        chunk_points = 1_000_000  # points per chunk
        first_pts = []
        last_pts = []

        cur_first = None
        prev = None

        while True:
            data = np.fromfile(f, dtype=np.float32, count=chunk_points * 3)
            if data.size == 0:
                break
            data = data.reshape((-1, 3))

            for p in data:
                # EOF marker
                if np.isinf(p[0]) and np.isinf(p[1]) and np.isinf(p[2]):
                    if cur_first is not None and prev is not None:
                        first_pts.append(cur_first)
                        last_pts.append(prev)
                    return (np.asarray(first_pts, dtype=np.float32),
                            np.asarray(last_pts, dtype=np.float32))

                # streamline delimiter
                if np.isnan(p[0]) and np.isnan(p[1]) and np.isnan(p[2]):
                    if cur_first is not None and prev is not None:
                        first_pts.append(cur_first)
                        last_pts.append(prev)
                    cur_first = None
                    prev = None
                    continue

                # regular point
                if cur_first is None:
                    cur_first = p.copy()
                prev = p.copy()

        # If the file ends without Inf marker (uncommon), finalize
        if cur_first is not None and prev is not None:
            first_pts.append(cur_first)
            last_pts.append(prev)

        return np.asarray(first_pts, dtype=np.float32), np.asarray(last_pts, dtype=np.float32)


def parse_header_labels(header: str):
    labels = [x.strip() for x in header.split(",") if x.strip()]
    if not labels:
        raise RuntimeError("Header is empty. Provide --header as comma-separated labels.")
    return labels


def enforce_atlas_labels_match_header(targets_nii_path: str, n_targets: int):
    """
    Enforce that the atlas label IDs are exactly {1..n_targets}, no gaps.
    This prevents silent column mislabeling.
    """
    timg = nib.load(targets_nii_path)
    tdat = timg.get_fdata()

    labs = np.asarray(np.rint(tdat), dtype=np.int64).ravel()
    labs = labs[labs > 0]
    if labs.size == 0:
        raise RuntimeError(f"Targets atlas has no positive labels: {targets_nii_path}")

    uniq = np.unique(labs)
    expected = np.arange(1, n_targets + 1, dtype=np.int64)

    if uniq.size != expected.size or not np.array_equal(uniq, expected):
        raise RuntimeError(
            "Targets atlas label IDs do not match header length/order.\n"
            f"  header n_targets = {n_targets}\n"
            f"  atlas unique labels (min..max) = {uniq[0]}..{uniq[-1]} (count {uniq.size})\n"
            "Expected atlas labels to be exactly 1..n_targets with no gaps.\n"
            "Fix: ensure your atlas is relabeled to 1..N OR update header+aggregation logic."
        )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed-tck", required=True)
    ap.add_argument("--seed-weights", required=True)
    ap.add_argument("--seed-mask", required=True)
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--header", required=True, help="Comma-separated target labels (column names)")
    ap.add_argument("--targets-nii", required=True, help="Targets atlas NIfTI used for label-ID validation")
    ap.add_argument("--out-weighted", required=True)
    ap.add_argument("--out-counts", required=True)
    ap.add_argument("--out-voxel-index", required=True)
    args = ap.parse_args()

    # Load seed mask
    seed_img = nib.load(args.seed_mask)
    seed_data = np.asarray(seed_img.dataobj)
    seed_aff = seed_img.affine
    inv_seed_aff = np.linalg.inv(seed_aff)

    seed_vox = np.array(np.argwhere(seed_data > 0), dtype=np.int32)  # (N,3) ijk
    if seed_vox.size == 0:
        raise RuntimeError("Seed mask is empty (no voxels > 0).")

    # Build row index mapping (voxel ijk -> row)
    I, J, K = seed_data.shape
    seed_hash = seed_vox[:, 0] + I * (seed_vox[:, 1] + J * seed_vox[:, 2])
    row_index = {int(h): idx for idx, h in enumerate(seed_hash)}
    n_seed = seed_vox.shape[0]

    # Header labels from config (Option A)
    target_labels = parse_header_labels(args.header)
    n_targets = len(target_labels)

    # HARD SAFETY CHECK: atlas labels must be exactly 1..n_targets
    enforce_atlas_labels_match_header(args.targets_nii, n_targets)

    # Endpoints from tck
    first_pts, last_pts = read_tck_endpoints(args.seed_tck)
    n_sl = first_pts.shape[0]

    # Weights (one per streamline)
    try:
        w = np.loadtxt(args.seed_weights, dtype=np.float64)
    except Exception as e:
        raise RuntimeError(
            f"Failed to read weights file: {args.seed_weights} ({type(e).__name__})"
        ) from None
    w = np.atleast_1d(w)

    # Assignments (node_i node_j per streamline)
    assign = np.loadtxt(args.assignments, dtype=np.int64)
    if assign.ndim == 1 and assign.size == 2:
        assign = assign.reshape((1, 2))

    # Hard sanity check: these should match exactly
    if w.shape[0] != n_sl:
        raise RuntimeError(f"Weight count mismatch: {w.shape[0]} weights vs {n_sl} streamlines")
    if assign.shape[0] != n_sl or assign.shape[1] < 2:
        raise RuntimeError(f"Assignment shape mismatch: {assign.shape} vs expected ({n_sl},2)")

    node_a = assign[:, 0]
    node_b = assign[:, 1]

    # Choose the nonzero node as target; if both set, this picks node_a
    # (for your use-case seed is NOT in the atlas, so exactly one endpoint should be labeled)
    target_node = np.where(node_a > 0, node_a, node_b)
    target_node = np.where(target_node > 0, target_node, np.maximum(node_a, node_b))

    # Convert endpoints to seed-mask voxel indices
    ones = np.ones((n_sl, 1), dtype=np.float64)
    first_h = np.hstack([first_pts.astype(np.float64), ones])
    last_h = np.hstack([last_pts.astype(np.float64), ones])

    first_ijk = (inv_seed_aff @ first_h.T).T[:, :3]
    last_ijk = (inv_seed_aff @ last_h.T).T[:, :3]

    # Robust-ish voxelization
    first_ijk = np.floor(first_ijk + 1e-6).astype(np.int64)
    last_ijk = np.floor(last_ijk + 1e-6).astype(np.int64)

    def in_bounds(ijk):
        return (
            (ijk[:, 0] >= 0) & (ijk[:, 0] < I) &
            (ijk[:, 1] >= 0) & (ijk[:, 1] < J) &
            (ijk[:, 2] >= 0) & (ijk[:, 2] < K)
        )

    f_ok = in_bounds(first_ijk)
    l_ok = in_bounds(last_ijk)

    f_in_seed = np.zeros(n_sl, dtype=bool)
    l_in_seed = np.zeros(n_sl, dtype=bool)

    f_idx = first_ijk[f_ok]
    l_idx = last_ijk[l_ok]

    f_in_seed[f_ok] = seed_data[f_idx[:, 0], f_idx[:, 1], f_idx[:, 2]] > 0
    l_in_seed[l_ok] = seed_data[l_idx[:, 0], l_idx[:, 1], l_idx[:, 2]] > 0

    seed_ijk = np.full((n_sl, 3), -1, dtype=np.int64)
    seed_ijk[f_in_seed] = first_ijk[f_in_seed]
    use_last = (~f_in_seed) & l_in_seed
    seed_ijk[use_last] = last_ijk[use_last]

    # Valid streamlines must have a seed voxel + a valid target node id
    valid = (seed_ijk[:, 0] >= 0) & (target_node > 0) & (target_node <= n_targets)

    weighted = np.zeros((n_seed, n_targets), dtype=np.float64)
    counts = np.zeros((n_seed, n_targets), dtype=np.int32)

    tgt_col = target_node[valid] - 1  # 0-based col
    s = seed_ijk[valid]
    s_hash = s[:, 0] + I * (s[:, 1] + J * s[:, 2])

    # --- vectorized mapping from s_hash -> row index ---

    # sort seed_hash once
    sort_idx = np.argsort(seed_hash)
    seed_hash_sorted = seed_hash[sort_idx]

    # find positions of s_hash in sorted seed_hash
    pos = np.searchsorted(seed_hash_sorted, s_hash)

    # check which are valid matches
    in_range = pos < seed_hash_sorted.size
    matched = np.zeros_like(in_range, dtype=bool)
    matched[in_range] = seed_hash_sorted[pos[in_range]] == s_hash[in_range]

    # map back to original row indices
    rows = np.full(s_hash.shape, -1, dtype=np.int64)
    rows[matched] = sort_idx[pos[matched]]

    # filter valid rows
    ok = rows >= 0

    rows = rows[ok]
    tgt_col = tgt_col[ok]
    wv = w[valid][ok]

    np.add.at(weighted, (rows, tgt_col), wv)
    np.add.at(counts, (rows, tgt_col), 1.0)

    # voxel index mapping
    vox_world = nib.affines.apply_affine(seed_aff, seed_vox.astype(np.float64))
    out_map = np.column_stack([
        np.arange(n_seed, dtype=np.int64),
        seed_vox.astype(np.int64),
        vox_world.astype(np.float64),
    ])
    np.savetxt(
        args.out_voxel_index,
        out_map,
        delimiter=",",
        header="row,i,j,k,x,y,z",
        comments=""
    )

    header = ",".join(target_labels)
    np.savetxt(args.out_weighted, weighted, delimiter=",", header=header, comments="")
    np.savetxt(args.out_counts, counts, delimiter=",", header=header, comments="")

    return 0


if __name__ == "__main__":
    sys.exit(main())
