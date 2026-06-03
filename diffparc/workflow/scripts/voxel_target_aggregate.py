#!/usr/bin/env python3

import argparse
import os
import sys
import json

import nibabel as nib
import numpy as np


def log(msg: str):
    print(msg, file=sys.stderr, flush=True)


def read_tck_endpoints(tck_path: str):
    """
    Minimal .tck reader returning first and last point per streamline.

    MRtrix .tck uses:
      NaN NaN NaN as streamline delimiter
      Inf Inf Inf as end-of-file marker
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

        chunk_points = 1_000_000
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
                if np.isinf(p[0]) and np.isinf(p[1]) and np.isinf(p[2]):
                    if cur_first is not None and prev is not None:
                        first_pts.append(cur_first)
                        last_pts.append(prev)

                    return (
                        np.asarray(first_pts, dtype=np.float32),
                        np.asarray(last_pts, dtype=np.float32),
                    )

                if np.isnan(p[0]) and np.isnan(p[1]) and np.isnan(p[2]):
                    if cur_first is not None and prev is not None:
                        first_pts.append(cur_first)
                        last_pts.append(prev)

                    cur_first = None
                    prev = None
                    continue

                if cur_first is None:
                    cur_first = p.copy()

                prev = p.copy()

        if cur_first is not None and prev is not None:
            first_pts.append(cur_first)
            last_pts.append(prev)

        return (
            np.asarray(first_pts, dtype=np.float32),
            np.asarray(last_pts, dtype=np.float32),
        )


def parse_header_labels(header: str):
    labels = [x.strip() for x in header.split(",") if x.strip()]

    if not labels:
        raise RuntimeError("Header is empty. Provide --header as comma-separated labels.")

    return labels


def enforce_atlas_labels_match_header(targets_nii_path: str, n_targets: int):
    """
    Require atlas labels to be exactly 1..N.

    This prevents silent column mislabeling.
    """
    timg = nib.load(targets_nii_path)
    tdat = np.asarray(timg.dataobj)

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
            f"  atlas unique labels min..max = {uniq[0]}..{uniq[-1]} "
            f"(count {uniq.size})\n"
            "Expected atlas labels to be exactly 1..n_targets with no gaps.\n"
            "Fix: relabel atlas to 1..N, or update header and aggregation logic."
        )


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--seed-tck", required=True)
    ap.add_argument("--seed-weights", required=True)
    ap.add_argument("--seed-mask", required=True)
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--targets-nii", required=True)
    ap.add_argument("--header", required=True)

    ap.add_argument("--out-weighted", required=True)
    ap.add_argument("--out-counts", required=True)
    ap.add_argument("--out-voxel-index", required=True)
    ap.add_argument("--out-qc", required=False, help="Output QC summary as JSON")

    ap.add_argument(
        "--keep-both-endpoints-in-seed",
        action="store_true",
        help=(
            "Keep streamlines whose BOTH endpoints fall inside the seed mask. "
        ),
    )

    args = ap.parse_args()

    target_labels = parse_header_labels(args.header)
    n_targets = len(target_labels)

    enforce_atlas_labels_match_header(args.targets_nii, n_targets)

    seed_img = nib.load(args.seed_mask)
    seed_data = np.asarray(seed_img.dataobj)
    seed_aff = seed_img.affine
    inv_seed_aff = np.linalg.inv(seed_aff)

    seed_vox = np.asarray(np.argwhere(seed_data > 0), dtype=np.int32)

    if seed_vox.size == 0:
        raise RuntimeError("Seed mask is empty. No voxels > 0.")

    I, J, K = seed_data.shape
    n_seed = seed_vox.shape[0]

    seed_hash = seed_vox[:, 0] + I * (seed_vox[:, 1] + J * seed_vox[:, 2])

    sort_idx = np.argsort(seed_hash)
    seed_hash_sorted = seed_hash[sort_idx]

    first_pts, last_pts = read_tck_endpoints(args.seed_tck)
    n_sl = first_pts.shape[0]

    try:
        weights = np.loadtxt(args.seed_weights, dtype=np.float64)
    except Exception as e:
        raise RuntimeError(
            f"Failed to read weights file: {args.seed_weights} ({type(e).__name__})"
        ) from None

    weights = np.atleast_1d(weights)

    assignments = np.loadtxt(args.assignments, dtype=np.int64)

    if assignments.ndim == 1 and assignments.size == 2:
        assignments = assignments.reshape((1, 2))

    if weights.shape[0] != n_sl:
        raise RuntimeError(
            f"Weight count mismatch: {weights.shape[0]} weights vs {n_sl} streamlines"
        )

    if assignments.shape[0] != n_sl or assignments.shape[1] < 2:
        raise RuntimeError(
            f"Assignment shape mismatch: {assignments.shape} vs expected ({n_sl}, 2)"
        )

    node_a = assignments[:, 0]
    node_b = assignments[:, 1]

    target_node = np.where(node_a > 0, node_a, node_b)
    target_node = np.where(target_node > 0, target_node, np.maximum(node_a, node_b))

    ones = np.ones((n_sl, 1), dtype=np.float64)

    first_h = np.hstack([first_pts.astype(np.float64), ones])
    last_h = np.hstack([last_pts.astype(np.float64), ones])

    first_ijk = (inv_seed_aff @ first_h.T).T[:, :3]
    last_ijk = (inv_seed_aff @ last_h.T).T[:, :3]

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

    both_in_seed = f_in_seed & l_in_seed
    if np.any(both_in_seed):
        log(
            f"WARNING: {both_in_seed.sum()} streamlines had both endpoints inside seed."
        )

    seed_ijk = np.full((n_sl, 3), -1, dtype=np.int64)

    seed_ijk[f_in_seed] = first_ijk[f_in_seed]

    use_last = (~f_in_seed) & l_in_seed
    seed_ijk[use_last] = last_ijk[use_last]

    valid = (
        (seed_ijk[:, 0] >= 0) &
        (target_node > 0) &
        (target_node <= n_targets)
    )

    if not args.keep_both_endpoints_in_seed:
        valid = valid & ~both_in_seed

    tgt_col = target_node[valid] - 1
    selected_seed_ijk = seed_ijk[valid]

    selected_hash = (
        selected_seed_ijk[:, 0] +
        I * (selected_seed_ijk[:, 1] + J * selected_seed_ijk[:, 2])
    )

    pos = np.searchsorted(seed_hash_sorted, selected_hash)

    in_range = pos < seed_hash_sorted.size
    matched = np.zeros_like(in_range, dtype=bool)
    matched[in_range] = seed_hash_sorted[pos[in_range]] == selected_hash[in_range]

    rows = np.full(selected_hash.shape, -1, dtype=np.int64)
    rows[matched] = sort_idx[pos[matched]]

    ok = rows >= 0

    rows = rows[ok]
    tgt_col = tgt_col[ok]
    selected_weights = weights[valid][ok]

    weighted = np.zeros((n_seed, n_targets), dtype=np.float64)
    counts = np.zeros((n_seed, n_targets), dtype=np.int32)

    np.add.at(weighted, (rows, tgt_col), selected_weights)
    np.add.at(counts, (rows, tgt_col), 1)

    vox_world = nib.affines.apply_affine(seed_aff, seed_vox.astype(np.float64))

    voxel_index = np.column_stack([
        np.arange(n_seed, dtype=np.int64),
        seed_vox.astype(np.int64),
        vox_world.astype(np.float64),
    ])

    os.makedirs(os.path.dirname(args.out_weighted), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_counts), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_voxel_index), exist_ok=True)

    header = ",".join(target_labels)

    np.savetxt(
        args.out_voxel_index,
        voxel_index,
        delimiter=",",
        header="row,i,j,k,x,y,z",
        comments="",
    )

    np.savetxt(
        args.out_weighted,
        weighted,
        delimiter=",",
        header=header,
        comments="",
    )

    np.savetxt(
        args.out_counts,
        counts,
        delimiter=",",
        header=header,
        comments="",
        fmt="%d",
    )
    
    if args.out_qc:
        os.makedirs(os.path.dirname(args.out_qc), exist_ok=True)

        qc = {
            "streamlines_read": int(n_sl),
            "weights_read": int(weights.shape[0]),
            "assignments_read": int(assignments.shape[0]),
            "seed_voxels": int(n_seed),
            "target_columns": int(n_targets),
            "endpoint_in_seed_first": int(f_in_seed.sum()),
            "endpoint_in_seed_last": int(l_in_seed.sum()),
            "both_endpoints_in_seed": int(both_in_seed.sum()),
            "drop_both_endpoints_in_seed": bool(not args.keep_both_endpoints_in_seed),
            "valid_streamlines_used_before_row_mapping": int(valid.sum()),
            "valid_streamlines_mapped_to_seed_rows": int(ok.sum()),
            "dropped_after_valid_check": int(n_sl - valid.sum()),
            "dropped_after_row_mapping": int(valid.sum() - ok.sum()),
        }

        with open(args.out_qc, "w") as f:
            json.dump(qc, f, indent=2)    

    log("QC summary:")
    log(f"  Streamlines read: {n_sl}")
    log(f"  Weights read: {weights.shape[0]}")
    log(f"  Assignments read: {assignments.shape[0]}")
    log(f"  Seed voxels: {n_seed}")
    log(f"  Target columns: {n_targets}")
    log(f"  Endpoint in seed, first endpoint: {f_in_seed.sum()}")
    log(f"  Endpoint in seed, last endpoint: {l_in_seed.sum()}")
    log(f"  Both endpoints in seed: {both_in_seed.sum()}")
    log(f"  Valid streamlines used before row mapping: {valid.sum()}")
    log(f"  Valid streamlines mapped to seed rows: {ok.sum()}")
    log(f"  Dropped after valid check: {n_sl - valid.sum()}")
    log(f"  Dropped after row mapping: {valid.sum() - ok.sum()}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
