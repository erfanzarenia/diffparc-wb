#!/usr/bin/env python3
"""Build a combined node-parcellation image for MRtrix ``tck2connectome``.

Every voxel of the (resliced) seed ROI becomes its *own* node, and the target
atlas labels are shifted to occupy node ids *above* the seed voxels. Running
``tck2connectome`` on this image therefore yields a connectome whose
[seed-voxel x atlas-target] block is exactly the voxel-wise connectivity
fingerprint we want (see ``slice_connectome_block.py``).

Node-id layout (atlas = 1..n_targets, seed = n_targets+1..n_targets+n_seed):

  * Seed voxels are written last and take priority, so every seed id is always
    present in the image. The maximum node id is therefore fixed at
    ``n_targets + n_seed`` regardless of how the streamlines fall, which means
    the connectome matrix has a fixed, predictable size and the atlas columns
    stay aligned to their label index even if an atlas parcel happens to be
    missing in native space (it just becomes an all-zero column).

The seed voxels are enumerated in C-order (``np.argwhere``), i.e. the same order
NumPy uses for boolean-mask assignment, so downstream code that does
``image[seed_mask > 0] = matrix[:, t]`` lines up row-for-row.
"""

import argparse
import csv
import os
import sys

import nibabel as nib
import nibabel.affines as naff
import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed-mask", required=True,
                    help="Binary seed mask, already resliced onto the node grid.")
    ap.add_argument("--targets-nii", required=True,
                    help="Target atlas dseg on the SAME grid as --seed-mask.")
    ap.add_argument("--n-targets", type=int, required=True,
                    help="Number of atlas targets (length of the label header).")
    ap.add_argument("--out-nodes", required=True,
                    help="Output integer node-parcellation image.")
    ap.add_argument("--out-voxel-index", required=True,
                    help="Output CSV mapping matrix rows -> seed voxel ijk/xyz.")
    args = ap.parse_args()

    seed_img = nib.load(args.seed_mask)
    atlas_img = nib.load(args.targets_nii)

    seed = np.asarray(seed_img.dataobj)
    atlas = np.rint(np.asarray(atlas_img.dataobj)).astype(np.int64)

    if seed.shape != atlas.shape or not np.allclose(
        seed_img.affine, atlas_img.affine, atol=1e-3
    ):
        sys.exit(
            "seed mask and atlas are not on the same voxel grid; "
            "reslice the seed onto the atlas grid first."
        )

    # C-order enumeration matches NumPy boolean-mask ordering used downstream.
    seed_vox = np.argwhere(seed > 0)
    n_seed = int(seed_vox.shape[0])
    if n_seed == 0:
        sys.exit("seed mask is empty (no voxels > 0).")

    nodes = np.zeros(seed.shape, dtype=np.int32)

    # Atlas occupies ids 1..n_targets (assumes labels are within 1..n_targets).
    amask = atlas > 0
    nodes[amask] = atlas[amask].astype(np.int32)

    # Seed voxels occupy ids n_targets+1..n_targets+n_seed and take priority.
    node_ids = np.arange(1, n_seed + 1, dtype=np.int32) + args.n_targets
    nodes[seed_vox[:, 0], seed_vox[:, 1], seed_vox[:, 2]] = node_ids

    # Fresh header from the affine (avoids inheriting scl_slope/dtype from the
    # uchar seed mask, which would corrupt the integer node labels).
    out_img = nib.Nifti1Image(nodes, seed_img.affine)
    out_img.set_data_dtype(np.int32)
    os.makedirs(os.path.dirname(args.out_nodes) or ".", exist_ok=True)
    out_img.to_filename(args.out_nodes)

    world = naff.apply_affine(seed_img.affine, seed_vox.astype(float))
    os.makedirs(os.path.dirname(args.out_voxel_index) or ".", exist_ok=True)
    with open(args.out_voxel_index, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["row", "i", "j", "k", "x", "y", "z", "node_id"])
        for r in range(n_seed):
            i, j, k = (int(v) for v in seed_vox[r])
            x, y, z = (float(v) for v in world[r])
            w.writerow([r, i, j, k, f"{x:.4f}", f"{y:.4f}", f"{z:.4f}",
                        int(node_ids[r])])

    print(f"build_voxel_nodes: n_seed={n_seed} n_targets={args.n_targets} "
          f"max_node_id={args.n_targets + n_seed}", file=sys.stderr)


if __name__ == "__main__":
    main()
