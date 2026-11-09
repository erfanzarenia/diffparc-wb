#!/usr/bin/env python
"""Ventral posterior (VP) thalamus mask for post-SIFT2 streamline exclusion.

On the grid of the FIRST thalamus partial-volume maps (the 5TT grid), every voxel with
FIRST thalamus partial volume > 0 takes the nearest Tian thalamic label (THA-*) of its own
hemisphere from the subject-space target atlas; the mask is the voxels whose label is
THA-VPm or THA-VPl. The fill covers the FIRST thalamus rim that the Tian atlas leaves
unlabelled, which would otherwise stay as terminal sub-cortical grey matter.
"""
import argparse

import nibabel as nib
import numpy as np
from scipy import ndimage

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--thal-pv-l", required=True)
p.add_argument("--thal-pv-r", required=True)
p.add_argument("--dseg", required=True, help="subject-space target atlas (labels 1..N)")
p.add_argument("--labels", required=True, help="comma-separated atlas label names, in label order")
p.add_argument("--out", required=True)
a = p.parse_args()

names = a.labels.split(",")
thal = {h: [i + 1 for i, n in enumerate(names) if n.startswith("THA-") and n.endswith(f"-{h.lower()}h")] for h in "LR"}
vp = [i + 1 for i, n in enumerate(names) if n.split("-")[:2] in (["THA", "VPm"], ["THA", "VPl"])]
assert all(thal[h] for h in "LR") and len(vp) == 4, "atlas needs Tian THA-* labels incl. THA-VPm/VPl in both hemispheres"

pv_img = {"L": nib.load(a.thal_pv_l), "R": nib.load(a.thal_pv_r)}
ref = pv_img["L"]
assert np.allclose(ref.affine, pv_img["R"].affine) and ref.shape == pv_img["R"].shape
shp = ref.shape[:3]

# atlas label at each grid voxel centre (nearest neighbour in world coordinates)
at = nib.load(a.dseg)
al = np.rint(np.asarray(at.dataobj)).astype(np.int32)
ijk = np.indices(shp).reshape(3, -1).T
xyz = ijk @ ref.affine[:3, :3].T + ref.affine[:3, 3]
inv = np.linalg.inv(at.affine)
v = np.rint(xyz @ inv[:3, :3].T + inv[:3, 3]).astype(int)
ok = np.all((v >= 0) & (v < np.array(al.shape[:3])), axis=1)
atl = np.zeros(len(v), np.int32)
atl[ok] = al[v[ok, 0], v[ok, 1], v[ok, 2]]
atl = atl.reshape(shp)

fill = np.zeros(shp, np.int32)
for h in "LR":
    th = np.asarray(pv_img[h].dataobj, dtype=np.float32) > 0
    lab = np.isin(atl, thal[h])
    _, idx = ndimage.distance_transform_edt(~lab, return_indices=True)
    nearest = atl[idx[0], idx[1], idx[2]]
    fill[th] = nearest[th]
    print(f"{h}: FIRST thalamus voxels {int(th.sum())}, unlabelled and filled {100 * (th & ~lab).sum() / max(th.sum(), 1):.1f}%")

mask = np.isin(fill, vp).astype(np.uint8)
vox = float(np.prod(ref.header.get_zooms()[:3]))
print(f"VP mask: {int(mask.sum())} voxels, {mask.sum() * vox:.0f} mm3")
out = nib.Nifti1Image(mask, ref.affine)
out.set_data_dtype(np.uint8)
out.to_filename(a.out)
