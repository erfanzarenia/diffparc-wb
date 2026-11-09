# Target atlas spaces and the one-step NLin6 → 2009c correction

This directory is named after the registration template (`MNI152NLin2009cAsym`), but **not every
atlas in it is in that space**. Sensitivity analysis
`/imaging/ezarenia/analyses/act_atlas_sensitivity_2026-11-03/REPORT.md`, Section 4, verified this.

## Which files are in which space

| Target key(s) | File(s) | Real space | How verified |
|---|---|---|---|
| `Yeo7TianS3`, `Yeo17TianS3` | `Yeo_Tian/Yeo{7,17}_TianS3.nii.gz` | **MNI152NLin6Asym** (FreeSurfer-conformed 256³ grid) | Cortex voxels equal `Yeo/Tight_Yeo{7,17}.nii.gz`. Every Tian voxel is the world-coordinate nearest neighbour of `Tian/Tian_Subcortex_S3_7T.nii.gz` + 7 (or + 17). No warp was ever applied. |
| `SynthSegTianS3` (subcortex half) | `SynthSeg_Tian/TianS3_subcortex_dseg.nii.gz` | **MNI152NLin6Asym** | Its 57,135 Tian voxels are identical to the Tian part of `Yeo7_TianS3`. |
| `Schaefer100TianS3`, `Schaefer500TianS3` | `Schaefer_Tian/Schaefer2018_*_Tian2020_*.nii.gz` | **MNI152NLin6Asym** (FSL 182×218×182 grid) | Schaefer labels 1..N are unchanged. Tian is N+1..N+54 and comes from `Tian/Tian_S3_Schaefer_Resampled.nii.gz`, which is a world-coordinate nearest-neighbour copy of the 7T S3 file. |
| `Yeo7`, `Yeo17` (cortex only) | `Yeo/Tight_Yeo*.nii.gz` | MNI152NLin6Asym (same Yeo source) | **Not changed on this branch.** No seed uses them. Add `template_space: MNI152NLin6Asym` before you use them. |
| `harvardoxford*` | `HarvardOxford/*.nii.gz` | Unverified (2009c 193×229×193 grid) | Not checked. No seed uses them. |
| CIT168 seeds | `CIT168*/` | 2009c | Unchanged on purpose. |

The production chain treated all of these as 2009c and applied only the subject ↔ 2009c warp.
Inside the brain the NLin6 → 2009c displacement has a median of 1.45 mm (95th percentile 2.67 mm,
maximum 5.21 mm), and about 0.7 mm at the SNc.

## Naming: Tian 7T Scale III

Every `*TianS3*` key, file and `desc-` tag here contains **Tian 7T Scale III** (54 parcels),
voxel-identical to `Tian/Tian_Subcortex_S3_7T.nii.gz`. It is neither 3T S4 nor 7T S4 (62 parcels).
A pure file rename (commit 378d95a) had introduced "S4" names; the original S3 names are restored, so
outputs produced in between carry `TianS4` for the same atlas.

## How the pipeline corrects it (single resampling)

A target entry may carry `template_space: <space>`. When that differs from `template`,
`transform_targets_to_subject` (and `warp_hybrid_subcortex_to_nodegrid` for the SynthSeg hybrid)
appends `template_space_xfms[<space>]` to the subject ← template chain:

```
antsApplyTransforms -n NearestNeighbor -i <NLin6 atlas> -r <subject node grid> \
  -t [subject→2009c affine,1] <subject→2009c invwarp> <2009c←NLin6 xfm .h5>
```

The atlas is resampled once, by nearest neighbour, through the concatenated transform. **Do not**
pre-warp an atlas into a 2009c file and then resample it again: that two-step route roughly doubles
the pallidal change (+4.17 pp against +2.16 pp) and inflates pallidal volume (report, Table 8).

The transform is TemplateFlow's `tpl-MNI152NLin2009cAsym_from-MNI152NLin6Asym_mode-image_xfm.h5`
(~200 MB, md5 `1f9b29f2d5cb898cc97243fe92ee9304`). It is git-ignored. Put it at
`xfm/tpl-MNI152NLin2009cAsym_from-MNI152NLin6Asym_mode-image_xfm.h5` in this directory, or point
`template_space_xfms.MNI152NLin6Asym` at an absolute path. Two ways to get it:

```bash
cp /imaging/ezarenia/analyses/act_atlas_sensitivity_2026-11-03/atlas/downloads/templateflow/tpl-MNI152NLin2009cAsym_from-MNI152NLin6Asym_mode-image_xfm.h5 diffparc/resources/tpl-MNI152NLin2009cAsym/xfm/
```

```bash
python -c "from templateflow import api; print(api.get('MNI152NLin2009cAsym', **{'from': 'MNI152NLin6Asym'}, suffix='xfm', extension='.h5'))"
```

If a target names a `template_space` with no transform configured, the workflow fails at parse time.

Verified on sub-1782 (2026-11-06): the `Yeo7TianS3` dseg from this rule matches
`work/sub-1782/prep/dseg_corrected_onestep.nii.gz` from the sensitivity analysis on every voxel
(0 of 4,651,920 differ). The Tian parts of `Yeo17TianS3` and `SynthSegTianS3` match it exactly too.

## Re-deriving existing outputs after this change

Only the target side changes. **Do not re-run tracking or SIFT2.**

Re-run, per subject:

1. `transform_targets_to_subject`, plus `warp_hybrid_subcortex_to_nodegrid` → `build_hybrid_target`
   for `SynthSegTianS3`, giving new `anat/sub-*_desc-<target>_dseg.nii.gz`.
2. `build_seed_nodes` → `voxelwise_connectivity` (`tck2connectome`) → `connectivity_variant`.
3. `microstructure_connectome` (FA/MD/RD), if enabled.
4. Mask sweep: `sweep_build_seed_nodes` → `sweep_voxelwise_connectivity` → `sweep_connectivity_variant`.
   `sweep_build_seed_nodes` reads the dseg through `ancient()`, so a new dseg **will not** trigger it
   on its own. Force it, e.g. with `--forcerun sweep_build_seed_nodes`.
5. Enrichment sweep, if enabled: `enrich_build_selected_nodes` → `enrich_voxelwise_connectivity*`.
6. Everything downstream of the matrices: maxprob/winner maps, fingerprints, parcellations (k
   selection), group tables.

Keep as they are: registration, FOD/5TT, `transform_seed_to_subject` (CIT168), WB and ROI `tckgen`,
`final_tractogram`, SIFT2 weights and µ, and `filter_tractogram` (it depends only on the seed mask).

Practical notes:

- Dry-run first (`-n`) with `--allowed-rules` limited to the rules above, and confirm that no
  `tckgen`/`tcksift2`/`run_sift2` job appears.
- Corrected outputs keep the same filenames as the baseline ones (`desc-Yeo7TianS3`). Write them to
  a new output root, or move the baseline files aside first, so the two cannot be mixed up.
- Existing SNIT derivatives use an older layout (`tracts/`, `Yeo7TianS3` names; see the run configs in
  `batch-*/config/`). Re-deriving those needs either the matching pipeline layout or a targeted
  script that applies the steps above to the existing tractograms and weights.
