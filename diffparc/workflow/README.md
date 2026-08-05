# Workflow

The pipeline is a set of Snakemake rule modules under [`rules/`](rules/), wired together by the
[`Snakefile`](Snakefile). The `Snakefile` defines the final targets (in `rule all`) and includes
the rule files in dependency order; the anatomical stages always run, and the diffusion stages run
unless `--no_vol` / `anat_only` is set.

## Pipeline stages

| Stage | Rule module | What it does |
|---|---|---|
| Shared helpers | `common.smk`, `container_check.smk` | Path helpers; container sanity check. |
| T1w preprocessing | `preproc_t1.smk` | Import T1w, SynthStrip brain mask, N4, masked T1w. |
| T1w → template | `reg_t1_to_template.smk` | Deformable (greedy) registration to `MNI152NLin2009cAsym`. |
| Segmentation | `synthseg.smk` | SynthSeg on subject and template; seed probability maps via shape injection. |
| Seeds & targets | `prop_seeds_targets.smk` | Warp template seeds and (template) target atlases into subject space; binarise/trim seed masks. |
| Hybrid target | `target_synthseg_tian.smk` | Assemble the hybrid `SynthSegTianS4` target: reslice native SynthSeg cortex + warp the Tian S4 subcortex + merge into one contiguous dseg. Only for `kind: hybrid_synthseg_tian` targets; template targets skip it. |
| DWI preprocessing | `preproc_dwi.smk` | Self-contained raw-DWI path: denoise → degibbs → topup/eddy or motion correction → masking → native bias correction. *(Skipped in import mode.)* |
| Motion correction | `motioncorr.smk` | Volume-wise rigid motion correction (used when eddy is off). |
| Brain masking | `masking_bet_from-b0.smk`, `masking_b0_to_template.smk` | DWI brain mask via BET-on-b0 or registration to the template. |
| DWI → T1w | `reg_dwi_to_t1.smk` | Import preprocessed DWI, or rigidly register + resample the native DWI onto the T1w grid. |
| MRtrix modelling | `mrtrix.smk` | DWI→MIF, MSMT-CSD response/FOD (+ mtnormalise), CSD fallback, tensor → FA/MD/RD. |
| ACT priors | `wb_prep.smk` | 5TT image and seed-mask preparation for anatomically-constrained tracking. |
| Whole-brain tractography | `mrtrix_wb.smk` | Chunked, ACT-constrained iFOD2 whole-brain tracking, merged. |
| ROI enrichment | `mrtrix_roi_enrichment.smk` | ROI-seeded enrichment tracking; merged with the whole-brain tractogram into the final tractogram. |
| SIFT2 | `sift2.smk` | SIFT2 weighting of the final tractogram (+ global µ). |
| Connectivity | `connectivity.smk` | Seed-filter streamlines, build per-voxel nodes, `tck2connectome` → seed→target matrices (+ variants) and tractography QC. |
| Microstructure | `microstructure_connectomes.smk` | Optional FA/MD/RD-weighted atlas connectomes from the final tractogram. |

## Conventions

- **Output paths** are built with Snakebids' `bids()` helper, so filenames follow BIDS entity
  conventions (`sub-…_hemi-…_label-…_desc-…_…`).
- **Output layout** — per subject, files are grouped by kind: `anat/`, `dwi/`, `warps/` (inputs and
  intermediates); `tracts/` (tractograms plus their SIFT2 weights); `connectivity/` (the
  seed→target and microstructure matrices — the deliverables); and `qc/` (the global SIFT2 µ), split into
  `qc/tractography/` (TDI, endpoints, DEC, lengths, tckinfo) and `qc/connectivity/` (assignments,
  seed-voxel index, per-matrix QC).
- **Containers** — every rule that calls a neuroimaging tool declares the `singularity` container;
  run with `--use-singularity`.
- **Grouping** — most rules use `group: "subj"` so a subject's work can be bundled into one
  cluster job (see `--group-components`).
- **Resources** — tractography-heavy rules read `threads`/`mem_mb`/`time` from the config
  `resources` block via a small `res()` helper.
- **Scripts** — non-trivial logic lives in [`scripts/`](scripts/) (voxel-node building,
  connectome slicing/labelling, tractography QC, hybrid-target merging, DWI preprocessing helpers).
- **Targets are generic** — every atlas in a seed's `targets:` list flows through the same
  `transform`/hybrid → `build_seed_nodes` → `tck2connectome` path, so adding one reuses the shared
  tractography and SIFT2 (see the [config guide](../config/README.md#seeds-and-targets)).

See the [configuration guide](../config/README.md) for the knobs that drive these rules.
