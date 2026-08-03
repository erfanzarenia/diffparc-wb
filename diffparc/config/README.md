# Configuration

All defaults live in [`snakebids.yml`](snakebids.yml). Any key can be overridden on the command
line (`--key value`), and a few are exposed as named CLI options via the `parse_args` section
(see the [main README](../../README.md#usage)). This guide covers the settings you are most
likely to change; the file itself is organised into labelled sections and commented inline.

## Execution environment

| Key | Purpose |
|---|---|
| `singularity.diffparc` | Absolute path to the **primary pipeline container** (`diffparc-deps_with_fsl.sif`), whose home is `resources/apptainer/`. Override per machine with `$DIFFPARC_CONTAINER` instead of editing this. |
| `tmp_dir` | Scratch directory for temporary tractograms. Override with `$DIFFPARC_TMPDIR`. **Must be a shared filesystem** on a cluster (not node-local), because grouped jobs share these files. |

## DWI source & preprocessing

| Key | Purpose |
|---|---|
| `in_snakedwi_dir` / `in_prepdwi_dir` | Import already-preprocessed DWI (recommended). Set via the matching CLI flags. |
| `dwi_biascorrect` | Native-space bias correction in the self-contained path. Keep `False` when snakedwi already did it. |
| `masking_method` | `b0_BET` or `b0_SyN`: brain masking for the self-contained path. |
| `resample_dwi` | Target grid for DWI in T1w space (default: single resample to 1.25 mm). Must match snakedwi. |

## Registration template

`template`, `template_t1w`, `template_mask`, `template_b0` define the reference space
(`MNI152NLin2009cAsym`, shipped in [`../resources`](../resources)). `ants` holds the
registration parameters.

## Seeds and targets

This is the core of what the pipeline computes: connectivity **from** each seed **to** the
parcels of a target atlas.

**`select_seeds`** lists the seed(s) to run. Each entry under **`seeds`** defines a deep-brain ROI:

| Field | Meaning |
|---|---|
| `template_probseg` | Template-space probability map for the seed (hemisphere-split). |
| `probseg_threshold` | Threshold used to binarise the seed mask. |
| `use_synthseg` | If `True`, derive the seed from SynthSeg + shape injection instead of the template prior. |
| `targets` | Which target atlas(es) to connect this seed to. |

Bundled seeds (midbrain dopaminergic nuclei, from CIT168):

| Seed | Region |
|---|---|
| `vtasncpbp` *(default)* | VTA + SNc + PBP (full dopaminergic midbrain) |
| `vtasnc` | VTA + SNc |
| `snc` | SNc only (nigrostriatal) |
| `vtapbp` | VTA + PBP (mesolimbic) |

**`targets`** is a menu of whole-brain atlases; each entry has a `template_dseg` and an ordered
`labels` list (the labels become the columns of the connectivity matrix). Bundled options include
`Yeo7TianS3` *(default)* and `Yeo17TianS3` (Yeo cortical networks + Tian S3 subcortical),
`Schaefer100TianS3` / `Schaefer500TianS3`, `Yeo7` / `Yeo17`, and several `harvardoxford*` variants.
To parcellate a seed against a different atlas, set that seed's `targets` accordingly.

## Tractography

| Key | Purpose |
|---|---|
| `fod_algorithm` | `msmt_csd` (default) or single-shell `csd`. |
| `mrtrix` | iFOD2 parameters: `step`, `angle`, `minlength`, `maxlength`, `cutoff`, and `target_search_radius` (mm) for endpoint→parcel assignment. |
| `wb_chunk_count`, `wb_chunk_streamlines` | Whole-brain tractography is generated in this many chunks (distinct RNG seeds) of this many streamlines each, then merged. |
| `mrtrix_rng_seed` | Base RNG seed for reproducibility. |
| `roi_enrichment` | `enabled` and `seeds_per_voxel` for ROI-targeted enrichment tracking. |

## Connectivity outputs

`connectivity_variants` toggles which matrices are written (all additive; `raw` and `mu_scaled`
are always produced):

| Variant | tck2connectome scaling |
|---|---|
| `volume_bias_corrected` | `-scale_invnodevol` |
| `length_bias_corrected` | `-scale_length` |
| `volume_length_bias_corrected` | both |
| `volume_length_bias_corrected_mu_scaled` | both, then × SIFT2 µ |

`microstructure_connectomes` (+ `_targets`, `_metrics`) enables FA/MD/RD-weighted atlas
connectomes. `tractography_qc` (+ `_nsubsample`, `_seed`) controls the visual QC outputs.

## Experiments (config-sweep branch)

Two additive parameter sweeps, each fully config-gated; see the
[README](../../README.md#experiments) for what they are for.

**`mask_sweep`** (seed-mask sensitivity sweep):

| Key | Purpose |
|---|---|
| `enabled` | Turn the seed-mask sweep on/off. |
| `thresholds` | Probseg thresholds to sweep (e.g. `0.1, 0.2, 0.3`). |
| `mask_sources` | Seed-mask sources to compare (e.g. `native`, `brainsteminjected`). |

**`enrichment_sweep`** (ROI-enrichment strength sweep):

| Key | Purpose |
|---|---|
| `enabled` | Turn the enrichment sweep on/off. |
| `seed`, `hemis` | Which ROI seed (per-hemi mask must exist) and hemispheres to sweep. |
| `levels` | Enrichment levels (seeds-per-voxel-equivalent); the max is generated, the rest are reproducible nested subsamples of it. |
| `subsample_seed` | RNG seed for the nested subsampling. |
| `targets` | Atlas(es) for the per-condition connectivity fingerprint. |
| `volume_bias_correction` | Also emit a volume-bias-corrected fingerprint per condition. |
| `sift2_propagated` | Also emit a second SIFT2 branch (weights propagated from the maximal tractogram) alongside the default re-fit branch. |

This branch also defines extra seeds (`snc`, `vtapbp`, …) beyond the core `vtasncpbp`.

## Resources

The `resources` block sets per-rule `threads` / `mem_mb` / `time_min` for the tractography-heavy
rules. Tune these to your hardware; on a cluster they feed directly into the scheduler's job
requests. Only rules that scale with tractography parameters are listed; everything else uses
sensible built-in defaults.
