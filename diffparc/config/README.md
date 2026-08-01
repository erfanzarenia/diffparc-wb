# Configuration

All defaults live in [`snakebids.yml`](snakebids.yml). Any key can be overridden on the command
line (`--key value`), and a few are exposed as named CLI options via the `parse_args` section
(see the [main README](../../README.md#usage)). This guide covers the settings you are most
likely to change; the file itself is organised into labelled sections and commented inline.

## Execution environment

| Key | Purpose |
|---|---|
| `singularity.diffparc` | Absolute path to the **primary pipeline container** — `diffparc-deps_with_fsl.sif`, whose home is `resources/apptainer/`. Override per machine with `$DIFFPARC_CONTAINER` instead of editing this. |
| `tmp_dir` | Scratch directory for temporary tractograms. Override with `$DIFFPARC_TMPDIR`. **Must be a shared filesystem** on a cluster (not node-local), because grouped jobs share these files. |

## DWI source & preprocessing

| Key | Purpose |
|---|---|
| `in_snakedwi_dir` / `in_prepdwi_dir` | Import already-preprocessed DWI (recommended). Set via the matching CLI flags. |
| `dwi_biascorrect` | Native-space bias correction in the self-contained path. Keep `False` when snakedwi already did it. |
| `masking_method` | `b0_BET` or `b0_SyN` — brain masking for the self-contained path. |
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
| `vtasncpbp` *(default)* | VTA + SNc + PBP — full dopaminergic midbrain |
| `vtasnc` | VTA + SNc |
| `snc` | SNc only (nigrostriatal) |
| `vtapbp` | VTA + PBP (mesolimbic) |

**`targets`** is a menu of whole-brain atlases; each entry has an ordered `labels` list (the labels
become the columns of the connectivity matrix). A seed is parcellated against **every** atlas in its
`targets:` list, so listing several produces several matrices from the *same* tractography: the
expensive tracking and SIFT2 stages are shared, and each extra atlas adds only a cheap warp plus
`tck2connectome`. Select all, or any subset, by editing that list.

This branch (`feature/target-sweep`) ships four target atlases, each pairing a cortical
parcellation with the **Tian S4** subcortex (54 parcels, 27 per hemisphere):

| Target | Cortex | Cortex source | Matrix columns |
|---|---|---|---|
| `Yeo7TianS4` *(default)* | Yeo 7 networks | template (warped from MNI) | 61 |
| `Yeo17TianS4` | Yeo 17 networks | template (warped from MNI) | 71 |
| `Schaefer100TianS4` | Schaefer 100 (17-network) | template (warped from MNI) | 154 |
| `SynthSegTianS4` | SynthSeg Desikan-Killiany (68) | **native subject space** (reslice, no warp) | 122 |

`SynthSegTianS4` is a **hybrid**: the cortex is a per-subject SynthSeg `--parc` segmentation (68
Desikan-Killiany parcels) resliced onto the tracking grid with no inter-subject warp, while the Tian
S4 subcortex is warped from MNI as usual (assembled by
[`target_synthseg_tian.smk`](../workflow/rules/target_synthseg_tian.smk)). It is modular: the cortical
columns are pinned by `resources/tpl-MNI152NLin2009cAsym/SynthSeg_Tian/synthseg_cortex_labels.tsv` (so
matrices stay aligned across subjects, with any absent parcel becoming a zero column), and pointing
`cortex_lut` / `subcortex_dseg` at other files swaps either half. Because its cortex is anatomical
(gyral), its columns are **not** comparable to the Yeo/Schaefer functional ones. Also bundled:
cortex-only `Yeo7` / `Yeo17` and several `harvardoxford*` variants.

> **Tian scale.** These subcortical parcels are Tian **Scale IV (S4)** (54 parcels), despite the
> upstream files historically carrying an "S3" suffix. The atlases, label lists and target keys were
> corrected to `...S4` so the `desc-` tag in every output filename reflects the real scale. The
> unused, incomplete `Schaefer500TianS3` entry is left untouched and flagged inline in the config.

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

## Resources

The `resources` block sets per-rule `threads` / `mem_mb` / `time_min` for the tractography-heavy
rules. Tune these to your hardware; on a cluster they feed directly into the scheduler's job
requests. Only rules that scale with tractography parameters are listed — everything else uses
sensible built-in defaults.
