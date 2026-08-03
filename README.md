# diffparc-wb

**Diffusion-based structural connectivity of subcortical nuclei.**

`diffparc-wb` is a [Snakebids](https://github.com/khanlab/snakebids) BIDS app that builds
subject-specific **structural connectivity fingerprints** for deep-brain seed regions,
by default the midbrain dopaminergic nuclei (VTA, SNc, PBP), and maps their connectivity
onto a whole-brain cortical + subcortical target atlas. Each seed voxel is characterised by its streamline connectivity to every target parcel.

The pipeline runs anatomically-constrained (ACT) whole-brain tractography with MRtrix3,
enriches it with ROI-targeted tracking from the seed, applies SIFT2 weighting, and computes
voxelwise seed→target connectivity matrices together with microstructure-weighted connectomes
and tractography QC.

> **Origin & attribution.** This project began as a fork of
> [khanlab/diffparc-surf](https://github.com/khanlab/diffparc-surf) (Ali Khan, Khan Lab).
> It has since evolved substantially. The surface/VBM and voxel-seeding paths were removed and
> replaced by a whole-brain ACT + ROI-enrichment tractography and connectivity workflow, while
> retaining the original as its foundation. See [Acknowledgements](#acknowledgements).

---

## Contents

- [Overview](#overview)
- [How it works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Inputs](#inputs)
- [Usage](#usage)
- [Configuration](#configuration)
- [Outputs](#outputs)
- [Repository structure](#repository-structure)
- [Acknowledgements](#acknowledgements)
- [Citation](#citation)
- [License](#license)

---

## Overview

Given a BIDS dataset with a T1w image and (preprocessed) diffusion data, `diffparc-wb`:

1. Prepares the T1w image and registers it to the `MNI152NLin2009cAsym` template.
2. Warps the **seed** priors (e.g. VTA / SNc / PBP from CIT168) and a whole-brain **target**
   atlas (e.g. Yeo-7 cortical + Tian S3 subcortical) into subject space.
3. Fits MSMT-CSD fibre orientation distributions and a diffusion tensor (FA/MD/RD).
4. Runs **whole-brain, ACT-constrained tractography**, enriched with **ROI-targeted**
   tracking seeded in the deep-brain ROI.
5. Applies **SIFT2** to weight streamlines by their fit to the diffusion signal.
6. Computes **voxelwise seed→target connectivity matrices** (plus geometry- and
   microstructure-informed variants) and visual tractography QC.

The primary output is, for each seed voxel, a connectivity fingerprint across the target
parcels.

**Key features**

- Anatomically-constrained (5TT/ACT) whole-brain tractography with reproducible, chunked seeding.
- ROI enrichment to boost streamline counts from small deep-brain seeds.
- SIFT2 weighting for quantitative connectivity.
- Configurable connectivity variants (raw, SIFT2-µ scaled, volume- and length-bias corrected).
- Optional FA/MD/RD-weighted microstructure connectomes.
- Fully containerised tools and BIDS-app CLI; scales from a laptop to a SLURM cluster.

## Requirements

- **Python** 3.8–3.10
- **Snakemake** 7.x and **Snakebids** 0.7.x (installed with the app; see below)
- **Singularity / Apptainer** to run the neuroimaging tools inside the pipeline container
- The **pipeline container** — `diffparc/resources/apptainer/diffparc-deps_with_fsl.sif` — a
  **required** ~12 GB image that provides every neuroimaging tool the rules call:
  FSL, MRtrix3, ANTs, Convert3D (`c3d`/`c4d`), `greedy`, NiftyReg (`reg_aladin`),
  and FreeSurfer's SynthSeg / SynthStrip. It is an FSL-enabled build of
  [`khanlab/diffparc-deps`](https://hub.docker.com/r/khanlab/diffparc-deps). Because of its size
  it is **not** committed to the repository (it is git-ignored), so a fresh clone does not include
  it — see [Installation](#installation) for where to place it.

Whole-brain tractography is compute- and memory-intensive. A multi-core machine (or an HPC
cluster) is recommended; per-rule resource requests are defined in the config under `resources`.

## Installation

```bash
git clone https://github.com/erfanzarenia/diffparc-wb.git
cd diffparc-wb
pip install .           # or: poetry install
```

This installs the `diffparc-wb` console command. You can also run the app directly with
`python diffparc/run.py …`.

The pipeline runs its tools inside `diffparc-deps_with_fsl.sif` (see [Requirements](#requirements)).
It is not part of the clone, so obtain it and place it at its expected home in the repo:

```
diffparc/resources/apptainer/diffparc-deps_with_fsl.sif
```

The container location is read from `singularity.diffparc` in the config, or from the
`DIFFPARC_CONTAINER` environment variable — which takes precedence and is the easiest way to point
at your own copy without editing tracked config:

```bash
export DIFFPARC_CONTAINER=/abs/path/to/diffparc-deps_with_fsl.sif
```

## Inputs

A [BIDS](https://bids.neuroimaging.org/) dataset containing, per subject:

- a **T1w** anatomical image, and
- **diffusion (DWI)** data.

Diffusion preprocessing can be provided in one of two ways:

| Mode | How | When to use |
|---|---|---|
| **Import (recommended)** | `--in_snakedwi_dir <dir>` (or `--in_prepdwi_dir`) | Reuse DWI that has already been denoised, distortion/eddy corrected, bias corrected and registered to T1w by [snakedwi](https://github.com/khanlab/snakedwi). |
| **Self-contained** | *(no import flag)* | Let `diffparc-wb` run its own raw-DWI preprocessing (denoise → degibbs → topup/eddy or motion correction → masking → bias correction). |

The typical workflow is to run `snakedwi` first, then point `diffparc-wb` at its output with
`--in_snakedwi_dir`.

## Usage

`diffparc-wb` follows the standard BIDS-app / Snakebids interface:

```
diffparc-wb <bids_dir> <output_dir> participant [options]
```

**Dry run** — print what would execute without running it:

```bash
diffparc-wb /data/bids /data/derivatives/diffparc-wb participant -np
```

**Run**, importing preprocessed DWI from snakedwi, using containers and all local cores:

```bash
diffparc-wb /data/bids /data/derivatives/diffparc-wb participant \
    --in_snakedwi_dir /data/derivatives/snakedwi \
    --use-singularity --cores all
```

**Selected subjects only:**

```bash
diffparc-wb /data/bids /data/derivatives/diffparc-wb participant \
    --participant_label 001 002 --use-singularity --cores 8
```

**Generate an HTML report** after a run:

```bash
diffparc-wb /data/bids /data/derivatives/diffparc-wb participant --report
```

App-specific options (in addition to all standard Snakemake flags):

| Option | Description |
|---|---|
| `--in_snakedwi_dir <dir>` | Import preprocessed DWI from a `snakedwi` output directory. |
| `--in_prepdwi_dir <dir>` | Import preprocessed DWI from a `prepdwi` output directory. |
| `--participant_label <ids…>` | Process only these subjects (no `sub-` prefix). |
| `--exclude_participant_label <ids…>` | Process all but these subjects. |
| `--masking_method {b0_BET,b0_SyN}` | Brain-masking method for self-contained DWI preprocessing. |
| `--b0_bet_frac <f>` | BET fractional-intensity threshold for b0 masking. |
| `--no_vol` | Disable the volume/tractography workflow (anatomical steps only). |

> Set `DIFFPARC_CONTAINER` (image path) and `DIFFPARC_TMPDIR` (a **shared** scratch filesystem)
> in the environment to override the container and temp locations without editing tracked config.

## Configuration

Defaults live in [`diffparc/config/snakebids.yml`](diffparc/config/snakebids.yml). The most
commonly adjusted settings:

- **Seeds** (`select_seeds`, `seeds`): which deep-brain ROI(s) to estimate the connectivity for.
- **Targets** (`targets`, per-seed `targets`): the target atlas the seed connectivity is estimated to.
- **Tractography** (`mrtrix`, `wb_chunk_*`, `roi_enrichment`): step/angle/length, streamline
  counts (per the whole bran), ROI enrichment streamline counts (per voxel).
- **Connectivity variants** (`connectivity_variants`): which corrected matrices to emit.
- **Microstructure / QC** (`microstructure_connectomes`, `tractography_qc`).
- **Resources** (`resources`): per-rule threads/memory/time for tuning on your hardware.

See the [configuration guide](diffparc/config/README.md) for a full walkthrough of the seed
and target menus and every knob.

## Outputs

Results are written under `<output_dir>` in BIDS-derivatives layout, per subject:

```
sub-<id>/
├── anat/          preprocessed T1w, brain mask, SynthSeg, subject-space seeds & target atlas
├── dwi/           MRtrix FODs, tensor and FA/MD/RD maps, brain mask
├── warps/         subject↔template affine and (inverse) warps
├── tracts/        tractograms (.tck) plus their SIFT2 weights and µ
├── connectivity/  the primary deliverable: seed→target connectivity matrices (+ microstructure connectomes)
└── qc/            qc/tractography/ (TDI, endpoints, DEC, length histograms, tckinfo) and
                   qc/connectivity/ (streamline→parcel assignments, seed-voxel index, per-matrix QC)
```

The primary deliverables are the per-seed connectivity matrices in `connectivity/`:

- `..._label-<seed>_desc-<targets>_connectivity_matrix.csv` — raw SIFT2-weighted
  **seed-voxel × target-parcel** connectivity (plus `meas-muscaled` and, when enabled,
  volume/length bias-corrected variants).
- `..._desc-<targets>_meas-{FA,MD,RD}_connectivity_matrix.csv` — microstructure-weighted
  atlas connectomes (optional).

## Repository structure

```
diffparc/
├── run.py                  BIDS-app entry point (the `diffparc-wb` command)
├── config/snakebids.yml    default configuration (see config/README.md)
├── workflow/
│   ├── Snakefile           targets and rule includes
│   ├── rules/*.smk         pipeline stages (see workflow/README.md)
│   ├── scripts/            helper scripts (node building, connectome slicing, QC, …)
│   └── report/             Snakemake HTML-report captions
└── resources/              MNI template, tissue priors, seed & target atlases, pipeline container (apptainer/)
Dockerfile                  builds the app (Snakebids) container image
pyproject.toml              packaging and console entry points
```

## Acknowledgements

This repository is derived from
[**khanlab/diffparc-surf**](https://github.com/khanlab/diffparc-surf) by **Ali Khan** and the
[Khan Lab](https://www.khanlab.ca/), built on the
[Snakebids](https://github.com/khanlab/snakebids) framework. The original diffusion-parcellation
concept and BIDS-app scaffolding are its foundation. This fork has since been substantially
re-engineered, most notably replacing the surface/VBM and voxel-seeding pathways with a
whole-brain ACT + ROI-enrichment tractography and connectivity workflow.

It also relies on FSL, MRtrix3, ANTs, Convert3D, greedy, NiftyReg, and FreeSurfer's
SynthSeg / SynthStrip; please cite those tools where appropriate.

## Citation

If you use this pipeline, please cite the underlying tools (MRtrix3, FSL, ANTs, SynthSeg) and
acknowledge the original [diffparc-surf](https://github.com/khanlab/diffparc-surf). A citation
for this work will be added here.

## License

Released under the MIT License — see [LICENSE](LICENSE). The original diffparc-surf is likewise
MIT licensed.
