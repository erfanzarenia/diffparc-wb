# Resources

Template and atlas files used by the pipeline. Paths here are referenced from
[`../config/snakebids.yml`](../config/snakebids.yml) (seed `template_probseg`, target
`template_dseg`, and the `template_*` keys).

## Contents

| Path | Contents |
|---|---|
| `apptainer/diffparc-deps_with_fsl.sif` | **Primary pipeline container** (required, ~12 GB, git-ignored): the FSL-enabled diffparc-deps image every rule runs inside. See [Installation](../../README.md#installation) for how to obtain and place it. |
| `tpl-MNI152NLin2009cAsym/` | Reference template: T1w, brain mask, GM/WM/CSF tissue priors, and a DWI-space b0. |
| `tpl-MNI152NLin2009cAsym/CIT168_prob/` | **Seed** probability maps (VTA / SNc / PBP and their combinations), hemisphere-split. |
| `tpl-MNI152NLin2009cAsym/{Yeo,Tian,Yeo_Tian,Schaefer,Schaefer_Tian,HarvardOxford}/` | **Target** atlases (cortical networks and subcortical parcellations) and their label lists. |
| `synthseg_simple_labels.tsv`, `synthseg_cortparc_labels.tsv` | SynthSeg label lookup tables. |

Seeds are defined in the template space and warped into each subject via the subject↔template
transforms (see [`workflow/rules/prop_seeds_targets.smk`](../workflow/rules/prop_seeds_targets.smk)).

## Provenance

These atlases are redistributed from their original sources: CIT168 (subcortical/midbrain
nuclei), Yeo (cortical networks), Tian (subcortical), Schaefer (cortical), and Harvard-Oxford,
in the `MNI152NLin2009cAsym` space. Please cite the corresponding atlas publications when using
the results.
