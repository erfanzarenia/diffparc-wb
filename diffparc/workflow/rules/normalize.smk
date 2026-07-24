# normalize.smk -- normalize DTI metric maps (z-score / percentile) to reduce
# site effects, sampling the normalization factors from WM only (via a tissue
# dseg) to avoid confounds from the amount of CSF etc.


rule resample_tissue_dseg_to_dwi:
    input:
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="tissue",
            method="synthseg",
            suffix="dseg.nii.gz"
        ),
        ref=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="tissue",
            resliced="dwi",
            method="synthseg",
            suffix="dseg.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "c3d -int 0 {input.ref} {input.dseg} -reslice-identity -o {output.dseg}"


rule zscore_norm:
    input:
        metric=bids(
            root=root, datatype="dwi", suffix="{metric}.nii.gz", **subj_wildcards
        ),
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="tissue",
            resliced="dwi",
            method="synthseg",
            suffix="dseg.nii.gz"
        ),
    params:
        dseg_label=config["vbm"]["tissue_lut"]["WM"],
    output:
        metric=bids(
            root=root, datatype="dwi", suffix="znorm{metric}.nii.gz", **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/normalize_by_zscore.py"


rule perc_norm:
    input:
        metric=bids(
            root=root, datatype="dwi", suffix="{metric}.nii.gz", **subj_wildcards
        ),
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="tissue",
            resliced="dwi",
            method="synthseg",
            suffix="dseg.nii.gz"
        ),
    params:
        dseg_label=config["vbm"]["tissue_lut"]["WM"],
        lower_perc=2,
        upper_perc=98,
    output:
        metric=bids(
            root=root, datatype="dwi", suffix="pnorm{metric}.nii.gz", **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/normalize_by_percentile.py"
