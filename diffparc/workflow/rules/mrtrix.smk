# mrtrix.smk -- MRtrix DWI modelling: DWI->MIF, optional upsampling, MSMT-CSD
# response/FOD estimation (+ mtnormalise), a single-shell CSD fallback, and
# diffusion-tensor fitting for FA/MD. Feeds the tractography/connectome rules.


def res(rule, key, default):
    return config.get("resources", {}).get(rule, {}).get(key, default)


# Diffusion-tensor (DTI) source depends on the pipeline mode:
#  * Self-contained (non-import): fit the tensor in NATIVE DWI space -- on the
#    bias-corrected DWI (dwibiascorrect_native, preproc_dwi.smk) when enabled --
#    then resample the scalar FA/MD to T1w. The native-fit rules below are gated
#    on this mode via _native_dti.
#  * Import (in_snakedwi_dir / in_prepdwi_dir): the DWI is imported already fully
#    preprocessed by the external pipeline, so the T1w-space tensor fit
#    (dwi_to_tensor -> tensor_to_metrics) is used as-is.
_native_dti = (not config.get("in_snakedwi_dir")) and (
    not config.get("in_prepdwi_dir")
)


def get_dwi_for_csd(wc):
    # nii2mif's output: already fully preprocessed (bias-corrected upstream when
    # enabled -- see _native_dti; import mode brings it in already corrected).
    return rules.nii2mif.output.dwi


rule nii2mif:
    input:
        dwi=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="preproc",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        bval=bids(
            root=root,
            suffix="dwi.bval",
            desc="preproc",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        bvec=bids(
            root=root,
            suffix="dwi.bvec",
            desc="preproc",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        mask=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        dwi=temp(
            bids(
                root=root,
                datatype="dwi",
                suffix="dwi.mif",
                **subj_wildcards,
            )
        ),
        mask=temp(
            bids(
                root=root,
                datatype="dwi",
                suffix="mask.mif",
                **subj_wildcards,
            )
        ),
    threads: lambda wc: res("nii2mif", "threads", 2)
    resources:
        mem_mb=lambda wc: res("nii2mif", "mem_mb", 4000),
        time=lambda wc: res("nii2mif", "time_min", 10)
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "mrconvert {input.dwi} {output.dwi} -fslgrad {input.bvec} {input.bval} -nthreads {threads} && "
        "mrconvert {input.mask} {output.mask} -nthreads {threads}"


rule upsample_dwi:
    input:
        dwi=get_dwi_for_csd,
    output:
        dwi=temp(
            bids(
                root=root,
                datatype="dwi",
                desc="upsampled",
                suffix="dwi.mif",
                **subj_wildcards,
            )
        ),
    params:
        vox=lambda wc: config.get("dwi_upsample_vox", 1.25),
    log:
        "logs/sub-{subject}/dwi/sub-{subject}_upsample_dwi.log",
    threads: lambda wc: res("upsample_dwi", "threads", 4)
    resources:
        mem_mb=lambda wc: res("upsample_dwi", "mem_mb", 16000),
        time=lambda wc: res("upsample_dwi", "time_min", 30),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.dwi}")"
        mkdir -p "$(dirname "{log}")"

        mrgrid -nthreads {threads} -force "{input.dwi}" regrid \
          -voxel {params.vox} "{output.dwi}" &> "{log}"
        """


rule upsample_dwi_mask:
    input:
        mask=rules.nii2mif.output.mask,
        template=rules.upsample_dwi.output.dwi,
    output:
        mask=bids(
            root=root,
            datatype="dwi",
            desc="upsampled",
            suffix="mask.mif",
            **subj_wildcards,
        ),
    log:
        "logs/sub-{subject}/dwi/sub-{subject}_upsample_dwi_mask.log",
    threads: lambda wc: res("upsample_dwi_mask", "threads", 1)
    resources:
        mem_mb=lambda wc: res("upsample_dwi_mask", "mem_mb", 4000),
        time=lambda wc: res("upsample_dwi_mask", "time_min", 10),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.mask}")"
        mkdir -p "$(dirname "{log}")"

        mrgrid -nthreads {threads} -force "{input.mask}" regrid \
          -template "{input.template}" -interp nearest "{output.mask}" &> "{log}"
        """


def get_dwi_for_fod(wildcards):
    # FOD/tracking/SIFT2 branch: upsampled DWI when enabled, else get_dwi_for_csd
    # (bias-corrected upstream in the self-contained pipeline; see _native_dti).
    # The DTI/tensor branch uses its own native fit (dwi2tensor_native, below).
    if config.get("dwi_upsample", False):
        return rules.upsample_dwi.output.dwi
    return get_dwi_for_csd(wildcards)


def get_mask_for_fod(wildcards):
    # Brain mask matched to the FOD grid (upsampled when enabled, else native).
    if config.get("dwi_upsample", False):
        return rules.upsample_dwi_mask.output.mask
    return rules.nii2mif.output.mask


# -----------------------------
# MSMT-CSD response + FOD (primary FOD path; config fod_algorithm: msmt_csd)
# -----------------------------
rule dwi2response_msmt:
    # Dhollander, T.; Mito, R.; Raffelt, D. & Connelly, A. Improved white matter response function estimation for 3-tissue constrained spherical deconvolution. Proc Intl Soc Mag Reson Med, 2019, 555
    input:
        dwi=get_dwi_for_fod,
        mask=get_mask_for_fod,
    output:
        wm_rf=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="wm",
            suffix="response.txt",
            **subj_wildcards,
        ),
        gm_rf=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="gm",
            suffix="response.txt",
            **subj_wildcards,
        ),
        csf_rf=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="csf",
            suffix="response.txt",
            **subj_wildcards,
        ),
    benchmark:
        "benchmarks/sub-{subject}/dwi/sub-{subject}_alg-msmt_desc-dwi2response.tsv"
    threads: lambda wc: res("dwi2response_msmt", "threads", 2)
    resources:
        mem_mb=lambda wc: res("dwi2response_msmt", "mem_mb", 4000),
        time=lambda wc: res("dwi2response_msmt", "time_min", 10)
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "dwi2response dhollander {input.dwi} {output.wm_rf} {output.gm_rf} {output.csf_rf}  -nthreads {threads} -mask {input.mask}"


rule dwi2fod_msmt:
    # Jeurissen, B; Tournier, J-D; Dhollander, T; Connelly, A & Sijbers, J. Multi-tissue constrained spherical deconvolution for improved analysis of multi-shell diffusion MRI data. NeuroImage, 2014, 103, 411-426
    input:
        dwi=get_dwi_for_fod,
        mask=get_mask_for_fod,
        wm_rf=rules.dwi2response_msmt.output.wm_rf,
        gm_rf=rules.dwi2response_msmt.output.gm_rf,
        csf_rf=rules.dwi2response_msmt.output.csf_rf,
    output:
        wm_fod=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="wm",
            suffix="fod.mif",
            **subj_wildcards,
        ),
        gm_fod=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="gm",
            suffix="fod.mif",
            **subj_wildcards,
        ),
        csf_fod=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="csf",
            suffix="fod.mif",
            **subj_wildcards,
        ),
    benchmark:
        "benchmarks/sub-{subject}/dwi/sub-{subject}_alg-msmt_desc-dwi2fod.tsv"
    threads: lambda wc: res("dwi2fod_msmt", "threads", 12)
    resources:
        mem_mb=lambda wc: res("dwi2fod_msmt", "mem_mb", 4000),
        time=lambda wc: res("dwi2fod_msmt", "time_min", 15)
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "dwi2fod -nthreads {threads} -mask {input.mask} msmt_csd {input.dwi} {input.wm_rf} {output.wm_fod} {input.gm_rf} {output.gm_fod} {input.csf_rf} {output.csf_fod} "


rule mtnormalise:
    # Raffelt, D.; Dhollander, T.; Tournier, J.-D.; Tabbara, R.; Smith, R. E.; Pierre, E. & Connelly, A. Bias Field Correction and Intensity Normalisation for Quantitative Analysis of Apparent Fibre Density. In Proc. ISMRM, 2017, 26, 3541
    # Dhollander, T.; Tabbara, R.; Rosnarho-Tornstrand, J.; Tournier, J.-D.; Raffelt, D. & Connelly, A. Multi-tissue log-domain intensity and inhomogeneity normalisation for quantitative apparent fibre density. In Proc. ISMRM, 2021, 29, 2472
    input:
        wm_fod=rules.dwi2fod_msmt.output.wm_fod,
        gm_fod=rules.dwi2fod_msmt.output.gm_fod,
        csf_fod=rules.dwi2fod_msmt.output.csf_fod,
        mask=get_mask_for_fod,
    output:
        wm_fod=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="wmnorm",
            suffix="fod.mif",
            **subj_wildcards,
        ),
        gm_fod=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="normalized",
            suffix="gm_fod.mif",
            **subj_wildcards,
        ),
        csf_fod=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="normalized",
            suffix="csf_fod.mif",
            **subj_wildcards,
        ),
    benchmark:
        "benchmarks/sub-{subject}/dwi/sub-{subject}_alg-msmt_desc-mtnormalise.tsv"
    threads: lambda wc: res("mtnormalise", "threads", 2)
    resources:
        mem_mb=lambda wc: res("mtnormalise", "mem_mb", 4000),
        time=lambda wc: res("mtnormalise", "time_min", 10)
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "mtnormalise -nthreads {threads} -mask {input.mask} {input.wm_fod} {output.wm_fod} {input.gm_fod} {output.gm_fod} {input.csf_fod} {output.csf_fod}"


# -----------------------------
# Single-shell CSD fallback (config fod_algorithm: csd)
# -----------------------------
rule dwi2response_csd:
    input:
        dwi=get_dwi_for_csd,
        mask=rules.nii2mif.output.mask,
    output:
        wm_rf=bids(
            root=root,
            datatype="dwi",
            alg="csd",
            desc="wm",
            suffix="response.txt",
            **subj_wildcards,
        ),
    threads: 8
    resources:
        mem_mb=32000,
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "dwi2response fa {input.dwi} {output.wm_rf} -nthreads {threads} -mask {input.mask}"


rule dwi2fod_csd:
    input:
        dwi=get_dwi_for_csd,
        mask=rules.nii2mif.output.mask,
        wm_rf=rules.dwi2response_csd.output.wm_rf,
    output:
        wm_fod=bids(
            root=root,
            datatype="dwi",
            alg="csd",
            desc="wm",
            suffix="fod.mif",
            **subj_wildcards,
        ),
    threads: 8
    resources:
        mem_mb=32000,
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "dwi2fod -nthreads {threads} -mask {input.mask} csd {input.dwi} {input.wm_rf} {output.wm_fod}  "


# -----------------------------
# Diffusion tensor -> FA/MD in T1w space (overridden by the native-space fit
# below when _native_dti; also feeds RD in microstructure_connectomes.smk)
# -----------------------------
rule dwi_to_tensor:
    input:
        dwi=get_dwi_for_csd,
        mask=rules.nii2mif.output.mask,
    output:
        tensor=bids(
            root=root,
            datatype="dwi",
            suffix="tensor.mif",
            **subj_wildcards,
        ),
    threads: 4
    resources:
        mem_mb=16000,
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "dwi2tensor -mask {input.mask} {input.dwi} {output.tensor} -nthreads {threads}"


rule tensor_to_metrics:
    input:
        tensor=rules.dwi_to_tensor.output.tensor,
        mask=rules.nii2mif.output.mask,
    output:
        fa=bids(
            root=root,
            datatype="dwi",
            suffix="FA.nii.gz",
            **subj_wildcards,
        ),
        md=bids(
            root=root,
            datatype="dwi",
            suffix="MD.nii.gz",
            **subj_wildcards,
        ),
    threads: 4
    resources:
        mem_mb=16000,
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "tensor2metric -mask {input.mask} -fa {output.fa} -adc {output.md} {input.tensor}"


# ---------------------------------------------------------------------------
# Native-space DTI (self-contained pipeline only; gated by _native_dti above).
#
# Fit the tensor in NATIVE DWI space, derive the rotation-invariant scalars
# (FA, MD) there, then resample ONLY those scalars to T1w. Fitting before
# resampling avoids interpolating the 4D directional signal, and because FA/MD
# are rotation invariant the scalar resample needs no tensor reorientation.
#
# Outputs reuse the FA.nii.gz / MD.nii.gz paths of the T1w-space fit
# (tensor_to_metrics); the ruleorder below makes this native-fit resample
# authoritative for them. Import mode keeps the T1w-space fit unchanged.
# ---------------------------------------------------------------------------

if _native_dti:

    def get_native_dwi_for_tensor(wildcards):
        # Same selection as reg_dwi_to_t1.smk's get_native_dwi_desc(): one
        # native DWI (bias-corrected when enabled, via rule
        # dwibiascorrect_native in preproc_dwi.smk) feeds both this tensor
        # fit and, via registration/resample, the FOD branch.
        desc = "biascorr" if config.get("dwi_biascorrect", True) else "preproc"
        prefix = bids(
            root=root, suffix="dwi", desc=desc, datatype="dwi", **subj_wildcards
        )
        return {
            "dwi": f"{prefix}.nii.gz",
            "bvec": f"{prefix}.bvec",
            "bval": f"{prefix}.bval",
        }

    rule dwi2tensor_native:
        input:
            unpack(get_native_dwi_for_tensor),
            mask=get_dwi_mask(),
        output:
            tensor=bids(
                root=root, datatype="dwi", proc="native",
                suffix="tensor.mif", **subj_wildcards
            ),
        threads: lambda wc: res("dwi2tensor_native", "threads", 4)
        resources:
            mem_mb=lambda wc: res("dwi2tensor_native", "mem_mb", 16000),
            time=lambda wc: res("dwi2tensor_native", "time_min", 10),
        group:
            "subj"
        container:
            config["singularity"]["diffparc"]
        shell:
            "dwi2tensor -fslgrad {input.bvec} {input.bval} "
            "-mask {input.mask} {input.dwi} {output.tensor} "
            "-nthreads {threads}"

    rule tensor2metric_native:
        input:
            tensor=rules.dwi2tensor_native.output.tensor,
            mask=get_dwi_mask(),
        output:
            fa=bids(
                root=root, datatype="dwi", proc="native",
                suffix="FA.nii.gz", **subj_wildcards
            ),
            md=bids(
                root=root, datatype="dwi", proc="native",
                suffix="MD.nii.gz", **subj_wildcards
            ),
        threads: lambda wc: res("tensor2metric_native", "threads", 2)
        resources:
            mem_mb=lambda wc: res("tensor2metric_native", "mem_mb", 8000),
            time=lambda wc: res("tensor2metric_native", "time_min", 10),
        group:
            "subj"
        container:
            config["singularity"]["diffparc"]
        shell:
            "tensor2metric -mask {input.mask} -fa {output.fa} -adc {output.md} "
            "{input.tensor} -nthreads {threads}"

    rule resample_dti_metric_to_t1w:
        # Scalar (rotation-invariant) resample of the native-fit metric onto the
        # T1w grid, via the dwi->T1w rigid transform. Trilinear keeps FA within
        # [0, 1] and MD non-negative (no overshoot); the reference grid is the
        # existing T1w brain mask, so the output grid matches the previous
        # (T1w-fit) FA/MD exactly.
        wildcard_constraints:
            metric="FA|MD",
        input:
            metric=bids(
                root=root, datatype="dwi", proc="native",
                suffix="{metric}.nii.gz", **subj_wildcards
            ),
            ref=bids(
                root=root, suffix="mask.nii.gz", desc="brain", space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi", **subj_wildcards
            ),
            xfm=bids(
                root=root, suffix="xfm.txt", from_="dwi", to="T1w",
                type_="itk", datatype="dwi", **subj_wildcards
            ),
        output:
            metric=bids(
                root=root, datatype="dwi",
                suffix="{metric}.nii.gz", **subj_wildcards
            ),
        params:
            interpolation="Linear",
        threads: lambda wc: res("resample_dti_metric_to_t1w", "threads", 1)
        resources:
            mem_mb=lambda wc: res("resample_dti_metric_to_t1w", "mem_mb", 4000),
            time=lambda wc: res("resample_dti_metric_to_t1w", "time_min", 10),
        group:
            "subj"
        container:
            config["singularity"]["diffparc"]
        shell:
            "antsApplyTransforms -d 3 --input-image-type 0 "
            "--input {input.metric} --reference-image {input.ref} "
            "--transform {input.xfm} --interpolation {params.interpolation} "
            "--output {output.metric} --verbose"

    # FA/MD now come from the native-space fit (resampled), not the T1w-space
    # tensor fit. Same output paths -> resolve the overlap in favour of the
    # native-fit resample.
    ruleorder: resample_dti_metric_to_t1w > tensor_to_metrics


def get_fod_for_tracking(wildcards):
    if config["fod_algorithm"] == "csd":
        return (
            bids(
                root=root,
                datatype="dwi",
                alg="csd",
                desc="wm",
                suffix="fod.mif",
                **subj_wildcards,
            ),
        )
    elif config["fod_algorithm"] == "msmt_csd":
        return (
            bids(
                root=root,
                datatype="dwi",
                desc="wmnorm",
                alg="msmt",
                suffix="fod.mif",
                **subj_wildcards,
            ),
        )
