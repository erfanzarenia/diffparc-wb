# mask_sweep.smk
# Experimental mask sweep.
# Compares multiple seed sources:
#   1. native: existing template-to-subject seed probseg
#   2. brainsteminjected: template SN/VTA prior warped along a template-brainstem to subject-brainstem warp
# Reuses existing final tractogram + SIFT2 weights.
# Does NOT overwrite main pipeline masks.

import os

SUBJECTS = sorted(set(inputs.input_zip_lists["T1w"]["subject"]))

VOXEL_AGG_SCRIPT = os.path.join(workflow.basedir, "scripts", "voxel_target_aggregate.py")

SWEEP_SEEDS = ["vtasnc", "snc", "vtapbp", "vtasncpbp"]
HEMIS = ["L", "R"]
TARGETS = ["Yeo7TianS3"]

MASK_SWEEP_ENABLED = config.get("mask_sweep", {}).get("enabled", False)

MASK_SOURCES = config.get("mask_sweep", {}).get(
    "mask_sources",
    ["native", "brainsteminjected"]
)

THRESHOLDS = [
    str(t) for t in config.get("mask_sweep", {}).get("thresholds", [])
]
THRESHOLD_TAGS = {
    str(t): str(t).replace(".", "p") for t in config["mask_sweep"]["thresholds"]
}


def threshold_tag_list():
    return [THRESHOLD_TAGS[t] for t in THRESHOLDS]


def template_prefix_str():
    return get_template_prefix(
        root=root,
        subj_wildcards=subj_wildcards,
        template=config["template"],
    )


def template_seed_probseg_path(wc):
    p = config["seeds"][wc.seed]["template_probseg"]

    if "{hemi}" in p:
        p = p.format(hemi=wc.hemi)

    return os.path.join(workflow.basedir, "..", p)


def sweep_seed_probseg(wc):
    if wc.mask_source == "native":
        return (
            f"sub-{wc.subject}/anat/"
            f"sub-{wc.subject}_hemi-{wc.hemi}_label-{wc.seed}_probseg.nii.gz"
        )

    if wc.mask_source == "brainsteminjected":
        return (
            f"sub-{wc.subject}/anat/mask_sweep/{wc.seed}/brainsteminjected/"
            f"sub-{wc.subject}_hemi-{wc.hemi}_label-{wc.seed}"
            f"_desc-brainsteminjected_probseg.nii.gz"
        )

    raise ValueError(f"Unknown mask_source: {wc.mask_source}")


rule all_mask_sweep:
    input:
        expand(
            "sub-{subject}/tracts/mask_sweep/{seed}/{mask_source}/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_connectivity_matrix.csv",
            subject=SUBJECTS,
            mask_source=MASK_SOURCES,
            hemi=HEMIS,
            seed=SWEEP_SEEDS,
            targets=TARGETS,
            thr_tag=threshold_tag_list(),
        ) if MASK_SWEEP_ENABLED else []
        
        
# -----------------------------
# Brainstem Injection seed source
# -----------------------------

rule sweep_extract_subject_brainstem:
    input:
        dseg="sub-{subject}/anat/sub-{subject}_desc-synthseg_dseg.nii.gz",
    output:
        mask=(
            "sub-{subject}/anat/mask_sweep/brainsteminjected/"
            "sub-{subject}_label-brainstem_desc-synthseg_mask.nii.gz"
        ),
    log:
        "logs/sub-{subject}/mask_sweep/brainsteminjected/sub-{subject}_extract_subject_brainstem.log",
    benchmark:
        "benchmarks/sub-{subject}/mask_sweep/brainsteminjected/sub-{subject}_extract_subject_brainstem.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        time=10
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.mask}")"
        mkdir -p "$(dirname "{log}")"

        c3d "{input.dseg}" \
          -retain-labels 16 \
          -binarize \
          -o "{output.mask}" \
          &> "{log}"
        """


rule sweep_extract_template_brainstem:
    input:
        dseg=get_template_prefix(
            root=root,
            subj_wildcards=subj_wildcards,
            template=config["template"],
        ) + "_desc-synthseg_dseg.nii.gz",
    output:
        mask=get_template_prefix(
            root=root,
            subj_wildcards=subj_wildcards,
            template=config["template"],
        ) + "_label-brainstem_desc-synthseg_mask.nii.gz",
    log:
        "logs/sub-{subject}/mask_sweep/brainsteminjected/sub-{subject}_extract_template_brainstem.log",
    benchmark:
        "benchmarks/sub-{subject}/mask_sweep/brainsteminjected/sub-{subject}_extract_template_brainstem.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        time=10
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.mask}")"
        mkdir -p "$(dirname "{log}")"

        c3d "{input.dseg}" \
          -retain-labels 16 \
          -binarize \
          -o "{output.mask}" \
          &> "{log}"
        """


rule sweep_brainstem_proxy_to_subject:
    input:
        fixed=rules.sweep_extract_subject_brainstem.output.mask,
        moving=rules.sweep_extract_template_brainstem.output.mask,
    output:
        affine=(
            "sub-{subject}/warps/mask_sweep/brainsteminjected/"
            "sub-{subject}_label-brainstem_desc-proxy_from-template_to-subject_affine.txt"
        ),
        warp=(
            "sub-{subject}/warps/mask_sweep/brainsteminjected/"
            "sub-{subject}_label-brainstem_desc-proxy_from-template_to-subject_warp.nii.gz"
        ),
        invwarp=(
            "sub-{subject}/warps/mask_sweep/brainsteminjected/"
            "sub-{subject}_label-brainstem_desc-proxy_from-template_to-subject_invwarp.nii.gz"
        ),
        warped=(
            "sub-{subject}/anat/mask_sweep/brainsteminjected/"
            "sub-{subject}_label-brainstem_desc-proxywarped_mask.nii.gz"
        ),
    log:
        "logs/sub-{subject}/mask_sweep/brainsteminjected/sub-{subject}_brainstem_proxywarp.log",
    benchmark:
        "benchmarks/sub-{subject}/mask_sweep/brainsteminjected/sub-{subject}_brainstem_proxywarp.tsv",
    threads: lambda wc: config.get("mrtrix", {}).get("proxy_threads", 8)
    resources:
        mem_mb=16000,
        time=60
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.affine}")"
        mkdir -p "$(dirname "{output.warped}")"
        mkdir -p "$(dirname "{log}")"

        greedy -d 3 -threads {threads} \
          -a \
          -moments 2 \
          -m SSD \
          -i "{input.fixed}" "{input.moving}" \
          -o "{output.affine}" \
          -n 100x50x10 \
          &> "{log}"

        greedy -d 3 -threads {threads} \
          -it "{output.affine}" \
          -m SSD \
          -i "{input.fixed}" "{input.moving}" \
          -o "{output.warp}" \
          -oinv "{output.invwarp}" \
          -n 100x50x10 \
          -s 1.732vox 0.707vox \
          -e 1.0 \
          &>> "{log}"

        greedy -d 3 -threads {threads} \
          -rf "{input.fixed}" \
          -rm "{input.moving}" "{output.warped}" \
          -r "{output.warp}" "{output.affine}" \
          &>> "{log}"
        """


rule sweep_warp_seed_with_brainstem_proxy:
    input:
        seed=template_seed_probseg_path,
        fixed=rules.sweep_extract_subject_brainstem.output.mask,
        affine=rules.sweep_brainstem_proxy_to_subject.output.affine,
        warp=rules.sweep_brainstem_proxy_to_subject.output.warp,
    output:
        probseg=(
            "sub-{subject}/anat/mask_sweep/{seed}/brainsteminjected/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-brainsteminjected_probseg.nii.gz"
        ),
    log:
        "logs/sub-{subject}/mask_sweep/{seed}/brainsteminjected/sub-{subject}_hemi-{hemi}_label-{seed}_brainsteminjected_seedwarp.log",
    benchmark:
        "benchmarks/sub-{subject}/mask_sweep/{seed}/brainsteminjected/sub-{subject}_hemi-{hemi}_label-{seed}_brainsteminjected_seedwarp.tsv",
    threads: lambda wc: config.get("mrtrix", {}).get("proxy_threads", 8)
    resources:
        mem_mb=16000,
        time=30
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.probseg}")"
        mkdir -p "$(dirname "{log}")"

        greedy -d 3 -threads {threads} \
          -rf "{input.fixed}" \
          -rm "{input.seed}" "{output.probseg}" \
          -r "{input.warp}" "{input.affine}" \
          &> "{log}"
        """


# -----------------------------
# Threshold seed source
# -----------------------------

rule sweep_binarize_seed:
    input:
        seed_prob=sweep_seed_probseg,
    output:
        seed_mask=(
            "sub-{subject}/anat/mask_sweep/{seed}/{mask_source}/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_mask.nii.gz"
        ),
    log:
        "logs/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_binarize.log",
    benchmark:
        "benchmarks/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_binarize.tsv",
    params:
        threshold=lambda wc: wc.thr_tag.replace("p", "."),
        resample_res=lambda wc: config["resample_seed_res"],
    threads: 1
    resources:
        mem_mb=4000,
        time=10
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.seed_mask}")"
        mkdir -p "$(dirname "{log}")"

        c3d -int 0 "{input.seed_prob}" \
          -threshold {params.threshold} inf 1 0 \
          -resample-mm {params.resample_res} \
          -trim 0vox \
          -type uchar \
          -o "{output.seed_mask}" \
          &> "{log}"
        """


# -----------------------------
# Filter existing tractogram by thresholded seed
# -----------------------------

rule sweep_filter_tractogram:
    input:
        tractogram=ancient(
            "sub-{subject}/tracts/sub-{subject}_desc-final_method-mrtrix_tractography.tck"
        ),
        sift2_weights=ancient(
            "sub-{subject}/tracts/sub-{subject}_desc-final_method-mrtrix_sift2_weights.txt"
        ),
        seed_mask=rules.sweep_binarize_seed.output.seed_mask,
    output:
        tractogram=temp(
            config["tmp_dir"]
            + "/sub-{subject}/tracts/mask_sweep/{seed}/{mask_source}/"
            + "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-seedfiltered_tractography.tck"
        ),
        sift2_weights=temp(
            config["tmp_dir"]
            + "/sub-{subject}/tracts/mask_sweep/{seed}/{mask_source}/"
            + "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-seedfiltered_sift2_weights.txt"
        ),
        tckinfo=(
            "sub-{subject}/qc/mask_sweep/{seed}/{mask_source}/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-seedfiltered_method-mrtrix_tractography_tckinfo.txt"
        ),
    log:
        "logs/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_filter_tractogram.log",
    benchmark:
        "benchmarks/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_filter_tractogram.tsv",
    threads: lambda wc: config.get("mrtrix", {}).get("filter_threads", 4)
    resources:
        mem_mb=16000,
        time=90
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tractogram}")"
        mkdir -p "$(dirname "{output.tckinfo}")"
        mkdir -p "$(dirname "{log}")"

        tckedit "{input.tractogram}" "{output.tractogram}" \
          -include "{input.seed_mask}" \
          -ends_only \
          -tck_weights_in "{input.sift2_weights}" \
          -tck_weights_out "{output.sift2_weights}" \
          &> "{log}"

        tckinfo "{output.tractogram}" > "{output.tckinfo}" 2>> "{log}"
        """

# -----------------------------
# Assign filtered streamlines to targets
# -----------------------------

rule sweep_targets_assignments:
    input:
        tractogram=rules.sweep_filter_tractogram.output.tractogram,
        dseg=ancient("sub-{subject}/anat/sub-{subject}_desc-{targets}_dseg.nii.gz"),
    output:
        assignments=(
            "sub-{subject}/qc/mask_sweep/{seed}/{mask_source}/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_method-mrtrix_assignments.txt"
        ),
    log:
        "logs/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_assignments.log",
    benchmark:
        "benchmarks/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_assignments.tsv",
    threads: lambda wc: config.get("mrtrix", {}).get("assign_threads", 4)
    resources:
        mem_mb=16000,
        time=30
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.assignments}")"
        mkdir -p "$(dirname "{log}")"

        tmpmat="$(mktemp)"
        trap 'rm -f "$tmpmat"' EXIT

        tck2connectome -nthreads {threads} -quiet -force \
          -assignment_radial_search 2 \
          -out_assignments "{output.assignments}" \
          "{input.tractogram}" "{input.dseg}" "$tmpmat" \
          &> "{log}"
        """

# -----------------------------
# Voxelwise aggregation
# -----------------------------

rule sweep_voxelwise_connectivity:
    input:
        tractogram=rules.sweep_filter_tractogram.output.tractogram,
        sift2_weights=rules.sweep_filter_tractogram.output.sift2_weights,
        seed_mask=rules.sweep_binarize_seed.output.seed_mask,
        assignments=rules.sweep_targets_assignments.output.assignments,
        targets_dseg=ancient("sub-{subject}/anat/sub-{subject}_desc-{targets}_dseg.nii.gz"),
    output:
        connectivity_matrix=(
            "sub-{subject}/tracts/mask_sweep/{seed}/{mask_source}/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_connectivity_matrix.csv"
        ),
        seed_voxel_index=(
            "sub-{subject}/qc/mask_sweep/{seed}/{mask_source}/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_seed_voxel_index.csv"
        ),
        qc_metrics=(
            "sub-{subject}/qc/mask_sweep/{seed}/{mask_source}/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_connectivity_qc.json"
        ),
        global_target_counts=temp(
            config["tmp_dir"]
            + "/sub-{subject}/qc/mask_sweep/{seed}/{mask_source}/"
            + "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_global_target_counts.csv"
        ),
    log:
        "logs/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_voxelwise_connectivity.log",
    benchmark:
        "benchmarks/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_voxelwise_connectivity.tsv",
    params:
        script=lambda wc: VOXEL_AGG_SCRIPT,
        target_labels=lambda wc: ",".join(config["targets"][wc.targets]["labels"]),
    threads: lambda wc: config.get("mrtrix", {}).get("voxelagg_threads", 4)
    resources:
        mem_mb=16000,
        time=90
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.connectivity_matrix}")"
        mkdir -p "$(dirname "{output.seed_voxel_index}")"
        mkdir -p "$(dirname "{output.qc_metrics}")"
        mkdir -p "$(dirname "{output.global_target_counts}")"
        mkdir -p "$(dirname "{log}")"

        python "{params.script}" \
          --seed-tck "{input.tractogram}" \
          --seed-weights "{input.sift2_weights}" \
          --seed-mask "{input.seed_mask}" \
          --assignments "{input.assignments}" \
          --targets-nii "{input.targets_dseg}" \
          --header "{params.target_labels}" \
          --out-weighted "{output.connectivity_matrix}" \
          --out-counts "{output.global_target_counts}" \
          --out-voxel-index "{output.seed_voxel_index}" \
          --out-qc "{output.qc_metrics}" \
          &> "{log}"
        """
