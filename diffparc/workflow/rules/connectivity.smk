import os

def res(rule, key, default):
    return config.get("resources", {}).get(rule, {}).get(key, default)


VOXEL_AGG_SCRIPT = os.path.join(
    workflow.basedir,
    "scripts",
    "voxel_target_aggregate.py",
)


# -----------------------------
# Filter final tractogram by seed endpoints
# -----------------------------
rule filter_tractogram:
    """
    Extract seed-connected streamlines from the final tractogram.

    Keeps streamlines where at least one endpoint falls inside the seed mask,
    and writes the matching subset of SIFT2 weights.
    """
    input:
        tractogram=rules.final_tractogram.output.tck,
        sift2_weights=rules.run_sift2.output.weights,
        seed_mask=bids(
            root=root,
            datatype="anat",
            hemi="{hemi}",
            label="{seed}",
            suffix="mask.nii.gz",
            **subj_wildcards,
        ),
    output:
        tractogram=bids(
            root=root,
            datatype="tracts",
            method="mrtrix",
            desc="seedfiltered",
            hemi="{hemi}",
            label="{seed}",
            suffix="tractography.tck",
            **subj_wildcards,
        ),
        sift2_weights=bids(
            root=root,
            datatype="tracts",
            method="mrtrix",
            desc="seedfiltered",
            hemi="{hemi}",
            label="{seed}",
            suffix="sift2_weights.txt",
            **subj_wildcards,
        ),
        tractogram_info=(
            "sub-{subject}/qc/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}"
            "_desc-seedfiltered_tractography_tckinfo.txt"
        ),
    log:
        (
            "logs/sub-{subject}/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}"
            "_desc-seedfiltered_filter_tractogram.log"
        ),
    benchmark:
        (
            "benchmarks/sub-{subject}/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}"
            "_desc-seedfiltered_filter_tractogram.tsv"
        ),
    threads: lambda wc: res("filter_tractogram", "threads", 2)
    resources:
        mem_mb=lambda wc: res("filter_tractogram", "mem_mb", 16000),
        time=lambda wc: res("filter_tractogram", "time_min", 30)
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.tractogram}")"
        mkdir -p "$(dirname "{output.tractogram_info}")"
        mkdir -p "$(dirname "{log}")"

        tckedit "{input.tractogram}" "{output.tractogram}" \
          -include "{input.seed_mask}" \
          -ends_only \
          -tck_weights_in "{input.sift2_weights}" \
          -tck_weights_out "{output.sift2_weights}" \
          &> "{log}"

        tckinfo "{output.tractogram}" > "{output.tractogram_info}" 2>> "{log}"
        """


# -----------------------------
# Assign seed-filtered streamlines to targets
# -----------------------------
rule targets_assignments:
    input:
        tractogram=rules.filter_tractogram.output.tractogram,
        targets_dseg_mif="sub-{subject}/anat/sub-{subject}_desc-{targets}_dseg.mif",
    output:
        assignments=(
            "sub-{subject}/qc/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_method-mrtrix_assignments.txt"
        ),
    log:
        (
            "logs/sub-{subject}/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_targets_assignments.log"
        ),
    benchmark:
        (
            "benchmarks/sub-{subject}/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_targets_assignments.tsv"
        ),
    threads: lambda wc: res("targets_assignments", "threads", 2)
    resources:
        mem_mb=lambda wc: res("targets_assignments", "mem_mb", 4000),
        time=lambda wc: res("targets_assignments", "time_min", 5)
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
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
          "{input.tractogram}" "{input.targets_dseg_mif}" "$tmpmat" \
          &> "{log}"
        """


# -----------------------------
# Voxelwise seed-to-target connectivity
# -----------------------------
rule voxelwise_connectivity:
    """
    Generate voxelwise seed-to-target connectivity matrices from the
    seed-filtered tractogram and corresponding SIFT2 weights.
    """
    input:
        tractogram=rules.filter_tractogram.output.tractogram,
        sift2_weights=rules.filter_tractogram.output.sift2_weights,
        seed_mask=bids(
            root=root,
            datatype="anat",
            hemi="{hemi}",
            label="{seed}",
            suffix="mask.nii.gz",
            **subj_wildcards,
        ),
        assignments=rules.targets_assignments.output.assignments,
        targets_dseg="sub-{subject}/anat/sub-{subject}_desc-{targets}_dseg.nii.gz",
    output:
        connectivity_matrix=(
            "sub-{subject}/tracts/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_connectivity_matrix.csv"
        ),
        seed_voxel_index=(
            "sub-{subject}/qc/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_seed_voxel_index.csv"
        ),
        qc_metrics=(
            "sub-{subject}/qc/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_connectivity_qc.json"
        ),
        global_target_counts=temp(
            config["tmp_dir"]
            + "/sub-{subject}/qc/connectivity/"
            + "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            + "_global_target_counts.csv"
        ),
    log:
        (
            "logs/sub-{subject}/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_voxelwise_connectivity.log"
        ),
    benchmark:
        (
            "benchmarks/sub-{subject}/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_voxelwise_connectivity.tsv"
        ),
    params:
        script=lambda wc: VOXEL_AGG_SCRIPT,
        target_labels=lambda wc: ",".join(config["targets"][wc.targets]["labels"]),
    threads: lambda wc: res("voxelwise_connectivity", "threads", 2)
    resources:
        mem_mb=lambda wc: res("voxelwise_connectivity", "mem_mb", 4000),
        time=lambda wc: res("voxelwise_connectivity", "time_min", 15)
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
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
