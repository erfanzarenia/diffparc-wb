# microstructure_connectomes.smk -- optional, additive module (MAIN pipeline).
#
# For the whole-brain (final) tractogram + its SIFT2 weights, sample diffusion
# microstructure (FA, MD, RD) along each streamline, then build SIFT2-weighted,
# mean-edge atlas connectomes:
#   1. tcksample -stat_tck mean            -> one mean metric per streamline
#   2. tck2connectome -scale_file <metric> -stat_edge mean -tck_weights_in <sift2>
#
# Reuses existing products: final_tractogram, run_sift2, and the tensor-derived
# FA/MD (dwi_to_tensor / tensor_to_metrics). Only RD is generated here, from the
# SAME tensor. Off by default (config microstructure_connectomes); outputs are
# preserved (not temp) for QC/reuse.

import os


def res(rule, key, default):
    return config.get("resources", {}).get(rule, {}).get(key, default)


LABEL_CONNECTOME_SCRIPT = os.path.join(
    workflow.basedir, "scripts", "label_connectome_csv.py"
)

MICRO_ENABLED = bool(config.get("microstructure_connectomes", False))
MICRO_TARGETS = config.get("microstructure_connectomes_targets", ["Yeo7TianS3"])
MICRO_METRICS = config.get("microstructure_connectomes_metrics", ["FA", "MD", "RD"])
MICRO_SUBJECTS = sorted(set(inputs.input_zip_lists["T1w"]["subject"]))


# `metric` (FA|MD|RD) is constrained per-rule below rather than globally, to
# keep the wildcard scoped to this module.


# -----------------------------
# (1) RD map (FA + MD already come from rules.tensor_to_metrics; reuse them).
#     Same tensor + mask as tensor_to_metrics -> MSMT-CSD-compatible. Preserved.
# -----------------------------
rule tensor_to_rd:
    input:
        tensor=rules.dwi_to_tensor.output.tensor,
        mask=rules.nii2mif.output.mask,
    output:
        rd=bids(root=root, datatype="dwi", suffix="RD.nii.gz", **subj_wildcards),
    log:
        "logs/sub-{subject}/dwi/sub-{subject}_tensor_to_rd.log",
    benchmark:
        "benchmarks/sub-{subject}/dwi/sub-{subject}_tensor_to_rd.tsv"
    threads: lambda wc: res("tensor_to_rd", "threads", 1)
    resources:
        mem_mb=lambda wc: res("tensor_to_rd", "mem_mb", 2000),
        time=lambda wc: res("tensor_to_rd", "time_min", 10),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.rd}")"
        mkdir -p "$(dirname "{log}")"

        tensor2metric -nthreads {threads} -mask "{input.mask}" \
          -rd "{output.rd}" "{input.tensor}" &> "{log}"
        """


# -----------------------------
# (2) Per-streamline mean metric along the whole-brain tractogram (preserved).
#     The {metric} map (FA/MD/RD) is selected by the wildcard; FA.nii.gz / MD.nii.gz
#     come from tensor_to_metrics, RD.nii.gz from tensor_to_rd above.
# -----------------------------
rule microstructure_sample:
    wildcard_constraints:
        metric="FA|MD|RD",
    input:
        tck=rules.final_tractogram.output.tck,
        metric_map=bids(
            root=root, datatype="dwi", suffix="{metric}.nii.gz", **subj_wildcards
        ),
    output:
        sample=(
            "sub-{subject}/tracts/"
            "sub-{subject}_desc-final_metric-{metric}_streamlinesample.csv"
        ),
    log:
        "logs/sub-{subject}/connectivity/sub-{subject}_desc-final_metric-{metric}_tcksample.log",
    benchmark:
        "benchmarks/sub-{subject}/connectivity/sub-{subject}_desc-final_metric-{metric}_tcksample.tsv"
    threads: lambda wc: res("microstructure_sample", "threads", 2)
    resources:
        mem_mb=lambda wc: res("microstructure_sample", "mem_mb", 4000),
        time=lambda wc: res("microstructure_sample", "time_min", 75),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.sample}")"
        mkdir -p "$(dirname "{log}")"

        tcksample -nthreads {threads} -force -stat_tck mean \
          "{input.tck}" "{input.metric_map}" "{output.sample}" &> "{log}"
        """


# -----------------------------
# (3) Mean-metric, SIFT2-weighted atlas connectome (clearly named .csv).
#     Edge value = SIFT2-weighted mean of the per-streamline metric over the
#     streamlines assigned to that parcel pair (-scale_file + -stat_edge mean +
#     -tck_weights_in). Same atlas dseg + radial assignment as the main pipeline.
# -----------------------------
rule microstructure_connectome:
    wildcard_constraints:
        metric="FA|MD|RD",
    input:
        tck=rules.final_tractogram.output.tck,
        sift2_weights=rules.run_sift2.output.weights,
        scale=rules.microstructure_sample.output.sample,
        nodes=rules.transform_targets_to_subject.output.targets,
    output:
        connectivity_matrix=(
            "sub-{subject}/tracts/"
            "sub-{subject}_desc-{targets}_meas-{metric}_connectivity_matrix.csv"
        ),
        connectome_raw=temp(
            config["tmp_dir"]
            + "/sub-{subject}/tracts/"
            + "sub-{subject}_desc-{targets}_meas-{metric}_connectome_raw.txt"
        ),
    params:
        label_script=lambda wc: LABEL_CONNECTOME_SCRIPT,
        target_labels=lambda wc: ",".join(config["targets"][wc.targets]["labels"]),
        radius=lambda wc: config.get("mrtrix", {}).get("target_search_radius", 4),
    log:
        "logs/sub-{subject}/connectivity/sub-{subject}_desc-{targets}_meas-{metric}_connectome.log",
    benchmark:
        "benchmarks/sub-{subject}/connectivity/sub-{subject}_desc-{targets}_meas-{metric}_connectome.tsv"
    threads: lambda wc: res("microstructure_connectome", "threads", 2)
    resources:
        mem_mb=lambda wc: res("microstructure_connectome", "mem_mb", 8000),
        time=lambda wc: res("microstructure_connectome", "time_min", 30),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.connectivity_matrix}")"
        mkdir -p "$(dirname "{output.connectome_raw}")"
        mkdir -p "$(dirname "{log}")"

        tck2connectome -nthreads {threads} -force -symmetric \
          -assignment_radial_search {params.radius} \
          -tck_weights_in "{input.sift2_weights}" \
          -scale_file "{input.scale}" \
          -stat_edge mean \
          "{input.tck}" "{input.nodes}" "{output.connectome_raw}" \
          &> "{log}"

        python "{params.label_script}" \
          --matrix "{output.connectome_raw}" \
          --labels "{params.target_labels}" \
          --out-csv "{output.connectivity_matrix}" \
          &>> "{log}"
        """


# -----------------------------
# Convenience aggregate target (also wired into `rule all` via the Snakefile helper)
# -----------------------------
rule all_microstructure_connectomes:
    input:
        expand(
            "sub-{subject}/tracts/"
            "sub-{subject}_desc-{targets}_meas-{metric}_connectivity_matrix.csv",
            subject=MICRO_SUBJECTS, targets=MICRO_TARGETS, metric=MICRO_METRICS,
        ) if MICRO_ENABLED else [],
