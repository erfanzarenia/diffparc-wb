# wb_prep.smk -- prerequisites for tracking: binarize/trim the subject seed
# masks, and build the 5TT image (+ GMWMI) that makes tracking and SIFT2
# anatomically constrained (ACT).


def res(rule, key, default):
    return config.get("resources", {}).get(rule, {}).get(key, default)


# -----------------------------
# Seed binarization + cleanup
# -----------------------------
rule binarize_trim_subject_seed:
    input:
        seed=get_subject_seed_probseg,
    params:
        threshold=lambda wc: config["seeds"][wc.seed]["probseg_threshold"],
        resample_res=lambda wc: config["resample_seed_res"],
    output:
        seed_thr=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            datatype="anat",
            suffix="mask.nii.gz"
        ),
    log:
        bids(
            root="logs",
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            suffix="binarizetrimsubjectseed.log"
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        "c3d -int 0 {input.seed} -threshold {params.threshold} inf 1 0 "
        "-dilate 1 3x3x3vox -resample-mm {params.resample_res} -trim 0vox "
        "-type uchar -o {output.seed_thr} &> {log}"


# -----------------------------
# ACT prerequisites: 5TT + GMWMI
# -----------------------------
rule act_5ttgen:
    """
    Build 5TT image for ACT from T1.
    NOTE: Assumes your preproc T1 path is:
          sub-{subject}/anat/sub-{subject}_desc-preproc_T1w.nii.gz
    """
    input:
        t1="sub-{subject}/anat/sub-{subject}_desc-preproc_T1w.nii.gz",
    output:
        five_tt="sub-{subject}/anat/sub-{subject}_desc-5tt_method-mrtrix.mif",
    log:
        "logs/sub-{subject}/act_5ttgen.log",
    threads: lambda wc: res("act_5ttgen", "threads", 4)
    resources:
        mem_mb=lambda wc: res("act_5ttgen", "mem_mb", 12000),
        time=lambda wc: res("act_5ttgen", "time_min", 60)
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.five_tt}")"
        mkdir -p "$(dirname "{log}")"

        5ttgen fsl -force "{input.t1}" "{output.five_tt}" &> "{log}"
        5ttcheck "{output.five_tt}" >> "{log}" 2>&1
        """


rule act_gmwmi:
    """
    Build GMWMI seed image for ACT seeding.
    """
    input:
        five_tt=rules.act_5ttgen.output.five_tt,
    output:
        gmwmi="sub-{subject}/anat/sub-{subject}_desc-gmwmi_method-mrtrix.mif",
    log:
        "logs/sub-{subject}/act_gmwmi.log",
    threads: lambda wc: res("act_gmwmi", "threads", 1)
    resources:
        mem_mb=lambda wc: res("act_gmwmi", "mem_mb", 4000),
        time=lambda wc: res("act_gmwmi", "time_min", 10)
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.gmwmi}")"
        mkdir -p "$(dirname "{log}")"

        5tt2gmwmi -force "{input.five_tt}" "{output.gmwmi}" &> "{log}"
        """


