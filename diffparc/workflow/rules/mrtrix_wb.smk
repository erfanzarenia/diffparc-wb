# mrtrix_wb.smk
# Whole-brain tractography + optional ROI enrichment + SIFT2 (WITH ACT)
# Then: endpoint-in-seed filtering + per-streamline target assignments + voxelwise seed->target aggregation (CSV headers from config)

import os

# Script path
VOXEL_AGG_SCRIPT = os.path.join(workflow.basedir, "scripts", "voxel_target_aggregate.py")

def res(rule, key, default):
    return config.get("resources", {}).get(rule, {}).get(key, default)

# ROI toggle + seeds list
def roi_enabled():
    return config.get("roi_enrichment", {}).get("enabled", False)

SEEDS = ["vtasnc"]

WB_CHUNKS = list(range(1, int(config.get("wb_chunk_count", 1)) + 1))

def wb_chunk_seed(wc):
    return int(config.get("mrtrix_rng_seed", 42)) + int(wc.chunk)

def wb_chunk_streamlines():
    return int(config.get("wb_chunk_streamlines", config.get("wb_sl_count", 20000000)))
 
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
        "-resample-mm {params.resample_res} -trim 0vox "
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
    threads: lambda wc: res("act_5ttgen", "threads", 8)
    resources:
        mem_mb=lambda wc: res("act_5ttgen", "mem_mb", 32000),
        time=lambda wc: res("act_5ttgen", "time_min", 180)
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
    threads: lambda wc: res("act_gmwmi", "threads", 2)
    resources:
        mem_mb=lambda wc: res("act_gmwmi", "mem_mb", 8000),
        time=lambda wc: res("act_gmwmi", "time_min", 30)
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


# -----------------------------
# Whole-brain tractography (ACT)
# -----------------------------

rule wb_tckgen_chunk:
    input:
        wm_fod=get_fod_for_tracking,
        mask=bids(root=root, datatype="dwi", suffix="mask.mif", **subj_wildcards),
        five_tt=rules.act_5ttgen.output.five_tt,
        gmwmi=rules.act_gmwmi.output.gmwmi,
    params:
        step=lambda wc: config.get("mrtrix", {}).get("step", 0.35),
        angle=lambda wc: config.get("mrtrix", {}).get("angle", 45),
        minlength=lambda wc: config.get("mrtrix", {}).get("minlength", 3),
        maxlength=lambda wc: config.get("mrtrix", {}).get("maxlength", 250),
        cutoff=lambda wc: config.get("mrtrix", {}).get("cutoff", 0.06),
        n_streamlines=lambda wc: wb_chunk_streamlines(),
        mrtrix_rng_seed=lambda wc: wb_chunk_seed(wc),
    output:
        tck=temp(
            "sub-{subject}/tracts/sub-{subject}_desc-wb_chunk-{chunk}_method-mrtrix_tractography.tck"
        ),
        tckinfo="sub-{subject}/tracts/sub-{subject}_desc-wb_chunk-{chunk}_method-mrtrix_tractography_tckinfo.txt",
    log:
        "logs/sub-{subject}/tracts/sub-{subject}_desc-wb_chunk-{chunk}_tckgen.log",
    benchmark:
        "benchmarks/sub-{subject}/tracts/sub-{subject}_desc-wb_chunk-{chunk}_tckgen.tsv"
    threads: lambda wc: res("wb_tckgen_chunk", "threads", 10)
    resources:
        mem_mb=lambda wc: res("wb_tckgen_chunk", "mem_mb", 40000),
        time=lambda wc: res("wb_tckgen_chunk", "time_min", 1200)
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tck}")"
        mkdir -p "$(dirname "{log}")"
        export MRTRIX_RNG_SEED={params.mrtrix_rng_seed}

        tckgen -nthreads {threads} -algorithm iFOD2 \
          -act "{input.five_tt}" \
          -backtrack \
          -crop_at_gmwmi \
          -seed_gmwmi "{input.gmwmi}" \
          -mask "{input.mask}" \
          -step {params.step} -angle {params.angle} \
          -minlength {params.minlength} -maxlength {params.maxlength} \
          -cutoff {params.cutoff} \
          -select {params.n_streamlines} \
          "{input.wm_fod}" "{output.tck}" \
          &> "{log}"

        tckinfo "{output.tck}" > "{output.tckinfo}" 2>> "{log}"
        """

rule wb_tckgen:
    input:
        tcks=lambda wc: expand(
            "sub-{subject}/tracts/sub-{subject}_desc-wb_chunk-{chunk}_method-mrtrix_tractography.tck",
            subject=wc.subject,
            chunk=WB_CHUNKS
        ),
        tckinfos=lambda wc: expand(
            "sub-{subject}/tracts/sub-{subject}_desc-wb_chunk-{chunk}_method-mrtrix_tractography_tckinfo.txt",
            subject=wc.subject,
            chunk=WB_CHUNKS
        ),
    output:
        tck=bids(
            root=root, datatype="tracts", method="mrtrix",
            desc="wb", suffix="tractography.tck", **subj_wildcards
        ),
        tckinfo=bids(
            root=root, datatype="tracts", method="mrtrix",
            desc="wb", suffix="tractography_tckinfo.txt", **subj_wildcards
        ),
    log:
        "logs/sub-{subject}/tracts/sub-{subject}_desc-wb_merge.log",
    threads: lambda wc: res("wb_tckgen_merge", "threads", 2)
    resources:
        mem_mb=lambda wc: res("wb_tckgen_merge", "mem_mb", 8000),
        time=lambda wc: res("wb_tckgen_merge", "time_min", 60)
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tck}")"
        mkdir -p "$(dirname "{log}")"

        tckedit -force {input.tcks} "{output.tck}" &> "{log}"
        tckinfo "{output.tck}" > "{output.tckinfo}" 2>> "{log}"
        """
        
# -----------------------------
# Optional ROI enrichment
# -----------------------------
rule roi_tckgen:
    """
    ROI-targeted enrichment tractography (UNIDIRECTIONAL) with ACT.
    Uses -seed_random_per_voxel as seeding attempts per voxel.
    """
    input:
        wm_fod=get_fod_for_tracking,
        mask=bids(root=root, datatype="dwi", suffix="mask.mif", **subj_wildcards),
        seed_mask=bids(
            root=root, **subj_wildcards, hemi="{hemi}", label="{seed}",
            datatype="anat", suffix="mask.nii.gz"
        ),
        five_tt=rules.act_5ttgen.output.five_tt,
    params:
        spv=lambda wc: config["roi_enrichment"]["seeds_per_voxel"],
        step=lambda wc: config.get("mrtrix", {}).get("step", 0.35),
        angle=lambda wc: config.get("mrtrix", {}).get("angle", 45),
        minlength=lambda wc: config.get("mrtrix", {}).get("minlength", 3),
        maxlength=lambda wc: config.get("mrtrix", {}).get("maxlength", 250),
        cutoff=lambda wc: config.get("mrtrix", {}).get("cutoff", 0.06),
        rng=lambda wc: {
            "L": int(config.get("mrtrix_rng_seed", 42)),
            "R": int(config.get("mrtrix_rng_seed", 42)) + 1000,
        }[wc.hemi],
    output:
        tck=temp(bids(
            root=config["tmp_dir"], datatype="tracts", method="mrtrix",
            desc="roi", hemi="{hemi}", label="{seed}",
            suffix="tractography.tck", **subj_wildcards
        )),
        tckinfo=bids(
            root=root, datatype="tracts", method="mrtrix",
            desc="roi", hemi="{hemi}", label="{seed}",
            suffix="tractography_tckinfo.txt", **subj_wildcards
        ),
    log:
        "logs/sub-{subject}/tracts/sub-{subject}_hemi-{hemi}_label-{seed}_roi_tckgen.log",
    benchmark:
        "benchmarks/sub-{subject}/tracts/sub-{subject}_hemi-{hemi}_label-{seed}_roi_tckgen.tsv"
    threads: lambda wc: res("roi_tckgen", "threads", 6)
    resources:
        mem_mb=lambda wc: res("roi_tckgen", "mem_mb", 24000),
        time=lambda wc: res("roi_tckgen", "time_min", 120)
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        export MRTRIX_RNG_SEED={params.rng}
        mkdir -p "$(dirname "{output.tck}")"
        mkdir -p "$(dirname "{output.tckinfo}")"
        mkdir -p "$(dirname "{log}")"
    
        tckgen -nthreads {threads} -algorithm iFOD2 \
          -act "{input.five_tt}" \
          -backtrack \
          -crop_at_gmwmi \
          -mask "{input.mask}" \
          -seed_unidirectional \
          -step {params.step} -angle {params.angle} \
          -minlength {params.minlength} -maxlength {params.maxlength} \
          -cutoff {params.cutoff} \
          -seed_random_per_voxel "{input.seed_mask}" {params.spv} \
          "{input.wm_fod}" "{output.tck}" \
          &> "{log}"

        tckinfo "{output.tck}" > "{output.tckinfo}" 2>> "{log}"
        """

rule roi_merge_all:
    input:
        lambda wc: expand(
            bids(
                root=config["tmp_dir"], datatype="tracts", method="mrtrix",
                desc="roi", hemi="{hemi}", label="{seed}",
                suffix="tractography.tck", **subj_wildcards
            ),
            subject=wc.subject,
            hemi=["L", "R"],
            seed=SEEDS
        ) if roi_enabled() else []
    output:
        tck=temp(bids(
            root=config["tmp_dir"], datatype="tracts", method="mrtrix",
            desc="roi", suffix="tractography_merged.tck", **subj_wildcards
        )),
        tckinfo=bids(
            root=root, datatype="tracts", method="mrtrix",
            desc="roi", suffix="tractography_merged_tckinfo.txt", **subj_wildcards
        ),
    log:
        "logs/sub-{subject}/tracts/sub-{subject}_desc-roi_merge_all.log",
    params:
        roi_on=lambda wc: 1 if roi_enabled() else 0
    threads: lambda wc: res("roi_merge_all", "threads", 2)
    resources:
        mem_mb=lambda wc: res("roi_merge_all", "mem_mb", 8000),
        time=lambda wc: res("roi_merge_all", "time_min", 30)
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tck}")"
        mkdir -p "$(dirname "{output.tckinfo}")"
        mkdir -p "$(dirname "{log}")"

        if [ "{params.roi_on}" -eq 1 ]; then
          inputs=""
          for f in {input}; do
            if [ -s "$f" ]; then
              inputs="$inputs $f"
            fi
          done

          if [ -z "$inputs" ]; then
            echo "roi_merge_all: ROI enabled but no non-empty input tcks found" >&2
            exit 1
          fi

          tckedit -force $inputs "{output.tck}" &> "{log}"
          tckinfo "{output.tck}" > "{output.tckinfo}" 2>> "{log}"
        else
          : > "{output.tck}"
          printf "ROI disabled\n" > "{output.tckinfo}"
          printf "ROI disabled\n" > "{log}"
        fi
        """


# -----------------------------
# Merge WB + ROI tractograms (or WB only)
# -----------------------------
rule tractogram_for_sift2:
    input:
        wb=rules.wb_tckgen.output.tck,
        roi=rules.roi_merge_all.output.tck
    output:
        tck=bids(
            root=root, datatype="tracts", method="mrtrix",
            desc="wbplusroi", suffix="tractography.tck", **subj_wildcards
        ),
        tckinfo=bids(
            root=root, datatype="tracts", method="mrtrix",
            desc="wbplusroi", suffix="tractography_tckinfo.txt", **subj_wildcards
        ),
    log:
        "logs/sub-{subject}/tracts/sub-{subject}_desc-wbplusroi_merge.log",
    threads: lambda wc: res("tractogram_for_sift2", "threads", 2)
    resources:
        mem_mb=lambda wc: res("tractogram_for_sift2", "mem_mb", 8000),
        time=lambda wc: res("tractogram_for_sift2", "time_min", 90)
    container: config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tck}")"
        mkdir -p "$(dirname "{log}")"

        if [ -s "{input.roi}" ]; then
          tckedit -force "{input.wb}" "{input.roi}" "{output.tck}" &> "{log}"
        else
          ln -sf "$(readlink -f "{input.wb}")" "{output.tck}" &> "{log}"
        fi

        tckinfo "{output.tck}" > "{output.tckinfo}" 2>> "{log}"
        """


# -----------------------------
# SIFT2 weights (ACT-aware)
# -----------------------------
rule wb_tcksift2:
    input:
        tck=rules.tractogram_for_sift2.output.tck,
        wm_fod=get_fod_for_tracking,
        five_tt=rules.act_5ttgen.output.five_tt,
    output:
        weights=bids(
            root=root, datatype="tracts", method="mrtrix",
            desc="wbplusroi", suffix="sift2_weights.txt", **subj_wildcards
        ),
        mu=bids(
            root=root, datatype="tracts", method="mrtrix",
            desc="wbplusroi", suffix="sift2_mu.txt", **subj_wildcards
        ),
    log:
        "logs/sub-{subject}/tracts/sub-{subject}_desc-wbplusroi_method-mrtrix_tcksift2.log",
    benchmark:
        "benchmarks/sub-{subject}/tracts/sub-{subject}_desc-wbplusroi_method-mrtrix_tcksift2.tsv",
    threads: lambda wc: res("wb_tcksift2", "threads", 10)
    resources:
        mem_mb=lambda wc: res("wb_tcksift2", "mem_mb", 40000),
        time=lambda wc: res("wb_tcksift2", "time_min", 300)
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:    
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"

        tcksift2 -nthreads {threads} -act "{input.five_tt}" -out_mu "{output.mu}" \
          "{input.tck}" "{input.wm_fod}" "{output.weights}" \
          &> "{log}"
        """


# -----------------------------
# Filter seed endpoints (and subset weights)
# -----------------------------
rule wb_filter_seed_tck:
    """
    Keep streamlines where AT LEAST ONE ENDPOINT is inside the seed mask.
    """
    input:
        wb_tck=rules.tractogram_for_sift2.output.tck,
        wb_weights=rules.wb_tcksift2.output.weights,
        seed_mask=bids(
            root=root, **subj_wildcards, hemi="{hemi}", label="{seed}",
            datatype="anat", suffix="mask.nii.gz"
        ),
    output:
        seed_tck=temp(bids(
            root=config["tmp_dir"], datatype="tracts", method="mrtrix",
            desc="wbplusroi", hemi="{hemi}", label="{seed}",
            suffix="seed_endpoints.tck", **subj_wildcards
        )),
        seed_weights=temp(bids(
            root=config["tmp_dir"], datatype="tracts", method="mrtrix",
            desc="wbplusroi", hemi="{hemi}", label="{seed}",
            suffix="seed_endpoints_weights.txt", **subj_wildcards
        )),
        seed_tckinfo=bids(
            root=root, datatype="tracts", hemi="{hemi}", label="{seed}",
            method="mrtrix", desc="wbplusroi",
            suffix="seed_endpoints_tckinfo.txt", **subj_wildcards
        ),
    log:
        "logs/sub-{subject}/tracts/sub-{subject}_hemi-{hemi}_label-{seed}_wb_filter_seed_tck.log",
    benchmark:
        "benchmarks/sub-{subject}/tracts/sub-{subject}_hemi-{hemi}_label-{seed}_wb_filter_seed_tck.tsv",
    threads: lambda wc: res("wb_filter_seed_tck", "threads", 4)
    resources:
        mem_mb=lambda wc: res("wb_filter_seed_tck", "mem_mb", 16000),
        time=lambda wc: res("wb_filter_seed_tck", "time_min", 90)
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.seed_tck}")"
        mkdir -p "$(dirname "{output.seed_tckinfo}")"
        mkdir -p "$(dirname "{log}")"

        tckedit "{input.wb_tck}" "{output.seed_tck}" \
          -include "{input.seed_mask}" \
          -ends_only \
          -tck_weights_in "{input.wb_weights}" \
          -tck_weights_out "{output.seed_weights}" \
          &> "{log}"

        tckinfo "{output.seed_tck}" > "{output.seed_tckinfo}" 2>> "{log}"
        """


# -----------------------------
# Per-streamline target assignments (MRtrix)
# -----------------------------
rule wb_targets_assignments:
    input:
        tck=rules.wb_filter_seed_tck.output.seed_tck,
        dseg_mif="sub-{subject}/anat/sub-{subject}_desc-{targets}_dseg.mif",
    output:
        assignments=bids(
            root=root,
            datatype="tracts",
            method="mrtrix",
            desc="{targets}",
            hemi="{hemi}",
            label="{seed}",
            suffix="assignments.txt",
            **subj_wildcards
        ),
    log:
        "logs/sub-{subject}/tracts/sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}_assignments.log",
    benchmark:
        "benchmarks/sub-{subject}/tracts/sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}_assignments.tsv",
    threads: lambda wc: res("wb_targets_assignments", "threads", 4)
    resources:
        mem_mb=lambda wc: res("wb_targets_assignments", "mem_mb", 16000),
        time=lambda wc: res("wb_targets_assignments", "time_min", 30)
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
          "{input.tck}" "{input.dseg_mif}" "$tmpmat" \
          &> "{log}"
        """


# -----------------------------
# Voxelwise seed->targets aggregation (Python)
# -----------------------------
rule wb_voxelwise_seed_to_targets_matrix:
    input:
        seed_tck=rules.wb_filter_seed_tck.output.seed_tck,
        seed_weights=rules.wb_filter_seed_tck.output.seed_weights,
        seed_mask="sub-{subject}/anat/sub-{subject}_hemi-{hemi}_label-{seed}_mask.nii.gz",
        assignments=rules.wb_targets_assignments.output.assignments,
        targets_nii="sub-{subject}/anat/sub-{subject}_desc-{targets}_dseg.nii.gz",
    output:
        weighted="sub-{subject}/tracts/sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}_sift2_conn.csv",
        voxel_index="sub-{subject}/tracts/sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}_seed_voxel_index.csv",
        qc="sub-{subject}/tracts/sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}_voxelagg_qc.json",
        counts=temp("{}".format(
            config["tmp_dir"] + "/sub-{subject}/tracts/sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}_conn_voxelwise_global_counts.csv"
        )),
    log:
        "logs/sub-{subject}/tracts/sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}_voxelagg.log",
    benchmark:
        "benchmarks/sub-{subject}/tracts/sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}_voxelagg.tsv",
    params:
        script=lambda wc: VOXEL_AGG_SCRIPT,
        header=lambda wc: ",".join(config["targets"][wc.targets]["labels"]),
    threads: lambda wc: res("wb_voxelwise_seed_to_targets_matrix", "threads", 4)
    resources:
        mem_mb=lambda wc: res("wb_voxelwise_seed_to_targets_matrix", "mem_mb", 16000),
        time=lambda wc: res("wb_voxelwise_seed_to_targets_matrix", "time_min", 90)
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.weighted}")"
        mkdir -p "$(dirname "{output.counts}")"
        mkdir -p "$(dirname "{log}")"
    
        python "{params.script}" \
          --seed-tck "{input.seed_tck}" \
          --seed-weights "{input.seed_weights}" \
          --seed-mask "{input.seed_mask}" \
          --assignments "{input.assignments}" \
          --targets-nii "{input.targets_nii}" \
          --header "{params.header}" \
          --out-weighted "{output.weighted}" \
          --out-counts "{output.counts}" \
          --out-voxel-index "{output.voxel_index}" \
          --out-qc "{output.qc}" \
          &> "{log}"
        """
