# enrichment_sweep.smk
# Experimental, fully ISOLATED ROI-enrichment sweep (additive; mirrors mask_sweep.smk).
#
# Sweeps ROI seeding strength ("enrichment level") WITHOUT changing any core rule:
#   1. Generate ONE maximal-enrichment ROI tractogram (voxelwise ROI seeding at
#      spv = max(levels)). This is the ONLY ROI tractogram generation.
#   2. Lower levels = reproducible RANDOM SUBSAMPLE of that maximal ROI tractogram
#      (no re-seeding, no re-running tckgen).
#   3. WB baseline = ONE whole-brain tractogram (no ROI), reused read-only from
#      the core pipeline (rules.wb_tckgen_merge). This is the reference condition.
#   4. Per enrichment level: combine fixed WB + subsampled ROI (no recompute).
#   5. SIFT2 run independently per condition (WB baseline + each level), identical
#      params; logs mu, cost, weights.
#   6. Safety check (tckinfo) applied unchanged, once on WB baseline + once/condition.
#
# Only entry points into the rest of the workflow: this include line and the
# get_enrichment_sweep_outputs() helper added to `rule all`. Off by default.

import os


def res(rule, key, default):
    return config.get("resources", {}).get(rule, {}).get(key, default)


SUBSAMPLE_SCRIPT = os.path.join(workflow.basedir, "scripts", "subsample_tractogram.py")
METRICS_SCRIPT = os.path.join(workflow.basedir, "scripts", "enrichment_condition_metrics.py")
# Reused from the main connectivity pipeline (connectivity.smk) -- no duplicate logic.
SLICE_BLOCK_SCRIPT = os.path.join(workflow.basedir, "scripts", "slice_connectome_block.py")

ENRICH_CFG = config.get("enrichment_sweep", {})
ENRICH_ENABLED = ENRICH_CFG.get("enabled", False)
ENRICH_SEED = ENRICH_CFG.get("seed", "vtasncpbp")
ENRICH_HEMIS = ENRICH_CFG.get("hemis", ["L", "R"])
ENRICH_LEVELS = [int(x) for x in ENRICH_CFG.get("levels", [])]
ENRICH_MAX_LEVEL = max(ENRICH_LEVELS) if ENRICH_LEVELS else 0
ENRICH_SUBSAMPLE_SEED = int(ENRICH_CFG.get("subsample_seed", 42))
ENRICH_TARGETS = ENRICH_CFG.get("targets", ["Yeo7TianS3"])
# Optional: also emit a volume-bias-corrected fingerprint (tck2connectome
# -scale_invnodevol) per condition, in addition to the uncorrected one.
ENRICH_VOLNORM = bool(ENRICH_CFG.get("volume_bias_correction", False))
# Conditions: WB baseline + each requested level (as strings).
ENRICH_CONDS = ["wb"] + [str(lv) for lv in ENRICH_LEVELS]

ENRICH_SUBJECTS = sorted(set(inputs.input_zip_lists["T1w"]["subject"]))

# Path templates (plain sub-{subject}, mirroring mask_sweep; no session entity).
_TMP = config["tmp_dir"] + "/sub-{subject}/tracts/enrichment_sweep"
_TRACTS = "sub-{subject}/tracts/enrichment_sweep"
_QC = "sub-{subject}/qc/enrichment_sweep"

ROI_MAX_HEMI = _TMP + "/sub-{subject}_hemi-{hemi}_label-{seed}_desc-roimax_tractography.tck"
ROI_MAX_MERGED = _TMP + "/sub-{subject}_label-{seed}_desc-roimax_tractography.tck"
ROI_SUB = _TMP + "/sub-{subject}_label-{seed}_level-{level}_desc-roi_tractography.tck"
COMBINED = _TMP + "/sub-{subject}_label-{seed}_level-{level}_desc-combined_tractography.tck"

# WB baseline tractogram (read-only reuse of the core pipeline output).
WB_TCK = rules.wb_tckgen_merge.output.tck


wildcard_constraints:
    cond="wb|[0-9]+",
    level="[0-9]+",


def enrich_hemi_rng(hemi):
    base = int(config.get("mrtrix_rng_seed", 42))
    return base if hemi == "L" else base + 1000


def enrich_cond_tck(wc):
    """Tractogram for a condition: WB baseline for 'wb', else the combined tck."""
    if wc.cond == "wb":
        return WB_TCK.format(subject=wc.subject)
    return COMBINED.format(subject=wc.subject, seed=wc.seed, level=wc.cond)


def enrich_cond_roi_dep(wc):
    """ROI dependency for a condition (none for the WB baseline)."""
    if wc.cond == "wb":
        return []
    return ROI_SUB.format(subject=wc.subject, seed=wc.seed, level=wc.cond)


def enrich_roi_arg(wc):
    """ROI tractogram path passed to the metrics script ('' for WB baseline)."""
    if wc.cond == "wb":
        return ""
    return ROI_SUB.format(subject=wc.subject, seed=wc.seed, level=wc.cond)


# -----------------------------
# (1) Maximal-enrichment ROI tractogram  -- the ONLY ROI generation
# -----------------------------
rule enrich_roi_tckgen_max:
    """
    Voxelwise ROI-seeded tractography at the MAXIMAL enrichment level.
    Same tckgen algorithm/params as the core roi_tckgen; only the seeds-per-voxel
    differs (= max(levels)). Generated once per hemi, then merged.
    """
    input:
        wm_fod=get_fod_for_tracking,
        mask=bids(root=root, datatype="dwi", suffix="mask.mif", **subj_wildcards),
        seed_mask=bids(
            root=root, **subj_wildcards, hemi="{hemi}", label="{seed}",
            datatype="anat", suffix="mask.nii.gz",
        ),
        five_tt=rules.act_5ttgen.output.five_tt,
    params:
        spv=lambda wc: ENRICH_MAX_LEVEL,
        step=lambda wc: config.get("mrtrix", {}).get("step", 0.35),
        angle=lambda wc: config.get("mrtrix", {}).get("angle", 45),
        minlength=lambda wc: config.get("mrtrix", {}).get("minlength", 3),
        maxlength=lambda wc: config.get("mrtrix", {}).get("maxlength", 250),
        cutoff=lambda wc: config.get("mrtrix", {}).get("cutoff", 0.06),
        rng=lambda wc: enrich_hemi_rng(wc.hemi),
    output:
        tck=temp(ROI_MAX_HEMI),
        tckinfo=_QC + "/sub-{subject}_hemi-{hemi}_label-{seed}_desc-roimax_tckinfo.txt",
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_roi_tckgen_max.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_roi_tckgen_max.tsv"
    threads: lambda wc: res("enrich_roi_tckgen_max", "threads", 6)
    resources:
        mem_mb=lambda wc: res("enrich_roi_tckgen_max", "mem_mb", 32000),
        time=lambda wc: res("enrich_roi_tckgen_max", "time_min", 240),
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


rule enrich_roi_merge_max:
    """Merge per-hemi maximal ROI tractograms into the single maximal ROI tractogram."""
    input:
        tcks=lambda wc: expand(ROI_MAX_HEMI, hemi=ENRICH_HEMIS, allow_missing=True),
    output:
        tck=temp(ROI_MAX_MERGED),
        tckinfo=_QC + "/sub-{subject}_label-{seed}_desc-roimax_merged_tckinfo.txt",
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_roi_merge_max.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_roi_merge_max.tsv"
    threads: lambda wc: res("enrich_roi_merge_max", "threads", 2)
    resources:
        mem_mb=lambda wc: res("enrich_roi_merge_max", "mem_mb", 2000),
        time=lambda wc: res("enrich_roi_merge_max", "time_min", 15),
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

        tckedit -force {input.tcks} "{output.tck}" &> "{log}"
        tckinfo "{output.tck}" > "{output.tckinfo}" 2>> "{log}"
        """


# -----------------------------
# (2) Lower levels = random subsample ONLY (no re-seeding)
# -----------------------------
rule enrich_roi_subsample:
    """Reproducible random subsample of the maximal ROI tractogram to this level."""
    input:
        tck=rules.enrich_roi_merge_max.output.tck,
    output:
        tck=temp(ROI_SUB),
        count=_QC + "/sub-{subject}_label-{seed}_level-{level}_roi_count.txt",
    params:
        script=lambda wc: SUBSAMPLE_SCRIPT,
        max_level=lambda wc: ENRICH_MAX_LEVEL,
        subsample_seed=lambda wc: ENRICH_SUBSAMPLE_SEED,
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_level-{level}_roi_subsample.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_level-{level}_roi_subsample.tsv"
    threads: lambda wc: res("enrich_roi_subsample", "threads", 1)
    resources:
        mem_mb=lambda wc: res("enrich_roi_subsample", "mem_mb", 32000),
        time=lambda wc: res("enrich_roi_subsample", "time_min", 60),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tck}")"
        mkdir -p "$(dirname "{output.count}")"
        mkdir -p "$(dirname "{log}")"

        python "{params.script}" \
          --in-tck "{input.tck}" \
          --out-tck "{output.tck}" \
          --level {wildcards.level} \
          --max-level {params.max_level} \
          --seed {params.subsample_seed} \
          --out-count "{output.count}" \
          &> "{log}"
        """


# -----------------------------
# (3) Per-level tractogram = fixed WB + subsampled ROI (no recompute)
# -----------------------------
rule enrich_combine:
    input:
        wb=rules.wb_tckgen_merge.output.tck,
        roi=rules.enrich_roi_subsample.output.tck,
    output:
        tck=temp(COMBINED),
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_level-{level}_combine.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_level-{level}_combine.tsv"
    threads: lambda wc: res("enrich_combine", "threads", 2)
    resources:
        mem_mb=lambda wc: res("enrich_combine", "mem_mb", 4000),
        time=lambda wc: res("enrich_combine", "time_min", 30),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tck}")"
        mkdir -p "$(dirname "{log}")"

        tckedit -force "{input.wb}" "{input.roi}" "{output.tck}" &> "{log}"
        """


# -----------------------------
# (4) Safety check (UNCHANGED tckinfo), once on WB baseline + once per condition
# -----------------------------
rule enrich_safety_check:
    input:
        tck=enrich_cond_tck,
    output:
        tckinfo=_QC + "/sub-{subject}_label-{seed}_cond-{cond}_desc-enrich_tckinfo.txt",
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_cond-{cond}_safety_check.log",
    threads: lambda wc: res("enrich_safety_check", "threads", 1)
    resources:
        mem_mb=lambda wc: res("enrich_safety_check", "mem_mb", 2000),
        time=lambda wc: res("enrich_safety_check", "time_min", 15),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tckinfo}")"
        mkdir -p "$(dirname "{log}")"

        tckinfo "{input.tck}" > "{output.tckinfo}" 2> "{log}"
        """


# -----------------------------
# (5) SIFT2 per condition (identical params; logs mu/cost/weights)
# -----------------------------
rule enrich_sift2:
    input:
        tck=enrich_cond_tck,
        wm_fod=get_fod_for_tracking,
        five_tt=rules.act_5ttgen.output.five_tt,
    output:
        weights=_TRACTS + "/sub-{subject}_label-{seed}_cond-{cond}_desc-enrich_sift2_weights.txt",
        mu=_QC + "/sub-{subject}_label-{seed}_cond-{cond}_sift2_mu.txt",
        sift2_csv=_QC + "/sub-{subject}_label-{seed}_cond-{cond}_sift2_stats.csv",
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_cond-{cond}_sift2.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_cond-{cond}_sift2.tsv"
    threads: lambda wc: res("enrich_sift2", "threads", 10)
    resources:
        mem_mb=lambda wc: res("enrich_sift2", "mem_mb", 40000),
        time=lambda wc: res("enrich_sift2", "time_min", 300),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.weights}")"
        mkdir -p "$(dirname "{output.mu}")"
        mkdir -p "$(dirname "{log}")"

        tcksift2 -nthreads {threads} -act "{input.five_tt}" \
          -out_mu "{output.mu}" \
          -csv "{output.sift2_csv}" \
          "{input.tck}" "{input.wm_fod}" "{output.weights}" \
          &> "{log}"
        """


# -----------------------------
# (6) Voxelwise connectivity fingerprint per condition (same logic as main pipeline)
#     Reuses connectivity.smk: reslice_seed_to_nodegrid (filter mask) and
#     build_seed_nodes (per-voxel node parcellation) -- no duplicate connectome logic.
# -----------------------------
rule enrich_filter_tractogram:
    """Filter a condition's tractogram by the (node-grid) seed, carrying its SIFT2 weights."""
    input:
        tractogram=enrich_cond_tck,
        sift2_weights=rules.enrich_sift2.output.weights,
        seed_mask=rules.reslice_seed_to_nodegrid.output.mask,
    output:
        tractogram=temp(
            _TMP + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_desc-seedfiltered_tractography.tck"
        ),
        sift2_weights=temp(
            _TMP + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_desc-seedfiltered_sift2_weights.txt"
        ),
        tckinfo=(
            _QC + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_desc-seedfiltered_tckinfo.txt"
        ),
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_filter_tractogram.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_filter_tractogram.tsv"
    threads: lambda wc: res("enrich_filter_tractogram", "threads", 2)
    resources:
        mem_mb=lambda wc: res("enrich_filter_tractogram", "mem_mb", 16000),
        time=lambda wc: res("enrich_filter_tractogram", "time_min", 60),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
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


rule enrich_voxelwise_connectivity:
    """
    Voxelwise seed-to-target connectivity for one condition, identical to the main
    pipeline: tck2connectome on the reused per-voxel node image (build_seed_nodes),
    then slice the [seed voxel x atlas target] block.
    """
    input:
        tractogram=rules.enrich_filter_tractogram.output.tractogram,
        sift2_weights=rules.enrich_filter_tractogram.output.sift2_weights,
        nodes=rules.build_seed_nodes.output.nodes,
        seed_voxel_index=rules.build_seed_nodes.output.seed_voxel_index,
    output:
        connectivity_matrix=(
            _TRACTS + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_desc-{targets}_connectivity_matrix.csv"
        ),
        assignments=(
            _QC + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_desc-{targets}_method-mrtrix_assignments.txt"
        ),
        qc_metrics=(
            _QC + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_desc-{targets}_connectivity_qc.json"
        ),
        connectome_full=temp(
            _TMP + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_desc-{targets}_connectome_full.txt"
        ),
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_desc-{targets}_voxelwise_connectivity.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_desc-{targets}_voxelwise_connectivity.tsv"
    params:
        slice_script=lambda wc: SLICE_BLOCK_SCRIPT,
        target_labels=lambda wc: ",".join(config["targets"][wc.targets]["labels"]),
        target_search_radius=lambda wc: config.get("mrtrix", {}).get("target_search_radius", 4),
    threads: lambda wc: res("enrich_voxelwise_connectivity", "threads", 2)
    resources:
        mem_mb=lambda wc: res("enrich_voxelwise_connectivity", "mem_mb", 8000),
        time=lambda wc: res("enrich_voxelwise_connectivity", "time_min", 30),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.connectivity_matrix}")"
        mkdir -p "$(dirname "{output.assignments}")"
        mkdir -p "$(dirname "{output.connectome_full}")"
        mkdir -p "$(dirname "{log}")"

        tck2connectome -nthreads {threads} -force -symmetric \
          -assignment_radial_search {params.target_search_radius} \
          -tck_weights_in "{input.sift2_weights}" \
          -out_assignments "{output.assignments}" \
          "{input.tractogram}" "{input.nodes}" "{output.connectome_full}" \
          &> "{log}"

        python "{params.slice_script}" \
          --connectome "{output.connectome_full}" \
          --voxel-index "{input.seed_voxel_index}" \
          --header "{params.target_labels}" \
          --out-matrix "{output.connectivity_matrix}" \
          --out-qc "{output.qc_metrics}" \
          &>> "{log}"
        """


# -----------------------------
# (6b) OPTIONAL volume-bias-corrected connectivity fingerprint
#      Identical to enrich_voxelwise_connectivity, but adds MRtrix's native
#      `-scale_invnodevol`, which scales every connectome edge by the inverse of
#      its two node volumes (2 / (vol_i + vol_j), volumes in voxels). Seed nodes
#      are single voxels (vol = 1), so for a seed-voxel x target edge this is
#      2 / (1 + vol_target) -- i.e. it normalizes for target-parcel volume.
#
#      Fully additive: emitted ONLY when enrichment_sweep.volume_bias_correction
#      is true (requested via the entry point / all_enrichment_sweep). The
#      `_scale-invnodevol` tag sits between the constrained `cond` and `desc`
#      entities so its path can never collide with the uncorrected rule above.
#      Nothing upstream or in connectivity.smk is touched.
# -----------------------------
rule enrich_voxelwise_connectivity_invnodevol:
    input:
        tractogram=rules.enrich_filter_tractogram.output.tractogram,
        sift2_weights=rules.enrich_filter_tractogram.output.sift2_weights,
        nodes=rules.build_seed_nodes.output.nodes,
        seed_voxel_index=rules.build_seed_nodes.output.seed_voxel_index,
    output:
        connectivity_matrix=(
            _TRACTS + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_scale-invnodevol_desc-{targets}_connectivity_matrix.csv"
        ),
        qc_metrics=(
            _QC + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_scale-invnodevol_desc-{targets}_connectivity_qc.json"
        ),
        connectome_full=temp(
            _TMP + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_scale-invnodevol_desc-{targets}_connectome_full.txt"
        ),
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_scale-invnodevol_desc-{targets}_voxelwise_connectivity.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_scale-invnodevol_desc-{targets}_voxelwise_connectivity.tsv"
    params:
        slice_script=lambda wc: SLICE_BLOCK_SCRIPT,
        target_labels=lambda wc: ",".join(config["targets"][wc.targets]["labels"]),
        target_search_radius=lambda wc: config.get("mrtrix", {}).get("target_search_radius", 4),
    threads: lambda wc: res("enrich_voxelwise_connectivity_invnodevol", "threads", 2)
    resources:
        mem_mb=lambda wc: res("enrich_voxelwise_connectivity_invnodevol", "mem_mb", 8000),
        time=lambda wc: res("enrich_voxelwise_connectivity_invnodevol", "time_min", 30),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.connectivity_matrix}")"
        mkdir -p "$(dirname "{output.qc_metrics}")"
        mkdir -p "$(dirname "{output.connectome_full}")"
        mkdir -p "$(dirname "{log}")"

        tck2connectome -nthreads {threads} -force -symmetric \
          -scale_invnodevol \
          -assignment_radial_search {params.target_search_radius} \
          -tck_weights_in "{input.sift2_weights}" \
          "{input.tractogram}" "{input.nodes}" "{output.connectome_full}" \
          &> "{log}"

        python "{params.slice_script}" \
          --connectome "{output.connectome_full}" \
          --voxel-index "{input.seed_voxel_index}" \
          --header "{params.target_labels}" \
          --out-matrix "{output.connectivity_matrix}" \
          --out-qc "{output.qc_metrics}" \
          &>> "{log}"
        """


# -----------------------------
# (7) Per-condition metrics row (logging only)
# -----------------------------
rule enrich_condition_row:
    input:
        mu=rules.enrich_sift2.output.mu,
        sift2_csv=rules.enrich_sift2.output.sift2_csv,
        safety=rules.enrich_safety_check.output.tckinfo,
        wb=rules.wb_tckgen_merge.output.tck,
        roi=enrich_cond_roi_dep,
    output:
        row=_QC + "/sub-{subject}_label-{seed}_cond-{cond}_metrics.csv",
    params:
        script=lambda wc: METRICS_SCRIPT,
        level_label=lambda wc: "baseline" if wc.cond == "wb" else wc.cond,
        roi_arg=enrich_roi_arg,
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_cond-{cond}_metrics.log",
    threads: 1
    resources:
        mem_mb=lambda wc: res("enrich_condition_row", "mem_mb", 2000),
        time=lambda wc: res("enrich_condition_row", "time_min", 10),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.row}")"
        mkdir -p "$(dirname "{log}")"

        python "{params.script}" \
          --subject "{wildcards.subject}" \
          --condition "{wildcards.cond}" \
          --enrichment-level "{params.level_label}" \
          --wb-tck "{input.wb}" \
          --roi-tck "{params.roi_arg}" \
          --mu "{input.mu}" \
          --sift2-csv "{input.sift2_csv}" \
          --out-row "{output.row}" \
          &> "{log}"
        """


# -----------------------------
# (8) Per-subject summary (ENTRY target) — concatenate condition rows
# -----------------------------
rule enrich_summary:
    input:
        rows=lambda wc: expand(
            _QC + "/sub-{subject}_label-{seed}_cond-{cond}_metrics.csv",
            subject=wc.subject, seed=wc.seed, cond=ENRICH_CONDS,
        ),
    output:
        summary=_QC + "/sub-{subject}_label-{seed}_enrichment_sweep_summary.csv",
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_summary.log",
    threads: 1
    resources:
        mem_mb=lambda wc: res("enrich_summary", "mem_mb", 1000),
        time=lambda wc: res("enrich_summary", "time_min", 10),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.summary}")"
        mkdir -p "$(dirname "{log}")"

        set -- {input.rows}
        head -n 1 "$1" > "{output.summary}"
        for f in "$@"; do
          tail -n +2 "$f" >> "{output.summary}"
        done
        """


# -----------------------------
# Convenience aggregate target (also wired into `rule all` via Snakefile helper)
# -----------------------------
rule all_enrichment_sweep:
    input:
        (
            expand(
                _QC + "/sub-{subject}_label-{seed}_enrichment_sweep_summary.csv",
                subject=ENRICH_SUBJECTS, seed=ENRICH_SEED,
            )
            + expand(
                _TRACTS + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
                "_desc-{targets}_connectivity_matrix.csv",
                subject=ENRICH_SUBJECTS, seed=ENRICH_SEED, hemi=ENRICH_HEMIS,
                cond=ENRICH_CONDS, targets=ENRICH_TARGETS,
            )
            + (
                expand(
                    _TRACTS + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
                    "_scale-invnodevol_desc-{targets}_connectivity_matrix.csv",
                    subject=ENRICH_SUBJECTS, seed=ENRICH_SEED, hemi=ENRICH_HEMIS,
                    cond=ENRICH_CONDS, targets=ENRICH_TARGETS,
                )
                if ENRICH_VOLNORM else []
            )
        ) if ENRICH_ENABLED else [],
