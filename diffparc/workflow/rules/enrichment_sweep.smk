# enrichment_sweep.smk
# Experimental, fully ISOLATED ROI-enrichment sweep (additive; mirrors mask_sweep.smk).
#
# Sweeps ROI seeding strength ("enrichment level") WITHOUT changing any core rule:
#   1. Generate ONE maximal-enrichment ROI tractogram (voxelwise ROI seeding at
#      spv = max(levels)). This is the ONLY ROI tractogram generation.
#   2. Lower levels = reproducible NESTED random subsample of that maximal ROI
#      tractogram (no re-seeding, no re-running tckgen). Nested = each higher
#      level is a strict superset of the lower ones (1000 ⊂ 3000 ⊂ 5000).
#   3. WB baseline = ONE whole-brain tractogram (no ROI), reused read-only from
#      the core pipeline (rules.wb_tckgen_merge). This is the reference condition.
#   4. Per enrichment level: combine fixed WB + subsampled ROI (no recompute).
#   5. SIFT2, default branch ("refit"): run independently per condition (WB
#      baseline + each level), identical params; logs mu, cost, weights.
#      Optional second branch ("propagated", config sift2_propagated): fit SIFT2
#      ONCE on the maximal combined tractogram and carry those weights down to
#      lower levels via the SAME nested subsample (streamlines + weights), so the
#      two branches differ ONLY in the SIFT2 weighting, never in the streamlines.
#   6. Safety check (tckinfo) applied unchanged, once on WB baseline + once/condition.
#
# Only entry points into the rest of the workflow: this include line and the
# get_enrichment_sweep_outputs() helper added to `rule all`. Off by default.

import os


def res(rule, key, default):
    return config.get("resources", {}).get(rule, {}).get(key, default)


SUBSAMPLE_SCRIPT = os.path.join(workflow.basedir, "scripts", "subsample_tractogram.py")
METRICS_SCRIPT = os.path.join(workflow.basedir, "scripts", "enrichment_condition_metrics.py")
PROPAGATE_SCRIPT = os.path.join(workflow.basedir, "scripts", "propagate_sift2_weights.py")
# Reused from the main connectivity pipeline (connectivity.smk) -- no duplicate logic.
SLICE_BLOCK_SCRIPT = os.path.join(workflow.basedir, "scripts", "slice_connectome_block.py")
# Same per-voxel node builder the main pipeline uses; called here on the selected
# max-volume sweep mask so node building stays byte-for-byte the same logic.
BUILD_NODES_SCRIPT = os.path.join(workflow.basedir, "scripts", "build_voxel_nodes.py")

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
# Optional: also emit the "propagated" SIFT2 branch (fit once on the maximal
# combined tractogram, carry weights down the nested subsample) alongside the
# default per-level "refit" branch.
ENRICH_SIFT2_PROP = bool(ENRICH_CFG.get("sift2_propagated", False))
# Conditions: WB baseline + each requested level (as strings).
ENRICH_CONDS = ["wb"] + [str(lv) for lv in ENRICH_LEVELS]

ENRICH_SUBJECTS = sorted(set(inputs.input_zip_lists["T1w"]["subject"]))

# -----------------------------------------------------------------------------
# Seed = LARGEST-volume mask_sweep mask (added)
# -----------------------------------------------------------------------------
# Instead of the canonical thresholded+dilated seed, the enrichment sweep seeds
# from the mask with the most voxels among the mask_sweep candidates, restricted
# to the NATIVE-source masks of ENRICH_SEED across all configured sweep
# thresholds. The selection is per (subject, hemi), purely by nonzero voxel
# count, deterministic with a sorted-filename tie-break (no RNG). The SAME
# selected mask drives the whole enrichment chain (seeding + filtering + node
# parcellation), so the seeded and analysed regions stay voxel-consistent.
# These candidate masks are produced by mask_sweep.smk:sweep_binarize_seed, whose
# rule is always defined; requesting them here pulls them into the DAG regardless
# of the mask_sweep.enabled (rule all) gate. mask_sweep.smk / connectivity.smk
# are left untouched.
ENRICH_SWEEP_MASK_SOURCE = "native"
ENRICH_SWEEP_THR_TAGS = [
    str(t).replace(".", "p")
    for t in config.get("mask_sweep", {}).get("thresholds", [])
]
if ENRICH_ENABLED and not ENRICH_SWEEP_THR_TAGS:
    raise ValueError(
        "enrichment_sweep is enabled and now selects its seed from the mask_sweep "
        "candidate masks, but config['mask_sweep']['thresholds'] is empty/missing. "
        "Define the sweep thresholds (the candidate masks) or disable enrichment_sweep."
    )


def enrich_sweep_mask_candidates(wc):
    """Native-source ENRICH_SEED sweep masks (one per configured threshold) for
    this (subject, hemi), in deterministic sorted-filename order. These are the
    candidates the max-volume selection chooses between."""
    return sorted(
        f"sub-{wc.subject}/anat/mask_sweep/{wc.seed}/{ENRICH_SWEEP_MASK_SOURCE}/"
        f"sub-{wc.subject}_hemi-{wc.hemi}_label-{wc.seed}_thr-{tag}_mask.nii.gz"
        for tag in ENRICH_SWEEP_THR_TAGS
    )

# Path templates (plain sub-{subject}, mirroring mask_sweep; no session entity).
_TMP = config["tmp_dir"] + "/sub-{subject}/tracts/enrichment_sweep"
_TRACTS = "sub-{subject}/tracts/enrichment_sweep"
_QC = "sub-{subject}/qc/enrichment_sweep"

ROI_MAX_HEMI = _TMP + "/sub-{subject}_hemi-{hemi}_label-{seed}_desc-roimax_tractography.tck"
ROI_MAX_MERGED = _TMP + "/sub-{subject}_label-{seed}_desc-roimax_tractography.tck"
ROI_SUB = _TMP + "/sub-{subject}_label-{seed}_level-{level}_desc-roi_tractography.tck"
COMBINED = _TMP + "/sub-{subject}_label-{seed}_level-{level}_desc-combined_tractography.tck"

# -----------------------------------------------------------------------------
# Output directory layout (ORGANIZATIONAL ONLY -- filenames/computations unchanged).
# Persistent outputs are grouped by enrichment level (cond-{cond}/), and within
# QC by type, so each level's results are self-contained and easy to compare:
#   tracts/enrichment_sweep/cond-{cond}/            primary (mu-scaled) connectomes + SIFT2 weights
#   qc/enrichment_sweep/cond-{cond}/connectome/     raw connectomes, qc json, assignments
#   qc/enrichment_sweep/cond-{cond}/sift2/          mu, sift2 stats
#   qc/enrichment_sweep/cond-{cond}/tractogram/     tckinfo + visual QC (tdi/endpoints/decmap/subset/lengths)
#   qc/enrichment_sweep/cond-{cond}/                per-condition metrics row
#   qc/enrichment_sweep/roi_max/                    level-independent maximal-ROI tractography QC
#   qc/enrichment_sweep/<summary.csv>               cross-level summary (top level)
_T_COND = _TRACTS + "/cond-{cond}"
_Q_COND = _QC + "/cond-{cond}"
_Q_CONN = _Q_COND + "/connectome"
_Q_SIFT2 = _Q_COND + "/sift2"
_Q_TRACT = _Q_COND + "/tractogram"
_Q_ROIMAX = _QC + "/roi_max"

# Max-volume sweep-mask selection (added). Selected mask + its selection report
# live in QC (auditable, persistent); the node-grid reslice and per-voxel node
# image are scratch/temp like their main-pipeline counterparts.
SEL_MASK = _QC + "/maxvol_mask/sub-{subject}_hemi-{hemi}_label-{seed}_desc-maxvolsweep_mask.nii.gz"
SEL_REPORT = _QC + "/maxvol_mask/sub-{subject}_hemi-{hemi}_label-{seed}_desc-maxvolsweep_selection.csv"
SEL_NODEGRID = _TMP + "/sub-{subject}_hemi-{hemi}_label-{seed}_desc-maxvolsweep_res-nodegrid_mask.nii.gz"
SEL_NODES = _TMP + "/sub-{subject}_hemi-{hemi}_label-{seed}_maxvolsweep_{targets}_nodes.nii.gz"
SEL_VOXIDX = _QC + "/maxvol_mask/sub-{subject}_hemi-{hemi}_label-{seed}_maxvolsweep_{targets}_seed_voxel_index.csv"

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


# --- Propagated-SIFT2 branch helpers (only used when sift2_propagated is set) ---
def enrich_weights_max(wc):
    """SIFT2 weights from the maximal combined tractogram -- the single fit the
    propagated branch reuses for every level."""
    return rules.enrich_sift2.output.weights.format(
        subject=wc.subject, seed=wc.seed, cond=str(ENRICH_MAX_LEVEL)
    )


def enrich_mu_max(wc):
    """SIFT2 mu (proportionality coefficient) from the MAXIMAL-level fit. This is
    the calibration scalar for the propagated branch, which has no per-condition
    SIFT2 run of its own: its weights are the max fit carried down, so mu_max puts
    the propagated connectome in the SAME fibre-density units as the refit branch,
    making the no-refit/subsample bias directly comparable."""
    return rules.enrich_sift2.output.mu.format(
        subject=wc.subject, seed=wc.seed, cond=str(ENRICH_MAX_LEVEL)
    )


def enrich_roi_max_meta(wc):
    """Subsample meta at the maximal level (supplies n_total = ROI_max count)."""
    return rules.enrich_roi_subsample.output.meta.format(
        subject=wc.subject, seed=wc.seed, level=str(ENRICH_MAX_LEVEL)
    )


def enrich_cond_indices(wc):
    """This level's subsample meta (nested ROI indices); none for the WB baseline."""
    if wc.cond == "wb":
        return []
    return rules.enrich_roi_subsample.output.meta.format(
        subject=wc.subject, seed=wc.seed, level=wc.cond
    )


def enrich_indices_path(wc):
    """This level's subsample meta path, or '' for the WB baseline (the script
    treats an empty --indices as 'no ROI streamlines')."""
    if wc.cond == "wb":
        return ""
    return rules.enrich_roi_subsample.output.meta.format(
        subject=wc.subject, seed=wc.seed, level=wc.cond
    )


# -----------------------------
# (0) Select the largest-volume mask_sweep mask  (drives the whole chain)
# -----------------------------
rule enrich_select_max_volume_mask:
    """Pick the largest-volume (nonzero voxel count) mask_sweep mask for this
    (subject, hemi) from the native-source ENRICH_SEED candidates.

    Deterministic: candidates arrive in sorted-filename order and the strict
    maximum is kept, so a tie falls to the first sorted filename -- no RNG.
    Selection is by voxel count (mrstats -output count over the mask's own
    nonzero voxels), never file size or header metadata. The per-candidate counts
    and the winner are written to a CSV report for provenance/QC."""
    input:
        candidates=enrich_sweep_mask_candidates,
    output:
        mask=SEL_MASK,
        report=SEL_REPORT,
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_select_max_volume_mask.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_select_max_volume_mask.tsv"
    threads: 1
    resources:
        mem_mb=lambda wc: res("enrich_select_max_volume_mask", "mem_mb", 2000),
        time=lambda wc: res("enrich_select_max_volume_mask", "time_min", 10),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.mask}")"
        mkdir -p "$(dirname "{log}")"
        : > "{log}"

        best=""
        best_n=-1
        printf 'candidate,voxel_count\n' > "{output.report}"
        # Candidates are in sorted-filename order; strict '>' keeps the FIRST
        # maximum => deterministic lexicographic tie-break. Count = the mask's own
        # nonzero voxels (mrstats restricted to the mask itself).
        for m in {input.candidates}; do
          n="$(mrstats "$m" -mask "$m" -output count 2>> "{log}" | tr -d '[:space:]')"
          n="${{n:-0}}"
          printf '%s,%s\n' "$m" "$n" >> "{output.report}"
          if [ "$n" -gt "$best_n" ]; then
            best_n="$n"
            best="$m"
          fi
        done

        if [ -z "$best" ]; then
          echo "enrich_select_max_volume_mask: no candidate masks found" >> "{log}"
          exit 1
        fi

        printf 'selected,%s,%s\n' "$best" "$best_n" >> "{output.report}"
        echo "enrich_select_max_volume_mask: selected $best (voxel_count=$best_n)" >> "{log}"
        cp -f "$best" "{output.mask}"
        """


# -----------------------------
# (0b) Reslice selected mask onto the node grid (== target-atlas grid)
#      Mirrors connectivity.smk:reslice_seed_to_nodegrid, on the selected mask.
#      The single resliced mask feeds BOTH filtering and node building below.
# -----------------------------
rule enrich_reslice_selected_to_nodegrid:
    input:
        seed_mask=rules.enrich_select_max_volume_mask.output.mask,
        ref=bids(
            root=root, suffix="mask.nii.gz", desc="brain", space="T1w",
            res="upsampled", datatype="dwi", **subj_wildcards,
        ),
    output:
        mask=temp(SEL_NODEGRID),
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_reslice_selected_to_nodegrid.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_reslice_selected_to_nodegrid.tsv"
    threads: 1
    resources:
        mem_mb=lambda wc: res("enrich_reslice_selected_to_nodegrid", "mem_mb", 4000),
        time=lambda wc: res("enrich_reslice_selected_to_nodegrid", "time_min", 10),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.mask}")"
        mkdir -p "$(dirname "{log}")"

        c3d -int 0 "{input.ref}" "{input.seed_mask}" \
          -reslice-identity \
          -threshold 0.5 inf 1 0 \
          -type uchar \
          -o "{output.mask}" \
          &> "{log}"
        """


# -----------------------------
# (0c) Per-voxel node parcellation from the selected mask
#      Same build_voxel_nodes.py logic as connectivity.smk:build_seed_nodes.
# -----------------------------
rule enrich_build_selected_nodes:
    input:
        seed_mask=rules.enrich_reslice_selected_to_nodegrid.output.mask,
        targets_dseg="sub-{subject}/anat/sub-{subject}_desc-{targets}_dseg.nii.gz",
    output:
        nodes=temp(SEL_NODES),
        seed_voxel_index=SEL_VOXIDX,
    params:
        script=lambda wc: BUILD_NODES_SCRIPT,
        n_targets=lambda wc: len(config["targets"][wc.targets]["labels"]),
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}_build_selected_nodes.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}_build_selected_nodes.tsv"
    threads: 1
    resources:
        mem_mb=lambda wc: res("enrich_build_selected_nodes", "mem_mb", 4000),
        time=lambda wc: res("enrich_build_selected_nodes", "time_min", 10),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.nodes}")"
        mkdir -p "$(dirname "{output.seed_voxel_index}")"
        mkdir -p "$(dirname "{log}")"

        python "{params.script}" \
          --seed-mask "{input.seed_mask}" \
          --targets-nii "{input.targets_dseg}" \
          --n-targets {params.n_targets} \
          --out-nodes "{output.nodes}" \
          --out-voxel-index "{output.seed_voxel_index}" \
          &> "{log}"
        """


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
        seed_mask=rules.enrich_select_max_volume_mask.output.mask,
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
        tckinfo=_Q_ROIMAX + "/sub-{subject}_hemi-{hemi}_label-{seed}_desc-roimax_tckinfo.txt",
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_roi_tckgen_max.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_roi_tckgen_max.tsv"
    threads: lambda wc: res("enrich_roi_tckgen_max", "threads", 8)
    resources:
        mem_mb=lambda wc: res("enrich_roi_tckgen_max", "mem_mb", 4000),
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
        tckinfo=_Q_ROIMAX + "/sub-{subject}_label-{seed}_desc-roimax_merged_tckinfo.txt",
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
    """Reproducible NESTED random subsample of the maximal ROI tractogram to this level.

    The `meta` JSON records n_total (ROI_max count) and the selected indices; it
    lets the optional propagated-SIFT2 branch carry max-fit weights down to this
    level without re-fitting (scripts/propagate_sift2_weights.py).
    """
    input:
        tck=rules.enrich_roi_merge_max.output.tck,
    output:
        tck=temp(ROI_SUB),
        roi_count=_QC + "/cond-{level}/tractogram/sub-{subject}_label-{seed}_level-{level}_roi_count.txt",
        meta=temp(
            _TMP + "/sub-{subject}_label-{seed}_level-{level}_desc-roi_subsample_meta.json"
        ),
    params:
        script=lambda wc: SUBSAMPLE_SCRIPT,
        max_level=lambda wc: ENRICH_MAX_LEVEL,
        subsample_seed=lambda wc: ENRICH_SUBSAMPLE_SEED,
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_level-{level}_roi_subsample.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_level-{level}_roi_subsample.tsv"
    threads: lambda wc: res("enrich_roi_subsample", "threads", 4)
    resources:
        mem_mb=lambda wc: res("enrich_roi_subsample", "mem_mb", 24000),
        time=lambda wc: res("enrich_roi_subsample", "time_min", 60),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tck}")"
        mkdir -p "$(dirname "{output.roi_count}")"
        mkdir -p "$(dirname "{log}")"

        python "{params.script}" \
          --in-tck "{input.tck}" \
          --out-tck "{output.tck}" \
          --level {wildcards.level} \
          --max-level {params.max_level} \
          --seed {params.subsample_seed} \
          --out-count "{output.roi_count}" \
          --out-meta "{output.meta}" \
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
        tckinfo=_Q_TRACT + "/sub-{subject}_label-{seed}_cond-{cond}_desc-enrich_tckinfo.txt",
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_cond-{cond}_safety_check.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_cond-{cond}_safety_check.tsv"
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
        weights=_T_COND + "/sub-{subject}_label-{seed}_cond-{cond}_desc-enrich_sift2_weights.txt",
        mu=_Q_SIFT2 + "/sub-{subject}_label-{seed}_cond-{cond}_sift2_mu.txt",
        sift2_csv=_Q_SIFT2 + "/sub-{subject}_label-{seed}_cond-{cond}_sift2_stats.csv",
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
        seed_mask=rules.enrich_reslice_selected_to_nodegrid.output.mask,
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
            _Q_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
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
        nodes=rules.enrich_build_selected_nodes.output.nodes,
        seed_voxel_index=rules.enrich_build_selected_nodes.output.seed_voxel_index,
        mu=rules.enrich_sift2.output.mu,
    output:
        connectivity_matrix=(
            _T_COND + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_desc-{targets}_connectivity_matrix.csv"
        ),
        connectivity_matrix_raw=(
            _Q_CONN + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_desc-{targets}_meas-raw_connectivity_matrix.csv"
        ),
        assignments=(
            _Q_CONN + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_desc-{targets}_method-mrtrix_assignments.txt"
        ),
        qc_metrics=(
            _Q_CONN + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
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

        if grep -q '[0-9]' "{input.sift2_weights}"; then
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
            --out-matrix-raw "{output.connectivity_matrix_raw}" \
            --mu-file "{input.mu}" \
            --out-qc "{output.qc_metrics}" \
            &>> "{log}"
        else
          echo "empty seed-filtered tractogram (0 streamlines): writing zero connectome" > "{log}"
          : > "{output.connectome_full}"
          : > "{output.assignments}"
          python "{params.slice_script}" \
            --connectome "{output.connectome_full}" \
            --voxel-index "{input.seed_voxel_index}" \
            --header "{params.target_labels}" \
            --out-matrix "{output.connectivity_matrix}" \
            --out-matrix-raw "{output.connectivity_matrix_raw}" \
            --mu-file "{input.mu}" \
            --out-qc "{output.qc_metrics}" \
            --empty \
            &>> "{log}"
        fi
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
        nodes=rules.enrich_build_selected_nodes.output.nodes,
        seed_voxel_index=rules.enrich_build_selected_nodes.output.seed_voxel_index,
        mu=rules.enrich_sift2.output.mu,
    output:
        connectivity_matrix=(
            _T_COND + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_scale-invnodevol_desc-{targets}_connectivity_matrix.csv"
        ),
        connectivity_matrix_raw=(
            _Q_CONN + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_scale-invnodevol_desc-{targets}_meas-raw_connectivity_matrix.csv"
        ),
        qc_metrics=(
            _Q_CONN + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
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

        if grep -q '[0-9]' "{input.sift2_weights}"; then
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
            --out-matrix-raw "{output.connectivity_matrix_raw}" \
            --mu-file "{input.mu}" \
            --out-qc "{output.qc_metrics}" \
            &>> "{log}"
        else
          echo "empty seed-filtered tractogram (0 streamlines): writing zero connectome" > "{log}"
          : > "{output.connectome_full}"
          python "{params.slice_script}" \
            --connectome "{output.connectome_full}" \
            --voxel-index "{input.seed_voxel_index}" \
            --header "{params.target_labels}" \
            --out-matrix "{output.connectivity_matrix}" \
            --out-matrix-raw "{output.connectivity_matrix_raw}" \
            --mu-file "{input.mu}" \
            --out-qc "{output.qc_metrics}" \
            --empty \
            &>> "{log}"
        fi
        """


# -----------------------------
# (6c) OPTIONAL propagated-SIFT2 branch -- per-condition propagated weights.
#      No SIFT2 re-fit: reuse the maximal-level SIFT2 weights (fit on WB ++ ROI_max)
#      and reindex them onto this condition's tractogram --
#        cond=wb    -> WB block only
#        cond=level -> WB block ++ ROI_max weights at the nested-subset indices
#      The result aligns 1:1 with enrich_cond_tck (WB, or WB ++ ROI_sub(level)),
#      which is the SAME tractogram the refit branch uses -- so the two branches
#      differ ONLY in the SIFT2 weights. Built only when sift2_propagated is set.
# -----------------------------
rule enrich_propagate_weights:
    input:
        weights_max=enrich_weights_max,
        roi_max_meta=enrich_roi_max_meta,
        indices=enrich_cond_indices,
    output:
        weights=temp(
            _TMP + "/sub-{subject}_label-{seed}_cond-{cond}_desc-propagated_sift2_weights.txt"
        ),
    params:
        script=lambda wc: PROPAGATE_SCRIPT,
        indices_path=enrich_indices_path,
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_cond-{cond}_propagate_weights.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_cond-{cond}_propagate_weights.tsv"
    threads: 1
    resources:
        mem_mb=lambda wc: res("enrich_propagate_weights", "mem_mb", 4000),
        time=lambda wc: res("enrich_propagate_weights", "time_min", 15),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.weights}")"
        mkdir -p "$(dirname "{log}")"

        python "{params.script}" \
          --weights-max "{input.weights_max}" \
          --roi-max-meta "{input.roi_max_meta}" \
          --indices "{params.indices_path}" \
          --out-weights "{output.weights}" \
          &> "{log}"
        """


# -----------------------------
# (6d) Propagated branch: filter by the seed (node grid), carrying propagated weights.
#      Identical to enrich_filter_tractogram except for the weights and the
#      `_sift2-propagated` tag.
# -----------------------------
rule enrich_filter_tractogram_propagated:
    input:
        tractogram=enrich_cond_tck,
        sift2_weights=rules.enrich_propagate_weights.output.weights,
        seed_mask=rules.enrich_reslice_selected_to_nodegrid.output.mask,
    output:
        tractogram=temp(
            _TMP + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_sift2-propagated_desc-seedfiltered_tractography.tck"
        ),
        sift2_weights=temp(
            _TMP + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_sift2-propagated_desc-seedfiltered_sift2_weights.txt"
        ),
        tckinfo=(
            _Q_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_sift2-propagated_desc-seedfiltered_tckinfo.txt"
        ),
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_sift2-propagated_filter_tractogram.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_sift2-propagated_filter_tractogram.tsv"
    threads: lambda wc: res("enrich_filter_tractogram_propagated", "threads", 2)
    resources:
        mem_mb=lambda wc: res("enrich_filter_tractogram_propagated", "mem_mb", 16000),
        time=lambda wc: res("enrich_filter_tractogram_propagated", "time_min", 60),
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


# -----------------------------
# (6e) Propagated branch: voxelwise connectivity fingerprint (uncorrected).
#      Same connectome logic as enrich_voxelwise_connectivity; `_sift2-propagated`
#      tag keeps its paths disjoint from the refit branch.
# -----------------------------
rule enrich_voxelwise_connectivity_propagated:
    input:
        tractogram=rules.enrich_filter_tractogram_propagated.output.tractogram,
        sift2_weights=rules.enrich_filter_tractogram_propagated.output.sift2_weights,
        nodes=rules.enrich_build_selected_nodes.output.nodes,
        seed_voxel_index=rules.enrich_build_selected_nodes.output.seed_voxel_index,
        mu=enrich_mu_max,
    output:
        connectivity_matrix=(
            _T_COND + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_sift2-propagated_desc-{targets}_connectivity_matrix.csv"
        ),
        connectivity_matrix_raw=(
            _Q_CONN + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_sift2-propagated_desc-{targets}_meas-raw_connectivity_matrix.csv"
        ),
        assignments=(
            _Q_CONN + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_sift2-propagated_desc-{targets}_method-mrtrix_assignments.txt"
        ),
        qc_metrics=(
            _Q_CONN + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_sift2-propagated_desc-{targets}_connectivity_qc.json"
        ),
        connectome_full=temp(
            _TMP + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_sift2-propagated_desc-{targets}_connectome_full.txt"
        ),
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_sift2-propagated_desc-{targets}_voxelwise_connectivity.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_sift2-propagated_desc-{targets}_voxelwise_connectivity.tsv"
    params:
        slice_script=lambda wc: SLICE_BLOCK_SCRIPT,
        target_labels=lambda wc: ",".join(config["targets"][wc.targets]["labels"]),
        target_search_radius=lambda wc: config.get("mrtrix", {}).get("target_search_radius", 4),
    threads: lambda wc: res("enrich_voxelwise_connectivity_propagated", "threads", 2)
    resources:
        mem_mb=lambda wc: res("enrich_voxelwise_connectivity_propagated", "mem_mb", 8000),
        time=lambda wc: res("enrich_voxelwise_connectivity_propagated", "time_min", 30),
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

        if grep -q '[0-9]' "{input.sift2_weights}"; then
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
            --out-matrix-raw "{output.connectivity_matrix_raw}" \
            --mu-file "{input.mu}" \
            --out-qc "{output.qc_metrics}" \
            &>> "{log}"
        else
          echo "empty seed-filtered tractogram (0 streamlines): writing zero connectome" > "{log}"
          : > "{output.connectome_full}"
          : > "{output.assignments}"
          python "{params.slice_script}" \
            --connectome "{output.connectome_full}" \
            --voxel-index "{input.seed_voxel_index}" \
            --header "{params.target_labels}" \
            --out-matrix "{output.connectivity_matrix}" \
            --out-matrix-raw "{output.connectivity_matrix_raw}" \
            --mu-file "{input.mu}" \
            --out-qc "{output.qc_metrics}" \
            --empty \
            &>> "{log}"
        fi
        """


# -----------------------------
# (6f) Propagated branch + volume-bias correction (only requested when BOTH
#      sift2_propagated AND volume_bias_correction are set). Mirrors
#      enrich_voxelwise_connectivity_invnodevol on the propagated weights.
# -----------------------------
rule enrich_voxelwise_connectivity_propagated_invnodevol:
    input:
        tractogram=rules.enrich_filter_tractogram_propagated.output.tractogram,
        sift2_weights=rules.enrich_filter_tractogram_propagated.output.sift2_weights,
        nodes=rules.enrich_build_selected_nodes.output.nodes,
        seed_voxel_index=rules.enrich_build_selected_nodes.output.seed_voxel_index,
        mu=enrich_mu_max,
    output:
        connectivity_matrix=(
            _T_COND + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_sift2-propagated_scale-invnodevol_desc-{targets}_connectivity_matrix.csv"
        ),
        connectivity_matrix_raw=(
            _Q_CONN + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_sift2-propagated_scale-invnodevol_desc-{targets}_meas-raw_connectivity_matrix.csv"
        ),
        qc_metrics=(
            _Q_CONN + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_sift2-propagated_scale-invnodevol_desc-{targets}_connectivity_qc.json"
        ),
        connectome_full=temp(
            _TMP + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            "_sift2-propagated_scale-invnodevol_desc-{targets}_connectome_full.txt"
        ),
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_sift2-propagated_scale-invnodevol_desc-{targets}_voxelwise_connectivity.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_sift2-propagated_scale-invnodevol_desc-{targets}_voxelwise_connectivity.tsv"
    params:
        slice_script=lambda wc: SLICE_BLOCK_SCRIPT,
        target_labels=lambda wc: ",".join(config["targets"][wc.targets]["labels"]),
        target_search_radius=lambda wc: config.get("mrtrix", {}).get("target_search_radius", 4),
    threads: lambda wc: res("enrich_voxelwise_connectivity_propagated_invnodevol", "threads", 2)
    resources:
        mem_mb=lambda wc: res("enrich_voxelwise_connectivity_propagated_invnodevol", "mem_mb", 8000),
        time=lambda wc: res("enrich_voxelwise_connectivity_propagated_invnodevol", "time_min", 30),
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

        if grep -q '[0-9]' "{input.sift2_weights}"; then
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
            --out-matrix-raw "{output.connectivity_matrix_raw}" \
            --mu-file "{input.mu}" \
            --out-qc "{output.qc_metrics}" \
            &>> "{log}"
        else
          echo "empty seed-filtered tractogram (0 streamlines): writing zero connectome" > "{log}"
          : > "{output.connectome_full}"
          python "{params.slice_script}" \
            --connectome "{output.connectome_full}" \
            --voxel-index "{input.seed_voxel_index}" \
            --header "{params.target_labels}" \
            --out-matrix "{output.connectivity_matrix}" \
            --out-matrix-raw "{output.connectivity_matrix_raw}" \
            --mu-file "{input.mu}" \
            --out-qc "{output.qc_metrics}" \
            --empty \
            &>> "{log}"
        fi
        """


# -----------------------------
# (6g) Tractography visual QC (gated by config tractography_qc).
#      Same artifacts as connectivity.smk (TDI, endpoints, track-weighted DEC,
#      random 200k subset, length CSV+PNG) via the shared scripts/tractogram_qc.py,
#      for the WB baseline tractogram and each condition's seed-filtered tractogram.
#      The FOD-based DEC background (rule fod2dec) is reused per subject from
#      connectivity.smk -- not redefined here. Constants TRACTOGRAM_QC_SCRIPT /
#      QC_NSUBSAMPLE / QC_SUBSAMPLE_SEED come from connectivity.smk (included first).
# -----------------------------
rule enrich_qc_tractogram_wb:
    """Lightweight QC for the whole-brain baseline tractogram (subset only).

    The full QC artifacts (TDI / endpoint / DEC density maps, full-tractogram length
    histogram, random reservoir subset) each stream the entire ~100M-streamline WB
    tractogram -- a single-core, multi-hour straggler that held the node idle at the
    end of the run. The WB-baseline QC here is not load-bearing (the main connectivity
    QC, rule qc_tractogram_wb, keeps the full treatment), so this rule keeps ONLY the
    one trivial artifact: a small leading-N subset for mrview inspection, written with
    `tckedit -number` (early-exit: reads ~N streamlines, NOT the whole tractogram).
    WB tckgen uses `-seed_dynamic`, so leading-N is a spatially uniform whole-brain
    sample. The expensive density maps + full-tractogram length histogram are dropped.
    """
    input:
        tck=rules.wb_tckgen_merge.output.tck,
    output:
        subset=_QC + "/cond-wb/tractogram/sub-{subject}_label-{seed}_desc-wb_subset.tck",
    params:
        n_subsample=QC_NSUBSAMPLE,
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_desc-wb_tractogram_qc.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_desc-wb_tractogram_qc.tsv"
    threads: lambda wc: res("enrich_qc_tractogram_wb", "threads", 2)
    resources:
        mem_mb=lambda wc: res("enrich_qc_tractogram_wb", "mem_mb", 2000),
        time=lambda wc: res("enrich_qc_tractogram_wb", "time_min", 20),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.subset}")"
        mkdir -p "$(dirname "{log}")"

        tckedit -force -nthreads {threads} -number {params.n_subsample} \
          "{input.tck}" "{output.subset}" &> "{log}"
        """


rule enrich_qc_tractogram_filtered:
    """Visual QC for one condition's seed-filtered tractogram + its SIFT2 weights."""
    input:
        tck=rules.enrich_filter_tractogram.output.tractogram,
        weights=rules.enrich_filter_tractogram.output.sift2_weights,
        template=bids(root=root, datatype="dwi", suffix="mask.mif", **subj_wildcards),
    output:
        tdi=_Q_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_desc-seedfiltered_tdi.mif",
        endpoints=_Q_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_desc-seedfiltered_endpoints.mif",
        decmap=_Q_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_desc-seedfiltered_decmap.mif",
        subset=_Q_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_desc-seedfiltered_subset.tck",
        lengths_csv=_Q_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_desc-seedfiltered_lengths.csv",
        lengths_png=_Q_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_desc-seedfiltered_lengths.png",
    params:
        script=lambda wc: TRACTOGRAM_QC_SCRIPT,
        n_subsample=QC_NSUBSAMPLE,
        seed=QC_SUBSAMPLE_SEED,
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_desc-seedfiltered_tractogram_qc.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_desc-seedfiltered_tractogram_qc.tsv"
    threads: lambda wc: res("enrich_qc_tractogram_filtered", "threads", 2)
    resources:
        mem_mb=lambda wc: res("enrich_qc_tractogram_filtered", "mem_mb", 2000),
        time=lambda wc: res("enrich_qc_tractogram_filtered", "time_min", 30),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tdi}")"
        mkdir -p "$(dirname "{log}")"

        python "{params.script}" \
          --tck "{input.tck}" \
          --template "{input.template}" \
          --weights "{input.weights}" \
          --out-tdi "{output.tdi}" \
          --out-endpoints "{output.endpoints}" \
          --out-dec "{output.decmap}" \
          --out-subset "{output.subset}" \
          --out-lengths-csv "{output.lengths_csv}" \
          --out-lengths-png "{output.lengths_png}" \
          --n-subsample {params.n_subsample} \
          --seed {params.seed} \
          --nthreads {threads} \
          &> "{log}"
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
        row=_Q_COND + "/sub-{subject}_label-{seed}_cond-{cond}_metrics.csv",
    params:
        script=lambda wc: METRICS_SCRIPT,
        level_label=lambda wc: "baseline" if wc.cond == "wb" else wc.cond,
        roi_arg=enrich_roi_arg,
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_cond-{cond}_metrics.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_cond-{cond}_metrics.tsv"
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
            _Q_COND + "/sub-{subject}_label-{seed}_cond-{cond}_metrics.csv",
            subject=wc.subject, seed=wc.seed, cond=ENRICH_CONDS,
        ),
    output:
        summary=_QC + "/sub-{subject}_label-{seed}_enrichment_sweep_summary.csv",
    log:
        "logs/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_summary.log",
    benchmark:
        "benchmarks/sub-{subject}/enrichment_sweep/sub-{subject}_label-{seed}_summary.tsv"
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
def _enrich_fingerprint_outputs():
    """Connectivity-fingerprint CSV targets across both SIFT2 branches and the
    optional volume-bias correction. Kept in sync with get_enrichment_sweep_outputs()
    in the Snakefile (that helper runs at parse time, before this module is included,
    so the two cannot share code)."""
    def fp(extra_tag):
        return expand(
            _T_COND + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}"
            + extra_tag + "_desc-{targets}_connectivity_matrix.csv",
            subject=ENRICH_SUBJECTS, seed=ENRICH_SEED, hemi=ENRICH_HEMIS,
            cond=ENRICH_CONDS, targets=ENRICH_TARGETS,
        )

    out = list(fp(""))                              # refit branch, uncorrected
    if ENRICH_VOLNORM:
        out += fp("_scale-invnodevol")              # refit branch, volume-bias-corrected
    if ENRICH_SIFT2_PROP:
        out += fp("_sift2-propagated")              # propagated branch, uncorrected
        if ENRICH_VOLNORM:
            out += fp("_sift2-propagated_scale-invnodevol")  # propagated branch, vbc
    return out


# Tractography visual-QC targets (gated by config tractography_qc): WB baseline +
# every condition's seed-filtered tractogram, plus the per-subject FOD DEC
# background (rule fod2dec, from connectivity.smk).
_ENRICH_QC_ARTS = ["tdi.mif", "endpoints.mif", "decmap.mif", "subset.tck",
                   "lengths.csv", "lengths.png"]
# WB baseline keeps ONLY the trivial subset (see rule enrich_qc_tractogram_wb): the
# full density maps + length histogram over the ~100M WB tractogram were dropped as
# an expensive single-core straggler. Filtered (small) tractograms keep all artifacts.
_ENRICH_QC_WB_ARTS = ["subset.tck"]


def _enrich_qc_outputs():
    if not bool(config.get("tractography_qc", True)):
        return []
    wb = expand(
        _QC + "/cond-wb/tractogram/sub-{subject}_label-{seed}_desc-wb_{art}",
        subject=ENRICH_SUBJECTS, seed=ENRICH_SEED, art=_ENRICH_QC_WB_ARTS,
    )
    filtered = expand(
        _Q_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_cond-{cond}_desc-seedfiltered_{art}",
        subject=ENRICH_SUBJECTS, seed=ENRICH_SEED, hemi=ENRICH_HEMIS,
        cond=ENRICH_CONDS, art=_ENRICH_QC_ARTS,
    )
    fod = expand(
        "sub-{subject}/qc/connectivity/sub-{subject}_desc-fod_decmap.mif",
        subject=ENRICH_SUBJECTS,
    )
    return wb + filtered + fod


rule all_enrichment_sweep:
    input:
        (
            expand(
                _QC + "/sub-{subject}_label-{seed}_enrichment_sweep_summary.csv",
                subject=ENRICH_SUBJECTS, seed=ENRICH_SEED,
            )
            + _enrich_fingerprint_outputs()
            + _enrich_qc_outputs()
        ) if ENRICH_ENABLED else [],
