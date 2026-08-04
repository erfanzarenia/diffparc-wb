# mask_sweep.smk
# Experimental mask sweep.
# Compares multiple seed sources:
#   1. native: existing template-to-subject seed probseg
#   2. brainsteminjected: template SN/VTA prior warped along a template-brainstem to subject-brainstem warp
# Reuses existing final tractogram + SIFT2 weights.
# Does NOT overwrite main pipeline masks.

import os

SUBJECTS = sorted(set(inputs.input_zip_lists["T1w"]["subject"]))

BUILD_NODES_SCRIPT = os.path.join(workflow.basedir, "scripts", "build_voxel_nodes.py")
SLICE_BLOCK_SCRIPT = os.path.join(workflow.basedir, "scripts", "slice_connectome_block.py")

SWEEP_SEEDS = ["vtasnc", "snc", "vtapbp", "vtasncpbp"]
HEMIS = ["L", "R"]
TARGETS = ["Yeo7TianS4"]

# QC directory layout (ORGANIZATIONAL). The sweep has many combos (variant x
# threshold x seed x mask_source x hemi), so QC is grouped by threshold and then by
# type -- mirroring rules/enrichment_sweep.smk's cond-{cond}/{connectome,tractogram}/
# hierarchy -- instead of cluttering one flat directory:
#   qc/mask_sweep/{seed}/{mask_source}/thr-{thr_tag}/connectome/   connectivity qc json,
#                                                                  assignments (per variant)
#   qc/mask_sweep/{seed}/{mask_source}/thr-{thr_tag}/tractogram/   seed-filtered tckinfo +
#                                                                  visual QC (tdi/endpoints/
#                                                                  decmap/subset/lengths)
_MS_QC_THR = "sub-{subject}/qc/mask_sweep/{seed}/{mask_source}/thr-{thr_tag}"
_MS_QC_CONN = _MS_QC_THR + "/connectome"
_MS_QC_TRACT = _MS_QC_THR + "/tractogram"

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


# snc and vtapbp touch at the winner-take-all boundary in template space, then
# are warped to subject space INDEPENDENTLY with linear interpolation, which
# smears each probseg ~1 voxel past its edge and makes the thresholded masks
# overlap. Re-imposing the argmax in subject space (keep snc where it out-probs
# vtapbp, and vice-versa) restores disjointness. Maps each split seed to its
# counterpart; seeds absent here (vtasnc, vtasncpbp) get no argmax.
SWEEP_ARGMAX_PARTNER = {"snc": "vtapbp", "vtapbp": "snc"}


def _seed_probseg_path(subject, hemi, seed, mask_source):
    if mask_source == "native":
        return (
            f"sub-{subject}/anat/"
            f"sub-{subject}_hemi-{hemi}_label-{seed}_probseg.nii.gz"
        )
    if mask_source == "brainsteminjected":
        return (
            f"sub-{subject}/anat/mask_sweep/{seed}/brainsteminjected/"
            f"sub-{subject}_hemi-{hemi}_label-{seed}"
            f"_desc-brainsteminjected_probseg.nii.gz"
        )
    raise ValueError(f"Unknown mask_source: {mask_source}")


def sweep_seed_probseg(wc):
    return _seed_probseg_path(wc.subject, wc.hemi, wc.seed, wc.mask_source)


def sweep_partner_probseg(wc):
    # Counterpart probseg for the subject-space argmax. For seeds without a
    # partner, return the seed's own probseg (placeholder; unused in the shell).
    partner = SWEEP_ARGMAX_PARTNER.get(wc.seed, wc.seed)
    return _seed_probseg_path(wc.subject, wc.hemi, partner, wc.mask_source)


def _mask_sweep_all_outputs():
    """Convenience-target output list for all_mask_sweep. Kept consistent with the
    Snakefile's get_mask_sweep_outputs(): baseline matrices + config-gated corrected
    variants + (when tractography_qc) the sweep visual QC."""
    if not MASK_SWEEP_ENABLED:
        return []

    common = dict(
        subject=SUBJECTS,
        mask_source=MASK_SOURCES,
        hemi=HEMIS,
        seed=SWEEP_SEEDS,
        targets=TARGETS,
        thr_tag=threshold_tag_list(),
    )

    out = expand(
        "sub-{subject}/connectivity/mask_sweep/{seed}/{mask_source}/"
        "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_connectivity_matrix.csv",
        **common,
    )

    meas = conn_variant_meas_enabled()
    if meas:
        out += expand(
            "sub-{subject}/connectivity/mask_sweep/{seed}/{mask_source}/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}"
            "_meas-{meas}_connectivity_matrix.csv",
            meas=meas,
            **common,
        )

    if config.get("tractography_qc", True):
        arts = ["tdi.mif", "endpoints.mif", "decmap.mif", "subset.tck",
                "lengths.csv", "lengths.png"]
        # Visual QC is variant-independent (one tractogram per threshold): no {targets}.
        out += expand(
            _MS_QC_TRACT + "/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-seedfiltered_{art}",
            subject=SUBJECTS,
            mask_source=MASK_SOURCES,
            hemi=HEMIS,
            seed=SWEEP_SEEDS,
            thr_tag=threshold_tag_list(),
            art=arts,
        )

    return out


rule all_mask_sweep:
    input:
        _mask_sweep_all_outputs()

        
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
        partner_prob=sweep_partner_probseg,
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
        argmax=lambda wc: 1 if wc.seed in SWEEP_ARGMAX_PARTNER else 0,
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

        if [ "{params.argmax}" -eq 1 ]; then
          # Subject-space winner-take-all. The split seeds touch in template space
          # and are warped here independently with linear interpolation, which
          # smears each probseg past its edge so the thresholded masks overlap.
          # Fix: keep a voxel only where this seed's prob is strictly greater than
          # its partner's. Strict-greater for both seeds makes overlap impossible.
          #
          # Step 1: reslice partner onto the seed grid (c3d order: REF MOVING),
          #         then win = (seed - partner) > 0.
          win="$(mktemp --suffix=.nii.gz)"
          trap 'rm -f "$win"' EXIT
          c3d -int 0 "{input.seed_prob}" "{input.partner_prob}" -reslice-identity \
            "{input.seed_prob}" -scale -1 -add \
            -threshold 0 inf 0 1 \
            -o "$win" \
            &> "{log}"

          # Step 2: thresholded seed AND win, then resample/trim as usual.
          c3d -int 0 "{input.seed_prob}" -threshold {params.threshold} inf 1 0 \
            "$win" -multiply \
            -resample-mm {params.resample_res} \
            -trim 0vox \
            -type uchar \
            -o "{output.seed_mask}" \
            &>> "{log}"
        else
          c3d -int 0 "{input.seed_prob}" \
            -threshold {params.threshold} inf 1 0 \
            -resample-mm {params.resample_res} \
            -trim 0vox \
            -type uchar \
            -o "{output.seed_mask}" \
            &> "{log}"
        fi
        """


# -----------------------------
# Paste thresholded seed onto the full-FOV node grid (== target-atlas grid)
# -----------------------------

rule sweep_reslice_seed_to_nodegrid:
    """
    Reslice the thresholded sweep seed onto the full-FOV upsampled reference grid
    (== target-atlas grid) with nearest-neighbour, used for BOTH filtering and
    node building so they are voxel-consistent. For the native mask source this
    is a lossless paste (same lattice); for brainsteminjected it is a genuine NN
    regrid, inherent to that experimental source.
    """
    input:
        seed_mask=rules.sweep_binarize_seed.output.seed_mask,
        ref=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res="upsampled",
            datatype="dwi",
            **subj_wildcards,
        ),
    output:
        mask=temp(
            config["tmp_dir"]
            + "/sub-{subject}/tracts/mask_sweep/{seed}/{mask_source}/"
            + "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_res-nodegrid_mask.nii.gz"
        ),
    log:
        "logs/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_reslice_seed_to_nodegrid.log",
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

        c3d -int 0 "{input.ref}" "{input.seed_mask}" \
          -reslice-identity \
          -threshold 0.5 inf 1 0 \
          -type uchar \
          -o "{output.mask}" \
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
        seed_mask=rules.sweep_reslice_seed_to_nodegrid.output.mask,
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
            _MS_QC_TRACT + "/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-seedfiltered_method-mrtrix_tractography_tckinfo.txt"
        ),
    log:
        "logs/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_filter_tractogram.log",
    benchmark:
        "benchmarks/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_filter_tractogram.tsv",
    threads: lambda wc: config.get("mrtrix", {}).get("filter_threads", 2)
    resources:
        mem_mb=8000,
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
# Per-voxel node parcellation (seed voxels + offset target atlas)
# -----------------------------

rule sweep_build_seed_nodes:
    input:
        seed_mask=rules.sweep_reslice_seed_to_nodegrid.output.mask,
        targets_dseg=ancient("sub-{subject}/anat/sub-{subject}_desc-{targets}_dseg.nii.gz"),
    output:
        nodes=temp(
            config["tmp_dir"]
            + "/sub-{subject}/tracts/mask_sweep/{seed}/{mask_source}/"
            + "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_nodes.nii.gz"
        ),
        seed_voxel_index=(
            "sub-{subject}/qc/mask_sweep/{seed}/{mask_source}/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_seed_voxel_index.csv"
        ),
    log:
        "logs/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_build_seed_nodes.log",
    benchmark:
        "benchmarks/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_build_seed_nodes.tsv",
    params:
        script=lambda wc: BUILD_NODES_SCRIPT,
        n_targets=lambda wc: len(config["targets"][wc.targets]["labels"]),
    threads: 1
    resources:
        mem_mb=4000,
        time=10
    container:
        config["singularity"]["diffparc"]
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
# Voxelwise seed-to-target connectivity (MRtrix tck2connectome)
# -----------------------------

rule sweep_voxelwise_connectivity:
    input:
        tractogram=rules.sweep_filter_tractogram.output.tractogram,
        sift2_weights=rules.sweep_filter_tractogram.output.sift2_weights,
        nodes=rules.sweep_build_seed_nodes.output.nodes,
        seed_voxel_index=rules.sweep_build_seed_nodes.output.seed_voxel_index,
        # Main-pipeline SIFT2 mu: the sweep filters the SAME final tractogram +
        # SIFT2 weights, so the same global mu calibrates these sub-blocks.
        # ancient() so a touched main SIFT2 fit doesn't force the sweep to re-run,
        # matching how the tractogram/weights are referenced above.
        mu=ancient(
            "sub-{subject}/tracts/sub-{subject}_desc-final_method-mrtrix_sift2_mu.txt"
        ),
    output:
        # MAIN deliverable = the RAW SIFT2 weight-sum block (untagged csv); the
        # SIFT2-mu-scaled block is preserved alongside as meas-muscaled. Same scheme
        # as connectivity.smk (mirrored naming).
        connectivity_matrix=(
            "sub-{subject}/connectivity/mask_sweep/{seed}/{mask_source}/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_connectivity_matrix.csv"
        ),
        connectivity_matrix_muscaled=(
            "sub-{subject}/connectivity/mask_sweep/{seed}/{mask_source}/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_meas-muscaled_connectivity_matrix.csv"
        ),
        assignments=(
            _MS_QC_CONN + "/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_method-mrtrix_assignments.txt"
        ),
        qc_metrics=(
            _MS_QC_CONN + "/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_connectivity_qc.json"
        ),
        connectome_full=temp(
            config["tmp_dir"]
            + "/sub-{subject}/qc/mask_sweep/{seed}/{mask_source}/thr-{thr_tag}/connectome/"
            + "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_connectome_full.txt"
        ),
    log:
        "logs/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_voxelwise_connectivity.log",
    benchmark:
        "benchmarks/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_voxelwise_connectivity.tsv",
    params:
        slice_script=lambda wc: SLICE_BLOCK_SCRIPT,
        target_labels=lambda wc: ",".join(config["targets"][wc.targets]["labels"]),
        target_search_radius=lambda wc: config.get("mrtrix", {}).get("target_search_radius", 4),
    threads: lambda wc: config.get("mrtrix", {}).get("voxelagg_threads", 1)
    resources:
        mem_mb=2000,
        time=90
    container:
        config["singularity"]["diffparc"]
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

        # --out-matrix is the script's "primary" (mu-scaled) -> meas-muscaled;
        # --out-matrix-raw is the RAW block -> the untagged main deliverable.
        python "{params.slice_script}" \
          --connectome "{output.connectome_full}" \
          --voxel-index "{input.seed_voxel_index}" \
          --header "{params.target_labels}" \
          --out-matrix "{output.connectivity_matrix_muscaled}" \
          --out-matrix-raw "{output.connectivity_matrix}" \
          --mu-file "{input.mu}" \
          --out-qc "{output.qc_metrics}" \
          &>> "{log}"
        """


# -----------------------------
# Configurable connectivity-matrix variants (config: connectivity_variants)
# Mirrors rules/connectivity.smk:connectivity_variant on the sweep seed-filtered
# tractogram, using the SAME meas-tag naming (thr-{thr_tag} is the only extra
# sweep entity). Additive to the baseline raw + mu-scaled matrices above.
# -----------------------------

rule sweep_connectivity_variant:
    wildcard_constraints:
        meas=CONN_VARIANT_MEAS_CONSTRAINT,
    input:
        tractogram=rules.sweep_filter_tractogram.output.tractogram,
        sift2_weights=rules.sweep_filter_tractogram.output.sift2_weights,
        nodes=rules.sweep_build_seed_nodes.output.nodes,
        seed_voxel_index=rules.sweep_build_seed_nodes.output.seed_voxel_index,
        # Same global main-pipeline SIFT2 mu the baseline sweep matrix uses.
        mu=ancient(
            "sub-{subject}/tracts/sub-{subject}_desc-final_method-mrtrix_sift2_mu.txt"
        ),
    output:
        connectivity_matrix=(
            "sub-{subject}/connectivity/mask_sweep/{seed}/{mask_source}/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}"
            "_meas-{meas}_connectivity_matrix.csv"
        ),
        qc_metrics=(
            _MS_QC_CONN + "/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}"
            "_meas-{meas}_connectivity_qc.json"
        ),
        connectome_full=temp(
            config["tmp_dir"]
            + "/sub-{subject}/qc/mask_sweep/{seed}/{mask_source}/thr-{thr_tag}/connectome/"
            + "sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}"
            + "_meas-{meas}_connectome_full.txt"
        ),
    log:
        "logs/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_meas-{meas}_connectivity_variant.log",
    benchmark:
        "benchmarks/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-{targets}_meas-{meas}_connectivity_variant.tsv",
    params:
        slice_script=lambda wc: SLICE_BLOCK_SCRIPT,
        target_labels=lambda wc: ",".join(config["targets"][wc.targets]["labels"]),
        target_search_radius=lambda wc: config.get("mrtrix", {}).get("target_search_radius", 4),
        scale_flags=lambda wc: CONN_VARIANT_RECIPES[wc.meas]["scale"],
        use_mu=lambda wc: CONN_VARIANT_RECIPES[wc.meas]["mu"],
    threads: lambda wc: config.get("mrtrix", {}).get("voxelagg_threads", 1)
    resources:
        mem_mb=2000,
        time=90
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.connectivity_matrix}")"
        mkdir -p "$(dirname "{output.qc_metrics}")"
        mkdir -p "$(dirname "{output.connectome_full}")"
        mkdir -p "$(dirname "{log}")"

        tck2connectome -nthreads {threads} -force -symmetric \
          {params.scale_flags} \
          -assignment_radial_search {params.target_search_radius} \
          -tck_weights_in "{input.sift2_weights}" \
          "{input.tractogram}" "{input.nodes}" "{output.connectome_full}" \
          &> "{log}"

        MU_ARG=""
        if [ "{params.use_mu}" -eq 1 ]; then MU_ARG="--mu-file {input.mu}"; fi

        python "{params.slice_script}" \
          --connectome "{output.connectome_full}" \
          --voxel-index "{input.seed_voxel_index}" \
          --header "{params.target_labels}" \
          --out-matrix "{output.connectivity_matrix}" \
          $MU_ARG \
          --out-qc "{output.qc_metrics}" \
          &>> "{log}"
        """


# -----------------------------
# Tractography visual QC for the sweep seed-filtered tractogram.
# Mirrors rules/connectivity.smk:qc_tractogram_filtered (same scripts/tractogram_qc.py
# artifacts: TDI, endpoint map, track-weighted DEC, random subset, length CSV+PNG),
# organized under thr-{thr_tag}/tractogram/. Reuses the shared TRACTOGRAM_QC_SCRIPT /
# QC_NSUBSAMPLE / QC_SUBSAMPLE_SEED from connectivity.smk. Gated by tractography_qc
# (via get_mask_sweep_outputs / all_mask_sweep).
# -----------------------------

rule sweep_qc_tractogram_filtered:
    input:
        tck=rules.sweep_filter_tractogram.output.tractogram,
        weights=rules.sweep_filter_tractogram.output.sift2_weights,
        template=bids(root=root, datatype="dwi", suffix="mask.mif", **subj_wildcards),
    output:
        tdi=_MS_QC_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-seedfiltered_tdi.mif",
        endpoints=_MS_QC_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-seedfiltered_endpoints.mif",
        decmap=_MS_QC_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-seedfiltered_decmap.mif",
        subset=_MS_QC_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-seedfiltered_subset.tck",
        lengths_csv=_MS_QC_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-seedfiltered_lengths.csv",
        lengths_png=_MS_QC_TRACT + "/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-seedfiltered_lengths.png",
    params:
        script=lambda wc: TRACTOGRAM_QC_SCRIPT,
        n_subsample=QC_NSUBSAMPLE,
        seed=QC_SUBSAMPLE_SEED,
    log:
        "logs/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-seedfiltered_tractogram_qc.log",
    benchmark:
        "benchmarks/sub-{subject}/mask_sweep/{seed}/{mask_source}/sub-{subject}_hemi-{hemi}_label-{seed}_thr-{thr_tag}_desc-seedfiltered_tractogram_qc.tsv",
    threads: 2
    resources:
        mem_mb=2000,
        time=30
    container:
        config["singularity"]["diffparc"]
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
