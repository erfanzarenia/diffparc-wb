import os

def res(rule, key, default):
    return config.get("resources", {}).get(rule, {}).get(key, default)


BUILD_NODES_SCRIPT = os.path.join(
    workflow.basedir, "scripts", "build_voxel_nodes.py",
)
SLICE_BLOCK_SCRIPT = os.path.join(
    workflow.basedir, "scripts", "slice_connectome_block.py",
)


# -----------------------------
# Paste seed mask onto the full-FOV node grid (== target-atlas grid)
# -----------------------------
rule reslice_seed_to_nodegrid:
    """
    Reslice the (trimmed) seed mask back onto the full-FOV upsampled reference
    grid, so it shares an *identical* voxel grid with the target atlas.

    The seed and the target atlas are both warped to this same reference grid
    (see prop_seeds_targets.smk), and the canonical seed mask is only a *cropped*
    (``-trim``) sub-region of that grid. Reslicing with nearest-neighbour
    (``-int 0``) onto the reference therefore just pastes the seed back into the
    full field of view -- no interpolation of seed values, no grid shift. The
    seed mask voxels are preserved exactly.

    This single resliced mask is used for BOTH streamline filtering and node
    building, guaranteeing they are voxel-for-voxel consistent.
    """
    input:
        seed_mask=bids(
            root=root,
            datatype="anat",
            hemi="{hemi}",
            label="{seed}",
            suffix="mask.nii.gz",
            **subj_wildcards,
        ),
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
        mask=bids(
            root=root,
            datatype="anat",
            hemi="{hemi}",
            label="{seed}",
            res="nodegrid",
            suffix="mask.nii.gz",
            **subj_wildcards,
        ),
    log:
        (
            "logs/sub-{subject}/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_reslice_seed_to_nodegrid.log"
        ),
    threads: 1
    resources:
        mem_mb=lambda wc: res("reslice_seed_to_nodegrid", "mem_mb", 4000),
        time=lambda wc: res("reslice_seed_to_nodegrid", "time_min", 10)
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.mask}")"
        mkdir -p "$(dirname "{log}")"

        # REF MOVING -reslice-identity: paste seed (MOVING) into ref grid (REF),
        # nearest-neighbour, no interpolation. Threshold is a binary safety net.
        c3d -int 0 "{input.ref}" "{input.seed_mask}" \
          -reslice-identity \
          -threshold 0.5 inf 1 0 \
          -type uchar \
          -o "{output.mask}" \
          &> "{log}"
        """


# -----------------------------
# Filter final tractogram by seed endpoints
# -----------------------------
rule filter_tractogram:
    """
    Extract seed-connected streamlines from the final tractogram.

    Keeps streamlines where at least one endpoint falls inside the (node-grid)
    seed mask, and writes the matching subset of SIFT2 weights. Using the same
    node-grid mask here as in node building guarantees that every retained seed
    endpoint lands in a labelled seed voxel, so tck2connectome assigns it to its
    own voxel (effective end-voxel assignment on the seed side).
    """
    input:
        tractogram=rules.final_tractogram.output.tck,
        sift2_weights=rules.run_sift2.output.weights,
        seed_mask=rules.reslice_seed_to_nodegrid.output.mask,
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
# Per-voxel node parcellation (seed voxels + offset target atlas)
# -----------------------------
rule build_seed_nodes:
    """
    Build the combined node image consumed by tck2connectome.

    The seed mask is already on the target-atlas grid (reslice_seed_to_nodegrid),
    so we just label every seed voxel as its own node and shift the atlas labels
    above the seed range. See build_voxel_nodes.py.
    """
    input:
        seed_mask=rules.reslice_seed_to_nodegrid.output.mask,
        targets_dseg="sub-{subject}/anat/sub-{subject}_desc-{targets}_dseg.nii.gz",
    output:
        nodes=temp(
            config["tmp_dir"]
            + "/sub-{subject}/tracts/"
            + "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}_nodes.nii.gz"
        ),
        seed_voxel_index=(
            "sub-{subject}/qc/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_seed_voxel_index.csv"
        ),
    log:
        (
            "logs/sub-{subject}/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_build_seed_nodes.log"
        ),
    benchmark:
        (
            "benchmarks/sub-{subject}/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_build_seed_nodes.tsv"
        ),
    params:
        script=lambda wc: BUILD_NODES_SCRIPT,
        n_targets=lambda wc: len(config["targets"][wc.targets]["labels"]),
    threads: lambda wc: res("build_seed_nodes", "threads", 1)
    resources:
        mem_mb=lambda wc: res("build_seed_nodes", "mem_mb", 4000),
        time=lambda wc: res("build_seed_nodes", "time_min", 10)
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
# Voxelwise seed-to-target connectivity (MRtrix tck2connectome)
# -----------------------------
rule voxelwise_connectivity:
    """
    Voxelwise seed-to-target connectivity from the seed-filtered tractogram and
    SIFT2 weights, computed natively with tck2connectome on the per-voxel node
    image, then sliced down to the [seed voxel x atlas target] block.

    Asymmetric endpoint assignment (native): the seed side is effectively
    end-voxel, because tck2connectome's radial search tests the endpoint's own
    voxel first and every seed endpoint sits inside a labelled seed voxel; the
    target side uses the full ``target_search_radius`` because distal endpoints
    stop in unlabelled white matter and must search outward to reach a parcel.
    """
    input:
        tractogram=rules.filter_tractogram.output.tractogram,
        sift2_weights=rules.filter_tractogram.output.sift2_weights,
        nodes=rules.build_seed_nodes.output.nodes,
        seed_voxel_index=rules.build_seed_nodes.output.seed_voxel_index,
    output:
        connectivity_matrix=(
            "sub-{subject}/tracts/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_connectivity_matrix.csv"
        ),
        assignments=(
            "sub-{subject}/qc/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_method-mrtrix_assignments.txt"
        ),
        qc_metrics=(
            "sub-{subject}/qc/connectivity/"
            "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            "_connectivity_qc.json"
        ),
        connectome_full=temp(
            config["tmp_dir"]
            + "/sub-{subject}/qc/connectivity/"
            + "sub-{subject}_hemi-{hemi}_label-{seed}_desc-{targets}"
            + "_connectome_full.txt"
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
        slice_script=lambda wc: SLICE_BLOCK_SCRIPT,
        target_labels=lambda wc: ",".join(config["targets"][wc.targets]["labels"]),
        target_search_radius=lambda wc: config.get("mrtrix", {}).get(
            "target_search_radius", 4
        ),
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
