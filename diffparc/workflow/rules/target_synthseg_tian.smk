# target_synthseg_tian.smk -- build a HYBRID target atlas (native SynthSeg cortex +
# template-warped subcortex, e.g. Tian S4) in subject node-grid space, then hand it
# to the SAME build_seed_nodes -> tck2connectome path as any template target.
#
# Nothing here touches tractography, SIFT2 or calibration: the hybrid only changes
# the target-side node image. Cortex is subject-native (SynthSeg --parc), so it needs
# a nearest-neighbour RESLICE onto the node grid, not an inter-subject warp; subcortex
# is warped in with the existing subject<->template transforms (identical mechanism to
# prop_seeds_targets.smk::transform_targets_to_subject).
#
# Hybrid targets are config entries with `kind: hybrid_synthseg_tian` (see snakebids.yml);
# the {targets} wildcard on these rules is constrained to HYBRID_TARGETS_RE (defined in
# the Snakefile) so it never collides with the template-target warp rule.

import os

HYBRID_TARGET_SCRIPT = os.path.join(workflow.basedir, "scripts", "build_hybrid_target.py")


rule reslice_synthseg_cortex_to_nodegrid:
    """Reslice the SynthSeg --parc cortical parcellation (subject-native, T1w space)
    onto the upsampled node grid with nearest-neighbour -- pure grid change, no
    registration (cortex is already in the subject's own space)."""
    input:
        cortex=bids(
            root=root,
            datatype="anat",
            desc="synthsegcortparc",
            suffix="dseg.nii.gz",
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
        cortex=temp(
            bids(
                root=root,
                datatype="anat",
                res="nodegrid",
                desc="synthsegcortparc",
                suffix="dseg.nii.gz",
                **subj_wildcards,
            )
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        "c3d -int 0 {input.ref} {input.cortex} -reslice-identity -o {output.cortex}"


rule warp_hybrid_subcortex_to_nodegrid:
    """Warp the hybrid's subcortical template atlas (e.g. Tian S4) into subject
    node-grid space (nearest-neighbour), reusing the existing subject<->template
    transforms. Same warp as transform_targets_to_subject, on the subcortex-only dseg."""
    wildcard_constraints:
        targets=HYBRID_TARGETS_RE,
    input:
        subcortex=lambda wildcards: os.path.join(
            workflow.basedir, "..", config["targets"][wildcards.targets]["subcortex_dseg"]
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
        inv_warp=bids(
            root=root,
            datatype="warps",
            suffix="invwarp.nii.gz",
            from_="subject",
            to=config["template"],
            **subj_wildcards,
        ),
        affine_xfm_itk=bids(
            root=root,
            datatype="warps",
            suffix="affine.txt",
            from_="subject",
            to=config["template"],
            desc="itk",
            **subj_wildcards,
        ),
    output:
        subcortex=temp(
            bids(
                root=root,
                datatype="anat",
                res="nodegrid",
                desc="{targets}",
                suffix="subcortexdseg.nii.gz",
                **subj_wildcards,
            )
        ),
    envmodules:
        "ants",
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    threads: 8
    shell:
        "antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.subcortex} "
        "-o {output.subcortex} -r {input.ref} -t [{input.affine_xfm_itk},1] {input.inv_warp}"


rule build_hybrid_target:
    """Merge native cortex + warped subcortex into one contiguous 1..N target dseg
    (scripts/build_hybrid_target.py). The output desc-{targets}_dseg is consumed
    unchanged by build_seed_nodes -> tck2connectome, exactly like a template target."""
    wildcard_constraints:
        targets=HYBRID_TARGETS_RE,
    input:
        cortex=rules.reslice_synthseg_cortex_to_nodegrid.output.cortex,
        subcortex=rules.warp_hybrid_subcortex_to_nodegrid.output.subcortex,
        lut=lambda wildcards: os.path.join(
            workflow.basedir, "..", config["targets"][wildcards.targets]["cortex_lut"]
        ),
    params:
        script=HYBRID_TARGET_SCRIPT,
        n_cortex=lambda wildcards: config["targets"][wildcards.targets]["n_cortex"],
    output:
        dseg=bids(
            root=root,
            datatype="anat",
            desc="{targets}",
            suffix="dseg.nii.gz",
            **subj_wildcards,
        ),
        qc=(
            "sub-{subject}/qc/connectivity/"
            "sub-{subject}_desc-{targets}_hybrid_target_qc.json"
        ),
    log:
        bids(
            root="logs",
            desc="{targets}",
            suffix="buildhybridtarget.log",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.dseg}")" "$(dirname "{output.qc}")" "$(dirname "{log}")"

        python "{params.script}" \
          --cortex-dseg "{input.cortex}" \
          --cortex-lut "{input.lut}" \
          --subcortex-dseg "{input.subcortex}" \
          --n-cortex {params.n_cortex} \
          --out-dseg "{output.dseg}" \
          --out-qc "{output.qc}" \
          &> "{log}"
        """
