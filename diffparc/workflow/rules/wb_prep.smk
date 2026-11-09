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
        # threshold, then dilate by one voxel (radius 1x1x1vox) before resampling/trimming
        "c3d -int 0 {input.seed} -threshold {params.threshold} inf 1 0 "
        "-dilate 1 1x1x1vox -resample-mm {params.resample_res} -trim 0vox "
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


# -----------------------------
# Optional: globus pallidus as WM (config act.gp_as_wm) -- FIRST partial volumes
# -----------------------------
# 5ttgen fsl sums the partial volumes of these ten FIRST structures into one sGM volume,
# so the pallidum's share is rebuilt from FIRST (same call, same T1) rather than read back.
FIRST_HEMIS = ["L", "R"]
FIRST_STRUCTS = ["Accu", "Caud", "Pall", "Puta", "Thal"]
FIRST_PV = "sub-{subject}/anat/first/sub-{subject}_hemi-{hemi}_label-{struct}_desc-first_probseg.nii.gz"


def first_pv(structs):
    return expand(FIRST_PV, subject="{subject}", hemi=FIRST_HEMIS, struct=structs)


rule act_first_pv:
    """
    FIRST partial-volume maps of the ten 5ttgen-fsl sub-cortical structures on the 5TT
    grid, replicating 5ttgen fsl: T1 with strides -1,+2,+3 (regridded to 1 mm for FIRST
    when coarser than 1.225 mm), run_first_all -m none, meshconvert -transform first2real,
    mesh2voxel. The 5TT is a crop of the T1 lattice, so the maps match its sGM voxels.
    """
    input:
        t1=rules.act_5ttgen.input.t1,
        five_tt=rules.act_5ttgen.output.five_tt,
    output:
        pv=first_pv(FIRST_STRUCTS),
    params:
        # same order as output.pv (expand: hemi outer, struct inner)
        structs=" ".join(f"{h}_{s}" for h in FIRST_HEMIS for s in FIRST_STRUCTS),
    log:
        "logs/sub-{subject}/act_first_pv.log",
    threads: 1
    resources:
        mem_mb=lambda wc: res("act_first_pv", "mem_mb", 8000),
        time=lambda wc: res("act_first_pv", "time_min", 90)
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.pv[0]}")" "$(dirname "{log}")"
        structs=({params.structs})
        outs=()
        for o in {output.pv}; do outs+=("$(realpath -m "$o")"); done
        log="$(realpath -m "{log}")"
        tmp="$(mktemp -d "${{TMPDIR:-/tmp}}/first_pv_XXXXXX")"
        trap 'rm -rf "$tmp"' EXIT

        mrconvert -coord 3 0 -axes 0,1,2 "{input.five_tt}" "$tmp/tpl3d.mif" &> "$log"
        mrconvert "{input.t1}" "$tmp/T1.nii" -strides -1,+2,+3 &>> "$log"
        cd "$tmp"
        first_in=T1.nii
        if mrinfo T1.nii -spacing | awk '{{exit !(($1*$2*$3)^(1/3) > 1.225)}}'; then
          mrgrid T1.nii regrid T1_1mm.nii -voxel 1.0 -interp sinc &>> "$log"
          first_in=T1_1mm.nii
        fi
        run_first_all -m none -s "$(IFS=,; echo "${{structs[*]}}")" -i "$first_in" -o first &>> "$log"
        # run_first_all may hand the fits to fsl_sub; wait for every mesh (up to 60 min)
        for i in $(seq 1 240); do
          [ "$(ls first-*_first.vtk 2>/dev/null | wc -l)" -ge "${{#structs[@]}}" ] && break
          sleep 15
        done
        for i in "${{!structs[@]}}"; do
          hs="${{structs[$i]}}"
          meshconvert "first-${{hs}}_first.vtk" "$hs.vtk" -transform first2real "$first_in" &>> "$log"
          mesh2voxel "$hs.vtk" tpl3d.mif "${{outs[$i]}}" &>> "$log"
        done
        """


rule act_5tt_gp_as_wm:
    """
    5TT with the FIRST pallidum moved from sGM to WM. 5ttgen fsl sets
    sGM = min(1 - CSF, sum of the ten FIRST partial volumes); the rule first checks that the
    rebuilt sum reproduces the original sGM, then sets sGM = min(original sGM, sum of the
    eight non-pallidal structures) and adds the removed fraction to WM. cGM, CSF and pathology
    are unchanged, so the tissue sum is preserved.
    """
    input:
        five_tt=rules.act_5ttgen.output.five_tt,
        pv_all=first_pv(FIRST_STRUCTS),
        pv_keep=first_pv([s for s in FIRST_STRUCTS if s != "Pall"]),
    output:
        five_tt="sub-{subject}/anat/sub-{subject}_desc-5ttgpwm_method-mrtrix.mif",
    log:
        "logs/sub-{subject}/act_5tt_gp_as_wm.log",
    threads: 1
    resources:
        mem_mb=lambda wc: res("act_5tt_gp_as_wm", "mem_mb", 4000),
        time=lambda wc: res("act_5tt_gp_as_wm", "time_min", 15)
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        tmp="$(mktemp -d "${{TMPDIR:-/tmp}}/gpwm_XXXXXX")"
        trap 'rm -rf "$tmp"' EXIT

        for i in 0 1 2 3 4; do
          mrconvert -coord 3 $i -axes 0,1,2 "{input.five_tt}" "$tmp/t$i.mif" &>> "{log}"
        done
        # QC: FIRST must reproduce the original sGM, min(1 - CSF, sum of all ten structures)
        mrmath {input.pv_all} sum - | mrcalc - 1 -min "$tmp/all.mif" &>> "{log}"
        mrcalc 1 "$tmp/t3.mif" -sub "$tmp/all.mif" -min "$tmp/t1.mif" -sub -abs "$tmp/d.mif" &>> "{log}"
        maxd="$(mrstats "$tmp/d.mif" -output max)"
        echo "max |rebuilt sGM - 5ttgen sGM| = $maxd" >> "{log}"
        if ! awk -v d="$maxd" 'BEGIN{{exit !(d <= 1e-4)}}'; then
          echo "ERROR: FIRST rerun does not reproduce the 5TT sGM (max diff $maxd)" >> "{log}"; exit 1
        fi
        # sGM without the pallidum; the removed fraction goes to WM
        mrmath {input.pv_keep} sum - | mrcalc - 1 -min "$tmp/keep.mif" &>> "{log}"
        mrcalc "$tmp/t1.mif" "$tmp/keep.mif" -min "$tmp/sgm.mif" &>> "{log}"
        mrcalc "$tmp/t2.mif" "$tmp/t1.mif" "$tmp/sgm.mif" -sub -add "$tmp/wm.mif" &>> "{log}"
        mrcat "$tmp/t0.mif" "$tmp/sgm.mif" "$tmp/wm.mif" "$tmp/t3.mif" "$tmp/t4.mif" - -axis 3 \
          | mrconvert - "{output.five_tt}" -strides +2,+3,+4,+1 -datatype float32 -force &>> "{log}"
        5ttcheck "{output.five_tt}" >> "{log}" 2>&1
        """


# 5TT used for tracking and SIFT2: the original, or the GP-as-WM version (config act.gp_as_wm).
ACT_5TT = (
    rules.act_5tt_gp_as_wm.output.five_tt
    if config.get("act", {}).get("gp_as_wm", False)
    else rules.act_5ttgen.output.five_tt
)


rule act_gmwmi:
    """
    Build GMWMI seed image for ACT seeding.
    """
    input:
        five_tt=ACT_5TT,
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


