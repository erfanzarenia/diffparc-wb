def res(rule, key, default):
    return config.get("resources", {}).get(rule, {}).get(key, default)    

# ROI toggle + seeds list
def roi_enabled():
    return config.get("roi_enrichment", {}).get("enabled", False)

SEEDS = ["vtasncpbp"]

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
        tckinfo="sub-{subject}/qc/sub-{subject}_hemi-{hemi}_label-{seed}_desc-roi_method-mrtrix_tractography_tckinfo.txt",
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

rule roi_tckgen_merge:
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
        tckinfo="sub-{subject}/qc/sub-{subject}_desc-roi_method-mrtrix_tractography_merged_tckinfo.txt",
    log:
        "logs/sub-{subject}/tracts/sub-{subject}_desc-roi_tckgen_merge.log",
    benchmark:
        "benchmarks/sub-{subject}/tracts/sub-{subject}_roi_tckgen_merge.tsv"
    params:
        roi_on=lambda wc: 1 if roi_enabled() else 0
    threads: lambda wc: res("roi_tckgen_merge", "threads", 2)
    resources:
        mem_mb=lambda wc: res("roi_tckgen_merge", "mem_mb", 1000),
        time=lambda wc: res("roi_tckgen_merge", "time_min", 10)
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
            echo "roi_tckgen_merge: ROI enabled but no non-empty input tcks found" >&2
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
rule final_tractogram:
    input:
        wb=rules.wb_tckgen_merge.output.tck,
        roi=rules.roi_tckgen_merge.output.tck
    output:
        tck=temp(bids(
            root=root, datatype="tracts", method="mrtrix",
            desc="final", suffix="tractography.tck", **subj_wildcards
        )),
        tckinfo="sub-{subject}/qc/sub-{subject}_desc-final_method-mrtrix_tractography_tckinfo.txt",
    log:
        "logs/sub-{subject}/tracts/sub-{subject}_desc-final_tractogram.log",
    benchmark:
        "benchmarks/sub-{subject}/tracts/sub-{subject}_desc-final_tractogram.tsv"
    threads: lambda wc: res("tractogram_for_sift2", "threads", 2)
    resources:
        mem_mb=lambda wc: res("final_tractogram", "mem_mb", 2000),
        time=lambda wc: res("final_tractogram", "time_min", 30)
    container: config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tck}")"
        mkdir -p "$(dirname "{output.tckinfo}")"
        mkdir -p "$(dirname "{log}")"

        if [ -s "{input.roi}" ]; then
          tckedit -force "{input.wb}" "{input.roi}" "{output.tck}" &> "{log}"
        else
          ln -sf "$(readlink -f "{input.wb}")" "{output.tck}" &> "{log}"
        fi

        tckinfo "{output.tck}" > "{output.tckinfo}" 2>> "{log}"
        """
