# mrtrix_wb.smk -- whole-brain, ACT-constrained tractography (iFOD2). Generated
# in parallel chunks with distinct RNG seeds, then merged into one tractogram.


def res(rule, key, default):
    return config.get("resources", {}).get(rule, {}).get(key, default)

WB_CHUNKS = list(range(1, int(config.get("wb_chunk_count", 1)) + 1))

def wb_chunk_seed(wc):
    return int(config.get("mrtrix_rng_seed", 42)) + int(wc.chunk)

def wb_chunk_streamlines():
    return int(config.get("wb_chunk_streamlines", config.get("wb_sl_count", 20000000)))

# -----------------------------
# Whole-brain tractography (ACT)
# -----------------------------

rule wb_tckgen_chunk:
    input:
        wm_fod=get_fod_for_tracking,
        mask=bids(root=root, datatype="dwi", suffix="mask.mif", **subj_wildcards),
        five_tt=rules.act_5ttgen.output.five_tt,
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
        tckinfo="sub-{subject}/qc/tractography/sub-{subject}_desc-wb_chunk-{chunk}_method-mrtrix_tractography_tckinfo.txt",
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
        mkdir -p "$(dirname "{output.tckinfo}")"
        mkdir -p "$(dirname "{log}")"
        export MRTRIX_RNG_SEED={params.mrtrix_rng_seed}

        # Whole-brain tracking is BIDIRECTIONAL (MRtrix default: no
        # -seed_unidirectional flag -> tracking proceeds both ways from each
        # dynamic seed) and uses -backtrack (ACT-aware backtracking).
        tckgen -nthreads {threads} -algorithm iFOD2 \
          -act "{input.five_tt}" \
          -backtrack \
          -crop_at_gmwmi \
          -seed_dynamic "{input.wm_fod}" \
          -mask "{input.mask}" \
          -step {params.step} -angle {params.angle} \
          -minlength {params.minlength} -maxlength {params.maxlength} \
          -cutoff {params.cutoff} \
          -select {params.n_streamlines} \
          "{input.wm_fod}" "{output.tck}" \
          &> "{log}"

        tckinfo "{output.tck}" > "{output.tckinfo}" 2>> "{log}"
        """

rule wb_tckgen_merge:
    input:
        tcks=lambda wc: expand(
            "sub-{subject}/tracts/sub-{subject}_desc-wb_chunk-{chunk}_method-mrtrix_tractography.tck",
            subject=wc.subject,
            chunk=WB_CHUNKS
        ),
        tckinfos=lambda wc: expand(
            "sub-{subject}/qc/tractography/sub-{subject}_desc-wb_chunk-{chunk}_method-mrtrix_tractography_tckinfo.txt",
            subject=wc.subject,
            chunk=WB_CHUNKS
        ),
    output:
        tck=bids(
            root=root, datatype="tracts", method="mrtrix",
            desc="wb", suffix="tractography.tck", **subj_wildcards
        ),
        tckinfo="sub-{subject}/qc/tractography/sub-{subject}_desc-wb_method-mrtrix_tractography_tckinfo.txt",
    log:
        "logs/sub-{subject}/tracts/sub-{subject}_desc-wb_tckgen_merge.log",
    benchmark:
        "benchmarks/sub-{subject}/tracts/sub-{subject}_desc-wb_tckgen_merge.tsv"
    threads: lambda wc: res("wb_tckgen_merge", "threads", 2)
    resources:
        mem_mb=lambda wc: res("wb_tckgen_merge", "mem_mb", 2000),
        time=lambda wc: res("wb_tckgen_merge", "time_min", 45)
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
