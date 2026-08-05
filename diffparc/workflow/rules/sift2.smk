# sift2.smk -- SIFT2 filtering of the final tractogram (ACT-aware), producing
# per-streamline weights and the global mu scaling factor.


def res(rule, key, default):
    return config.get("resources", {}).get(rule, {}).get(key, default)


rule run_sift2:
    input:
        tck=rules.final_tractogram.output.tck,
        wm_fod=get_fod_for_tracking,
        five_tt=rules.act_5ttgen.output.five_tt,
    output:
        weights=bids(
            root=root, datatype="tracts", method="mrtrix",
            desc="final", suffix="sift2_weights.txt", **subj_wildcards
        ),
        mu="sub-{subject}/qc/sub-{subject}_desc-final_method-mrtrix_sift2_mu.txt",
    log:
        "logs/sub-{subject}/tracts/sub-{subject}_desc-sift2_method-mrtrix.log",
    benchmark:
        "benchmarks/sub-{subject}/tracts/sub-{subject}_desc-sift2_method-mrtrix.tsv",
    threads: lambda wc: res("sift2", "threads", 10)
    resources:
        mem_mb=lambda wc: res("sift2", "mem_mb", 40000),
        time=lambda wc: res("sift2", "time_min", 300)
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        mkdir -p "$(dirname "{output.mu}")"

        tcksift2 -nthreads {threads} -act "{input.five_tt}" -out_mu "{output.mu}" \
          "{input.tck}" "{input.wm_fod}" "{output.weights}" \
          &> "{log}"
        """

