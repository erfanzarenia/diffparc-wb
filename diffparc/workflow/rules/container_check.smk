rule check_diffparc_container:
    output:
        "logs/container_check_fsl.ok"
    singularity:
        config["singularity"]["diffparc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p logs

        echo "=== PATH ==="
        echo "$PATH"

        echo "=== FSL ==="
        command -v fast flirt bet fslmaths fslstats run_first_all
        fslversion

        echo "=== MRTRIX ==="
        command -v 5ttgen tckgen tcksift2 mrconvert

        echo "=== PYTHON ==="
        command -v python3
        python3 --version
        python3 - <<'PY'
try:
    import imp
    print("imp: OK")
except Exception as e:
    print("imp: FAIL ->", e)
PY

        touch {output}
        """

