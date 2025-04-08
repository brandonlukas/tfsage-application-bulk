rule filter_motifs:
    input:
        motif2factors="resources/motif_databases/{motif_db}.motif2factors.txt",
        pfm="resources/motif_databases/{motif_db}.pfm",
    output:
        motif2factors="resources/motif_databases_filtered/{motif_db}_filtered.motif2factors.txt",
        pfm="resources/motif_databases_filtered/{motif_db}_filtered.pfm",
    params:
        factor_list=config["motif_scan"]["factors"],
    script:
        "../scripts/motif_scan/filter_motifs.py"


rule motif_scan:
    input:
        config["downloads_dir"] + "{threshold}/{experiment}.bed",
        pfm="resources/motif_databases_filtered/{motif_db}_filtered.pfm",
    output:
        config["results_dir"] + "motif_scan/{threshold}/{motif_db}/{experiment}.bed",
    params:
        opts=lambda w, input: f"-g GRCh38 -p {input.pfm}",
    threads: 24
    resources:
        mem_mb=24000,
    conda:
        "gimme_env"
    shell:
        """
        gimme scan {input[0]} {params.opts} -N {threads} -b > {output}
        """
