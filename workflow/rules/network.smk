rule find_targets_tfsage:
    input:
        predictions=config["results_dir"]
        + "generate/{threshold}/{query_id}_{factor}.parquet",
        query_file=config["downloads_dir"] + "{threshold}/{query_id}.bed",
    output:
        config["results_dir"]
        + "network/{threshold}/{linkage}/target_genes/tfsage/{query_id}_{factor}.json",
    params:
        method_class="tfsage",
        method_name=f"TFSage-{config['tfsage']['n']}",
    script:
        "../scripts/network/find_targets.py"


rule find_targets_motif_scan:
    input:
        predictions=config["results_dir"]
        + "motif_scan/{threshold}/{motif_db}/{query_id}.bed",
        query_file=config["downloads_dir"] + "{threshold}/{query_id}.bed",
        m2f_path="resources/motif_databases_filtered/{motif_db}_filtered.motif2factors.txt",
    output:
        config["results_dir"]
        + "network/{threshold}/{linkage}/target_genes/motif_scan/{motif_db}/{query_id}_{factor}.json",
    params:
        method_class="motif scan",
        method_name="{motif_db}",
    script:
        "../scripts/network/find_targets.py"


rule enrichment_analysis:
    input:
        tfsage=expand(
            config["results_dir"]
            + "network/{threshold}/{linkage}/target_genes/tfsage/{query_id}_{factor}.json",
            threshold=config["thresholds"],
            linkage=config["linkages"],
            query_id=config["query_ids"],
            factor=factor_list,
        ),
        motif_scan=expand(
            config["results_dir"]
            + "network/{threshold}/{linkage}/target_genes/motif_scan/{motif_db}/{query_id}_{factor}.json",
            threshold=config["thresholds"],
            linkage=config["linkages"],
            motif_db=config["motif_scan"]["motif_dbs"],
            query_id=config["query_ids"],
            factor=config["motif_scan"]["factors"],
        ),
    output:
        config["results_dir"] + "network/enrichment_analysis.csv",
    params:
        deseq2_results=config["deseq2_results"],
        padj_threshold=0.05,
        log2fc_threshold=1,
    script:
        "../scripts/network/enrichment_analysis.py"
