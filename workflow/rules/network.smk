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


rule compile_targets_tfsage:
    input:
        expand(
            config["results_dir"]
            + "network/{threshold}/{linkage}/target_genes/tfsage/{query_id}_{factor}.json",
            factor=factor_list,
            allow_missing=True,
        ),
    output:
        config["results_dir"]
        + "network/{threshold}/{linkage}/networks/tfsage/{query_id}.parquet",
    script:
        "../scripts/network/compile_targets.py"


rule compile_targets_motif_scan:
    input:
        expand(
            config["results_dir"]
            + "network/{threshold}/{linkage}/target_genes/motif_scan/{motif_db}/{query_id}_{factor}.json",
            factor=config["motif_scan"]["factors"],
            allow_missing=True,
        ),
    output:
        config["results_dir"]
        + "network/{threshold}/{linkage}/networks/motif_scan/{motif_db}/{query_id}.parquet",
    script:
        "../scripts/network/compile_targets.py"


rule compute_enrichment_tfsage:
    input:
        config["results_dir"]
        + "network/{threshold}/{linkage}/networks/tfsage/{query_id}.parquet",
    output:
        config["results_dir"]
        + "network/{threshold}/{linkage}/enrichment/tfsage/{query_id}.json",
    params:
        method_class="tfsage",
        method_name=f"TFSage-{config['tfsage']['n']}",
        de_res=config["enrichment_analysis"]["de_res"],
        padj_threshold=config["enrichment_analysis"]["padj_threshold"],
        log2fc_threshold=config["enrichment_analysis"]["log2fc_threshold"],
    script:
        "../scripts/network/compute_enrichment.py"


rule compute_enrichment_motif_scan:
    input:
        config["results_dir"]
        + "network/{threshold}/{linkage}/networks/motif_scan/{motif_db}/{query_id}.parquet",
    output:
        config["results_dir"]
        + "network/{threshold}/{linkage}/enrichment/motif_scan/{motif_db}/{query_id}.json",
    params:
        method_class="motif scan",
        method_name="{motif_db}",
        de_res=config["enrichment_analysis"]["de_res"],
        padj_threshold=config["enrichment_analysis"]["padj_threshold"],
        log2fc_threshold=config["enrichment_analysis"]["log2fc_threshold"],
    script:
        "../scripts/network/compute_enrichment.py"


rule compute_enrichment_decoupler:
    input:
        "resources/decoupler/{dc_network}.csv",
    output:
        config["results_dir"] + "network/curated/enrichment/decoupler/{dc_network}.json",
    params:
        method_class="decoupler",
        method_name="{dc_network}",
        de_res=config["enrichment_analysis"]["de_res"],
        padj_threshold=config["enrichment_analysis"]["padj_threshold"],
        log2fc_threshold=config["enrichment_analysis"]["log2fc_threshold"],
    script:
        "../scripts/network/compute_enrichment.py"


rule aggregate_enrichment:
    input:
        tfsage=expand(
            config["results_dir"]
            + "network/{threshold}/{linkage}/enrichment/tfsage/{query_id}.json",
            threshold=config["thresholds"],
            linkage=config["linkages"],
            query_id=config["query_ids"],
        ),
        motif_scan=expand(
            config["results_dir"]
            + "network/{threshold}/{linkage}/enrichment/motif_scan/{motif_db}/{query_id}.json",
            threshold=config["thresholds"],
            linkage=config["linkages"],
            motif_db=config["motif_scan"]["motif_dbs"],
            query_id=config["query_ids"],
        ),
        decoupler=expand(
            config["results_dir"]
            + "network/curated/enrichment/decoupler/{dc_network}.json",
            dc_network=config["decoupler"]["dc_networks"],
        ),
    output:
        config["results_dir"] + "network/aggregate_enrichment.csv",
    script:
        "../scripts/network/aggregate_enrichment.py"
