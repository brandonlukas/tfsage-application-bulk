rule find_targets_tfsage:
    input:
        predictions=config["results_dir"]
        + "generate/{threshold}/tfsage_experiments/head_{n}/{query_id}_{factor}.parquet",
        query_file=config["results_dir"] + "downloads/{threshold}/{query_id}.bed",
    output:
        config["results_dir"]
        + "network/{threshold}/target_genes/{linkage_name}/tfsage/head_{n}/{query_id}_{factor}.json",
    params:
        genome="hg38",
        linkage_config=lambda w: config["network"].get(w.linkage_name),
        method_class="tfsage",
        method_name=lambda w: f"TFSage-{w.n}",
    script:
        "../scripts/network/find_targets.py"


rule find_targets_motif_scan:
    input:
        predictions=config["results_dir"]
        + "motif_scan/{threshold}/{motif_db}/{query_id}.bed",
        m2f_path="resources/motif_databases_filtered/{motif_db}_filtered.motif2factors.txt",
        query_file=config["results_dir"] + "downloads/{threshold}/{query_id}.bed",
    output:
        config["results_dir"]
        + "network/{threshold}/target_genes/{linkage_name}/motif_scan/{motif_db}/{query_id}_{factor}.json",
    params:
        genome="hg38",
        linkage_config=lambda w: config["network"].get(w.linkage_name),
        method_class="motif scan",
        method_name="{motif_db}",
    script:
        "../scripts/network/find_targets.py"


rule enrichment_analysis:
    input:
        tfsage=expand(
            config["results_dir"]
            + "network/{threshold}/target_genes/{linkage_name}/tfsage/head_{n}/{query_id}_{factor}.json",
            threshold=config["thresholds"],
            linkage_name=config["network"].keys(),
            n=config["tfsage_n"],
            query_id=query_set,
            factor=factor_list,
        ),
        motif_scan=expand(
            config["results_dir"]
            + "network/{threshold}/target_genes/{linkage_name}/motif_scan/{motif_db}/{query_id}_{factor}.json",
            threshold=config["thresholds"],
            linkage_name=config["network"].keys(),
            motif_db=config["motif_dbs"],
            query_id=query_set,
            factor=factor_list,
        ),
    output:
        config["results_dir"] + "network/enrichment_analysis.csv",
    params:
        de_path=config["deseq2_report"],
        padj_threshold=0.05,
        log2fc_threshold=1,
    script:
        "../scripts/network/enrichment_analysis.py"
