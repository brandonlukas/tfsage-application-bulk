rule find_target_genes_tfsage:
    input:
        predictions=config["results_dir"]
        + "generate/{mode}/{threshold}/tfsage_experiments/head_{n}/{query_id}_{factor}.parquet",
    output:
        config["results_dir"]
        + "network/{threshold}/target_genes/tfsage/{mode}/head_{n}/{query_id}_{factor}.json",
    params:
        genome="hg38",
        network_config=config["network"],
        method_class="tfsage",
        method_name=lambda w: f"TFSage-{w.n}",
    script:
        "../scripts/network/find_target_genes.py"


rule find_target_genes_motif_scan:
    input:
        predictions=config["results_dir"]
        + "motif_scan/{threshold}/{motif_db}/{query_id}.bed",
        m2f_path="resources/motif_databases_filtered/{motif_db}_filtered.motif2factors.txt",
    output:
        config["results_dir"]
        + "network/{threshold}/target_genes/motif_scan/{motif_db}/{query_id}_{factor}.json",
    params:
        genome="hg38",
        network_config=config["network"],
        method_class="motif scan",
        method_name="{motif_db}",
    script:
        "../scripts/network/find_target_genes.py"


rule aggregate_target_genes_tfsage:
    input:
        expand(
            config["results_dir"]
            + "network/{threshold}/target_genes/tfsage/{mode}/head_{n}/{query_id}_{factor}.json",
            factor=factor_list,
            allow_missing=True,
        ),
    output:
        config["results_dir"]
        + "network/{threshold}/networks/tfsage/{mode}/head_{n}/{query_id}.parquet",
    script:
        "../scripts/network/aggregate_target_genes.py"


rule aggregate_target_genes_motif_scan:
    input:
        expand(
            config["results_dir"]
            + "network/{threshold}/target_genes/motif_scan/{motif_db}/{query_id}_{factor}.json",
            factor=factor_list,
            allow_missing=True,
        ),
    output:
        config["results_dir"]
        + "network/{threshold}/networks/motif_scan/{motif_db}/{query_id}.parquet",
    script:
        "../scripts/network/aggregate_target_genes.py"
