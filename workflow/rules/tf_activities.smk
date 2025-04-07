rule infer_activities_tfsage:
    input:
        expression="data/GSE128242_counts.csv",
        network=config["results_dir"]
        + "network/{threshold}/networks/tfsage/head_{n}/{query_id}.parquet",
    output:
        config["results_dir"]
        + "tf_activities/{threshold}/tfsage/head_{n}/{query_id}.json",
    params:
        method_class="tfsage",
        method_name=lambda w: f"TFSage-{w.n}",
    conda:
        "decoupler_env"
    script:
        "../scripts/tf_activities/infer_activities.py"


rule infer_activities_motif_scan:
    input:
        expression="data/GSE128242_counts.csv",
        network=config["results_dir"]
        + "network/{threshold}/networks/motif_scan/{motif_db}/{query_id}.parquet",
    output:
        config["results_dir"]
        + "tf_activities/{threshold}/motif_scan/{motif_db}/{query_id}.json",
    params:
        method_class="motif scan",
        method_name="{motif_db}",
    conda:
        "decoupler_env"
    script:
        "../scripts/tf_activities/infer_activities.py"


rule infer_activities_curated:
    input:
        expression="data/GSE128242_counts.csv",
    output:
        config["results_dir"] + "tf_activities/curated/{curated_db}/{query_id}.json",
    params:
        method_class="curated",
        method_name="{curated_db}",
        factor_list=factor_list,
    conda:
        "decoupler_env"
    script:
        "../scripts/tf_activities/infer_activities.py"


rule aggregate_activities:
    input:
        tfsage=expand(
            config["results_dir"]
            + "tf_activities/{threshold}/tfsage/head_{n}/{query_id}.json",
            query_id=query_set,
            n=config["tfsage_n"],
            threshold=config["thresholds"],
        ),
        motif_scan=expand(
            config["results_dir"]
            + "tf_activities/{threshold}/motif_scan/{motif_db}/{query_id}.json",
            query_id=query_set,
            motif_db=config["motif_dbs"],
            threshold=config["thresholds"],
        ),
        curated=expand(
            config["results_dir"]
            + "tf_activities/curated/{curated_db}/{query_id}.json",
            query_id=query_set,
            curated_db=config["curated_dbs"],
        ),
    output:
        config["results_dir"] + "tf_activities/results.parquet",
    script:
        "../scripts/tf_activities/aggregate_activities.py"
