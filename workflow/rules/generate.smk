rule presearch:
    input:
        distances=config["results_dir"] + "data/{threshold}/distances.parquet",
        metadata=config["metadata"],
    output:
        config["results_dir"] + "presearch/{threshold}/results.parquet",
    params:
        query_ids=config["query_ids"],
        n=config["tfsage"]["n"],
    script:
        "../scripts/generate/presearch.py"


rule synthesize_experiments:
    input:
        config["results_dir"] + "presearch/{threshold}/results.parquet",
    output:
        config["results_dir"] + "generate/{threshold}/{query_id}_{factor}.parquet",
    params:
        data_dir=config["downloads_dir"] + "{threshold}/",
    script:
        "../scripts/generate/synthesize_experiments.py"
