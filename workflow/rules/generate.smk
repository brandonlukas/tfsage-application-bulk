rule synthesize_experiments:
    input:
        distances=config["results_dir"] + "data/{threshold}/distances.parquet",
        metadata=config["metadata"],
    output:
        config["results_dir"] + "generate/{threshold}/{query_id}_{factor}.parquet",
    params:
        data_dir=config["downloads_dir"] + "{threshold}/",
        n=config["tfsage"]["n"],
    script:
        "../scripts/generate/synthesize_experiments.py"
