rule synthesize_experiments:
    input:
        distances=config["results_dir"] + "data/{threshold}/distances.parquet",
        metadata=config["chip_metadata"],
    output:
        config["results_dir"]
        + "generate/{threshold}/tfsage_experiments/head_{n}/{query_id}_{factor}.parquet",
    params:
        data_dir=lambda w: config["results_dir"] + f"downloads/{w.threshold}/",
    script:
        "../scripts/generate/synthesize_experiments.py"
