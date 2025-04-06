import pandas as pd

benchmark_df = pd.read_csv(config["benchmark_set"])


rule synthesize_experiments:
    input:
        distances=config["results_dir"]
        + "data/{threshold}/distances/exclude_holdout.parquet",
        metadata=config["chip_metadata"],
    output:
        config["results_dir"]
        + "generate/{threshold}/tfsage_experiments/head_{n}/{query_id}_{factor}.parquet",
    params:
        data_dir=lambda w: config["results_dir"] + f"downloads/{w.threshold}/",
    script:
        "../scripts/generate/synthesize_experiments.py"


rule generate_test_set:
    input:
        query_file=config["results_dir"] + "downloads/{threshold}/{query_id}.bed",
        target_file=config["results_dir"] + "downloads/{threshold}/{target_id}.bed",
    output:
        config["results_dir"]
        + "generate/{threshold}/test_sets/{test_set_name}/{query_id}_{target_id}.parquet",
    params:
        config_test_set=lambda w: config["test_sets"][w.test_set_name],
    script:
        "../scripts/generate/generate_test_set.py"


rule compute_metrics_tfsage:
    input:
        predictions=config["results_dir"]
        + "generate/{threshold}/tfsage_experiments/head_{n}/{query_id}_{factor}.parquet",
        test_set=config["results_dir"]
        + "generate/{threshold}/test_sets/{test_set_name}/{query_id}_{target_id}.parquet",
    output:
        config["results_dir"]
        + "generate/{threshold}/metrics/{test_set_name}/tfsage/head_{n}/{query_id}_{factor}_{target_id}.json",
    params:
        method_class="tfsage",
        method_name=lambda w: f"TFSage-{w.n}",
        query_assay=lambda w: benchmark_df.query("ID_query == @w.query_id").iloc[0][
            "TrackType_query"
        ],
        cell_type=lambda w: benchmark_df.query("ID_query == @w.query_id").iloc[0][
            "CellType"
        ],
    script:
        "../scripts/generate/compute_classification_metrics.py"


rule compute_metrics_motif_scan:
    input:
        predictions=config["results_dir"]
        + "motif_scan/{threshold}/{motif_db}/{query_id}.bed",
        test_set=config["results_dir"]
        + "generate/{threshold}/test_sets/{test_set_name}/{query_id}_{target_id}.parquet",
        m2f_path="resources/motif_databases_filtered/{motif_db}_filtered.motif2factors.txt",
    output:
        config["results_dir"]
        + "generate/{threshold}/metrics/{test_set_name}/motif_scan/{motif_db}/{query_id}_{factor}_{target_id}.json",
    params:
        method_class="motif scan",
        method_name="{motif_db}",
        query_assay=lambda w: benchmark_df.query("ID_query == @w.query_id").iloc[0][
            "TrackType_query"
        ],
        cell_type=lambda w: benchmark_df.query("ID_query == @w.query_id").iloc[0][
            "CellType"
        ],
    script:
        "../scripts/generate/compute_classification_metrics.py"


rule aggregate_metrics_generate:
    input:
        all_tfsage=expand(
            expand(
                config["results_dir"]
                + "generate/{threshold}/metrics/{test_set_name}/tfsage/head_{n}/{{query_id}}_{{factor}}_{{target_id}}.json",
                n=config["tfsage_n"],
                test_set_name=list(config["test_sets"].keys()),
                threshold=config["thresholds"],
            ),
            zip,
            query_id=benchmark_df["ID_query"],
            target_id=benchmark_df["ID_target"],
            factor=benchmark_df["TrackType_target"],
        ),
        all_motif_scan=expand(
            expand(
                config["results_dir"]
                + "generate/{threshold}/metrics/{test_set_name}/motif_scan/{motif_db}/{{query_id}}_{{factor}}_{{target_id}}.json",
                motif_db=config["motif_dbs"],
                test_set_name=list(config["test_sets"].keys()),
                threshold=config["thresholds"],
            ),
            zip,
            query_id=benchmark_df["ID_query"],
            target_id=benchmark_df["ID_target"],
            factor=benchmark_df["TrackType_target"],
        ),
    output:
        config["results_dir"] + "generate/benchmark_metrics.csv",
    script:
        "../scripts/generate/aggregate_metrics.py"
