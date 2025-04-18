import pandas as pd
from tfsage.download import download_chip_atlas
from tfsage.search import compute_distances


rule download:
    output:
        config["downloads_dir"] + "{threshold}/{experiment}.bed",
    retries: 5
    run:
        download_chip_atlas(wildcards.experiment, output[0], wildcards.threshold)


def all_experiments(wildcards):
    experiments = pd.read_csv(config["metadata"])["ID"].tolist()
    template = config["downloads_dir"] + f"{wildcards.threshold}/{{experiment}}.bed"
    return expand(template, experiment=experiments)


rule features_tss:
    input:
        all_experiments,
    output:
        config["results_dir"] + "data/{threshold}/rp_matrix/tss.parquet",
    threads: 48
    resources:
        mem_mb=96000,
    script:
        "../scripts/prepare_data/features_tss.py"


rule features_gene:
    input:
        config["results_dir"] + "data/{threshold}/rp_matrix/tss.parquet",
    output:
        config["results_dir"] + "data/{threshold}/rp_matrix/gene.parquet",
    script:
        "../scripts/prepare_data/features_gene.py"


rule embeddings:
    input:
        rp_matrix=config["results_dir"] + "data/{threshold}/rp_matrix/gene.parquet",
        metadata=config["metadata"],
    output:
        config["results_dir"] + "data/{threshold}/embeddings.parquet",
    params:
        align_key="TrackClass",
    script:
        "../scripts/prepare_data/embeddings.py"


rule distances:
    input:
        config["results_dir"] + "data/{threshold}/embeddings.parquet",
    output:
        config["results_dir"] + "data/{threshold}/distances.parquet",
    run:
        embeddings = pd.read_parquet(input[0])
        distances = compute_distances(embeddings)
        distances.to_parquet(output[0])
