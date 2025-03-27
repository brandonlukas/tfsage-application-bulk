import pandas as pd
from tfsage.download import download_chip_atlas
from tfsage.search import compute_distances


rule download:
    output:
        config["results_dir"] + "downloads/{threshold}/{experiment}.bed",
    retries: 5
    params:
        genome="hg38",
    run:
        download_chip_atlas(
            wildcards.experiment, output[0], wildcards.threshold, params.genome
        )


def all_experiments(wildcards):
    experiments = pd.read_csv(config["chip_metadata"])["ID"].tolist()
    template = (
        config["results_dir"] + f"downloads/{wildcards.threshold}/{{experiment}}.bed"
    )
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
        metadata=config["chip_metadata"],
    output:
        config["results_dir"] + "data/{threshold}/embeddings/{mode}.parquet",
    params:
        align_key="TrackClass",
        target_factors=config["target_factors"],
    script:
        "../scripts/prepare_data/embeddings.py"


rule distances:
    input:
        config["results_dir"] + "data/{threshold}/embeddings/{mode}.parquet",
    output:
        config["results_dir"] + "data/{threshold}/distances/{mode}.parquet",
    run:
        embeddings = pd.read_parquet(input[0])
        distances = compute_distances(embeddings)
        distances.to_parquet(output[0])
