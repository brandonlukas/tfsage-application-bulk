import pandas as pd
from snakemake.script import snakemake


def aggregate_genes(input_file, output_file):
    # Aggregate gene-level data
    df = pd.read_parquet(input_file)
    df["gene"] = df.index.str.split(":").str[1]
    df = df.groupby("gene").mean()
    df.index.name = None

    # Save aggregated data
    df.to_parquet(output_file)


aggregate_genes(snakemake.input[0], snakemake.output[0])
