import pandas as pd
import json
from snakemake.script import snakemake


def load_json_targets(file):
    with open(file) as f:
        data = json.load(f)
    res = (data["factor"], data["target_genes"])
    return res


def aggregate_target_genes(input_files, output_file):
    # Step 1: Load all factors and their target genes
    factor_to_genes = {}
    for file in input_files:
        factor, genes = load_json_targets(file)
        factor_to_genes[factor] = genes

    # Step 2: Get a unique list of all genes
    all_genes = set(gene for genes in factor_to_genes.values() for gene in genes)

    # Step 3: Initialize DataFrame (rows=genes, cols=factors) with 0
    df = pd.DataFrame(
        0,
        index=sorted(all_genes),
        columns=sorted(
            factor_to_genes.keys(),
        ),
    )

    # Step 4: Fill in 1s where gene is a target of a factor
    for factor, genes in factor_to_genes.items():
        df.loc[genes, factor] = 1

    df.to_parquet(output_file)


aggregate_target_genes(snakemake.input, snakemake.output[0])
