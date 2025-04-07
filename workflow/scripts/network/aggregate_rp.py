import pandas as pd
import json
from snakemake.script import snakemake


def load_json_targets(file):
    with open(file) as f:
        data = json.load(f)

    res = (data["factor"], data["target_genes"])
    return res


def aggregate_rp(input_files, output_file):
    # Step 1: Load all factors and their target genes
    factor_to_genes = {}
    for file in input_files:
        factor, genes = load_json_targets(file)
        factor_to_genes[factor] = genes

    # Step 2: Make dataframe
    df = pd.DataFrame(factor_to_genes)

    # Step 3: Save to parquet
    df.to_parquet(output_file)


aggregate_rp(snakemake.input, snakemake.output[0])
