import pandas as pd
import json
from scipy.stats import fisher_exact
from snakemake.script import snakemake


def compute_enrichment(input_file, output_file, params, wildcards):
    # --- 1. Load DESeq2 Report ---
    all_genes, de_genes = get_gene_sets(params)

    # --- 2. Loop through each source and test ---
    if input_file.endswith(".parquet"):
        net = pd.read_parquet(input_file)
    elif input_file.endswith(".csv"):
        net = pd.read_csv(input_file)
    else:
        raise ValueError(f"Unsupported file type: {input_file}")

    records = {}
    for source in net["source"].unique():
        target_genes = net.query("source == @source & weight != 0")["target"].tolist()

        enrichment = test_enrichment(target_genes, all_genes, de_genes)
        records[source] = enrichment

    # --- 3. Combine to final results dataframe ---
    results = {
        "method_class": params.method_class,
        "method_name": params.method_name,
        "threshold": wildcards.get("threshold", None),
        "linkage": wildcards.get("linkage", None),
        "query_id": wildcards.get("query_id", None),
        "enrichment": records,
    }

    # Save results to file
    with open(output_file, "w") as f:
        json.dump(results, f, indent=4)


def get_gene_sets(params):
    df = pd.read_csv(params.de_res)

    # Clean DE table: remove rows with NA in padj
    df["is_de"] = (
        (~df["padj"].isna())
        & (df["padj"] < params.padj_threshold)
        & (df["log2FoldChange"].abs() > params.log2fc_threshold)
    )

    # Make gene_name lookup set
    all_genes = set(df["gene_name"])
    de_genes = set(df.loc[df["is_de"], "gene_name"])
    return all_genes, de_genes


def test_enrichment(target_genes, all_genes, de_genes):
    target_set = set(target_genes) & all_genes
    non_target_set = all_genes - target_set

    a = len(target_set & de_genes)  # target & DE
    b = len(target_set) - a  # target & not DE
    c = len(non_target_set & de_genes)  # non-target & DE
    d = len(non_target_set) - c  # non-target & not DE

    table = [[a, b], [c, d]]
    oddsratio, pvalue = fisher_exact(table, alternative="greater")

    return {
        "num_targets": len(target_set),
        "num_de_genes": len(de_genes),
        "num_de_targets": a,
        "odds_ratio": oddsratio,
        "pvalue": pvalue,
    }


compute_enrichment(
    snakemake.input[0], snakemake.output[0], snakemake.params, snakemake.wildcards
)
