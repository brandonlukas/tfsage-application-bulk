import pandas as pd
import numpy as np
import json
from scipy.stats import fisher_exact
from snakemake.script import snakemake


def enrichment_analysis(json_files, output_file, params):
    # --- 1. Load DESeq2 Report ---
    de_path = params.get("deseq2_results")
    de_df = pd.read_csv(de_path)

    # Clean DE table: remove rows with NA in padj
    de_df["is_de"] = (
        (~de_df["padj"].isna())
        & (de_df["padj"] < params.padj_threshold)
        & (de_df["log2FoldChange"].abs() > params.log2fc_threshold)
    )

    # Make gene_name lookup set
    all_genes = set(de_df["gene_name"])
    de_genes = set(de_df.loc[de_df["is_de"], "gene_name"])

    # --- 2. Loop through JSONs and test ---
    records = []
    for path in json_files:
        with open(path) as f:
            js = json.load(f)

        method_class = js["method_class"]
        method_name = js["method_name"]
        threshold = js.get("threshold")
        query_id = js["query_id"]
        factor = js["factor"]
        linkage_name = js["linkage_name"]
        target_genes = list(js["target_genes"].keys())

        enrich = test_enrichment(target_genes, all_genes, de_genes)
        record = {
            "method_class": method_class,
            "method_name": method_name,
            "threshold": threshold,
            "query_id": query_id,
            "factor": factor,
            "linkage_name": linkage_name,
            **enrich,
        }

        records.append(record)

    # --- 3. Combine to final results dataframe ---
    results_df = pd.DataFrame.from_records(records)

    # Save to disk
    results_df.to_csv(output_file, index=False)


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


enrichment_analysis(
    snakemake.input,
    snakemake.output[0],
    snakemake.params,
)
