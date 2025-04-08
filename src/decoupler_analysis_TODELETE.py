import pandas as pd
import decoupler as dc
from glob import glob
import json
from tqdm import tqdm

with open("resources/query_set.txt") as f:
    query_set = [line.strip() for line in f.readlines()]

thresholds = ["05", "10", "20"]


def parse_json(file):
    with open(file) as f:
        data = json.load(f)

    source = data["factor"]
    target_dict = data["target_genes"]
    return source, target_dict


def build_network(files):
    data = []
    for file in files:
        source, target_dict = parse_json(file)
        for target, weight in target_dict.items():
            data.append((source, target, weight))

    # Build long-form DataFrame
    net = pd.DataFrame(data, columns=["source", "target", "weight"])
    return net


def estimate_activities(files, mat):
    net = build_network(files)
    (estimate, *_) = dc.run_aucell(mat, net)
    df = (
        estimate.reset_index()
        .melt(id_vars="index")
        .rename(columns={"index": "sample_id", "value": "aucell"})
    )
    return df


counts = pd.read_csv("data/GSE128242_counts.csv")
mat = counts.drop(columns="gene_id").groupby("gene_name").sum().T

# Use DE statistic instead of gene expression
res = pd.read_csv("data/GSE128242_deseq2_report.csv")
mat = res[["gene_name", "stat"]].dropna().groupby("gene_name").mean().T
# de_genes = res.query("padj < 0.05 & abs(log2FoldChange) > 1")["gene_name"].unique()
# mat = mat[de_genes]

pattern = "/Users/brandonlukas/Library/CloudStorage/Box-Box/data/tfsage/application-bulk/results/network/{threshold}/target_genes/linkage1/tfsage/head_15/{query_id}_*.json"
df_list = []
for query_id in tqdm(query_set):
    for threshold in thresholds:
        files = glob(pattern.format(query_id=query_id, threshold=threshold))
        df = estimate_activities(files, mat).assign(
            method_class="tfsage",
            method_name="TFSage-15",
            threshold=threshold,
            query_id=query_id,
            linkage_name="linkage1",
        )
        df_list.append(df)

df = pd.concat(df_list, ignore_index=True)
df.to_csv("test.csv", index=False)
