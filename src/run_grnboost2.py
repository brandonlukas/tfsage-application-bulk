import pandas as pd
from arboreto.algo import grnboost2


def get_factor_list():
    import dask.dataframe as dd

    factor_list = (
        dd.read_csv("resources/metadata.csv")
        .query("TrackClass == 'TFs and others'")["TrackType"]
        .unique()
        .compute()
        .tolist()
    )
    return factor_list


df = pd.read_csv("data/GSE128242_counts.csv")
factor_list = get_factor_list()

network = grnboost2(
    expression_data=df.drop(columns=["gene_id", "gene_name"]).values.T,
    gene_names=df.gene_name,
    tf_names=factor_list,
    seed=42,
    verbose=True,
)

network.rename(
    columns={
        "TF": "source",
        "importance": "weight",
    }
).to_csv("data/GSE128242_grnboost2.parquet", index=False)
