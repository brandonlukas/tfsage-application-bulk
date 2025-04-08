import dask.dataframe as dd
import pandas as pd
import numpy as np
from snakemake.script import snakemake


def presearch(distances_file, metadata_file, output_file, params):
    metadata = pd.read_csv(metadata_file)
    ddf = dd.read_parquet(distances_file)

    # Compute scoring function
    p = ddf.melt().value.var().compute()
    scoring_function = lambda x: np.exp(-(x**2) / p)

    df = (
        ddf[params.query_ids]
        .reset_index()
        .melt(id_vars="index")
        .astype({"variable": "string[pyarrow]"})
        .merge(metadata, left_on="index", right_on="ID")
        .query("TrackClass == 'TFs and others'")
        .rename(
            columns={
                "index": "experiment_id",
                "variable": "query_id",
                "value": "distance",
            },
        )
        .groupby(["query_id", "TrackType"])
        .apply(lambda df: df.nsmallest(params.n, "distance"))
        .assign(score=lambda df: scoring_function(df.distance))
        .reset_index(drop=True)
        .compute()
    )

    df.to_parquet(output_file, index=False)


presearch(
    snakemake.input.distances,
    snakemake.input.metadata,
    snakemake.output[0],
    snakemake.params,
)
