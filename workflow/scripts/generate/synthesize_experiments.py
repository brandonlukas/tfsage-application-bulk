import pandas as pd
import numpy as np
from tfsage.search import generate_ranked_list
from tfsage import generate
from snakemake.script import snakemake


def synthesize_experiments(
    distances_file, metadata_file, output_file, params, wildcards
):
    # Read distances
    distances_df = pd.read_parquet(distances_file)

    # Generate scoring function
    p = distances_df.to_numpy().var()
    scoring_function = lambda x: np.exp(-(x**2) / p)

    # Read metadata
    metadata = pd.read_csv(metadata_file).set_index("ID")

    # Generate ranked list for the query, removing the cell type
    ranked_list = (
        generate_ranked_list(
            distances_df,
            wildcards.query_id,
            metadata=metadata,
            scoring_function=scoring_function,
        )
        .query("TrackType == @wildcards.factor")
        .filter(["score", "TrackType", "CellClass", "CellType"], axis=1)
        .sort_values("score", ascending=False)
        .head(params.n)
    )

    bed_files = [params.data_dir + x + ".bed" for x in ranked_list.index]
    weights = ranked_list["score"].tolist()
    result = generate_result(bed_files, weights)
    result.to_parquet(output_file)


def generate_result(bed_files, weights):
    if len(bed_files) == 1:
        # "Hack" to make the function work with only one bed file
        bed_files = [bed_files[0], bed_files[0]]
        df = generate.synthesize_experiments(bed_files, weights).drop("file_1", axis=1)
    else:
        df = generate.synthesize_experiments(bed_files, weights)
    return df


synthesize_experiments(
    snakemake.input.distances,
    snakemake.input.metadata,
    snakemake.output[0],
    snakemake.params,
    snakemake.wildcards,
)
