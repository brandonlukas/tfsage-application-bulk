import pandas as pd
from tfsage.embedding import generate_embeddings
from snakemake.script import snakemake


def embeddings(input, output, params, wildcards):
    rp_matrix = pd.read_parquet(input.rp_matrix)
    metadata = pd.read_csv(input.metadata).set_index("ID")

    if wildcards.mode == "exclude_holdout":
        # Build regex pattern: "LEIO-factor1|MYO-factor1|LEIO-factor2|MYO-factor2|..."
        pattern = "|".join(
            [f"LEIO-{factor}|MYO-{factor}" for factor in params.holdout_factors]
        )

        # Create boolean masks
        cond1 = metadata["TrackType"].isin(params.holdout_factors)
        cond2 = metadata["Title"].str.contains(pattern, na=False)

        # Apply filter
        metadata = metadata[~(cond1 & cond2)]
        rp_matrix = rp_matrix[metadata.index]

    embeddings = generate_embeddings(
        rp_matrix_df=rp_matrix,
        metadata_df=metadata,
        align_key=params.align_key,
    )
    embeddings.to_parquet(output[0])


embeddings(snakemake.input, snakemake.output, snakemake.params, snakemake.wildcards)
