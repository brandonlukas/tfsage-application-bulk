import pandas as pd
from tfsage.embedding import generate_embeddings
from snakemake.script import snakemake


def embeddings(input, output, params):
    rp_matrix = pd.read_parquet(input.rp_matrix)
    metadata = pd.read_csv(input.metadata).set_index("ID")

    embeddings = generate_embeddings(
        rp_matrix_df=rp_matrix,
        metadata_df=metadata,
        align_key=params.align_key,
    )
    embeddings.to_parquet(output[0])


embeddings(snakemake.input, snakemake.output, snakemake.params)
