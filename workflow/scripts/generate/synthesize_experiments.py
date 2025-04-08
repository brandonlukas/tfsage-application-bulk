import dask.dataframe as dd
from tfsage import generate
from snakemake.script import snakemake


def synthesize_experiments(input_file, output_file, params, wildcards):
    # Fetch ranked list from presearch data
    ranked_list = (
        dd.read_parquet(input_file)
        .query("query_id == @wildcards.query_id")
        .query("TrackType == @wildcards.factor")
        .compute()
    )

    bed_files = [params.data_dir + f"{x}.bed" for x in ranked_list["experiment_id"]]
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
    snakemake.input[0], snakemake.output[0], snakemake.params, snakemake.wildcards
)
