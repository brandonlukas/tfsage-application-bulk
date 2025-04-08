# Download file and refer to table schema:
# https://github.com/inutano/chip-atlas/wiki#experimentList_schema
import dask.dataframe as dd


def read_file(file_path, max_cols=70):
    """
    Read the experimentList.tab file with Dask and return a DataFrame.
    The file is tab-separated and has no header. We will set the column names
    to a list of integers from 0 to max_cols-1. We will also set the dtype to str
    to avoid issues with parsing. We will set on_bad_lines to "warn" to avoid
    dropping lines and assume_missing to True to treat all columns as potentially missing.
    """
    column_names = list(range(max_cols))
    ddf = dd.read_csv(
        file_path,
        sep="\t",
        header=None,
        names=column_names,
        dtype=str,
        on_bad_lines="warn",
        assume_missing=True,
    )

    # Basic filtering
    ddf = ddf[ddf[1] == "hg38"]
    ddf = ddf[ddf[2].isin(["TFs and others", "Histone", "ATAC-Seq"])]
    ddf = ddf[(ddf[2] != "Histone") | (ddf[3] == "H3K27ac")]
    df = ddf.compute()
    return df


def match_chip_atlas_schema(
    df,
    schema_columns=[
        "ID",
        "Assembly",
        "TrackClass",
        "TrackType",
        "CellClass",
        "CellType",
        "Description",
        "ProcessingLogs",
        "Title",
        "MetaData",
    ],
):
    """
    Rename columns to match the schema.
    Schema-based column names (from the Chip-Atlas wiki)
    """
    # Make sure we have at least 9 columns before merging
    num_cols = df.shape[1]

    if num_cols > 9:
        # Join columns 9+ into a single MetaData string
        metadata_cols = df.columns[9:]
        df["MetaData"] = (
            df[metadata_cols]
            .astype(str)
            .apply(
                lambda row: "\t".join(
                    [val for val in row if val not in ["", "nan", "<NA>", "None"]]
                ),
                axis=1,
            )
        )
        # Drop the original metadata columns
        df = df.drop(columns=metadata_cols)
    else:
        # If there's no extra metadata, just make sure there's a MetaData column
        df["MetaData"] = ""

    # Now rename columns
    df.columns = schema_columns
    return df


df = read_file("resources/experimentList.tab")

# Include relevant H3K27ac and ATAC-seq experiments for Uterus
condition1 = ((df[2] == "Histone") | (df[2] == "ATAC-Seq")) & (df[4] == "Uterus")

# Include almost all TFs (exclude CTCF and Epitope tags?)
condition2 = df[2] == "TFs and others"  # & (~df[3].isin(["CTCF", "Epitope tags"]))

# Filter the DataFrame based on the conditions
df = df[condition1 | condition2]

# Drop columns that are completely empty and reset the index
df = df.dropna(axis=1, how="all").reset_index(drop=True)

# Rename columns to match the schema
df = match_chip_atlas_schema(df)

# Save with header!
df.to_csv("resources/metadata.csv", index=False)
