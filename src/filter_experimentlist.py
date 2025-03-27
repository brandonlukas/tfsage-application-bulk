import pandas as pd
import dask.dataframe as dd

# (Download) and table schema
# https://github.com/inutano/chip-atlas/wiki#experimentList_schema

max_cols = 70  # pick a safe upper bound
column_names = list(range(max_cols))

ddf = dd.read_csv(
    "resources/experimentList.tab",
    sep="\t",
    header=None,
    names=column_names,
    dtype=str,
    on_bad_lines="warn",  # will still warn, but won't drop the lines
    assume_missing=True,  # treats all columns as potentially missing
)

ddf = ddf[ddf[1] == "hg38"]
ddf = ddf[ddf[2].isin(["TFs and others", "Histone", "ATAC-Seq"])]
ddf = ddf[(ddf[2] != "Histone") | (ddf[3] == "H3K27ac")]
df = ddf.compute()

# Define the conditions
condition1 = ((df[2] == "Histone") | (df[2] == "ATAC-Seq")) & (df[4] == "Uterus")

factor_list = ["MED12", "FOS", "JUN", "CDK8", "PGR", "ESR1"]
condition2 = (df[2] == "TFs and others") & (df[3].isin(factor_list))

# Combine them
filtered_df = df[condition1 | condition2]

# Drop columns that are completely empty
filtered_df = filtered_df.dropna(axis=1, how="all")

filtered_df.reset_index(drop=True, inplace=True)

# Schema-based column names (from the Chip-Atlas wiki)
schema_columns = [
    "ID",  # 0
    "Assembly",  # 1
    "TrackClass",  # 2
    "TrackType",  # 3
    "CellClass",  # 4
    "CellType",  # 5
    "Description",  # 6
    "ProcessingLogs",  # 7
    "Title",  # 8
    "MetaData",  # 9 — everything beyond column 8 gets joined into this
]

# Make sure we have at least 9 columns before merging
num_cols = filtered_df.shape[1]

if num_cols > 9:
    # Join columns 9+ into a single MetaData string
    metadata_cols = filtered_df.columns[9:]
    filtered_df["MetaData"] = (
        filtered_df[metadata_cols]
        .astype(str)
        .apply(
            lambda row: "\t".join(
                [val for val in row if val not in ["", "nan", "<NA>", "None"]]
            ),
            axis=1,
        )
    )
    # Drop the original metadata columns
    filtered_df = filtered_df.drop(columns=metadata_cols)
else:
    # If there's no extra metadata, just make sure there's a MetaData column
    filtered_df["MetaData"] = ""

# Now rename columns
filtered_df.columns = schema_columns

# Save with header!
filtered_df.to_csv("resources/metadata.csv", index=False)
