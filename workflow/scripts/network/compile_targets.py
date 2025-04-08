import pandas as pd
import json
from snakemake.script import snakemake


def compile_targets(json_files, output_file):
    df_list = [json_to_net(f) for f in json_files]
    df = pd.concat(df_list, ignore_index=True)
    df.to_parquet(output_file, index=False)


def json_to_net(json_file):
    with open(json_file) as f:
        js = json.load(f)

    # Select target genes with weight != 0
    target_genes = {k: v for k, v in js["target_genes"].items() if v != 0}
    df = (
        pd.Series(target_genes)
        .to_frame(name="weight")
        .reset_index(names="target")
        .assign(
            method_class=js["method_class"],
            method_name=js["method_name"],
            threshold=js["threshold"],
            linkage_name=js["linkage_name"],
            query_id=js["query_id"],
            source=js["factor"],
        )
    )

    str_columns = df.dtypes[df.dtypes == "object"].index
    df[str_columns] = df[str_columns].astype("string")

    column_order = ["source", "target", "weight", *df.columns]
    column_order = list(dict.fromkeys(column_order))
    df = df.filter(column_order, axis=1)
    return df


compile_targets(snakemake.input, snakemake.output[0])
