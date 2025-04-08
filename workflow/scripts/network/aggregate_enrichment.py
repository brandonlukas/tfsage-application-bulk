import pandas as pd
import json
from snakemake.script import snakemake


def json_to_dataframe(json_file):
    with open(json_file) as f:
        js = json.load(f)

    df = (
        pd.DataFrame.from_records(js["enrichment"])
        .T.reset_index(names="factor")
        .assign(
            method_class=js["method_class"],
            method_name=js["method_name"],
            threshold=js["threshold"],
            linkage=js["linkage"],
            query_id=js["query_id"],
        )
    )
    return df


def aggregate_enrichment(json_files, output_file):
    df_list = [json_to_dataframe(f) for f in json_files]
    df = pd.concat(df_list, ignore_index=True)
    df.to_csv(output_file, index=False)


aggregate_enrichment(snakemake.input, snakemake.output[0])
