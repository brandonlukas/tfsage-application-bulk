import json
import pandas as pd
from snakemake.script import snakemake


def aggregate_activities(input_files, output_file):
    records = []

    for file in input_files:
        with open(file) as f:
            data = json.load(f)

        tf_activities = data.get("tf_activities")
        if tf_activities is None:
            continue

        meta = {
            "method_class": data.get("method_class"),
            "method_name": data.get("method_name"),
            "query_id": data.get("query_id"),
            "threshold": data.get("threshold"),
        }

        for sample_id, tf_scores in tf_activities.items():
            for tf_name, score in tf_scores.items():
                row = meta.copy()
                row.update({"sample_id": sample_id, "factor": tf_name, "score": score})
                records.append(row)

    df = pd.DataFrame(records)
    df.to_parquet(output_file, index=False)


aggregate_activities(snakemake.input, snakemake.output[0])
