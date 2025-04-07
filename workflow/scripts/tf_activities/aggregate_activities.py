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

        # Extract metadata
        meta = {
            "method_class": data.get("method_class"),
            "method_name": data.get("method_name"),
            "mode": data.get("mode"),
            "threshold": data.get("threshold"),
            "sample_id": data.get("sample_id"),
        }

        # Extract TF activities and create one row per TF
        for tf_name, activity in tf_activities.items():
            row = meta.copy()
            row["factor"] = tf_name
            row["ulm_score"] = activity
            records.append(row)

    # Create the final DataFrame
    df = pd.DataFrame(records)
    df.to_parquet(output_file, index=False)


aggregate_activities(snakemake.input, snakemake.output[0])
