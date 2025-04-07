import decoupler as dc
import pandas as pd
import json
from snakemake.script import snakemake


def infer_activities(input_expression, input_network, output_file, params, wildcards):
    method_class = params.get("method_class", None)
    method_name = params.get("method_name", None)
    threshold = wildcards.get("threshold", None)

    net = load_network(input_network, method_name, params)
    mat = (
        pd.read_csv(input_expression)
        .drop(columns=["gene_id"])
        .query("gene_name in @net['target'].unique()")
        .groupby("gene_name")
        .max()
        .transpose()
    )

    try:
        acts, *_ = dc.run_ulm(mat, net)
    except ValueError:
        acts = None

    # Combine everything into one dictionary
    results = {
        "method_class": method_class,
        "method_name": method_name,
        "threshold": threshold,
        "tf_activities": acts.T.to_dict() if acts is not None else None,
    }

    # Save results to file
    with open(output_file, "w") as f:
        json.dump(results, f, indent=4)


def load_network(input_network, method_name, params):
    if input_network is not None:
        net = (
            pd.read_parquet(input_network)
            .reset_index()
            .melt(id_vars="index")
            .rename(
                {"index": "target", "variable": "source", "value": "weight"}, axis=1
            )
            .filter(["source", "target", "weight"], axis=1)
            .query("weight > 0")
        )

    else:
        if method_name == "dorothea":
            net = dc.get_dorothea()
        elif method_name == "collectri":
            net = dc.get_collectri()

        factor_list = params.get("factor_list")
        net = net.query("source in @factor_list")

    return net


infer_activities(
    input_expression=snakemake.input["expression"],
    input_network=getattr(snakemake.input, "network", None),
    output_file=snakemake.output[0],
    params=snakemake.params,
    wildcards=snakemake.wildcards,
)
