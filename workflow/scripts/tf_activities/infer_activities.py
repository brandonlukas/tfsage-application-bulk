import decoupler as dc
import pandas as pd
import json
from snakemake.script import snakemake


def infer_activities(
    input_gene_expression, input_network, output_file, params, wildcards
):
    sample_id = params.get("sample_id", None)
    method_class = params.get("method_class", None)
    method_name = params.get("method_name", None)
    mode = wildcards.get("mode", None)
    threshold = wildcards.get("threshold", None)

    net = load_network(input_network, method_name, params)
    mat = (
        pd.read_csv(input_gene_expression)
        .filter(["gene_name", sample_id], axis=1)
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
        "mode": mode,
        "threshold": threshold,
        "sample_id": sample_id,
        "tf_activities": acts.T.to_dict()[sample_id] if acts is not None else None,
    }

    # Save results to file
    with open(output_file, "w") as f:
        json.dump(results, f, indent=4)


def load_network(input_network, method_name, params):
    if input_network is not None:
        network = pd.read_parquet(input_network)
        net = (
            network.reset_index()
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
    input_gene_expression=snakemake.input["gene_expression"],
    input_network=getattr(snakemake.input, "network", None),
    output_file=snakemake.output[0],
    params=snakemake.params,
    wildcards=snakemake.wildcards,
)
