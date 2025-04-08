import pandas as pd
import dask.dataframe as dd
import json
import pybedtools
from tfsage.features import load_region_set, extract_features
from tfsage.utils import load_features_bed, closest_features
from snakemake.script import snakemake


def find_targets(
    predictions_file: str,
    query_file: str,
    output_file: str,
    params: dict,
    wildcards: dict,
    m2f_path: str | None = None,
):
    # Load predictions
    if params.method_class == "motif scan":
        predictions = load_predictions_motif_scan(
            predictions_file, wildcards.factor, m2f_path
        )
    elif params.method_class == "tfsage":
        predictions = pd.read_parquet(predictions_file).assign(score=lambda x: x["sum"])
    else:
        raise ValueError(f"Unknown method class: {params.method_class}")

    target_genes = {}
    if predictions is not None:
        predictions = predictions.query("score > 0")
        if not predictions.empty:
            predictions = intersect_with_query(predictions, query_file)
            if wildcards.linkage == "rp-grn":
                gene_loc_set = load_region_set()
                target_genes = rp_grn_target_genes(predictions, gene_loc_set)
            elif wildcards.linkage == "cellranger-arc":
                features_bed = load_features_bed()
                target_genes = cellranger_arc_target_genes(predictions, features_bed)

            # Remove NaN keys
            target_genes = remove_nan_keys(target_genes)

    # Combine everything into one dictionary
    results = {
        "method_class": params.method_class,
        "method_name": params.method_name,
        "threshold": wildcards.threshold,
        "linkage": wildcards.linkage,
        "query_id": wildcards.query_id,
        "factor": wildcards.factor,
        "target_genes": target_genes,
    }

    # Save results to file
    with open(output_file, "w") as f:
        json.dump(results, f, indent=4)


def remove_nan_keys(a_dict):
    # Remove keys that are NaN or not a string
    a_dict = {k: v for k, v in a_dict.items() if isinstance(k, str) and k != "NaN"}
    return a_dict


def load_predictions_motif_scan(
    predictions_file, factor, m2f_path
) -> pd.DataFrame | None:
    motif_list = (
        pd.read_csv(m2f_path, sep="\t").query("Factor == @factor")["Motif"].tolist()
    )

    if not len(motif_list):
        return None

    ddf = dd.read_csv(predictions_file, sep="\t", header=None, comment="#")
    ddf = ddf[ddf[3].isin(motif_list)]
    predictions = ddf.compute()
    predictions.columns = ["chrom", "start", "end", "motif", "score", "strand"]
    return predictions


def intersect_with_query(predictions, query_file):
    bed_file = pybedtools.BedTool.from_dataframe(predictions)
    query_bed = pybedtools.BedTool(query_file)
    bed_file = bed_file.intersect(query_bed, u=True)

    predictions = bed_file.to_dataframe(
        disable_auto_names=True, names=predictions.columns
    )
    return predictions


def cellranger_arc_target_genes(predictions, features_bed):
    bed_file = pybedtools.BedTool.from_dataframe(predictions).sort()
    df = closest_features(bed_file, features_bed)
    df = df.query("-1000 < distance < 100")
    df["gene"] = df["name"].str.split(":").str[1]
    target_genes = {gene: 1.0 for gene in df["gene"].unique()}
    return target_genes


def rp_grn_target_genes(predictions, gene_loc_set):
    bed_file = pybedtools.BedTool.from_dataframe(predictions)
    df = extract_features(bed_file.fn, gene_loc_set).to_frame()
    df["gene"] = df.index.str.split(":").str[1]
    df = df.groupby("gene").mean()
    target_genes = df.to_dict()[0]
    return target_genes


find_targets(
    snakemake.input.predictions,
    snakemake.input.query_file,
    snakemake.output[0],
    snakemake.params,
    snakemake.wildcards,
    m2f_path=getattr(snakemake.input, "m2f_path", None),
)
