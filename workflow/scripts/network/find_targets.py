import pandas as pd
import dask.dataframe as dd
import json
import pybedtools
from tfsage.features import load_region_set, extract_features
from tfsage.utils import load_features_bed, closest_features
from snakemake.script import snakemake


def find_targets(
    input_file: str,
    output_file: str,
    params: dict,
    wildcards: dict,
    genome: str = "hg38",
    score_threshold: float = 0,
    score_threshold_q: float | None = None,
    m2f_path: str | None = None,
    query_file: str | None = None,
):
    method_class = params.get("method_class")
    method_name = params.get("method_name")
    threshold = wildcards.get("threshold")
    query_id = wildcards.get("query_id")
    factor = wildcards.get("factor")
    linkage_name = wildcards.get("linkage_name")

    genome = params.get("genome", genome)
    score_threshold = params["linkage_config"].get("score_threshold", score_threshold)
    score_threshold_q = params["linkage_config"].get(
        "score_threshold_q", score_threshold_q
    )

    predictions = load_predictions(input_file, method_class, factor, m2f_path)
    if predictions is not None:
        df_predictions = filter_predictions(
            predictions, score_threshold, score_threshold_q
        )

        if df_predictions.empty:
            target_genes = {}
        else:
            if params["linkage_config"].get("intersect_query"):
                df_predictions = intersect_with_query(df_predictions, query_file)

            if params["linkage_config"].get("use_rp"):
                gene_loc_set = load_region_set(genome)
                target_genes = rp_target_genes(df_predictions, gene_loc_set)
            else:
                features_bed = load_features_bed(genome)
                target_genes = binary_target_genes(df_predictions, features_bed)
    else:
        target_genes = {}

    # Combine everything into one dictionary
    results = {
        "method_class": method_class,
        "method_name": method_name,
        "threshold": threshold,
        "query_id": query_id,
        "factor": factor,
        "linkage_name": linkage_name,
        "target_genes": target_genes,
    }

    # Save results to file
    with open(output_file, "w") as f:
        json.dump(results, f, indent=4)


def load_predictions(
    predictions_path,
    method_class,
    factor: str | None = None,
    m2f_path: str | None = None,
) -> pd.DataFrame | None:
    if method_class == "tfsage":
        predictions = pd.read_parquet(predictions_path).assign(score=lambda x: x["sum"])
    elif method_class == "motif scan":
        predictions = load_predictions_motif_scan(predictions_path, factor, m2f_path)
    else:
        raise ValueError(f"Unknown method class: {method_class}")
    return predictions


def load_predictions_motif_scan(
    predictions_path, factor, m2f_path
) -> pd.DataFrame | None:
    factor = factor.split("-")[0]
    motif_list = (
        pd.read_csv(m2f_path, sep="\t").query("Factor == @factor")["Motif"].tolist()
    )
    if len(motif_list) == 0:
        return None

    ddf = dd.read_csv(predictions_path, sep="\t", header=None, comment="#")
    ddf = ddf[ddf[3].isin(motif_list)]
    predictions = ddf.compute()
    predictions.columns = ["chrom", "start", "end", "motif", "score", "strand"]
    return predictions


def filter_predictions(predictions, score_threshold, score_threshold_q):
    if score_threshold_q is not None:
        score_threshold = predictions["score"].quantile(score_threshold_q)

    df_predictions = predictions.query("score >= @score_threshold")
    return df_predictions


def intersect_with_query(df_predictions, query_file):
    bed_file = pybedtools.BedTool.from_dataframe(df_predictions)
    if query_file is not None:
        query_bed = pybedtools.BedTool(query_file)
        bed_file = bed_file.intersect(query_bed, u=True)

    df_predictions = bed_file.to_dataframe(
        disable_auto_names=True, names=df_predictions.columns
    )
    return df_predictions


def rp_target_genes(df_predictions, gene_loc_set):
    bed_file = pybedtools.BedTool.from_dataframe(df_predictions)
    df = extract_features(bed_file.fn, gene_loc_set).to_frame()
    df["gene"] = df.index.str.split(":").str[1]
    df = df.groupby("gene").mean()
    target_genes = df.to_dict()[0]
    return target_genes


def binary_target_genes(df_predictions, features_bed):
    bed_file = pybedtools.BedTool.from_dataframe(df_predictions).sort()
    df = closest_features(bed_file, features_bed)
    df = df.query("-1000 < distance < 100")
    df["gene"] = df["name"].str.split(":").str[1]
    target_genes = {gene: 1.0 for gene in df["gene"].unique()}
    return target_genes


find_targets(
    snakemake.input["predictions"],
    snakemake.output[0],
    snakemake.params,
    snakemake.wildcards,
    m2f_path=getattr(snakemake.input, "m2f_path", None),
    query_file=getattr(snakemake.input, "query_file", None),
)
