import pandas as pd
import dask.dataframe as dd
import json
import pybedtools
from tfsage.features import load_region_set, extract_features
from snakemake.script import snakemake


def compute_rp(
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
    mode = wildcards.get("mode", None)
    threshold = wildcards.get("threshold", None)
    query_id = wildcards.get("query_id", None)
    factor = wildcards.get("factor", None)

    genome = params.get("genome", genome)
    score_threshold = params["network_config"].get("score_threshold", score_threshold)
    score_threshold_q = params["network_config"].get(
        "score_threshold_q", score_threshold_q
    )
    method_class = params.get("method_class", None)
    method_name = params.get("method_name", None)

    gene_loc_set = load_region_set(genome)
    predictions = load_predictions(input_file, method_class, factor, m2f_path)

    if predictions is not None:
        if score_threshold_q is not None:
            score_threshold = predictions["score"].quantile(score_threshold_q)

        df_predictions = predictions.query("score >= @score_threshold")
        if df_predictions.empty:
            target_genes = {}
        else:
            bed_file = pybedtools.BedTool.from_dataframe(df_predictions)
            if query_file is not None:
                # Intersect with query file before extracting features
                query_bed = pybedtools.BedTool(query_file)
                bed_file = bed_file.intersect(query_bed, u=True)

            df = extract_features(bed_file.fn, gene_loc_set).to_frame()
            df["gene"] = df.index.str.split(":").str[1]
            df = df.groupby("gene").mean()
            target_genes = df.to_dict()[0]
    else:
        target_genes = {}

    # Combine everything into one dictionary
    results = {
        "method_class": method_class,
        "method_name": method_name,
        "mode": mode,
        "threshold": threshold,
        "query_id": query_id,
        "factor": factor,
        "score_threshold": score_threshold,
        "score_threshold_q": score_threshold_q,
        "genome": genome,
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


compute_rp(
    snakemake.input["predictions"],
    snakemake.output[0],
    snakemake.params,
    snakemake.wildcards,
    m2f_path=getattr(snakemake.input, "m2f_path", None),
    query_file=getattr(snakemake.input, "query_file", None),
)
