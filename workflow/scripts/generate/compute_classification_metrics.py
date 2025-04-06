import pandas as pd
import dask.dataframe as dd
import json
from tfsage import utils
from snakemake.script import snakemake


def compute_classification_metrics(
    predictions_path,
    test_set_path,
    output_file,
    method_class,
    method_name,
    query_id,
    target_id,
    factor,
    m2f_path,
    test_set_name,
    query_assay,
    cell_type,
    threshold,
):
    predictions = load_predictions(predictions_path, method_class, factor, m2f_path)
    test_set = pd.read_parquet(test_set_path)

    n_samples = len(test_set)
    p_positive = test_set["positive"].mean()

    if predictions is not None:
        # Intersect predictions with test set
        df = utils.intersect_predictions_with_test_set(predictions, test_set)
        df["score"] = utils.sanitize_scores(df["score"])

        # Get true labels and scores
        y_true = df["positive"]
        y_score = df["score"]
        y_pred = y_score > 0

        # Compute metrics
        metrics = utils.compute_classification_metrics(y_true, y_score, y_pred)
    else:
        metrics = None

    # Combine everything into one dictionary
    results = {
        "method_class": method_class,
        "method_name": method_name,
        "query_id": query_id,
        "query_assay": query_assay,
        "cell_type": cell_type,
        "factor": factor,
        "target_id": target_id,
        "test_set": test_set_name,
        "n_samples": n_samples,
        "p_positive": p_positive,
        "metrics": metrics,
        "threshold": threshold,
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


compute_classification_metrics(
    snakemake.input.predictions,
    snakemake.input.test_set,
    snakemake.output[0],
    snakemake.params.method_class,
    snakemake.params.method_name,
    snakemake.wildcards.query_id,
    snakemake.wildcards.target_id,
    snakemake.wildcards.factor,
    getattr(snakemake.input, "m2f_path", None),
    snakemake.wildcards.test_set_name,
    snakemake.params.query_assay,
    snakemake.params.cell_type,
    snakemake.wildcards.threshold,
)
