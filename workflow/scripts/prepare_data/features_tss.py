from pathlib import Path
from tfsage.features import extract_features_parallel, load_region_set
from snakemake.script import snakemake


def extract_features(bed_files, output_file, max_workers=None):
    # Load reference gene locations
    gene_loc_set = load_region_set("hg38")

    # Extract features
    df = extract_features_parallel(bed_files, gene_loc_set, max_workers=max_workers)
    df.columns = [f"{Path(f).stem}" for f in bed_files]

    # Save results
    df.to_parquet(output_file)


extract_features(snakemake.input, snakemake.output[0], snakemake.threads)
