import pandas as pd
from tqdm import tqdm
from snakemake.script import snakemake


def filter_motifs(m2f_infile, pfm_infile, m2f_outfile, pfm_outfile, tf_list):
    # Filter motif2factors
    m2f_filtered = pd.read_csv(m2f_infile, sep="\t").query("Factor in @tf_list")
    motif_list = m2f_filtered["Motif"].unique()

    # Filter PFM file
    pfm_filtered = filter_pfm(pfm_infile, motif_list)

    # Save filtered motif2factors and PFM
    m2f_filtered.to_csv(m2f_outfile, sep="\t", index=False)
    with open(pfm_outfile, "w") as f:
        f.writelines(pfm_filtered)


def filter_pfm(pfm_file, motif_list):
    with open(pfm_file, "r") as f:
        pfm_lines = f.readlines()

    include = False
    lines = []
    for line in tqdm(pfm_lines):
        if line.startswith(">"):
            motif = line.replace(">", "").strip()
            include = motif in motif_list
        if include or line.startswith("#"):
            lines.append(line)

    return lines


filter_motifs(
    snakemake.input.motif2factors,
    snakemake.input.pfm,
    snakemake.output.motif2factors,
    snakemake.output.pfm,
    snakemake.params.tf_list,
)
