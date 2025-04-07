library(DESeq2)
library(tidyverse)

expr <- read_csv("data/GSE128242_counts.csv")
mat <- expr %>%
  column_to_rownames("gene_id") %>%
  select(where(is.numeric))

coldata <- tibble(
  sample_id = colnames(mat),
  condition = ifelse(grepl("AP-1-KO", sample_id), "AP1-KO", "WT")
) %>%
  column_to_rownames("sample_id")

dds <- DESeqDataSetFromMatrix(
  countData = mat,
  colData = coldata,
  design = ~condition
)

dds$condition <- factor(dds$condition, levels = c("WT", "AP1-KO"))

dds <- DESeq(dds)
res <- results(dds)
res

df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  left_join(
    expr %>%
      select(gene_id, gene_name),
    by = "gene_id"
  ) %>%
  relocate(gene_name, .after = gene_id)

write_csv(df, "data/GSE128242_deseq2_report.csv")
