# ==============================================================================
# Pipeline: Targeted Heatmap of LOY and Immune Evasion Genes (GSE290114)
# Author: Narges Shayesteh
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
})

cat("[1/4] Loading raw counts, metadata, and DEG results...\n")
counts_raw <- read.table("data/raw/GSE290114_gene_count.txt.gz", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
samplesheet <- read.csv("metadata/samplesheet.csv", stringsAsFactors = FALSE)

# Explicitly select only the sample columns for count matrix
rownames(counts_raw) <- counts_raw$gene_id
count_cols <- samplesheet$sample_id
count_mat <- as.matrix(counts_raw[, count_cols])
mode(count_mat) <- "integer"

rownames(samplesheet) <- samplesheet$sample_id

cat("[2/4] Normalizing counts with DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = samplesheet, design = ~ condition)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

# Select key genes (Y-chromosome and Immune response)
key_genes <- c(
  "Ddx3y", "Eif2s3y", "Kdm5d", "Uty",
  "Cxcl1", "Cxcl2", "Cxcl10", "Ccl2", "Ccl7",
  "Csf2", "Tnf", "Txnip"
)

gene_annot <- counts_raw[, c("gene_id", "gene_name")]
matched_genes <- gene_annot[gene_annot$gene_name %in% key_genes, ]
matched_genes <- matched_genes[!duplicated(matched_genes$gene_name), ]

sub_counts <- norm_counts[matched_genes$gene_id, ]
rownames(sub_counts) <- matched_genes$gene_name

# Add log2 transformation (log2(normalized_counts + 1))
sub_log2 <- log2(sub_counts + 1)

# Sort genes for clear biological grouping
y_genes <- intersect(c("Ddx3y", "Eif2s3y", "Kdm5d", "Uty"), rownames(sub_log2))
immune_genes <- setdiff(rownames(sub_log2), y_genes)
sub_log2 <- sub_log2[c(y_genes, immune_genes), ]

# Column annotation
anno_col <- data.frame(
  Condition = factor(samplesheet$condition, levels = c("Y_Scr", "Y_KO")),
  row.names = samplesheet$sample_id
)
anno_colors <- list(
  Condition = c(Y_Scr = "#2B83BA", Y_KO = "#D7191C")
)

cat(sprintf("[3/4] Plotting heatmap for %d target genes...\n", nrow(sub_log2)))
png("results/figures/key_genes_heatmap.png", width = 2100, height = 2400, res = 300)
pheatmap(
  sub_log2,
  scale = "row",
  annotation_col = anno_col,
  annotation_colors = anno_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  fontsize = 11,
  fontsize_row = 11,
  fontsize_col = 11,
  angle_col = 0,
  main = "Expression of Y-Chromosome & Repressed Immune Genes (GSE290114)"
)
dev.off()

cat("[4/4] Heatmap saved to results/figures/key_genes_heatmap.png\n")
