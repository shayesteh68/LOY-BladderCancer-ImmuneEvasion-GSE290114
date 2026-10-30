# ============================================================================
# Pipeline: Differential expression analysis for GSE290114
# Contrast: Y_KO versus Y_Scr in MB49 mouse bladder cancer cells
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

cat("[1/5] Loading raw counts and sample metadata...\n")
raw_data <- read.delim(
  "data/raw/GSE290114_gene_count.txt.gz",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata <- read.csv("metadata/samplesheet.csv", stringsAsFactors = FALSE)
required_metadata <- c("sample_id", "condition")
if (!all(required_metadata %in% names(metadata))) {
  stop("Sample sheet must contain sample_id and condition columns.")
}
if (anyDuplicated(metadata$sample_id) || anyDuplicated(raw_data$gene_id)) {
  stop("Sample IDs and gene IDs must be unique.")
}
if (!all(c("gene_id", "gene_name", "gene_chr", "gene_biotype", "gene_description") %in% names(raw_data))) {
  stop("Raw count matrix is missing required gene annotation columns.")
}
count_cols <- metadata$sample_id
if (!all(count_cols %in% names(raw_data))) {
  stop("One or more sample IDs are absent from the raw count matrix.")
}
if (!setequal(unique(metadata$condition), c("Y_Scr", "Y_KO"))) {
  stop("Expected exactly the Y_Scr and Y_KO conditions.")
}
metadata$condition <- factor(metadata$condition, levels = c("Y_Scr", "Y_KO"))
replicate_counts <- table(metadata$condition)
if (any(replicate_counts < 2L)) {
  stop("At least two samples per condition are required for this model.")
}

counts_matrix <- as.matrix(raw_data[, count_cols, drop = FALSE])
if (!is.numeric(counts_matrix) || anyNA(counts_matrix) || any(!is.finite(counts_matrix)) ||
    any(counts_matrix < 0) || any(counts_matrix != floor(counts_matrix))) {
  stop("Raw counts must be finite, non-negative integers.")
}
storage.mode(counts_matrix) <- "integer"
rownames(counts_matrix) <- raw_data$gene_id
rownames(metadata) <- metadata$sample_id
gene_annotations <- raw_data[, c("gene_id", "gene_name", "gene_chr", "gene_biotype", "gene_description"), drop = FALSE]

if (!identical(colnames(counts_matrix), metadata$sample_id)) {
  stop("Count columns and sample metadata are not aligned.")
}
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cat("[2/5] Constructing DESeq2 dataset and pre-filtering...\n")
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = metadata,
  design = ~ condition
)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(sprintf("Retained %d of %d genes after the count >= 10 filter.\n", sum(keep), nrow(counts_matrix)))

cat("[3/5] Performing VST and PCA...\n")
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c("Y_Scr" = "#2b5c8f", "Y_KO" = "#d95f02"), drop = FALSE) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA: MB49 Y-KO versus Y-Scr",
    subtitle = "GSE290114 bulk RNA-seq",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Condition"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
ggsave("results/figures/pca_plot.png", plot = pca_plot, width = 7, height = 5, dpi = 300)

cat("[4/5] Running DESeq2 differential analysis...\n")
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Y_KO", "Y_Scr"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
annotation_index <- match(res_df$gene_id, gene_annotations$gene_id)
if (anyNA(annotation_index)) {
  stop("Could not map every tested gene back to its raw annotation row.")
}
annotation_columns <- setdiff(names(gene_annotations), "gene_id")
res_df <- cbind(res_df, gene_annotations[annotation_index, annotation_columns, drop = FALSE])
res_df <- res_df[order(res_df$padj, res_df$gene_id, na.last = TRUE), , drop = FALSE]
write.csv(res_df, "results/tables/deseq2_deg_all.csv", row.names = FALSE)

sig_df <- res_df[
  !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1.0,
  ,
  drop = FALSE
]
write.csv(sig_df, "results/tables/deseq2_deg_significant.csv", row.names = FALSE)
n_up <- sum(sig_df$log2FoldChange > 0)
n_down <- sum(sig_df$log2FoldChange < 0)
cat(sprintf(
  "Significant DEGs (padj < 0.05, |log2FC| >= 1): %d total (%d up, %d down in Y-KO).\n",
  nrow(sig_df), n_up, n_down
))

cat("[5/5] Generating the volcano plot and run summary...\n")
plot_df <- res_df[
  is.finite(res_df$log2FoldChange) & is.finite(res_df$padj) & res_df$padj > 0,
  ,
  drop = FALSE
]
plot_df$Significance <- "Not Significant"
plot_df$Significance[plot_df$padj < 0.05 & plot_df$log2FoldChange >= 1.0] <- "Up-regulated in Y-KO"
plot_df$Significance[plot_df$padj < 0.05 & plot_df$log2FoldChange <= -1.0] <- "Down-regulated in Y-KO"
plot_df$Significance <- factor(
  plot_df$Significance,
  levels = c("Up-regulated in Y-KO", "Down-regulated in Y-KO", "Not Significant")
)
cat(sprintf("Omitted %d rows with non-finite or zero adjusted p-values from the volcano plot.\n", nrow(res_df) - nrow(plot_df)))
volcano_plot <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c(
    "Up-regulated in Y-KO" = "#d95f02",
    "Down-regulated in Y-KO" = "#2b5c8f",
    "Not Significant" = "gray70"
  ), drop = FALSE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.6) +
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano plot: MB49 Y-KO versus Y-Scr",
    subtitle = "BH adjusted p-value < 0.05 and |log2FC| >= 1",
    x = "Log2 fold change",
    y = "-Log10 adjusted p-value",
    color = "Significance"
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
ggsave("results/figures/volcano_plot.png", plot = volcano_plot, width = 7, height = 6, dpi = 300)

summary_lines <- c(
  "GSE290114 DESeq2 run summary",
  sprintf("Samples: %d (%s)", nrow(metadata), paste(names(replicate_counts), as.integer(replicate_counts), sep = "=", collapse = "; ")),
  "Design: ~ condition; contrast: Y_KO versus Y_Scr",
  "Count filter: row sum >= 10 across samples",
  "DEG thresholds: BH adjusted p-value < 0.05 and absolute log2 fold change >= 1",
  sprintf("Genes tested: %d", nrow(res_df)),
  sprintf("Significant DEGs: %d (up in Y_KO: %d; down in Y_KO: %d)", nrow(sig_df), n_up, n_down),
  sprintf("R: %s", R.version.string),
  sprintf("DESeq2: %s", as.character(packageVersion("DESeq2"))),
  sprintf("ggplot2: %s", as.character(packageVersion("ggplot2")))
)
writeLines(summary_lines, "results/tables/summary.txt")
cat("Done. Figures and tables are saved under results/.\n")
