# ==============================================================================
# Pipeline: Differential Expression Analysis for GSE290114 (MB49 Y-KO vs Y-Scr)
# Author: Narges Shayesteh
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

cat("[1/5] Loading raw count data and sample metadata...\n")
raw_data <- read.delim("data/raw/GSE290114_gene_count.txt.gz", sep = "\t", stringsAsFactors = FALSE)
metadata <- read.csv("metadata/samplesheet.csv", stringsAsFactors = FALSE)

# Extract counts and annotations
count_cols <- metadata$sample_id
counts_matrix <- as.matrix(raw_data[, count_cols])
rownames(counts_matrix) <- raw_data$gene_id

gene_annotations <- raw_data[, c("gene_id", "gene_name", "gene_chr", "gene_biotype", "gene_description")]

# Ensure sample alignment
stopifnot(all(colnames(counts_matrix) == metadata$sample_id))
rownames(metadata) <- metadata$sample_id

cat("[2/5] Constructing DESeq2 dataset and pre-filtering...\n")
metadata$condition <- factor(metadata$condition, levels = c("Y_Scr", "Y_KO"))

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = metadata,
  design = ~ condition
)

# Filter low counts (keep genes with at least 10 counts across all samples)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(sprintf("Retained %d genes after filtering low counts.\n", sum(keep)))

cat("[3/5] Performing Variance Stabilizing Transformation (VST) & PCA...\n")
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c("Y_Scr" = "#2b5c8f", "Y_KO" = "#d95f02")) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA: MB49 Bladder Cancer (Y-KO vs Y-Scr)",
    subtitle = "GSE290114 (Bulk RNA-seq)",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Condition"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("results/figures/pca_plot.png", plot = pca_plot, width = 7, height = 5, dpi = 300)

cat("[4/5] Running DESeq2 differential analysis...\n")
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Y_KO", "Y_Scr"))

# Merge with annotations
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- merge(res_df, gene_annotations, by = "gene_id", all.x = TRUE)
res_df <- res_df[order(res_df$padj), ]

# Export all results
write.csv(res_df, "results/tables/deseq2_deg_all.csv", row.names = FALSE)

# Significant DEGs (padj < 0.05 & |log2FC| >= 1.0)
sig_df <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1.0, ]

write.csv(sig_df, "results/tables/deseq2_deg_significant.csv", row.names = FALSE)

n_up <- sum(sig_df$log2FoldChange > 0)
n_down <- sum(sig_df$log2FoldChange < 0)
cat(sprintf("Significant DEGs (padj < 0.05, |log2FC| >= 1.0): Total = %d (Up in Y-KO: %d, Down in Y-KO: %d)\n",
            nrow(sig_df), n_up, n_down))

cat("[5/5] Generating Volcano Plot...\n")
res_df$Significance <- "Not Significant"
res_df$Significance[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >= 1.0] <- "Up-regulated in Y-KO"
res_df$Significance[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange <= -1.0] <- "Down-regulated in Y-KO"
res_df$Significance <- factor(res_df$Significance, levels = c("Up-regulated in Y-KO", "Down-regulated in Y-KO", "Not Significant"))

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c(
    "Up-regulated in Y-KO" = "#d95f02",
    "Down-regulated in Y-KO" = "#2b5c8f",
    "Not Significant" = "gray70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.6) +
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot: MB49 Y-KO vs Y-Scr",
    subtitle = "FDR < 0.05, |log2FC| >= 1.0",
    x = "log2 Fold Change",
    y = "-log10 adjusted P-value"
  ) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("results/figures/volcano_plot.png", plot = volcano_plot, width = 7, height = 6, dpi = 300)

cat("Done! Figures saved in results/figures/ and tables in results/tables/\n")
