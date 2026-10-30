# ============================================================================
# Pipeline: Targeted expression heatmap for GSE290114
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
})

cat("[1/4] Loading raw counts and sample metadata...\n")
counts_raw <- read.delim(
  "data/raw/GSE290114_gene_count.txt.gz",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
samplesheet <- read.csv("metadata/samplesheet.csv", stringsAsFactors = FALSE)
required_metadata <- c("sample_id", "condition")
if (!all(required_metadata %in% names(samplesheet))) {
  stop("Sample sheet must contain sample_id and condition columns.")
}
if (anyDuplicated(samplesheet$sample_id) || anyDuplicated(counts_raw$gene_id)) {
  stop("Sample IDs and gene IDs must be unique.")
}
count_cols <- samplesheet$sample_id
if (!all(count_cols %in% names(counts_raw))) {
  stop("One or more sample IDs are absent from the raw count matrix.")
}
if (!all(c("gene_id", "gene_name") %in% names(counts_raw))) {
  stop("Raw count matrix must contain gene_id and gene_name columns.")
}
if (!setequal(unique(samplesheet$condition), c("Y_Scr", "Y_KO"))) {
  stop("Expected exactly the Y_Scr and Y_KO conditions.")
}
samplesheet$condition <- factor(samplesheet$condition, levels = c("Y_Scr", "Y_KO"))
if (any(table(samplesheet$condition) < 2L)) {
  stop("At least two biological samples per condition are required.")
}

count_mat <- as.matrix(counts_raw[, count_cols, drop = FALSE])
if (!is.numeric(count_mat) || anyNA(count_mat) || any(!is.finite(count_mat)) ||
    any(count_mat < 0) || any(count_mat != floor(count_mat))) {
  stop("Raw counts must be finite, non-negative integers.")
}
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- counts_raw$gene_id
rownames(samplesheet) <- samplesheet$sample_id

cat("[2/4] Normalizing counts with DESeq2...\n")
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = samplesheet,
  design = ~ condition
)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

key_genes <- c(
  "Ddx3y", "Eif2s3y", "Kdm5d", "Uty",
  "Cxcl1", "Cxcl2", "Cxcl10", "Ccl2", "Ccl7", "Csf2", "Tnf", "Txnip"
)
gene_annot <- counts_raw[, c("gene_id", "gene_name"), drop = FALSE]
matched_genes <- gene_annot[
  !is.na(gene_annot$gene_name) & gene_annot$gene_name %in% key_genes,
  ,
  drop = FALSE
]
matched_genes <- matched_genes[!duplicated(matched_genes$gene_name), , drop = FALSE]

coverage <- data.frame(
  target_gene = key_genes,
  matched = key_genes %in% matched_genes$gene_name,
  gene_id = matched_genes$gene_id[match(key_genes, matched_genes$gene_name)],
  status = ifelse(key_genes %in% matched_genes$gene_name, "matched", "not_found"),
  stringsAsFactors = FALSE
)
if (any(!coverage$matched)) {
  warning("Some requested target genes were not found; see the coverage table.")
}

sub_counts <- norm_counts[matched_genes$gene_id, , drop = FALSE]
rownames(sub_counts) <- matched_genes$gene_name
sub_log2 <- log2(sub_counts + 1)
row_sd <- apply(sub_log2, 1, sd)
variable_rows <- is.finite(row_sd) & row_sd > 0
coverage$status[coverage$matched & !coverage$target_gene %in% rownames(sub_log2)[variable_rows]] <- "zero_variance_not_plotted"
coverage$status[coverage$target_gene %in% rownames(sub_log2)[variable_rows]] <- "plotted"
sub_log2 <- sub_log2[variable_rows, , drop = FALSE]
if (!nrow(sub_log2)) {
  stop("No variable target genes are available for the heatmap.")
}

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
write.csv(coverage, "results/tables/heatmap_target_gene_coverage.csv", row.names = FALSE)

y_genes <- intersect(c("Ddx3y", "Eif2s3y", "Kdm5d", "Uty"), rownames(sub_log2))
immune_genes <- setdiff(rownames(sub_log2), y_genes)
sub_log2 <- sub_log2[c(y_genes, immune_genes), , drop = FALSE]

anno_col <- data.frame(
  Condition = samplesheet$condition,
  row.names = samplesheet$sample_id
)
anno_colors <- list(Condition = c(Y_Scr = "#2B83BA", Y_KO = "#D7191C"))

cat(sprintf("[3/4] Plotting heatmap for %d variable target genes...\n", nrow(sub_log2)))
png("results/figures/key_genes_heatmap.png", width = 2700, height = 2400, res = 300)
pheatmap(
  sub_log2,
  scale = "row",
  annotation_col = anno_col,
  annotation_colors = anno_colors,
  annotation_names_col = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 0,
  main = "GSE290114: Y-chromosome and immune targets"
)
dev.off()
cat("[4/4] Heatmap and target coverage table saved.\n")
