# ============================================================================
# Pipeline: TRRUST over-representation analysis for GSE290114
# ============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ggplot2)
})

cat("[1/4] Loading TRRUST interactions and differential-expression results...\n")
trrust <- read.delim(
  "data/annotation/trrust_rawdata.mouse.tsv",
  header = FALSE,
  stringsAsFactors = FALSE
)
if (ncol(trrust) < 4L) {
  stop("TRRUST input must have TF, target, mode, and PMID columns.")
}
colnames(trrust)[1:4] <- c("TF", "Target", "Mode", "PMID")
term2gene <- unique(trrust[, c("TF", "Target")])
term2gene$TF <- trimws(term2gene$TF)
term2gene$Target <- trimws(term2gene$Target)
term2gene <- term2gene[
  !is.na(term2gene$TF) & nzchar(term2gene$TF) &
    !is.na(term2gene$Target) & nzchar(term2gene$Target),
  ,
  drop = FALSE
]

deg_sig <- read.csv("results/tables/deseq2_deg_significant.csv", stringsAsFactors = FALSE)
deg_all <- read.csv("results/tables/deseq2_deg_all.csv", stringsAsFactors = FALSE)
if (!all(c("gene_name", "log2FoldChange", "padj") %in% names(deg_sig)) ||
    !"gene_name" %in% names(deg_all)) {
  stop("DE result files are missing the required gene and statistic columns.")
}
down_candidates <- unique(na.omit(deg_sig$gene_name[
  deg_sig$log2FoldChange <= -1.0 & deg_sig$padj < 0.05
]))
universe_candidates <- unique(na.omit(deg_all$gene_name))
# The enrichment background is the tested genes represented in TRRUST.
universe_genes <- intersect(universe_candidates, unique(term2gene$Target))
down_genes <- intersect(down_candidates, universe_genes)
cat(sprintf(
  "-> Downregulated target genes represented in TRRUST: %d\n-> Tested genes represented in TRRUST: %d\n",
  length(down_genes), length(universe_genes)
))

if (length(down_genes) < 5L || length(universe_genes) < 5L) {
  stop("Too few TRRUST-mapped genes to run the enrichment test.")
}

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
result_path <- "results/tables/trrust_tf_enrichment_downregulated.csv"
summary_path <- "results/tables/trrust_tf_enrichment_summary.csv"
figure_path <- "results/figures/tf_enrichment_downregulated.png"

cat("[2/4] Performing over-representation analysis...\n")
tf_enrich <- enricher(
  gene = down_genes,
  universe = universe_genes,
  TERM2GENE = term2gene,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 5,
  maxGSSize = 500
)
tf_res <- as.data.frame(tf_enrich)
if (is.null(tf_enrich) || nrow(tf_res) == 0L) {
  unlink(c(result_path, figure_path), force = TRUE)
  summary <- data.frame(
    downregulated_trrust_targets = length(down_genes),
    tested_trrust_targets = length(universe_genes),
    significant_transcription_factors = 0L,
    status = "no_significant_factors",
    analysis_date = as.character(Sys.Date()),
    stringsAsFactors = FALSE
  )
  write.csv(summary, summary_path, row.names = FALSE)
  cat("No TRRUST factors passed the adjusted p-value threshold.\n")
  quit(save = "no", status = 0)
}

write.csv(tf_res, result_path, row.names = FALSE)
summary <- data.frame(
  downregulated_trrust_targets = length(down_genes),
  tested_trrust_targets = length(universe_genes),
  significant_transcription_factors = nrow(tf_res),
  status = "significant_factors_found",
  analysis_date = as.character(Sys.Date()),
  stringsAsFactors = FALSE
)
write.csv(summary, summary_path, row.names = FALSE)
cat(sprintf("-> Found %d significantly enriched TRRUST factors.\n", nrow(tf_res)))

cat("[3/4] Plotting the numeric gene ratios for the top factors...\n")
top_n <- min(15L, nrow(tf_res))
plot_df <- head(tf_res, top_n)
ratio_parts <- strsplit(as.character(plot_df$GeneRatio), "/", fixed = TRUE)
plot_df$GeneRatioNumeric <- vapply(
  ratio_parts,
  function(x) {
    if (length(x) != 2L) return(NA_real_)
    numerator <- suppressWarnings(as.numeric(x[1]))
    denominator <- suppressWarnings(as.numeric(x[2]))
    if (!is.finite(numerator) || !is.finite(denominator) || denominator <= 0) {
      return(NA_real_)
    }
    numerator / denominator
  },
  numeric(1)
)
if (anyNA(plot_df$GeneRatioNumeric)) {
  stop("Could not convert one or more TRRUST GeneRatio values to numeric ratios.")
}
plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))

p <- ggplot(plot_df, aes(x = GeneRatioNumeric, y = Description, size = Count, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "#E41A1C", high = "#377EB8") +
  theme_bw(base_size = 12) +
  labs(
    title = "TRRUST factors linked to genes downregulated in Y-KO",
    subtitle = "Candidate regulatory associations; hypergeometric ORA with BH adjustment",
    x = "Gene ratio among annotated downregulated targets",
    y = "Transcription factor",
    color = "Adjusted p-value",
    size = "Target count"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 10)
  )
ggsave(figure_path, plot = p, width = 8, height = 6, dpi = 300)
cat("[4/4] TRRUST table, summary, and figure saved.\n")
