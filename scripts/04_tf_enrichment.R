# ==============================================================================
# Pipeline: Upstream Transcription Factor Enrichment Analysis (TRRUST v2)
# Dataset: GSE290114 (MB49 Y-KO vs Y-Scr)
# Author: Narges Shayesteh
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ggplot2)
})

cat("[1/4] Loading TRRUST transcription factor database...\n")
trrust <- read.table("data/annotation/trrust_rawdata.mouse.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(trrust) <- c("TF", "Target", "Mode", "PMID")

# Prepare TERM2GENE data frame for clusterProfiler: (Term = TF, Gene = Target)
term2gene <- trrust[, c("TF", "Target")]

cat("[2/4] Loading differential expression results...\n")
deg_sig <- read.csv("results/tables/deseq2_deg_significant.csv", stringsAsFactors = FALSE)
deg_all <- read.csv("results/tables/deseq2_deg_all.csv", stringsAsFactors = FALSE)

# Extract down-regulated genes (repressed by LOY)
down_genes <- deg_sig$gene_name[deg_sig$log2FoldChange <= -1.0 & deg_sig$padj < 0.05]
universe_genes <- unique(deg_all$gene_name)

cat(sprintf("-> Target down-regulated genes: %d\n", length(down_genes)))
cat(sprintf("-> Tested universe genes: %d\n", length(universe_genes)))

cat("[3/4] Performing Over-Representation Analysis for TFs...\n")
tf_enrich <- enricher(
  gene = down_genes,
  universe = universe_genes,
  TERM2GENE = term2gene,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 5,
  maxGSSize = 500
)

if (is.null(tf_enrich) || nrow(as.data.frame(tf_enrich)) == 0) {
  stop("No significant TFs enriched at p.adjust < 0.05. Check gene symbol consistency.")
}

tf_res <- as.data.frame(tf_enrich)
cat(sprintf("-> Found %d significantly enriched transcription factors.\n", nrow(tf_res)))

# Save full results table
write.csv(tf_res, "results/tables/trrust_tf_enrichment_downregulated.csv", row.names = FALSE)
cat("-> Saved results to results/tables/trrust_tf_enrichment_downregulated.csv\n")

cat("[4/4] Generating Dotplot for top enriched transcription factors...\n")
top_n <- min(15, nrow(tf_res))
plot_df <- head(tf_res, top_n)
plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))

p <- ggplot(plot_df, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "#E41A1C", high = "#377EB8") +
  theme_bw(base_size = 12) +
  labs(
    title = "Top Upstream Transcription Factors Repressed in Y-KO (GSE290114)",
    subtitle = "TRRUST v2 Database | Hypergeometric Test (BH adjusted)",
    x = "Gene Ratio",
    y = "Transcription Factor",
    color = "p.adjust",
    size = "Target Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 10)
  )

ggsave("results/figures/tf_enrichment_downregulated.png", plot = p, width = 8, height = 6, dpi = 300)
cat("-> Figure saved to results/figures/tf_enrichment_downregulated.png\n")
