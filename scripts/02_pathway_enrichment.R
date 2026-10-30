# ============================================================================
# Pipeline: KEGG over-representation analysis for GSE290114
# ============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ggplot2)
})

cat("[1/4] Loading differential-expression results and gene mapping...\n")
sig_df <- read.csv("results/tables/deseq2_deg_significant.csv", stringsAsFactors = FALSE)
deg_all <- read.csv("results/tables/deseq2_deg_all.csv", stringsAsFactors = FALSE)
gene_map <- read.delim("data/reference/mouse_gene2entrez.tsv", stringsAsFactors = FALSE)
if (!all(c("gene_name", "log2FoldChange", "padj") %in% names(sig_df)) ||
    !"gene_name" %in% names(deg_all) ||
    !all(c("Symbol", "GeneID") %in% names(gene_map))) {
  stop("Input tables are missing required gene, fold-change, or mapping columns.")
}

gene_map <- unique(gene_map[, c("Symbol", "GeneID")])
all_symbols <- unique(na.omit(deg_all$gene_name))
universe_map <- merge(
  data.frame(gene_name = all_symbols, stringsAsFactors = FALSE),
  gene_map,
  by.x = "gene_name",
  by.y = "Symbol"
)
universe_entrez <- unique(as.character(universe_map$GeneID[!is.na(universe_map$GeneID)]))
if (length(universe_entrez) < 100L) {
  stop("Too few tested genes map to Entrez IDs to define an enrichment background.")
}

sig_mapped <- merge(
  sig_df,
  gene_map,
  by.x = "gene_name",
  by.y = "Symbol"
)
down_entrez <- unique(as.character(sig_mapped$GeneID[
  !is.na(sig_mapped$GeneID) & sig_mapped$log2FoldChange <= -1.0
]))
up_entrez <- unique(as.character(sig_mapped$GeneID[
  !is.na(sig_mapped$GeneID) & sig_mapped$log2FoldChange >= 1.0
]))
cat(sprintf(
  "Mapped significant genes to Entrez IDs: %d down, %d up; tested background: %d IDs.\n",
  length(down_entrez), length(up_entrez), length(universe_entrez)
))

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

run_kegg <- function(gene_ids, direction) {
  table_path <- sprintf("results/tables/kegg_enrichment_%s.csv", direction)
  figure_path <- sprintf("results/figures/kegg_enrichment_%s.png", direction)
  if (!length(gene_ids)) {
    unlink(c(table_path, figure_path), force = TRUE)
    return(data.frame(
      direction = direction,
      mapped_gene_ids = 0L,
      background_gene_ids = length(universe_entrez),
      significant_pathways = 0L,
      status = "no_mapped_genes",
      stringsAsFactors = FALSE
    ))
  }

  result <- enrichKEGG(
    gene = gene_ids,
    universe = universe_entrez,
    organism = "mmu",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  result_df <- as.data.frame(result)
  if (!is.null(result) && nrow(result_df) > 0L) {
    write.csv(result_df, table_path, row.names = FALSE)
    plot_df <- head(result_df, 15L)
    plot_df$GeneRatioNumeric <- vapply(
      strsplit(as.character(plot_df$GeneRatio), "/", fixed = TRUE),
      function(x) as.numeric(x[1]) / as.numeric(x[2]),
      numeric(1)
    )
    plot_df$Description <- sub(" - Mus musculus.*$", "", plot_df$Description)
    plot_df$Description <- vapply(
      strsplit(plot_df$Description, " ", fixed = TRUE),
      function(words) paste(strwrap(paste(words, collapse = " "), width = 38), collapse = "\n"),
      character(1)
    )
    plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
    direction_label <- if (direction == "downregulated") "downregulated" else "upregulated"
    plot <- ggplot(plot_df, aes(x = GeneRatioNumeric, y = Description, size = Count, color = p.adjust)) +
      geom_point() +
      scale_color_gradient(low = "#D95F02", high = "#2166AC") +
      scale_size(range = c(3, 10)) +
      theme_bw(base_size = 11) +
      labs(
        title = sprintf("KEGG pathways: genes %s in Y-KO", direction_label),
        x = "Gene ratio",
        y = NULL,
        color = "Adjusted p-value",
        size = "Gene count"
      ) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 9)
      )
    ggsave(figure_path, plot = plot, width = 10, height = 8, dpi = 300)
    status <- "significant_pathways_found"
  } else {
    unlink(c(table_path, figure_path), force = TRUE)
    status <- "no_significant_pathways"
  }
  data.frame(
    direction = direction,
    mapped_gene_ids = length(gene_ids),
    background_gene_ids = length(universe_entrez),
    significant_pathways = nrow(result_df),
    status = status,
    stringsAsFactors = FALSE
  )
}

cat("[2/4] Testing pathways for genes down in Y-KO...\n")
down_summary <- run_kegg(down_entrez, "downregulated")
cat("[3/4] Testing pathways for genes up in Y-KO...\n")
up_summary <- run_kegg(up_entrez, "upregulated")
summary <- rbind(down_summary, up_summary)
summary$analysis_date <- as.character(Sys.Date())
summary$annotation_source <- "KEGG REST API (mouse/mmu)"
write.csv(summary, "results/tables/kegg_enrichment_summary.csv", row.names = FALSE)
print(summary[, c("direction", "significant_pathways", "status")], row.names = FALSE)
cat("[4/4] KEGG analysis completed. Annotation is retrieved online at run time.\n")
