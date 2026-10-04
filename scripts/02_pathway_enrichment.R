# ==============================================================================
# Pipeline: KEGG Pathway Enrichment Analysis for GSE290114 (Y-KO vs Y-Scr)
# Author: Narges Shayesteh
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ggplot2)
})

cat("[1/4] Loading significant DEGs and gene mapping...\n")
sig_df <- read.csv("results/tables/deseq2_deg_significant.csv", stringsAsFactors = FALSE)
gene_map <- read.delim("data/reference/mouse_gene2entrez.tsv", stringsAsFactors = FALSE)

# Merge DEGs with Entrez GeneID
sig_mapped <- merge(sig_df, gene_map, by.x = "gene_name", by.y = "Symbol")

down_entrez <- as.character(unique(sig_mapped$GeneID[sig_mapped$log2FoldChange <= -1.0]))
up_entrez   <- as.character(unique(sig_mapped$GeneID[sig_mapped$log2FoldChange >= 1.0]))

cat(sprintf("Mapped to Entrez IDs: %d Down-regulated, %d Up-regulated\n", 
            length(down_entrez), length(up_entrez)))

# --- KEGG Enrichment for Down-regulated genes ---
cat("[2/4] Running KEGG Pathway Enrichment for Down-regulated genes...\n")
kegg_down <- enrichKEGG(
  gene         = down_entrez,
  organism     = "mmu",
  pvalueCutoff = 0.05,
  pAdjustMethod= "BH",
  qvalueCutoff = 0.05
)

if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
  kegg_down_df <- as.data.frame(kegg_down)
  write.csv(kegg_down_df, "results/tables/kegg_enrichment_downregulated.csv", row.names = FALSE)
  
  p_kegg_down <- dotplot(kegg_down, showCategory = 15, title = "KEGG Pathways: Down-regulated in Y-KO") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave("results/figures/kegg_enrichment_downregulated.png", plot = p_kegg_down, width = 9, height = 6, dpi = 300)
  cat(sprintf("Found %d significant KEGG pathways for down-regulated genes.\n", nrow(kegg_down_df)))
} else {
  cat("No significant KEGG pathways found for down-regulated genes at FDR < 0.05.\n")
}

# --- KEGG Enrichment for Up-regulated genes ---
cat("[3/4] Running KEGG Pathway Enrichment for Up-regulated genes...\n")
kegg_up <- enrichKEGG(
  gene         = up_entrez,
  organism     = "mmu",
  pvalueCutoff = 0.05,
  pAdjustMethod= "BH",
  qvalueCutoff = 0.05
)

if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
  kegg_up_df <- as.data.frame(kegg_up)
  write.csv(kegg_up_df, "results/tables/kegg_enrichment_upregulated.csv", row.names = FALSE)
  
  p_kegg_up <- dotplot(kegg_up, showCategory = 15, title = "KEGG Pathways: Up-regulated in Y-KO") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave("results/figures/kegg_enrichment_upregulated.png", plot = p_kegg_up, width = 9, height = 6, dpi = 300)
  cat(sprintf("Found %d significant KEGG pathways for up-regulated genes.\n", nrow(kegg_up_df)))
} else {
  cat("No significant KEGG pathways found for up-regulated genes at FDR < 0.05.\n")
}

cat("[4/4] Completed successfully.\n")
