# End-to-End Transcriptomic Dissection of LOY-Driven Immune Evasion in Bladder Cancer

[![Domain: Bioinformatics | Oncogenomics](https://img.shields.io/badge/Domain-Bioinformatics%20%7C%20Oncogenomics-blue.svg)](https://github.com/shayesteh68)
[![Pipeline: DESeq2 | KEGG | TRRUST](https://img.shields.io/badge/Pipeline-DESeq2%20%7C%20KEGG%20%7C%20TRRUST-brightgreen.svg)](https://github.com/shayesteh68)
[![GEO: GSE290114](https://img.shields.io/badge/GEO-GSE290114-orange.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290114)

## 📌 Scientific Context & Objectives

Loss of the Y chromosome (LOY) is one of the most frequent somatic chromosomal alterations in male cancer and has been associated with tumour aggressiveness and with altered anti-tumour immunity, including impaired response to immune checkpoint blockade.

This repository holds the computational pipeline and the committed analysis outputs for a comparative RNA-seq analysis of **CRISPR-mediated whole-Y-chromosome knockout MB49 urothelial (bladder) carcinoma cells (`Y_KO`)** versus **scrambled-control MB49 cells (`Y_Scr`)**. All analyses use the mouse RNA-seq dataset **NCBI GEO accession [GSE290114](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290114)**.

The pipeline is organised as four sequential analyses:

1. **Differential expression analysis (DESeq2)** of `Y_KO` versus `Y_Scr`, including sample-level quality control (VST-transformed PCA) and a volcano plot.
2. **KEGG pathway enrichment** (`clusterProfiler::enrichKEGG`, `organism = "mmu"`) of the differentially expressed gene sets.
3. **Targeted heatmap** of Y-chromosome-linked genes together with immune / chemokine-related genes.
4. **Transcription factor (TF) enrichment** of the down-regulated gene set against the **TRRUST v2** mouse regulatory interaction database.

## 🔬 Dataset Overview

- **Accession:** NCBI GEO [GSE290114](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290114)
- **Reference article:** Chen X. *et al.*, *Concurrent loss of the Y chromosome in cancer and T cells impacts outcome*, ***Nature*** **642**, 1041–1050 (2025). DOI: [10.1038/s41586-025-09071-2](https://doi.org/10.1038/s41586-025-09071-2)
- **Organism:** *Mus musculus*
- **Cell model:** MB49 murine bladder cancer cell line
- **Experimental groups:**
  - `Y_Scr` — scrambled control (Y chromosome intact)
  - `Y_KO` — CRISPR whole-Y-chromosome knockout
- **Tested contrast:** `Y_KO` versus `Y_Scr`
- **Significance thresholds applied in the pipeline:** `padj < 0.05` and `|log2FC| >= 1.0`

## 🗂️ Repository Structure

```text
LOY-BladderCancer-ImmuneEvasion-GSE290114/
├── README.md
├── data/
│   └── annotation/
│       └── trrust_rawdata.mouse.tsv
├── scripts/
│   ├── 01_deseq2_analysis.R        # DESeq2: QC, VST/PCA, differential expression
│   ├── 02_pathway_enrichment.R     # KEGG pathway enrichment of DEG sets
│   ├── 03_plot_heatmap.R           # Targeted heatmap (Y-chromosome + immune genes)
│   └── 04_tf_enrichment.R          # TRRUST v2 transcription factor enrichment
└── results/
    ├── tables/
    │   ├── deseq2_deg_all.csv
    │   ├── deseq2_deg_significant.csv
    │   ├── kegg_enrichment_downregulated.csv
    │   └── trrust_tf_enrichment_downregulated.csv
    └── figures/
        ├── pca_plot.png
        ├── volcano_plot.png
        ├── key_genes_heatmap.png
        ├── kegg_enrichment_downregulated.png
        └── tf_enrichment_downregulated.png
```

## 📥 Required Input Files (not included in this repository)

The scripts read their inputs from paths that are deliberately kept outside version control. Before running the pipeline, place the following files in the expected locations:

| Expected path | Description |
| --- | --- |
| `data/raw/GSE290114_gene_count.txt.gz` | Raw gene-level count matrix for GSE290114 (download from GEO). |
| `metadata/samplesheet.csv` | Sample sheet with the columns `sample_id` and `condition`; the `condition` column must contain the levels `Y_Scr` and `Y_KO`. |
| `data/reference/mouse_gene2entrez.tsv` | Mouse gene symbol to Entrez ID mapping used for KEGG enrichment. |
| `data/annotation/trrust_rawdata.mouse.tsv` | TRRUST v2 mouse regulatory interactions (already included in this repository). |

Raw sequencing data and the count matrix are not redistributed here; they are available from the GEO accession above.

## ⚙️ Environment & Dependencies

- **R:** `>= 4.3.0`
- **Key Bioconductor packages:** `DESeq2` (differential expression), `clusterProfiler` and `org.Mm.eg.db` (KEGG enrichment with the mouse `mmu` annotation), `pheatmap` (heatmap rendering)
- **Key CRAN packages:** `ggplot2` (PCA and volcano plots)

The complete list of loaded packages is declared at the top of each script, and the conda environment is captured in `environment.yml`. A `sessionInfo()` snapshot is not currently committed to this repository.

## 🚀 Pipeline Execution

```bash
# Clone the repository
git clone https://github.com/shayesteh68/LOY-BladderCancer-ImmuneEvasion-GSE290114.git
cd LOY-BladderCancer-ImmuneEvasion-GSE290114

# Step 1: Differential expression analysis (QC, PCA, volcano plot)
Rscript scripts/01_deseq2_analysis.R

# Step 2: KEGG pathway enrichment
Rscript scripts/02_pathway_enrichment.R

# Step 3: Targeted heatmap
Rscript scripts/03_plot_heatmap.R

# Step 4: TRRUST transcription factor enrichment
Rscript scripts/04_tf_enrichment.R
```

All scripts assume they are launched from the repository root directory and write their outputs into `results/tables/` and `results/figures/`.

## 📊 Committed Outputs

**Tables (`results/tables/`)**

- `deseq2_deg_all.csv` — DESeq2 results for all tested genes.
- `deseq2_deg_significant.csv` — genes passing the significance thresholds.
- `kegg_enrichment_downregulated.csv` — KEGG terms enriched in down-regulated genes.
- `trrust_tf_enrichment_downregulated.csv` — TRRUST v2 transcription factors enriched in down-regulated genes.

**Figures (`results/figures/`)**

- `pca_plot.png` — PCA of VST-transformed samples.
- `volcano_plot.png` — volcano plot of the `Y_KO` versus `Y_Scr` contrast.
- `key_genes_heatmap.png` — targeted heatmap of Y-chromosome and immune-related genes.
- `kegg_enrichment_downregulated.png` — KEGG enrichment bar/dot plot (down-regulated genes).
- `tf_enrichment_downregulated.png` — TRRUST TF enrichment plot (down-regulated genes).

## ⚠️ Scope and Limitations

- Enrichment results are reported for the **down-regulated** gene set only. Up-regulated KEGG enrichment outputs (the corresponding table and figure) are **not** part of this repository in its current state.
- The raw count matrix and the sample sheet are not versioned here and must be obtained from GEO (see *Required Input Files*).
- This repository focuses on bulk RNA-seq differential expression, pathway and transcription-factor enrichment; no single-cell or deconvolution analysis is included.

## 👤 Author & Contact

**Narges Shayesteh**
*Bioinformatics Specialist | Computational Transcriptomics*

- **GitHub:** [@shayesteh68](https://github.com/shayesteh68)
- **LinkedIn:** [Narges Shayesteh](https://www.linkedin.com/in/narges-shayesteh)
