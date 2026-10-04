# End-to-End Transcriptomic Dissection of LOY-Driven Immune Evasion in Bladder Cancer

[![Bioinformatics](https://img.shields.io/badge/Domain-Bioinformatics%20%7C%20Oncogenomics-blue.svg)](https://github.com/shayesteh68)
[![Pipeline](https://img.shields.io/badge/Pipeline-DESeq2%20%7C%20fgsea%20%7C%20Seurat-brightgreen.svg)](https://github.com/shayesteh68)
[![Dataset](https://img.shields.io/badge/GEO-GSE290113-orange.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290113)

## 📌 Scientific Context & Objectives
Loss of the Y chromosome (LOY) in malignant epithelial cells and infiltrating tumor-infiltrating lymphocytes (TILs) is a major somatic driver of aggressive tumor behavior and immune checkpoint resistance in male bladder carcinomas.

This project delivers a robust, publication-grade computational framework to dissect:
1. **Differential Expression Analysis (DESeq2):** Contrasting CRISPR-mediated Y-chromosome knockout (`CRISPR Y-KO`) against control scrambled models (`CRISPR Y-Scr`) in MB49 urothelial cancer cells (Accession: **GSE290113**).
2. **Pathway & Gene Set Enrichment Analysis (fgsea):** Mapping down-regulated antigen presentation machinery, interferon-gamma response depletion, and active immunosuppressive signaling cascades.
3. **Tumor Microenvironment Deconvolution:** Identifying transcriptomic proxies of CD8+ T-cell exhaustion associated with LOY.

---

## 🔬 Dataset Overview
- **Accession:** NCBI GEO [GSE290113](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290113)
- **Reference Article:** *Concurrent loss of the Y chromosome in cancer and T cells impacts outcome* (*Nature*, 2025).
- **Organism:** *Mus musculus*
- **Experimental Groups:**
  - `CRISPR Y-Scr` (Control / Y-intact MB49 bladder cancer line)
  - `CRISPR Y-KO` (CRISPR-deleted entire Y chromosome MB49 line)

---

## 🗂️ Repository Structure

```text
Theodorescu-LOY-Immune-Evasion-Pipeline/
├── data/
│   ├── raw/                   # Raw counts / GEO supplementary matrices
│   └── metadata/              # Sample phenotype annotations & batch info
├── scripts/
│   ├── 01_qc_normalization.R  # Count QC, filtering, and VST/rlog transformation
│   ├── 02_deseq2_differential.R # Negative binomial GLM fitting (Y-KO vs Y-Scr)
│   ├── 03_gsea_enrichment.R   # MSigDB Hallmarks & Reactome pathway analysis
│   └── 04_publication_figures.R # High-resolution Volcano, Heatmaps, and PCA plots
├── results/
│   ├── tables/                # DEG lists (padj < 0.05, |log2FC| > 1.0)
│   └── figures/               # Publication-ready vector & raster exports (300 DPI)
├── environment.yml            # Conda environment specifications
└── README.md                  # Project documentation
```

---

## ⚙️ Environment Setup & Reproducibility

### Conda Environment Specification
```bash
# Clone the repository
git clone https://github.com/shayesteh68/Theodorescu-LOY-Immune-Evasion-Pipeline.git
cd Theodorescu-LOY-Immune-Evasion-Pipeline

# Create and activate environment
conda env create -f environment.yml
conda activate loy-rnaseq-env
```

### Core Dependencies
- **R version:** `>= 4.3.0`
- **Bioconductor Packages:** `DESeq2 (v1.40.0)`, `fgsea (v1.26.0)`, `ComplexHeatmap (v2.16.0)`, `EnhancedVolcano (v1.18.0)`, `org.Mm.eg.db (v3.17.0)`
- **CRAN Packages:** `tidyverse (v2.0.0)`, `pheatmap (v1.0.12)`, `patchwork (v1.1.3)`

---

## 🚀 Pipeline Execution

```bash
# Step 1: Preprocessing and Quality Control
Rscript scripts/01_qc_normalization.R

# Step 2: Differential Expression Modeling
Rscript scripts/02_deseq2_differential.R

# Step 3: Fast Gene Set Enrichment Analysis
Rscript scripts/03_gsea_enrichment.R

# Step 4: Generate Publication Figures
Rscript scripts/04_publication_figures.R
```

---

## 👤 Author & Contact
**Narges Shayesteh**  
*Bioinformatics Specialist | Computational Transcriptomics*  
- **GitHub:** [@shayesteh68](https://github.com/shayesteh68)  
- **LinkedIn:** [Narges Shayesteh](https://www.linkedin.com/in/narges-shayesteh)
