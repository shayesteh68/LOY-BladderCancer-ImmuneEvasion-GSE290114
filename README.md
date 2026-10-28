# LOY bladder cancer RNA-seq analysis (GSE290114)

This repository compares CRISPR Y-chromosome knockout (`Y_KO`) with scrambled-control (`Y_Scr`) MB49 mouse bladder cancer cells. NCBI GEO lists six bulk RNA-seq samples: three per group. The scripts use DESeq2 for differential expression, then run KEGG enrichment, targeted heatmaps, and TRRUST transcription-factor enrichment.

## Inputs and outputs

- `data/raw/GSE290114_gene_count.txt.gz` - GEO gene-count matrix.
- `metadata/samplesheet.csv` - sample IDs and condition labels.
- `data/reference/mouse_gene2entrez.tsv` - mouse gene-to-Entrez mapping.
- `data/annotation/trrust_rawdata.mouse.tsv` - TRRUST mouse interactions.
- `scripts/01_deseq2_analysis.R` through `scripts/04_tf_enrichment.R` - analysis steps in order.
- `results/` - previously committed outputs and logs. These files were not regenerated during this audit.

## Setup and run

The Conda environment is named `loy-rnaseq`:

```bash
conda env create -f environment.yml
conda activate loy-rnaseq
Rscript scripts/01_deseq2_analysis.R
Rscript scripts/02_pathway_enrichment.R
Rscript scripts/03_plot_heatmap.R
Rscript scripts/04_tf_enrichment.R
```

The scripts expect the input files above and write analysis outputs under `results/`.

## Restoring missing inputs

If inputs are absent, `restore_loy_inputs.sh` can retrieve the GEO count matrix and rebuild the sample sheet and gene mapping. It requires Bash, `curl`, `gzip`, `awk`, `sort`, `md5sum`, and `git`; it does not run the R analysis. From the repository root, use:

```bash
bash restore_loy_inputs.sh --repo .
```

The restore script checks that `results/` remains unchanged and reports Git status. It does not stage, commit, or push files.

## Interpretation

The DESeq2 contrast is `Y_KO` versus `Y_Scr`; positive log2 fold change means higher expression in `Y_KO`. These results describe a six-sample experiment. Existing result files are historical snapshots and should be regenerated before being treated as verified by the current scripts.

Code in this repository is licensed under the MIT License (see `LICENSE`). GEO inputs and TRRUST resources remain subject to their respective source terms.
