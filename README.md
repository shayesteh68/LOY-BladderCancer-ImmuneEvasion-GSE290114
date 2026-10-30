# LOY bladder cancer RNA-seq analysis (GSE290114)

This repository analyzes the six bulk RNA-seq samples from the MB49 mouse bladder cancer experiment in GEO series GSE290114: three Y-chromosome knockout (`Y_KO`) replicates and three scrambled-control (`Y_Scr`) replicates.

## Workflow and current outputs

Run the scripts in order from the repository root:

```bash
Rscript scripts/01_deseq2_analysis.R
Rscript scripts/02_pathway_enrichment.R
Rscript scripts/03_plot_heatmap.R
Rscript scripts/04_tf_enrichment.R
```

The current committed run contains 765 DEGs at BH adjusted p-value < 0.05 and |log2 fold change| >= 1 (181 higher and 584 lower in Y_KO). The KEGG step found 48 significant pathways for genes lower in Y_KO and none for genes higher in Y_KO. TRRUST returned nine enriched candidate transcription factors among the downregulated targets. The heatmap covers all 12 requested targets; `results/tables/heatmap_target_gene_coverage.csv` records which genes were plotted.

The scripts write DESeq2 tables and summaries to `results/tables/`, and figures to `results/figures/`. KEGG annotations are retrieved from the KEGG REST service at run time, so the live database can change; the run summary records the date and annotation source. Raw GEO counts and annotation resources are not tracked in Git.

## Inputs and restoring them

- `data/raw/GSE290114_gene_count.txt.gz` - GEO gene-count matrix.
- `metadata/samplesheet.csv` - sample IDs and condition labels; this file is tracked.
- `data/reference/mouse_gene2entrez.tsv` - mouse symbol-to-Entrez mapping.
- `data/annotation/trrust_rawdata.mouse.tsv` - TRRUST mouse transcription-factor interactions.

If the ignored data files are absent, restore them from the repository root:

```bash
bash restore_loy_inputs.sh --repo .
```

The restore script checks downloaded inputs and reports Git status. It does not run the R analysis or stage, commit, or push files. GEO and TRRUST resources remain subject to their source terms.

## Environment

The scripts run on Linux or WSL with Conda. For the audited Linux x86_64 package set:

```bash
conda create --name loy-rnaseq --file environment-linux-64.lock.txt
conda activate loy-rnaseq
```

`environment.yml` contains the pinned direct dependencies for other compatible Conda platforms:

```bash
conda env create --file environment.yml
conda activate loy-rnaseq
```

The lock records R 4.4.3, DESeq2 1.46.0, clusterProfiler 4.14.0, GO.db 3.20.0, ggplot2 4.0.3, and pheatmap 1.0.13.

## Interpretation

The DESeq2 contrast is `Y_KO` versus `Y_Scr`; positive log2 fold change means higher expression in Y_KO. The model compares three replicates per condition in this MB49 experiment and does not support generalization to other cell lines or tumors.

The KEGG over-representation background is the set of tested genes that map to Entrez IDs. The TRRUST background is the set of tested genes represented in TRRUST. Enriched TRRUST factors are candidate regulators whose recorded targets overlap the downregulated gene set; this analysis does not measure transcription-factor expression or activity and does not establish causality.

The heatmap shows log2 normalized counts with each gene scaled across samples. It is a descriptive visualization, not a separate statistical test.

## License

Code in this repository is licensed under the MIT License; see [LICENSE](LICENSE). GEO and TRRUST inputs retain their respective source terms.
