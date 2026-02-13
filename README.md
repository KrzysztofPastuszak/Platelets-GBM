# Platelets-GBM

Platelet RNA profiling for glioblastoma: GSEA, TDEA, and Elastic Net analysis.

This repository contains code and data accompanying the manuscript:

> **Deepening the Understanding of Platelet RNA as Biomarker for Glioblastoma through GSEA, TDEA and Elastic Net Regularization Analysis**
>
> Anna Giczewska\*, Krzysztof Pastuszak\*, Leon van Hout, Edward Post, Jip Ramaker, Ken Zwaan, Anna Supernat\*\*, Bart A. Westerman\*\*

## Overview

Blood platelets alter their RNA content in response to tumor signals, a phenomenon referred to as tumor-educated platelets (TEPs). This project analyzes platelet RNA-seq data from 84 glioblastoma (GBM) patients and 319 healthy controls to investigate the molecular processes underlying platelet reprogramming in GBM. Two complementary analytical strategies are applied:

1. **GSEA** (Gene Set Enrichment Analysis) on the full expressed transcriptome, identifying broad pathway-level shifts between GBM and healthy platelets.
2. **TDEA with Elastic Net regularization** (Threshold-Based Differential Expression Analysis), a more targeted approach that first filters genes by differential expression thresholds, then fits a regularized logistic regression to select the most informative gene panel.

The GSEA arm highlights platelet regulatory processes such as coagulation, wound healing, and cytoskeletal organization. The Elastic Net arm surfaces immune activation, bacterial/DAMP response, and cell signaling pathways. Both point to systemic platelet involvement in the GBM tumor microenvironment.

## Repository structure

```
Platelets-GBM/
├── programs/
│   └── basic_model_analysis.R        # Main analysis script (TDEA + Elastic Net + enrichment)
├── data/
│   ├── rawdata/
│   │   ├── countsNormalized.tsv.zip  # DESeq2-normalized counts (zipped due to GitHub file size limits; unzip before use)
│   │   ├── sampleInfo.tsv            # Sample metadata (group, age, sex, isolation site, split)
│   │   ├── train_ids_final.tsv       # Training set sample IDs (60%)
│   │   └── testId_ids_final.tsv      # Test set sample IDs (40%)
│   ├── dgeGenesEnsembl75.RData       # Gene annotations (Ensembl 75, HGNC symbols + descriptions)
│   ├── c2.cp.reactome.v2022.1.Hs.symbols.gmt  # Reactome gene sets (MSigDB)
│   ├── enet/
│   │   └── allModelCoeffs_GBM.csv    # Elastic Net coefficients for the GBM classifier
│   └── pathways/
│       ├── GSEA/                     # GSEA results (HC vs GBM)
│       └── ORA_Enet/                 # GO ORA results for the Elastic Net gene panel
├── BM/                               # Brain metastasis analysis (supplementary)
│   ├── allModelCoeffs_BM.csv         # Elastic Net coefficients for the BM classifier
│   ├── BM_GO_all_pathways_final_model_train_split_60.xlsx
│   └── gene_set_enrichment_all_hc_vs_bm_ref_hc2.xlsx
├── Report/
│   ├── Confusion_matrix_test.pdf     # Confusion matrix heatmap (test set)
│   ├── reactome_top_10.pdf           # Top 10 Reactome pathways (bar chart)
│   ├── gene_ontology_top_10.pdf      # Top 10 GO Biological Process terms (bar chart)
│   ├── GO_BP_venn_adjusted.pdf       # Venn: significant GO BP terms (FDR-adjusted p < 0.05)
│   └── GO_BP_venn_unadjusted.pdf     # Venn: significant GO BP terms (unadjusted p < 0.05)
├── Platelets-GBM.Rproj              # RStudio project file
└── LICENSE                           # GPL-3.0
```

## Analysis pipeline

The `programs/basic_model_analysis.R` script runs the complete TDEA + Elastic Net pipeline:

1. **Load data.** DESeq2-normalized counts and sample metadata are read in. One IDH-mutant sample (reclassified as non-GBM under updated WHO criteria) and all NKI-site controls are excluded, leaving 84 GBM patients and 319 healthy controls.

2. **Train/test split.** The cohort is split 60/40 into training (51 GBM, 192 HC) and test (33 GBM, 127 HC) sets. The split used in the manuscript is stored in `train_ids_final.tsv` and `testId_ids_final.tsv` for reproducibility.

3. **Differential expression (training set only).** For each gene, a Wilcoxon rank-sum test compares expression between GBM and HC in the training set, followed by FDR correction. Genes meeting all three TDEA thresholds are retained as candidates:
   - FDR q-value < 0.1
   - Mean fold change > 1.3 (or < 1/1.3)
   - Median expression > 3 in at least one group

4. **Elastic Net classification.** A glmnet model is trained on the candidate genes using 5x5 repeated cross-validation with random hyperparameter search. Two rounds of modeling are performed: first on all candidate genes, then on the reduced panel of genes whose absolute coefficient exceeds 0.05.

5. **Pathway enrichment.** The final gene panel is tested for enrichment in Reactome pathways (ORA via `enricher`) and Gene Ontology Biological Process terms (ORA via `enrichGO`). Top 10 results are plotted as grouped bar charts showing the -log10 p-value alongside the percentage of genes selected vs. the background rate.

6. **Three-way GO comparison.** GO Biological Process enrichment is run separately on all three nested gene sets (TDEA candidates, non-zero Elastic Net coefficients, and the final panel). Venn diagrams visualize the overlap of significant terms at both FDR-adjusted and unadjusted p < 0.05 thresholds.

## Data

The RNA-seq data reanalyzes publicly available platelet transcriptome profiles from [Sol et al. (2020)](https://doi.org/10.1016/j.xcrm.2020.100101), GEO accession [GSE156050](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156050). Counts were normalized using the DESeq2 variance-stabilizing transformation.

## Requirements

R >= 4.1.2 with the following packages:

- DESeq2, edgeR
- caret, glmnet
- matrixStats, cvAUC
- clusterProfiler, org.Hs.eg.db
- ggplot2, dplyr
- openxlsx
- party
- VennDiagram

## Running the analysis

```bash
# First, unzip the normalized counts (compressed due to GitHub file size limits)
cd data/rawdata
unzip countsNormalized.tsv.zip
cd ../../
```

```r
# Open the RStudio project
# Set working directory to programs/
setwd("programs")

# Run the full pipeline
source("basic_model_analysis.R")
```

Note that the Elastic Net training step involves repeated cross-validation and may take several minutes depending on hardware. Outputs are written to `Report/` (PDF figures) and `programs/` (enrichment tables and intermediate RData files).

## Key results

The Elastic Net classifier achieved an AUC of 0.91 on the held-out test set. The final gene panel includes genes involved in immune activation (B2M, LGALS2), PI3K signaling (PIK3IP1), calcium-dependent kinase activity (CAMK2D), and chromatin regulation (CBX7), among others. Pathway analysis of this panel highlighted immune response, bacterial/DAMP sensing, and protein folding pathways, complementing the broader coagulation and wound healing signatures identified by GSEA.

## Citation

If you use this code or data, please cite:

> Giczewska A\*, Pastuszak K\*, van Hout L, Post E, Ramaker J, Zwaan K, Supernat A\*\*, Westerman BA\*\*. Deepening the Understanding of Platelet RNA as Biomarker for Glioblastoma through GSEA, TDEA and Elastic Net Regularization Analysis.

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).
