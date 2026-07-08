# scGeneExpressionExplorer

<div align="center">

**An interactive Shiny application for downstream single-cell and single-nucleus RNA-seq analysis**

Differential Expression • GSEA • CellChat • Gene Expression Explorer • UMAP Visualization • SCC Integration

</div>

---

## Overview

**scGeneExpressionExplorer** is an R/Shiny application designed to simplify downstream analysis of single-cell (scRNA-seq) and single-nucleus (snRNA-seq) transcriptomic data.

Instead of requiring users to manually write R scripts for each analysis, the application provides an intuitive graphical interface that integrates several commonly used downstream analyses into a single workflow.

The application currently supports:

- Differential gene expression analysis
- Gene Set Enrichment Analysis (GSEA)
- Gene expression visualization
- Cell-cell communication analysis using CellChat
- Interactive visualization of UMAPs, marker genes, heatmaps, volcano plots, and communication networks
- Automated CellChat execution on Boston University's Shared Computing Cluster (SCC)

The project is intended for researchers who wish to perform reproducible downstream analyses on processed Seurat objects without extensive programming experience while still maintaining compatibility with standard R workflows.

---

## Tutorial

A detailed tutorial describing the graphical interface and example workflows is available here:

**[How to use scGeneExpressionExplorer.pdf](How%20to%20use%20scGeneExpressionExplorer.pdf)**

The tutorial includes:

- Data preparation
- Differential expression analysis
- Gene Set Enrichment Analysis (GSEA)
- Cell-cell communication analysis
- Gene expression visualization
- Example screenshots of the Shiny interface

# Features

## Differential Expression Analysis

- Compare any user-selected sample and cluster combinations
- SEGEX-compatible differential expression output
- Export complete DEG tables
- Interactive volcano plots
- Multiple comparison support

---

## Gene Set Enrichment Analysis (GSEA)

- Pre-ranked GSEA
- Multiple ranking strategies
- GO Biological Process
- Reactome
- Hallmark
- MSigDB integration via **msigdbr**
- Interactive enrichment visualization

---

## Gene Expression Explorer

Visualize gene expression across:

- Cell clusters
- Sample groups
- Experimental conditions

Generate summary tables including:

- Mean expression
- Relative expression
- Percentage of expressing cells

---

## Cell-Cell Communication Analysis

Integrated CellChat workflow supporting:

- Human and mouse databases
- Single dataset analysis
- CONTROL vs TREATMENT comparison
- Pathway-specific visualization
- Ligand-receptor exploration
- Network centrality analysis
- Communication heatmaps
- Bubble plots
- Circle plots
- Chord diagrams
- Pathway contribution analysis

---

## UMAP Visualization

Built-in visualization utilities include:

- UMAP
- DotPlot
- Marker heatmap
- Cluster summary tables
- Automatic cell-type annotation
- Marker-based cluster prediction

---

## Boston University SCC Integration

Large CellChat analyses can be submitted directly to BU SCC through the graphical interface.

The application automatically:

- Generates parameter files
- Creates qsub scripts
- Copies required input files
- Submits CellChat jobs
- Collects results
- Exports communication tables
- Generates summary reports

---

# Workflow

The typical workflow is illustrated below.

```text
Seurat Object
      │
      ▼
Load into Shiny
      │
      ├──────── Differential Expression
      │                 │
      │                 ▼
      │            DEG Table
      │                 │
      │                 ▼
      │               GSEA
      │
      ├──────── Gene Expression Explorer
      │
      ├──────── UMAP Visualization
      │
      └──────── CellChat Analysis
                        │
                        ▼
          SCC (optional for large datasets)
                        │
                        ▼
        Communication Networks & Figures
```

---

# Tutorial

A detailed user guide is provided in:

**📄 How to use scGeneExpressionExplorer.pdf**

The tutorial includes:

- Data preparation
- Loading Seurat objects
- Differential expression analysis
- GSEA
- CellChat analysis
- Gene expression visualization
- Example screenshots
- Typical analysis workflow

---

# Repository Structure

```
scGeneExpressionExplorer/
│
├── app.R
├── CCC.R
├── cellchat_runner.R
├── DE_analysis.R
├── GSEA_analysis.R
├── gene_expression_output.R
├── UmapPlot.R
├── renv.lock
├── renv/
├── tutorial/
│     └── How to use scGeneExpressionExplorer.pdf
│
└── README.md
```

---

# Repository Contents

| File | Description |
|------|-------------|
| `app.R` | Main Shiny application containing both the user interface and server logic. |
| `CCC.R` | CellChat helper functions for generating SCC jobs and parameter files. |
| `cellchat_runner.R` | Command-line CellChat workflow executed on SCC. |
| `DE_analysis.R` | Differential expression analysis using the SEGEX-compatible workflow. |
| `GSEA_analysis.R` | Gene Set Enrichment Analysis utilities based on DEG output. |
| `gene_expression_output.R` | Gene expression summary table generation. |
| `UmapPlot.R` | UMAP visualization, marker heatmaps, cell-type annotation, and plotting utilities. |
| `renv.lock` | Locked package versions for reproducible installation. |

---

# Installation

## Requirements

- R ≥ 4.4
- Bioconductor ≥ 3.20
- Git
- renv

---

## Clone the repository

```bash
git clone https://github.com/your_repository/scGeneExpressionExplorer.git

cd scGeneExpressionExplorer
```

---

## Install renv

```r
install.packages("renv")
```

---

## Restore the project environment

```r
renv::restore()
```

This will install all package versions recorded in **renv.lock**.

---

## Install the customized CellChat package

This project uses a customized version of CellChat instead of the original package.

Install it using

```r
renv::install("Ye1203/CellChat")
```

or

```r
devtools::install_github("Ye1203/CellChat")
```

---

## Install additional dependencies

Some helper packages were originally developed by Max and may need to be installed separately.

```r
renv::install("mpyatkov/FindMarkersLoupe")

renv::install("mpyatkov/NotationConverter")
```

If `NotationConverter` is only available locally on SCC, install it from the local project directory instead.

---

# Quick Start

Launch the application with

```r
shiny::runApp("app.R")
```

The graphical interface allows users to

1. Load Seurat objects
2. Configure sample metadata
3. Perform differential expression analysis
4. Run GSEA
5. Visualize gene expression
6. Submit CellChat analyses
7. Explore CellChat results interactively

---
# CellChat Workflow

Cell-cell communication analysis in **scGeneExpressionExplorer** is implemented as a two-stage workflow.

The graphical interface is responsible for collecting user parameters and preparing the analysis, while the computationally intensive CellChat workflow is executed as a standalone command-line script on Boston University's Shared Computing Cluster (SCC).

```
Shiny App
    │
    ▼
Generate parameter files
    │
    ▼
Generate qsub script
    │
    ▼
Submit SCC Job
    │
    ▼
cellchat_runner.R
    │
    ▼
CellChat Analysis
    │
    ▼
Save Results
    │
    ▼
Interactive Visualization
```

The CellChat workflow supports:

- Single dataset analysis
- CONTROL vs TREATMENT comparison
- Human and mouse databases
- Automatic pathway filtering
- Multiple communication probability estimation methods
- Parallel execution
- Excel export
- Communication network visualization

---

# Customized CellChat Package

This project uses a customized version of CellChat rather than the original package.

Repository:

https://github.com/Ye1203/CellChat

The customized package contains numerous improvements that were developed specifically for this application.

# Application-level CellChat Modifications

Besides modifications to the CellChat package itself, the Shiny application contains substantial workflow extensions implemented in

- `app.R`
- `CCC.R`
- `cellchat_runner.R`

These additions include

## Interactive Parameter Selection

The graphical interface allows users to configure

- species
- pathway categories
- probability estimation method
- minimum cell threshold
- sample grouping
- comparison settings

without writing any R code.

---

## Automated SCC Submission

The application automatically

- saves analysis parameters
- copies the Seurat object
- generates a reproducible qsub script
- submits the analysis
- monitors execution
- stores all outputs

This eliminates manual command-line preparation.

---

## Reproducible Analysis

Each CellChat analysis automatically records

- input Seurat object
- sample information
- selected parameters
- runtime settings
- CPU allocation
- generated qsub script

making analyses fully reproducible.

---

## Automatic Metadata Processing

The workflow additionally performs several preprocessing steps automatically.

Examples include

- CONTROL/TREATMENT splitting
- cluster label validation
- automatic "C" prefix for numeric clusters
- pathway database filtering
- merged CellChat generation

These preprocessing steps reduce common user errors while preserving compatibility with CellChat.

---

## Diagnostic Utilities

The application also provides helper functions for selecting an appropriate communication probability summary method.

Users can compare

- triMean
- truncatedMean (0.1)
- truncatedMean (0.05)

according to ligand-receptor signal retention across clusters before running CellChat.

---

# Differential Expression Analysis

Differential expression analysis is implemented in `DE_analysis.R`.

The workflow is adapted from Max's original **FindMarkersLoupe** implementation and generates SEGEX-compatible differential expression tables.

Compared with standard Seurat workflows, this implementation emphasizes compatibility with downstream SEGEX analysis pipelines.

For conventional scRNA-seq studies, users may instead choose to perform differential expression using

```r
Seurat::FindMarkers()
```

before importing results into downstream analyses.

---

# Gene Set Enrichment Analysis

The GSEA module performs pre-ranked enrichment analysis using differential expression output.

Supported ranking methods include

- fc
- log2fc
- fc × -log10(p-value)
- log2fc × -log10(p-value)

Supported databases include

- GO Biological Process
- Reactome
- Hallmark

through the **msigdbr** package.

For single-cell datasets without biological replicates, interpretation should primarily focus on the **Normalized Enrichment Score (NES)** rather than statistical significance.

---

# Running on Boston University SCC

Large CellChat analyses are intended to be executed on BU SCC.

The application automatically generates a command similar to

```bash
module load R/4.4.3

Rscript cellchat_runner.R \
    --seurat_file input.rds \
    --sample_info_file sample_info.rds \
    --cluster_column seurat_clusters \
    --species mouse \
    --pathway_type "Secreted Signaling|Cell-Cell Contact" \
    --prob_type triMean \
    --min_cells 10 \
    --cores 8 \
    --save_path output \
    --email your_email
```

The output directory typically contains

- CellChat objects
- merged CellChat objects
- Excel communication tables
- parameter files
- qsub scripts
- analysis summaries
- log files

---

# Code Sources and Acknowledgements

This project incorporates and extends several existing open-source tools.

## Max's Analysis Code

Parts of the differential expression workflow were adapted from Max's internal analysis code.

These components include

- FindMarkersLoupe workflow
- SEGEX-compatible DEG export
- UMAP helper utilities
- Marker visualization functions

Original code references are preserved within the corresponding source files whenever applicable.

---

## Open-source Packages

This project is built upon numerous outstanding open-source R packages, including

- Seurat
- CellChat
- clusterProfiler
- msigdbr
- ComplexHeatmap
- ggplot2
- patchwork
- openxlsx
- future
- shiny

The developers of these packages are gratefully acknowledged.

---

# Reproducibility

This project uses **renv** for package management.

The complete software environment is recorded in

```
renv.lock
```

To restore the environment

```r
install.packages("renv")

renv::restore()
```

Whenever package versions are updated, synchronize the lock file using

```r
renv::snapshot()
```

---

# Frequently Asked Questions

### Why use a customized CellChat package?

The customized version contains visualization improvements, workflow extensions, and bug fixes that are required by the Shiny application.

---

### Can I use the original CellChat package?

The application may still run, but some visualization modules and comparison plots may not behave as expected. The customized package is therefore recommended.

---

### Can I run CellChat locally?

Yes.

However, large datasets are recommended to be analyzed on BU SCC because CellChat can be computationally intensive.

---

### Is this application limited to mouse data?

No.

Both mouse and human CellChat databases are supported.

---

# Citation

If you use this application in your research, please cite the original publications of the software packages used in your analyses, including

- Seurat
- CellChat
- clusterProfiler
- msigdbr

Please also acknowledge this repository whenever appropriate.

---

# Contact

**Bingtian Ye**

Graduate Student  
Boston University

Email:

- btye@bu.edu
- biangtian@icloud.com

GitHub:

https://github.com/Ye1203

---

# License

This repository contains original code together with modifications of several open-source projects.

Please follow the licenses of the corresponding upstream projects when redistributing or extending this software.