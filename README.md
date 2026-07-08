# scGeneExpressionExplorer README

This folder contains the Shiny application and analysis helper scripts for single-cell gene expression exploration, differential expression analysis, GSEA, UMAP/cell-type marker visualization, gene expression summary export, and CellChat-based cell-cell communication analysis.

## Code sources and acknowledgements

Part of the original analysis code was adapted from Max's codebase.

- `DE_analysis.R` uses Max's `FindMarkersLoupe` workflow and notes the original Max code location as `/projectnb/wax-dk/max/RSRC/G190`.
- GitHub dependency from Max: `mpyatkov/FindMarkersLoupe`
  - Install command: `devtools::install_github("mpyatkov/FindMarkersLoupe")`
- `app.R` also uses `NotationConverter`, with the commented install source:
  - `devtools::install_github("mpyatkov/NotationConverter")`

The CellChat code used by this project is a modified fork of CellChat:

- Modified CellChat repository: <https://github.com/Ye1203/CellChat/tree/415849a378014b220e8f82e04c3835f346ed92ad>
- This fork is based on `jinworks/CellChat` and should be used instead of installing the default upstream CellChat when reproducing this app environment.

## Files

| File | Purpose |
|---|---|
| `app.R` | Main Shiny application. Provides UI and server logic for loading Seurat, DEG, GSEA, and CellChat objects; plotting UMAPs, dot plots, heatmaps, volcano plots, GSEA results, and CellChat visualizations. |
| `DE_analysis.R` | Differential expression helper code. Uses Max's `FindMarkersLoupe`/SEGEX-style workflow to compare selected sample and cluster groups and export DEG tables. |
| `GSEA_analysis.R` | Pre-ranked GSEA helper code using DEG output. Supports ranking by `fc`, `log2fc`, `fc-pvalue`, or `log2fc-pvalue`, and gene-set databases from MSigDB through `msigdbr`. |
| `UmapPlot.R` | UMAP, dot plot, cluster summary table, marker heatmap, and cell-type marker helper functions. Includes marker lists and cluster annotation utilities. |
| `gene_expression_output.R` | Creates gene expression summary tables by one or two metadata variables. Outputs mean expression, relative expression, and non-zero expression percentage. |
| `CCC.R` | Front-end/helper functions for CellChat analysis submission. Writes parameter files, copies the Seurat object, generates a `qsub` script, submits SCC jobs, and provides a diagnostic function for choosing the CellChat expression-summary method. |
| `cellchat_runner.R` | Command-line CellChat runner used by SCC jobs. Loads the Seurat object, splits control/treatment groups if needed, runs CellChat, saves `.rds` results, exports Excel communication tables, writes summaries, and optionally sends completion/error emails. |
| `renv.lock` | Locked R package environment for reproducibility. |

## Main functionality

### 1. Shiny app

Run the main app with:

```r
shiny::runApp("app.R")
```

The app supports:

- Seurat object loading
- metadata-driven sample/cluster comparison setup
- differential expression analysis
- GSEA visualization
- UMAP, dot plot, heatmap, and expression summary outputs
- CellChat result loading and visualization
- CellChat analysis submission to SCC using `qsub`

### 2. Differential expression

`DE_analysis.R` defines `compute_de()`, which compares selected sample/cluster groups by creating temporary identities named `GROUP1` and `GROUP2`. It then runs `FindMarkersLoupe()` and converts output to SEGEX-compatible format.

Important note: this workflow follows a pseudo-bulk/SEGEX-compatible style. For general scRNA-seq or snRNA-seq differential expression, especially when not constrained by SEGEX compatibility, `Seurat::FindMarkers()` may be more appropriate.

### 3. GSEA

`GSEA_analysis.R` defines `gsea_analysis_from_tsv()`, which reads DEG output, builds a ranked gene list, retrieves gene sets using `msigdbr`, and runs `clusterProfiler::GSEA()`.

Supported ranking methods:

- `fc`
- `log2fc`
- `fc-pvalue`
- `log2fc-pvalue`

Recommended interpretation: in single-cell analysis without true biological replicates, focus mainly on NES direction and magnitude rather than over-interpreting p-values.

### 4. Gene expression summary

`gene_expression_output.R` defines `create_expression_summary()`. For each selected gene and metadata group, it reports:

- mean expression
- relative expression, normalized by the maximum mean expression for that gene across groups
- non-zero proportion (%)

This is useful for exporting compact expression tables for marker genes or selected genes.

### 5. UMAP and cell-type marker visualization

`UmapPlot.R` provides helper functions for:

- UMAP plotting
- dot plots
- cluster summary tables
- marker-based cell-type prediction
- marker heatmaps
- combined UMAP/table/dotplot/heatmap layouts
- optional UMAP reclustering utilities

### 6. CellChat analysis

CellChat analysis is split into two layers:

- `CCC.R`: prepares and submits jobs from the Shiny app.
- `cellchat_runner.R`: runs CellChat on SCC as a command-line script.

The CellChat runner supports:

- single combined analysis (`All Together`)
- CONTROL vs TREATMENT comparison
- species selection: human or mouse
- pathway type filtering, such as `Secreted Signaling`, `Cell-Cell Contact`, `ECM-Receptor`, and `Non-protein Signaling`
- probability method selection: `triMean`, `truncatedMean_0.1`, or `truncatedMean_0.05`
- minimum cell filtering
- merged CellChat output for comparison visualization
- Excel export of `net` and `netP` communication tables
- run summary text files
- optional email notification

## CellChat code modification

This project should use the modified CellChat fork here:

<https://github.com/Ye1203/CellChat>

Use this version when restoring or rebuilding the environment:

```r
renv::install("Ye1203/CellChat")
renv::snapshot()
```

The app-side CellChat code was also customized in `app.R`, `CCC.R`, and `cellchat_runner.R` to support this workflow. Main project-level changes include:

- Shiny UI for selecting CellChat input parameters.
- SCC job submission through generated `qsub` scripts.
- Saving a reproducible parameter file for each CellChat run.
- Automatic CONTROL/TREATMENT splitting based on selected sample metadata.
- Automatic addition of a `C` prefix to numeric cluster labels containing `0`, avoiding CellChat or plotting issues with numeric/zero cluster IDs.
- Pathway-type filtering before CellChat probability calculation.
- Support for multiple probability summary methods: `triMean`, `truncatedMean(0.1)`, and `truncatedMean(0.05)`.
- Diagnostic function to recommend a probability method based on ligand/receptor signal retention across clusters.
- Export of CellChat communication tables to Excel.
- Support for merged CellChat objects and comparison visualizations in the Shiny app.

## Environment setup with renv

The project uses `renv` for reproducible R package management. The provided `renv.lock` was generated for:

- R version: `4.4.3`
- Bioconductor version: `3.20`
- CRAN repository: `https://cloud.r-project.org`

### 1. Start R in the project folder

```bash
cd /path/to/scGeneExpressionExplorer
module load R/4.4.3   # on BU SCC, if applicable
R
```

### 2. Install renv if needed

```r
install.packages("renv")
```

### 3. Restore the locked environment

```r
renv::restore()
```

This installs package versions from `renv.lock` into the project-local `renv/library`.

### 4. Install the modified CellChat fork

If `renv::restore()` installs upstream CellChat or an older CellChat commit, reinstall the modified fork:

```r
renv::install("Ye1203/CellChat")
renv::snapshot()
```

### 5. Install Max-related GitHub dependencies if missing

```r
renv::install("mpyatkov/FindMarkersLoupe")
renv::install("mpyatkov/NotationConverter")
renv::snapshot()
```

If `NotationConverter` is only available locally on SCC, install it from the local path used by the project/team instead of GitHub.

## Running on SCC

The CellChat analysis workflow is designed for SCC job submission.

`CCC.R` generates a job script similar to:

```bash
module load R/4.4.3
Rscript cellchat_runner.R \
  --seurat_file=/path/to/input.rds \
  --sample_info_file=/path/to/cellchat_sample_info.rds \
  --cluster_column=seurat_clusters \
  --species=mouse \
  --pathway_type='Secreted Signaling|Cell-Cell Contact' \
  --prob_type=triMean \
  --min_cells=10 \
  --cores=8 \
  --save_path=/path/to/output \
  --email=NA
```

The generated output folder may include:

- `cellchat_parameters.txt`
- `cellchat_sample_info.rds`
- `cellchat_analysis.qsub`
- `ccc.log`
- `cellchat_object.rds` for all-together analysis
- `cellchat_compare.rds` for CONTROL/TREATMENT analysis
- `cellchat_tables.xlsx`
- `analysis_summary.txt`

## Recommended project structure

```text
scGeneExpressionExplorer/
├── app.R
├── DE_analysis.R
├── GSEA_analysis.R
├── UmapPlot.R
├── gene_expression_output.R
├── CCC.R
├── cellchat_runner.R
├── renv.lock
├── renv/
└── README.md
```

## Notes and caveats

- Keep `app.R`, helper scripts, and `renv.lock` in the same project folder unless paths are updated manually.
- `cellchat_runner.R` currently contains a hard-coded `project_root`; update it if the project is moved to a different SCC directory.
- For large Seurat objects and CellChat jobs, run on SCC rather than a local laptop.
- CellChat outputs can be large. Store results in a dedicated output directory for each run.
- After changing package versions, run `renv::snapshot()` so the lockfile stays synchronized with the working environment.
