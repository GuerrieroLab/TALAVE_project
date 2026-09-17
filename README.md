# TALAVE Project

This repository contains the end-to-end TALAVE analysis workflow across clinical metadata, NanoString IO360, DSP protein, CyCIF, and integrated cross-platform analyses.

## Repository structure

- `Rmd/`: the seven figure-generating notebooks shipped in this Zenodo supplement (IO360, DSP protein, CyCIF, integrated)
- `functions/`: shared bootstrap, I/O, analysis, QC, and plotting helpers

This is a figure-generation-only supplement: the upstream data-loading, QC, cell-type
calling, and batch-correction notebooks that produced the checkpoint objects each notebook
consumes are not included here. `data/` (precomputed checkpoint objects each notebook loads
from) and `figures/` (generated outputs) are not part of this GitHub repository — they are
distributed via the accompanying Zenodo deposit, one archive per modality folder under
`data/obj/` (see "Downloading the code and data" below).

## Downloading the code and data

Both the code and the data checkpoints are on Zenodo: **[DOI: 10.5281/zenodo.21815064]**
(placeholder -- replace with the actual DOI/link once published).

The Zenodo record has five files: `TALAVE_project_code.zip` (this repository's tracked
contents) and four data archives, one per modality: `clinical.zip`, `io360.zip`,
`protein.zip`, `cycif.zip`.

```bash
# 1. Code: extract to wherever you want the repository to live
unzip TALAVE_project_code.zip -d TALAVE_project
cd TALAVE_project

# 2. Data: each archive's contents go directly into data/obj/<modality>/
mkdir -p data/obj/clinical data/obj/io360 data/obj/protein data/obj/cycif
unzip /path/to/clinical.zip -d data/obj/clinical
unzip /path/to/io360.zip    -d data/obj/io360
unzip /path/to/protein.zip  -d data/obj/protein
unzip /path/to/cycif.zip    -d data/obj/cycif
```

The result should be a `data/obj/` layout matching exactly what's described in "Shared
inputs and metadata" below and what `functions/config_paths.R` expects (`clinical_obj_dir`,
`io360_obj_dir`, `protein_obj_dir`, `cycif_obj_dir` are each `file.path(root_dir, "data",
"obj", "<modality>")`) -- no notebook reads from anywhere else under `data/`. Then follow
"Setup" below to point `root_dir` at wherever you extracted the code.

## Analysis workflow

The notebooks are organized as a staged pipeline:

1. `Rmd/1_io360_expression_analysis.Rmd`
   Per-gene mixed-effects modeling and all gene/pathway figures — cGAS-STING, DDR, volcano plots, FGSEA results, cross-platform heatmaps.
2. `Rmd/2_dsp_protein_expression_analysis.Rmd`
   Models protein-level longitudinal changes by segment and produces contrast summaries and protein heatmaps.
3. `Rmd/3_cycif_bulk_celltype_analysis.Rmd`
   Uses batch-corrected CyCIF objects to quantify cell-type frequencies, tumor frequencies, marker-positive populations, double-positive states, T-cell stratifications, and tumor/stroma classification.
4. `Rmd/4_cycif_umap_analysis.Rmd`
   Runs UMAP-based cell-state analysis (All Cells, Tumor Cells, Macrophages), producing the macrophage clustering later merged with `5`'s TME clustering.
5. `Rmd/5_cycif_spatial_frnn_tme.Rmd`
   Loads the spatial FRNN checkpoint, generates [EDF2C] cell-type hierarchy, [Fig 2C] tumor cell density, and performs TME-level Leiden clustering with its UMAP/cluster-frequency figures.
6. `Rmd/6_cycif_spatial_continued.Rmd`
   Continues the spatial pipeline: adjacent-cell and tumor-border distances, tumor density, macrophage and T-cell spatial programs, TS ratios, and marker-density regressions.
7. `Rmd/7_integrated_cross_platform_analysis.Rmd`
   Harmonizes shared markers across IO360 and DSP protein data, subsets shared clinical cohorts, computes cross-platform correlation summaries, compares CyCIF with other technologies, and fits integrated CD8 signal models.

## R objects used for the CyCIF project

[`cycif_object_reference.md`](cycif_object_reference.md) is a reference glossary for the R
objects passed between the CyCIF notebooks (`3`-`6`) -- what each object is, which notebook
creates it, which checkpoint file it's saved to/loaded from, and its row/column structure.
Not part of the run sequence itself (no code chunks); a lookup for when an object name in
those notebooks isn't self-explanatory.

## Shared inputs and metadata

All inputs are pre-compiled `.rda`/`.rds` checkpoints organized one folder per modality
under `data/obj/`. No notebook reads a raw xlsx/csv/tsv file. On Zenodo, each modality
folder below is distributed as its own archive (see "Downloading the code and data" above).

**`data/obj/clinical/`**
- `clinical_metadata.rds`: harmonized clinical metadata, loaded via `load_clinical_metadata()` (`functions/io_clinical.R`) in notebook `7`. A pre-specified patient-exclusion criterion has already been applied upstream in every checkpoint object -- no notebook filters on patient ID.
- `patient_colors.rds`: shared patient color palette, loaded via `load_pt_colors()` (`functions/plotting.R`) in notebooks `1`-`6`.

**`data/obj/io360/`** (notebooks `1`, `7`)
- `io360_expression_input_data.rda`: IO360 expression/pathway checkpoint.
- `sting_from_literature.rds`: cGAS-STING literature gene signature.
- `msigdb_v7.5.1_almac.rds`: third-party MSigDB reference database, used for FGSEA.
- `io360_gene_contrast_results.rds`: precomputed per-gene LMER contrast results (p-values/FDR) feeding the volcano plots -- notebook `1` regenerates and saves this if missing, but it's slow to rebuild, so include it.

**`data/obj/protein/`** (notebooks `2`, `7`)
- `dsp_protein_expression_input_data.rda`: DSP protein expression checkpoint, loaded via `get_dsp_protein_dirs()` (`functions/io_dsp_protein.R`).
- `protein_gene_lookup.rds`: gene/protein name lookup table, loaded via `load_dsp_protein_lookup()`.

**`data/obj/cycif/`** (notebooks `3`-`7`)
- `cs5.rds` (notebooks `3`, `5`, `7`), `cs7.rds` (notebooks `3`-`7`): the 5-cycle and 7-cycle
  `CycifStack` objects, split into separate files so a notebook that only needs one (`3`, `5`,
  `7` need both; `4`, `6` need only `cs7`) doesn't pay the load cost of the other -- `cs5`
  alone is ~1.8GB in memory.
- `cycif_spatial_frnn_input_data.rda`: the remaining spatial FRNN checkpoint objects
  (`frnn7.cmb`, `frnn7.exp`, `frnn7.rcn`, `is.sele.7`), loaded by notebooks `4`, `5`, `6`.
- `cycif5_neighborhood_profiles.rda`, `cycif7_neighborhood_profiles.rda`: 5-cycle/7-cycle recurrent-cellular-neighborhood composition + expression checkpoints.
- `cycif7_macrophage_clusters.rds`: `frnn7.cmb` with `cluster_mac` populated, saved by
  notebook `4`.
- `cycif7_tme_clusters.rds`, `cycif7_tme_and_macrophage_clusters.rds`: `frnn7.cmb` snapshots
  saved by notebook `5` -- the former has just `cluster_tme`, the latter also merges in `4`'s
  `cluster_mac` and is what notebook `6` actually loads.
- `cycif7_tme_colors_and_selection.rda`, `cycif7_tme_cluster_assignments.rds`: TME-clustering outputs, saved by notebook `5` and consumed by notebook `6`.
- `cycif7_umap_tumor_and_macrophage.rda`: saved UMAP embedding coordinates for tumor/macrophage cell populations.
- `cycif7_tcell_density_bayesian_fits.rda`: cached Bayesian (brms) model fits for the CD4T/CD8T/T_other hierarchical density model in notebook `6` -- required (the fitting code that would regenerate it is disabled).
- `cycif7_tumor_border_distances_by_sample.rds`, `cycif7_neighborhood_composition_by_sample.rds`: per-sample distance and neighborhood-composition objects, one list entry per sample.
- `pdl1_cc3_boundary_coords.rds`: manually defined PDL1/cCasp3 classification boundary, used to reconstruct `classify_pdl1_cc3()` in notebook `6`.
- `cell_type_definitions.rds`: cell-type hierarchy definitions for [EDF2C], loaded in notebook `5`.

## Setup

Every notebook's bootstrap chunk calls `initialize_project()` (`functions/bootstrap.R`), which
sources `functions/config_paths.R`. That file requires a single variable, `root_dir`, to
already exist -- every other path (`data_dir`, `fig_dir`, `source_dir`, and every
modality-specific subdirectory such as `cycif_obj_dir`, `io360_obj_dir`, etc.) is derived
from `root_dir` alone via `file.path(root_dir, ...)`. Nothing else needs to be set manually.

`root_dir` is expected to come from an R startup file, `.Rprofile`, placed at the repository
root (not tracked in this repo -- create your own). Its entire contents should be:

```r
root_dir <- "/absolute/path/to/this/repository"
```

RStudio/R sources `.Rprofile` automatically on startup when the working directory is the
repository root (e.g. opening `TALAVE_project.Rproj`). No `.Renviron` is required to run any
notebook -- that file, if present, only holds machine-specific settings (custom library
paths, memory limits) unrelated to this pipeline.

Notebooks `3`, `4`, `5`, `6`, and `7` also require the `CycifAnalyzeR` package, which is not on
CRAN/Bioconductor. Each of those notebooks' bootstrap chunks checks for it and installs it
from GitHub automatically if missing (`remotes::install_github("<org>/CycifAnalyzeR")`) --
the `remotes` package must already be installed for that to work.

## Session info

The analyses in this repository were run with the following R version and package versions
(`sessionInfo()`, captured from the environment used to produce the shipped figures):

```
R version 4.3.3 (2024-02-29)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS 26.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
  [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods
  [9] base

other attached packages:
  [1] deldir_2.0-4          dplyr_1.1.4           matrixStats_1.4.1
  [4] forcats_1.0.0         purrr_1.0.2           tidyr_1.3.1
  [7] tibble_3.2.1          data.table_1.16.0     RANN_2.6.2
  [10] igraph_2.0.3          leidenbase_0.1.32     FNN_1.1.4
  [13] rlang_1.1.4           reshape2_1.4.4        uwot_0.2.2
  [16] brms_2.22.0           Rcpp_1.0.13           broom.mixed_0.2.9.6
  [19] emmeans_1.10.6        lmerTest_3.1-3        lme4_1.1-35.5
  [22] Matrix_1.6-5          broom_1.0.7           org.Hs.eg.db_3.18.0
  [25] AnnotationDbi_1.64.1  IRanges_2.36.0        S4Vectors_0.40.2
  [28] Biobase_2.62.0        BiocGenerics_0.48.1   circlize_0.4.16
  [31] ComplexHeatmap_2.18.0 RColorBrewer_1.1-3    viridis_0.6.5
  [34] viridisLite_0.4.2     patchwork_1.2.0       ggrepel_0.9.6
  [37] ggplot2_3.5.1         cowplot_1.1.3         stringr_1.5.1
  [40] fgsea_1.28.0          MASS_7.3-60.0.1
```

## Running notebooks

Each notebook begins by sourcing `functions/bootstrap.R` and calling `initialize_project()`, then loading modality-specific helper scripts from `functions/` and a precomputed checkpoint from `data/`. Notebooks are intended to be run in numeric order (`1` -> `7`).

**Note on `3`-`6`:** `3` (bulk cell-type analysis) and `4` (UMAP analysis) each only need the
base `cs5`/`cs7`/`frnn7.cmb` checkpoint and have no dependency on each other or on `5`. `5`
(TME clustering)'s own analysis is likewise independent, but its final export step merges in
`4`'s `cluster_mac` output to produce the combined checkpoint that `6` loads. Run them in
order `3` -> `4` -> `5` -> `6` for a single top-to-bottom pass; each is independently
re-runnable from a clean R session via explicit `.rda`/`.rds` file checkpoints (no shared
live-session state is assumed).
