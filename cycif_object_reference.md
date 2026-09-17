# CyCIF object reference (notebooks 3-6)

Reference glossary for the key R objects passed between the four CyCIF notebooks
(`3_cycif_bulk_celltype_analysis.Rmd`, `4_cycif_umap_analysis.Rmd`,
`5_cycif_spatial_frnn_tme.Rmd`, `6_cycif_spatial_continued.Rmd`). Not part of the run
sequence -- no code chunks, nothing to execute. Intended as a lookup when an object name in
one of those notebooks isn't self-explanatory from its column names or usage.

`frnn7.cmb`, `frnn7.exp`, `frnn7.rcn`, and `is.sele.7` are loaded, not created, by any
notebook in this repo -- they come from `cycif_spatial_frnn_input_data.rda` /
`cycif7_neighborhood_profiles.rda`, checkpoints produced upstream of this Zenodo supplement
(see `README.md`, "Shared inputs and metadata"). Their descriptions below are inferred from
how they're used in notebooks `4`-`6`, not from their generating code.

## Core per-cell tables

**`cs5`, `cs7`**: the 5-cycle and 7-cycle `CycifStack` objects (from the `CycifAnalyzeR`
package), loaded from `cs5.rds`/`cs7.rds`. The underlying single-cell CyCIF data (marker
expression, coordinates, cell-type calls) for every sample; almost everything else in this
reference is either derived from these or joined against them.

**`pd5`, `pd7`**: `pData(cs5)` / `pData(cs7)` -- one row per sample, clinical/experimental
metadata (BRCA status, TimePoint, Patient.ID, etc).

**`cts5`, `cts7`, `cts7d`**: `cell_types(cs5)` / `cell_types(cs7)` -- one row per cell, the
called cell-type label (`cell_types` column) plus supporting per-cell classification fields.
`cts7d` is `cts7` with `NA`-cell-type rows dropped.

**`frnn7.cmb`** (loaded from `cycif_spatial_frnn_input_data.rda`; re-saved with additional
columns at successive pipeline stages -- see "Checkpoint chain" below): the main combined
per-7-cycle-cell table. One row per cell. Columns include `sample`, `cell_types`,
marker `.log` expression columns (e.g. `PDL1.log`, `cCaspase3.log`, `gH2AX.log`,
`pTBK1.log`), and the 11 `uniq.cts` neighborhood-composition columns (`Tumor_panCK`,
`Tumor_neg`, `Fibro`, `CD8T`, `CD4T`, `T_other`, `Mac_CD68`, `Mac_CD163`, `Mac_CD68_CD163`,
`Immune_other`, `all_other`) embedded directly (the same composition data also lives
standalone in `frnn7.rcn`). Downstream notebooks add `cluster_tme`, `cluster_mac`, `dist`,
`Tumor` (= `Tumor_panCK + Tumor_neg`), `predicted_prob_stroma`, `is_pred_stroma`, and (in
notebook `6`) per-TME `mac_neighbor_count_TME_*`/`tme_neighbor_label` columns. Row order is
assumed aligned across the whole checkpoint chain and with `cts7`/`cs7`/`lst.dists7`/
`lst.frnn7`'s per-sample ordering -- notebook `6` itself flags this as an assumption worth
verifying before trusting ROI-level output that depends on it.

**`frnn7.rcn`** ("Recurrent Cellular Neighborhood" composition): one row per cell, same 11
`uniq.cts` columns as embedded in `frnn7.cmb` -- a raw per-cell-type neighbor *count* within
some spatial radius (not a fraction; `rowSums(frnn7.rcn) >= 10` is used as an "enough
neighbors sampled" QC threshold in several places, e.g. `is.ge10.7 <- rowSums(frnn7.rcn) >=
10`). This is the standalone version of the composition data that's also embedded directly in
`frnn7.cmb`.

**`frnn7.exp`**: the paired neighborhood mean-expression profile, one row per cell. Used as
`mat.exp <- frnn7.exp` and concatenated onto the normalized RCN matrix
(`cbind(mat_row_norm, mat.exp %>% select(-1))`) to build the TME UMAP/Leiden feature space in
notebook `5`. The `select(-1)` implies its first column is an id/index rather than a feature
column, but this isn't stated directly anywhere -- treat as inferred, not confirmed.

**`lst.dists7`** (loaded from `cycif7_tumor_border_distances_by_sample.rds`): a named list,
one element per sample. Each element has `$dist$dist` (numeric vector, per-cell signed
distance to the tumor border -- used to build `frnn7.cmb$dist`) and `$borders` (a list of
border polygon geometries per sample, later intersected with `sf` squares for the ROI
figures).

**`lst.frnn7`** (loaded from `cycif7_neighborhood_composition_by_sample.rds`): a named list,
one element per sample, each an object with slot `@nn@id` -- a list of neighbor-cell-index
vectors, one per cell, giving the indices (into that sample's cells) of its spatial
neighbors. This is the fixed-radius nearest-neighbor (FRNN) search result that `frnn7.rcn`/
`frnn7.exp` are presumably summarized from upstream.

## Selection / QC vectors

All are boolean, length = the population they subset.

**`is.sele.7`**: loaded (not created) from `cycif_spatial_frnn_input_data.rda`; a baseline
per-cell QC/eligibility flag used everywhere as a starting filter (`filter(is.sele.7)`,
`n>=10 & is.sele.7`). Its exact construction criterion isn't shown in any of notebooks
`3`-`6` -- inferred to be some upstream inclusion flag (e.g. not-out-of-ROI), not confirmed.

**`is.sele.tme`** (notebook `5`): `n>=10 & is.sele.7`, where `n` is each cell's total
neighbor count computed from `frnn7.cmb`'s embedded RCN columns. Cells with enough sampled
neighbors *and* passing baseline QC -- the subsample (roughly 1-3% of cells per cell type) on
which the actual TME Leiden clustering runs; `cluster_tme` is assigned directly only for
these cells, then propagated to everyone else (see below).

**`is.sele.tum`** (notebook `4`): a per-sample-balanced down-sample of Tumor cells
(`Tumor_panCK`+`Tumor_neg`), capped at the 20th percentile of per-sample tumor-cell counts so
no single high-cellularity sample dominates the tumor UMAP.

**`is.sele.mac`** (notebook `4`): the analogous balanced down-sample restricted to
macrophage subtypes (`Mac_CD68`, `Mac_CD68_CD163`, `Mac_CD163`).

**`is.sele.full`** (notebook `4`): `is.sele.mac` (defined only over macrophage cells)
scattered back into a full-length boolean vector over every cell with a called type -- i.e.
`is.sele.mac` re-indexed into `frnn7.cmb`'s row space, used to assign
`frnn7.cmb$cluster_mac[is.sele.full] <- ...`.

## UMAP embedding tables

**`dt.tme.1`** (notebook `5`; saved to `cycif7_tme_cluster_assignments.rds`): one row per
`is.sele.tme` cell. `frnn7.cmb`'s columns for that row, plus `umap_x`/`umap_y` (TME UMAP
coordinates), metadata joined from `pd7`, and (after Leiden clustering) `cluster` -- a
`"TME_1"`..`"TME_17"` factor.

**`dt.tum.1`** (notebook `4`): one row per `is.sele.tum` tumor cell. `umap_x`/`umap_y`
(tumor-cell UMAP), `sample`, `cell_types`, marker expression, and metadata joined from `pd7`.

**`dt.mac.1`, `dt.mac.2`** (notebook `4`): two parallel tables, one row per `is.sele.mac`
macrophage cell, built the same way as `dt.tum.1`. By analogy with the earlier `exp.1`/
`exp.2` pattern in this notebook (continuous log expression vs. `>0.5`-thresholded binary
expression), `dt.mac.2` is likely the thresholded variant -- this isn't stated directly in a
comment, so treat as plausible rather than confirmed. After clustering, both get a `cluster`
factor column, `"Mac_1"`..`"Mac_10"`.

## Cluster-label columns (merged onto `frnn7.cmb`)

**`cluster_tme`** (notebook `5`): `"TME_1"`..`"TME_17"` (17 Leiden clusters,
`resolution_parameter = 2e-3`). Assigned directly only for `is.sele.tme` cells
(`frnn7.cmb$cluster_tme[is.sele.tme] <- as.character(dt.tme.1$cluster)`), then **explicitly
propagated** to the rest of the population via 20-NN majority vote in the same scaled
RCN+expression feature space (`mat.cmb.all`, scaled using the `is.sele.tme` population's own
center/scale) -- see notebook `5`'s "Propagate cluster_tme to the Full Cell Population"
section. Coverage after propagation is printed; any remaining `NA`s are cells where
`stats::complete.cases(mat.cmb.all)` is `FALSE` (missing expression data).

**`cluster_mac`** (notebook `4`): `"Mac_1"`..`"Mac_10"` (Leiden clusters on macrophage marker
expression, `resolution_parameter = 5e-4`). Assigned only for `is.sele.full` cells
(`frnn7.cmb$cluster_mac[is.sele.full] <- paste0("Mac_",cluster_labels)`) -- **there is no
propagation step for `cluster_mac`** anywhere in notebooks `4`, `5`, or `6`, unlike
`cluster_tme`. It is `NA` for every macrophage not in the down-sampled `is.sele.mac`/
`is.sele.full` set (a small minority of all macrophages). The macrophage-cluster-enrichment
analyses in notebook `6` (Fig 5D/5I/5J, EDF8G) work around this by joining `dt.mac.1$cluster`
directly onto a fresh macrophage-filtered subset, rather than relying on `frnn7.cmb$cluster_mac`.

## Marker taxonomies and palettes

**`uniq.cts`** (`get_cycif_uniq_cts()`): the fixed 11-element cell-type taxonomy used for RCN
composition -- `Tumor_panCK`, `Tumor_neg`, `Fibro`, `CD8T`, `CD4T`, `T_other`, `Mac_CD68`,
`Mac_CD163`, `Mac_CD68_CD163`, `Immune_other`, `all_other`.

**`cst7.abs`** (`colnames(cs7@cell_types$default@cell_state_def)`): the antibody/marker names
in the 7-cycle cell-state classification scheme (e.g. `PD1`, `PDL1`, `cCaspase3`, `gH2AX`,
`pTBK1`, `BCLXL`, `MCL1`, `pAKT`, `pERK`). Used to iterate marker-expression UMAP overlays and
heatmaps.

**`ct.cols7`** (`get_cycif_ct_colors()`): named hex-color vector keyed by the `uniq.cts`
labels.

**`ct.cols17`** (notebook `5`; saved in `cycif7_tme_colors_and_selection.rda`): named
hex-color vector keyed by `"TME_1"`..`"TME_17"`, built from a Spectral-palette permutation
chosen for visual distinctiveness between adjacent cluster numbers.

## Manual marker-classification boundaries

**`bound_coords`** (notebook `6`; loaded from `pdl1_cc3_boundary_coords.rds`): a manually
defined PDL1/cCaspase3 decision boundary, set by interactive clicking (`plot_scatter_locator_base()`/`split_by_clicks()` in `functions/plotting_cycif.R`) and reused as a fixed
constant here rather than re-clicked. A `list` with `$vline$x` (a single numeric PDL1
threshold) and `$curve` (a data.frame of `x,y` points defining a piecewise-linear
cCaspase3-vs-PDL1 boundary, interpolated via `findInterval()` in `classify_pdl1_cc3()`).

**`boundaries2`** (notebook `6`): `list(v=c(0.06300655), h=c(-0.09366919))` -- a hardcoded
gH2AX/pTBK1 threshold pair defining a simple quadrant split (not a curve), used in
`classify_gh2ax_ptbk1()`.

**`classification`** (notebook `6`, via `classify_pdl1_cc3()`): factor with levels
`"PDL1+CC3+"`, `"PDL1+CC3-"`, `"PDL1-"`. Assigned onto `frnn7.cmb.tum`/`frnn7.cmb.mac` (the
tumor/macrophage cell-type-filtered, `is.sele.tum`/`is.sele.mac`-restricted subsets built in
notebook `6`), not the full `frnn7.cmb`.

**`classification_ptbk`** (notebook `6`, via `classify_gh2ax_ptbk1()`): factor with 4 levels
(`gH2AX-/pTBK1-`, `gH2AX+/pTBK1-`, `gH2AX-/pTBK1+`, `gH2AX+/pTBK1+`), assigned the same way.

## Predicted tumor/stroma classification

**`is_pred_stroma` / `predicted_prob_stroma`** (notebook `6`): a logistic regression
predicting stroma membership from local cell-type composition (`frnn7.cmb`'s own `uniq.cts`
columns), trained on a size- and distance-stratified subsample of baseline-sample cells, then
applied to every cell in `frnn7.cmb`. `predicted_prob_stroma` is the 0-1 predicted
probability (`glm(is_stroma ~ Tumor + Fibro + CD8T + CD4T + T_other + Mac_CD68 + Mac_CD163 +
Mac_CD68_CD163 + Immune_other, family=binomial)`); `is_pred_stroma` is the thresholded binary
label (`predicted_prob_stroma > 0.5`). The training label itself (`is_stroma`, not saved onto
`frnn7.cmb`) comes from `dist > 0.5` on baseline-sample cells passing the neighbor-count QC,
stratified per sample to balance near-border (`dist` in `[-250, 250]`) vs. far cells.

## Checkpoint chain (which notebook produces what)

1. `cycif_spatial_frnn_input_data.rda` (upstream, not produced here): `frnn7.cmb`,
   `frnn7.exp`, `frnn7.rcn`, `is.sele.7`.
2. Notebook `4` adds `cluster_mac` (down-sampled subset only) -> saves
   `cycif7_macrophage_clusters.rds`.
3. Notebook `5` adds `cluster_tme` (fully propagated) -> saves `cycif7_tme_clusters.rds`,
   then merges in `4`'s `cluster_mac` -> saves `cycif7_tme_and_macrophage_clusters.rds` (also
   saves `dt.tme.1` to `cycif7_tme_cluster_assignments.rds` and `ct.cols17`/`is.sele.tme` to
   `cycif7_tme_colors_and_selection.rda`).
4. Notebook `6` loads `cycif7_tme_and_macrophage_clusters.rds` as its working `frnn7.cmb`.
