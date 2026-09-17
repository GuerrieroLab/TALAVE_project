## ---- Clinical ----

load_clinical_metadata <- function(filename, input_dir = clinical_obj_dir) {
  readRDS(file.path(input_dir, filename))
}

## ---- IO360 ----

load_io360_expression_data <- function(
  filename = "io360_expression_input_data.rda",
  input_dir = io360_obj_dir
) {
  load(file.path(input_dir, filename), envir = parent.frame())
}

## ---- DSP Protein ----

load_dsp_protein_lookup <- function(input_dir = protein_obj_dir, filename = "protein_gene_lookup.rds") {
  readRDS(file.path(input_dir, filename))
}

## ---- CyCIF ----

load_cycif_cell_type_definitions <- function(
  filename = file.path(cycif_obj_dir, "cell_type_definitions.rds")
) {
  readRDS(filename)
}

get_cycif_ct_colors <- function() {
  c(
    Tumor_panCK = "#A3E635",
    Tumor_neg = "#4B0082",
    Fibro = "#F4D03F",
    CD8T = "#9E0142",
    CD4T = "#F88D51",
    T_other = "#D7BDE2",
    Mac_CD68 = "#40E0D0",
    Mac_CD163 = "#27AE60",
    Mac_CD68_CD163 = "#3498DB",
    Immune_other = "#EAF69E",
    all_other = "grey70"
  )
}

get_cycif_ct_colors_white_blend <- function() {
  blend_with_white_palette(get_cycif_ct_colors(), blend_ratio = 0.5)
}

get_cycif_uniq_cts <- function() {
  c("Tumor_panCK","Tumor_neg","Fibro","CD8T","CD4T","T_other",
    "Mac_CD68","Mac_CD163","Mac_CD68_CD163","Immune_other","all_other")
}
