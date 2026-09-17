if(!exists("root_dir")) {
  stop("root_dir not set in .Rprofile")
}

source_dir <- file.path(root_dir, "functions")
data_dir <- file.path(root_dir, "data")
fig_dir <- file.path(root_dir, "figures")
obj_dir <- file.path(data_dir, "obj")

clinical_obj_dir <- file.path(obj_dir, "clinical")
cycif_obj_dir <- file.path(obj_dir, "cycif")
protein_obj_dir <- file.path(obj_dir, "protein")
io360_obj_dir <- file.path(obj_dir, "io360")

cycif_fig_dir <- file.path(fig_dir, "cycif_1")
protein_fig_dir <- file.path(fig_dir, "protein")
io360_fig_dir <- file.path(fig_dir, "io360")
integrated_fig_dir <- file.path(fig_dir, "integrated")
