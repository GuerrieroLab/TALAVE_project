initialize_project <- function(load_packages = TRUE) {
  source(file.path(root_dir, "functions", "config_paths.R"))

  if (isTRUE(load_packages)) {
    source(file.path(source_dir, "load_packages.R"))
  }

  invisible(
    list(
      root_dir = root_dir,
      source_dir = source_dir,
      data_dir = data_dir,
      fig_dir = fig_dir
    )
  )
}
