
export_cycif_slide_celltype_plots <- function(cs5, cs7, pd5, ct_cols, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pd5 <- pd5 %>%
    dplyr::mutate(ttl = paste(id, BRCA, BOR, paste0(PFS_days, "d"), paste0("LSP", sample), sep = ","))
  ttls5 <- pd5$ttl
  pts <- pd5$id
  names(pts) <- names(ttls5) <- pd5$id

  for (smpl in names(cs5)) {
    cat(smpl, "..\n")

    png(file.path(out_dir, paste0(pts[[smpl]], "_cy5.png")), width = 700, height = 700)
    slidePlot(
      cs5[[smpl]],
      plot_type = "cell_type",
      ncells = 1e5,
      legend = TRUE,
      mar = c(3, 3, 3, 10),
      uniq.cols = ct_cols,
      ttl = paste0(ttls5[smpl], "(cycle 5)")
    )
    dev.off()

    if (smpl %in% names(cs7)) {
      png(file.path(out_dir, paste0(pts[[smpl]], "_cy7.png")), width = 700, height = 700)
      slidePlot(
        cs7[[smpl]],
        plot_type = "cell_type",
        ncells = 1e5,
        legend = TRUE,
        mar = c(3, 3, 3, 10),
        uniq.cols = ct_cols,
        ttl = paste0(ttls5[smpl], "(cycle 7)")
      )
      dev.off()
    }
  }
}

