umap_styling_num <- function() {
  list(
    theme(
      legend.title = element_text(size = 28, margin = margin(b = 10)),
      legend.text = element_text(size = 24),
      legend.key.size = unit(1.2, "cm"),
      legend.spacing.y = unit(1, "cm"),
      legend.key.height = unit(1.2, "cm"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  )
}

common_scale <- function(max_limit = 15) {
  scale_color_distiller(
    palette = "YlGnBu",
    limits = c(0, max_limit),
    direction = 1,
    oob = scales::squish,
    name = "Density"
  )
}

umap_styling <- function() {
  list(
    theme(
      plot.title = element_text(size = 60, hjust = 0.5),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 24),
      legend.key.size = unit(1.2, "cm"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  )
}

umap_styling_1 <- function(
  title_size = 60,
  legend_title_margin = TRUE,
  legend_key_height = 1.2,
  legend_title_size = 39
) {
  legend_margin <- if (legend_title_margin) margin(b = 20) else margin()

  list(
    theme(
      plot.title = element_text(size = title_size, hjust = 0.5),
      legend.text = element_text(size = 36),
      legend.title = element_text(size = legend_title_size, margin = legend_margin),
      legend.key.size = unit(1.2, "cm"),
      legend.spacing.y = unit(1, "cm"),
      legend.key.height = unit(legend_key_height, "cm"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  )
}

umap_base <- function(data, mapping, point_size = 0.3) {
  lims <- list(
    xlim = range(data$umap_x, na.rm = TRUE),
    ylim = range(data$umap_y, na.rm = TRUE)
  )

  ggplot(data, mapping) +
    geom_point(size = point_size) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    xlim(lims$xlim) +
    ylim(lims$ylim) +
    umap_styling() +
    guides(color = guide_legend(override.aes = list(size = 6)))
}

umap_base_num <- function(data, mapping, point_size = 0.3) {
  lims <- list(
    xlim = range(data$umap_x, na.rm = TRUE),
    ylim = range(data$umap_y, na.rm = TRUE)
  )

  ggplot(data, mapping) +
    geom_point(size = point_size) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    xlim(lims$xlim) +
    ylim(lims$ylim) +
    umap_styling()
}

save_legend_pdf <- function(plot_with_legend, file, width = 4, height = 6) {
  legend <- get_legend(plot_with_legend)
  pdf(file, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(legend)
  dev.off()
}

blend_with_white_palette <- function(colors, blend_ratio = 0.5) {
  sapply(colors, function(col) {
    ramp <- colorRampPalette(c(col, "white"))
    n <- 100
    ramp(n)[round(blend_ratio * (n - 1)) + 1]
  }, USE.NAMES = TRUE)
}

plot_marker_tme_heatmap <- function(
  marker,
  tme_levels,
  frnn_meta,
  frnn_exp,
  bt3_levels = NULL,
  min_cells = 10,
  collapse_bt = TRUE,
  cap_quantile = 0.95,
  file = NULL,
  width = 6,
  height = 6,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  na_col = "grey90"
) {
  mat_tme <- build_marker_tme_matrix(
    marker = marker,
    tme_levels = tme_levels,
    frnn_meta = frnn_meta,
    frnn_exp = frnn_exp,
    bt3_levels = bt3_levels,
    min_cells = min_cells,
    collapse_bt = collapse_bt
  )

  mat_tme <- mat_tme[
    rowSums(!is.na(mat_tme)) > 0,
    colSums(!is.na(mat_tme)) > 0,
    drop = FALSE
  ]

  thr <- quantile(abs(as.vector(mat_tme)), probs = cap_quantile, na.rm = TRUE)
  if (!is.finite(thr) || thr == 0) {
    thr <- max(abs(mat_tme), na.rm = TRUE)
  }
  if (!is.finite(thr) || thr == 0) {
    thr <- 1
  }

  mat_plot <- pmin(pmax(mat_tme, -thr), thr)

  ht <- Heatmap(
    mat_plot,
    name = marker,
    col = colorRamp2(c(-thr, 0, thr), c("blue", "white", "red")),
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    na_col = na_col
  )

  if (!is.null(file)) {
    pdf(file, width = width, height = height)
    ComplexHeatmap::draw(ht)
    dev.off()
  }

  invisible(list(mat = mat_tme, mat_plot = mat_plot, ht = ht, thr = thr))
}

plot_scatter_locator_base <- function(
  df, xcol, ycol,
  compute_density = TRUE, dens_col = "density", k = 20,
  n_clicks = 10, number_clicks = TRUE,
  main_title = NULL, xlab = NULL, ylab = NULL,
  add_colorbar = TRUE,
  save_file = NULL, save_width = 1600, save_height = 1400, save_res = 200
) {
  stopifnot(xcol %in% names(df), ycol %in% names(df))

  ok <- is.finite(df[[xcol]]) & is.finite(df[[ycol]])
  df <- df[ok, , drop = FALSE]

  if (compute_density) {
    if (!(dens_col %in% names(df))) {
      coords <- cbind(df[[xcol]], df[[ycol]])
      nn <- FNN::get.knn(coords, k = min(k, nrow(coords) - 1L))
      dbar <- rowMeans(nn$nn.dist, na.rm = TRUE)
      dens <- 1 / pmax(dbar, 1e-8)
    } else {
      coords <- cbind(df[[xcol]], df[[ycol]])
      nn <- FNN::get.knn(coords, k = min(k, nrow(coords) - 1L))
      idx_mat <- cbind(seq_len(nrow(df)), nn$nn.index)
      dens <- apply(idx_mat, 1, function(ix) mean(df[[dens_col]][ix], na.rm = TRUE))
    }
  } else {
    stopifnot(dens_col %in% names(df))
    dens <- df[[dens_col]]
  }

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Please install.packages('RColorBrewer')")
  }
  pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))
  cols <- pal(256)
  rng <- range(dens[is.finite(dens)], na.rm = TRUE)
  if (diff(rng) == 0) rng <- rng + c(-0.5, 0.5)
  col_idx <- 1 + floor((dens - rng[1]) / diff(rng) * 255)
  col_idx[!is.finite(col_idx)] <- 1
  pt_col <- cols[pmin(pmax(col_idx, 1), 256)]

  if (is.null(main_title)) main_title <- ""
  if (is.null(xlab)) xlab <- paste0(xcol, " (log)")
  if (is.null(ylab)) ylab <- paste0(ycol, " (log)")

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(4, 4, 3, if (add_colorbar) 5 else 2) + 0.1)

  plot(df[[xcol]], df[[ycol]],
       xlab = xlab, ylab = ylab, main = main_title,
       pch = 16, cex = 0.7, col = pt_col)
  grid(col = "grey85")

  if (add_colorbar && requireNamespace("fields", quietly = TRUE)) {
    fields::image.plot(
      legend.only = TRUE, zlim = rng, col = cols,
      legend.args = list(text = "mean\ndensity", side = 4, line = 2.5, cex = 0.9)
    )
  }

  cat(sprintf("Click up to %d points; press ESC/Stop to finish early.\n", n_clicks))
  pts <- locator(n = n_clicks)

  if (!is.null(pts) && length(pts$x)) {
    points(pts$x, pts$y, pch = 1, cex = 1.3, lwd = 2, col = "red")
    if (number_clicks) {
      text(pts$x, pts$y, labels = seq_along(pts$x), col = "red", pos = 3, cex = 0.9)
    }
  } else {
    pts <- list(x = numeric(0), y = numeric(0))
  }

  if (!is.null(save_file)) {
    png(save_file, width = save_width, height = save_height, res = save_res)
    par(mar = c(4, 4, 3, if (add_colorbar) 5 else 2) + 0.1)
    plot(df[[xcol]], df[[ycol]],
         xlab = xlab, ylab = ylab, main = main_title,
         pch = 16, cex = 0.7, col = pt_col)
    grid(col = "grey85")
    if (add_colorbar && requireNamespace("fields", quietly = TRUE)) {
      fields::image.plot(
        legend.only = TRUE, zlim = rng, col = cols,
        legend.args = list(text = "mean\ndensity", side = 4, line = 2.5, cex = 0.9)
      )
    }
    if (length(pts$x)) {
      points(pts$x, pts$y, pch = 1, cex = 1.3, lwd = 2, col = "red")
      if (number_clicks) {
        text(pts$x, pts$y, labels = seq_along(pts$x), col = "red", pos = 3, cex = 0.9)
      }
    }
    dev.off()
  }

  invisible(list(
    clicked = data.frame(x = pts$x, y = pts$y),
    range_density = rng
  ))
}

split_by_clicks <- function(df, clicks, xcol, ycol,
                            idx_v = 2, idx_h = 10, idx_curve = 2:7,
                            draw = TRUE) {
  stopifnot(all(c("x", "y") %in% names(clicks)))
  x0 <- clicks$x[idx_v]
  y0 <- clicks$y[idx_h]

  x <- df[[xcol]]
  y <- df[[ycol]]
  quadrant <- ifelse(x <= x0 & y <= y0, "BL",
              ifelse(x > x0 & y <= y0, "BR",
              ifelse(x <= x0 & y > y0, "TL", "TR")))
  quadrant <- factor(quadrant, levels = c("BL", "BR", "TL", "TR"))

  curve_pts <- clicks[idx_curve, , drop = FALSE]
  ord <- order(curve_pts$x, curve_pts$y)
  curve_pts <- curve_pts[ord, , drop = FALSE]
  dupx <- duplicated(curve_pts$x)
  if (any(dupx)) curve_pts <- curve_pts[!dupx, , drop = FALSE]

  y_on_curve <- function(xq) {
    xi <- curve_pts$x
    yi <- curve_pts$y
    minx <- min(xi)
    maxx <- max(xi)
    xq2 <- pmin(pmax(xq, minx), maxx)
    i <- findInterval(xq2, xi, left.open = FALSE, rightmost.closed = TRUE)
    i <- pmin(pmax(i, 1), length(xi) - 1)
    x1 <- xi[i]
    x2 <- xi[i + 1]
    y1 <- yi[i]
    y2 <- yi[i + 1]
    y1 + (xq2 - x1) * (y2 - y1) / pmax(x2 - x1, 1e-12)
  }

  is_tr <- quadrant == "TR"
  tr_side <- rep(NA_character_, length(x))
  if (nrow(curve_pts) >= 2 && any(is_tr)) {
    y_curve <- y_on_curve(x[is_tr])
    tr_side[is_tr] <- ifelse(y[is_tr] > y_curve, "above", "below")
  }
  tr_side <- factor(tr_side, levels = c("below", "above"))

  if (draw) {
    abline(v = x0, h = y0, col = "grey30", lwd = 2, lty = 2)
    lines(curve_pts$x, curve_pts$y, col = "black", lwd = 2)
    points(clicks$x[2:7], clicks$y[2:7], pch = 19, col = "black")
    text(clicks$x[2], clicks$y[2], "#2", pos = 3)
    text(clicks$x[10], clicks$y[10], "#10", pos = 3)
  }

  out <- df
  out$quadrant <- quadrant
  out$tr_side <- tr_side
  out
}

make_viridis_legend <- function(limits = c(0, 5),
                                name = "Density",
                                file = NULL,
                                width = 4,
                                height = 6) {
  df_dummy <- data.frame(
    x = 1:100,
    y = 1,
    density = seq(limits[1], limits[2], length.out = 100)
  )

  p <- ggplot(df_dummy, aes(x = x, y = y, color = density)) +
    geom_point() +
    scale_color_viridis_c(
      option = "viridis",
      limits = limits,
      name = name,
      oob = scales::squish
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 16),
      legend.key.height = unit(2, "cm"),
      legend.key.width = unit(0.6, "cm")
    )

  legend <- cowplot::get_legend(p)

  if (!is.null(file)) {
    pdf(file, width = width, height = height)
    grid::grid.newpage()
    grid::grid.draw(legend)
    dev.off()
  }

  legend
}
