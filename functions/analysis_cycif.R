
prepare_sample_frequency_df <- function(freq_mat, pd) {
  t(freq_mat) %>%
    as.data.frame() %>%
    rownames_to_column("id") %>%
    left_join(pd, by = "id") %>%
    arrange(bt3) %>%
    group_by(Patient.ID) %>%
    mutate(
      max_bn = max(as.numeric(TimePoint)),
      text = ifelse(as.numeric(TimePoint) == max_bn, as.character(Patient.ID), "")
    ) %>%
    ungroup()
}

prepare_tumor_frequency_df <- function(freq_mat, pd) {
  freq_mat <- as.matrix(freq_mat)
  stopifnot(all(c("Tumor_panCK", "Tumor_neg") %in% rownames(freq_mat)))

  freq_df <- prepare_sample_frequency_df(freq_mat, pd)
  freq_df %>% mutate(Tumor = Tumor_panCK + Tumor_neg)
}

build_augmented_celltype_frequency <- function(cs) {
  ctf <- cellTypeFrequency(cs, simple = FALSE)
  rn.all <- rownames(ctf[[1]])
  tum.only <- ctf[[1]][grep("Tumor", rn.all), ]
  tum.all <- colSums(tum.only)
  tum.rel <- apply(tum.only, 2, function(x) x / sum(x))
  rownames(tum.rel) <- c("Tumor_panCK:Tumor_all", "Tumor_neg:Tumor_all")

  ctf[[1]] <- rbind(ctf[[1]], "Tumor_all:all" = tum.all)
  ctf <- c(ctf, Tumor = list(tum.rel))

  n.ctf <- do.call(rbind, ctf)
  rownames(n.ctf) <- sub(":", ".", rownames(n.ctf))
  n.ctf
}

compute_signed_logp <- function(df) {
  df %>%
    mutate(signed_logp = -log10(p.value) * sign(estimate)) %>%
    select(contrast, signed_logp)
}

compute_pvalue <- function(df) {
  df %>% select(contrast, p.value)
}

## BH correction that excludes NA from the family size (base p.adjust() otherwise counts
## NA positions toward n, inflating the correction), used per-column below since each
## column is one BRCA x timepoint contrast, corrected across its cell-type rows.
bh_ignore_na <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  out[ok] <- p.adjust(p[ok], method = "BH")
  out
}

## Refits the same model plot_ct() fits internally (lmer(y_var ~ bt3 + (1|Patient.ID)))
## and returns its residuals -- plot_ct() itself only returns list(plot, stat), not the
## model object, so this exists purely for the pooled normality diagnostic (mirrors
## get_gene_lmer_residuals() in functions/analysis_io360.R). Used across every cell-type/
## marker frequency variable analyzed via plot_ct() in this project, including the ones
## feeding the Fig.5A/EDF7A heatmaps (conds.bt3's underlying fits).
get_ct_freq_lmer_residuals <- function(df, y_var) {
  d <- df
  if (n_distinct(d$bt3) < 2 || n_distinct(d$Patient.ID) < 2) {
    return(NULL)
  }
  mod <- tryCatch(
    lmer(as.formula(paste(y_var, "~ bt3 + (1 | Patient.ID)")), data = d),
    error = function(e) NULL
  )
  if (is.null(mod)) return(NULL)
  as.numeric(residuals(mod))
}

add_cycif_umap_metadata <- function(df, pd, include_bor = FALSE) {
  join_cols <- c("id", "Patient.ID", "BRCA", "TimePoint", "bt3")
  if (include_bor) {
    join_cols <- c(join_cols, "BOR")
  }

  df %>%
    left_join(
      pd %>%
        rename(lspid = "sample") %>%
        mutate(
          batch = case_when(
            substr(lspid, 1, 2) == "13" ~ 1,
            substr(lspid, 1, 2) == "14" ~ 2,
            substr(lspid, 1, 2) == "22" ~ 3,
            TRUE ~ NA_real_
          ),
          batch = factor(as.numeric(batch))
        ) %>%
        select(all_of(join_cols), batch),
      by = c("sample" = "id")
    )
}

build_umap_embedding <- function(df, features, n_neighbors = 20, seed = 123, scale = TRUE) {
  mat <- df %>% select(all_of(features)) %>% as.matrix()
  set.seed(seed)
  umap <- uwot::umap(mat, scale = scale, n_neighbors = n_neighbors, verbose = TRUE)
  colnames(umap) <- c("umap_x", "umap_y")
  umap
}

bind_umap_to_tables <- function(df_list, umap, pd, include_bor = FALSE) {
  lapply(df_list, function(df) {
    df %>%
      cbind(umap) %>%
      add_cycif_umap_metadata(pd = pd, include_bor = include_bor)
  })
}

prepare_signed_logp_heatmap_inputs <- function(conds_bt3, sig_cut = -log10(0.05), cap = 3.2) {
  signed_logps_list <- lapply(conds_bt3, compute_signed_logp)
  signed_logps_df <- Reduce(function(x, y) merge(x, y, by = "contrast", all = TRUE), signed_logps_list)
  colnames(signed_logps_df)[-1] <- names(conds_bt3)

  signed_logps_matrix <- as.matrix(signed_logps_df[, -1])
  rownames(signed_logps_matrix) <- signed_logps_df$contrast
  colnames(signed_logps_matrix) <- sub("\\.", ":", colnames(signed_logps_matrix))

  cns <- colnames(signed_logps_matrix)
  removed <- cns[which(grepl("other", cns) & !grepl("T_other", cns))]
  cns <- cns[!cns %in% removed]

  signed_logps_matrix <- signed_logps_matrix[c(1, 2, 4, 7:9, 3, 5, 6), cns]
  signed_logps_matrix[7:9, ] <- -signed_logps_matrix[7:9, ]
  rownames(signed_logps_matrix)[7:9] <- c("WT_BS vs MUT_BS", "WT_BX2 vs MUT_BX2", "WT_BX3 vs MUT_BX3")

  ps1 <- t(signed_logps_matrix)
  div_by_all <- grepl(":all$", rownames(ps1))
  ps <- ps1[div_by_all, c(1, 2, 4, 5), drop = FALSE]
  rownames(ps) <- sub(":all", "", rownames(ps))
  colnames(ps) <- sub("(MUT_BS|WT_BS) vs ", "", colnames(ps))

  ## Parallel p-value matrix, put through the identical row-select/reorder/column-filter
  ## steps as signed_logps_matrix above (no sign flip needed -- p-values are unsigned).
  pvals_list <- lapply(conds_bt3, compute_pvalue)
  pvals_df <- Reduce(function(x, y) merge(x, y, by = "contrast", all = TRUE), pvals_list)
  colnames(pvals_df)[-1] <- names(conds_bt3)

  pvals_matrix <- as.matrix(pvals_df[, -1])
  rownames(pvals_matrix) <- pvals_df$contrast
  colnames(pvals_matrix) <- sub("\\.", ":", colnames(pvals_matrix))
  pvals_matrix <- pvals_matrix[c(1, 2, 4, 7:9, 3, 5, 6), cns]
  rownames(pvals_matrix)[7:9] <- c("WT_BS vs MUT_BS", "WT_BX2 vs MUT_BX2", "WT_BX3 vs MUT_BX3")

  pvals_t <- t(pvals_matrix)
  pvals_ps <- pvals_t[div_by_all, c(1, 2, 4, 5), drop = FALSE]
  rownames(pvals_ps) <- sub(":all", "", rownames(pvals_ps))
  colnames(pvals_ps) <- sub("(MUT_BS|WT_BS) vs ", "", colnames(pvals_ps))
  pvals_ps1 <- pvals_t

  ## Two-tier significance, same convention as IO360/DSP Protein: BH-adjusted per column
  ## (i.e. within each BRCA x timepoint contrast, across the cell-type rows tested against
  ## it) -- "fdr" if padj<sig_cut, "nominal" if only the raw p clears sig_cut, else "none".
  p_adj_ps <- apply(pvals_ps, 2, bh_ignore_na)
  dimnames(p_adj_ps) <- dimnames(pvals_ps)
  p_adj_ps1 <- apply(pvals_ps1, 2, bh_ignore_na)
  dimnames(p_adj_ps1) <- dimnames(pvals_ps1)

  make_tier <- function(pval_mat, padj_mat, cut) {
    tier <- matrix("none", nrow = nrow(pval_mat), ncol = ncol(pval_mat), dimnames = dimnames(pval_mat))
    tier[!is.na(pval_mat) & pval_mat < cut] <- "nominal"
    tier[!is.na(padj_mat) & padj_mat < cut] <- "fdr"
    tier
  }
  tier_ps <- make_tier(pvals_ps, p_adj_ps, 0.05)
  tier_ps1 <- make_tier(pvals_ps1, p_adj_ps1, 0.05)

  mat_val_filt <- ps
  mat_val_filt1 <- ps1
  mat_val_filt[mat_val_filt > cap] <- cap
  mat_val_filt[mat_val_filt < -cap] <- -cap
  mat_val_filt1[mat_val_filt1 > cap] <- cap
  mat_val_filt1[mat_val_filt1 < -cap] <- -cap

  rownames(mat_val_filt) <- gsub("\\.", "-", rownames(mat_val_filt))
  rownames(mat_val_filt1) <- gsub("\\.", "-", rownames(mat_val_filt1))
  rownames(tier_ps) <- gsub("\\.", "-", rownames(tier_ps))
  rownames(tier_ps1) <- gsub("\\.", "-", rownames(tier_ps1))
  tier_ps <- tier_ps[rownames(mat_val_filt), colnames(mat_val_filt)]
  tier_ps1 <- tier_ps1[rownames(mat_val_filt1), colnames(mat_val_filt1)]

  col_fun <- circlize::colorRamp2(c(-cap, 0, cap), c("#2166ac", "white", "#b2182b"))
  ann_colors <- list(
    BRCA = c(MUT = "#1b9e77", WT = "#7570b3", `MUT vs WT` = "#d95f02"),
    Time = c(BX2 = "#fc8d62", BX3 = "#8da0cb", `BX2 vs BX3` = "#66a61e")
  )
  col_annot <- data.frame(
    BRCA = c("MUT", "MUT", "WT", "WT"),
    Time = c("BX2", "BX3", "BX2", "BX3"),
    stringsAsFactors = FALSE,
    row.names = colnames(mat_val_filt)
  )
  col_annot1 <- data.frame(
    BRCA = c("MUT", "MUT", "MUT", "WT", "WT", "WT", "MUT vs WT", "MUT vs WT", "MUT vs WT"),
    Time = c("BX2", "BX3", "BX2 vs BX3", "BX2", "BX3", "BX2 vs BX3", "BX2", "BX3", "BX2 vs BX3"),
    stringsAsFactors = FALSE,
    row.names = colnames(mat_val_filt1)
  )

  list(
    mat_val_filt = mat_val_filt,
    mat_val_filt1 = mat_val_filt1,
    tier_ps = tier_ps,
    tier_ps1 = tier_ps1,
    col_fun = col_fun,
    ann_colors = ann_colors,
    col_annot = col_annot,
    col_annot1 = col_annot1
  )
}

build_signed_logp_heatmap <- function(mat, tier_mat, annotation_df, ann_colors, col_fun) {
  ha_top <- HeatmapAnnotation(
    df = annotation_df[, c("BRCA", "Time"), drop = FALSE],
    col = ann_colors[c("BRCA", "Time")],
    which = "column",
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10),
    simple_anno_size = unit(3.0, "mm")
  )

  Heatmap(
    mat,
    name = "Signed logP",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    clustering_distance_columns = "euclidean",
    row_order = 1:nrow(mat),
    show_row_names = TRUE,
    show_column_names = TRUE,
    top_annotation = ha_top,
    rect_gp = gpar(col = NA),
    cell_fun = function(j, i, x, y, w, h, fill) {
      ## Same two-tier convention as IO360/DSP Protein: nominal p<0.05 -> grey box (drawn
      ## first), FDR<0.05 -> black box on top (so it's never obscured by the grey layer).
      tier <- tier_mat[i, j]
      if (tier %in% c("nominal", "fdr")) {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(fill = NA, col = "grey50", lwd = 1.5))
      }
      if (tier == "fdr") {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(fill = NA, col = "black", lwd = 2))
      }
    }
  )
}

knn_impute_dt <- function(
  dt,
  k = 20,
  weighted = TRUE,
  rebuild_per_col = FALSE,
  fallback = c("median", "zero")
) {
  fallback <- match.arg(fallback)

  stopifnot(requireNamespace("RANN", quietly = TRUE))
  stopifnot(requireNamespace("matrixStats", quietly = TRUE))

  num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
  X <- as.matrix(dt[, ..num_cols])
  n <- nrow(X)
  p <- ncol(X)
  if (n < 2) return(dt)

  col_meds <- matrixStats::colMedians(X, na.rm = TRUE)
  X0 <- X
  for (j in seq_len(p)) {
    nas <- is.na(X0[, j])
    if (any(nas)) X0[nas, j] <- if (is.finite(col_meds[j])) col_meds[j] else 0
  }

  get_nn <- function(Xsrc, Xref = Xsrc, kk = k) {
    kk_eff <- min(kk + 1, n)
    nn <- RANN::nn2(Xref, Xsrc, k = kk_eff)
    idx <- nn$nn.idx
    d <- nn$nn.dists
    if (ncol(idx) > 1) {
      idx <- idx[, -1, drop = FALSE]
      d <- d[, -1, drop = FALSE]
    }
    list(idx = idx, d = d)
  }

  if (!rebuild_per_col) {
    nn_global <- get_nn(X0)
  }

  for (j in seq_len(p)) {
    rows_na <- which(is.na(X[, j]))
    if (!length(rows_na)) next

    if (rebuild_per_col) {
      feat <- X0[, -j, drop = FALSE]
      nn <- get_nn(feat)
    } else {
      nn <- nn_global
    }

    for (i in rows_na) {
      neigh <- nn$idx[i, ]
      dists <- nn$d[i, ]
      keep <- !is.na(X[neigh, j])

      if (!any(keep)) {
        X[i, j] <- if (fallback == "median" && is.finite(col_meds[j])) col_meds[j] else 0
        next
      }

      neigh <- neigh[keep]
      dists <- dists[keep]

      if (!weighted) {
        X[i, j] <- mean(X[neigh, j])
      } else {
        w <- 1 / pmax(dists, 1e-8)
        w <- w / sum(w)
        X[i, j] <- sum(w * X[neigh, j])
      }
    }
  }

  out <- copy(dt)
  out[, (num_cols) := as.data.frame(X)]
  out
}

make_square <- function(cx, cy, h) {
  coords <- matrix(
    c(
      cx - h, cy - h,
      cx + h, cy - h,
      cx + h, cy + h,
      cx - h, cy + h,
      cx - h, cy - h
    ),
    ncol = 2,
    byrow = TRUE
  )
  st_polygon(list(coords))
}

get_marker_means <- function(
  tme_cluster,
  marker,
  frnn_meta,
  frnn_exp,
  bt3_levels = NULL,
  collapse_bt = TRUE
) {
  cells <- frnn_meta %>%
    filter(cluster_tme == tme_cluster) %>%
    select(sample, cell_id, BRCA, TimePoint)

  exp_tbl <- cells %>%
    left_join(frnn_exp, by = c("sample", "cell_id"))

  if (!collapse_bt) {
    out <- exp_tbl %>%
      group_by(cell_types, BRCA, TimePoint) %>%
      summarize(
        marker_mean = mean(.data[[marker]], na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )

    if (!is.null(bt3_levels)) {
      out <- out %>%
        mutate(bt3 = factor(paste0(BRCA, "_", TimePoint), levels = bt3_levels))
    }

    out
  } else {
    exp_tbl %>%
      group_by(cell_types) %>%
      summarize(
        marker_mean = mean(.data[[marker]], na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
  }
}

build_marker_tme_matrix <- function(
  marker,
  tme_levels,
  frnn_meta,
  frnn_exp,
  bt3_levels = NULL,
  min_cells = 10,
  collapse_bt = TRUE
) {
  tme_list <- lapply(tme_levels, function(tme) {
    get_marker_means(
      tme_cluster = tme,
      marker = marker,
      frnn_meta = frnn_meta,
      frnn_exp = frnn_exp,
      bt3_levels = bt3_levels,
      collapse_bt = collapse_bt
    )
  })
  names(tme_list) <- tme_levels

  if (!collapse_bt) {
    tme_collapsed <- bind_rows(lapply(names(tme_list), function(tme) {
      tme_list[[tme]] %>%
        filter(n >= min_cells) %>%
        group_by(cell_types) %>%
        summarize(
          marker_mean = mean(marker_mean, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(TME = tme)
    }))
  } else {
    tme_collapsed <- bind_rows(lapply(names(tme_list), function(tme) {
      tme_list[[tme]] %>%
        filter(n >= min_cells) %>%
        mutate(TME = tme)
    }))
  }

  mat_tme <- tme_collapsed %>%
    select(TME, cell_types, marker_mean) %>%
    pivot_wider(names_from = cell_types, values_from = marker_mean) %>%
    column_to_rownames("TME") %>%
    as.matrix()

  cts_order <- intersect(get_cycif_uniq_cts(), colnames(mat_tme))
  mat_tme <- mat_tme[, cts_order, drop = FALSE]

  mat_tme[tme_levels, , drop = FALSE]
}

compute_cluster_stats <- function(df, ab1, clusters = c("TME_5", "TME_11", "TME_12", "TME_16")) {
  df_summary <- df %>%
    group_by(Patient.ID, cluster) %>%
    summarize(ab_mean = mean(.data[[ab1]], na.rm = TRUE), .groups = "drop")

  df_wide <- df_summary %>%
    pivot_wider(names_from = cluster, values_from = ab_mean)

  results <- lapply(clusters, function(cl) {
    if (!("others" %in% colnames(df_wide)) || !(cl %in% colnames(df_wide))) {
      return(data.frame(cluster = cl, p = NA, dir = NA))
    }

    test <- wilcox.test(df_wide[[cl]], df_wide[["others"]], paired = TRUE, exact = FALSE)
    diff <- median(df_wide[[cl]] - df_wide[["others"]], na.rm = TRUE)

    data.frame(
      cluster = cl,
      p = test$p.value,
      dir = ifelse(diff > 0, "up", "down")
    )
  })

  do.call(rbind, results) %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      sig = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ ""
      ),
      sig_col = case_when(
        dir == "up" ~ "#ff0000",
        dir == "down" ~ "#0072bc",
        TRUE ~ "black"
      )
    )
}

prepare_marker_positive_frequency_tables <- function(cts, exp_th, exp_log, pd, ab, this_cts) {
  cts <- cts %>%
    mutate(cell_types = fct_collapse(cell_types, Tumor = c("Tumor_panCK", "Tumor_neg")))

  if (any(grepl("Tumor", this_cts))) {
    this_cts <- unique(sub("Tumor_.+", "Tumor", this_cts))
    this_cts <- this_cts[!this_cts %in% c("Immune_other", "all_other")]
  }

  tmp <- cts %>%
    mutate(ab = exp_th[[ab]]) %>%
    mutate(is.pos = ab > 0.5) %>%
    mutate(sample = factor(sample, levels = as.character(pd$id))) %>%
    filter(cell_types %in% this_cts) %>%
    group_by(sample, cell_types) %>%
    summarize(pos_freq = mean(is.pos, na.rm = TRUE), .groups = "drop") %>%
    ungroup() %>%
    spread(cell_types, pos_freq)

  tmp1 <- cts %>%
    mutate(ab_log = exp_log[[ab]]) %>%
    mutate(sample = factor(sample, levels = as.character(pd$id))) %>%
    filter(cell_types %in% this_cts) %>%
    group_by(sample, cell_types) %>%
    summarize(mean_exp = mean(ab_log, na.rm = TRUE), .groups = "drop") %>%
    ungroup() %>%
    spread(cell_types, mean_exp)

  mat_fr <- as.matrix(tmp[-1])
  rownames(mat_fr) <- tmp$sample

  list(tmp = tmp, tmp1 = tmp1, mat_fr = mat_fr, selected_cts = this_cts)
}

prepare_density_coexpression_df <- function(df, x_var = "gH2AX", y_var = "pTBK1", palette = plasma(100)) {
  expr_mat <- df
  coex_mat <- cor(expr_mat, method = "spearman", use = "pairwise.complete.obs")
  dens <- kde2d(expr_mat[[x_var]], expr_mat[[y_var]], n = 200)
  ix <- findInterval(expr_mat[[x_var]], dens$x)
  iy <- findInterval(expr_mat[[y_var]], dens$y)
  d <- dens$z[cbind(ix, iy)]
  d_scaled <- log1p(d^0.5)
  expr_mat$density <- palette[as.numeric(cut(d_scaled, breaks = 100))]

  spearman <- cor.test(expr_mat[[x_var]], expr_mat[[y_var]], method = "spearman")
  pearson <- cor.test(expr_mat[[x_var]], expr_mat[[y_var]], method = "pearson")

  list(
    expr_mat = expr_mat,
    coex_mat = coex_mat,
    rho = unname(spearman$estimate),
    r = unname(pearson$estimate)
  )
}

prepare_sample_marker_means_df <- function(df, markers, group_cols = c("sample", "bt3", "Patient.ID", "TimePoint", "BOR")) {
  summary_names <- paste0(markers, "_mean")

  df %>%
    select(all_of(c(markers, group_cols))) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      across(all_of(markers), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"),
      .groups = "drop"
    ) %>%
    select(all_of(c(group_cols, summary_names)))
}

prepare_double_positive_frequency_df <- function(cts, exp_th, pd, ab1, ab2, selected_cts = "Tumor") {
  cts %>%
    mutate(cell_types = fct_collapse(cell_types, Tumor = c("Tumor_panCK", "Tumor_neg"))) %>%
    mutate(ab1 = exp_th[[ab1]]) %>%
    mutate(ab2 = exp_th[[ab2]]) %>%
    mutate(is.pos = ab1 > 0.5 & ab2 > 0.5) %>%
    mutate(sample = factor(sample, levels = as.character(pd$id))) %>%
    filter(cell_types %in% selected_cts) %>%
    group_by(sample, cell_types) %>%
    summarize(pos_frac = mean(is.pos, na.rm = TRUE), .groups = "drop") %>%
    left_join(pd %>% select(id, BRCA, TimePoint, PFS_days, Patient.ID, BOR, bt3), by = c("sample" = "id"))
}
