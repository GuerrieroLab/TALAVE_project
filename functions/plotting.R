load_pt_colors <- function(filename = "patient_colors.rds", input_dir = clinical_obj_dir) {
  readRDS(file.path(input_dir, filename))
}

plot_ct <- function(df, y_var, y_label, title=NULL, pt.cols=pt.cols, ymin=NULL) {
  required_cols <- c("Patient.ID", "TimePoint", "bt3","BOR",y_var)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required column(s):", paste(missing_cols, collapse = ", ")))
  }

  df <- df %>%
    as.data.frame() %>%
    group_by(Patient.ID) %>%
    mutate(
      max_bn = max(as.numeric(TimePoint)),
      text = ifelse(as.numeric(TimePoint) == max_bn, as.character(Patient.ID), "")
    ) %>%
    ungroup() %>%
    mutate(bt3 = factor(bt3, levels = c("MUT_BS", "MUT_BX2", "MUT_BX3", "WT_BS", "WT_BX2", "WT_BX3")))

  if (is.null(ymin)) {
    ymin <- min(df[[y_var]], na.rm = TRUE)
  } else if (ymin==0){
    ymin <- min(0, min(df[[y_var]], na.rm = TRUE))
  }

  y_data_max <- max(df[[y_var]], na.rm = TRUE)
  y_range <- y_data_max - ymin
  y_mag <- y_data_max * 1.05

  conds <- rbind(c(1,2),c(2,3),c(1,3),c(4,5),c(5,6),c(4,6),c(1,4),c(2,5),c(3,6))
  paired <- c(T,T,T,T,T,T,F,F,F)

  y_offset_levels <- c(1,1,2,1,1,2,3,4,5)

  y_offset_levels <- y_offset_levels + 1

  ys <- (y_offset_levels * 0.1 + 0.9) * y_range + ymin

  df.conds <- data.frame(
    idx = seq_along(paired),
    cond = conds,
    paired = paired,
    p.value = rep(NA, 9),
    x = rowMeans(conds),
    y = ys,
    x1 = conds[,1] + 0.03,
    x2 = conds[,2] - 0.03,
    y1 = ys,
    y2 = ys
  )

  model <- lmer(as.formula(paste(y_var, "~ bt3 + (1 | Patient.ID)")), data = df)
  em <- emmeans(model, ~ bt3)
  contrast_matrix <- list(
    "MUT_BS vs MUT_BX2" = c(1, -1,  0,  0,  0,  0),
    "MUT_BX2 vs MUT_BX3" = c(0,  1, -1,  0,  0,  0),
    "MUT_BS vs MUT_BX3" = c(1,  0, -1,  0,  0,  0),
    "WT_BS vs WT_BX2" = c(0,  0,  0,  1, -1,  0),
    "WT_BX2 vs WT_BX3" = c(0,  0,  0,  0,  1, -1),
    "WT_BS vs WT_BX3" = c(0,  0,  0,  1,  0, -1),
    "MUT_BS vs WT_BS" = c(1,  0,  0, -1,  0,  0),
    "MUT_BX2 vs WT_BX2" = c(0,  1,  0,  0, -1,  0),
    "MUT_BX3 vs WT_BX3" = c(0,  0,  1,  0,  0, -1)
  )
  contrast_results <- contrast(em, contrast_matrix, adjust = "none") %>%
    as.data.frame() %>%
    mutate(estimate = -estimate) %>%
    cbind(df.conds %>% select(x, y, x1, x2, y1, y2))

  ymax <- max(max(df[[y_var]], na.rm = TRUE), max(ys)) * 1.05

  p <- ggplot(df, aes(x = bt3, y = .data[[y_var]], fill = Patient.ID)) +
    geom_path(aes(group = Patient.ID, linetype = BOR), color = "grey60", size = 0.4) +
    geom_point(size = 3, alpha = 1, shape = 21, color = "black", stroke = 0.6) +
    geom_text_repel(aes(label = text), size = 5, color = "black") +
    scale_fill_manual(values = pt.cols) +
    xlab("") + ylab(y_label) + ggtitle(title) +
    theme_bw(base_size = 18) +
    theme(
      text = element_text(color = "black", face = "plain"),
      axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16,color = "black"),
      axis.text.y  = element_text(size = 15,color = "black"),
      axis.title.x = element_text(size = 17, margin = margin(r = 4),color = "black"),
      axis.title.y = element_text(size = 17, margin = margin(r = 4),color = "black"),
      plot.title   = element_text(size = 20, hjust = 0.5,color = "black"),
      legend.text  = element_text(size = 14,color = "black"),
      legend.title = element_text(size = 15,color = "black"),
      strip.text   = element_text(size = 16,color = "black"),
      legend.position = "none",
      panel.border = element_rect(color = "black", size = 1),
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black"),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1)
    ) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    geom_text(
      data = contrast_results,
      aes(
        x = x, y = y,
        label = ifelse(
          p.value < 0.001,
          paste0(formatC(p.value, format = "e", digits = 1)),
          paste0(signif(p.value, 2))
        ),
        color = ifelse(p.value < 0.05, ifelse(estimate > 0, "#B22222", "#0072B2"), "black")
      ),
      vjust = -0.5, size = 4.5,inherit.aes = FALSE
    ) +
    scale_color_identity() +
    geom_segment(
      data = contrast_results,
      aes(x = x1, xend = x2, y = y1, yend = y2),
      linetype = 1, color = "black", inherit.aes = FALSE
    )

  return(list(plot = p, stat = contrast_results))
}

plot_ct_prot <- function(df, y_var, y_label, title=NULL, pt.cols=pt.cols, ymin=NULL) {
  require(stringr)
  required_cols <- c("Patient.ID", "TimePoint", "bt3","BOR","segment",y_var)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required column(s):", paste(missing_cols, collapse = ", ")))
  }
  df <- df %>%
    mutate(
      bt4 = stringr::str_replace(as.character(bt3),
                                 "(.+)_(.+)",
                                 paste0("\\1_", segment, "_\\2"))
    ) %>%
    mutate(bt4 = factor(
      bt4,
      levels = c(
        "MUT_Tumor_BS", "MUT_Tumor_BX2", "MUT_Tumor_BX3",
        "MUT_Stroma_BS", "MUT_Stroma_BX2", "MUT_Stroma_BX3",
        "WT_Tumor_BS", "WT_Tumor_BX2", "WT_Tumor_BX3",
        "WT_Stroma_BS", "WT_Stroma_BX2", "WT_Stroma_BX3"
      )
    ))
  df.mean <- df %>%
    as.data.frame() %>%
    group_by(!!!syms(c(required_cols[1:5], "bt4"))) %>%
    summarise(mean_val = mean(.data[[y_var]], na.rm = TRUE), .groups = "drop") %>%
    group_by(Patient.ID, segment) %>%
    mutate(
      max_bn = max(as.numeric(TimePoint), na.rm = TRUE),
      text = ifelse(as.numeric(TimePoint) == max_bn, as.character(Patient.ID), "")
    ) %>%
    ungroup()
  
  if (is.null(ymin)) {
    ymin <- min(df[[y_var]], na.rm = TRUE)
  } else if (ymin==0){
    ymin <- min(0, min(df[[y_var]], na.rm = TRUE))
  }
  
  y_data_max <- max(df[[y_var]], na.rm = TRUE)
  y_range <- y_data_max - ymin
  y_mag <- y_data_max * 1.05
  
  conds <- rbind(c(1,2),c(2,3),c(1,3),c(4,5),c(5,6),c(4,6),c(1,4),c(2,5),c(3,6))
  paired <- rep(T,9)
  
  conds1 <- conds + 6
  paired1 <- rep(T,9)
  
  conds2 <- rbind(c(1,7),c(2,8),c(3,9),c(4,10),c(5,11),c(6,12))
  paired2 <- rep(F,6)

  conds <- rbind(conds, conds1, conds2)
  paired <- c(paired, paired1, paired2)
  
  y_offset_levels <- c(rep(c(1,1,2,1,1,2,3,4,5),2), 6,7,8,9,10,11)
  y_offset_levels <- y_offset_levels + 1
  
  ys <- (y_offset_levels * 0.07 + 0.8) * y_range + ymin
  
  df.conds <- data.frame(
    idx = seq_along(paired),
    cond = conds,
    paired = paired,
    p.value = rep(NA, 24),
    x = rowMeans(conds),
    y = ys,
    x1 = conds[,1] + 0.03,
    x2 = conds[,2] - 0.03,
    y1 = ys,
    y2 = ys
  )
  
  model <- lmer(as.formula(paste(y_var, "~ bt4 + (1 | Patient.ID)")), data = df)
  em <- emmeans(model, ~ bt4)
  contrast_matrix <- list(
    "MUT_Tumor_BS vs MUT_Tumor_BX2" = c(1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),
    "MUT_Tumor_BX2 vs MUT_Tumor_BX3" = c(0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0),
    "MUT_Tumor_BS vs MUT_Tumor_BX3" = c(1,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0),
    "MUT_Stroma_BS vs MUT_Stroma_BX2" = c(0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0),
    "MUT_Stroma_BX2 vs MUT_Stroma_BX3" = c(0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0),
    "MUT_Stroma_BS vs MUT_Stroma_BX3" = c(0,  0,  0,  1,  0, -1,  0,  0,  0,  0,  0,  0),
    "MUT_Tumor_BS vs MUT_Stroma_BX2" = c(1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0),
    "MUT_Tumor_BX2 vs MUT_Stroma_BX3" = c(0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0),
    "MUT_Tumor_BS vs MUT_Stroma_BX3" = c(0,  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0),
    
    "WT_Tumor_BS vs WT_Tumor_BX2" = c(0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0),
    "WT_Tumor_BX2 vs WT_Tumor_BX3" = c(0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0),
    "WT_Tumor_BS vs WT_Tumor_BX3" = c(0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0),
    "WT_Stroma_BS vs WT_Stroma_BX2" = c(0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  0),
    "WT_Stroma_BX2 vs WT_Stroma_BX3" = c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1),
    "WT_Stroma_BS vs WT_Stroma_BX3" = c(0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0, -1),
    "WT_Tumor_BS vs WT_Stroma_BX2" = c(0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0),
    "WT_Tumor_BX2 vs WT_Stroma_BX3" = c(0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0),
    "WT_Tumor_BS vs WT_Stroma_BX3" = c(0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1),
  
    "MUT_Tumor_BS vs WT_Tumor_BS" = c(1,  0,  0,  0,  0,  0,  -1,  0,  0,  0,  0,  0),
    "MUT_Tumor_BX2 vs WT_Tumor_BX2" = c(0,  1,  0,  0,  0,  0,  0,  -1,  0,  0,  0,  0),
    "MUT_Tumor_BX3 vs WT_Tumor_BX3" = c(0,  0,  1,  0,  0,  0,  0,  0,  -1,  0,  0,  0),
    "MUT_Stroma_BS vs WT_Stroma_BS" = c(0,  0,  0,  1,  0,  0,  0,  0,  0,  -1,  0,  0),
    "MUT_Stroma_BX2 vs WT_Stroma_BX2" = c(0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  -1,  0),
    "MUT_Stroma_BX3 vs WT_Stroma_BX3" = c(0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  -1)
  )
  
  contrast_results <- contrast(em, contrast_matrix, adjust = "none") %>%
    as.data.frame() %>%
    mutate(estimate = -estimate) %>%
    cbind(df.conds %>% select(x, y, x1, x2, y1, y2))
  
  ymax <- max(max(df[[y_var]], na.rm = TRUE), max(ys)) * 1.01
  
  p <- ggplot(df.mean, aes(x = bt4, y = mean_val, fill = Patient.ID)) +
    geom_path(
      aes(group = interaction(Patient.ID, segment), linetype = BOR),
      color = "grey60",
      size = 0.4
    ) +
    geom_point(size = 3, alpha = 1, shape = 21, color = "black", stroke = 0.6) +
    geom_text_repel(aes(label = text), size = 5, color = "black") +
    scale_fill_manual(values = pt.cols) +
    scale_x_discrete(labels = rep(c("BS", "BX2", "BX3"), 4)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    xlab("") + ylab(y_label) + ggtitle(title) +
    theme_bw(base_size = 18) +
    theme(
      text = element_text(color = "black", face = "plain"),
      axis.text.x  = element_text(angle = 0, vjust = 1, hjust = 0.5, size = 16,color = "black"),
      axis.text.y  = element_text(size = 15,color = "black"),
      axis.title.x = element_text(size = 17, margin = margin(r = 4),color = "black"),
      axis.title.y = element_text(size = 17, margin = margin(r = 4),color = "black"),
      plot.title   = element_text(size = 20, hjust = 0.5,color = "black"),
      legend.text  = element_text(size = 14,color = "black"),
      legend.title = element_text(size = 15,color = "black"),
      strip.text   = element_text(size = 16,color = "black"),
      legend.position = "none",
      panel.border = element_rect(color = "black", size = 1),
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black"),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1)
    ) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    geom_text(
      data = contrast_results,
      aes(
        x = x, y = y,
        label = ifelse(
          p.value < 0.001,
          paste0(formatC(p.value, format = "e", digits = 1)),
          paste0(signif(p.value, 2))
        ),
        color = ifelse(p.value < 0.05, ifelse(estimate > 0, "#B22222", "#0072B2"), "black")
      ),
      vjust = -0.2, size = 4.5,inherit.aes = FALSE
    ) +
    scale_color_identity() +
    geom_segment(
      data = contrast_results,
      aes(x = x1, xend = x2, y = y1, yend = y2),
      linetype = 1, color = "black", inherit.aes = FALSE
    )

  mg_seg_label <- 0.23
  seg_line_down <- 0

  mg_grp_label <- 0.37
  grp_line_up  <- 0.02
  
  segment_labels <- data.frame(
    x = c(2, 5, 8, 11),
    y = ymin - mg_seg_label * y_range,
    label = c("Tumor", "Stroma", "Tumor", "Stroma")
  )
  
  segment_lines <- data.frame(
    x1 = c(1, 4, 7, 10),
    x2 = c(3, 6, 9, 12),
    y  = ymin - (mg_seg_label + seg_line_down) * y_range
  )
  
  mut_labels <- data.frame(
    x = c(3.5, 9.5),
    y = ymin - mg_grp_label * y_range,
    label = c("MUT", "WT")
  )
  
  mut_lines <- data.frame(
    x1 = c(1, 7),
    x2 = c(6, 12),
    y  = ymin - (mg_grp_label - grp_line_up) * y_range
  )

  p <- p +
    geom_segment(data = segment_lines,
                 aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, color = "black", size = 0.6) +
    geom_text(data = segment_labels,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 6, vjust = 1.2) +
    geom_segment(data = mut_lines,
                 aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, color = "black", size = 0.6) +
    geom_text(data = mut_labels,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 7, vjust = 1.2) +
    coord_cartesian(ylim = c(ymin, ymax), clip = "off") +
    theme(plot.margin = margin(t = 5, r = 5, b = 50, l = 5))
  
  return(list(plot = p, stat = contrast_results))
}

make_volcano <- function(df, title = "Volcano",
                         p_cut = 0.05,
                         label_genes = NULL,
                         highlight_genes = NULL,
                         metric = c("p.value", "FDR"),
                         label_p_cut = 0.05,
                         n_label = 10) {
  metric <- match.arg(metric)
  stopifnot(all(c("gene", "BRCA", "estimate", "p.value") %in% colnames(df)))
  if (metric == "FDR") stopifnot("FDR" %in% colnames(df))

  ## metric drives the y-axis value, the dashed reference line, and the gene-labeling
  ## cutoff -- so "the FDR volcano" and "the p-value volcano" are two fully separate,
  ## internally consistent plots, not just a color toggle on top of an always-nominal-p
  ## y-axis. Point color is a separate, three-tier significance indicator when plotting
  ## against nominal p (metric = "p.value") and FDR is available: FDR<p_cut (red),
  ## nominal p<p_cut only (black), neither (grey80). When plotting against FDR directly,
  ## there's no separate "stronger" tier to show, so it stays a simple two-color scheme.
  ## highlight_genes are marked by shape (filled triangle), not color, so they don't collide with
  ## the red = FDR-significant convention -- a highlighted gene's color still reflects its
  ## own significance tier.
  has_fdr <- "FDR" %in% colnames(df)

  df2 <- df %>%
    mutate(
      metric_val = if (metric == "FDR") FDR else p.value,
      neglog10p = -log10(metric_val),
      sig_tier = if (metric == "p.value" && has_fdr) {
        dplyr::case_when(
          FDR < p_cut ~ "fdr",
          p.value < p_cut ~ "nominal",
          TRUE ~ "ns"
        )
      } else {
        ifelse(metric_val < p_cut, "fdr", "ns")
      },
      is_highlight = if (!is.null(highlight_genes)) gene %in% highlight_genes else FALSE
    )

  y_lab <- if (metric == "FDR") "-log10(FDR)" else "-log10(p)"
  show_legend <- metric == "p.value" && has_fdr

  p <- ggplot(df2, aes(x = estimate, y = neglog10p)) +
    geom_point(aes(color = sig_tier), size = 1.6, alpha = 0.95) +
    scale_color_manual(
      values = c("fdr" = "red", "nominal" = "black", "ns" = "grey80"),
      breaks = c("fdr", "nominal"),
      labels = c(paste0("FDR < ", p_cut), paste0("nominal p < ", p_cut)),
      name = "Significance"
    ) +
    geom_point(
      data = dplyr::filter(df2, is_highlight),
      aes(color = sig_tier),
      shape = 17,
      size = 2.6,
      alpha = 1,
      show.legend = FALSE
    ) +
    facet_wrap(~ BRCA, nrow = 1) +
    theme_bw() +
    geom_vline(xintercept = 0, linetype = "solid") +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed") +
    labs(x = "Effect (Δlog expr)", y = y_lab, title = title) +
    guides(color = if (show_legend) "legend" else "none")

  if (!is.null(label_genes)) {
    df_lab <- df2 %>%
      filter(gene %in% label_genes, metric_val < label_p_cut,
             !is.na(neglog10p), !is.na(estimate))
  } else {
    df_lab <- df2 %>%
      filter(metric_val < label_p_cut) %>%
      group_by(BRCA) %>%
      slice_min(order_by = metric_val, n = n_label, with_ties = FALSE) %>%
      ungroup()
  }

  if (!is.null(highlight_genes)) {
    df_hlab <- df2 %>%
      filter(gene %in% highlight_genes,
             !is.na(neglog10p), !is.na(estimate))
    df_lab <- bind_rows(df_lab, df_hlab) %>%
      distinct(gene, BRCA, .keep_all = TRUE)
  }

  if (nrow(df_lab) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = df_lab,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.2
    )
  }

  p
}

