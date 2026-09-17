
## Residual extraction mirroring plot_ct_prot()'s exact model (lmer(y_var ~ bt4 +
## (1|Patient.ID)), bt4 = BRCA x Segment x TimePoint, AOI-level data) -- for the A8
## normality check. Unlike IO360 (separate model per BRCA group), this is DSP Protein's
## existing, long-standing design: one joint model per protein spanning all of
## BRCA/Segment/Time. Not changed here, just diagnosed.
get_protein_lmer_residuals <- function(df, y_var) {
  required_cols <- c("Patient.ID", "TimePoint", "bt3", "BOR", "segment", y_var)
  if (!all(required_cols %in% names(df))) return(NULL)

  df2 <- df %>%
    mutate(
      bt4 = stringr::str_replace(as.character(bt3), "(.+)_(.+)", paste0("\\1_", segment, "_\\2"))
    ) %>%
    mutate(bt4 = factor(bt4, levels = c(
      "MUT_Tumor_BS", "MUT_Tumor_BX2", "MUT_Tumor_BX3",
      "MUT_Stroma_BS", "MUT_Stroma_BX2", "MUT_Stroma_BX3",
      "WT_Tumor_BS", "WT_Tumor_BX2", "WT_Tumor_BX3",
      "WT_Stroma_BS", "WT_Stroma_BX2", "WT_Stroma_BX3"
    )))

  mod <- tryCatch(
    lmer(as.formula(paste(y_var, "~ bt4 + (1 | Patient.ID)")), data = df2),
    error = function(e) NULL
  )
  if (is.null(mod)) return(NULL)

  as.numeric(residuals(mod))
}

build_dsp_protein_contrast_table <- function(prot.stats,
                                             keep_idx = c(1, 3, 4, 6, 10, 12, 13, 15)) {
  prot_long <- lapply(names(prot.stats), function(protein_name) {
    stat_df <- prot.stats[[protein_name]][keep_idx, ]

    data.frame(
      Protein = protein_name,
      Contrast = stat_df$contrast,
      estimate = stat_df$estimate,
      p.value = stat_df$p.value,
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows() %>%
    group_by(Contrast) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    mutate(
      signed_logp = sign(estimate) * -log10(p.value),
      signed_logq = sign(estimate) * -log10(p.adj)
    )

  prot_long
}

prepare_dsp_protein_heatmap_data <- function(prot_long, cap = 4, sig_cut = 0.05) {
  ## Same two-tier significance convention as IO360's heatmaps: color reflects nominal p
  ## (signed_logp), box border reflects a 3-level tier -- solid black for FDR<sig_cut,
  ## grey for nominal p<sig_cut only, no box for neither.
  prot_long2 <- prot_long %>%
    mutate(tier = case_when(
      p.adj < sig_cut ~ "fdr",
      p.value < sig_cut ~ "nominal",
      TRUE ~ "none"
    ))

  keep_prots <- prot_long2 %>%
    group_by(Protein) %>%
    summarise(any_sig = any(tier != "none", na.rm = TRUE), .groups = "drop") %>%
    filter(any_sig) %>%
    pull(Protein)

  heat_mat <- prot_long2 %>%
    filter(Protein %in% keep_prots) %>%
    select(Protein, Contrast, signed_logp) %>%
    pivot_wider(names_from = Contrast, values_from = signed_logp) %>%
    column_to_rownames("Protein") %>%
    as.matrix()

  tier_mat <- prot_long2 %>%
    filter(Protein %in% keep_prots) %>%
    select(Protein, Contrast, tier) %>%
    pivot_wider(names_from = Contrast, values_from = tier) %>%
    column_to_rownames("Protein") %>%
    as.matrix()

  tier_mat <- tier_mat[rownames(heat_mat), colnames(heat_mat)]
  tier_mat[is.na(tier_mat)] <- "none"

  heat_mat[heat_mat > cap] <- cap
  heat_mat[heat_mat < -cap] <- -cap

  rownames(heat_mat) <- gsub("\\.", "-", rownames(heat_mat))
  colnames(heat_mat) <- colnames(heat_mat) |>
    gsub(".* vs ", "", x = _) |>
    gsub("\\.", "-", x = _)

  row_dist <- stats::as.dist(1 - stats::cor(t(heat_mat), use = "pairwise.complete.obs"))
  row_hc <- stats::hclust(row_dist, method = "average")

  col_meta <- data.frame(
    BRCA = sub("_.*", "", colnames(heat_mat)),
    Segment = sub("^[^_]+_([^_]+)_.*", "\\1", colnames(heat_mat)),
    Time = sub(".*_", "", colnames(heat_mat)),
    stringsAsFactors = FALSE
  )

  ann_colors <- list(
    BRCA = c(MUT = "#1b9e77", WT = "#7570b3"),
    Segment = c(Tumor = "#666666", Stroma = "#bdbdbd"),
    Time = c(BX2 = "#fc8d62", BX3 = "#8da0cb")
  )

  list(
    prot_long = prot_long2,
    keep_prots = keep_prots,
    heat_mat = heat_mat,
    tier_mat = tier_mat,
    row_hc = row_hc,
    cap = cap,
    col_meta = col_meta,
    ann_colors = ann_colors
  )
}
