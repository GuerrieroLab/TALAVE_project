prep_tmp <- function(df) {
  df %>%
    filter(!is.na(BRCA), !is.na(TimePoint), !is.na(Patient.ID), !is.na(expr))
}

fit_lmer_contrasts <- function(df_one_gene, brca_level = c("WT", "MUT")) {
  brca_level <- match.arg(brca_level)

  d <- df_one_gene %>% filter(BRCA == brca_level)

  if (n_distinct(d$TimePoint) < 2 || n_distinct(d$Patient.ID) < 2) {
    return(tibble(
      BRCA = brca_level,
      contrast = c("BX2 - BS", "BX3 - BS"),
      estimate = NA_real_,
      SE = NA_real_,
      df = NA_real_,
      t.ratio = NA_real_,
      p.value = NA_real_
    ))
  }

  mod <- lmer(expr ~ TimePoint + (1 | Patient.ID), data = d, REML = FALSE)
  em <- emmeans(mod, ~ TimePoint)

  ct <- contrast(em, method = list(
    "BX2 - BS" = c(-1, 1, 0),
    "BX3 - BS" = c(-1, 0, 1)
  ))

  as_tibble(summary(ct, infer = c(TRUE, TRUE))) %>%
    transmute(
      BRCA = brca_level,
      contrast = contrast,
      estimate = estimate,
      SE = SE,
      df = df,
      t.ratio = t.ratio,
      p.value = p.value
    )
}

get_gene_lmer_residuals <- function(df_one_gene, brca_level = c("WT", "MUT")) {
  brca_level <- match.arg(brca_level)

  d <- df_one_gene %>% filter(BRCA == brca_level)

  if (n_distinct(d$TimePoint) < 2 || n_distinct(d$Patient.ID) < 2) {
    return(NULL)
  }

  mod <- lmer(expr ~ TimePoint + (1 | Patient.ID), data = d, REML = FALSE)
  as.numeric(residuals(mod))
}

## Joint-model versions: one lmer(expr ~ bt3 + (1|Patient.ID)) per gene, bt3 combining
## BRCA x TimePoint into a single 6-level factor -- same scheme as plot_ct() in
## functions/plotting.R. Unlike fit_lmer_contrasts()/get_gene_lmer_residuals() (which fit
## two independent models, one per BRCA group, and so can never test MUT vs WT directly),
## this pools variance estimation across both groups and additionally makes MUT-vs-WT
## contrasts possible. Output shape matches fit_lmer_contrasts() exactly (BRCA, contrast,
## estimate, SE, df, t.ratio, p.value) so it's a drop-in replacement wherever that shape is
## expected (e.g. all_res).
bt3_levels <- c("MUT_BS", "MUT_BX2", "MUT_BX3", "WT_BS", "WT_BX2", "WT_BX3")

## ---- GSEA / pathway-enrichment helpers ----
## Moved here from a local definition in TALAVE_IO360/11_io360_expression_analysis.Rmd
## (organizational cleanup; not reused elsewhere in the project, but analysis logic belongs
## here rather than in the notebook per this project's existing convention).

make_rank_signed_logp <- function(df, brca, use = c("p.value","FDR"), eps = 1e-300) {
  use <- match.arg(use)
  stopifnot(all(c("gene","BRCA","estimate",use) %in% colnames(df)))

  df %>%
    filter(BRCA == brca, !is.na(gene), gene != "", !is.na(estimate), !is.na(.data[[use]])) %>%
    mutate(
      p = pmax(.data[[use]], eps),
      signed_logp = sign(estimate) * (-log10(p))
    ) %>%
    group_by(gene) %>%
    summarise(score = max(signed_logp, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(score)) %>%
    tibble::deframe()
}

run_fgsea_rank <- function(rank_vec, pathways, nperm = 10000, minSize = 3, maxSize = 500) {
  fgsea(pathways = pathways, stats = rank_vec,
        nperm = nperm, minSize = minSize, maxSize = maxSize) %>%
    as_tibble() %>%
    arrange(padj)
}

classify_change <- function(mut_nes, mut_padj, wt_nes, wt_padj, padj_cut = 0.05) {
  case_when(
    mut_padj < padj_cut & wt_padj >= padj_cut & mut_nes > 0 ~ "Mut.up",
    mut_padj < padj_cut & wt_padj >= padj_cut & mut_nes < 0 ~ "Mut.dn",
    wt_padj  < padj_cut & mut_padj >= padj_cut & wt_nes  > 0 ~ "Wt.up",
    wt_padj  < padj_cut & mut_padj >= padj_cut & wt_nes  < 0 ~ "Wt.dn",
    mut_padj < padj_cut & wt_padj  < padj_cut ~ "Both.sig",
    TRUE ~ "Others"
  )
}

combine_fgsea <- function(fg_mut, fg_wt, padj_cut = 0.05) {
  fg_mut %>%
    select(pathway, mut_NES = NES, mut_pval = pval, mut_padj = padj) %>%
    left_join(
      fg_wt %>% select(pathway, wt_NES = NES, wt_pval = pval, wt_padj = padj),
      by = "pathway"
    ) %>%
    mutate(
      mut_slp = sign(mut_NES) * -log10(pmax(mut_pval, 1e-300)),
      wt_slp  = sign(wt_NES)  * -log10(pmax(wt_pval,  1e-300)),
      Change  = classify_change(mut_NES, mut_padj, wt_NES, wt_padj, padj_cut = padj_cut)
    )
}

format_pathway_label_ascii <- function(x) {
  source_tag <- dplyr::case_when(
    grepl("^HALLMARK_", x) ~ "(H)",
    grepl("^NANOSTRING_", x) ~ "(N)",
    grepl("^CUSTOM_", x) ~ "(C)",
    TRUE ~ ""
  )

  clean <- x
  clean <- gsub("^HALLMARK_", "", clean)
  clean <- gsub("^NANOSTRING_", "", clean)
  clean <- gsub("^CUSTOM_", "", clean)
  clean <- gsub("_", " ", clean)
  clean <- gsub("\\.", " ", clean)
  clean <- tolower(clean)

  clean <- gsub("fa brca hr", "Fanconi anemia-BRCA HR", clean)
  clean <- gsub("\\bnejh\\b", "NEJH", clean)
  clean <- gsub("\\bt cells\\b", "T cells", clean)
  clean <- gsub("\\bb cells\\b", "B cells", clean)
  clean <- gsub("\\be2f\\b", "E2F", clean)
  clean <- gsub("\\bmtorc1\\b", "MTORC1", clean)
  clean <- gsub("\\bg2m\\b", "G2M", clean)
  clean <- gsub("\\bmyc\\b", "MYC", clean)
  clean <- gsub("\\bkras\\b", "KRAS", clean)

  clean <- gsub("il6 jak stat3 signaling", "IL6-JAK-STAT3 signaling", clean)
  clean <- gsub("il2 stat5 signaling", "IL2-STAT5 signaling", clean)
  clean <- gsub("cgas sting", "cGAS-STING", clean)
  clean <- gsub("tnfa signaling via nfkb", "TNF-alpha signaling via NF-kB", clean)
  clean <- gsub("ifn downstream", "IFN downstream", clean)
  clean <- gsub("tgf beta signaling", "TGF-beta signaling", clean)

  clean <- gsub("epithelial mesenchymal transition", "epithelial-mesenchymal transition", clean)

  clean <- ifelse(
    grepl("^[A-Z]", clean),
    clean,
    paste0(toupper(substr(clean, 1, 1)), substr(clean, 2, nchar(clean)))
  )

  paste0(clean, " ", source_tag)
}

extract_signed_logp <- function(fg, pathway_name, condition_label, use_adj = TRUE, eps = 1e-300) {
  stopifnot(all(c("pathway","NES","pval","padj") %in% colnames(fg)))
  pcol <- if (use_adj) "padj" else "pval"

  fg %>%
    filter(pathway == pathway_name) %>%
    transmute(
      Condition = condition_label,
      pathway = pathway_name,
      NES = NES,
      p = .data[[pcol]],
      signed_logp = sign(NES) * (-log10(pmax(p, eps)))
    )
}
