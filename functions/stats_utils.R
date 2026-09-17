## General-purpose statistics helpers shared across platforms (not specific to any one
## assay), as opposed to functions/analysis_{io360,dsp_protein,cycif}.R.

## Repeated-subsample Shapiro-Wilk normality check. shapiro.test() is capped at n=5000 and
## becomes hyper-sensitive (rejects for trivial deviations) at large n regardless; a single
## subsample also has its own problem (the result depends on which draw you happen to get).
## Address both by drawing n_draws independent subsamples of size draw_size and looking at
## the distribution of W (the bounded, sample-size-robust effect-size-like statistic) across
## draws, rather than trusting a single p-value from a single draw. Used identically by the
## normality-diagnostic sections in IO360, DSP Protein, and CyCIF.
run_shapiro_draws <- function(x, draw_size, n_draws = 20) {
  draws <- replicate(
    n_draws,
    shapiro.test(sample(x, min(length(x), draw_size))),
    simplify = FALSE
  )
  tibble(
    W = sapply(draws, function(d) unname(d$statistic)),
    p = sapply(draws, function(d) d$p.value)
  )
}
