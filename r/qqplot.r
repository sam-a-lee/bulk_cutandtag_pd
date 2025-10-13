qq_p_gg <- function(p, ci = 0.95) {
  
  p <- p[is.finite(p) & !is.na(p)]
  p[p <= 0] <- .Machine$double.xmin
  p[p > 1] <- 1
  
  o <- sort(p)
  n <- length(o); i <- seq_len(n)
  e <- (i - 0.5) / n
  alpha <- (1 - ci) / 2
  lo <- qbeta(alpha, i, n - i + 1)
  hi <- qbeta(1 - alpha, i, n - i + 1)
  
  df <- data.frame(
    exp = -log10(e),
    obs = -log10(o),
    lower = -log10(lo),
    upper = -log10(hi)
  )
  
  ggplot(df, aes(exp, obs)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_point(size = 0.8, alpha = 0.6, color = "#00033f") +
    coord_equal() +
    labs(x = "Expected -log10(p)", y = "Observed -log10(p)") +
    theme_minimal()
}

lambda_from_p <- function(p, df = 1, na.rm = TRUE) {
  p <- p[is.finite(p)]
  if (na.rm) p <- p[!is.na(p)]
  chisq_obs <- qchisq(1 - p, df = df)
  median(chisq_obs) / qchisq(0.5, df = df)
}
