# heat + scree with ALL PCs in the heat map (PC1 included),
# now with a SECOND heat map of p-values under the coefficient map.
# - Loadings: PC scores per sample (e.g., p$x from prcomp)
# - Importance: proportion of variance per PC (e.g., (p$sdev^2)/sum(p$sdev^2))
# - Categorical: data.frame of factors (2-level -> r + p from Pearson test; >2-level -> eta^2 + p from ANOVA)
# - Continuous: data.frame of numeric (cor.test; Spearman by default)
heat_scree_plot <- function(Loadings, Importance, Num, Order,
                            Categorical = NULL, Continuous = NULL,
                            continuous_method = c("spearman", "pearson")) {
  continuous_method <- match.arg(continuous_method)
  
  # ---- bounds & scree (PC1 shown) ----
  Num <- min(Num, length(Importance), ncol(Loadings))
  scree_df <- data.frame(variance = Importance[1:Num],
                         PC       = seq_len(Num))
  scree <- ggplot2::ggplot(scree_df, ggplot2::aes(PC, variance)) +
    ggplot2::geom_bar(stat = "identity", color = "black", fill = "grey") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 10),
                   axis.title = ggplot2::element_text(size = 10),
                   plot.margin = grid::unit(c(0, 2.25, 0.5, 1), "cm")) +
    ggplot2::ylab("Variance") +
    ggplot2::scale_x_continuous(breaks = seq_len(Num))
  
  # ---- compute association coefficients AND p-values for ALL PCs ----
  # we'll build two matrices with the same shape: coef_mat (values), p_mat (p-values)
  coef_mat <- NULL
  p_mat    <- NULL
  
  # ----- Categorical -----
  if (!is.null(Categorical)) {
    cat_coef_rows <- list()
    cat_p_rows    <- list()
    for (covar in seq_len(ncol(Categorical))) {
      x <- Categorical[, covar]
      if (is.logical(x) || is.character(x)) x <- as.factor(x)
      x <- droplevels(as.factor(x))
      k <- nlevels(x)
      
      if (k == 2L) {
        # 2-level: point-biserial r (Pearson on 0/1) + p-value from cor.test
        xb <- as.numeric(x) - 1L
        coef_vec <- numeric(ncol(Loadings))
        p_vec    <- numeric(ncol(Loadings))
        for (PC in seq_len(ncol(Loadings))) {
          res <- suppressWarnings(stats::cor.test(Loadings[, PC], xb,
                                                  use = "pairwise.complete.obs",
                                                  method = "pearson"))
          coef_vec[PC] <- unname(res$estimate)
          p_vec[PC]    <- res$p.value
        }
      } else {
        # >2 levels: eta^2 from ANOVA as effect size + p-value from ANOVA F-test
        coef_vec <- numeric(ncol(Loadings))
        p_vec    <- numeric(ncol(Loadings))
        for (PC in seq_len(ncol(Loadings))) {
          fit <- stats::aov(Loadings[, PC] ~ x)
          sm  <- summary(fit)[[1]]
          ss  <- sm[, "Sum Sq"]
          coef_vec[PC] <- as.numeric(ss[1] / sum(ss))    # eta^2 in [0,1]
          p_vec[PC]    <- sm[1, "Pr(>F)"]
        }
      }
      cat_coef_rows[[covar]] <- coef_vec
      cat_p_rows[[covar]]    <- p_vec
    }
    coef_cat <- do.call(rbind, cat_coef_rows); rownames(coef_cat) <- colnames(Categorical)
    p_cat    <- do.call(rbind, cat_p_rows);    rownames(p_cat)    <- colnames(Categorical)
    
    coef_mat <- coef_cat
    p_mat    <- p_cat
  }
  
  # ----- Continuous -----
  if (!is.null(Continuous)) {
    cont_coef_rows <- list()
    cont_p_rows    <- list()
    for (covar in seq_len(ncol(Continuous))) {
      y <- as.numeric(Continuous[, covar])
      coef_vec <- numeric(ncol(Loadings))
      p_vec    <- numeric(ncol(Loadings))
      for (PC in seq_len(ncol(Loadings))) {
        res <- suppressWarnings(stats::cor.test(Loadings[, PC], y,
                                                use = "pairwise.complete.obs",
                                                method = continuous_method,
                                                exact = FALSE))
        coef_vec[PC] <- unname(res$estimate)
        p_vec[PC]    <- res$p.value
      }
      cont_coef_rows[[covar]] <- coef_vec
      cont_p_rows[[covar]]    <- p_vec
    }
    coef_cont <- do.call(rbind, cont_coef_rows); rownames(coef_cont) <- colnames(Continuous)
    p_cont    <- do.call(rbind, cont_p_rows);    rownames(p_cont)    <- colnames(Continuous)
    
    if (is.null(coef_mat)) {
      coef_mat <- coef_cont; p_mat <- p_cont
    } else {
      coef_mat <- rbind(coef_mat, coef_cont)
      p_mat    <- rbind(p_mat,    p_cont)
    }
  }
  
  if (is.null(coef_mat) || is.null(p_mat)) stop("Provide Categorical and/or Continuous metadata.")
  
  # Keep only PCs 1..Num and set column names to PC1..PCNum
  coef_mat <- as.data.frame(coef_mat[, 1:Num, drop = FALSE])
  p_mat    <- as.data.frame(p_mat[,    1:Num, drop = FALSE])
  colnames(coef_mat) <- paste0("PC", seq_len(Num))
  colnames(p_mat)    <- paste0("PC", seq_len(Num))
  
  # ---- to long format, with identical ordering ----
  coef_mat$meta <- rownames(coef_mat)
  p_mat$meta    <- rownames(p_mat)
  
  coef_long <- reshape2::melt(coef_mat, id.vars = "meta",
                              variable.name = "PC", value.name = "coef")
  p_long    <- reshape2::melt(p_mat,    id.vars = "meta",
                              variable.name = "PC", value.name = "p")
  
  # Row order by user-specified 'Order'
  ord <- Order
  row_levels <- unique(coef_long$meta)[rev(ord)]
  coef_long$meta <- factor(coef_long$meta, levels = row_levels)
  p_long$meta    <- factor(p_long$meta,    levels = row_levels)
  
  # Ensure identical column order (PC1..PCNum)
  coef_long$PC <- factor(coef_long$PC, levels = paste0("PC", seq_len(Num)))
  p_long$PC    <- factor(p_long$PC,    levels = paste0("PC", seq_len(Num)))
  
  # ---- heat map: coefficients ----
  heat_coef <- ggplot2::ggplot(coef_long, ggplot2::aes(PC, meta, fill = coef)) +
    ggplot2::geom_tile(color = "black", size = 0.5) +
    ggplot2::theme_gray(8) +
    ggplot2::scale_fill_gradient2(low = "#0571b0", mid = "white", high = "#ca0020",
                                  midpoint = 0, limits = c(-1, 1),
                                  oob = scales::squish, name = "Coefficient") +
    ggplot2::theme(axis.text  = ggplot2::element_text(size = 10, color = "black"),
                   axis.text.x = ggplot2::element_text(),
                   axis.title = ggplot2::element_text(size = 12),
                   legend.text  = ggplot2::element_text(size = 8),
                   legend.title = ggplot2::element_text(size = 10),
                   legend.position = "bottom",
                   legend.key.height = unit(0.25, "cm"),
                   legend.key.width = unit(0.25, "cm"),
                   plot.margin = grid::unit(c(0, 2.25, 0.5, 1), "cm")) +
    ggplot2::xlab("Principal Component") + ggplot2::ylab(NULL)
  
  # ---- heat map: p-values (same grid, stacked below) ----
  heat_p <- ggplot2::ggplot(p_long, ggplot2::aes(PC, meta, fill = p)) +
    ggplot2::geom_tile(color = "black", size = 0.5) +
    ggplot2::theme_gray(8) +
    ggplot2::scale_fill_gradient(limits = c(0, 1),
                                 low = "#08306b", high = "#ffffff",
                                 name = "p-value",
                                 labels = function(x) format(x, digits = 2, scientific = TRUE)) +
    ggplot2::theme(axis.text  = ggplot2::element_text(size = 10, color = "black"),
                   axis.text.x = ggplot2::element_text(),
                   axis.title = ggplot2::element_text(size = 12),
                   legend.text  = ggplot2::element_text(size = 8),
                   legend.title = ggplot2::element_text(size = 10),
                   legend.position = "bottom",
                   legend.key.height = unit(0.25, "cm"),
                   legend.key.width = unit(0.25, "cm"),
                   plot.margin = grid::unit(c(0, 2.25, 0.5, 1), "cm")) +
    ggplot2::xlab("Principal Component") + ggplot2::ylab(NULL)
  
  # ---- layout: scree (top) + coefficient heatmap + p-value heatmap ----
  cowplot::plot_grid(
    scree,
    heat_coef,
   # heat_p,
    ncol = 1,
    rel_heights = c(0.5, 0.5)
  )
}
