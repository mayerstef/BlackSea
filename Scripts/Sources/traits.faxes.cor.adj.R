traits.faxes.cor.adj<-function (sp_tr, sp_faxes_coord, tr_nm = NULL, faxes_nm = NULL, 
          plot = FALSE, name_file = NULL, color_signif = "darkblue", 
          color_non_signif = "gray80", stop_if_NA = TRUE) 
{
  if (missing(sp_tr)) {
    stop("Argument 'sp_tr' is mandatory.")
  }
  if (missing(sp_faxes_coord)) {
    stop("Argument 'sp_faxes_coord' is mandatory.")
  }
  check.sp.tr(sp_tr, stop_if_NA = stop_if_NA)
  check.sp.faxes.coord(sp_faxes_coord)
  if (!identical(sort(rownames(sp_tr)), sort(rownames(sp_faxes_coord)))) {
    stop("Species names mismatch between 'sp_tr' and 'sp_faxes_coord'.")
  }
  if (is.null(tr_nm)) {
    tr_nm <- names(sp_tr)
  }
  else {
    if (any(!(tr_nm %in% names(sp_tr)))) {
      stop("Trait names should be as in 'sp_tr'.")
    }
  }
  if (is.null(faxes_nm)) {
    faxes_nm <- colnames(sp_faxes_coord)
  }
  else {
    if (any(!(faxes_nm %in% colnames(sp_faxes_coord)))) {
      stop("Axes names should be as in 'sp_tr'.")
    }
  }
  if (plot) {
    if (length(faxes_nm) > 10) {
      stop("Number of axes to plot should be < 11.")
    }
    if (length(tr_nm) > 10) {
      stop("Number of traits to plot should be < 11.")
    }
  }
  res <- as.data.frame(matrix(NA, length(tr_nm) * length(faxes_nm), 
                              7, dimnames = list(NULL, c("trait", "axis", "test", "stat", 
                                                         "value", "p.value", "adj.p.value"))))
  flag <- 0
  for (i in tr_nm) {
    for (j in faxes_nm) {
      flag <- flag + 1
      trait <- NULL
      axis <- NULL
      data_ij <- data.frame(trait = sp_tr[, i], axis = sp_faxes_coord[, 
                                                                      j])
      if (is.numeric(data_ij$trait)) {
        test_ij <- "Linear Model"
        lm_ij <- summary(stats::lm(trait ~ axis, data = data_ij))
        
        # Store unadjusted p-value in a separate variable
        raw_p_value <- lm_ij$coefficients[2, 4]
        
        # Adjust p-value using Benjamini-Hochberg procedure
        adj_p_value <- p.adjust(raw_p_value, method = "BH")
        
        res[flag, c("trait", "axis", "test", "stat")] <- c(i, 
                                                           j, test_ij, "r2")
        res[flag, c("value", "p.value", "adj.p.value")] <- c(round(lm_ij$r.squared, 
                                                                   3), round(raw_p_value, 4), round(adj_p_value, 4))
      }
      else {
        test_ij <- "Kruskal-Wallis"
        kw_ij <- stats::kruskal.test(axis ~ trait, data = data_ij)
        kw_ij_eta2 <- rstatix::kruskal_effsize(data = data_ij, 
                                               axis ~ trait)
        
        # Store unadjusted p-value in a separate variable
        raw_p_value <- kw_ij$p.value
        
        # Adjust p-value using Benjamini-Hochberg procedure
        adj_p_value <- p.adjust(raw_p_value, method = "BH")
        
        res[flag, c("trait", "axis", "test", "stat")] <- c(i, 
                                                           j, test_ij, "eta2")
        res[flag, c("value", "p.value", "adj.p.value")] <- c(round(kw_ij_eta2$effsize, 
                                                                   3), round(raw_p_value, 4), round(adj_p_value, 4))
      }
      if (plot) {
        x_lab_ij <- NULL
        if (j == faxes_nm[length(faxes_nm)]) {
          x_lab_ij <- i
        }
        y_lab_ij <- NULL
        if (i == tr_nm[1]) {
          y_lab_ij <- j
        }
        if (res[flag, "adj.p.value"] < 0.05) {  # use adj.p.value instead of p.value
          col_cor <- color_signif
        }
        else {
          col_cor <- color_non_signif
        }
        gg_ij <- ggplot2::ggplot(data_ij, ggplot2::aes(trait, 
                                                       axis)) + ggplot2::xlab(x_lab_ij) + ggplot2::ylab(y_lab_ij) + 
          ggplot2::theme_bw()
        if (res[flag, "stat"] == "r2") {
          gg_ij <- gg_ij + ggplot2::geom_point(size = 2, 
                                               col = col_cor)
        }
        else {
          gg_ij <- gg_ij + ggplot2::geom_boxplot(colour = col_cor) + 
            ggplot2::geom_jitter(colour = col_cor, width = 0.3, 
                                 size = 1.5)
        }
        if (j == faxes_nm[1]) {
          col_j <- gg_ij
        }
        else {
          col_j <- col_j/gg_ij
        }
      }
    }
    if (plot) {
      if (i == tr_nm[1]) {
        tr_faxes_plot <- col_j
      }
      else {
        tr_faxes_plot <- tr_faxes_plot | col_j
      }
    }
  }
  if (!plot) {
    return(tr_faxes_stat = res)
  }
  else {
    tr_faxes_plot <- tr_faxes_plot + patchwork::plot_annotation(title = "Relation between traits and PCoA axes", 
                                                                caption = "Made with mFD package")
    if (is.null(name_file)) {
      return(list(tr_faxes_stat = res, tr_faxes_plot = tr_faxes_plot))
    }
    else {
      ggplot2::ggsave(filename = paste0(name_file, ".jpeg"), 
                      plot = tr_faxes_plot, device = "jpeg", width = length(tr_nm) * 
                        3, height = length(faxes_nm) * 2.33, units = "in", 
                      dpi = 300)
      return(tr_faxes_stat = res)
    }
  }
}
