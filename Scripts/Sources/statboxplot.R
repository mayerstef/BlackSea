
.centrality_ggrepel <- function(plot,
                                data,
                                x,
                                y,
                                centrality.path = FALSE,
                                centrality.path.args = list(linewidth = 1, color = "red", alpha = 0.5),
                                centrality.point.args = list(size = 5, color = "darkred"),
                                centrality.label.args = list(size = 3, nudge_x = 0.4, segment.linetype = 4),
                                ...) {
  centrality_df <- suppressWarnings(centrality_description(data, {{ x }}, {{ y }}, ...))
  
  # if there should be lines connecting mean values across groups
  if (isTRUE(centrality.path)) {
    plot <- plot +
      exec(
        geom_path,
        data = centrality_df,
        mapping = aes({{ x }}, {{ y }}, group = 1L),
        inherit.aes = FALSE,
        !!!centrality.path.args
      )
  }
  
  plot + # highlight the mean of each group
    exec(
      geom_point,
      mapping = aes({{ x }}, {{ y }}),
      data = centrality_df,
      inherit.aes = FALSE,
      !!!centrality.point.args
    ) + # attach the labels with means to the plot
    exec(
      ggrepel::geom_label_repel,
      data = centrality_df,
      mapping = aes({{ x }}, {{ y }}, label = expression),
      inherit.aes = FALSE,
      parse = TRUE,
      !!!centrality.label.args
    ) + # adding sample size labels to the x axes
    scale_x_discrete(labels = unique(centrality_df$n.expression))
}

.ggsignif_adder <- function(plot,
                            data,
                            x,
                            y,
                            mpc_df,
                            pairwise.display = "significant",
                            ggsignif.args = list(textsize = 3, tip_length = 0.01, na.rm = TRUE),
                            ...) {
  # creating a column for group combinations
  mpc_df %<>% mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c))
  
  # for Bayes Factor, there will be no "p.value" column
  if ("p.value" %in% names(mpc_df)) {
    if (startsWith(pairwise.display, "s")) mpc_df %<>% filter(p.value < 0.05) # sig
    if (startsWith(pairwise.display, "n")) mpc_df %<>% filter(p.value >= 0.05) # non-sig
    
    # proceed only if there are any significant comparisons to display
    if (nrow(mpc_df) == 0L) {
      return(plot)
    }
  }
  
  # arrange the data frame so that annotations are properly aligned
  mpc_df %<>% arrange(group1, group2)
  
  # adding ggsignif comparisons to the plot
  plot +
    exec(
      ggsignif::geom_signif,
      comparisons = mpc_df$groups,
      map_signif_level = TRUE,
      y_position = .ggsignif_xy(pull(data, {{ x }}), pull(data, {{ y }})),
      annotations = as.character(mpc_df$expression),
      test = NULL,
      parse = TRUE,
      !!!ggsignif.args
    )
}


.ggsignif_xy <- function(x, y) {
  # number of comparisons and size of each step
  n_comps <- length(utils::combn(x = unique(x), m = 2L, simplify = FALSE))
  step_length <- (max(y, na.rm = TRUE) - min(y, na.rm = TRUE)) / 20
  
  # start and end position on `y`-axis for the `ggsignif` lines
  y_start <- max(y, na.rm = TRUE) * (1 + 0.025)
  y_end <- y_start + (step_length * n_comps)
  
  # creating a vector of positions for the `ggsignif` lines
  seq(y_start, y_end, length.out = n_comps)
}

.pairwise_seclabel <- function(test.description, pairwise.display = "significant") {
  # single quote (') needs to be escaped inside glue expressions
  test <- sub("'", "\\'", test.description, fixed = TRUE)
  
  # which comparisons were displayed?
  display <- case_when(
    substr(pairwise.display, 1L, 1L) == "s" ~ "significant",
    substr(pairwise.display, 1L, 1L) == "n" ~ "non-significant",
    TRUE ~ "all"
  )
  
  parse(text = glue("list('Pairwise test:'~bold('{test}'), 'Bars shown:'~bold('{display}'))"))
}


.aesthetic_addon <- function(plot,
                             x,
                             fill_var,
                             xlab = NULL,
                             ylab = NULL,
                             title = NULL,
                             subtitle = NULL,
                             caption = NULL,
                             seclabel = NULL,
                             ggtheme = ggstatsplot::theme_ggstatsplot(),
                             package = "RColorBrewer",
                             palette = "Dark2",
                             ggplot.component = NULL,
                             ...) {
  # if no. of factor levels is greater than the default palette color count
  .palette_message(package, palette, length(unique(levels(x)))[[1L]])
  
  plot +
    labs(
      x        = xlab,
      y        = ylab,
      title    = title,
      subtitle = subtitle,
      caption  = caption,
      color    = xlab
    ) +
    ggtheme +
    # no matter the theme, the following ought to be part of a ggstatsplot plot
    theme(legend.position = "none") +
    paletteer::scale_color_paletteer_d(paste0(package, "::", palette)) +
    scale_y_continuous(sec.axis = dup_axis(name = seclabel, breaks = NULL, labels = NULL)) +
    # this is the hail mary way for users to override these defaults
    ggplot.component
}

.f_switch <- function(test) ifelse(test == "t", two_sample_test, oneway_anova)

.onAttach <- function(lib, pkg) {
  packageStartupMessage(
    "You can cite this package as:
     Patil, I. (2021). Visualizations with statistical details: The 'ggstatsplot' approach.
     Journal of Open Source Software, 6(61), 3167, doi:10.21105/joss.03167"
  )
}

.histo_labeller <- function(plot, x, centrality.line.args, ...) {
  # compute centrality measure (with a temporary data frame)
  df_central <- suppressWarnings(centrality_description(tibble(.x = ".x", "var" = x), .x, var, ...))
  
  # adding a vertical line corresponding to centrality parameter
  plot +
    exec(geom_vline, xintercept = df_central$var, !!!centrality.line.args) +
    scale_x_continuous(
      sec.axis = sec_axis(
        trans = ~.,
        name = NULL,
        labels = parse(text = df_central$expression),
        breaks = df_central$var
      )
    )
}

.eval_f <- function(.f, ...) {
  tryCatch(
    suppressWarnings(suppressMessages(exec(.f, ...))),
    error = function(e) NULL
  )
}

.palette_message <- function(package, palette, min_length) {
  palette_length <- paletteer::palettes_d_names %>%
    filter(package == !!package, palette == !!palette) %$%
    length[[1L]]
  
  are_enough_colors_available <- palette_length > min_length
  
  if (!are_enough_colors_available) {
    rlang::inform(c(
      "Number of labels is greater than default palette color count.",
      "Select another color `palette` (and/or `package`)."
    ))
  }
  
  return(are_enough_colors_available)
}

ggbetweenstats <- function(data,
                           x,
                           y,
                           fill_var = NULL,
                           type = "parametric",
                           pairwise.comparisons = TRUE,
                           pairwise.display = "significant",
                           p.adjust.method = "holm",
                           effsize.type = "unbiased",
                           bf.prior = 0.707,
                           bf.message = TRUE,
                           results.subtitle = TRUE,
                           xlab = NULL,
                           ylab = NULL,
                           caption = NULL,
                           title = NULL,
                           subtitle = NULL,
                           k = 2L,
                           var.equal = FALSE,
                           conf.level = 0.95,
                           nboot = 100L,
                           tr = 0.2,
                           centrality.plotting = TRUE,
                           centrality.type = type,
                           centrality.point.args = list(size = 5, color = "darkred"),
                           centrality.label.args = list(
                             size = 3,
                             nudge_x = 0.4,
                             segment.linetype = 4,
                             min.segment.length = 0
                           ),
                           point.args = list(
                             position = ggplot2::position_jitterdodge(dodge.width = 0.60),
                             alpha = 0.4,
                             size = 3,
                             stroke = 0,
                             na.rm = TRUE
                           ),
                           boxplot.args = list(width = 0.3, alpha = 0.2, na.rm = TRUE),
                           violin.args = list(width = 0.5, alpha = 0.2, na.rm = TRUE),
                           ggsignif.args = list(textsize = 3, tip_length = 0.01, na.rm = TRUE),
                           ggtheme = ggstatsplot::theme_ggstatsplot(),
                           package = "RColorBrewer",
                           palette = "Dark2",
                           ggplot.component = NULL,
                           ...) {
  
  # data -----------------------------------
  
  # convert entered stats type to a standard notation
  type <- stats_type_switch(type)
  
  # make sure both quoted and unquoted arguments are allowed
  c(x, y) %<-% c(ensym(x), ensym(y))
  
  data %<>%
    select({{ x }}, {{ y }}) %>%
    tidyr::drop_na() %>%
    mutate({{ x }} := droplevels(as.factor({{ x }})))
  
  # statistical analysis ------------------------------------------
  
  # test to run; depends on the no. of levels of the independent variable
  test <- ifelse(nlevels(data %>% pull({{ x }})) < 3L, "t", "anova")
  
  if (results.subtitle) {
    # relevant arguments for statistical tests
    .f.args <- list(
      data         = data,
      x            = as_string(x),
      y            = as_string(y),
      effsize.type = effsize.type,
      conf.level   = conf.level,
      var.equal    = var.equal,
      k            = k,
      tr           = tr,
      paired       = FALSE,
      bf.prior     = bf.prior,
      nboot        = nboot
    )
    
    .f <- .f_switch(test)
    subtitle_df <- .eval_f(.f, !!!.f.args, type = type)
    subtitle <- if (!is.null(subtitle_df)) subtitle_df$expression[[1L]]
    
    if (type == "parametric" && bf.message) {
      caption_df <- .eval_f(.f, !!!.f.args, type = "bayes")
      caption <- if (!is.null(caption_df)) caption_df$expression[[1L]]
    }
  }
  
  # plot -----------------------------------
   # Adding fill_var to the ggplot aesthetic if it's not NULL
  if (!is.null(fill_var)) {
    fill_var <- enexpr(fill_var)
    plot <- ggplot(data, mapping = aes({{ x }}, {{ y }}, fill = {{ fill_var }})) +
      exec(geom_point, aes(color = {{ x }}), !!!point.args) +
      exec(geom_boxplot, aes(fill = {{ fill_var }}), outlier.shape = NA, !!!boxplot.args) +
      exec(geom_violin, !!!violin.args)
  } else {
    plot <- ggplot(data, mapping = aes({{ x }}, {{ y }})) +
      exec(geom_point, aes(color = {{ x }}), !!!point.args) +
      exec(geom_boxplot, outlier.shape = NA, !!!boxplot.args) +
      exec(geom_violin, !!!violin.args)
  }

  # centrality tagging -------------------------------------
  
  if (isTRUE(centrality.plotting)) {
    plot <- suppressWarnings(.centrality_ggrepel(
      plot                  = plot,
      data                  = data,
      x                     = {{ x }},
      y                     = {{ y }},
      k                     = k,
      type                  = stats_type_switch(centrality.type),
      tr                    = tr,
      centrality.point.args = centrality.point.args,
      centrality.label.args = centrality.label.args
    ))
  }
  
  # ggsignif labels -------------------------------------
  
  seclabel <- NULL
  
  if (isTRUE(pairwise.comparisons) && test == "anova") {
    mpc_df <- pairwise_comparisons(
      data            = data,
      x               = {{ x }},
      y               = {{ y }},
      type            = type,
      tr              = tr,
      paired          = FALSE,
      var.equal       = var.equal,
      p.adjust.method = p.adjust.method,
      k               = k
    )
    
    # adding the layer for pairwise comparisons
    plot <- .ggsignif_adder(
      plot             = plot,
      mpc_df           = mpc_df,
      data             = data,
      x                = {{ x }},
      y                = {{ y }},
      pairwise.display = pairwise.display,
      ggsignif.args    = ggsignif.args
    )
    
    # preparing the secondary label axis to give pairwise comparisons test details
    seclabel <- .pairwise_seclabel(
      unique(mpc_df$test),
      ifelse(type == "bayes", "all", pairwise.display)
    )
  }
  
  # annotations ------------------------
  
  .aesthetic_addon(
    plot             = plot,
    x                = data %>% pull({{ x }}),
    xlab             = xlab %||% as_name(x),
    ylab             = ylab %||% as_name(y),
    title            = title,
    subtitle         = subtitle,
    caption          = caption,
    seclabel         = seclabel,
    ggtheme          = ggtheme,
    package          = package,
    palette          = palette,
    ggplot.component = ggplot.component
  )
}


#' @title Violin plots for group or condition comparisons in between-subjects
#'   designs repeated across all levels of a grouping variable.
#' @name grouped_ggbetweenstats
#'
#' @description
#'
#' Helper function for `ggstatsplot::ggbetweenstats` to apply this function
#' across multiple levels of a given factor and combining the resulting plots
#' using `ggstatsplot::combine_plots`.
#'
#' @inheritParams ggbetweenstats
#' @inheritParams .grouped_list
#' @inheritParams combine_plots
#' @inheritDotParams ggbetweenstats -title
#'
#' @seealso \code{\link{ggbetweenstats}}, \code{\link{ggwithinstats}},
#'  \code{\link{grouped_ggwithinstats}}
#'
#' @inherit ggbetweenstats return references
#'
#' @examplesIf identical(Sys.getenv("NOT_CRAN"), "true") && requireNamespace("PMCMRplus", quietly = TRUE)
#' # for reproducibility
#' set.seed(123)
#' library(PMCMRplus) # for pairwise comparisons
#' library(dplyr, warn.conflicts = FALSE)
#' library(ggplot2)
#'
#' # the most basic function call
#' grouped_ggbetweenstats(
#'   data = filter(ggplot2::mpg, drv != "4"),
#'   x = year,
#'   y = hwy,
#'   grouping.var = drv
#' )
#'
#' # modifying individual plots using `ggplot.component` argument
#' grouped_ggbetweenstats(
#'   data = filter(
#'     movies_long,
#'     genre %in% c("Action", "Comedy"),
#'     mpaa %in% c("R", "PG")
#'   ),
#'   x = genre,
#'   y = rating,
#'   grouping.var = mpaa,
#'   ggplot.component = scale_y_continuous(
#'     breaks = seq(1, 9, 1),
#'     limits = (c(1, 9))
#'   )
#' )
#' @export
grouped_ggbetweenstats <- function(data,
                                   ...,
                                   grouping.var,
                                   plotgrid.args = list(),
                                   annotation.args = list()) {
  purrr::pmap(
    .l = .grouped_list(data, {{ grouping.var }}),
    .f = ggbetweenstats,
    ...
  ) %>%
    combine_plots(plotgrid.args, annotation.args)
}
