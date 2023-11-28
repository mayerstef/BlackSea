library(ggplot2)
library(ggthemes)

bio_indices <- bio_indices %>%
  mutate(code = case_when(
    str_detect(Bioregions, "Turkey_bs1") ~ "TR (W)",
    str_detect(Bioregions, "Turkey_bs2") ~ "TR (E)",
    str_detect(Bioregions, "Turkey_marmara") ~ "Marmara (TR)",
    str_detect(Bioregions, "Ukraine_azov") ~ "Azov (UA)",
    str_detect(Bioregions, "Ukraine") ~ "UA",
    str_detect(Bioregions, "Russia_azov") ~ "Azov (RU)",
    str_detect(Bioregions, "Russia") ~ "RU",
    str_detect(Bioregions, "Georgia") ~ "GE",
    str_detect(Bioregions, "Bulgaria") ~ "BG",
    str_detect(Bioregions, "Romania") ~ "RO"
  ))


selected_vars1 <- c(
  "fricSES",
  "fdisSES",
  "fmpdSES",
  "fnndSES",
  "feveSES",
  "fdivSES",
  "foriSES",
  "fspeSES"
)

selected_vars2 <- c("NTI",
                    "NRI",
                    "PD.SES",
                    "NTI.SES",
                    "NRI.SES",
                    "MPD.SES",
                    "MNTD.SES",
                    "VPD.SES",
                    "VNTD.SES",
                    "zscore_ts",
                    "zscore_rw")

# Create a folder for selected_vars1 plots
dir.create("selected_vars1_plots")

# Loop over selected_vars1
for (var in selected_vars1) {
  indice <- bio_indices[[var]]
  region <- bio_indices$code
  kruskal <- kruskal.test(indice ~ factor(region), data = bio_indices)
  posthoc_dunn <- dunnTest(indice ~ factor(region), data = bio_indices, method = "bonferroni")
  
  # Create boxplot
  plot1 <- ggplot(bio_indices, aes(x = region, y = indice, fill = Countries)) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Greens") +
    theme_ggstatsplot() +
    labs(title = paste("Functional diversity", var, "among Bioregions"),
         x = "Bioregions",
         y = paste("Functional diversity", var),
         fill = "Bioregions") +
    theme(legend.position = "none")
  
  # Save the plot
  plot_path <- file.path("selected_vars1_plots", paste(var, ".png", sep = ""))
  ggsave(plot_path, plot = plot1, width = 6, height = 4, dpi = 300)
}

# Create a folder for selected_vars2 plots
dir.create("selected_vars2_plots")

# Loop over selected_vars2
for (var in selected_vars2) {
  indice <- bio_indices[[var]]
  region <- bio_indices$code
  kruskal <- kruskal.test(indice ~ factor(region), data = bio_indices)
  posthoc_dunn <- dunnTest(indice ~ factor(region), data = bio_indices, method = "bonferroni")
  
  # Create boxplot
  plot2 <- ggplot(bio_indices, aes(x = region, y = indice, fill = Countries)) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Reds") +
    theme_ggstatsplot()+
    labs(title = paste("Phylogenetic Diversity", var, "among Bioregions"),
         x = "Bioregions",
         y = paste("Phylogenetic Diversity", var),
         fill = "Bioregions") +
    theme(legend.position = "none")
  
  # Save the plot
  plot_path <- file.path("selected_vars2_plots", paste(var, ".png", sep = ""))
  ggsave(plot_path, plot = plot2, width = 6, height = 4, dpi = 300)
}


# Create a folder for selected_vars1 plots
dir.create("selected_vars3_plots")

# Loop over selected_vars1
for (var in selected_vars1) {
  indice <- bio_indices[[var]]
  region <- bio_indices$code
  kruskal <- kruskal.test(indice ~ factor(region), data = bio_indices)
  posthoc_dunn <- dunnTest(indice ~ factor(region), data = bio_indices, method = "bonferroni")
  
  # Create a column for significance level based on the adjusted p-value
  posthoc_dunn$res$significance <- ifelse(posthoc_dunn$res$`P.adj` < 0.05, "***",
                                          ifelse(posthoc_dunn$res$`P.adj` < 0.01, "**",
                                                 ifelse(posthoc_dunn$res$`P.adj` < 0.1, "*", "ns")))
  
  # Create bar plot for the Z-values with significance level as labels
  plot1 <- ggplot(posthoc_dunn$res, aes(x = Z, y = reorder(Comparison, -Z), fill = significance)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = significance), vjust = -0.5, size = 3.5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Z-value") +
    ylab("Comparison") +
    labs(fill = "Significance") +
    scale_fill_brewer(palette = "Greens") +
    theme_ggstatsplot()
  
  # Save the plot
  plot_path <- file.path("selected_vars3_plots", paste(var, ".png", sep = ""))
  ggsave(plot_path, plot = plot1, width = 6, height = 8, dpi = 300)
}

# Create a folder for selected_vars2 plots
dir.create("selected_vars4_plots")

# Loop over selected_vars2
for (var in selected_vars2) {
  indice <- bio_indices[[var]]
  region <- bio_indices$code
  kruskal <- kruskal.test(indice ~ factor(region), data = bio_indices)
  posthoc_dunn <- dunnTest(indice ~ factor(region), data = bio_indices, method = "bonferroni")
  
  # Create a column for significance level based on the adjusted p-value
  posthoc_dunn$res$significance <- ifelse(posthoc_dunn$res$`P.adj` < 0.05, "***",
                                          ifelse(posthoc_dunn$res$`P.adj` < 0.01, "**",
                                                 ifelse(posthoc_dunn$res$`P.adj` < 0.1, "*", "ns")))
  
  # Create bar plot for the Z-values with significance level as labels
  plot2 <- ggplot(posthoc_dunn$res, aes(x = Z, y = reorder(Comparison, -Z), fill = significance)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = significance), vjust = -0.5, size = 3.5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Z-value") +
    ylab("Comparison") +
    labs(fill = "Significance") +
    scale_fill_brewer(palette = "Reds") +
    theme_ggstatsplot()
  
  # Save the plot
  plot_path <- file.path("selected_vars4_plots", paste(var, ".png", sep = ""))
  ggsave(plot_path, plot = plot2, width = 6, height = 8, dpi = 300)
}



library(FSA)

load("./Data/bio_indices.RData")

# Nonparametric tests
region <- bio_indices$Bioregions

# Function to perform Kruskal-Wallis and Dunn's tests for a variable
perform_tests <- function(variable) {
  indice <- bio_indices[[variable]]
  
  cat("Variable:", variable, "\n")
  
  # Kruskal-Wallis test
  cat("Kruskal-Wallis Test:\n")
  kruskal <- kruskal.test(indice ~ factor(region), data = bio_indices)
  print(kruskal)
  
  # Post hoc Dunn's test
  cat("Dunn's Test (Bonferroni-adjusted):\n")
  posthoc_dunn <- dunnTest(indice ~ factor(region), data = bio_indices, method = "bonferroni")
  print(posthoc_dunn)
  
  cat("\n")
  
  # Store results in a list
  results <- list(
    kruskal = kruskal,
    posthoc_dunn = posthoc_dunn
  )
  
  return(results)
}

# Create an empty list to store results
results_list <- list()

# Iterate over selected variables and perform tests
for (var in selected_vars) {
  results <- perform_tests(var)
  results_list[[var]] <- results
}

# Access and print results for a specific variable
variable <- "SR"
print(results_list[[variable]])



library(FSA)
library(ggpubr)
library(ggstatsplot)
library(statsExpressions)
library(stringr)
library(rlang)
library(glue)

load("./Data/bio_indices.RData")

# Function to generate boxplot and save the plot
generate_boxplot <- function(variable) {
  cat("Generating plot for variable:", variable, "\n")
  
  plot_title <- paste("Boxplot -", variable)
  
  bioregions <- ggbetweenstats(
    data = bio_indices,
    x = Bioregions,
    y = !!sym(variable),
    fill = bio_indices$Countries,
    plot.type = "box",
    type = "non-parametric",
    ggtheme = ggpubr::theme_pubclean(),
    boxplot.args = list(color = "black", alpha = 1, width = 0.7),
    pairwise.comparisons = TRUE,
    pairwise.display = "ns",
    p.adjust.method = "bonferroni",
    xlab = "Bioregions",
    ylab = variable,
    point.args = list(alpha = 0),
    violin.args = list(alpha = 0, width = 0),
    centrality.type = "parametric",
    centrality.label.args = list(
      size = 2.5,
      fill = "white",
      label.r = .5,
      nudge_x = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, -.3, -.3),
      nudge_y = c(0, 2, 0, 0, -1, 0, 2, 0, -2, -2),
      segment.linetype = 4,
      min.segment.length = 0
    ),
    centrality.point.args = list(size = 3, color = "darkred"),
    ggsignif.args = list(textsize = 3, tip_length = 0.01, na.rm = TRUE),
    messages = FALSE,
    palette="Blues",
  )
  
  # Create directory if it doesn't exist
  if (!file.exists("./AnovaPlots")) {
    dir.create("./AnovaPlots")
  }
  
  # Save the plot
  plot_filename <- paste0("./AnovaPlots/", variable, "_boxplot_pubclean.png")
  ggsave(plot_filename, bioregions, width = 10, height = 6.5, dpi = 300, bg = "transparent")
  
  cat("Plot saved as:", plot_filename, "\n\n")
}

# Iterate over selected variables and generate boxplot for each variable
for (var in selected_vars) {
  generate_boxplot(var)
}

