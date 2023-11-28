library(ggplot2)
library(classInt)

ggplot_swatch <- function(df, aes_fill, countries, palette = "Blue-Red 3", breaks_style = "quantile", n_col = 10) {
  
  # Compute the colors
  cols <- hcl.colors(n = n_col, palette = palette, rev = FALSE)
  
  # Compute the breaks
  breaks <- classInt::classIntervals(df[[aes_fill]], n = n_col, style = breaks_style)$brks
  
  # Cut the fill variable according to the breaks
  df$fill_factor <- cut(df[[aes_fill]], breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  # Generate the plot
  p <- ggplot() +
    geom_tile(data = df, aes(x = X, y = Y, fill = as.factor(fill_factor))) +
    geom_sf(data = countries) +
    scale_y_continuous(breaks = seq(from = 40, to = 48, by = 2)) +
    scale_fill_manual(values = cols, name = aes_fill) +
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      pad_x = unit(1.5, "cm"),
      pad_y = unit(1, "cm"),
      style = north_arrow_fancy_orienteering,
      height = unit(1.5, "cm"),
      width = unit(1.5, "cm")
    ) +
    annotation_scale(
      location = "bl",
      width_hint = 0.3,
      text_cex = 0.6,
      pad_y = unit(0.08,'cm'),
      pad_x = unit(0.5,'cm')
    ) +
    theme_bw() +
    theme(axis.title = element_blank(), legend.title = element_text(size = 10)) +
    ggtitle("Distribution of Weighted Endemism of Fish in the Black Sea")
  
  return(p)
}

# Example usage
 ggplot_swatch(df, "WE", countries)
