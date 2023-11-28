indice <- bio_indices$PD
(map_PD_SES <- ggplot() + 
    geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
    geom_sf(data = st_geometry(Bioreg), fill = NA) +
    geom_sf_label(aes(label = code), 
                  data = Bioreg, 
                  color = "#001219ff",
                  label.size = NA, 
                  label.r = unit(0.5, "lines"), 
                  alpha = 0.5) +
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
    geom_sf(data = countries, colour = "black", fill = "#D7D7D7") + 
    scale_fill_gradientn(colours = gradient, 
                         #breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "PD") +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "br") +
    theme_pubr() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(), 
          legend.title = element_text(size = 10, colour = "black"),
          legend.position = c(0.94, .85), 
          legend.background = element_rect(fill = "white",
                                           size = 0.5, linetype = "solid", 
                                           colour = "black")))
ggsave("./Maps/PD.png", map_PD_SES, width = 10, height = 6.5, dpi = 300, bg = "transparent")


indice <- bio_indices$PD.SES
(map_PD_SES <- ggplot() + 
    geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
    geom_sf(data = st_geometry(Bioreg), fill = NA) +
    geom_sf_label(aes(label = code), 
                  data = Bioreg, 
                  color = "#001219ff",
                  label.size = NA, 
                  label.r = unit(0.5, "lines"), 
                  alpha = 0.5) +
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
    geom_sf(data = countries, colour = "black", fill = "#D7D7D7") + 
    scale_fill_gradientn(colours = gradient, 
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "PD (SES TS)") +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "br") +
    theme_pubr() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(), 
          legend.title = element_text(size = 10, colour = "black"),
          legend.position = c(0.94, .85), 
          legend.background = element_rect(fill = "white",
                                           size = 0.5, linetype = "solid", 
                                           colour = "black")))
ggsave("./Maps/PD_SES.png", map_PD_SES, width = 10, height = 6.5, dpi = 300, bg = "transparent")

# Map for "NTI"
indice <- bio_indices$NTI
(map_NTI <- ggplot() + 
    geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
    geom_sf(data = st_geometry(Bioreg), fill = NA) +
    geom_sf_label(aes(label = code), 
                  data = Bioreg, 
                  color = "#001219ff",
                  label.size = NA, 
                  label.r = unit(0.5, "lines"), 
                  alpha = 0.5) +
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
    geom_sf(data = countries, colour = "black", fill = "#D7D7D7") + 
    scale_fill_gradientn(colours = gradient_pal(200), 
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "NTI") +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "br") +
    theme_pubr() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(), 
          legend.title = element_text(size = 10, colour = "black"),
          legend.position = c(0.94, .85), 
          legend.background = element_rect(fill = "white",
                                           size = 0.5, linetype = "solid", 
                                           colour = "black")))
ggsave("./Maps/NTI.png", map_NTI, width = 10, height = 6.5, dpi = 300, bg = "transparent")

# Repeat the same structure for the remaining variables...

# Map for "NTI.SES"
indice <- bio_indices$NTI.SES
(map_NTI_SES <- ggplot() + 
    geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
    geom_sf(data = st_geometry(Bioreg), fill = NA) +
    geom_sf_label(aes(label = code), 
                  data = Bioreg, 
                  color = "#001219ff",
                  label.size = NA, 
                  label.r = unit(0.5, "lines"), 
                  alpha = 0.5) +
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
    geom_sf(data = countries, colour = "black", fill = "#D7D7D7") + 
    scale_fill_gradientn(colours = gradient, #gradient_pal(200), 
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "NTI (SES TS)") +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "br") +
    theme_pubr() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(), 
          legend.title = element_text(size = 10, colour = "black"),
          legend.position = c(0.94, .85), 
          legend.background = element_rect(fill = "white",
                                           size = 0.5, linetype = "solid", 
                                           colour = "black")))
ggsave("./Maps/NTI_SES.png", map_NTI_SES, width = 10, height = 6.5, dpi = 300, bg = "transparent")

# Repeat the same structure for the remaining variables...

pal <- c(
  "#001219ff", 
  "#003946",
  "#005f73ff", 
  "#057985",
  "#0e9396",
  "#94d2bdff",
  "#C1E5D9",
  "#f1f1de",
  "#e9d8a6",
  "#ecba53",
  "#ee9b00ff", 
  "#ca6702",
  "#bb3e03ff",
  "#ae2012ff",
  "#9b2226",
  "#340A00")


# Map for "NRI"
indice <- bio_indices$NRI
(map_NRI <- ggplot() + 
    geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
    geom_sf(data = st_geometry(Bioreg), fill = NA) +
    geom_sf_label(aes(label = code), 
                  data = Bioreg, 
                  color = "#001219ff",
                  label.size = NA, 
                  label.r = unit(0.5, "lines"), 
                  alpha = 0.5) +
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
    geom_sf(data = countries, colour = "black", fill = "#D7D7D7") + 
    scale_fill_gradientn(colours = pal, 
                         #colours = gradient, 
                         breaks = c(-5, -1.96, -1.28, 0, 1.28, 1.96),
                         name = "NRI") +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "br") +
    theme_pubr() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(), 
          legend.title = element_text(size = 10, colour = "black"),
          legend.position = c(0.94, .85), 
          legend.background = element_rect(fill = "white",
                                           size = 0.5, linetype = "solid", 
                                           colour = "black")))
ggsave("./Maps/NRI.png", map_NRI, width = 10, height = 6.5, dpi = 300, bg = "transparent")

pal <- c(
  "#001219ff", 
  "#003946",
  "#005f73ff", 
  "#057985",
  "#0e9396",
  "#94d2bdff",
  "#C1E5D9",
  "#f1f1de",
  "#e9d8a6",
  "#ecba53",
  "#ee9b00ff", 
  #"#ca6702",
  #"#bb3e03ff",
  #"#ae2012ff",
  "#9b2226")
 # "#340A00")

# Map for "NRI.SES"
indice <- bio_indices$NRI.SES
(map_NRI_SES <- ggplot() + 
    geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
    geom_sf(data = st_geometry(Bioreg), fill = NA) +
    geom_sf_label(aes(label = code), 
                  data = Bioreg, 
                  color = "#001219ff",
                  label.size = NA, 
                  label.r = unit(0.5, "lines"), 
                  alpha = 0.5) +
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
    geom_sf(data = countries, colour = "black", fill = "#D7D7D7") + 
    scale_fill_gradientn(colours = pal, #colours = gradient, 
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "NRI (SES TS)") +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "br") +
    theme_pubr() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(), 
          legend.title = element_text(size = 10, colour = "black"),
          legend.position = c(0.94, .85), 
          legend.background = element_rect(fill = "white",
                                           size = 0.5, linetype = "solid", 
                                           colour = "black")))
ggsave("./Maps/NRI_SES.png", map_NRI_SES, width = 10, height = 6.5, dpi = 300, bg = "transparent")

# Map for "MPD.SES"
indice <- bio_indices$MPD.SES
(map_MPD_SES <- ggplot() + 
    geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
    geom_sf(data = st_geometry(Bioreg), fill = NA) +
    geom_sf_label(aes(label = code), 
                  data = Bioreg, 
                  color = "#001219ff",
                  label.size = NA, 
                  label.r = unit(0.5, "lines"), 
                  alpha = 0.5) +
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
    geom_sf(data = countries, colour = "black", fill = "#D7D7D7") + 
    scale_fill_gradientn(colours = gradient_pal(200), 
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "MPD (SES)") +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "br") +
    theme_pubr() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(), 
          legend.title = element_text(size = 10, colour = "black"),
          legend.position = c(0.94, .85), 
          legend.background = element_rect(fill = "white",
                                           size = 0.5, linetype = "solid", 
                                           colour = "black")))
ggsave("./Maps/MPD_SES.png", map_MPD_SES, width = 10, height = 6.5, dpi = 300, bg = "transparent")

# Repeat the same structure for the remaining variables...

# Map for "MNTD.SES"
indice <- bio_indices$MNTD.SES
(map_MNTD_SES <- ggplot() + 
    geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
    geom_sf(data = st_geometry(Bioreg), fill = NA) +
    geom_sf_label(aes(label = code), 
                  data = Bioreg, 
                  color = "#001219ff",
                  label.size = NA, 
                  label.r = unit(0.5, "lines"), 
                  alpha = 0.5) +
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
    geom_sf(data = countries, colour = "black", fill = "#D7D7D7") + 
    scale_fill_gradientn(colours = gradient_pal(200), 
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "MNTD (SES)") +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "br") +
    theme_pubr() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(), 
          legend.title = element_text(size = 10, colour = "black"),
          legend.position = c(0.94, .85), 
          legend.background = element_rect(fill = "white",
                                           size = 0.5, linetype = "solid", 
                                           colour = "black")))
ggsave("./Maps/MNTD_SES.png", map_MNTD_SES, width = 10, height = 6.5, dpi = 300, bg = "transparent")

# Repeat the same structure for the remaining variables...
pal <- c(
  "#001219ff", 
  "#003946",
  "#005f73ff", 
  "#057985",
  "#0e9396",
  "#94d2bdff",
  "#C1E5D9",
  "#f1f1de",
  "#e9d8a6",
  "#ecba53",
  "#ee9b00ff", 
  #"#ca6702",
  "#bb3e03ff",
  #"#ae2012ff",
  "#9b2226",
"#340A00")
# Map for "VPD.SES"
indice <- bio_indices$VPD.SES
(map_VPD_SES <- ggplot() + 
    geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
    geom_sf(data = st_geometry(Bioreg), fill = NA) +
    geom_sf_label(aes(label = code), 
                  data = Bioreg, 
                  color = "#001219ff",
                  label.size = NA, 
                  label.r = unit(0.5, "lines"), 
                  alpha = 0.5) +
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
    geom_sf(data = countries, colour = "black", fill = "#D7D7D7") + 
    scale_fill_gradientn(colours = pal, #colours = gradient, 
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.645, 1.96),
                         name = "VPD (SES)") +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "br") +
    theme_pubr() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(), 
          legend.title = element_text(size = 10, colour = "black"),
          legend.position = c(0.94, .85), 
          legend.background = element_rect(fill = "white",
                                           size = 0.5, linetype = "solid", 
                                           colour = "black")))
ggsave("./Maps/VPD_SES.png", map_VPD_SES, width = 10, height = 6.5, dpi = 300, bg = "transparent")

# Repeat the same structure for the remaining variables...

pal <- c(
  "#001219ff", 
  #"#003946",
  # "#005f73ff", 
  "#057985",
  #"#0e9396",
  "#94d2bdff",
  #"#C1E5D9",
  "#f1f1de",
  "#e9d8a6",
  "#ecba53",
  "#ee9b00ff", 
  "#ca6702",
  "#bb3e03ff",
  "#ae2012ff",
  "#9b2226",
  "#340A00")
# Map for "VNTD.SES"
indice <- bio_indices$VNTD.SES
(map_VNTD_SES <- ggplot() + 
    geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
    geom_sf(data = st_geometry(Bioreg), fill = NA) +
    geom_sf_label(aes(label = code), 
                  data = Bioreg, 
                  color = "#001219ff",
                  label.size = NA, 
                  label.r = unit(0.5, "lines"), 
                  alpha = 0.5) +
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
    geom_sf(data = countries, colour = "black", fill = "#D7D7D7") + 
    scale_fill_gradientn(colours = pal,#gradient, 
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "VNTD (SES)") +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "br") +
    theme_pubr() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(), 
          legend.title = element_text(size = 10, colour = "black"),
          legend.position = c(0.94, .85), 
          legend.background = element_rect(fill = "white",
                                           size = 0.5, linetype = "solid", 
                                           colour = "black")))
ggsave("./Maps/VNTD_SES.png", map_VNTD_SES, width = 10, height = 6.5, dpi = 300, bg = "transparent")

# Repeat the same structure for the remaining variables...
pal <- c(
 # "turquoise",
  "#001219ff", 
  "#001219ff", 
  #"#001219ff", 
  "#003946",
  #"#003946",
  "#003946",
  #"#005f73ff", 
  "#005f73ff", 
  #"#057985",
  "#0e9396",
  "#0e9396",
 "#0e9396",
 "#94d2bdff",
 #"#94d2bdff",
 # "#94d2bdff",
  "#C1E5D9",
  "#f1f1de",
 "#f1f1de",
 "#f1f1de",
  "#e9d8a6",
  "#ecba53",
  "#ee9b00ff", 
  "#ca6702",
  #"#bb3e03ff",
  "#ae2012ff",
  #"#9b2226",
  "#340A00")
  #"#340A00")
# Map for "zscore_ts"
indice <- bio_indices$zscore_ts
(map_zscore_ts <- ggplot() + 
    geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
    geom_sf(data = st_geometry(Bioreg), fill = NA) +
    geom_sf_label(aes(label = code), 
                  data = Bioreg, 
                  color = "#001219ff",
                  label.size = NA, 
                  label.r = unit(0.5, "lines"), 
                  alpha = 0.5) +
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
    #geom_sf(data = countries, colour = "black", fill = "#D7D7D7") + 
    scale_fill_gradientn(colours = pal,#gradient, 
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "PD (SES ts)") +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "br") +
    theme_pubr() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(), 
          legend.title = element_text(size = 10, colour = "black"),
          legend.position = c(0.94, .85), 
          legend.background = element_rect(fill = "white",
                                           size = 0.5, linetype = "solid", 
                                           colour = "black")))
ggsave("./Maps2/zscore_ts.png", map_zscore_ts, width = 10, height = 6.5, dpi = 300, bg = "transparent")

# Repeat the same structure for the remaining variables...

# Map for "zscore_rw"
indice <- bio_indices$zscore_rw
(map_zscore_rw <- ggplot() + 
    geom_sf(data = countries, colour = "black", fill = "#D7D7D7") + 
    geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
    geom_sf(data = st_geometry(Bioreg), fill = NA) +
    geom_sf_label(aes(label = code), 
                  data = Bioreg, 
                  color = "#001219ff",
                  label.size = NA, 
                  label.r = unit(0.5, "lines"), 
                  alpha = 0.5) +
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
      scale_fill_gradientn(colours = gradient, 
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "PD (SES rw") +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "br") +
    theme_pubr() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(), 
          legend.title = element_text(size = 10, colour = "black"),
          legend.position = c(0.94, .85), 
          legend.background = element_rect(fill = "white",
                                           size = 0.5, linetype = "solid", 
                                           colour = "black")))
ggsave("./Maps/zscore_rw.png", map_zscore_rw, width = 10, height = 6.5, dpi = 300, bg = "transparent")
