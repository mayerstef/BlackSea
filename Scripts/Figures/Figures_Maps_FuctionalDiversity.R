#Select Indice to Map
indice<-bio_indices$SR
(map <- ggplot() + 
    geom_tile(data = bio_indices, aes(x = X, y = Y, fill = indice))+
    geom_sf(data = st_geometry(Bioreg), fill = NA) +
    geom_sf_label(aes(label = code), 
                  data=Bioreg, 
                  color="#001219ff",
                  label.size  = NA, 
                  label.r = unit(0.5, "lines"), 
                  alpha = 0.5)+
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
    #geom_sf(data = countries, colour = "black", fill="#D7D7D7")+ 
    scale_fill_gradientn(colours = gradient, 
                         breaks = c(5, 50, 95, 140), #this is for looking at SR
                         name="SR") + #Don't forget to change the legend name
    annotation_north_arrow(
      location = "tl",
      which_north = "true",
      style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "br") +
    theme_pubr()+
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(), 
          legend.title = element_text(size = 10, colour = "black"),
          legend.position = c(0.94, .85), 
          legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black")))
ggsave("./Maps3/SR.png", map, width = 10, height = 6.5, dpi = 300, bg = "transparent")

################################################################################
# Map for "fdis"
indice <- bio_indices$fric
(map_fric <- ggplot() + 
  geom_tile(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
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
  scale_fill_gradientn(colours = gradient, name = "FRic") +
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
ggsave("./Maps3/FRic.png", map_fric, width = 10, height = 6.5, dpi = 300, bg = "transparent")

indice<-bio_indices$fricSES
(map_fric <- ggplot() + 
    geom_tile(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
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
    scale_fill_gradientn(colours = gradient, 
                         breaks = c(-1.96, -1.28, 0, 1.28,1.96), #this is for looking at SR
                         name = "FRic (SES)") +
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
ggsave("./Maps3/FRicSES.png", map_fric, width = 10, height = 6.5, dpi = 300, bg = "transparent")



indice<-bio_indices$fricSES
(map_fric <- ggplot() + 
    geom_tile(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
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
    scale_fill_gradientn(colours = gradient, 
                         breaks = c(-1.96, -1.28, 0, 1.28,1.96), #this is for looking at SR
                         name = "FRic (SES)") +
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
ggsave("./Maps3/FRicSES.png", map_fric, width = 10, height = 6.5, dpi = 300, bg = "transparent")

pal <- c(
  "#001219ff", 
  "#003946",
  "#005f73ff", 
  #"#057985",
# "#0e9396",
 "#0e9396",
"#94d2bdff",
  "#94d2bdff",
  "#C1E5D9",
  "#f1f1de",
  "#e9d8a6",
  "#ecba53",
  "#ee9b00ff", 
  "#ca6702",
  "#bb3e03ff",
  #"#ae2012ff",
  "#9b2226",
  "#340A00")
# Map for "fdisSES"
indice <- bio_indices$fdisSES
(map_fdisSES <- ggplot() + 
    geom_tile(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
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
                         name = "FDis (SES)") +
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
map_fdisSES
ggsave("./Maps3/fdisSES.png", map_fdisSES, width = 10, height = 6.5, dpi = 300, bg = "transparent")



pal <- c(
  "#001219ff", 
  #"#003946",
  "#005f73ff", 
  #"#057985",
  # "#0e9396",
  "#0e9396",
  "#94d2bdff",
  "#C1E5D9",
  "#f1f1de",
  "#e9d8a6",
  "#ecba53",
  "#ee9b00ff", 
 # "#ca6702",
  "#bb3e03ff",
  "#ae2012ff",
  #"#9b2226",
  "#340A00",
  "#340A00")
indice <- bio_indices$fmpdSES
(map_fmpdSES <- ggplot() + 
    geom_tile(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
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
                         name = "FMPD (SES)") +
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
ggsave("./Maps3/fmpdSES.png", map_fmpdSES, width = 10, height = 6.5, dpi = 300, bg = "transparent")


# Map for "fnndSES"
indice <- bio_indices$fnndSES
(map_fnndSES <- ggplot() + 
    geom_tile(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
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
    scale_fill_gradientn(colours = gradient, 
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "FNND (SES)") +
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
ggsave("./Maps3/fnndSES.png", map_fnndSES, width = 10, height = 6.5, dpi = 300, bg = "transparent")

# Repeat the same structure for the remaining variables...

# Map for "feveSES"
indice <- bio_indices$feveSES
(map_feveSES <- ggplot() + 
    geom_tile(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
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
    scale_fill_gradientn(colours = gradient, 
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "FEve (SES)") +
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
ggsave("./Maps3/feveSES.png", map_feveSES, width = 10, height = 6.5, dpi = 300, bg = "transparent")

# Repeat the same structure for the remaining variables...

# Map for "fdivSES"
indice <- bio_indices$fdivSES
(map_fdivSES <- ggplot() + 
    geom_tile(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
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
    scale_fill_gradientn(colours = gradient, 
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "FDiv (SES)") +
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
ggsave("./Maps3/fdivSES.png", map_fdivSES, width = 10, height = 6.5, dpi = 300, bg = "transparent")



# Map for "foriSES"
indice <- bio_indices$foriSES
(map_foriSES <- ggplot() + 
    geom_tile(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
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
                         name = "FOri (SES)") +
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
map_foriSES
ggsave("./Maps3/foriSES.png", map_foriSES, width = 10, height = 6.5, dpi = 300, bg = "transparent")


pal <- c(
  "#001219ff", 
  "#003946",
  #"#005f73ff", 
  #"#057985",
  # "#0e9396",
  "#0e9396",
  "#0e9396",
  "#94d2bdff",
  "#C1E5D9",
  "#f1f1de",
  "#e9d8a6",
  "#ecba53",
  "#ee9b00ff", 
  # "#ca6702",
  "#bb3e03ff",
  "#ae2012ff",
  "#340A00",
  #"#9b2226",
  "#340A00",
  "#340A00")

indice <- bio_indices$fspeSES
(map_fspeSES <- ggplot() +
    geom_tile(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
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
    scale_fill_gradientn(colours = pal,#gradient_pal(200),
                         breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                         name = "FSpe (SES)") +
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
map_fspeSES
ggsave("./Maps3/fspeSES.png", map_fspeSES, width = 10, height = 6.5, dpi = 300, bg = "transparent")


indice <- bio_indices$fmpdSES
(map_fmpd <- ggplot() + 
  geom_tile(data = bio_indices, aes(x = X, y = Y, fill = indice)) +
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
  scale_fill_gradientn(colours = gradient, 
                       breaks = c(-1.96, -1.28, 0, 1.28, 1.96),
                       name = "FMPD (SES)") +
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
ggsave("./Maps3/fmpd.png", map_fmpd, width = 10, height = 6.5, dpi = 300, bg = "transparent")



