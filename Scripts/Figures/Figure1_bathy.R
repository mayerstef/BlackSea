Bioreg<-st_read("./Data/GIS/10bioregions_blacksea.shp") #10 features, see script "BlackSea_Bioregions" for shapefile creation
Bioregion_fished <- st_read("Data/GIS/fishedEEZ_blacksea.shp")
NOFISH_row <- Bioregion_fished[Bioregion_fished$UNION == "NOFISH",]
Bioreg2 <- bind_rows(Bioreg, NOFISH_row)
Bioreg_NOFISH <- st_difference(Bioreg2, NOFISH_row)
Bioreg_NOFISH<-Bioreg_NOFISH[,-2]


Bioreg <- Bioreg %>%
  mutate(code = case_when(
    str_detect(UNION, "Turkey_bs1") ~ "TR (West)",
    str_detect(UNION, "Turkey_bs2") ~ "TR (East)",
    str_detect(UNION, "Turkey_marmara") ~ "TR (Marmara)",
    str_detect(UNION, "Ukraine_azov") ~ "UA (Azov)",
    str_detect(UNION, "Ukraine") ~ "UA",
    str_detect(UNION, "Russia_azov") ~ "RU (Azov)",
    str_detect(UNION, "Russia") ~ "RU",
    str_detect(UNION, "Georgia") ~ "GE",
    str_detect(UNION, "Bulgaria") ~ "BG",
    str_detect(UNION, "Romania") ~ "RO",
    TRUE ~ UNION # if none of the conditions are met, use UNION value
  ))


Bioreg_NOFISH <- Bioreg_NOFISH %>%
  mutate(code = case_when(
    str_detect(UNION, "Turkey_bs1") ~ "TR (West)",
    str_detect(UNION, "Turkey_bs2") ~ "TR (East)",
    str_detect(UNION, "Turkey_marmara") ~ "TR (Marmara)",
    str_detect(UNION, "Ukraine_azov") ~ "UA (Azov)",
    str_detect(UNION, "Ukraine") ~ "UA",
    str_detect(UNION, "Russia_azov") ~ "RU (Azov)",
    str_detect(UNION, "Russia") ~ "RU",
    str_detect(UNION, "Georgia") ~ "GE",
    str_detect(UNION, "Bulgaria") ~ "BG",
    str_detect(UNION, "Romania") ~ "RO",
    TRUE ~ UNION # if none of the condit
))

library(ggOceanMaps)
dt <- data.frame(lon = c(xmin, xmin, xmax,xmax), lat = c(ymin, ymax, ymin, ymax))
test<-basemap(data = dt, bathymetry = TRUE, #bathy.style = "contour_blues", 
              land.border.col = "transparent", land.col = "transparent",  grid.col = "transparent") + 
  #  theme(panel.background = element_rect(fill = "grey"),
  #       panel.ontop = FALSE) +
  scale_color_viridis_d(option = "G", end =0.7)
test
test+geom_sf(data = st_geometry(Bioreg), fill = NA, lwd=.2, color="black")+
  geom_sf_label(aes(label = code), 
                data=Bioreg, 
                color="#001219ff",
                size =3,
                #label.size  = NA,  
                label.r = unit(0.5, "lines"), 
                alpha = 1)+
  ggspatial::annotation_scale(location = "br") + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true") 

test+geom_sf(data = Bioreg2, fill = NA) +
  geom_sf(data = NOFISH_row, fill = "#ae2012ff", color = NA) +
  geom_sf_label(aes(label = code), 
                data=Bioreg_NOFISH, 
                color="#001219ff",
                size =3,
                #label.size  = NA, 
                label.r = unit(0.5, "lines"), 
                alpha = 1)+
 
  ggspatial::annotation_scale(location = "br") + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true") 

test+geom_sf(data = st_geometry(Bioreg2), fill = NA, lwd=.2, color="black")+
  geom_sf(data = NOFISH_row, fill = "#ae2012ff", color = NA, alpha =0.2) +
  geom_sf_label(aes(label = code), 
                data=Bioreg2, 
                color="#001219ff",
                size =3,
                #label.size  = NA, 
                label.r = unit(0.5, "lines"), 
                alpha = 1)+
  ggspatial::annotation_scale(location = "br") + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true") 


test + 
  geom_sf(data = st_geometry(Bioreg2), fill = NA, lwd = 0.2, color = "black") +
  geom_sf_label(
    aes(label = code), 
    data = Bioreg2, 
    color = "#001219ff",
    size = 3,
    label.r = unit(0.5, "lines"), 
    alpha = 1
  ) +
  geom_sf_pattern(
    data = NOFISH_row,
    pattern = "stripe",
    pattern_angle = 45,
    pattern_density = 0.01,
    #pattern_spacing = 0.1,
    pattern_key_scale_factor = 0.2, 
    pattern_alpha=0.7
  ) +
  geom_sf_label(
    aes(label = code), 
    data = Bioreg2, 
    color = "#001219ff",
    size = 3,
    label.r = unit(0.5, "lines"), 
    alpha = 1
  ) +
  ggspatial::annotation_scale(location = "br") + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true")
