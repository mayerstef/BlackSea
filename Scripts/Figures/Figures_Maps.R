################################################################################
####                                                                        ####
####         Master thesis:LINKING PATTERNS IN PHYLOGENY, TRAITS,           ####
####            ABIOTIC VARIABlES AND SPACE: BLACK SEA FISHES               ####
####                                                                        ####
####   Stefanie Mayer                                                       ####
####                                                                        ####
####        calbouy@ethz.ch                                                 ####
################################################################################

library(ggstatsplot)
library(ggplot2)
library(ggpubr)

### Projection definition
proj84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") 
sf::sf_use_s2(F) #global option for using the S2 geometry library to FALSE


# Define the extent
xmin <- 26.223404293
xmax <- 41.923404293
ymin <- 39.879949463
ymax <- 47.379949463
cell_size<-.1

ncol <- ceiling((xmax - xmin) / cell_size)
nrow <- ceiling((ymax - ymin) / cell_size)

# Calculate the smallest centroid coordinates for each dimension
centre_x <- xmin + cell_size / 2
centre_y <- ymin + cell_size / 2

# Create the grid topology
grd <- sp::GridTopology(c(centre_x, centre_y), c(cell_size, cell_size), c(ncol, nrow))
grd <- raster(grd)
grd[] <- 0
grd@crs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

countries <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  st_geometry() %>%
  st_crop(grd)
countries <- st_as_sf(countries, crs = proj84) 

Bioreg<-st_read("./Data/GIS/10bioregions_blacksea.shp") #10 features, see script "BlackSea_Bioregions" for shapefile creation
Bioreg <- st_transform(Bioreg, proj84)

Bioreg <- Bioreg %>%
  mutate(code = case_when(
    grepl("Turkey", UNION) ~ "TR",
    grepl("Ukraine", UNION) ~ "UA",
    grepl("Russia", UNION) ~ "RU",
    grepl("Georgia", UNION) ~ "GE",
    grepl("Bulgaria", UNION) ~ "BG",
    grepl("Romania", UNION) ~ "RO",
    TRUE ~ NA_character_
  ))

#EEZ<-st_read("./Data/GIS/EEZ/EEZ_Land_v3_202030.shp")
#EEZ<- EEZ %>% 
#  st_geometry()%>%
#  st_crop(grd)



#Adjust colors selected based on diversity indice plotted
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

# Convert the color palette to a gradient
gradient_pal <- colorRampPalette(pal)
num_colors <- 50
gradient <- gradient_pal(num_colors)
print(gradient)

#Select Indice to Map
indice<-bio_indices$SR
(map <- ggplot() + 
  geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice))+
  geom_sf(data = st_geometry(Bioreg), fill = NA) +
  geom_sf_label(aes(label = code), 
                data=Bioreg, 
                color="#001219ff",
                label.size  = NA, 
                label.r = unit(0.5, "lines"), 
                alpha = 0.5)+
  scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
  geom_sf(data = countries, colour = "black", fill="#D7D7D7")+ 
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

map

ggsave("./Maps/SR.png", map, width = 10, height = 6.5, dpi = 300, bg = "transparent")

#For the maps below, EEZ and Bioregion labs are not plotted 
indice<-bio_indices$SES_FRic3
(map2 <- ggplot() + 
    geom_raster(data = bio_indices, aes(x = X, y = Y, fill = indice))+
    #geom_sf(data = st_geometry(Bioreg), fill = NA) +
    #geom_sf_label(aes(label = code), 
    #              data=Bioreg, 
    #              color="#001219ff",
    #              label.size  = NA, 
    #              label.r = unit(0.5, "lines"), 
    #              alpha = 0.5)+
    scale_y_continuous(breaks = seq(from = 39, to = 48, by = 2), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(from = 26, to = 42, by = 2), expand = c(0,0)) +
    geom_sf(data = countries, colour = "black", fill="#D7D7D7")+ 
    scale_fill_gradientn(colours = pal, 
                         #breaks = c(5, 50, 95, 140),
                         name="FRic (SES)") +
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
map2

hist(SES_FRic2)

