library(phyloregion)
library(sf)
library(sp)
library(dplyr)

createSPDF <- function(region_name, data, spdf){
  spd <- sf::as_Spatial(st_geometry(data[data$UNION == region_name, ]), IDs = as.character(1:nrow(data[data$UNION == region_name, ])))
  df <- data[data$UNION == region_name, ]
  df$geometry <- NULL
  df <- as.data.frame(df)
  row.names(df) <- sapply(slot(spd, "polygons"), function(x) slot(x, "ID")) # Ensure row names of df match IDs of polygons
  spd <- sp::SpatialPolygonsDataFrame(spd, data = df)
  proj4string(spd) <- proj84
  spdf_region <- selectbylocation(spdf, spd)
  return(spdf_region)
}
# And then you can call this function for each of your regions:
spdf_tmar <- createSPDF("Turkey_marmara", Bioreg_bf, spdf)
spdf_ua <- createSPDF("Ukraine_azov", Bioreg_bf, spdf)
spdf_ra <- createSPDF("Russia_azov", Bioreg_bf, spdf)
spdf_G <- createSPDF("Georgia", Bioreg_bf, spdf)
spdf_B <- createSPDF("Bulgaria", Bioreg_bf, spdf)
spdf_Rom <- createSPDF("Romania", Bioreg_bf, spdf)
spdf_U <- createSPDF("Ukraine", Bioreg_bf, spdf)
spdf_Rus <- createSPDF("Russia", Bioreg_bf, spdf)
spdf_T_bs1 <- createSPDF("Turkey_bs1", Bioreg_bf, spdf)
spdf_T_bs2 <- createSPDF("Turkey_bs2", Bioreg_bf, spdf)

B_rows<-rownames(spdf_B@data) %>% as.data.frame()
B_rows$Region<-"Bulgaria"
Tmar_rows <- rownames(spdf_tmar@data) %>% as.data.frame()
Tmar_rows$Region <- "Turkey_marmara"
Ua_rows <- rownames(spdf_ua@data) %>% as.data.frame()
Ua_rows$Region <- "Ukraine_azov"
Ra_rows <- rownames(spdf_ra@data) %>% as.data.frame()
Ra_rows$Region <- "Russia_azov"
G_rows <- rownames(spdf_G@data) %>% as.data.frame()
G_rows$Region <- "Georgia"
Rom_rows <- rownames(spdf_Rom@data) %>% as.data.frame()
Rom_rows$Region <- "Romania"
U_rows <- rownames(spdf_U@data) %>% as.data.frame()
U_rows$Region <- "Ukraine"
Rus_rows <- rownames(spdf_Rus@data) %>% as.data.frame()
Rus_rows$Region <- "Russia"
T_bs1_rows <- rownames(spdf_T_bs1@data) %>% as.data.frame()
T_bs1_rows$Region <- "Turkey_bs1"
T_bs2_rows <- rownames(spdf_T_bs2@data) %>% as.data.frame()
T_bs2_rows$Region <- "Turkey_bs2"

Bioreg_rows <- rbind(B_rows, Tmar_rows, Ua_rows, Ra_rows, G_rows, Rom_rows, 
                     U_rows, Rus_rows, T_bs1_rows, T_bs2_rows)
colnames(Bioreg_rows)[1] <- "grids"
Bioreg_rows <- Bioreg_rows[!duplicated(Bioreg_rows$grids), ]

BSRegion1<-long2sparse(Bioreg_rows, grids = "grids", species = "Region") %>%
  as.matrix()%>%
  as.data.frame() %>% 
  tibble::rownames_to_column("grids")


# And then you can call this function for each of your regions:
spdf_tmar <- createSPDF("Turkey_marmara", Bioreg, spdf)
spdf_ua <- createSPDF("Ukraine_azov", Bioreg, spdf)
spdf_ra <- createSPDF("Russia_azov", Bioreg, spdf)
spdf_G <- createSPDF("Georgia", Bioreg, spdf)
spdf_B <- createSPDF("Bulgaria", Bioreg, spdf)
spdf_Rom <- createSPDF("Romania", Bioreg, spdf)
spdf_U <- createSPDF("Ukraine", Bioreg, spdf)
spdf_Rus <- createSPDF("Russia", Bioreg, spdf)
spdf_T_bs1 <- createSPDF("Turkey_bs1", Bioreg, spdf)
spdf_T_bs2 <- createSPDF("Turkey_bs2", Bioreg, spdf)

#plot(spdf)
#plot(spdf_tmar, col = "red",add=T)
#plot(spdf_ua, col = "blue", add = TRUE)
#plot(spdf_ra, col = "green", add = TRUE)
#plot(spdf_G, col = "yellow", add = TRUE)
#plot(spdf_B, col = "orange", add = TRUE)
#plot(spdf_Rom, col = "purple", add = TRUE)
#plot(spdf_U, col = "cyan", add = TRUE)
#plot(spdf_Rus, col = "brown", add = TRUE)
#plot(spdf_T_bs1, col = "pink", add = TRUE)
#plot(spdf_T_bs2, col = "gray", add = TRUE)


B_rows<-rownames(spdf_B@data) %>% as.data.frame()
B_rows$Region<-"Bulgaria"
Tmar_rows <- rownames(spdf_tmar@data) %>% as.data.frame()
Tmar_rows$Region <- "Turkey_marmara"
Ua_rows <- rownames(spdf_ua@data) %>% as.data.frame()
Ua_rows$Region <- "Ukraine_azov"
Ra_rows <- rownames(spdf_ra@data) %>% as.data.frame()
Ra_rows$Region <- "Russia_azov"
G_rows <- rownames(spdf_G@data) %>% as.data.frame()
G_rows$Region <- "Georgia"
Rom_rows <- rownames(spdf_Rom@data) %>% as.data.frame()
Rom_rows$Region <- "Romania"
U_rows <- rownames(spdf_U@data) %>% as.data.frame()
U_rows$Region <- "Ukraine"
Rus_rows <- rownames(spdf_Rus@data) %>% as.data.frame()
Rus_rows$Region <- "Russia"
T_bs1_rows <- rownames(spdf_T_bs1@data) %>% as.data.frame()
T_bs1_rows$Region <- "Turkey_bs1"
T_bs2_rows <- rownames(spdf_T_bs2@data) %>% as.data.frame()
T_bs2_rows$Region <- "Turkey_bs2"

Bioreg_rows2 <- rbind(B_rows, Tmar_rows, Ua_rows, Ra_rows, G_rows, Rom_rows, 
                     U_rows, Rus_rows, T_bs1_rows, T_bs2_rows)
colnames(Bioreg_rows2)[1] <- "grids"
#duplicates <- Bioreg_rows$grids[duplicated(Bioreg_rows$grids)]

rm(B_rows, Tmar_rows, Ua_rows, Ra_rows, G_rows, Rom_rows, 
    U_rows, Rus_rows, T_bs1_rows, T_bs2_rows)
rm(spdf_tmar, spdf_ua, spdf_ra, spdf_G, spdf_B, spdf_Rom, spdf_U, spdf_Rus, spdf_T_bs1, spdf_T_bs2)

BSRegion2<-long2sparse(Bioreg_rows2, grids = "grids", species = "Region") %>%
  as.matrix()%>%
  as.data.frame() %>% 
  tibble::rownames_to_column("grids")

BSRegion3<-anti_join(BSRegion1, BSRegion2, by= "grids")
BSRegions<-rbind(BSRegion2, BSRegion3)
rm(Bioreg_rows, Bioreg_rows2, BSRegion1, BSRegion2, BSRegion3)
