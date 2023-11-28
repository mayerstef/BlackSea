# Load packages (leaflet allows to load google maps)
library(sdmpredictors)
library(leaflet)
library(dplyr)


xmin <- 26.223404293
xmax <- 41.923404293
ymin <- 39.879949463
ymax <- 47.379949463


# Generate a data.frame with the sites of interest
my.sites <- data.frame(Name="Black Sea", Lon=c(xmin, xmax) , Lat=c(ymin, ymax) )
my.sites

# Visualise sites of interest in google maps
m <- leaflet()
m <- addTiles(m)
m <- addMarkers(m, lng=my.sites$Lon, lat=my.sites$Lat, popup=my.sites$Name)
m


# List layers avaialble in Bio-ORACLE v2
layers.bio2 <- list_layers( datasets="Bio-ORACLE" )

# Download environmental data layers (Max. Temperature, Min. Salinity and Min. Nitrates at the sea bottom)
env.layers <- load_layers( layer_codes <- c("BO22_chlorange_bdmean",
                                                    "BO22_chlorange_ss",
                                                    "BO22_chlomean_bdmean",
                                                    "BO22_chlomean_ss",
                                                    "BO22_dissoxmean_bdmean",
                                                    "BO22_dissoxmean_ss",
                                                    "BO22_dissoxrange_bdmean",
                                                    "BO22_dissoxrange_ss",
                                                    "BO22_ironmean_bdmean",
                                                    "BO22_ironmean_ss",
                                                    "BO22_ironrange_bdmean",
                                                    "BO22_ironrange_ss",
                                                    "BO22_nitratemean_bdmean",
                                                    "BO22_nitratemean_ss",
                                                    "BO22_nitraterange_bdmean",
                                                    "BO22_nitraterange_ss",
                                                    "BO22_phosphatemean_bdmean",
                                                    "BO22_phosphatemean_ss",
                                                    "BO22_phosphaterange_bdmean",
                                                    "BO22_phosphaterange_ss",
                                                    "BO22_ppmean_bdmean",
                                                    "BO22_ppmean_ss",
                                                    "BO22_pprange_bdmean",
                                                    "BO22_pprange_ss",
                                                    "BO22_salinitymean_bdmean",
                                                    "BO22_salinitymean_ss",
                                                    "BO22_salinityrange_bdmean",
                                                    "BO22_salinityrange_ss",
                                                    "BO22_silicatemean_bdmean",
                                                    "BO22_silicatemean_ss",
                                                    "BO22_silicaterange_bdmean",
                                                    "BO22_silicaterange_ss",
                                                    "BO22_tempmean_bdmean",
                                                    "BO22_tempmean_ss",
                                                    "BO22_temprange_bdmean",
                                                    "BO22_temprange_ss",
                                                    "BO22_ph"),
                                   equalarea=FALSE, rasterstack=TRUE)


cropped_stack <- crop(env.layers, extent(xmin, xmax, ymin, ymax))
cropped_stack <- projectRaster(cropped_stack, crs = proj84)
env_df <- as.data.frame(cropped_stack, xy = TRUE)
env_df$x <- round(env_df$x, 1)
env_df$y <- round(env_df$y, 1)
env_df$grids <- paste(env_df$x, env_df$y, sep = "_")
colnames(env_df)[3:39] <- gsub("^BO22_", "", colnames(env_df)[3:39])
env_df <- env_df %>%
  rename_with(toupper, c(1, 2)) %>% 
  relocate(grids, .before = 1) %>%
  na.omit()
#save(cropped_stack, file="Env_BiOracle.RData")
#save(env_df, file="env_df.RData")
env_df<-inner_join(env_df, grids)
env_df<-aggregate(.~grids, data=env_df, FUN=mean)
# Get the common columns that should be in all dataframes
common_cols <- c("grids", "X", "Y", "Bioregions")
VIF_analysis <- function(x){
  x <- as.data.frame(x)
  
  varname <- vector()
  Rsquared <- vector()
  VIF <- vector()
  
  for(i in 1:ncol(x)){
    varname <- c(varname, colnames(x)[i])
    mod <- lm(data=x[,-i], x[,i]~.)
    
    R2 <- summary(mod)$r.squared
    Rsquared <- c(Rsquared,R2)
    
    VIF <- c(VIF,1/(1-R2))
  }
  output <- data.frame(variable=varname, Rsquared=Rsquared, VIF=VIF)
}

selected_vars <- c("dissoxrange_ss",
                   "dissoxrange_bdmean", 
                   "chlorange_bdmean", 
                   "ironrange_ss",
                   "nitraterange_ss", 
                   "salinityrange_ss", 
                   "silicaterange_ss", 
                   "temprange_ss", 
                   "phosphatemean_bdmean",
                   "ppmean_bdmean", #this one was the last one added and can be removed 
                   "ph")

selected_vars2 <- c("dissoxrange_ss",
                   "dissoxrange_bdmean", 
                   "chlorange_bdmean", 
                   "ironrange_ss",
                   "nitraterange_ss", 
                   "salinityrange_ss", 
                   "silicaterange_ss", 
                   "temprange_ss", 
                   "phosphatemean_bdmean",
                   "ppmean_bdmean",
                   "tempmean_bdmean",
                   "ph")

env <- env_df[, selected_vars]
vif_results<-VIF_analysis(env)

library(vegan)
env_stand <- scale(env) 
env.pca <- rda(env_stand)
par(mfrow = c(1,1))
screeplot(env.pca, bstick=TRUE)
summary(env.pca) #Cumulative Proportion 0.348 0.5487
par(mfrow = c(1,2))
cleanplot.pca(env.pca, scaling = 1, plot.sites = F, label.sites = F)
cleanplot.pca(env.pca, scaling = 2, plot.sites = F, label.sites = F)
# Create a dataframe for each suffix
env_df<-tibble::column_to_rownames(env_df, var="grids")
#range_ss
e1<- env_df[, c(grep("range_ss$", names(env_df), value = TRUE))]
#range_bdmean 
e2 <- env_df[, c(grep("range_bdmean$", names(env_df), value = TRUE))]
#mean_ss 
e3<- env_df[, c(grep("mean_ss$", names(env_df), value = TRUE))]
#mean_bdmean
e4<- env_df[, c(grep("mean_bdmean$", names(env_df), value = TRUE))]
















