################################################################################
####                                                                        ####
####         Master thesis:LINKING PATTERNS IN PHYLOGENY, TRAITS,           ####
####            ABIOTIC VARIABlES AND SPACE: BLACK SEA FISHES               ####
####                                                                        ####
####   Stefanie Mayer                                                       ####
####                                                                        ####
####        calbouy@ethz.ch                                                 ####
################################################################################

#================== Workspace preparation & Data loading ======================#
### Library loading
# Data Manipulation and Analysis
lib_general <- c("dplyr", "reshape2","readxl", "ggplot2", "Matrix", "stringr")

# Spatial Data and Geospatial Analysis
lib_gis <- c("raster", "rgeos", "rgdal", "sf", "sp", "shape", "geometry","geosphere", 
             "lwgeom", "fasterize", "terra","rnaturalearthdata", "rnaturalearth", 
             "shape", "ggspatial", "ggpubr", "prettymapr")

# Phylogenetic Analysis
lib_phylo <- c("ape", "PhyloMeasures", "tidytree", "TreeTools", "phyloregion")

#Choose which libraries you would like to load
sapply(lib_general,library,character.only=TRUE)
sapply(lib_gis,library,character.only=TRUE)
sapply(lib_phylo,library,character.only=TRUE)

### Set Working directory
setwd("MAYER_MasterThesis")
setwd("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis")

### Projection definition
proj84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") 
sf::sf_use_s2(F) #global option for using the S2 geometry library to FALSE

### Load Data
world <- ne_countries(scale = "medium", returnclass = "sf")
coast <- st_read("./Data/GIS/coastline_blacksea.shp") #only used for maps, not so important
bathy <- raster("./Data/GIS/gebco_blacksea.asc") 
Bioreg<-st_read("./Data/GIS/10bioregions_blacksea.shp") #10 features, see script "BlackSea_Bioregions" for shapefile creation
Bioreg <- st_transform(Bioreg, proj84)
Bioreg_bf <- st_buffer(Bioreg, dist = .1)

#Load finalized traits extracted from Fishbase (see "FishBase Extraction" script)
traits_sp <- read_excel("./Data/traits_sp.xlsx", 
                        col_types = c("numeric", "text", "text", 
                                      "text", "text", "text", "text", "text", 
                                      "text", "numeric", "numeric", "text", 
                                      "text", "text", "text", "text", "text", 
                                      "numeric", "numeric", "text", "text", 
                                      "text", "text", "text", "text"))
traits_sp$X <- gsub(" ", "_", traits_sp$X, fixed=TRUE)


#==================== PART 1 - Manual Distribution Maps =======================#
source("Scripts/Sources/manual_distribution_rasters.R")

# Create list of CSV file names in Data/Data_occurences/ directory
listfile <- list.files("Data/Occurences/input/", pattern = "\\.csv$")
listfile <- gsub("\\.csv$", "", listfile)

for (rep_sp in listfile[1:length(listfile)]) {
  print(rep_sp)
  # Read occurrence data
  Occ <- read.csv(paste0("Data/Occurences/input/", rep_sp, ".csv"),
                  stringsAsFactors = FALSE, header = TRUE)
  
  # Get species name and bathymetric range
  Species_name <- rep_sp
  ind <- which(traits_sp[,"X"] == Species_name)
  
  # Extract bathymetric range from traits_sp data frame
  bathy_sp <- c(traits_sp[ind, "Depth_min"], traits_sp[ind, "Depth_max"])
  
  # Correct bathymetric range: 
  # set min depth to -10 if environment is freshwater ( for Anadromous/Catadromous species)
  # check if Marine, then pelagic, then see if it`s
  # benthopelagic 
  # pelagic-neritic: restricted to continental shelf (200m)
  # oceanic-pelagic: no min depth 
  if (traits_sp$Env_4[ind] == "Freshwater") { 
    bathy_sp[1] <- -10
  }
  if (str_detect(traits_sp$Env_2[ind], "pelagic") && str_detect(traits_sp$Env_1[ind], "Marine")) {
    if (traits_sp$Env_3[ind] != "benthopelagic" && traits_sp$Env_3[ind] != "pelagic-neritic") {
      bathy_sp[2] <- 11000  
    }
    if (traits_sp$Env_3[ind] == "pelagic-neritic" && bathy_sp[2] < 200) {
      bathy_sp[2] <- 200  
    }
    if (traits_sp$Env_3[ind] == "oceanic-pelagic" && bathy_sp[2] < 5000) {
      bathy_sp[2] <- 5000
    }
  }  
  
  # Select necessary columns and group by longitude, latitude, and species
  Occ <- Occ[, c("decimalLongitude", "decimalLatitude", "species")]
  Occ <- Occ %>%  
    dplyr::rename(longitude = decimalLongitude, latitude = decimalLatitude) %>% 
    dplyr::group_by(longitude, latitude, species) %>%
    dplyr::summarise(n_sp = n_distinct(species)) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-n_sp) %>% 
    dplyr::filter(!(longitude == 35.0 & latitude == 43.0)) #remove outlier (from museum in Occurence data)
  
  print(dim(Occ))
  
  # Filter occurrences based on distance and/or random sampling
  # Filter out locations with more than 80,000 occurrences
  if (dim(Occ)[1] > 80000) {
    Occ <- Occ %>%
      dplyr::select(longitude, latitude, species)
    # Compute minimum distance between occurrence points and filter out locations beyond the 90th percentile
  } else if (dim(Occ)[1] >= 20) {
    matrice_dist <- spDists(coordinates(Occ[, c("longitude", "latitude")]))
    matrice_dist[matrice_dist == 0] <- NA
    dist_min_3 <- apply(matrice_dist, 2, find.5.min)
    Occ$dist_min <- apply(dist_min_3, 2, mean)
    limit_distance <- quantile(dist_min_3, 0.9)
    Occ <- Occ %>%
      dplyr::filter(dist_min <= limit_distance) 
    print(dim(Occ))
    # Keep locations with less than 20 occurrences
  } else {
    Occ <- Occ %>%
      dplyr::select(longitude, latitude, species)
  }
  
  # Sample at most 130,000 occurrence points
  # Randomly sample occurrences if too many
  if (nrow(Occ) > 130000) {
    Occ <- Occ[sample(1:nrow(Occ), 130000),]
  }
  
  # Apply functions to occurrence data
  Rast_sp <- get_gbif2rast_data(Occ = Occ %>% 
                                  select(longitude, latitude, species),
                                proj = proj84,
                                Bioreg = Bioreg,
                                Bioreg_bf = Bioreg_bf,
                                Species_name = Species_name,
                                name_shp = "UNION",
                                bathy_sp = bathy_sp,
                                corrected_bathy = TRUE,
                                bathy = bathy)
  
  if( !is.character(Rast_sp)){
    writeRaster(Rast_sp, filename= paste0("Rasters/manual/", Species_name,".asc"), overwrite= T)
    write.csv(Occ, file = paste0("Data/Occurences/output", Species_name,".csv"))
    
  }
  
} #end 

#===================== PART 2 - IUCN Distribution Maps ========================#

#Load IUCN polygons which were downloaded from IUCN Red List website 
# see "IUCN process" script for details on how they were edited
IUCN<-st_read("GIS/IUCN_polygons.shp")

#Rename species to match traits_sp
IUCN <- IUCN %>% 
  dplyr::mutate(X = case_when(
    X == "Mycteroperca caninus" ~ "Epinephelus_caninus",
    X == "Mycteroperca costae" ~ "Epinephelus_costae",
    X == "Mycteroperca marginatus" ~ "Epinephelus_marginatus",
    X == "Chelon ramada" ~ "Liza_ramada",
    X == "Chelon saliens" ~ "Liza_saliens",
    X == "Mullus barbatus" ~ "Mullus_barbatus_barbatus",
    X == "Auxis rochei" ~ "Auxis_rochei_rochei",
    X == "Gobius ophiocephalus" ~ "Zosterisessor_ophiocephalus",
    TRUE ~ X
  ))
IUCN$X<-gsub(" ", "_", IUCN$X, fixed=TRUE)

select<-intersect(IUCN$X, traits_sp$X) 
(remove<-setdiff(traits_sp$X, IUCN$X))
subset_IUCN <- IUCN[IUCN$X %in% select, ]

# Select Species that only have freshwater distribution in IUCN and buffer them
selected_species <- c(
  "Alosa fallax",
  "Alosa immaculata",
  "Alosa maeotica",
  "Alosa tanaica",
  "Aphanius fasciatus", 
  "Babka gymnotrachelus",
  "Caspiosoma caspium",
  "Clupeonella cultriventris",
  "Cyprinus carpio",
  "Gymnocephalus cernua",
  "Liza ramada",
  "Neogobius fluviatilis",
  "Pelecus cultratus",
  "Ponticola kessleri",
  "Ponticola syrman",
  "Proterorhinus semilunaris",
  "Pungitius platygaster",
  "Salmo labrax",
  "Salmo trutta",
  "Squalius cephalus"
) 
selected_species<-gsub(" ", "_", selected_species)

species_to_buffer <- subset_IUCN[subset_IUCN$X %in% selected_species, ]
rest<-subset_IUCN[!subset_IUCN$X %in% selected_species, ]
buffered <- st_buffer(species_to_buffer, dist = .1) 
final_IUCN<-rbind(rest, buffered)
sf_species <- st_as_sf(final_IUCN) # convert shp_species to an sf object
sf_species$sci_name <-sf_species$X

IUCNcrop <- function(sf_species, traits_sp, bathy) {
  # create empty list to store output rasters
  output_rasters <- list()
  # loop through each species in the sf_species data frame
  for (i in seq_len(nrow(sf_species))) {
    print(sf_species$sci_name[i])
    # Get species name and bathymetric range from traits_sp data frame
    Species_name <- sf_species$sci_name[i]
    ind <- which(traits_sp$X == Species_name)
    # Check if ind is empty
    if (length(ind) == 0) {
      warning(paste0("No match found for:", Species_name, " :/ "))
      next
    }
    # Extract bathymetric range from traits_sp data frame
    bathy_sp <- c(traits_sp[ind, "Depth_min"], traits_sp[ind, "Depth_max"])
    # Correct bathymetric range based on environment
    if (traits_sp$Env_4[ind] == "Freshwater") {
      bathy_sp[1] <- -10
    }
    if (str_detect(traits_sp$Env_2[ind], "pelagic") && str_detect(traits_sp$Env_1[ind], "Marine")) {
      if (traits_sp$Env_3[ind] != "benthopelagic" && traits_sp$Env_3[ind] != "pelagic-neritic") {
        bathy_sp[2] <- 11000  
      }
      if (traits_sp$Env_3[ind] == "pelagic-neritic" && bathy_sp[2] < 200) {
        bathy_sp[2] <- 200  
      }
      if (traits_sp$Env_3[ind] == "oceanic-pelagic" && bathy_sp[2] < 5000) {
        bathy_sp[2] <- 5000
      }
    }  
    
    # Crop the bathymetry raster to the species polygon and rasterize polygons
    row_data <- sf_species[i,]
    r1 <- crop(bathy, row_data) 
    raster_sp <- fasterize(row_data, r1)
    # Extract raster values from the bathymetry raster and correct for species range
    data_sp <- values(raster_sp)
    data_r1 <- values(r1)
    data_r1[which(is.na(data_sp))] <- NA 
    data_r1[which(data_r1 > -bathy_sp$Depth_min)] <- NA
    data_r1[which(data_r1 <= -bathy_sp$Depth_max)] <- NA
    values(raster_sp) <- data_r1 # replace species raster values with corrected bathymetry values
    values(raster_sp) <- ifelse(is.na(data_r1), 0, 1)
    # Save output raster to list
    output_rasters[[i]] <- raster_sp
    print(paste0("Output raster saved for species ", sf_species$sci_name[i]))
  }
  # Write the output rasters to file
  filenames <- paste0("Rasters/IUCN/", sf_species$sci_name, ".asc")
  for (i in seq_along(output_rasters)) {
    writeRaster(output_rasters[[i]], filename = filenames[i], overwrite = TRUE)
  }
  return(output_rasters)
}

# Call IUCNcrop function
IUCN_rasters <- IUCNcrop(sf_species, traits_sp, bathy)


#===================== PART 3 - FAO Distribution Maps =========================#
setwd("Data/FAO")
zip_files <- list.files(pattern = "*.zip")
# loop through each file in the zip_files vector and extract its contents
for (file in zip_files) {
  unzip(file, exdir = ".")
}

FAOshp <- list.files(pattern = ".shp", recursive = TRUE)
FAOall<- do.call(rbind, lapply(FAOshp, st_read))
FAOall <- left_join(FAOall, species %>% select(code, scientific_name), by = c("ALPHACODE" = "code")) %>% 
  st_transform(., proj84)

setwd("MAYER_MasterThesis")
names(FAOall)[names(FAOall) == "scientific_name"] <- "X"
sf_species <- st_as_sf(FAOall) # convert shp_species to an sf object
sf_species$sci_name <- gsub(" ", "_", sf_species$X)

FAOcrop <- function(sf_species, traits_sp, bathy) {
  # create empty list to store output rasters
  output_rasters <- list()
  # loop through each species in the sf_species data frame
  for (i in seq_len(nrow(sf_species))) {
    print(sf_species$sci_name[i])
    # Get species name and bathymetric range from traits_sp data frame
    Species_name <- sf_species$sci_name[i]
    ind <- which(traits_sp$X == Species_name)
    # Check if ind is empty
    if (length(ind) == 0) {
      warning(paste0("No match found for:", Species_name, " :/ "))
      next
    }
    # Extract bathymetric range from traits_sp data frame
    bathy_sp <- c(traits_sp[ind, "Depth_min"], traits_sp[ind, "Depth_max"])
    # Correct bathymetric range based on environment
    if (traits_sp$Env_4[ind] == "Freshwater") {
      bathy_sp[1] <- -10
    }
    if (str_detect(traits_sp$Env_2[ind], "pelagic") && str_detect(traits_sp$Env_1[ind], "Marine")) {
      if (traits_sp$Env_3[ind] != "benthopelagic" && traits_sp$Env_3[ind] != "pelagic-neritic") {
        bathy_sp[2] <- 11000  
      }
      if (traits_sp$Env_3[ind] == "pelagic-neritic" && bathy_sp[2] < 200) {
        bathy_sp[2] <- 200  
      }
      if (traits_sp$Env_3[ind] == "oceanic-pelagic" && bathy_sp[2] < 5000) {
        bathy_sp[2] <- 5000
      }
    }  
    
    # Crop the bathymetry raster to the species polygon and rasterize polygons
    row_data <- sf_species[i,]
    r1 <- crop(bathy, row_data) 
    raster_sp <- fasterize(row_data, r1)
    # Extract raster values from the bathymetry raster and correct for species range
    data_sp <- values(raster_sp)
    data_r1 <- values(r1)
    data_r1[which(is.na(data_sp))] <- NA 
    data_r1[which(data_r1 > -bathy_sp$Depth_min)] <- NA
    data_r1[which(data_r1 <= -bathy_sp$Depth_max)] <- NA
    values(raster_sp) <- data_r1 # replace species raster values with corrected bathymetry values
    values(raster_sp) <- ifelse(is.na(data_r1), 0, 1)
    # Save output raster to list
    output_rasters[[i]] <- raster_sp
    print(paste0("Output raster saved for species ", sf_species$sci_name[i]))
  }
  # Write the output rasters to file
  filenames <- paste0("Rasters/FAO", sf_species$sci_name, ".asc")
  for (i in seq_along(output_rasters)) {
    if (!is.null(output_rasters[[i]])) {
      writeRaster(output_rasters[[i]], filename = filenames[i], overwrite = TRUE)
    } else {
      warning(paste0("Skipping saving NULL raster for species ", sf_species$sci_name[i]))
    }
  }
}

# Call IUCNcrop function
FAO_rasters <- FAOcrop(sf_species, traits_sp, bathy)


#==================== PART 4 - Fishery Distribution Maps ======================#

Bioregion_fished <- st_read("Data/GIS/fishedEEZ_blacksea.shp")
fisheries<-read_xlsx("Data/fisheries.xlsx")
unique_species <- unique(fisheries$Species)
num_species <- length(unique_species)
species_depth <- fisheries %>%
  group_by(Species) %>%
  summarise(Depth_min = unique(Depth_min),
            Depth_max = unique(Depth_max))


process_species <- function(species_name) {
  ind <- which(species_depth$Species == species_name)
  bathy_sp <- c(species_depth[ind, "Depth_min"], species_depth[ind, "Depth_max"])
  
  if (fisheries$Env_4[ind] == "Freshwater") {
    bathy_sp[1] <- -10
    print(paste("Setting bathy_sp[1] to -10. New bathy_sp:", bathy_sp))
    flush.console()
  }
  if (str_detect(fisheries$Env_2[ind], "pelagic") && str_detect(fisheries$Env_1[ind], "Marine")) {
    if (fisheries$Env_3[ind] != "benthopelagic" && fisheries$Env_3[ind] != "pelagic-neritic") {
      bathy_sp[2] <- 11000
      print(paste("Setting bathy_sp[2] to 11000. New bathy_sp:", bathy_sp))
      flush.console()
    }
    if (fisheries$Env_3[ind] == "pelagic-neritic" && bathy_sp[2] < 200) {
      bathy_sp[2] <- 200
      print(paste("Setting bathy_sp[2] to 200. New bathy_sp:", bathy_sp))
      flush.console()
    }
    if (fisheries$Env_3[ind] == "oceanic-pelagic") {
      if (is.na(bathy_sp[2])) {
        bathy_sp[2] <- Inf
        print(paste("Setting bathy_sp[2] to Inf. New bathy_sp:", bathy_sp))
        flush.console()
      }
    }
    # add condition to set min depth to -10 if environment is freshwater
  }
  species_eez <- fisheries  %>% filter(Species == species_name ) %>% pull(EEZ)
  if ("USSR" %in% species_eez) {
    species_eez <- species_eez[species_eez != "USSR"]
    species_eez <- c(species_eez, "Russia", "Ukraine", "Russia_azov", "Ukraine_azov")
  }
  if ("USSR_azov" %in% species_eez) {
    species_eez <- species_eez[species_eez != "USSR_azov"]
    species_eez <- c(species_eez, "Russia_azov", "Ukraine_azov")
  }
  eez_subset <- Bioregion_fished  %>% filter(UNION %in% species_eez)
  r1 <- crop(bathy, eez_subset)
  raster_sp <- fasterize(eez_subset, r1)
  
  data_sp <- values(raster_sp)
  data_r1 <- values(r1)
  data_r1[which(is.na(data_sp))] <- NA
  data_r1[which(data_r1 > -bathy_sp$Depth_min)] <- NA
  data_r1[which(data_r1<= -bathy_sp$Depth_max)] <- NA
  values(raster_sp) <- data_r1
  
  val <- values(raster_sp)
  val[!is.na(val)] <- 1
  val[is.na(val)] <- 0
  values(raster_sp) <- val
  
  writeRaster(raster_sp, filename=paste0("Rasters/fishery/", species_name, ".asc"), overwrite=TRUE, format="ascii")
  print(paste0("Output raster saved for species ", species_name))
}

# Process each unique species
for (i in seq_len(num_species)) {
  print(unique_species[i])
  process_species(unique_species[i])
}

######################################################################################
decisiontree <- read.csv("Data/decisiontree.csv")
decisiontree$X<-gsub(" ", "_", decisiontree$X)
ExcludedSpecies <- decisiontree[rowSums(decisiontree[,2:5] == 0) == 4,]
SpeciesToCopy <- decisiontree[rowSums(decisiontree[,2:5] == 1) == 1,]
decisiontree <- decisiontree %>% anti_join(ExcludedSpecies) 
coast_bf<-Bioreg_bf  %>%
  summarise()
coast_bf<-st_transform(coast_bf, proj84)

xmin <- 26.223404293 #lon
xmax <- 41.923404293 #lon
ymin <- 39.879949463 #lat
ymax <- 47.379949463 #lat


PAMashup <- function(df, spatialmask, xmin, xmax, ymin, ymax) {
  # Create the final folder if it doesn't exist
  dir.create("Rasters/final", showWarnings = FALSE)
  new_extent <- raster::extent(xmin, xmax, ymin, ymax)
  
  process_species <- function(i) {
    species_name <- df$X[i]
    folder_flags <- df[i, 2:5] # Extract binary values from columns 2 to 5
    
    # Check if all folder flags are 0; if so, skip this species
    if (sum(folder_flags) == 0) {
      print(paste0("Skipped species: ", species_name))
      return(NULL)
    }
    
    # Get the folder names where the binary value is 1
    active_folders <- colnames(df)[2:5][folder_flags == 1]
    print(paste0("Processing species: ", species_name))
    
    # Process raster files based on the active folders
    rasters_to_process <- lapply(active_folders, function(folder) {
      raster_to_crop <- raster(file.path(folder, paste0(species_name, ".asc")))
      raster_to_mask <- crop(raster_to_crop, new_extent)
      masked_raster <- mask(raster_to_mask, spatialmask, updatevalue=0)
      
      return(masked_raster) # Add the missing return statement
    })
    
    if (length(rasters_to_process) == 1) {
      processed_raster <- rasters_to_process[[1]]
      print(paste0("Copied file for ", species_name))
    } else {
      processed_raster <- do.call(raster::mosaic, c(rasters_to_process, fun = max))
      print(paste0("Merged rasters and saved file for ", species_name))
    }
    
    writeRaster(processed_raster, file.path("Rasters/final/", paste0(species_name, ".asc")), overwrite = TRUE)
    
  }
  
  lapply(seq_len(nrow(df)), process_species)
}

PAMashup(decisiontree,coast_bf, xmin, xmax, ymin, ymax)

########################################################################################
#this one is proper masking
listfile <- list.files("Rasters/final/")

rgb <- col2rgb("grey70")
col <- rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255, alpha = 250)

lapply(listfile, function(i){
  if (!dir.exists("Rastes/Maps")) {
    dir.create("Rastes/Maps")
  }
  Rast_sp <- raster(paste0("Rasters/final/", i))
  Species_name <- str_sub(i, start = 1, end = -5)
  # Set the map coordinates
  Min_x <- 25.4073
  Min_y <- 38.2849
  Max_x <- 43.6943
  Max_y <- 48.000
  # Load world countries
  countries <- ne_countries(scale = "medium", returnclass = "sf")
  countries <- st_as_sf(countries, crs = proj84)
  # Plot the map
  png(paste0("maps/", Species_name, ".png"), width = 2000, height = 1400, units = "px", pointsize = 10, res = 300)
  plot(rnorm(1000), col = "white", type = "n", axes = FALSE, 
       xlab = "", ylab = "", 
       xlim = c(Min_x, Max_x), ylim = c(Min_y, Max_y), main = as.expression(bquote(paste("Species Distribution of ", italic(.(Species_name))))), 
       mar = c(5, 5, 4, 2) + 0.1)
  image(Rast_sp, col = c("transparent", "royalblue1"), add = TRUE)
  plot(coast, col = "black", add = TRUE, lwd = 0.5)
  plot(countries, col = "gray80", add = TRUE)
  #Add the x and y axis ticks
  axis(side = 1, line = 0, cex.axis = 0.6, lwd = 0.5, tcl = -0.25, bg = "white", labels = c("25°E", "29°E", "33°E", "38°E", "44°E"), at = c(25, 29, 33, 38, 44), mgp = c(3, 0.25, 0))
  axis(side = 2, line = 0, las = 2, cex.axis = 0.6, lwd = 0.5, tcl = -0.25, bg = "white", labels = c("38°N", "42°N", "46°N", "50°N"), at = c(38, 42, 46, 50), mgp = c(3, 0.5, 0))
  addnortharrow(
    pos = "topright",
    padin = c(0.15, 0.15),
    scale = 0.5,
    lwd = 1,
    border = "black",
    cols = c("white", "black"),
    text.col = "black"
  )
  addscalebar(
    plotunit = NULL,
    plotepsg = NULL,
    widthhint = 0.15,
    unitcategory = "metric",
    htin = 0.1,
    padin = c(0.15, 0.15),
    style = "bar",
    bar.cols = c("black", "white"),
    lwd = 1,
    linecol = "black",
    tick.cex = 0.7,
    labelpadin = 0.08,
    label.cex = 0.8,
    label.col = "black",
    pos = "bottomright"
  )
  dev.off()
})

##################################################################################

# Intersection function 2 load
get_vect_species <- function(grd=grd, sp=Sp1){
  Sp1 <- raster(sp)
  pts_d <- coordinates(Sp1)[which(values(Sp1) == 1), ]
  rast_cel <- unique(cellFromXY(grd, pts_d))
  a <- grd
  a[] <- 0
  a[][rast_cel] <- 1
  a[]
} 

list_species <- list.files("Rasters/final/", pattern = "\\.asc$")

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

# Transform your grid and extract the grid coordinates
grd_coord <- as(SpatialGrid(grd), "SpatialPolygons")
coord <- do.call(rbind, lapply(1:length(grd_coord@polygons), function(x) {
  grd_coord@polygons[[x]]@labpt
}))
#grille <- data.frame(ID = seq(1, dim(coord)[1], 1), X = coord[,1], Y = coord[,2])

#Decimal place to round coordinates (depending on resolution chosen)
dp <- 1 #1 for 0.1 or 0.5 resolution 

grille <- data.frame(ID = paste(round(coord[,1], dp), round(coord[,2], dp), sep = "_"),
                     X = round(coord[,1], dp),
                     Y = round(coord[,2], dp))

# Transform your grid in raster
grd <- raster(grd)
grd[] <- 0
grd@crs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Intersection between the grid and each species raster
res <- lapply(1:length(list_species), function(x){
  cat("Processing: ", list_species[x], "\n")
  get_vect_species(sp = list_species[x], grd = grd)
})

list_species <- gsub("\\.asc$", "", list_species)

# Matrix P/A  
mat_PA <- cbind(grille, do.call(cbind, res)) 
mat_PA <- mat_PA%>%
  tibble::column_to_rownames(var = "ID") #
colnames(mat_PA)[-c(1:2)] <- list_species

matPAfull<-as.data.frame(mat_PA)
row_names <- rownames(matPAfull)
matPAfull$grids<-rownames(matPAfull)
for (x in 1:11775) { #dont forget to change the total number of row numbers
  grd_coord@polygons[[x]]@ID <- row_names[x]
}
spdf <- SpatialPolygonsDataFrame(grd_coord, matPAfull)


richness <- mat_PA %>%
  mutate(sum = rowSums(dplyr::select(., -c(X, Y)))) %>%
  dplyr::select(X, Y, sum)

countries <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  st_geometry() %>%
  st_crop(grd)
countries <- st_as_sf(countries, crs = proj84) 

gridmap <- ggplot() + 
  geom_raster(data = SES_PD1, aes(x = X, y = Y, fill = ))+
  geom_sf(data = countries) + 
  scale_fill_viridis_c(
    option = 'turbo', 
    na.value = alpha("transparent", 0), 
    alpha = 1, 
    name = 'Species Richness', 
    direction=1, 
    breaks = c(5, 50, 95, 140)
  ) +
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
  theme(axis.title = element_blank(), legend.title = element_text(size = 10))
gridmap

#gridmap %>% ggsave('gridBS001x001.png',.,dpi = 320, device='png',width=18,height=12,units='cm')
#saveRDS(mat_PA, file = "mat_PA05.RDS")

###========= CLEANING COMMUNITY MATRIX FOR BIODIVERSITY ANALYSES ============###
mat_PA<-matPAfull[, -190]
mat_PA <- mat_PA[rowSums(mat_PA[,-c(1:2)]) != 0,]
sums <- colSums(mat_PA)
zero_cols <- sums == 0
if (any(zero_cols)) {
  cat("The following columns sum to 0:", paste(names(mat_PA)[zero_cols], collapse = ", "), "\n")
} else {
  cat("No columns sum to 0.\n")
}
mat_PA$Gobius_bucchichi <- NULL #remove species

#These need to be changed based on the species names used in the phylogenetic trees:
colnames(mat_PA)[colnames(mat_PA) == "Alosa_caspia"] <- "Alosa_caspia_caspia"
colnames(mat_PA)[colnames(mat_PA) == "Chelon_auratus"] <- "Liza_aurata"
colnames(mat_PA)[colnames(mat_PA) == "Planiliza_haematocheilus"] <- "Liza_haematocheila"
matPA<-mat_PA[ , -c(1:2)]
matPA <- matPA[, order(colnames(matPA))]
fishPA<-t(matPA)
matPA<-as.matrix(matPA)
matPA <- Matrix(data = matPA, sparse = TRUE)
save(matPA, mat_PA, fishPA, matPAfull, spdf, file="res01x01.RData")

###========================= LOAD AND PROCESS TREES =========================###
#Load and Process Trees
treesBS <- readRDS("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis/Phylo/Trees/Trees.Raboski.Actino.Stein.Chondry.203Sp.100.tr.RDS")
tips<-rownames(fishPA)
pruned<-lapply(treesBS,function(treesBS,tips)
  drop.tip(treesBS,setdiff(treesBS$tip.label,tips)),tips=tips)
class(pruned)<-"multiPhylo"
class(pruned)
unique_trees<-unique(pruned, incomparables = FALSE,
                     use.edge.length = FALSE,
                     use.tip.label = TRUE)
pruned_fish<-AllTipLabels(unique_trees)
SpToRemove <- setdiff(unique_trees[[1]]$tip.label, rownames(fishPA)) # species list to remove from the trees because not in the table.--> should be empty
goodtrees <- lapply(seq(1,length(unique_trees)), function(i){
  drop.tip(unique_trees[[i]], SpToRemove) #removing of the species of the trees that are not in the datas fishPA
})
rownames(fishPA[order(rownames(fishPA)), ] ) ==pruned_fish[order(pruned_fish)]
all(rownames(fishPA[order(rownames(fishPA)), ]) == pruned_fish[order(pruned_fish)])
t<-phytools::consensus.edges(goodtrees,method="least.squares") #use this one
#rm(treesBS, pruned, unique_trees)
#rm(fishPA, goodtrees, mat_PA, matPAfull)

###======================= LOAD AND PROCESS IUCN DATA =======================###
load("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis/Data/IUCNinfo.RData")
IUCNinfo <- species_df[,-c(18:24,27:30)]
all(IUCNinfo$X[order(IUCNinfo$X)] == pruned_fish[order(pruned_fish)]) #IF FALSE, run NEXT LINE
IUCNinfo <- IUCNinfo %>%
  mutate(X = case_when(
    X == "Alosa_caspia" ~ "Alosa_caspia_caspia",
    X == "Chelon_auratus" ~ "Liza_aurata",
    X == "Planiliza_haematocheilus" ~ "Liza_haematocheila",
    TRUE ~ X
  ))
all(IUCNinfo$X[order(IUCNinfo$X)] == pruned_fish[order(pruned_fish)])
IUCNbs <- IUCNinfo[, c("X", "category")]

#=== === === === === === === SAVE PART 1 DATA=== === === === === === === === ===#
source("Scripts/Sources/CreateSpatialPolygonBioregions.R")
BSRegions<-tibble::column_to_rownames(BSRegions, var = "grids")
BSRegions<-as.matrix(BSRegions)
BSRegions<-Matrix(BSRegions)
BSrows<-sparse2long(BSRegions)
duplicated_rows <- duplicated(BSrows$grids)
BSrows<- BSrows[!duplicated_rows, ]

save(BSRegions, spdf, matPA, t, goodtrees, fishPA, mat_PA, file="Part1.RData")
save(BSrows, spdf, matPA, t, goodtrees, file="part1.RData")


