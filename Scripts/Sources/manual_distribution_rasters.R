################################################################################
###                          Range mapping script                            ###
###           Stefanie Mayer, modified from Camille Albouy's script          ###
###                               06/2023                                    ###
################################################################################

### Function: conv_function(): Create a convex hull polygon from a set of coordinates
conv_function  <- function (x, polygons = TRUE, proj = NULL, g = NULL, buffer = 1) {
  
  # If there are less than or equal to three unique points, create a buffer around the points
  if (nrow(unique(round(x, 2))) <= 3) {  
    sp_coord <- SpatialPoints(x, proj4string = proj)
    return(gBuffer(sp_coord, width = buffer)@polygons[[1]]) 
  } else {
    # Otherwise, create a convex hull polygon from the points
    convex <- x[unique(as.vector(convhulln(x, option = "FA")$hull)),] #coordinates of convex points
    convex <- apply(convex, 2, as.numeric)
    # Arrange the points in clockwise order
    coord_conv <- convex[order(-1 * atan2(convex[,2] - mean(range(convex[,2])),
                                          convex[,1] - mean(range(convex[,1])))),] 
    # If polygons is TRUE, return the convex hull as a SpatialPolygons object
    if (polygons == TRUE) {
      # Add the first point to the end to close the polygon
      coord_conv <- rbind(coord_conv, coord_conv[1,])
      # Create a Polygons object from the convex hull coordinates
      P1 <- Polygons(srl=list(Polygon(coord_conv, hole = FALSE)), ID = "PolygA")
      # Create a SpatialPolygons object from the Polygons object
      P1 <- SpatialPolygons(Srl=list(P1), proj4string = proj)
      cat("### End of intersection for region : ", g, "\n") 
      # Buffer the resulting polygon if desired and return the first polygon in the buffer
      return(gBuffer(P1, width = buffer)@polygons[[1]]) #Enlarged polygons of 1°
    } else {
      # If polygons is FALSE, return the convex hull coordinates
      cat("### End of intersection for region with a<4 : ", g, "\n") 
      return(coord_conv)  #polygons coordinates
    } # end of ifelse
  }# end of ifeslse
} # end of conv_function  

#this is the function (before is internal )
get_gbif2rast_data <- function(Occ = Occ, proj = proj84, Bioreg = Bioreg, Bioreg_bf = Bioreg_bf, 
                               Species_name = "Capros_aper", name_shp = "UNION", bathy_sp = c(-10,120), 
                               corrected_bathy = TRUE, bathy = bathy) {
  # Print some information
  cat("########### Computation for species: ", Species_name, " #######################", "\n")
  cat("### Cleaning occurrences ...", "\n")   
  # Select coordinates columns and rename them
  occ_coord <- data.frame(Occ, stringsAsFactors = FALSE)[, c(1:2)]
  colnames(occ_coord) <- c("X", "Y")
  # Return if the number of occurrences is less than 5
  if (nrow(occ_coord) < 5) {
    return("Species in Gbif but not enough occurrences")
  } else {
    # Remove coordinates (0, 0)
    cond <- which(occ_coord[, 2] == 0 & occ_coord[, 1] == 0)
    if (length(cond) != 0) {
      occ_coord <- occ_coord[-which(occ_coord[, 2] == 0 & occ_coord[, 1] == 0), ]
    }
    # Transform coordinates into sf object
    occ_points <- st_as_sf(SpatialPoints(coords = occ_coord)) 
    occ_points$id <- seq(1, nrow(occ_coord), 1)
    sf::sf_use_s2(F) 
    st_crs(occ_points) <- proj 
    # Print information
    cat("### Projection adjustment for bioregion shapefile...", "\n") 
    # Crop bioregion shapefiles to the extent of occurrences
    BioregC_bf <- st_crop(st_as_sf(Bioreg_bf), extent(occ_points)) # corresponds to the right zone
    BioregC <- st_crop(st_as_sf(Bioreg), extent(occ_points))  
    # Print and remove unused variables
    print(get("BioregC"))
    print(name_shp)
    uniqu <- unique(get("BioregC")[[name_shp]])  #rename unique to uniqu (to not have it twice)
    gc() 
    # Print information
    cat("### Interscetion between occurences of ",Species_name," and bioregions ...", "\n") 
    #Using multiple bioregions will make computation time faster
    # SP_dist() function: Calculate intersection between occurrences and bioregions
    SP_dist <- lapply(1:length(uniqu), function(g) {
      cat("g=", g, '\n')
      tmp <- BioregC[get("BioregC")[[name_shp]] == uniqu[g],]
      if (st_geometry_type(tmp) != "MULTIPOLYGON"){ 
        tmp <- st_collection_extract(tmp, type="POLYGON") 
      }
      tmp <- as_Spatial(tmp)
      occ <- as(occ_points,"Spatial") 
      proj4string(tmp)<- proj
      #buffer the regions 
      tmp_bf <- BioregC_bf[get("BioregC_bf")[[name_shp]] == uniqu[g],] #inshore margin of 0.5° to include estuaries, ...
      # if the buffer is not a MULTIPOLYGON object, convert it to POLYGON
      if (st_geometry_type(tmp_bf) != "MULTIPOLYGON"){ 
        tmp_bf <- st_collection_extract(tmp_bf, type="POLYGON") 
      }
      # convert the buffer to a Spatial object and assign a projection
      tmp_bf <- as_Spatial(tmp_bf)        
      proj4string(tmp_bf)<- proj
      # find the occurrence data points within the buffer
      cel2get <- which(!is.na(over(occ,tmp_bf)[,1])) #takes points with 0.5 margin
      rm(occ) ; rm(tmp_bf) #remove objects to free up memory and prevent memory errors.
      # extract the occurrence data coordinates and round them to three decimal places
      a <- data.frame(round(occ_coord[cel2get,],3))
      rm(cel2get) #remove object
      # if the region is the Pacific ocean, do the following
      if (uniqu[g]=="South Pacific Ocean" |  uniqu[g]=="North Pacific Ocean"){
        cat("g=", g, " = Pacific ocean", "\n")
        b <- a # copy the occurrence data coordinates into b
        # adjust coordinates that are in both the West and East Pacific
        if (length(table(b[,1] < 0)) > 1) {
          b[which(b[,1] < 0),1] <- b[which(b[,1] < 0),1] + 360
        }
        # if there are no points in b, return NA
        if (nrow(b) == 0) {
          cat("g=",g, " ; a=0",'\n')
          return(NA)
        } #end of if
        # if there are 1 or 2 points in b, create a buffer around them
        if (nrow(b) == 1 | nrow(b) == 2) {
          cat("b= 1 or 2",'\n')
          bf <- gBuffer(SpatialPoints(coords=data.frame(b),proj4string=proj),width=1)
          df <- data.frame(ID=character(), stringsAsFactors=FALSE )
          for(i in bf@polygons ){ df <- rbind(df, data.frame(ID=i@ID, stringsAsFactors=FALSE)) }
          row.names(df) <- df$ID 
          bf_df <- SpatialPolygonsDataFrame(bf, df)
          return(bf_df)
        } #end of if 
        # if there are more than 2 points in b, create a polygon and divide it into sub-polygons
        if (nrow(b) > 2) {
          # create a polygon from the set of points b
          Poly <- SpatialPolygons(Srl=list(conv_function(b,proj=proj,g=g, buffer=0.000001)),proj4string=proj84)
          # create a line that goes from the left to the right of the map (180 degrees east to 180 degrees west, at the equator)
          line <- SpatialLines(list(Lines(list(Line(cbind(c(180,180),c(90,-90)))), ID="line")), proj4string=proj84)
          # check if there is an intersection between the polygon and the line
          if (class(gIntersection(Poly, line)) == "NULL"){  
            # if there is no intersection, buffer the polygon to create a buffer zone around it
            Poly <- gBuffer(Poly,width=1) #buffer in two step to avoid intersection with the line if longitude occurrences ~179.5 ...
            # return the intersection between the buffer and the input point set tmp
            return( as(st_intersection(st_as_sf(Poly),st_as_sf(tmp)), "Spatial") )
          } 
          # if there is an intersection between the polygon and the line, divide the polygon into sub-polygons
          cat("Occurrences in both side of pacific", "\n")
          # buffer the polygon to create a buffer zone around it
          Poly <- gBuffer(Poly,width=1) # avoid intersection with the line if longitude occurrences ~179.5 ...
          lpi <- gIntersection(Poly, line)        # intersect your line with the polygon
          blpi <- gBuffer(lpi, width = 0.000001)    # create a small buffer around the line intersection
          dpi <- gDifference(Poly, blpi) # subtract the small buffer from the polygon to split it into sub-polygons
          Pol_list <- list() # empty list 
          # Loop through each sub-polygon and extract the intersection with the input point set tmp
          for (i in 1:length(dpi@polygons[[1]]@Polygons)) {
            cat("i=", i, "\n") # extract the i-th sub-polygon
            Pol <- SpatialPolygons(list(Polygons(list(dpi@polygons[[1]]@Polygons[[i]]), as.character(i)))) 
            # if the sub-polygon crosses the antimeridian (180 degrees east/west), shift its coordinates by 360 degrees
            if (sum(Pol@polygons[[1]]@Polygons[[1]]@coords[,1] > 180) > 0) {
              Pol@polygons[[1]]@Polygons[[1]]@coords[,1] <- Pol@polygons[[1]]@Polygons[[1]]@coords[,1] - 360
              Pol@bbox <- as.matrix(extent(Pol@polygons[[1]]@Polygons[[1]]@coords), byrow)
            }
            proj4string(Pol) <- proj
            # extract the intersection between the sub-polygon and the input point set tmp
            Res <- st_intersection(st_as_sf(Pol), st_as_sf(tmp))
            if (nrow(Res) == 0) {
              return(NA) #or return NA if no intersection
            }
            if (str_detect(st_geometry_type(Res), "POLYGON") == F) {
              Res <- st_collection_extract(Res, type = "POLYGON")
            }
            Pol_list[[i]] <- as(Res, "Spatial")
          }
          return(Pol_list)
        } #end of if(nrow(b)>2 )
      }   else {
        
        # If there are no points within the buffer, return NA
        if (nrow(a) == 0) {
          cat("g=",g, " ; a=0",'\n') 
          return(NA)
        } 
        # If there is only 1 or 2 points within the buffer
        if (nrow(a) == 1 | nrow(a) == 2) { 
          cat("a= 1 or 2",'\n') 
          # Create a buffer around the points and return as a SpatialPolygonsDataFrame
          bf <- gBuffer(SpatialPoints(coords=data.frame(a), proj4string=proj), width=1)
          df <- data.frame(ID=character(), stringsAsFactors=FALSE )
          for(i in bf@polygons ){ 
            df <- rbind(df, data.frame(ID=i@ID, stringsAsFactors=FALSE))  
          }
          row.names(df) <- df$ID
          bf_df <- SpatialPolygonsDataFrame(bf, df)
          return(bf_df)
        } 
        # If more than 2 points,
        if (nrow(a) > 2) { 
          cat("g=",g, " ; a>2",'\n')
          # Creates a polygon using the conv_function() function and a buffer of 1 degree.
          Poly <- SpatialPolygons(Srl=list(conv_function(a, proj=proj, g=g, buffer=1)), proj4string=proj84)
          Res <- st_intersection(st_as_sf(Poly), st_as_sf(tmp)) # makes intersection with the precise coast shapefile
          if (nrow(Res)==0) { 
            return(NA) 
          } #in cases points are taken in the buffer but not in tmp
          #If  the intersection is not a polygon, the polygon is extracted from the intersection 
          if (str_detect(st_geometry_type(Res), "POLYGON")==F){ 
            Res <- st_collection_extract(Res, type="POLYGON") 
          }
          Res <- as(Res, "Spatial") # and converted into a  SpatialPolygons 
          return(Res)
        } # end of if a>2
      }  # end of else
    })  # end of SP_dist function
    cat("######## End of intersection ###########", "\n")
    L <- unlist(SP_dist[!is.na(SP_dist)])
    rm(SP_dist)
    gc()
    if(length(L) == 0){
      return("Species not in the considered area")
    }
    # Create spatial polygons from L and set the ID of each polygon in L to its respective name
    names_poly <- paste0("Poly_", seq(1, length(L), 1))
    for(i in 1:length(L)){
      L[[i]]@polygons[[1]]@ID <- names_poly[i]
    }
    # Create spatial polygons from L
    shp_species <- SpatialPolygons(Srl = lapply(L, function(y) {y@polygons[[1]]}), proj4string = proj) 
    
    ############ Bathymetry correction ###############
    cat("### Bathymetry correction for species: ", Species_name, " ###", "\n")
    sf_species <- st_as_sf(shp_species) # convert shp_species to an sf object
    r1 <- crop(bathy, shp_species) # crop the bathymetry raster to the species polygon
    raster_sp <- fasterize(sf_species, r1) # rasterize polygons (high performance function)
    if (corrected_bathy == TRUE){
      data_sp <- values(raster_sp) # extract raster values from the species raster
      data_r1 <- values(r1) # extract raster values from the bathymetry raster
      data_r1[which(is.na(data_sp))] <- NA # remove bathymetry values that are out of the species polygon
      data_r1[which(data_r1 > -bathy_sp$Depth_min)] <- NA # remove values that are out of the species bathymetric range
      data_r1[which(data_r1 <= -bathy_sp$Depth_max)] <- NA 
      values(raster_sp) <- data_r1 # replace species raster values with corrected bathymetry values
    }
    # can change 15 to something else
    #raster_sp <- aggregate(raster_sp, fact = 2.5) # Species range raster of 0,1° resolution -> need lots of RAM for global species
    val <- values(raster_sp)
    val[!is.na(val)] <- 1
    val[is.na(val)] <- 0
    values(raster_sp) <- val # presence/absence raster
    cat("########### End of computation for species: ", Species_name, " #######################", "\n")
    return(raster_sp)
  } 
} # end of function get_gbif2rast_data 

################

find.5.min <- function(i){i[order(i, decreasing=TRUE)][1:4]}
`%ni%` <- Negate(`%in%`)