library(knitr)

# Define your libraries here
my_libraries <- c("dplyr", "tidyr", "reshape2", "readxl", "ggplot2", "Matrix", "stringr", "ggpubr", 
                                 "raster", "rgeos", "rgdal", "sf", "sp", "shape", "geometry", "geosphere", 
                                 "lwgeom", "fasterize", "terra", "rnaturalearthdata", "rnaturalearth", 
                                 "ggspatial", "prettymapr", 
                                 "ape", "PhyloMeasures", "tidytree", "TreeTools", "phyloregion", "picante", 
                                 "mFD", "hillR", "readr", "magrittr", "rstatix", "FD", "vegan")

knitr::write_bib(my_libraries)
