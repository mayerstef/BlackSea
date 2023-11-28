################################################################################
####                                                                        ####
####         Master thesis:LINKING PATTERNS IN PHYLOGENY, TRAITS,           ####
####            ABIOTIC VARIABlES AND SPACE: BLACK SEA FISHES               ####
####                                                                        ####
####                               PART 2                                   ####
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

# Phylogenetic Diversity Analysis
lib_phylo <- c("ape", "PhyloMeasures", "tidytree", "TreeTools", "phyloregion")

#Functional Diversity Analysis
lib_fun <- c("mFD", "hillR", "readr", "magrittr", "vegan")

#Choose which libraries you would like to load
sapply(lib_general,library,character.only=TRUE)
sapply(lib_gis,library,character.only=TRUE)
sapply(lib_phylo,library,character.only=TRUE)
sapply(lib_fun,library,character.only=TRUE)

#Load Part 1 Data
load("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis/BACKUPs/Part1.RData")

#========================== Diversity Pipelines ===============================#

#For Analyses with 0.5x0.5 resolution:
source("./Scripts/Pipelines/TaxonomicDiversity.R")
taxodiversity<-taxonomic_diversity(tab_eDNA = fishPA, abundance = F, has.replicates = F)

fishPA<-t(as.matrix(matPA))
source("./Scripts/Pipelines/functionstoLOAD.R")
source("./Scripts/Pipelines/PhyloDiversity.R")
PhylogeneticDiversity01<-phylo_diversity(trees = goodtrees, 
                                         tab_eDNA = fishPA,
                                         savepathindic = "./Phylo/outputs01x01/",  
                                         abundance=FALSE,
                                         has.replicates = F)
# Set the folder path
folder_path <- "./Phylo/outputs01x01/"
file_paths <- paste0(folder_path, "Index_results_Tree_", 1:91, ".rds")
data_list <- list()
for (file_path in file_paths) {
  data <- readRDS(file_path)
  data_list[[file_path]] <- data
}
# Calculate the mean and min of the data
mean_data <- Reduce("+", data_list) / length(data_list)
min_data <- Reduce(pmin, data_list)
write_rds(mean_data, file="./Phylo/outputs01x01/Index_results_Tree_mean.rds")
write_rds(min_data, file="./Phylo/outputs01x01/Index_results_Tree_min.rds")

#========================= Phylogenetic Diversity =============================#
#========================== Phyloregion Package ===============================#
#Phylogenetic diversity 
mypd <- PD(matPA, t)
M <- merge(spdf, data.frame(grids=names(mypd), pd=mypd), by="grids")
M <- M[!is.na(M@data$pd),]

#Weighted endemism 
Endm <- weighted_endemism(matPA)
m <- merge(spdf, data.frame(grids=names(Endm), WE=Endm), by="grids")
m <- m[!is.na(m@data$WE),]

#Phylogenetic endemism 
pe <- phylo_endemism(matPA, t)
mx <- merge(spdf, data.frame(grids=names(pe), pe=pe), by="grids")
mx <- mx[!is.na(mx$pe),]

#Evolutionary Distinctiveness and Global Endangerment
source("./Scripts/Sources/EDGE_updated.R")
unique(IUCNbs$category)
EDGE <- EDGE(IUCNbs, t, Redlist = "category", species="X")
y <- map_trait(matPA, EDGE, FUN = sd, shp=spdf)

#Phyloregions show evolutionary affinities among disjunct assemblages (function plot.phyloregion). 
pb <- phylobeta(matPA, t)
phyreg <- phyloregion(pb[[1]], shp=spdf)

#Ordination of phyloregions in NMDS space shows that different phyloregions differ strongly in evolutionary uniqueness (function plot_NMDS). 
par(mar=rep(4,4))
plot_NMDS(phyreg, cex=3)
text_NMDS(phyreg)

#===================== CHECK FOR HOTSPOTS/COLDSPOTS ========================###
##Endm 
C <- coldspots(Endm) # coldspots
H <- hotspots(Endm) # hotspots
DF <- data.frame(grids=names(C), cold=C, hot=H)
HC <- merge(spdf, DF, by = "grids", all = TRUE)
plot(spdf, border = "grey", col = "lightgrey",
     main = "Weighted Endemism Hotspots and Coldspots")
subset_cold <- !is.na(HC@data$cold) & HC@data$cold == 1
cold_data <- HC[subset_cold, ]
plot(cold_data, col = "#057985", add = TRUE, border = NA)
subset_hot <- !is.na(HC@data$hot) & HC@data$hot == 1
hot_data <- HC[subset_hot, ]
plot(hot_data, col = "#ae2012ff",, add = TRUE, border = NA)
legend("bottomleft", fill = c("#057985", "#ae2012ff")
       legend = c("coldspots", "hotspots"), bty = "n", inset = .092)
hotcoldWE<-hot_data@data


#Generate a sparse community matrix as input for clustering regions based on the similairity of
longPA<-sparse2long(matPA)
longPA <- longPA %>% select(ncol(longPA), everything())
BSRegions_M<- BSRegions %>% 
  tibble::column_to_rownames("grids") %>% 
  as.matrix() %>% 
  Matrix()
BSRegions_long<-sparse2long(BSRegions_M)

#Phylogenetic diversity standardized for species richness
SES_PD1<-PD_ses(matPA, t, model = "tipshuffle", reps = 1000)
SES_PD2<-PD_ses(matPA, t, model = "rowwise", reps = 1000)
# this one doesnt work gives same as PD 
#SES_PD3<-PD_ses(matPA, t, model = "colwise", reps = 1000) 
SES_PD12<-full_join(SES_PD1, SES_PD2, by = c("grids", "richness", "PD_obs", "reps"), 
suffix = c("_ts", "_rw")) %>% relocate("reps", .after = everything(.))

SES_PD <- full_join(BSRegions, SES_PD12) %>%
  na.omit()
 
SES_PD$X <- sub("_.*", "", SES_PD$grids)  # Extract the X values before the underscore
SES_PD$Y <- sub(".*_", "", SES_PD$grids)  # Extract the Y values after the underscore
SES_PD$X <- as.numeric(SES_PD$X)
SES_PD$Y <- as.numeric(SES_PD$Y)
#write.csv(SES_PD, file="SES_PD.csv")

########################### BETA DIVERSITY #####################################
#tax1 <- beta_core(matPA)
#tax2 <- beta_diss(matPA, index.family = "sorensen")
#tax3 <- beta_diss(matPA, index.family = "jaccard")
#save(tax1, tax2, tax3, file="betatax.RData")

#SES_Pbeta1<-phylobeta_ses(matPA, t, model = "tipshuffle", reps = 1000)
#SES_Pbeta2<-phylobeta_ses(matPA, t, model = "rowwise", reps = 1000)
#SES_Pbeta3<-phylobeta_ses(matPA, t, model = "colwise", reps = 1000)


#=========================== FUNCTIONAL DIVERSITY =============================##
#=========================== using the mFD package =============================#
#See "FishBase_traits.R" for the creation of Traits186 dataframe
#incase you're jumping to this part 
load("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis/Data/Part1.RData")
rm(BSRegions, fishPA, mat_PA, goodtrees,   Traits186)
fishPA<-as.matrix(matPA)
fishPA <- fishPA[rowSums(fishPA) >= 6,]#to be able to calculate FRic in 5D 
rm(matPA)

Traits186 <- read_csv("Data/Traits186.csv", col_types = cols(...1 = col_skip()))
fish_traits<-Traits186[, -c(2:8, 13,16, 20:35)]
column_names <- colnames(fish_traits[,-1])
column_classes <- sapply(fish_traits[,-1], class)

fish_traits_cat <- data.frame(
  trait_name = column_names,
  trait_type = ifelse(toupper(substr(column_classes, 1, 1)) == "N", "Q", "N"),
  stringsAsFactors = FALSE
)
fish_traits <- tibble::column_to_rownames(fish_traits, "X")

fish_traits$position <- as.factor(fish_traits$position)
fish_traits$shape <- as.factor(fish_traits$shape)
#fish_traits$morph <- as.factor(fish_traits$morph)
fish_traits$fert <- as.factor(fish_traits$fert)
fish_traits$repro <- as.factor(fish_traits$repro)
fish_traits$parent <- as.factor(fish_traits$parent)

Classif <- Traits186[, c(1:6)]
Classif <- Classif %>% dplyr::rename(., Taxon = X) %>% 
  relocate(Taxon, .after=Genus)
row.names(Classif)<-Classif$Taxon

# Display the table:
knitr::kable(head(fish_traits),
             caption = "Species x traits data frame")
# Display the table:
knitr::kable(head(fish_traits_cat), 
             caption = "Traits types based on **fish traits** dataset")

# Species traits summary:
fish_traits_summ <- mFD::sp.tr.summary(
  tr_cat     =fish_traits_cat,   
  sp_tr      =fish_traits, 
  stop_if_NA = TRUE)

fish_traits_summ
fish_traits_summ$"tr_types"                    
fish_traits_summ$"mod_list"   
fish_traits_summ$"tr_summary_list"   

# Summary of the assemblages * species dataframe:
asb_sp_fish_summ <- mFD::asb.sp.summary(asb_sp_w = fishPA)


head(asb_sp_fish_summ$"asb_sp_occ", 3)        # Species occurrences for the first 3 assemblages
asb_sp_fish_occ <- asb_sp_fish_summ$"asb_sp_occ"
asb_sp_fish_summ$"sp_tot_w"              # Species total biomass in all assemblages
asb_sp_fish_summ$"asb_tot_w"             # Total biomass per assemblage
asb_sp_fish_summ$"asb_sp_richn"           # Species richness per assemblage
asb_sp_fish_summ$"asb_sp_nm"[[1]]             # Names of species present in the first assemblage

sp_dist_fish <- mFD::funct.dist(
  sp_tr         =fish_traits,
  tr_cat        =fish_traits_cat,
  metric        = "gower", #non-continuous traits
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

fun_dist<-as.matrix(sp_dist_fish)
summary(as.matrix(sp_dist_fish))
#ranges from 
min(sp_dist_fish) #0.0002788234
max(sp_dist_fish) #0.8791891

#Looking at the quality of the functional space
fspaces_quality_fish <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fish,
  fdendro             = "average",
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute", "squared"),
  fdist_scaling       = c(TRUE, FALSE))

# display the table gathering quality metrics:
fspaces_quality_fish$"quality_fspaces"
round(fspaces_quality_fish$"quality_fspaces", 3)     
# retrieve the functional space associated with minimal quality metric: 
apply(fspaces_quality_fish$quality_fspaces, 2, which.min)

round(fspaces_quality_fish$"quality_fspaces", 2) 
#if you round to 2, 5 is best in every category 
# the more parsimonious (i.e., simpler) model is often preferred, following the principle of Occam's Razor.


fspaces_quality_fish$"quality_fspaces" %>%
  tibble::as_tibble(rownames = "Funct.space") %>%
  tidyr::pivot_longer(cols =! Funct.space, names_to = "quality_metric", values_to = "Quality") %>%
  ggplot2::ggplot(ggplot2::aes(x = Funct.space, y = Quality, 
                               color = quality_metric, shape = quality_metric)) +
  ggplot2::geom_point() 

#Honestly, from this 3 is still under average, so 3, 4, or 5 is probably ok 

mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_fish,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")


#this shows 5 looks best, but 4 seems okay as well. 
#many of the distances between pairs of species on the dendrogram are more than 0.3 different 
#from the distances based on their traits. This means that some species that have similar trait 
#values actually have large distances on the dendrogram. The dendrogram has a binary structure, 
#which means that many pairs of species have the same distance. In particular, all pairs of species
#that are on different sides of the root of the tree have the maximum distance.

# pick 3 random fish 
cherry<-"Serranus_hepatus"
lime<-"Dentex_macrophthalmus"
lemon<-"Acipenser_ruthenus"
fish_traits[c(cherry, lime, lemon), ]

fspaces_quality_fish$"details_fspaces"$"pairsp_fspaces_dist" %>%
  dplyr::filter(sp.x %in% c(cherry, lime, lemon) & 
                  sp.y %in% c(cherry, lime, lemon)) %>%
  dplyr::select(sp.x, sp.y, tr, pcoa_5d, tree_average) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), round, 2))

#sp.x                  sp.y   tr pcoa_5d tree_average
#1 Dentex_macrophthalmus    Acipenser_ruthenus 0.27    0.31         0.27
#2      Serranus_hepatus    Acipenser_ruthenus 0.35    0.38         0.29
#3      Serranus_hepatus Dentex_macrophthalmus 0.26    0.28         0.29\

# ===> because tree_average is lower than pcoa_5d, it shows that the dendogram
#      overestimates distance between some pairs of species having actually similar trait values.

#=============================================================================
#5. Test correlation between functional axes and traits

sp_faxes_coord_fish <- fspaces_quality_fish$"details_fspaces"$"sp_pc_coord"

source("./Scripts/Sources/check_inputs.R")
source("./Scripts/Sources/traits.faxes.cor.adj.R")
fish_tr_faxes <- traits.faxes.cor.adj(
  color_signif= "#005f73ff",
  color_non_signif = "#94d2bdff",
  sp_tr          =fish_traits, 
  sp_faxes_coord = sp_faxes_coord_fish[ , c("PC1", "PC2", "PC3", "PC4", "PC5")], 
  plot           = TRUE)

# Print traits with significant effect:
fish_tr_faxes$"tr_faxes_stat"[which(fish_tr_faxes$"tr_faxes_stat"$"adj.p.value" < 0.05), ]
fish_tr_faxes$"tr_faxes_stat"
# Return plots:
fish_tr_faxes$"tr_faxes_plot"
#--> after looking at this I think we only need the first 4, because 5 only add an 
#    explanation for "max", but it's already in PC1 (and better) 


#===============================================================================
#6. Plot functional space
sp_faxes_coord_fish <- fspaces_quality_fish$"details_fspaces"$"sp_pc_coord"

big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord_fish,
  faxes           = NULL,
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkgreen",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)


# Plot the graph with all pairs of axes:
big_plot$patchwork

sp_faxes_coord_fish[ , c("PC1", "PC2", "PC3", "PC4", "PC5")]

# doesnt work with FRic 
alpha_fd_indices_fish <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_fish[ , c("PC1", "PC2", "PC3", "PC4", "PC5")],
  asb_sp_w         = fishPA,
  ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fdiv", "fori", "fric", 
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

save(alpha_fd_indices_fish, file="alpha_fun.RData")

fd_ind_values_fish <- alpha_fd_indices_fish$"functional_diversity_indices"

details_list_fish <- alpha_fd_indices_fish$"details"
details_list_fish

# If not loaded from above:
asb_sp_fish_summ <- mFD::asb.sp.summary(asb_sp_w = fishPA)
asb_sp_fish_occ <- asb_sp_fish_summ$"asb_sp_occ"
sp_dist_fish <- mFD::funct.dist(
  sp_tr         =fish_traits,
  tr_cat        =fish_traits_cat,
  metric        = "gower", #non-continuous traits
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Compute alpha FD Hill with q = 0:
fish_FD0mean <- mFD::alpha.fd.hill(
  asb_sp_w = asb_sp_fish_occ, 
  sp_dist  = sp_dist_fish, 
  tau      = "mean", 
  q        = 0)

addFuncHill<-round(fish_FD0mean$"asb_FD_Hill", 2)
addFuncHill<- addFuncHill %>% 
  as.data.frame() %>% 
  rownames_to_column(var="grids")

fd_ind_values_fish <-fd_ind_values_fish %>% rownames_to_column(var="grids")
ALL<-full_join(fd_ind_values_fish, SES_PD)
ALL<-full_join(ALL, addFuncHill)

#write.csv(ALL, file = "ALL.csv", row.names = F)
################################################################################

IUCN<-Traits186[,c(1, 20)]
lookup <- c(0.001, 0.01, 0.1, 0.67, 0.999, 1, 1, NA_real_,NA_real_)
names(lookup) <- c("LC", "NT", "VU", "EN", "CR" , "EW", "EX", "DD", "NE")
IUCN$num <- lookup[IUCN$category]
IUCN<-IUCN %>% tibble::column_to_rownames(var="X")
fish_traits$IUCN_num<-IUCN$num
FUSE<-fuse(sp_dist_fish, as.matrix(sp_faxes_coord_fish), nb_NN = 5, fish_traits$IUCN_num, standGE = T)
FUSE<-FUSE %>% tibble::rownames_to_column("X")
IUCN<-IUCN %>% tibble::rownames_to_column("X")

EDGE <- EDGE(IUCN, t, Redlist = "category", species="X")
EDGE <- EDGE %>% as.data.frame() %>% 
  tibble::rownames_to_column("X") 
colnames(EDGE)[2] <- "EDGE"


FUDGE<-full_join(FUSE, EDGE) %>% na.omit()
#FUSE = 2  EDGE = 7

i<-2
x <- FUDGE[,i]
names(x)<-FUDGE[,1]
y_fuse <- map_trait(matPA, x, FUN = sd, shp=spdf)
fuse_df<-y_fuse@data

i<-7
x <- FUDGE[,i]
names(x)<-FUDGE[,1]
y_edge<- map_trait(matPA, x, FUN = sd, shp=spdf)
edge_df<-y_edge@data

#### go to Figures. Maps to plot
map<-ggplot() + 
  geom_tile(data = edge_df, aes(x = X, y = Y, fill = traits))+
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
                       #breaks = c(5, 50, 95, 140), #this is for looking at SR
                       name="EDGE") + #Don't forget to change the legend name
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
                                         colour ="black"))
map
ggsave("./Maps/EDGE2.png", map, width = 10, height = 6.5, dpi = 300, bg = "transparent")


map<-ggplot() + 
  geom_tile(data = fuse_df, aes(x = X, y = Y, fill = traits))+
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
                       #breaks = c(5, 50, 95, 140), #this is for looking at SR
                       name="FUSE") + #Don't forget to change the legend name
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
                                         colour ="black"))
map
ggsave("./Maps/FUSE2.png", map, width = 10, height = 6.5, dpi = 300, bg = "transparent")


map<-ggplot() + 
  geom_raster(data = hotcoldWE, aes(x = X, y = Y, fill = hot))+
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
                       #breaks = c(5, 50, 95, 140), #this is for looking at SR
                       name="FUSE") + #Don't forget to change the legend name
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
                                         colour ="black"))
map


matPA<-as.matrix(matPA)
dist <- vegdist(matPA, binary=T, method = "raup")
rm(unifrac_dist)
opt20<-optimal_phyloregion(dist, method="ward.D2", k=20)


#Phyloregions show evolutionary affinities among disjunct assemblages (function plot.phyloregion). 
pb <- phylobeta(matPA, t)
phyreg <- phyloregion(pb[[1]], shp=spdf)

#Ordination of phyloregions in NMDS space shows that different phyloregions differ strongly in evolutionary uniqueness (function plot_NMDS). 
par(mar=rep(4,4))
plot_NMDS(phyreg, cex=3)
text_NMDS(phyreg)
,



C <- coldspots(Endm) # coldspots
H <- hotspots(Endm) # hotspots
DF <- data.frame(grids=names(C), cold=C, hot=H)
HC <- merge(spdf, DF, by = "grids", all = TRUE)
plot(spdf, border = "grey", col = "lightgrey",
     main = "Weighted Endemism Hotspots and Coldspots")
subset_cold <- !is.na(HC@data$cold) & HC@data$cold == 1
cold_data <- HC[subset_cold, ]
plot(cold_data, col = "#057985", add = TRUE, border = NA)
subset_hot <- !is.na(HC@data$hot) & HC@data$hot == 1
hot_data <- HC[subset_hot, ]
plot(hot_data, col = "#ae2012ff",, add = TRUE, border = NA)
legend("bottomleft", fill = c("#057985", "#ae2012ff")
       legend = c("coldspots", "hotspots"), bty = "n", inset = .092)
hotcoldWE<-hot_data@data

