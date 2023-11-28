library(picante)
library(vegan)
library(FD)
library(dplyr)

fish_traits<-fish_traits[order(row.names(fish_traits)),]
species<-species[,order(colnames(species))]
null_com<-picante::randomizeMatrix(species, null.model = "trialswap", iterations=1000)

#calculate the funcitional diversity using the mFD package, 
#choosing 5 axes, but usingthe null_com instead of the fishPA

load("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis/NULLalpha_fun.RData")
load("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis/BACKUPs/alpha_fun.RData")

expectedFRic <- mean(nullFD$FRic)
null.sd<- sd(nullFD$FRic)
obsFRic<-obsFD$FRic

fun_indices<-alpha_fd_indices_fish[["functional_diversity_indices"]]
null_indices<-NULLalpha_fd_indices_fish[["functional_diversity_indices"]]

expected.fric <- mean(null_indices$fric)
null.sd<- sd(null_indices$fric)
obs.fric<-fun_indices$fric
fricSES <- (obs.fric- expected.fric) / null.sd
summary(fricSES)
fun_indices$fricSES <- fricSES

# Calculate fdisSES
expected_fdis <- mean(null_indices$fdis)
null.sd <- sd(null_indices$fdis)
obs_fdis <- fun_indices$fdis
fdisSES <- (obs_fdis - expected_fdis) / null.sd
fun_indices$fdisSES <- fdisSES

# Calculate fmpdSES
expected_fmpd <- mean(null_indices$fmpd)
null.sd <- sd(null_indices$fmpd)
obs_fmpd <- fun_indices$fmpd
fmpdSES <- (obs_fmpd - expected_fmpd) / null.sd
fun_indices$fmpdSES <- fmpdSES

# Calculate fnndSES
expected_fnnd <- mean(null_indices$fnnd)
null.sd <- sd(null_indices$fnnd)
obs_fnnd <- fun_indices$fnnd
fnndSES <- (obs_fnnd - expected_fnnd) / null.sd
fun_indices$fnndSES <- fnndSES

# Calculate feveSES
expected_feve <- mean(null_indices$feve)
null.sd <- sd(null_indices$feve)
obs_feve <- fun_indices$feve
feveSES <- (obs_feve - expected_feve) / null.sd
fun_indices$feveSES <- feveSES

# Calculate fdivSES
expected_fdiv <- mean(null_indices$fdiv)
null.sd <- sd(null_indices$fdiv)
obs_fdiv <- fun_indices$fdiv
fdivSES <- (obs_fdiv - expected_fdiv) / null.sd
fun_indices$fdivSES <- fdivSES

# Calculate foriSES
expected_fori <- mean(null_indices$fori)
null.sd <- sd(null_indices$fori)
obs_fori <- fun_indices$fori
foriSES <- (obs_fori - expected_fori) / null.sd
fun_indices$foriSES <- foriSES

# Calculate fspeSES
expected_fspe <- mean(null_indices$fspe)
null.sd <- sd(null_indices$fspe)
obs_fspe <- fun_indices$fspe
fspeSES <- (obs_fspe - expected_fspe) / null.sd
fun_indices$fspeSES <- fspeSES

# Calculate fide_PC1SES
expected_fide_PC1 <- mean(null_indices$fide_PC1)
null.sd <- sd(null_indices$fide_PC1)
obs_fide_PC1 <- fun_indices$fide_PC1
fide_PC1SES <- (obs_fide_PC1 - expected_fide_PC1) / null.sd
fun_indices$fide_PC1SES <- fide_PC1SES

# Calculate fide_PC2SES
expected_fide_PC2 <- mean(null_indices$fide_PC2)
null.sd <- sd(null_indices$fide_PC2)
obs_fide_PC2 <- fun_indices$fide_PC2
fide_PC2SES <- (obs_fide_PC2 - expected_fide_PC2) / null.sd
fun_indices$fide_PC2SES <- fide_PC2SES

# Calculate fide_PC3SES
expected_fide_PC3 <- mean(null_indices$fide_PC3)
null.sd <- sd(null_indices$fide_PC3)
obs_fide_PC3 <- fun_indices$fide_PC3
fide_PC3SES <- (obs_fide_PC3 - expected_fide_PC3) / null.sd
fun_indices$fide_PC3SES <- fide_PC3SES

# Calculate fide_PC4SES
expected_fide_PC4 <- mean(null_indices$fide_PC4)
null.sd <- sd(null_indices$fide_PC4)
obs_fide_PC4 <- fun_indices$fide_PC4
fide_PC4SES <- (obs_fide_PC4 - expected_fide_PC4) / null.sd
fun_indices$fide_PC4SES <- fide_PC4SES

# Calculate fide_PC5SES
expected_fide_PC5 <- mean(null_indices$fide_PC5)
null.sd <- sd(null_indices$fide_PC5)
obs_fide_PC5 <- fun_indices$fide_PC5
fide_PC5SES <- (obs_fide_PC5 - expected_fide_PC5) / null.sd
fun_indices$fide_PC5SES <- fide_PC5SES

# Calculate and print summary for each indice
# Calculate and print summary for each indice
cat("Summary for fdisSES:\n", summary(fun_indices$fdisSES), "\n")
cat("Summary for fmpdSES:\n", summary(fun_indices$fmpdSES), "\n")
cat("Summary for fnndSES:\n", summary(fun_indices$fnndSES), "\n")
cat("Summary for feveSES:\n", summary(fun_indices$feveSES), "\n")
cat("Summary for fdivSES:\n", summary(fun_indices$fdivSES), "\n")
cat("Summary for foriSES:\n", summary(fun_indices$foriSES), "\n")
cat("Summary for fspeSES:\n", summary(fun_indices$fspeSES), "\n")
cat("Summary for fide_PC1SES:\n", summary(fun_indices$fide_PC1SES), "\n")
cat("Summary for fide_PC2SES:\n", summary(fun_indices$fide_PC2SES), "\n")
cat("Summary for fide_PC3SES:\n", summary(fun_indices$fide_PC3SES), "\n")
cat("Summary for fide_PC4SES:\n", summary(fun_indices$fide_PC4SES), "\n")
cat("Summary for fide_PC5SES:\n", summary(fun_indices$fide_PC5SES), "\n")

save(alpha_fd_indices_fish, fun_indices, null_indices, NULLalpha_fd_indices_fish, file="FunctionalDiversityIndicesmFDwithSES.RData")

###################################################################################
#Load Part 1 Data
load("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis/BACKUPs/Part1.RData")
rm(BSRegions, fishPA, mat_PA, spdf, t)
null_com<-picante::randomizeMatrix(matPA, null.model = "trialswap", iterations=1000)

fishPA<-t(as.matrix(null_com))
source("./Scripts/Pipelines/functionstoLOAD.R")
source("./Scripts/Pipelines/PhyloDiversity.R")
PhylogeneticDiversity01<-phylo_diversity(trees = goodtrees, 
                                         tab_eDNA = fishPA,
                                         savepathindic = "./Phylo/nulloutputs01x01/",  
                                         abundance=FALSE,
                                         has.replicates = F)
# Set the folder path
folder_path <- "./Phylo/nulloutputs01x01/"
file_paths <- paste0(folder_path, "Index_results_Tree_", 1:91, ".rds")
data_list <- list()
for (file_path in file_paths) {
  data <- readRDS(file_path)
  data_list[[file_path]] <- data
}
# Calculate the mean and min of the data
mean_nulldata <- Reduce("+", data_list) / length(data_list)
min_nulldata <- Reduce(pmin, data_list)

library(readr)
write_rds(mean_nulldata, file="./Phylo/nulloutputs01x01/Index_results_Tree_mean.rds")
write_rds(min_nulldata, file="./Phylo/nulloutputs01x01/Index_results_Tree_min.rds")

obs_phylo<- readRDS("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis/Phylo/outputs01x01/Index_results_Tree_mean.rds")
sd_nulldata <- sqrt(Reduce("+", lapply(data_list, function(x) (x - mean_nulldata)^2)) / (length(data_list) - 1))

null_indices<-mean_nulldata
phylo_indices<-obs_phylo

expected_PD <- mean(null_indices$PD)
null_sd_PD <- sd(null_indices$PD)
phylo_indices$PD.SES <- (phylo_indices$PD - expected_PD) / null_sd_PD

expected_NTI <- mean(null_indices$NTI)
null_sd_NTI <- sd(null_indices$NTI)
phylo_indices$NTI.SES <- (phylo_indices$NTI - expected_NTI) / null_sd_NTI

expected_NRI <- mean(null_indices$NRI)
null_sd_NRI <- sd(null_indices$NRI)
phylo_indices$NRI.SES <- (phylo_indices$NRI - expected_NRI) / null_sd_NRI

expected_MPD <- mean(null_indices$MPD)
null_sd_MPD <- sd(null_indices$MPD)
phylo_indices$MPD.SES <- (phylo_indices$MPD - expected_MPD) / null_sd_MPD

expected_MNTD <- mean(null_indices$MNTD)
null_sd_MNTD <- sd(null_indices$MNTD)
phylo_indices$MNTD.SES <- (phylo_indices$MNTD - expected_MNTD) / null_sd_MNTD

expected_VPD <- mean(null_indices$VPD)
null_sd_VPD <- sd(null_indices$VPD)
phylo_indices$VPD.SES <- (phylo_indices$VPD - expected_VPD) / null_sd_VPD

expected_VNTD <- mean(null_indices$VNTD)
null_sd_VNTD <- sd(null_indices$VNTD)
phylo_indices$VNTD.SES <- (phylo_indices$VNTD - expected_VNTD) / null_sd_VNTD

# Calculate and print summary for each indice in phylo_indices
cat("Summary for PD.SES:\n", summary(phylo_indices$PD.SES), "\n")
cat("Summary for NTI.SES:\n", summary(phylo_indices$NTI.SES), "\n")
cat("Summary for NRI.SES:\n", summary(phylo_indices$NRI.SES), "\n")
cat("Summary for MPD.SES:\n", summary(phylo_indices$MPD.SES), "\n")
cat("Summary for MNTD.SES:\n", summary(phylo_indices$MNTD.SES), "\n")
cat("Summary for VPD.SES:\n", summary(phylo_indices$VPD.SES), "\n")
cat("Summary for VNTD.SES:\n", summary(phylo_indices$VNTD.SES), "\n")

cat("Summary for MPD.SES:\n", summary(phylo_indices$MPD.SES), "\n")
cat("Summary for MNTD.SES:\n", summary(phylo_indices$MNTD.SES), "\n")
cat("Summary for NTI.SES:\n", summary(phylo_indices$NTI), "\n")
cat("Summary for NRI.SES:\n", summary(phylo_indices$NRI), "\n")

save(phylo_indices, null_indices, file="PhylogeneticDiversityIndicesCalculatedSES.RData")
