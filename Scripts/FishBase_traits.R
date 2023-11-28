#load libraries 
my_libraries <- c("tidyverse", "readxl", "RCurl", "XML", "dplyr", "reshape2", 
                  "writexl", "rfishbase")
sapply(my_libraries, library, character.only = TRUE)

species<-rownames(fishPA)
species<-gsub("_", " ", species)

#rfishbase --> validate names before using function
fish<- validate_names(species)

#Get data using rfishbase package:
bsea_est<-estimate(fish)
#bsea_est <- bsea_est[, colSums(is.na(bsea_est)) == 0]
bsea_spec<-species(fish)
#bsea_spec <- bsea_spec[, colSums(is.na(bsea_spec)) == 0]
bsea_repro<-reproduction(fish)
bsea_eco<-ecology(fish)
bsea_ecos<-morphology(fish)
bsea_sys<-ecosystem(fish)
bsea_diet3<-fooditems(fish)
bsea_extra<-fishbase(fish)
#write_xlsx(bsea_repro,".\\rfishbase_repro.xlsx")
#write_xlsx(bsea_eco,".\\rfishbase_ecology.xlsx")
#write_xlsx(bsea_est,".\\rfishbase_estimate.xlsx")
#write_xlsx(bsea_ecos,".\\rfishbase_morphology.xlsx")
#write_xlsx(bsea_spec,".\\rfishbase_extrafishdata.xlsx")
#write_xlsx(bsea_sys,".\\rfishbase_ecosystem.xlsx")
#write_xlsx(bsea_diet3,".\\rfishbase_diet.xlsx")

#Get FishBase data using manual function
specieslist<-gsub(" ", "-", fish, fixed=TRUE) #change species name so they fit into the function 
specieslist

#FIRST LOAD FUNCTION, then continue 
source("./Scripts/Sources/FishBase.R")
### Example of how to apply the function
Data <- lapply(specieslist,get_fishbase_data)
Data_end <- do.call(rbind,Data)
write.csv(Data_end, file = "Fishbase_extract.csv")

######################### TRAITS_SP WAS MANUALLY MODIFIED ######################
##########################  TO ADD MISSING INFORMATION  ########################
#######################  WHICH WAS NOT PRESENT ON FISHBASE #####################

#We are going to take a few more columns from the rfishbase function, just incase
est_add<-bsea_est[,c(1, 3, 5:6, 8:11, 22)]
spec_add<-bsea_spec[,c(1,38, 58:59)]%>%
  mutate_at(vars(-1), ~ifelse(is.na(.), "unknown", .))
toadd<-full_join(est_add, spec_add) %>% rename(X = Species)
toadd$X[toadd$X == "Auxis rochei"] <- "Auxis rochei rochei"
toadd$X[toadd$X == "Diplodus sargus"] <- "Diplodus sargus sargus"
toadd$X[toadd$X == "Chelon ramada"] <- "Liza ramada"
toadd$X[toadd$X == "Chelon saliens"] <- "Liza saliens"
toadd$X[toadd$X == "Mullus barbatus"] <- "Mullus barbatus barbatus"
toadd$X[toadd$X == "Mullus ponticus"] <- "Mullus barbatus ponticus"
toadd$X[toadd$X == "Oblada melanurus"] <- "Oblada melanura"

#==============================================================================#
#======= FINAL FUNCTIONAL TRAIT MATRIX, ADDING UPDATED INFORMATION ============#
#==============================================================================#


#Get species list from matPA
species<-rownames(fishPA)%>% as.data.frame()
colnames(species)[1] <- "X"
#Change the names back to match traits_sp
species$X[species$X == "Alosa_caspia_caspia"] <- "Alosa_caspia"
species$X[species$X == "Liza_aurata"] <- "Chelon_auratus"
species$X[species$X == "Liza_haematocheila"] <- "Planiliza_haematocheilus"
species$X<-gsub("_", " ", species$X) 

#new dataframe with just the 186 species 
traits<-semi_join(traits_sp, species)

#Double check that all the traits are the same as the original traits you made 
Data_end<-Data_end %>% rownames_to_column("X")
Data_end$X<-gsub("-", " ", Data_end$X)
Data_end$X[Data_end$X == "Auxis rochei"] <- "Auxis rochei rochei"
Data_end$X[Data_end$X == "Diplodus sargus"] <- "Diplodus sargus sargus"
Data_end$X[Data_end$X == "Chelon ramada"] <- "Liza ramada"
Data_end$X[Data_end$X == "Chelon saliens"] <- "Liza saliens"
Data_end$X[Data_end$X == "Mullus barbatus"] <- "Mullus barbatus barbatus"
Data_end$X[Data_end$X == "Mullus ponticus"] <- "Mullus barbatus ponticus"
Data_end$X[Data_end$X == "Oblada melanurus"] <- "Oblada melanura"
Data_end2<-Data_end[,-c(6:8,10:15,20:21, 24:25, 27:34 )]
traits<-traits[,-c(1:8, 12:14, 16)]
pray<-full_join(traits, Data_end2, by="X", suffix=c("", "cry"))
pray <- pray[, order(colnames(pray))]
#everything looks in order, other than the changes that were manually made 

#========================BLACK SEA CHECKLIST =================================#
#Compare and contrast list of fish species in the Black Sea from the the
# http://www.blacksea-commission.org/
#NOTE: they also reffer to FishBase for fish in the Black Sea

#We are most interested in the origins of these fish which they have prepared in 
# Black Sea Check List
BlackSeaCheckList<-read_excel("Data/BlackSeaCheckList.xlsx")
#Origins were found manually for  "checklist_add"
checklist_add<-read_excel("Data/BlackSeaCheckList.xlsx", 
                          sheet = "missing from BlackSeaCheckList")

BlackSeaCheckList$X[BlackSeaCheckList$X == "Alburnus alburnus"] <- "Alburnus chalcoides"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Alosa caspia paleostomi"] <- "Alosa tanaica"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Atherina bonapartii"] <- "Atherina boyeri"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Belone belone euxini"] <- "Belone belone"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Chelidonichthys lucernus"] <- "Chelidonichthys lucerna"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Engraulis encrasicolus maeoticus"] <- "Engraulis encrasicolus"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Gymnammodytes cicerellus"] <- "Gymnammodytes cicerelus"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Gymnocephalus acerina"] <- "Gymnocephalus cernua"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Lipophrys adriaticus"] <- "Labrus mixtus"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Liza haematocheila"] <- "Planiliza haematocheilus"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Liza ramada"] <- "Chelon ramada"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Liza saliens"] <- "Chelon saliens"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Liza aurata"] <- "Chelon auratus"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Percarina demidoffi"] <- "Percarina demidoffii"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Neogobius eurycephalus"] <- "Ponticola eurycephalus"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Neogobius kessleri"] <- "Ponticola kessleri"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Neogobius platyrostris"] <- "Ponticola platyrostris"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Neogobius ratan"] <- "Ponticola ratan"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Neogobius syrman"] <- "Ponticola syrman"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Trachurus mediterraneus ponticus"] <- "Trachurus mediterraneus"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Tripterygion tripteronotus"] <- "Tripterygion tripteronotum"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Apletodon dentatus"] <- "Apletodon bacescui"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Oblada melanura"] <- "Oblada melanurus"
BlackSeaCheckList$X[BlackSeaCheckList$X == "Parablenniuis sanguinolentus"] <- "Parablennius sanguinolentus" 

check<-bind_rows(BlackSeaCheckList, checklist_add)
check$X[check$X == "Auxis rochei"] <- "Auxis rochei rochei"
check$X[check$X == "Diplodus sargus"] <- "Diplodus sargus sargus"
check$X[check$X == "Chelon ramada"] <- "Liza ramada"
check$X[check$X == "Chelon saliens"] <- "Liza saliens"
check$X[check$X == "Mullus barbatus"] <- "Mullus barbatus barbatus"
check$X[check$X == "Mullus ponticus"] <- "Mullus barbatus ponticus"
check$X[check$X == "Oblada melanurus"] <- "Oblada melanura"
Check<-semi_join(check, traits, by = "X")

(y<-setdiff(Check$X, toadd$X))

toadd<-full_join(Check, toadd)

#===================== GET UPDATED IUCN STATUSES ===============================#

load("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis/Data/IUCNinfo.RData")
IUCNadd <- species_df[,c(13,15,31)]%>% relocate(X, .before = 1)
IUCNadd$X<-gsub("_", " ", IUCNadd$X)
(y<-setdiff(IUCNadd$X, toadd$X))

toadd<-full_join(toadd, IUCNadd)
traits<-traits[,-c(1, 6:7, 12:14, 16)]%>% relocate(X, .before = common)

(y<-setdiff(traits$X, toadd$X))

traits_186<-full_join(traits, toadd)
traits_186<-traits_186 %>% relocate(category, .before = resilience)
traits_186$stat<-NULL
traits_186<-traits_186 %>% relocate(Common, .before = Depth_min)
traits_186$common <- ifelse(is.na(traits_186$common), traits_186$Common, traits_186$common)
traits_186$Common<-NULL


species<-rownames(fishPA)%>% as.data.frame()
colnames(species)[1] <- "X"
traits_186$X<-gsub(" ", "_", traits_186$X, fixed = T ) 

(y<-setdiff(traits_186$X, species$X))

traits_186 <- traits_186 %>%
  mutate(X = case_when(
    X == "Alosa_caspia" ~ "Alosa_caspia_caspia",
    X == "Chelon_auratus" ~ "Liza_aurata",
    X == "Planiliza_haematocheilus" ~ "Liza_haematocheila",
    TRUE ~ X
  ))
(y<-setdiff(traits_186$X, species$X))

library(janitor)
traits_186 <-clean_names(traits_186)

traits_186$depth_ave <- (traits_186$depth_min + traits_186$depth_max) / 2 
traits_186 <- traits_186 %>% relocate(depth_ave, .before = 9)
traits_186$depth_range <- traits_186$depth_max - traits_186$depth_min
traits_186 <- traits_186 %>% relocate(depth_range, .before = 10)
traits_186 <- traits_186 %>% rename(trend=population_trend)  %>% rename(position=env_3) %>%  rename(shape=shape1)
traits_186$resilience <- tolower(traits_186$resilience)
traits_186$environment <- tolower(traits_186$environment)
traits_186$origin <- tolower(traits_186$origin)
traits_186$trend <- tolower(traits_186$trend)
colnames(traits_186)[1:6] <- sapply(colnames(traits_186)[1:6], function(x) paste(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)), sep = ""))


load("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis/Data/IUCNinfo.RData")
IUCNinfo <- species_df[,-c(1:5,7:30)]
IUCNinfo <- IUCNinfo %>%
  mutate(X = case_when(
    X == "Alosa_caspia" ~ "Alosa_caspia_caspia",
    X == "Chelon_auratus" ~ "Liza_aurata",
    X == "Planiliza_haematocheilus" ~ "Liza_haematocheila",
    TRUE ~ X
  ))
addOrder<-full_join(IUCNinfo, traits_186)
addOrder$order[is.na(addOrder$order)] <- addOrder$Order[is.na(addOrder$order)]
addOrder$order <- toupper(addOrder$order)
addOrder <- addOrder %>%
  rename(morph = order) %>%
  relocate(morph, .after = shape)
traits_186<-addOrder
write.csv(traits_186, file="./Data/Traits186.csv")
