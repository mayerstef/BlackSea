library(readxl)
library(dplyr)
library(tidyr)

codes<-read.csv("Data/speciescode.csv", row.names = 1,  header= TRUE)
SAU <- read.csv("Data/SAU EEZ 794,966,804,647,100,268,642 v50-0.csv")
SAU <- SAU %>% 
  group_by(area_name, scientific_name, common_name, functional_group, commercial_group) %>% 
  summarise(Tot=sum(tonnes)) %>%
  filter(commercial_group != "Crustaceans", 
         commercial_group != "Molluscs",
         functional_group != "Cephalopods",
         functional_group != "Jellyfish") %>% 
  mutate(area_name = ifelse(area_name == "Russia (Black Sea)", "Russia", 
                            ifelse(area_name == "Turkey (Black Sea)", "Turkey",
                                   ifelse(area_name == "Turkey (Marmara Sea)", "Turkey_marmara", area_name))))

fishstat<-read.csv("Data/BlackSeaCountriesFishStatJRegionalCaptures.csv") 
colnames(fishstat)[1:4] <- c("country", "common_name", "scientific_name", "area")
fishstat <- fishstat %>%
  mutate(country = ifelse(country == "T\xfcrkiye", "Turkey", 
                          ifelse(country == "Russian Federation", "Russia",
                                 ifelse(country == "Un. Sov. Soc. Rep.", "USSR", country)))) %>%
  mutate(area = case_when(
    area == "Tunas (GFCM area)" ~ "",
    area == "Marmara Sea" ~ "_marmara",
    area == "Black Sea" ~ "",
    area == "Azov Sea" ~ "_azov")) %>% 
  filter(!(country %in% "Totals - Tonnes - live weight")) %>% 
  filter(!(country %in% "FAO-GFCM. 2022. Fishery and Aquaculture Statistics. GFCM capture production 1970-2020 (FishStatJ). In: FAO Fisheries and Aquaculture Division [online]. Rome. Updated 2022. www.fao.org/fishery/statistics/software/fishstatj/en"))

non_num<-sapply(fishstat, function(x) !is.numeric(x))
df<-fishstat[, !non_num] %>% 
  transmute(Tot = rowSums(.))
fishstat <- bind_cols(fishstat[,1:4], df) 
fishstat$area_name <- paste(fishstat$country, fishstat$area, sep = "")
fishstat$common_name <- gsub("\\(.*?\\)", "", fishstat$common_name)
fishstat<- fishstat %>% select(area_name, scientific_name, common_name, Tot)
fisheries<-full_join(SAU, fishstat, by=c("scientific_name", "area_name"),  suffix = c("SAU", "FISHSTAT")) # %>% 
fisheries <- fisheries %>%
  group_by(scientific_name) %>%
  fill(functional_group, commercial_group, common_nameSAU, common_nameFISHSTAT) %>%
  ungroup() %>%
  mutate(common_nameSAU = coalesce(common_nameSAU, common_nameFISHSTAT)) %>% 
  rename(X = scientific_name)

# Edit Discrepancies 
fisheries$X <- sub("Gadiformes", "Gaidropsarus mediterraneus", fisheries$X)
fisheries$X <- sub("Gaidropsarus spp", "Gaidropsarus mediterraneus", fisheries$X)
fisheries$X <- sub("Percarina demidoffi", "Percarina_demidoffii", fisheries$X)
fisheries$X <- sub("Mugil soiuy", "Planiliza_haematocheilus", fisheries$X)
fisheries$X <- sub("Liza haematocheilus", "Planiliza_haematocheilus", fisheries$X)
# some other edits may have been made (manually) but were not included in the
# final community matrix  

fisheries$X <- gsub(" ", "_", fisheries$X, fixed=TRUE)
fisheries<-inner_join(fisheries, traits_sp, by = "X")

fisheries<-fisheries %>% 
  rename(EEZ = area_name) %>% 
  rename(Species = X)
fisheries <- fisheries[!is.na(fisheries$EEZ),]
species_depth <- fisheries %>%
  group_by(Species) %>%
  summarise(Depth_min = unique(Depth_min),
            Depth_max = unique(Depth_max))

writexl::write_xlsx(fisheries, "Data/fisheries.xlsx")  
