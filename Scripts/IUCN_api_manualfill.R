library(rredlist)
library(purrr)
df<-traits_sp
df$X<-gsub("_", " ", df$X)
df <- df %>% 
  dplyr::mutate(X = case_when(
    X == "Epinephelus caninus" ~ "Mycteroperca caninus",
    X == "Epinephelus costae" ~ "Mycteroperca costae",
    X == "Epinephelus marginatus" ~ "Mycteroperca marginatus",
    X == "Liza ramada" ~ "Chelon ramada",
    X == "Liza saliens" ~ "Chelon saliens",
    X == "Mullus barbatus barbatus" ~ "Mullus barbatus",
    X == "Auxis rochei rochei" ~ "Auxis rochei",
    X == "Zosterisessor ophiocephalus" ~ "Gobius ophiocephalus",
    TRUE ~ X
  ))

search_species_list <- function(species_list, key) {
  results <- lapply(species_list, function(species) {
    rl_search(name = species, key = key)
  })
  return(results)
}
api_key <- "4139b1d7b6b9ca94388151e78bb55ad629df6d7bd0c9fbc301f710b96a4324c1"
species_names <- df$X
search_results <- search_species_list(species_list = species_names, key = api_key)
result_list <- map(search_results, pluck, "result") 
species_df <- map_df(result_list, ~ as.data.frame(.x))
species_df$X<-species_df$scientific_name

unique_species <- unique(species_df$X)
missing_species <- setdiff(species_names, unique_species)
print(missing_species)

species_df$X<-gsub(" ", "_", species_df$X)
species_df <- species_df %>% 
  dplyr::mutate(X = case_when(
    X == "Mycteroperca_caninus" ~ "Epinephelus_caninus",
    X == "Mycteroperca_costae" ~ "Epinephelus_costae",
    X == "Mycteroperca_marginatus" ~ "Epinephelus_marginatus",
    X == "Chelon_ramada" ~ "Liza_ramada",
    X == "Chelon_saliens" ~ "Liza_saliens",
    X == "Mullus_barbatus" ~ "Mullus_barbatus_barbatus",
    X == "Auxis_rochei" ~ "Auxis_rochei_rochei",
    X == "Gobius_ophiocephalus" ~ "Zosterisessor_ophiocephalus",
    TRUE ~ X
  ))

y<-colnames(mat_PA[,-c(1:2)])
x<-species_df$X
add<-as.data.frame(setdiff(y, x))
add
colnames(add)<-"X"
species_df<-full_join(species_df, add)
species_df <- species_df %>% 
  filter(X %in% y)

row<-rl_search("Aphia minuta", key= api_key, region = "europe")
species_df[179,c(1:30)] <- row$result
row<-rl_search("Diplodus sargus", key= api_key)
species_df[180,c(1:30)] <- row$result
row<-rl_search("Gaidropsarus mediterraneus", key= api_key, region = "europe")
species_df[181,c(1:30)] <- row$result
row<-rl_search("Gobius cobitis", key= api_key, region = "europe")
species_df[182,c(1:30)] <- row$result

species_df$category[species_df$X == "Mullus_barbatus_ponticus"] <- "NE"
species_df$category[species_df$X == "Pegusa_nasuta"] <- "NE"
species_df$category[species_df$X == "Planiliza_haematocheilus"] <- "NE"
species_df$category[species_df$X == "Scophthalmus_maeoticus"] <- "NE"
species_df$population_trend[is.na(species_df$population_trend)] <- "Unknown"
IUCNinfo <- species_df[,-c(18:24,27:30)]

save(species_df, file = "Data/IUCNinfo.RData")
