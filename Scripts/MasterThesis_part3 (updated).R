################################################################################
####                                                                        ####
####         Master thesis:LINKING PATTERNS IN PHYLOGENY, TRAITS,           ####
####            ABIOTIC VARIABlES AND SPACE: BLACK SEA FISHES               ####
####                                                                        ####
####                               PART 3                                   ####
####   Stefanie Mayer                                                       ####
####                                                                        ####
####        calbouy@ethz.ch                                                 ####
################################################################################
#================== Workspace preparation & Data loading ======================#
### Library loading
# Data Manipulation and Analysis
lib_general <- c("dplyr", "tidyr","reshape2","readxl", "ggplot2", "Matrix", "stringr", "ggpubr")

# Spatial Data and Geospatial Analysis
lib_gis <- c("raster", "rgeos", "rgdal", "sf", "sp", "shape", "geometry","geosphere", 
             "lwgeom", "fasterize", "terra","rnaturalearthdata", "rnaturalearth", 
             "shape", "ggspatial", "prettymapr")

# Phylogenetic Diversity Analysis
lib_phylo <- c("ape", "PhyloMeasures", "tidytree", "TreeTools", "phyloregion", "picante")

#Functional Diversity Analysis
lib_fun <- c("mFD", "hillR", "readr", "magrittr", "rstatix","FD", "vegan")



#Choose which libraries you would like to load
sapply(lib_general,library,character.only=TRUE)
sapply(lib_gis,library,character.only=TRUE)
sapply(lib_phylo,library,character.only=TRUE)
sapply(lib_fun,library,character.only=TRUE)

#Load part 1 data
load("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis/Data/part1.RData")
rm(t, goodtrees, matPA, spdf)
#Load Part 2 Data
#phylo_indices<-
load("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis/Data/PhylogeneticDiversityIndicesCalculatedSES.RData")
#functional_indices<-
load("C:/Users/stefa/OneDrive - UGent/IMBRSea/Thesis/FINAL/MAYER_MasterThesis/Data/FinalFunctionalIndicesSES.RData")
fun_indices<-tibble::rownames_to_column(fun_indices, var = "grids")
phylo_indices<-tibble::rownames_to_column(phylo_indices, var = "grids")
BSrows<-BSrows %>% dplyr::rename(Bioregions = species)
bio_indices<-inner_join(fun_indices, phylo_indices)
bio_indices<-inner_join(bio_indices, BSrows)
bio_indices<- bio_indices %>% 
  dplyr::relocate(Bioregions, .before = sp_richn)
bio_indices$Countries <- gsub("_.*", "", bio_indices$Bioregions)
#bio_indices<-bio_indices %>% tibble::column_to_rownames("grids")

#Results from phyloregion package + mfd package
SES_PD <- read_csv("Data/SES_PD.csv", col_types = cols(...1 = col_skip(), 
                                                  Bulgaria = col_skip(), Georgia = col_skip(), 
                                                  Romania = col_skip(), Russia = col_skip(), 
                                                  Russia_azov = col_skip(), Turkey_bs1 = col_skip(), 
                                                  Turkey_bs2 = col_skip(), Turkey_marmara = col_skip(), 
                                                  Ukraine = col_skip(),Ukraine_azov = col_skip())) 
bio_indices<-inner_join(bio_indices, SES_PD, by="grids")

#save(bio_indices, file="bio_indices.RData")

result <- bio_indices %>%
  group_by(Bioregions) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE))
result_min <- bio_indices %>%
  group_by(Bioregions) %>%
  summarize(across(where(is.numeric), min, na.rm = TRUE))
result_max <- bio_indices %>%
  group_by(Bioregions) %>%
  summarize(across(where(is.numeric), max, na.rm = TRUE))

#==============================================================================#
spc_tbl <- read_csv("Data/Traits186.csv", 
                    col_types = cols(...1 = col_skip()))

fishPA<-as.matrix(matPA)
fishPA <- fishPA[rowSums(fishPA) >= 6,]

mat_trait<-as.data.frame(bio_indices$grids)
colnames(mat_trait)<-"grids"
row.names(mat_trait)<-mat_trait$grids
mat_trait<-inner_join(mat_trait, BSrows)

# Create a list of variables that you want to calculate mean for each site
vars <- c("depth_ave", "depth_range", "troph", "max", "max_length_tl", "troph_2", "se_troph", 
          "a", "sd_log10a", "b", "sd_b", "max_length_sl", "vulnerability")

# Initialize a list to store all the new trait matrices
trait_continuous<-mat_trait
trait_matrices <- list()
for (var in vars) {
  var_vector <- spc_tbl[[var]]
  names(var_vector) <- spc_tbl$X
  var_matrix <- fishPA
  for (i in 1:ncol(fishPA)) {
    var_matrix[,i] <- var_vector[colnames(fishPA)[i]] * fishPA[,i]
  }
  trait_matrices[[var]] <- var_matrix
  trait_continuous[[paste0(var, "_mean")]] <- rowMeans(var_matrix, na.rm = TRUE)
}
# Replace 'var_name' with the name of the variable you want to plot
var_name <- "vulnerability"  # example
selected_df <- trait_continuous %>% dplyr::select(Bioregions, paste0(var_name, "_mean"))
colnames(selected_df)[2] <- "Mean"

# Plot
ggplot(selected_df, aes(x = Bioregions, y = Mean)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Bioregion", y = "Mean") +
  ggtitle(var_name)


# Create a list of categorical variables
cat_vars <- c("position", "shape", "repro", "fert", "parent", "category", 
              "resilience", "origin", "importance", "price_categ", "trend")

trait_cat<-mat_trait
# Loop over all categorical variables
for (cat_var in cat_vars) {
  levels <- unique(spc_tbl[[cat_var]])
  for (level in levels) {
    level_vector <- ifelse(spc_tbl[[cat_var]] == level, 1, 0)
    names(level_vector) <- spc_tbl$X
    level_matrix <- fishPA
    for (i in 1:ncol(fishPA)) {
      level_matrix[,i] <- level_vector[colnames(fishPA)[i]] * fishPA[,i]
    }
    trait_cat[[paste0(cat_var, "_", level, "_prop")]] <- rowSums(level_matrix, na.rm = TRUE) / rowSums(fishPA)
  }
}

mat_trait3<-as.data.frame(bio_indices$grids)
colnames(mat_trait3)<-"grids"
row.names(mat_trait3)<-mat_trait3$grids
BSrows[[2]]<-"Bioregions"
mat_trait3<-inner_join(mat_trait3, BSrows)

# Loop over all categorical variables
for (cat_var in cat_vars) {
  levels <- unique(spc_tbl[[cat_var]])
  for (level in levels) {
    level_vector <- ifelse(spc_tbl[[cat_var]] == level, 1, 0)
    names(level_vector) <- spc_tbl$X
    level_matrix <- fishPA
    for (i in 1:ncol(fishPA)) {
      level_matrix[,i] <- level_vector[colnames(fishPA)[i]] * fishPA[,i]
    }
    mat_trait3[[paste0(cat_var, "_", level, "_mean")]] <- rowMeans(level_matrix, na.rm = TRUE)
  }
}

# Select the trait
trait <- "position"  # Replace this with any category from cat_vars
selected_cols <- grep(paste0("^", trait, "_"), names(mat_trait3), value = TRUE)
selected_cols <- c(selected_cols, "Bioregions")
mean_df <- mat_trait3[selected_cols] %>% aggregate(. ~ Bioregions, ., mean)
se_df <- mat_trait3[selected_cols] %>% aggregate(. ~ Bioregions, ., function(x) sd(x) / sqrt(length(x)))
mean_long_df <- mean_df %>% pivot_longer(-Bioregions, names_to = "Trait", values_to = "Mean")
se_long_df <- se_df %>% pivot_longer(-Bioregions, names_to = "Trait", values_to = "SE")
long_df <- full_join(mean_long_df, se_long_df, by = c("Bioregions", "Trait"))
long_df$Trait <- sub(paste0("^", trait, "_"), "", long_df$Trait)
long_df$Trait <- sub("_mean", "", long_df$Trait)


# Plot
ggplot(long_df, aes(x = Trait, y = Mean, fill = Bioregions)) +
  geom_bar(stat = "identity", position = "stack") +
  #geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Bioregion", y = "Mean Proportion", fill = paste0(trait, " Type")) +
  scale_fill_viridis_d()



# Create a new list to hold the binary matrices
binary_mat_trait <- list()
for (cat_var in cat_vars) {
  levels <- unique(spc_tbl[[cat_var]])
  for (level in levels) {
    level_vector <- ifelse(spc_tbl[[cat_var]] == level, 1, 0)
    names(level_vector) <- spc_tbl$X
    level_matrix <- fishPA
    for (i in 1:ncol(fishPA)) {
      level_matrix[,i] <- level_vector[colnames(fishPA)[i]] * fishPA[,i]
    }
    binary_mat_trait[[paste0(cat_var, "_", level, "_bin")]] <- level_matrix
  }
}

#==============================================================================#
mat_trait$Countries<-bio_indices$Countries
mat_trait2$Countries<-bio_indices$Countries
mat_trait3$Countries<-bio_indices$Countries

long_data <- mat_trait %>%
  tidyr::pivot_longer(cols = c("depth_ave_mean", "depth_range_mean", "troph_mean", "max_mean", "vulnerability_mean"),
               names_to = "variable",
               values_to = "value")

ggplot(long_data, aes(x = value, y = Bioregions, fill=Countries)) +
  geom_boxplot() +
  scale_fill_brewer(palette=3)+
  facet_wrap(~variable, scales = "free_x") +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") +
  ylab("Bioregions") 


#==============================================================================#
install.packages("FSA")
library(car)
library(lmtest)
library(gvlma)
library(FSA)
model.data<-bio_indices
model.lm<-lm(sp_richn~factor(Bioregions),data=model.data)    #Model
model.data$res<-residuals(model.lm)		      #calculate residuals
model.data$res2<-model.data$res*model.data$res		      #calculate squared residuals
out<-lm(res2~factor(Bioregions),data=model.data)	      #model squared residuals again predictors
summary(out)			      #Overall F test
nobs(out)*summary(out)$r.squared		      #LM statistic
bptest(model.lm)                       #obtain LM statistic direct from model output	
nobs(out)*summary(out)$r.squared
gvlma(model.lm) 

# Check normality of residuals
qqnorm(resid(model.lm))
qqline(resid(model.lm))
shapiro.test(resid(model.lm)) #data too large 
car::leveneTest(resid(model.lm) ~ Bioregions, data =bio_indices)
par(mfrow = c(2, 2))
plot(model)

## === === === === === === use non parametric methods === === === === === === ##

#nonparametric tests 
#Comparing Species richness among Bioregions using a non-parametric test:
richness_kruskal <- kruskal.test(sp_richn ~ Bioregions, data = bio_indices)
richness_kruskal

#Comparing Phylogenetic diversity among Bioregions using a non-parametric test:
PD_kruskal <- kruskal.test(PD ~ Bioregions, data = bio_indices)
PD_kruskal

#Comparing Functional diversity (let's say, fric as an example) among Bioregions using a non-parametric test:
fric_kruskal <- kruskal.test(fric ~ Bioregions, data = bio_indices)
fric_kruskal

# Perform Dunn's test for post hoc pairwise comparisons
posthoc_dunn <- dunnTest(fric ~ factor(Bioregions), data = bio_indices, method = "bonferroni")
print(posthoc_dunn)
"red" #put this in APPEMDIX

# Create a box plot
box_plot <- ggplot(bio_indices, aes(x = Countries, y = sp_richn)) +
  geom_boxplot() +
  labs(title = "Box Plot of Species Richness by Bioregions")
box_plot

# Create a strip chart
strip_chart <- ggplot(bio_indices, aes(x = Bioregions, y = sp_richn)) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +
  labs(title = "Strip Chart of Response Variable by Group")

# Print the plots
print(box_plot)
print(strip_chart)



###############################################################################

duplicated_rows <- duplicated(BSrows$grids)
BSrows<- BSrows[!duplicated_rows, ]
df<-inner_join(SES_PD, BSrows)

df2<-df[,-c(1, 14:16)]
test <- aggregate(. ~ species, data = df2, FUN = mean)

df<-left_join(BSrows, SES_PD)
df2<-df[,-c(1, 15:17)]
test2 <- aggregate(. ~ species, data = df2, FUN = mean)

freq_table <- table(BSrows$species)
freq_table


###############################################################################

save(mat_trait, mat_trait2, mat_trait3, bio_indices, file="part3.RData")










#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4


see <- bio_indices%>%
  group_by(Bioregions) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)
see[,1:2]


newgrids<-bio_indices[,1:2]
df<-as.data.frame(fishPA) %>% tibble::rownames_to_column("grids")
df1<-inner_join(newgrids, df)
df <- df1 %>%
  group_by(Bioregions) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE)
df[,2:187] <- ifelse(df[,2:187] != 0, 1, 0)
df$sum <- rowSums(df[, 2:187])
df <- df[, c(1, 188, 2:187)]
df[,1:2] 
df2<- df1 %>%
  group_by(Bioregions) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE)
df2$sum <- rowSums(df2[, 2:187])
df2 <- df2[, c(1, 188, 2:187)]
df2[,1:2] 


df3<-as.data.frame(t(df)) %>% tibble::rownames_to_column("X") 
colnames(df3)<-df3[1,]
df3<-df3[-1,]
df3$X<-df3$Bioregions
df4<-full_join(df3, spc_tbl)
df4<-df4[-1,-1]

# Melt the DataFrame
df_long <- melt(df4, id.vars=c("X", "Class", "Order", "Fam", "Genus", "Common", "depth_min", "depth_max", "depth_ave", "depth_range", "position", "shape", "morph", "troph", "max", "maxmeas", "repro", "fert", "parent", "category", "resilience", "environment", "origin", "max_length_tl", "troph_2", "se_troph", "a", "sd_log10a", "b", "sd_b", "max_length_sl", "vulnerability", "importance", "price_categ", "trend"),
                measure.vars=c("Bulgaria", "Georgia", "Romania", "Russia", "Russia_azov", "Turkey_bs1", "Turkey_bs2", "Turkey_marmara", "Ukraine", "Ukraine_azov"),
                variable.name="Bioregions",
                value.name="Value")

# Group the data by taxa and regions and calculate species richness
species_richness <- df_long %>%
  group_by(Class, Order, Fam, Genus, Bioregions) %>%
  summarize(SpeciesRichness = n_distinct(Common))

# View the resulting species richness table
print(species_richness)

df1 <- df1 %>%
  mutate(sum = rowSums(across(where(is.numeric))))%>%
  separate(grids, into = c("X", "Y"), sep = "_")

dfs <- split(df1, df1$Bioregions)
BU<-dfs$Bulgaria
GE<-dfs$Georgia
RO<-dfs$Romania
RU<-dfs$Russia
RUa<-dfs$Russia_azov
UA<-dfs$Ukraine
UAa<-dfs$Ukraine_azov
TR1<-dfs$Turkey_bs1
TR2<-dfs$Turkey_bs2
TRm<-dfs$Turkey_marmara

test<-t(RUa)
colnames(test) <- paste(test[1, ], test[2, ], sep = "_")
test <- test[-c(1, 2,3), ]
test<-as.data.frame(test)
test <- test %>%
  mutate_all(as.numeric)
test$sum<-rowSums(test)
ughAzov1<-test
ughAzov2<-test
rm(test)


# Select the top 20 rows from each data frame based on the 'sum' column
top20_ughAzov1 <- head(ughAzov1[order(ughAzov1$sum, decreasing = TRUE), ], 20)
top20_ughAzov2 <- head(ughAzov2[order(ughAzov2$sum, decreasing = TRUE), ], 20)

# Create a new data frame 'ughAzov' with columns 'sum1', 'sum2', 'colname1', 'colname2'
ughAzov <- data.frame(sum1 = top20_ughAzov1$sum,
                      sum2 = top20_ughAzov2$sum,
                      rowname1 = row.names(top20_ughAzov1),
                      rowname2 = row.names(top20_ughAzov2))
print(ughAzov)
both<-intersect(ughAzov$rowname1, ughAzov$rowname2)
View(both)

df1.2 <- df1 %>%
  mutate(Bioregions = case_when(
    Bioregions == "Russia_azov" | Bioregions == "Ukraine_azov" ~ "azov",
    TRUE ~ Bioregions
  ))

dfs <- split(df1.2, df1.2$Bioregions)
Azov<-dfs$azov
Azov <- Azov %>%
  mutate(row_names = paste(X, Y, sep = "_")) 
row.names(Azov) <- Azov$row_names
Azov <- Azov[, -which(names(Azov) == "row_names")]
Azov<-Azov[,-c(1,2,3)]
Azov$sum<-NULL
pca <- rda(Azov)
print(pca)
simper_result <- simper(Azov)
print(simper_result)



######################################################################################



count_list <- list()
for (cat_var in cat_vars) {
  # Get the unique categories for this variable
  levels <- unique(spc_tbl[[cat_var]])
  for (level in levels) {
    # Create a vector marking species in the current category with 1, others with 0
    level_vector <- ifelse(spc_tbl[[cat_var]] == level, 1, 0)
    names(level_vector) <- spc_tbl$X
    # Create a copy of the community matrix
    level_matrix <- fishPA
    # Multiply the community matrix by the category vector. This zeros out species not in the category.
    for (i in 1:ncol(fishPA)) {
      level_matrix[,i] <- level_vector[colnames(fishPA)[i]] * fishPA[,i]
    }
    # calculate the sum (count) of each row (site)
    count_list[[paste0(cat_var, "_", level, "_count")]] <- rowSums(level_matrix, na.rm = TRUE)
  }
}


count_df <- do.call(cbind, count_list)
final_df <- cbind(mat_trait, count_df)

Num_Sites <- table(mat_trait$Bioregion)
Num_Sites_df <- data.frame(Bioregions = names(Num_Sites), Num_Sites = as.integer(Num_Sites))
final_df <- merge(final_df, Num_Sites_df, by = "Bioregions")
for (count_col in grep("_count$", names(final_df), value = TRUE)) {
  final_df[[paste0(count_col, "_per_site")]] <- final_df[[count_col]] / final_df$Num_Sites
}

library(tidyverse)
final_long_df <- final_df %>% 
  pivot_longer(cols = contains("_per_site"),
               names_to = "category_level",
               values_to = "count") %>% 
  separate(category_level, into = c("category", "level"), sep = "_", extra = "drop") 

final_long_df <- final_long_df %>%
  filter(category != "price") %>% 
  filter(category != "")

test<-aggregate(final_long_df, count~Bioregions+category+level, FUN=mean)

library(ggpubr)
ggplot(test, aes(x = level, y = count, fill = Bioregions)) +
  scale_fill_viridis_d(option="G")+
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ category, scales = "free_x", dir="v") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Level") +
  ylab("Count") +
  labs(fill = "Bioregions") +
  theme_pubr(legend = "right")

ggplot(final_long_df, aes(x = level, y = count, fill = Bioregions)) +
  geom_bar(stat = "identity") +
  facet_grid(~ category, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Level") +
  ylab("Count") +
  labs(fill = "Bioregions") +
  theme_bw()

# Calculate the proportion of each category's level for each bioregion
final_long_df2 <- final_long_df %>%
  group_by(category, level) %>%
  mutate(prop = count / sum(count))

final_long_df <- final_long_df %>%
  filter(category != "price") %>% 
  filter(category != "resilience") %>% 
  filter(category != "trend") %>% 
  filter(category != "importance") %>% 
  filter(category != "origin")



# Plot
ggplot(final_long_df2, aes(x = prop, y = level, fill = Bioregions)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(option="G")+
  facet_wrap(~ category, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Level") +
  ylab("Proportion") +
  labs(fill = "Bioregions") +
  scale_x_continuous(labels = scales::percent) +
  theme_pubr(legend="right")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4

# Plotting
ggplot(df_long2, aes(x = interaction(cat, variable), y = y_values, fill = bioregion)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(option = "G", direction = -1) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Categories", y = "Proportion") +
  ggtitle("Proportions per Category") +
  facet_wrap(~ variable, scales = "free_x", dir = "v") +
  theme_classic2() +
  scale_x_discrete(labels = function(x) sub("\\..*", "", x))




ggplot(df_long2, aes(x = interaction(cat, variable), y = y_values, fill = bioregion)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(option = "G", direction = 1, name = "Bioregions") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "", y = "Proportion (%)",
       title = "") +
  facet_wrap(~ variable, scales = "free_x", dir = "v") +
  theme_classic2(base_family = "Helvetica")+
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text=element_text(size=12),
        plot.title = element_text(size=16,face="bold"),
        legend.title = element_text(size=14),
        strip.text = element_text(size = 12)) +
  scale_x_discrete(labels = function(x) sub("\\..*", "", x))



