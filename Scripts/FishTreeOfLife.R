library(readxl)
library(rfishbase)
library(ape)
library(TreeTools)

# Read the species list from a file
specieslist <- read_excel("specieslist.xlsx")

# Validate species names using rfishbase
fish <- rfishbase::validate_names(species_list = specieslist)

# Convert species names to match the format in the tree
tips <- gsub(" ", "_", fish, fixed = TRUE)

# Read the tree from the file
tree <- read.tree(file = "./full.trees") # Data from Rabosky 2018 Fish Tree of Life

# Prune the tree based on the tips
pruned <- lapply(tree, function(tree, tips) drop.tip(tree, setdiff(tree$tip.label, tips)), tips = tips)
class(pruned) <- "multiPhylo"

# Remove duplicate trees
unique_trees <- unique(pruned, incomparables = FALSE, use.edge.length = FALSE, use.tip.label = TRUE)

# Get the number of unique tips
num_unique_tips <- Ntip(unique_trees)

# Get all the tip labels from the pruned trees
pruned_fish <- AllTipLabels(pruned)
pruned_fish_df <- as.data.frame(pruned_fish)

# Save the pruned species list to a file
write.csv(pruned_fish_df, file = "specieslist_new")

# Get the species that don't have a tree
species_without_tree <- setdiff(tips, pruned_fish)
species_without_tree

##==== EXPERT REQUIRED FOR ADDING MISSING SPECIES TO BACKBONE OF TREES =======##

