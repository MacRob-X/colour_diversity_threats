# Match IUCN and Jetz taxonomies for threat data

# Clear environment
rm(list=ls())

# Load libraries ----
library(dplyr)

## EDITABLE CODE ##
# Use latest IUCN assessment data or use most recent assessment data pre- specified cutoff year?
latest <- FALSE
# If not using latest assessment data, specify a cutoff year. Set to NULL if using latest.
cutoff_year <- 2023

# Load data ----

# Load threat data
# Set CSV filename and path
if(latest == TRUE){
  threat_filename <- paste0("iucn_threat_matrix_latest_", Sys.Date(), ".csv")
} else if(latest == FALSE){
  threat_filename <- paste("iucn_threat_matrix", cutoff_year, "cutoff_year.csv", sep = "_")
}
threat_matrix <- read.csv(
  here::here(
    "03_output_data", threat_filename
  )
)

# load colour pattern space
# load colour pattern space (created in Chapter 1 - patch-pipeline)
colspace_path <- "G:/My Drive/patch-pipeline/2_Patches/3_OutputData/Neognaths/2_PCA_ColourPattern_spaces/1_Raw_PCA/Neognaths.matchedsex.patches.250716.PCAcolspaces.rds"
colour_space <- readRDS(
  colspace_path
  )[["lab"]][["x"]]

# load Jetz taxonomic data
taxo_master <- read.csv(
  here::here(
    "01_input_data", "BLIOCPhyloMasterTax_2019_10_28.csv"
  )
)

# load AVONET BirdLife-BirdTree crosswalk 
avonet_crosswalk <- read.csv(
  "G:/My Drive/patch-pipeline/4_SharedInputData/avonet/avonet_v7_birdlife-birdtree-crosswalk.csv"
)


# Workflow ----

# Extract patch species names
patch_species <- unique(sapply(strsplit(rownames(colour_space), split = "-"), "[", 1))
jetz_species <- data.frame(
  jetz_species = patch_species,
  jetz_id = 1:length(patch_species)
)

# Clean AVONET crosswalk names
avonet_crosswalk <- avonet_crosswalk %>% 
  rename(
    species_birdlife = Species1,
    species_birdtree = Species3,
    match_type = Match.type,
    match_notes = Match.notes
  ) %>% 
  mutate(
    species_birdlife = sub(" ", "_", species_birdlife),
    species_birdtree = sub(" ", "_", species_birdtree),
    match_type = snakecase::to_snake_case(match_type)
  ) %>% 
  filter( # remove blank rows
    !if_all(everything(), function(x) x == "")
  )

# we can filter out newly described species as these are all species that don't exist in the 
# birdtree taxonomy
# # There are 4 species in the BT taxonomy that are marked as "invalid taxon" and don't have 
# a BL taxonomy name. These are Anthus_longicaudatus, Hypositta_perdita, Lophura_hatinhensis, 
# and Phyllastrephus_leucolepis. None of these species are in the patch data, so we can remove them
# UPDATE 2025-07-18 - these species are still not present in the 250716 patch data,so we can continue
# to filter them out
# There are also 143 extinct species which have no BT name. Keep these in, as we may get data on 
# extinct species in the future
avonet_crosswalk <- avonet_crosswalk %>% 
  filter(
    match_type != "newly_described_species",
    match_type != "invalid_taxon"
  )

# Get species in threat data
threat_species <- threat_matrix |> 
  distinct(binomial_name, .keep_all = T) |> 
  select(
    binomial_name, assessment_year, notes
  ) |> 
  rename(
    threat_species = binomial_name
  ) |> 
  mutate(
    threat_id = 1:n()
  )



# need to check that all species for which we have data appear in the AVONET BirdTree taxonomy
# identify any species which don't appear
patch_species[which(!(patch_species %in% avonet_crosswalk$species_birdtree))]
# Campylopterus_curvipennis doesn't appear. This is because it has since been renamed
# "Pampa_curvipennis" (as per IUCN Red List website, accessed 04/09/2024)


# Let's change the BirdTree name from "Pampa_curvipennis" to "Campylopterus_curvipennis"
avonet_crosswalk <- avonet_crosswalk %>% 
  mutate(
    species_birdtree = if_else(species_birdtree == "Pampa_curvipennis", "Campylopterus_curvipennis", species_birdtree)
  )


# Match Jetz species to AVONET(using BirdTree = Jetz)
jetz_avonet <- jetz_species |> 
  left_join(
    avonet_crosswalk,
    by = join_by("jetz_species" == "species_birdtree")
  )

# Now match this to the threat data
jetz_avonet_threat <- jetz_avonet |> 
  left_join(
    threat_species,
    by = join_by("species_birdlife" == "threat_species")
  )

# All the species that are 1-1 matched are fine, we don't need to do anything
# there are 6783 of these
one_to_one_matched <- jetz_avonet_threat |> 
  filter(
    match_type == "1_bl_to_1_bt"
  ) |> 
  mutate(
    threat_assign = "direct_assign"
  )

# All the species that are a single BL species to many Jetz species are also fine
# we can just use the threat data for the single BL species for many Jetz species
# there are 167 of these
one_bl_to_many_bt <- jetz_avonet_threat |> 
  filter(
    match_type == "1_bl_to_many_bt"
  ) |> 
  mutate(
    threat_assign = "direct_assign"
  )

# For the species that are many BL species to a single Jetz species, it's a bit more complicated
# I need to try to use the threat data for the single BL species that corresponds to the nominate
# Jetz subspecies, since we mostly have colour data for the nominate subspecies
# There are 720 of these
many_bl_to_one_bt <- jetz_avonet_threat |> 
  filter(
    match_type == "many_bl_to_1_bt"
  ) |> 
  mutate(
    threat_assign = "assign_nominate"
  )

# Check how many on these many_bl_to_1_bt have an obvious nominate subspecies candidate
nom_subspecies <- many_bl_to_one_bt |> 
  filter(
    jetz_species == species_birdlife
  ) |> 
  mutate(
    nominate = T
  )
# 593 of these, so most of them

# Check how many don't have an obvious nominate subspecies candidate
no_nom_subspecies <- many_bl_to_one_bt |> 
  filter(
    jetz_species != species_birdlife
  ) |> 
  filter(
    !(jetz_species %in% nom_subspecies$jetz_species)
  )
# there are 127 Jetz species in this category

# Get just the species part of the binomial and match based on this - species with matching
# species part of binomial are likely to be the nominate subspecies
no_nom_subspecies <- no_nom_subspecies |> 
  mutate(
    jetz_spec_only = sapply(strsplit(jetz_species, split = "_"), "[", 2),
    bl_spec_only = sapply(strsplit(species_birdlife, split = "_"), "[", 2)
  ) |> 
  mutate(
    nominate = ifelse(
      jetz_spec_only == bl_spec_only,
      T,
      F
    )
  )
unique(length(no_nom_subspecies$jetz_species[no_nom_subspecies$nominate == T]))
unique(length(no_nom_subspecies$jetz_species[no_nom_subspecies$nominate == F]))
# This directly matches another 94 Jetz species to a nominate subspecies, leaving 212 Jetz species unmatched to nominate subspecies


# Some of these will just be a case of an -us changing to an -a or vice versa - we can catch these

