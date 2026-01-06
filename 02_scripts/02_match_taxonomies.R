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
colspace_path <- "G:/My Drive/patch-pipeline/2_Patches/3_OutputData/Aves/2_PCA_ColourPattern_spaces/1_Raw_PCA/Aves.matchedsex.patches.250716.PCAcolspaces.rds"
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

# Check what match types exist
match_types <- unique(jetz_avonet_threat$match_type)

# All the species that are 1-1 matched are fine, we don't need to do anything
# there are 6806 of these
one_to_one_matched <- jetz_avonet_threat |> 
  filter(
    match_type == "1_bl_to_1_bt"
  ) |> 
  mutate(
    threat_assign = "direct_assign"
  )

# All the species that are a single BL species to many Jetz species are also fine
# we can just use the threat data for the single BL species for many Jetz species
# there are 168 of these
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
# 595 of these, so most of them

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
length(unique(no_nom_subspecies$jetz_species[no_nom_subspecies$nominate == T]))
length(unique(no_nom_subspecies$jetz_species)) - length(unique(no_nom_subspecies$jetz_species[no_nom_subspecies$nominate == T]))
# This directly matches another 94 Jetz species to a nominate subspecies, leaving 33 Jetz species unmatched to nominate subspecies

# Let's split into the species we've already got a nominate subspecies for and those we don't
fixed_no_nom_subspecies <- no_nom_subspecies |> 
  filter(
    nominate == TRUE
  )
still_no_nom_subspecies <- no_nom_subspecies |> 
  filter(
    !(jetz_species %in% unique(no_nom_subspecies$jetz_species[no_nom_subspecies$nominate == T]))
  )

# Some of these will just be a case of an -us changing to an -a or vice versa - we can catch these
# first ones where the Jetz species name ends in -a and the BL species name ends in -us
fixed_a_us_no_nom_subspecies <- still_no_nom_subspecies |> 
  mutate(
    short_jetz_spec = ifelse(grepl("a$", jetz_spec_only), stringr::str_sub(jetz_spec_only, end = -2), jetz_spec_only),
    short_bl_spec = ifelse(grepl("us$", bl_spec_only), stringr::str_sub(bl_spec_only, end = -3), bl_spec_only)
  ) |> 
  filter(
    short_jetz_spec == short_bl_spec
  ) |> 
  mutate(
    nominate = TRUE
  ) |> 
  select(
    -short_jetz_spec,
    -short_bl_spec
  )
# now ones where the Jetz species name ends in -us and the BL species name ends in -a
fixed_us_a_no_nom_subspecies <- still_no_nom_subspecies |> 
  mutate(
    short_jetz_spec = ifelse(grepl("us$", jetz_spec_only), stringr::str_sub(jetz_spec_only, end = -3), jetz_spec_only),
    short_bl_spec = ifelse(grepl("a$", bl_spec_only), stringr::str_sub(bl_spec_only, end = -2), bl_spec_only)
  ) |> 
  filter(
    short_jetz_spec == short_bl_spec
  ) |> 
  mutate(
    nominate = TRUE
  ) |> 
  select(
    -short_jetz_spec,
    -short_bl_spec
  )

# combine into the fixed no-nom subspecies df
fixed_no_nom_subspecies <- fixed_no_nom_subspecies |> 
  bind_rows(
    fixed_a_us_no_nom_subspecies,
    fixed_us_a_no_nom_subspecies
  )

# Now get the ones that STILL don't have a nominate subspecies assigned
still_no_nom_subspecies <- still_no_nom_subspecies |> 
  filter(
    !(jetz_species %in% fixed_no_nom_subspecies$jetz_species)
  )
length(unique(still_no_nom_subspecies$jetz_species))
# there are only 11 species left with no match - I can just manually check these using the taxonomy 
# section on the IUCN website and we're good to go - accessed 2025-12-09

# Jetz: Arses_telescophthalmus
# Nominate BL: Arses_telescopthalmus
still_no_nom_subspecies <- still_no_nom_subspecies |> 
  mutate(
    nominate = ifelse(
      jetz_species == "Arses_telescophthalmus" & species_birdlife == "Arses_telescopthalmus",
      TRUE,
      nominate
    )
  )
# Jetz: Cinclidium_leucurum
# Nominate BL: Myiomela_leucura
still_no_nom_subspecies <- still_no_nom_subspecies |> 
  mutate(
    nominate = ifelse(
      jetz_species == "Cinclidium_leucurum" & species_birdlife == "Myiomela_leucura",
      TRUE,
      nominate
    )
  )
# Jetz: Coracina_tenuirostris
# Nominate BL: Edolisoma_tenuirostre
# still_no_nom_subspecies <- still_no_nom_subspecies |> 
still_no_nom_subspecies <- still_no_nom_subspecies |> 
  mutate(
  nominate = ifelse(
    jetz_species == "Coracina_tenuirostris" & species_birdlife == "Edolisoma_tenuirostre",
    TRUE,
    nominate
  )
)
# Jetz: Monarcha_castaneiventris
# NOTE: this one is an imperfect match - from IUCN website 2025-12-19: "Monarcha castaneiventris 
# and M. erythrostictus (Sibley and Monroe [1990, 1993]) have been lumped and split into 
# M. castaneiventris, M. megarhynchus and M. ugiensis following del Hoyo and Collar (2016)."
# Since we already have M. castaneiventris BT matched to M. castaneiventris and we don't have
# colour information for M. erythrostictus, we can just use the existing match and ignore
# the other BL species (Monarcha_megarhynchus and Monarcha_ugiensis)

# Jetz: Nectarinia_afra
# Nominate BL: Cinnyris_afer
still_no_nom_subspecies <- still_no_nom_subspecies |> 
  mutate(
    nominate = ifelse(
      jetz_species == "Nectarinia_afra" & species_birdlife == "Cinnyris_afer",
      TRUE,
      nominate
    )
)
# Jetz: Phylloscopus_poliocephalus
# NOTE: imperfect match - from IUCN website 2025-12-19: "Phylloscopus poliocephalus and 
# P. makirensis (Sibley and Monroe [1990, 1993]) have been lumped and subsequently split into 
# P. poliocephalus, P. misoriensis and P. maforensis following del Hoyo and Collar (2016)."
# Since we already have P. poliocephalus and P. makirensis BT matched to P. poliocephalus BL, 
# we can ignore the other BL species (Phylloscopus_maforensis and Phylloscopus_misoriensis)

# Jetz: Picus_mentalis
# Nominate BL: Chrysophlegma_mentale
still_no_nom_subspecies <- still_no_nom_subspecies |> 
  mutate(
    nominate = ifelse(
      jetz_species == "Picus_mentalis" & species_birdlife == "Chrysophlegma_mentale",
      TRUE,
      nominate
    )
)
# Jetz: Stachyris_erythroptera
# Nominate BL: Cyanoderma_erythropterum
still_no_nom_subspecies <- still_no_nom_subspecies |> 
  mutate(
    nominate = ifelse(
      jetz_species == "Stachyris_erythroptera" & species_birdlife == "Cyanoderma_erythropterum",
      TRUE,
      nominate
    )
)
# Jetz: Tangara_cyanoptera
# This one is a strange one: Jetz taxonomy contains both Tangara_cyanoptera AND Thraupis_cyanoptera,
# even though these are considered synonyms by the IUCN. Thraupis_cyanoptera (BT) is 1-1 matched
# to Tangara_cyanoptera (BL), but Tangara_cyanoptera (BT) is matched to both Tangara_argentea and 
# Tangara_whitelyi (BL).
# From IUCN 2025-12-19: "Tangara argentea and T. whiteleyi (del Hoyo and Collar 2016) were previously lumped and listed as T. cyanoptera following SACC (2005 & updates); Sibley & Monroe (1990, 1993); Stotz et al. (1996)."
# I will use the illustrations in birdsoftheworld.org to visually inspect which species we have 
# images of
# Ok - it's pretty clear that: 
# Thraupis_cyanoptera (BT) = Tangara_cyanoptera (BL) : blue all over
# Tangara_cyanoptera (BT) = Tangara_argentea (BL) : Black head, yellow body, blue & black wings
# We don't have Tangara_whitelyi (BL) in our image set : Black head, off-white body
# Nominate BL: Tangara_argentea
still_no_nom_subspecies <- still_no_nom_subspecies |> 
  mutate(
    nominate = ifelse(
      jetz_species == "Tangara_cyanoptera" & species_birdlife == "Tangara_argentea",
      TRUE,
      nominate
    )
)
# Jetz: Tephrodornis_gularis
# There is no obvious nominate subspecies for this, as there are two BL matches (Tephrodornis_sylvicola
# and Tephrodornis_virgatus). We therefore follow Stewart et al 2025 and select one of these at random
# (see Stewart et al 2025 Nat Ecol Evol Supplementary Information - Reconciling BirdLife and BirdTree 
# taxonomies and including all BirdLife synonyms)
# Nominate subspecies: Tephrodornis_sylvicola
still_no_nom_subspecies <- still_no_nom_subspecies |> 
  mutate(
    nominate = ifelse(
      jetz_species == "Tephrodornis_gularis" & species_birdlife == "Tephrodornis_sylvicola",
      TRUE,
      nominate
    )
)
# Jetz: Zosterops_palpebrosus
# NOTE: imperfect match - from IUCN website 2025-12-19: "Oriental White-eye Zosterops palpebrosus 
# has been split into Indian White-eye Z. palpebrosus, Hume's White-eye Z. auriventer and Sangkar 
# White-eye Z. melanurus on the basis of thorough morphological comparisons (Wells et al. 2017a, b) 
# and genetic differentiation, morphology and vocalisations (Round et al. 2017, Lim et al. 2019)."
# Since we already have Zosterops_palpebrosus (BT) matched to Zosterops_palpebrosus (BL), we can ignore
# these other BL species (Zosterops_auriventer and Zosterops_melanurus)

# Now add these species into the fixed no-nominate dataset and we're basically done
fixed_no_nom_subspecies <- fixed_no_nom_subspecies |> 
  bind_rows(
    still_no_nom_subspecies[still_no_nom_subspecies$nominate == T, ]
  ) |> 
  select(
   -jetz_spec_only, -bl_spec_only, -nominate 
  )

# Put all the different match types together
nom_subspecies <- nom_subspecies |> 
  select(
    -nominate
  )
final_matched_data <- one_to_one_matched |> 
  bind_rows(
    one_bl_to_many_bt
  ) |> 
  bind_rows(
    nom_subspecies
  ) |> 
  bind_rows(
    fixed_no_nom_subspecies
  ) |> 
  select(
    jetz_species, species_birdlife
  )

# And finally, add the threat data to get the finished, matched threat dataset for the Jetz taxonomy
final_jetz_threat_data <- final_matched_data |> 
  left_join(
    threat_matrix,
    join_by("species_birdlife" == "binomial_name")
  )

# Write to CSV
write.csv(
  final_jetz_threat_data,
  file = here::here(
    "03_output_data", paste("jetz_threat_matrix", cutoff_year, "cutoff_year.csv", sep = "_")
  ), 
  row.names = F
)
