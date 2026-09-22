# Match IUCN and Jetz taxonomies for threat data

# Clear environment
rm(list=ls())

# Load libraries ----
library(dplyr)

## EDITABLE CODE ##
# Use latest IUCN assessment data or use most recent assessment data pre- specified cutoff year?
latest <- TRUE
# If not using latest assessment data, specify a cutoff year. Set to NULL if using latest.
cutoff_year <- NULL

# Load data ----

# Load threat data
# Set CSV filename and path
if(latest == TRUE){
  threat_filename <- paste0("iucn_threat_matrix_latest_2026-01-07.csv")
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

# load HBW/BirdLife taxonomic checklist v10
bl_10 <- readr::read_csv(
  here::here(
    "01_input_data", "birdlife_taxonomic_checklists", "HBW-BirdLife_Checklist_v10_Oct25",
    "simplified_rxm_digital_checklist_v10.csv"
  ), 
  skip_empty_rows = TRUE
)

# load HBW/BirdLife taxonomic checklist v5.0
bl_5 <- readr::read_csv(
  here::here(
    "01_input_data", "birdlife_taxonomic_checklists", "HBW-BirdLife_Checklist_v5_Dec20",
    "simplified_rxm_digital_checklist_v5.csv"
  ),
  skip_empty_rows = TRUE
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

# Check for species in the threat data which don't appear in the AVONET crosswalk
missing_threat_spp <- threat_species %>% 
  filter(
    !(threat_species %in% avonet_crosswalk$species_birdlife)
  ) %>% 
  pull(threat_species)

# Some of these are likely to be because AVONET uses HBW/BirdLife digital checklist v5.0 and 
# IUCN uses a later version 
# I believe it's version 10 (the latest version as of 2026-09-22, released 2025)

# I can therefore match the names between v5 and v10
# Just need to replace the IUCN 2023 species names (v10) with their v5.0 names, where possible

# Imported code
##################################################################################

# v10 - discard rows referring to subspecies
bl_10 <- bl_10 %>% 
  filter(
    is.na(sub_spp_id) # filter out subspecies
  ) |> 
  select(
    -c(authority, alt_common_names, sub_spp_id)
  ) |> 
  mutate(
    species_10 = sub(" ", "_", scientific_name)
  )


# v5.0 - filter out not recognised complexes
bl_5 <- bl_5 %>% 
  mutate(
    species_5 = sub(" ", "_", scientific_name)
  ) %>% 
  filter(
    iucn_cat_5 != "NR",
    iucn_cat_5 != "UR"
  ) %>% 
  select(
    species_5, sis_rec_id, iucn_cat_5
  )

# check if there are any species which are in the threat data but not the BL v10 checklist
extra_threat_spp <- threat_species %>% 
  filter(
    !(threat_species %in% bl_10$species_10)
  ) %>% 
  pull(threat_species)
# just 1 species in the threat data but not BL v10 - "Pyrrhulagra_portoricensis"

# This is now known as 'Melopyrrha_portoricensis' in the IUCN Red List (checked on iucnredlist.org,
# accessed 2026-09-22)
# It's 'Melopyrrha_portoricensis' in bl v10, so I can just remove 'Pyrrhulagra_portoricensis' row
# in the threat data as the actual threat data is under 'Melopyrrha_portoricensis'
threat_species <- threat_species %>% 
  filter(
    threat_species != "Pyrrhulagra_portoricensis"
  )
threat_matrix <- threat_matrix %>% 
  filter(
    binomial_name != "Pyrrhulagra_portoricensis"
  )

# check if there are any species which are in the avonet crosswalk BL list but not the BL v5.0 checklist
avonet_crosswalk$species_birdlife[which(!(avonet_crosswalk$species_birdlife %in% bl_5$species_5))]
# nope, all good - the versions match

# check if there are any species which are in the v5 checklist but only in the v10 checklist as
# NR species AND are 1_bl_to_1_bt mapping AND do not appear in the threat data - these will not have threat
# data to assign
# Note that most of these will NOT cause a problem since this list will include species 
# which do have threat data but whose scientific name has changed (e.g. Buettikoferella_bivittata/
# Cincloramphus_bivittatus)
# I could use the 2020 category I guess?
problem_species <- avonet_crosswalk[(avonet_crosswalk$species_birdlife %in% bl_5$species_5[!(bl_5$sis_rec_id %in% bl_10$sis_rec_id[bl_10$iucn_cat_10 != "NR"]) & !(bl_5$species_5 %in% threat_species$threat_species)]) & avonet_crosswalk$match_type == "1_bl_to_1_bt", ]
# there aren't any of these
rm(problem_species)


# join the birdlife v10 checklist to the v5.0 checklist (using SISRecID) and 
# filter to only the species which are missing from the AVONET crosswalk data
bl_checklist_match <- bl_5 %>% 
  full_join(
    bl_10, by = "sis_rec_id"
  ) %>% 
  # filter(
  #   species_10 %in% missing_birdlife_species | species_50 %in% missing_birdlife_species
  # ) %>% 
  select(
    -c(seq, order, family_name, family, subfamily, tribe, common_name, scientific_name, synonyms, taxonomic_sources, spc_rec_id)
  )

# check which 5.0 species don't have an 10 name equivalent
bl_checklist_match %>% 
  filter(
    is.na(species_10)
  )

# All 5.0 species have a 10 equivalent

# identify rows which have duplicated v10 species names
dupes_sp10 <- bl_checklist_match[bl_checklist_match$species_10 %in% unique(bl_checklist_match$species_10[duplicated(bl_checklist_match$species_10)]), ]

# remove duplicate bl_10 species which have NA bl_5
duplicated_spp <- unique(bl_checklist_match$species_10[duplicated(bl_checklist_match$species_10)])

for(sp in duplicated_spp){
  dupe_rows <- bl_checklist_match %>% 
    filter(
      species_10 == sp
    )
  # remove row if there's a row with corresponding BL v5 name AND rows with no corresponding BL v5 name
  if(nrow(dupe_rows[!is.na(dupe_rows$species_5), ]) != 0){
    bl_checklist_match <- bl_checklist_match %>% 
      filter(
        !(species_10 == sp & is.na(species_5))
      )
  }
}

# check if any duplicate rows left
dupes_sp10 <- bl_checklist_match[bl_checklist_match$species_10 %in% unique(bl_checklist_match$species_10[duplicated(bl_checklist_match$species_10)]), ]

# can remove the remaining duplicates if their v10 IUCN cat is NR and they have no v5 species name
# (as the NR species won't appear in my IUCN threat data)
duplicated_spp <- unique(bl_checklist_match$species_10[duplicated(bl_checklist_match$species_10)])

for(sp in duplicated_spp){
  dupe_rows <- bl_checklist_match %>% 
    filter(
      species_10 == sp
    )
  # remove rows if v10 IUCN cat is NR
  bl_checklist_match <- bl_checklist_match %>% 
    filter(
      !(species_10 == sp & iucn_cat_10 == "NR")
    )
  
}

# check if any duplicate rows left
dupes_sp10 <- bl_checklist_match[bl_checklist_match$species_10 %in% unique(bl_checklist_match$species_10[duplicated(bl_checklist_match$species_10)]), ]
# no duplicated rows left - all good. Let's remove the variables
rm(dupes_sp10, dupe_rows, duplicated_spp, sp)


# now use this data to replace the 2026 IUCN threat data species names with the v5.0 species names
# (don't replace if the v5.0 species name is NA)
threat_species_v5 <- threat_species %>% 
  left_join(
    bl_checklist_match,
    by = join_by(threat_species == species_10)
  ) %>% 
  rename(
    og_threat_species = threat_species
  ) %>% 
  mutate(
    species_birdlife = ifelse(!is.na(species_5), species_5, og_threat_species)
  ) %>% 
  select(
    -species_5
  )
# "species_birdlife" column now has the v5.0 names - original 2023 IUCN names preserved in
# "species_iucn_2023" column

# add v5 species name column to actual threat matrix
threat_matrix <- threat_species_v5 %>% 
  select(
    og_threat_species,
    species_birdlife
  ) %>% 
  right_join(
    threat_matrix,
    join_by("og_threat_species" == "binomial_name")
  ) %>% 
  rename(
    threat_binomial = og_threat_species
  )

##################################################################################


# Match Jetz species to AVONET(using BirdTree = Jetz)
jetz_avonet <- jetz_species |> 
  left_join(
    avonet_crosswalk,
    by = join_by("jetz_species" == "species_birdtree")
  )

# Now match this to the threat data
jetz_avonet_threat <- jetz_avonet |> 
  left_join(
    threat_species_v5,
    by = "species_birdlife"
  )

# Check what match types exist
match_types <- unique(jetz_avonet_threat$match_type)

# All the species that are 1-1 matched are fine, we don't need to do anything
# there are 6807 of these
one_to_one_matched <- jetz_avonet_threat |> 
  filter(
    match_type == "1_bl_to_1_bt"
  ) |> 
  mutate(
    threat_assign = "direct_assign"
  )

# All the species that are a single BL species to many Jetz species are also fine
# we can just use the threat data for the single BL species for many Jetz species
# there are 168 BT species of these, corresponding to 94 BL species
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
# There are 722 BT species of these, corresponding to 1727 BL species
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
    "species_birdlife"
  )

# Write to CSV
if(latest == TRUE){
  jetz_threat_filename <- paste0("jetz_threat_matrix_latest_2026-01-07.csv")
} else if(latest == FALSE){
  jetz_threat_filename <- paste("jetz_threat_matrix", cutoff_year, "cutoff_year.csv", sep = "_")
}

write.csv(
  final_jetz_threat_data,
  file = here::here(
    "03_output_data", jetz_threat_filename
  ), 
  row.names = F
)





