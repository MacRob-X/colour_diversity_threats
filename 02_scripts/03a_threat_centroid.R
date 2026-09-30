# Generate matrix of threats and colour space centroid distances
# 2026-09-29
# Robert MacDonald

# Clear environment
rm(list=ls())

# Load custom functions
source(
  here::here(
    "02_scripts", "R", 
    "03_colour_threats_functions.R"
  )
)

# Load libraries ----
library(dplyr)
library(ggplot2)

## EDITABLE CODE ##
# Use latest IUCN assessment data or use most recent assessment data pre- specified cutoff year?
latest <- TRUE
# If not using latest assessment data, specify a cutoff year. Set to NULL if using latest.
cutoff_year <- NULL
# Clade to focus on ("Aves", "Neognaths", "Neoaves", "Passeriformes")
clade <- "Aves"
# Exclude past threats?
# If TRUE, generated threat matrix will exclude threats classified as "Past, Unlikely to Return"
# leaving only ongoing, future, unknown, and past but likely to return threats
exclude_past_threats <-  TRUE


# Load data ----

# Load threat data (Jetz taxonomy version)
if(latest == TRUE){
  jetz_threat_filename <- paste0("jetz_threat_matrix_latest_2026-01-07.csv")
} else if(latest == FALSE){
  jetz_threat_filename <- paste("jetz_threat_matrix", cutoff_year, "cutoff_year.csv", sep = "_")
}
threat_matrix <- read.csv(
  file = here::here(
    "03_output_data", jetz_threat_filename
  )
)

# load colour pattern space (created in Chapter 1 - patch-pipeline)
colspace_path <- paste0("G:/My Drive/patch-pipeline/2_Patches/3_OutputData/", clade, "/2_PCA_ColourPattern_spaces/1_Raw_PCA/", clade, ".matchedsex.patches.250716.PCAcolspaces.rds")
colour_space <- readRDS(colspace_path)[["lab"]][["x"]]


# Data preparation ----

# exclude past threats, if desired
if(exclude_past_threats == TRUE) {
  threat_matrix <- threat_matrix %>% 
    filter(
      timing != "Past, Unlikely to Return"
    )
  past_threat_par <- "exclude_past_threats_"
} else {
  past_threat_par <- NULL
}

# add second-order threat codes to threat matrix (derived from third-order codes)
threat_matrix$second_ord_code <- stringr::str_extract(threat_matrix$code, "[^_]*_[^_]*")

# inspect species with missing threat data
missing_data_spp <- threat_matrix[which(is.na(threat_matrix$second_ord_code) & threat_matrix$notes != "no_threats"),]
# no missing data species
rm(missing_data_spp)

# Remove problem taxon (Cecropis hyperythra) 
# Hirundo_daurica is in there twice, because it corresponds to two BirdLife species
# but is an 'imperfect match' according to AVONET. It's a slightly complicated one also involving
# Hirundo_striolata
# Jetz spp H. striolata and H. daurica correspond to BL spp Cecropis daurica and C. hyperythra
# but in non-straightforward ways
# Essentially C. hyperythra only exists in Sri Lanka so is a subpopulation of the daurica/striolata complex
# The upshot I think is that I just remove the profile for C. hyperythra as in Jetz it is part of 
# daurica/striolata and these are both threatened by the same threat
# Each BL species has different threat profile (Cecropis_hyperythra has no threats, Cecropis_daurica
# is threatened by invasive species)
# ACTION TAKEN : remove the threat profile of C. hyperythra 
# This means I'm using the profile of the nominate one, so is in keeping with my approach to taxonomy matching
# I do this earlier in the pipeline
threat_matrix <- threat_matrix %>% 
  filter(
    species_birdlife != "Cecropis_hyperythra"
  )

# assign second-order IUCN threat types to grouped 'driver of extinction' categories
# same system as Stewart et al 2025 Nat Ecol Evol - grouping is provided in Supplementary Dataset 1
# of that paper
# There are many threats that aren't assigned to one of these groups - this is because these
# threats were non-significant in predicting IUCN threat level in Stewart et al 2025
# We assign these as 'FLAG' in case we want to do anything with them later
threat_matrix <- assign_ext_drivers(threat_matrix)

# For now, let's make all the flagged threats (i.e. those which are not significant predictors
# of extinction risk) have "no_sig_threats", as we might want to use them later
threat_matrix <- threat_matrix |> 
  mutate(
    ex_driver = ifelse(
      ex_driver == "FLAG",
      "no_sig_threats",
      ex_driver
    )
  )

## NOT RUN
# DECISION: let's also make ALL threats for non-threatened (i.e., LC) species NA, as we're not
# interested in threats to LC species
# threat_matrix <- threat_matrix |>
#   mutate(
#     ex_driver = ifelse(
#       iucn_cat == "LC",
#       NA,
#       ex_driver
#     )
#   )

# Pivot wider to get binary extinction driver variables
wide_threat_matrix <- threat_matrix %>% 
  select(
    jetz_species,
    iucn_cat,
    notes,
    ex_driver
  ) %>% 
  distinct() %>% 
  mutate(dummy = 1) %>% 
  tidyr::pivot_wider(
    names_from = ex_driver, values_from = dummy, values_fill = 0
  )

# calculate distance to centroid from colourspace
centr_dists <- dispRity::dispRity(colour_space, metric = dispRity::centroids)$disparity[[1]][[1]]
centr_dists <- data.frame(
  species = sapply(strsplit(rownames(colour_space), split = "-"), "[", 1),
  sex = sapply(strsplit(rownames(colour_space), split = "-"), "[", 2),
  centr_dist = centr_dists
)

# add distance to centroid onto threat matrix
threat_centr <- wide_threat_matrix |> 
  inner_join(centr_dists, by = join_by("jetz_species" == "species"))
# this throws a warning but it's just because we have male and female centroid distance data together
# - it's not a problem

# also add the raw PC values
colour_space <- as.data.frame(colour_space)
colour_space$jetz_species <- sapply(strsplit(rownames(colour_space), split = "-"), "[", 1)
colour_space$sex <- sapply(strsplit(rownames(colour_space), split = "-"), "[", 2)
threat_centr <- threat_centr %>% 
  left_join(
    colour_space,
    by = c("jetz_species", "sex")
  )

# I no longer need the 'notes' column as it only contains info about species with no threats, which is now in
# the 'no_threats' column
threat_centr <- threat_centr %>% 
  select(
    -notes
  )

# reorder columns
threat_centr <- threat_centr %>% 
  relocate(sex, .after = jetz_species) %>% 
  relocate(centr_dist, .after = sex)

# Save as CSV
if(latest == TRUE){
  filename <- paste0(past_threat_par, "centroid_jetz_threat_matrix_latest_2026-01-07.csv")
} else if(latest == FALSE){
  filename <- paste0(past_threat_par, "centroid_jetz_threat_matrix_", cutoff_year, "_cutoff_year.csv")
}
write.csv(
  threat_centr,
  file = here::here(
    "03_output_data", "03a_threat_centroid",
    filename
  ), row.names = FALSE
)
