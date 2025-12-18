# Obtain detailed information on threats to each avian species from IUCN Red List
# Based on approach of Stewart et al 2025 Nat Ecol Evol

# Clear environment
rm(list=ls())

# Load functions ----
source(
  here::here(
    "02_scripts", "R", "01_get_threat_info_functions.R"
  )
)

## EDITABLE CODE ##
# Use latest IUCN assessment data or use most recent assessment data pre- specified cutoff year?
latest <- TRUE
# If not using latest assessment data, specify a cutoff year. Set to NULL if using latest.
cutoff_year <- NULL

# Workflow ----

# set iucn API key
api <- "zNXM3miBBknHAGv5mKZraLuhPg2sUTs5c2Ay"

# Get full list of birds in IUCN data
# Takes ~20 minutes but only needs to be run once
# iucn_aves = rredlist::rl_class(class = "aves", key = api)

# Save to RDS
# saveRDS(
#   iucn_aves,
#   file = here::here(
#     "3_output_data", "iucn_data_class_aves.RDS"
#     
#   )
# )

# Read in full Aves IUCN data
iucn_aves <- readRDS(
  here::here(
    "03_output_data", "iucn_data_class_aves.RDS"
  )
)

# Get assessment data
aves_data <- iucn_aves$assessments

# Get unique species names with underscore
iucn_species <- unique(gsub(" ", "_", aves_data$taxon_scientific_name))
# Keep only full species - no subspecies
iucn_species <- iucn_species[which(stringr::str_count(iucn_species, "_") == 1)]


# Set up threat skeleton for species which have no threats
# use a species which has threats, but set them to NA
# let's use the assessment ID for the 2023 assessment of Aburria aburri
threats_skeleton <- make_threat_skeleton(binomial_name = "Aburria_aburri", assessment_year = "2023", api_key = api)


# Set CSV filename and path to write to
if(latest == TRUE){
  filename <- paste0("iucn_threat_matrix_latest_", Sys.Date(), ".csv")
} else if(latest == FALSE){
  filename <- paste("iucn_threat_matrix", cutoff_year, "cutoff_year.csv", sep = "_")
}
writepath <- here::here(
  "03_output_data", filename
)

# Get threat data for all species
# this will write to a new CSV if the specified filename and path doesn't exist yet, or append if
# it already exists
allspec_threat_data <- lapply(iucn_species, get_threat_data, latest = latest, year_cutoff = cutoff_year, csv_writepath = writepath, threats_skeleton = threats_skeleton)

# Convert data to matrix format (if you haven't written to a CSV, which you probably should always do)
threat_matrix <- do.call(rbind, allspec_threat_data)