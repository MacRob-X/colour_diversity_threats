# Effects of abatement of threat groups on colour pattern diversity
# Robert MacDonald

# Clear environment
rm(list=ls())

# Load libraries ----
library(dplyr)


## EDITABLE CODE ##
# Use latest IUCN assessment data or use most recent assessment data pre- specified cutoff year?
latest <- TRUE
# If not using latest assessment data, specify a cutoff year. Set to NULL if using latest.
cutoff_year <- NULL
# Clade to focus on ("Aves", "Neognaths", "Neoaves", "Passeriformes")
clade <- "Aves"

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

# Load raw species extinction matrices (1=extant, 0=extinct) for 100% abatement scope
extinctions = read.csv(
  here::here("01_input_data", "all_pamat_all_threats_1000_17122024.csv"), 
  row.names = 1
  ) # 100% scope


# Analysis ----

# calculate distance to centroid from colourspace
centr_dists <- dispRity::dispRity(colour_space, metric = dispRity::centroids)$disparity[[1]][[1]]
centr_dists <- data.frame(
  species = sapply(strsplit(rownames(colour_space), split = "-"), "[", 1),
  sex = sapply(strsplit(rownames(colour_space), split = "-"), "[", 2),
  centr_dists = centr_dists
)

# get species and sex with colourspace values
colour_space_sppsex <- data.frame(species = sapply(strsplit(rownames(colour_space), split = "-"), "[", 1),
                                  sex = sapply(strsplit(rownames(colour_space), split = "-"), "[", 2), 
                                  colour_space)

# trim extinctions data to only species inluded in colour dataset
ext <- extinctions[, colnames(extinctions) %in% colour_space_sppsex$species]

# trim colour space data to only species included in extinctions dataset
colour_space_sppsex <- colour_space_sppsex[colour_space_sppsex$species %in% colnames(ext), ]

# trim centroid distance data to only species included in extincitons dataset
centr_dists <- centr_dists[centr_dists$species %in% colnames(ext), ]

# put centroid distance rows in the same order as extinction matrix columns
ext <- ext[, order(colnames(ext))]
centr_dists <- centr_dists[order(centr_dists$species), ]

# get extinction matrix with centroid distances as values, with NA if species are extinct
# separately for males and females
ext_centr_matrix_m <- as.data.frame(t(t(as.matrix(ext)) * as.vector(centr_dists[centr_dists$sex == "M", "centr_dists"])))
ext[ext == 0] <- NA
# I can now use rowMeans to get the mean distance to centroid in
# each scenario, which I can compare to the baseline





############################################
# based on Kerry's approach

# get species and sex with colourspace values
colour_space_sppsex <- data.frame(species = sapply(strsplit(rownames(colour_space), split = "-"), "[", 1),
                                  sex = sapply(strsplit(rownames(colour_space), split = "-"), "[", 2), 
                                  colour_space)

# trim extinctions data to only species inluded in colour dataset
ext <- extinctions[, colnames(extinctions) %in% colour_space_sppsex$species]

# trim colour space data to only species included in extinctions dataset
colour_space_sppsex <- colour_space_sppsex[colour_space_sppsex$species %in% colnames(ext), ]

# Remove the 'full' (extant) row from the extinction matrix
ext = ext[rownames(ext) != "full", ]

# Calculate **Species Loss Avoided** (tots) for 100% scope
# sp_loss = (Total Extant Species - Scenario Species Remaining)
tots = as.data.frame(cbind(
  rownames(ext), # Scenario column
  stringr::str_split_i(as.character(rownames(ext)), "_", 1), # Extract the scenario/driver Code
  (rowSums(unique(extinctions[rownames(extinctions) == "full", -1])) -  rowSums(ext[, -1])) # # Sum of all extant species - Sum of species remaining in this scenario
))

# Assign column names and scope
colnames(tots) = c("Scenario", "Code", "sp_loss")
tots$scope = "100%"

# Rename the extinction scenario from 'none' (in file name) to 'comp' (complete extinction)
tots$Code = gsub("none", "comp", tots$Code)
tots$sp_loss = as.numeric(tots$sp_loss)
# Store the raw Species Loss Avoided for later use
tots$sp_loss_raw = tots$sp_loss

# Group and summarize the raw Species Loss Avoided
tots %>% group_by(Code) |> 
  subset(Scenario != "full") |> 
  summarise(
    mean_sp_loss = mean(sp_loss),
    sd_sp_loss = sd(sp_loss)
  )

# Calculate **Species Extinctions Averted** (the metric used in the paper's final analysis)
# Averted = Baseline Loss ('comp') - Scenario Loss Avoided
# The 'times = 21' accounts for 7 codes (cli, dis, hab, exp, inv, pol, all) * 3 scopes (100%, 50%, 10%).
# Note: 'sp_loss' here now represents the number of species **averted** (saved).
tots[tots$Code != "comp", "sp_loss"] <- rep(tots[tots$Code == "comp", "sp_loss"], times = 6) - (tots[tots$Code != "comp", "sp_loss"])

# Create a separate data frame for percentage averted *relative to the Baseline Loss*
tots_perc_comp=tots

# Calculate the **Proportion of Species Extinctions Averted** relative to the total projected loss ('comp').
tots_perc_comp[tots_perc_comp$Code != "comp", "sp_loss"] = (tots_perc_comp[tots_perc_comp$Code != "comp", "sp_loss"])/rep(tots_perc_comp[tots_perc_comp$Code == "comp", "sp_loss"], times = 6)

# Group and summarize the Averted Species data (Absolute number averted)
tots %>%
  dplyr::group_by(Code) %>%
  dplyr::summarise(SD = sd(sp_loss),
                   Value = mean(as.numeric(sp_loss)))


# Summarize Averted Species data (Absolute number averted and raw loss avoided)
tots_sum = tots %>%
  dplyr::group_by(Code) %>%
  dplyr::summarise(SD = sd(sp_loss),
                   Value = mean(as.numeric(sp_loss)),
                   SD_withres = sd(sp_loss_raw),
                   Value_withres = mean(sp_loss_raw))
