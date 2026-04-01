# Troubleshoot where the gradient in Kerry's simulation data is coming from

# Comparison between SR/colour loss in the real simulated data, and that when we randomise
# centroid distance across species

# Clear environment
rm(list=ls())

# Load libraries ----
library(dplyr)
library(ggplot2)

# Load custom functions ----
source(
  here::here(
    "02_scripts", "R", "04_threat_abatements_functions.R"
  )
)


## EDITABLE CODE ##
# Use latest IUCN assessment data or use most recent assessment data pre- specified cutoff year?
latest <- TRUE
# If not using latest assessment data, specify a cutoff year. Set to NULL if using latest.
cutoff_year <- NULL
# Clade to focus on ("Aves", "Neognaths", "Neoaves", "Passeriformes")
clade <- "Aves"
# Colour space to use (e.g. "lab" [CIELAB], "srgb" [sRGB], "xyz" [TCSxyz], "jndxyzlum" [JND with xyz and luminance])
space <- "lab"

# Load data ----

# load colour pattern space (created in Chapter 1 - patch-pipeline)
colspace_path <- paste0("G:/My Drive/patch-pipeline/2_Patches/3_OutputData/", clade, "/2_PCA_ColourPattern_spaces/1_Raw_PCA/", clade, ".matchedsex.patches.250716.PCAcolspaces.rds")
colour_space <- readRDS(colspace_path)[[space]][["x"]]

# Load raw species extinction matrices (1=extant, 0=extinct) for 100% abatement scope
extinctions = read.csv(
  here::here("01_input_data", "all_pamat_all_threats_1000.csv"), 
  row.names = 1
) # 100% scope

# Load IUCN threat categories
iucn <- read.csv(
  "G:/My Drive/patch-pipeline/4_SharedInputData/iucn_data/aves_iucn_2025_nominate.csv"
)

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

# convert to true matrix
ext <- as.matrix(ext)

# Remove the 'full' (extant) row from the extinction matrix
# ext = ext[rownames(ext) != "full", ]

# get extinction matrix with centroid distances as values, with NA if species are extinct
# separately for males and females
cd_vector <- centr_dists[centr_dists$sex == "M", "centr_dists"]
# multiple extinction matrix by centroid distance vector (this is an optimised way to run this
# matrix/vector calculation)
ext_centr_matrix_m <- sweep(ext, 2, cd_vector, `*`)
ext_centr_matrix_m[ext_centr_matrix_m == 0] <- NA
# I can now use rowMeans to get the mean distance to centroid in
# each scenario, which I can compare to the baseline
ext_means <- rowMeans(ext_centr_matrix_m, na.rm = T)
names(ext_means) <- rownames(ext_centr_matrix_m)

# Now do the same, but with randomised centroid distances (to break the link between 
# threat status and distance to centroid)
# first get list of species for which we have both extinction and colour data
spp_list <- colnames(ext)
cd_vector_randomised <- sample(centr_dists[centr_dists$sex == "M", "centr_dists"], size = length(spp_list), replace = FALSE)
# double check this contains exactly the same values as cd_vector but in a different order
identical(cd_vector, cd_vector_randomised) # expect FALSE
identical(cd_vector[order(cd_vector)], cd_vector_randomised[order(cd_vector_randomised)]) # expect TRUE
# multiple extinction matrix by centroid distance vector (this is an optimised way to run this
# matrix/vector calculation)
ext_centr_matrix_m_random <- sweep(ext, 2, cd_vector_randomised, `*`)
ext_centr_matrix_m_random[ext_centr_matrix_m_random == 0] <- NA
# I can now use rowMeans to get the mean distance to centroid in
# each scenario, which I can compare to the baseline
ext_means_random <- rowMeans(ext_centr_matrix_m_random, na.rm = T)
names(ext_means_random) <- rownames(ext_centr_matrix_m_random)

# check the random centroid distance have the same distribution as the true centroid distances
# (since the values should be exactly the same, jsut assigned to different species)
plot(density(cd_vector))
plot(density(cd_vector_randomised))
# all good

# check how the distribution of the means differs
plot(density(ext_means))
plot(density(ext_means_random))

# Do the same for species richness
ext_sr_values <- rowSums(ext)
names(ext_sr_values) <- rownames(ext)


# Calculate mean distance to centroid loss (as a percentage of the full assemblage mean
# distance to centroid)
# first sanity check that there are no extinctions in the 'full' scenario
all(extinctions["full", ] == 1)
all(ext["full", ] == 1)
# all good

# first get the mean distance to centroid and species richness of the full extant assemblage
ext_cd_mean_full <- ext_means["full"]
ext_sr_full <- ext_sr_values["full"]

# now calculate the loss percentage as the mean distance to centroid for each abatement simulation,
# as a percentage of full extant assemblage centroid distance
# compared to EXTANT ASSEMBLAGE - so this is actual loss, not loss avoided
mean_cd_loss_abs = ext_cd_mean_full - ext_means
random_cd_loss_abs = ext_cd_mean_full - ext_means_random
sr_loss_abs <- ext_sr_full - ext_sr_values
# create dataframe with results
centr_dist_sims <- data.frame(
  code = gsub('(.*)_\\w+', '\\1', names(ext_means)),
  mean_centr_dist = ext_means,
  mean_centr_dist_random = ext_means_random,
  mean_cd_loss_abs,
  random_cd_loss_abs,
  sr = ext_sr_values,
  sr_loss_abs
)
# remove the 'full' row, as there's by definition zero loss
centr_dist_sims <- centr_dist_sims[centr_dist_sims$code != "full", ]

# Plot absolute species richness loss vs absolute species richness loss with Kerry's 
# data and the real centroid distances
centr_dist_sims |> 
  filter(
    code != "none"
  ) |> 
  ggplot(aes(x = sr_loss_abs, y = mean_cd_loss_abs, colour = code)) + 
  geom_point() + 
  geom_smooth(method = "lm")

# Plot absolute species richness loss vs absolute species richness loss with Kerry's 
# data and the RANDOM centroid distances
centr_dist_sims |> 
  filter(
    code != "none"
  ) |> 
  ggplot(aes(x = sr_loss_abs, y = random_cd_loss_abs, colour = code)) + 
  geom_point() + 
  geom_smooth(method = "lm")
