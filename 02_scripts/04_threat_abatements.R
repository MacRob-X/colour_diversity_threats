# Effects of abatement of threat groups on colour pattern diversity
# Robert MacDonald

# Clear environment
rm(list=ls())

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

# Remove the 'full' (extant) row from the extinction matrix
# ext = ext[rownames(ext) != "full", ]

# get extinction matrix with centroid distances as values, with NA if species are extinct
# separately for males and females
ext_centr_matrix_m <- as.data.frame(t(t(as.matrix(ext)) * as.vector(centr_dists[centr_dists$sex == "M", "centr_dists"])))
ext_centr_matrix_m[ext_centr_matrix_m == 0] <- NA
# I can now use rowMeans to get the mean distance to centroid in
# each scenario, which I can compare to the baseline
ext_means <- rowMeans(ext_centr_matrix_m, na.rm = T)
names(ext_means) <- rownames(ext_centr_matrix_m)

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
cd_pc_loss <- 100 * (ext_cd_mean_full - ext_means) / ext_cd_mean_full
sr_loss_abs <- ext_sr_full - ext_sr_values
ext_sr_loss_pc <- 100 * (ext_sr_full - ext_sr_values) / ext_sr_full
centr_dist_sims <- data.frame(
  code = gsub('(.*)_\\w+', '\\1', names(ext_means)),
  mean_centr_dist = ext_means,
  mean_cd_loss_abs,
  cd_pc_loss = cd_pc_loss,
  cd_loss_avoided_abs = NA,
  cd_loss_avoided_pc = NA,
  sr = ext_sr_values,
  sr_loss_abs,
  sr_pc_loss = ext_sr_loss_pc,
  sr_loss_avoided_abs = NA,
  sr_loss_avoided_pc = NA
)
# remove the 'full' row, as there's by definition zero loss
centr_dist_sims <- centr_dist_sims[centr_dist_sims$code != "full", ]


# Plot absolute species richness loss vs absolute species richness loss
centr_dist_sims |> 
  filter(
    code != "none"
  ) |> 
  ggplot(aes(x = sr_loss_abs, y = mean_cd_loss_abs, colour = code)) + 
  geom_point() + 
  geom_smooth(method = "lm")
# we can see that the line for 'exp' - hunting & collection - is lower than the others, which means
# colour pattern diversity loss per number of species lost is lower when hunting & collection is
# abated

# Calculate absolute mean distance to centroid loss avoided compared with no abatement scenario
# positive number means higher mean distance to centroid under specific abatement than no abatement
centr_dist_sims[centr_dist_sims$code != "none", "cd_loss_avoided_abs"] = rep(centr_dist_sims[centr_dist_sims$code == "none", "mean_cd_loss_abs"], times = 7) - centr_dist_sims[centr_dist_sims$code != "none", "mean_cd_loss_abs"]

# Calculate absolute species richness loss AVOIDED compared with no abatement scenario
# positive number mean higher species richness under specific abatement compared with no abatement
centr_dist_sims[centr_dist_sims$code != "none", "sr_loss_avoided_abs"] = rep(centr_dist_sims[centr_dist_sims$code == "none", "sr_loss_abs"], times = 7) - centr_dist_sims[centr_dist_sims$code != "none", "sr_loss_abs"]

# plot absolute centroid distance loss AVOIDED vs absolute species richness loss AVOIDED
centr_dist_sims |> 
  filter(
    code != "none"
  ) |> 
  ggplot(aes(x = sr_loss_avoided_abs, y = cd_loss_avoided_abs, colour = code)) + 
  geom_point() + 
  geom_smooth(method = "lm")
# we can see that the line for 'exp' (Hunting & collection) is above the other lines - so abatement of
# hunting and collection avoids disproportionately more loss of colour pattern diversity accounting for
# decrease in species richness than abatement of other specific threats

# Now get as percentages for easier comparison
# THIS IS WHAT WE WILL ACTUALLY USE IN THE MANUSCRIPT
centr_dist_sims[centr_dist_sims$code != "none", "cd_loss_avoided_pc"] <- rep(centr_dist_sims[centr_dist_sims$code == "none", "cd_pc_loss"], times = 7) - centr_dist_sims[centr_dist_sims$code != "none", "cd_pc_loss"]

# And the same for species richness loss
centr_dist_sims[centr_dist_sims$code != "none", "sr_loss_avoided_pc"] <- rep(centr_dist_sims[centr_dist_sims$code == "none", "sr_pc_loss"], times = 7) - centr_dist_sims[centr_dist_sims$code != "none", "sr_pc_loss"]

# plot percentage centroid distance loss avoided vs percentage species richness loss avoided
centr_dist_sims |> 
  filter(
    code != "none"
  ) |> 
  ggplot(aes(x = sr_loss_avoided_pc, y = cd_loss_avoided_pc, colour = code)) + 
  geom_point() + 
  geom_smooth(method = "lm")
# we can see that the line for 'exp' (Hunting & collection) is above the other lines - so abatement of
# hunting and collection avoids disproportionately more loss of colour pattern diversity accounting for
# decrease in species richness than abatement of other specific threats


# FOR DOTTED LINE SHOWING MEAN COLOUR PATTERN DIVERSITY LOSS AVOIDANCE PER SPECIES RICHNESS UNDER FULL ABATEMENT

# Get mean % species loss avoided in full abatement ('all') scenario
n_sp <- ncol(ext)
sp_loss_av_all <- centr_dist_sims |> 
  filter(
    code == "all"
  ) |> 
  pull(
    sr_loss_avoided_pc
  ) |> 
  mean()
# Get mean % colour pattern diversity loss avoided in full abatement ('all') scenario
cd_loss_av_all <- centr_dist_sims |> 
  filter(
    code == "all"
  ) |> 
  pull(
    cd_loss_avoided_pc
  ) |> 
  mean()
# Use to calculate gradient of line
all_ab_grad <- cd_loss_av_all / sp_loss_av_all

# Now plot the same (percentage centroid distance loss avoided vs percentage species richness loss avoided)
# but with means instead of full simulation distributions 
centr_dist_sims |> 
  filter(
    code != "none",
  #  code != "all"
  ) |> 
  group_by(code) |> 
  summarise(
    mean_cd_loss_avoided = mean(cd_loss_avoided_pc),
    sd_cd_loss_avoided = sd(cd_loss_avoided_pc),
    mean_sr_loss_avoided = mean(sr_loss_avoided_pc),
    sd_sr_loss_avoided = sd(sr_loss_avoided_pc)
  ) |> 
  ggplot(aes(x = mean_sr_loss_avoided, y = mean_cd_loss_avoided, fill = code)) + 
  geom_errorbar(aes(xmin = mean_sr_loss_avoided - 0.5*sd_sr_loss_avoided, xmax = mean_sr_loss_avoided + 0.5*sd_sr_loss_avoided), colour = "black", width = 0.005) + 
  geom_errorbar(aes(ymin = mean_cd_loss_avoided - 0.5*sd_cd_loss_avoided, ymax = mean_cd_loss_avoided + 0.5*sd_cd_loss_avoided), colour = "black", width = 0.02) + 
  geom_point(shape = 21, size = 5) + 
  #### DIAGONAL LINE - must go through intercept
   geom_abline(intercept = 0, slope = all_ab_grad, lty = 2) +
  xlab("Avoided species richness loss (%)") + ylab("Avoided colour pattern diversity loss (%)") + 
  labs(fill = "Driver-specific abatement") + 
  scale_fill_discrete(labels = c("All", "Climate", "Disturbance", "Hunting", "Habitat", "Invasive", "Pollution")) + 
#  coord_fixed() + # to fix the aspect ratio
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 2025-01-12
# It's strange that my percentage avoided species loss figure fro habitat loss abatement
# is about half the figure in Kerry's paper
# OK - having gone through my figures and Kerry's figures, I suspect it's because of the 
# relatively limited and possibly skewed species sampling of our colour dataset
# We are working with 7609 species (that we have both colour and extinction data for) compared with
# 9873 that Kerry has extinction data for
# If the missing ~2000 species tend to be those at higher risk of extinction, which is a fair
# assumption as these are probably less likely to be present in the NHM collection, then it could
# affect the absolute loss values we see with species extinctions
# However, this shouldn't be a problem as I am only interested in the colour pattern diversity loss relative
# to the species richness loss - not the absolute values
# It does mean I can't directly compare figures with Kerry's paper, but I can compare the patterns
# (and actually it looks like the patterns are very similar, comparing the above plot with
# Fig. 1b in Stewart et al 2025)
# What is still strange though is that the difference only seems to be present for the 'hab'-threatened
# species (so the species loss avoided under the 'hab' abatement scenario is about half that seen in 
# Kerry's paper, but not for the other ones)
# Have just used Kerry's code to calculate SR loss, but limited to species in the colour data - this
# confirms that it is the limited species sampling that is producing the discrepancy, but not why the discrepancy
# is much higher for 'hab' - it presumably must be because habitat loss is such a big threat, most of 
# the species which are both endangered and missing from the colour dataset are threatened by habitat loss
# Again, shouldn't be a problem as I'm only interested in the ratio of colour diversity loss to 
# species richness loss












# THE JUICE
# Calculate colour pattern diversity loss AVOIDED under each abatement scenario (excluding the 'full' 
# row, which is the entire extant assemblage, and the no abatement scenario - we calculate the
# loss avoided relative to the no abatement scenario)

centr_dist_sims[centr_dist_sims$code != "none", "cd_loss_avoided_pc"] <- rep(centr_dist_sims[centr_dist_sims$code == "none", "cd_pc_loss"], times = 6) - centr_dist_sims[centr_dist_sims$code != "none", "cd_pc_loss"]

# And the same for species richness loss
centr_dist_sims[centr_dist_sims$code != "none", "sr_loss_avoided_pc"] <- rep(centr_dist_sims[centr_dist_sims$code == "none", "sr_pc_loss"], times = 6) - centr_dist_sims[centr_dist_sims$code != "none", "sr_pc_loss"]


# plot colour pattern diversity loss vs species richness loss
centr_dist_sims |> 
  group_by(code) |> 
  summarise(
    mean_cd_loss = mean(cd_pc_loss),
    sd_cd_loss = sd(cd_pc_loss),
    mean_sr_loss = mean(sr_pc_loss),
    sd_sr_loss = sd(sr_pc_loss)
  ) |> 
  ggplot(aes(x = mean_sr_loss, y = mean_cd_loss, colour = code, fill = code)) + 
  geom_errorbar(aes(xmin = mean_sr_loss - sd_sr_loss, xmax = mean_sr_loss + sd_sr_loss), colour = "black") + 
  geom_errorbar(aes(ymin = mean_cd_loss - sd_cd_loss, ymax = mean_cd_loss + sd_cd_loss), colour = "black") + 
  geom_point(shape = 21, size = 5)

# plot avoided colour pattern diversity loss vs avoided species richness loss
centr_dist_sims |> 
  group_by(code) |> 
  summarise(
    mean_cd_loss_avd = mean(cd_loss_avoided_pc),
    sd_cd_loss_avd = sd(cd_loss_avoided_pc),
    mean_sr_loss_avd = mean(sr_loss_avoided_pc),
    sd_sr_loss_avd = sd(sr_loss_avoided_pc)
  ) |> 
  ggplot(aes(x = mean_sr_loss_avd, y = mean_cd_loss_avd, colour = code, fill = code)) + 
  geom_errorbar(aes(xmin = mean_sr_loss_avd - sd_sr_loss_avd, xmax = mean_sr_loss_avd + sd_sr_loss_avd), colour = "black") + 
  geom_errorbar(aes(ymin = mean_cd_loss_avd - sd_cd_loss_avd, ymax = mean_cd_loss_avd + sd_cd_loss_avd), colour = "black") + 
  geom_point(shape = 21, size = 5)

