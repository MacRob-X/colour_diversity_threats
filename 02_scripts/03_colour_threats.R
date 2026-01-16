# Compare mean colour pattern diversity for species threatened by different IUCN threat groups
# Robert MacDonald

# Clear environment
rm(list=ls())

# Load libraries ----
library(dplyr)
library(ggplot2)
library(ggridges)

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

# Analysis ----

# add second-order threat codes to threat matrix (derived from third-order codes)
threat_matrix$second_ord_code <- stringr::str_extract(threat_matrix$code, "[^_]*_[^_]*")


# assign second-order IUCN threat types to grouped 'driver of extinction' categories
# same system as Stewart et al 2025 Nat Ecol Evol - grouping is provided in Supplementary Dataset 1
# of that paper
# There are many threats that aren't assigned to one of these groups - this is because these
# threats were non-significant in predicting IUCN threat level in Stewart et al 2025
# We assign these as 'FLAG' in case we want to do anything with them later

# Groups
# Accidental mortality and disturbance
acc_mort_codes <- c("4_2", "5_4", "6_3")
# Climate change and severe weather
clim_chan_codes <- c("11_1", "11_4")
# Habitat loss and degradation
hab_loss_codes <- c("1_2", "1_3", "2_1", "2_2", "2_3", "5_3", "7_1", "7_2")
# Hunting and collecting
hunt_col_codes <- "5_1"
# Invasive species and disease
invas_spec_codes <- c("8_1", "8_2")
# Other ["Threats that affected ten or fewer species were grouped with other threats"]
other_codes <- c("10_1", "10_2", "10_3", "12_1")
# Pollution
pollut_codes <- "9_3"

threat_matrix <- threat_matrix |> 
  mutate(
    ex_driver = second_ord_code
  ) |> 
  mutate( # surely there's a more elegant way to do this
    ex_driver = ifelse(
      ex_driver %in% acc_mort_codes,
      "acc_mort",
      ifelse(
        ex_driver %in% clim_chan_codes,
        "clim_chan",
        ifelse(
          ex_driver %in% hab_loss_codes,
          "hab_loss",
          ifelse(
            ex_driver %in% hunt_col_codes,
            "hunt_col",
            ifelse(
              ex_driver %in% invas_spec_codes,
              "invas_spec",
              ifelse(
                ex_driver %in% other_codes,
                "other",
                ifelse(
                  ex_driver %in% pollut_codes,
                  "pollut",
                  ifelse(
                    is.na(ex_driver),
                    NA,
                    "FLAG"
                  )
                )
              )
            )
          )
        )
      )
    )
  )


# For now, let's make all the flagged threats (i.e. those which are not significant predictors
# of extinction risk) NA, as we're not interested them
threat_matrix <- threat_matrix |> 
  mutate(
    ex_driver = ifelse(
      ex_driver == "FLAG",
      NA,
      ex_driver
    )
  )

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

# Plot bar chart of threat types 
threat_matrix |>
  ggplot(aes(x = ex_driver)) + 
  geom_bar()


# calculate distance to centroid from colourspace
centr_dists <- dispRity::dispRity(colour_space, metric = dispRity::centroids)$disparity[[1]][[1]]
centr_dists <- data.frame(
  species = sapply(strsplit(rownames(colour_space), split = "-"), "[", 1),
  sex = sapply(strsplit(rownames(colour_space), split = "-"), "[", 2),
  centr_dists = centr_dists
)

# add distance to centroid onto threat matrix
threat_centr <- threat_matrix |> 
  inner_join(centr_dists, by = join_by("jetz_species" == "species"))
# this throws a warning but it's just because we have male and female centroid distance data together
# - it's not a problem

# remove duplicates based on second-order code
threat_centr_clean <- threat_centr |> 
  distinct(second_ord_code, jetz_species, sex, .keep_all = TRUE)

# boxplot of centroid distances by threat type - exclude LC species
threat_centr_clean |> 
  filter(
    !is.na(sex),
    !is.na(ex_driver)
    #   iucn_cat != "LC"
  ) |> 
  ggplot(aes(x = ex_driver, y = centr_dists, fill = ex_driver)) + 
  geom_boxplot(outliers = F) + 
  facet_grid(~ sex)

# ANOVA to check if mean centroid distances of each extinction driver are different
an_m <- aov(centr_dists ~ ex_driver, data = threat_centr_clean[threat_centr_clean$sex == "M", ])
summary(an_m)
# post-hoc Tukey test
lm_m <- lm(centr_dists ~ ex_driver, data = threat_centr_clean[threat_centr_clean$sex == "M", ])
av_m <- aov(lm_m)
tukey_centr_dists_m <- TukeyHSD(av_m)

an_f <- aov(centr_dists ~ ex_driver, data = threat_centr_clean[threat_centr_clean$sex == "F", ])
summary(an_f)
# post-hoc Tukey test
lm_f <- lm(centr_dists ~ ex_driver, data = threat_centr_clean[threat_centr_clean$sex == "F", ])
av_f <- aov(lm_f)
tukey_centr_dists_f <- TukeyHSD(av_f)

# Individual PC axes ----

# Plot distribution of centroid distances for each threat type

# Density plot of all species' centroid distance (excluding LC) vs species threatened by each threat type
threat_centr_clean |> 
  bind_rows(
    threat_centr_clean |> 
      mutate(
        ex_driver = "all"
      )
  ) |> 
  filter(
    iucn_cat != "LC",
    !(ex_driver %in% c(NA, "other"))
  ) |> 
  filter(
    ex_driver %in% c("all", "hunt_col")
  ) |> 
  ggplot(aes(x = centr_dists, fill = ex_driver)) + 
  geom_density(alpha = 0.4) + 
  facet_wrap(~ sex)

# Get threat data with PC axes
colour_space_sppsex <- data.frame(species = sapply(strsplit(rownames(colour_space), split = "-"), "[", 1),
                                  sex = sapply(strsplit(rownames(colour_space), split = "-"), "[", 2), 
                                  colour_space)
threat_colour <- threat_matrix |> 
  inner_join(colour_space_sppsex, by = join_by("jetz_species" == "species"))
# this throws a warning but it's just because we have male and female data together
# - it's not a problem

# remove duplicates based on second-order code
threat_colour_clean <- threat_colour |> 
  distinct(second_ord_code, jetz_species, sex, .keep_all = TRUE)

# Pivot longer, so we can facet by PC axis
threat_colour_long <- threat_colour_clean |> 
  tidyr::pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "PC_value"
  ) 

# Density plot of all species (excluding LC) vs species threatened by each threat type,
# faceted by PC axis (first 7 PCs only)
threat_colour_long |> 
  bind_rows(
    threat_colour_long |> 
      mutate(
        ex_driver = "all"
      )
  ) |> 
  filter(
    iucn_cat != "LC",
    !(ex_driver %in% c(NA, "other")),
    PC %in% paste0("PC", 4)
  ) |> 
  mutate(
    ex_driver = factor(ex_driver, levels = c("all", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
    sex = factor(sex, levels = c("M", "F"))
  ) |> 
  filter(
    ex_driver %in% c("all", "hunt_col")
  ) |> 
  ggplot(aes(x = PC_value, fill = ex_driver)) + 
  geom_density(alpha = 0.4) + 
  geom_segment(x = 0, xend = 0, y = 0, yend = 0.04, lwd = 0.3, colour = "grey30") + # Vertical line at 0
  facet_grid(rows = vars(PC), cols = vars(sex)) + 
  theme_minimal()


# Same plot, but with absolute PC values - this is effectively distance to centroid of each PC
threat_colour_long |> 
  bind_rows(
    threat_colour_long |> 
      mutate(
        ex_driver = "all"
      )
  ) |> 
  filter(
    iucn_cat != "LC",
    !(ex_driver %in% c(NA, "other")),
    PC %in% paste0("PC", 3:3)
  ) |> 
  mutate(
    ex_driver = factor(ex_driver, levels = c("all", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
    sex = factor(sex, levels = c("M", "F"))
  ) |> 
  filter(
    ex_driver %in% c("all", "hunt_col")
  ) |> 
  mutate(
    PC_value = abs(PC_value)
  ) |> 
  ggplot(aes(x = PC_value, fill = ex_driver)) + 
  geom_density(alpha = 0.4) + 
  geom_segment(x = 0, xend = 0, y = 0, yend = 0.04, lwd = 0.3, colour = "grey30") + # Vertical line at 0
  facet_grid(rows = vars(PC), cols = vars(sex)) + 
  theme_minimal()


# I can do a Levene's test to test for equality of variances between two distributions - this will
# tell me if one distribution is more clustered around the mean than another
# This is useful for the PC axes (e.g. PC1) where it looks visually like the hunting and collection-threatened
# species are generally more dispersed, rather than being concentrated towards one end of the axes 
# (so it doesn't look like there's a difference in mean, just variance)
# Levene's test is robust to non-normal distributions
pc_stats_df <- threat_colour_clean |> 
  bind_rows(
    threat_colour_clean |> 
      mutate(
        ex_driver = "all"
      )
    ) |> 
  filter(
    iucn_cat != "LC",
    !(ex_driver %in% c(NA, "other"))
  )

levene_M <- car::leveneTest(PC1 ~ ex_driver, data = pc_stats_df[pc_stats_df$ex_driver %in% c("all", "hunt_col") & pc_stats_df$sex == "M", ])
levene_F <- car::leveneTest(PC1 ~ ex_driver, data = pc_stats_df[pc_stats_df$ex_driver %in% c("all", "hunt_col") & pc_stats_df$sex == "F", ])
levene_M
levene_F

# And we can also run a simple ANOVA to test if there's difference in the means of the two distributions
aov_m <- aov(PC7 ~ ex_driver, data = pc_stats_df[pc_stats_df$ex_driver %in% c("all", "hunt_col") & pc_stats_df$sex == "M", ])
aov_f <- aov(PC7 ~ ex_driver, data = pc_stats_df[pc_stats_df$ex_driver %in% c("all", "hunt_col") & pc_stats_df$sex == "F", ])
summary(aov_m)
summary(aov_f)

# Loop over the first 7 axes, performing ANOVA and Levene's test each time, for each pairing of 
# all-specific driver
threat_types <- c(hab_loss = "hab_loss", hunt_col = "hunt_col", invas_spec = "invas_spec", clim_chan = "clim_chan", pollut = "pollut", acc_mort = "acc_mort")
pc_res <- vector("list", length = length(threat_types))
names(pc_res) <- threat_types
lev_res <- rep(list(pc_res), times = 7)
str(lev_res)
names(lev_res) <- paste0("PC", 1:7)

lev_res <- lapply(
  paste0("PC", 1:7),
  function(pc_axis){
    
    pc_res <- lapply(
      threat_types,
      function(threat){
        
        form <- formula(get(pc_axis) ~ ex_driver)
        
        aov_m <- aov(form, data = pc_stats_df[pc_stats_df$ex_driver %in% c("all", threat) & pc_stats_df$sex == "M", ])
        aov_f <- aov(form, data = pc_stats_df[pc_stats_df$ex_driver %in% c("all", threat) & pc_stats_df$sex == "F", ])
        levene_m <- car::leveneTest(PC1 ~ ex_driver, data = pc_stats_df[pc_stats_df$ex_driver %in% c("all", "hunt_col") & pc_stats_df$sex == "M", ])
        levene_f <- car::leveneTest(PC1 ~ ex_driver, data = pc_stats_df[pc_stats_df$ex_driver %in% c("all", "hunt_col") & pc_stats_df$sex == "F", ])
        
        res <- list(
          aov_m, aov_f, levene_m, levene_f
        )
        names(res) <- c("aov_m", "aov_f", "levene_m", "levene_f")
        
        return(res)
        
      }
    )
    names(pc_res) <- threat_types
    return(pc_res)
  }
)
names(lev_res) <- paste0("PC", 1:7)

# Can then access the individual PC/threat type results
# E.g. for PC1, to look at the difference in distributions in species threatened by all threat types
# vs only species threatened by hunting and collection
summary(lev_res$PC1$hunt_col$aov_m)
summary(lev_res$PC1$hunt_col$aov_f)
lev_res$PC1$hunt_col$levene_m
lev_res$PC1$hunt_col$levene_f
# So, we can see that the species which are threatened by hunting and collection tend to 
# be further out towards the edges of the first PC (as shown by the Levene's test), but not in either 
# direction in particular (as shown by ANOVA)
# So we interpret this to mean that hunting and collection disproportionately threatens both light
# and dark birds
# for PC3, look at the difference in distributions in species threatened by all threat types
# vs only species threatened by hunting and collection
summary(lev_res$PC3$hunt_col$aov_m)
summary(lev_res$PC3$hunt_col$aov_f)
lev_res$PC3$hunt_col$levene_m
lev_res$PC3$hunt_col$levene_f
# for PC3 (countershading vs reverse countershading), we see that species threatened by hunting and 
# collection tend to be further from the centroid, and tend to have more positive PC3 values - 
# meaning that hunting and collection tends to threaten reverse-countershaded species more
# than countershaded species
