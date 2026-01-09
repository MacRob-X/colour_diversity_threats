# Compare mean colour pattern diversity for species threatened by different IUCN threat groups
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

# Choose axes to loop over - use 7, as parallel analysis indicates this many statistically significant dimensions
axes <- paste0("PC", 1:7)

# get species and sex with colourspace values
colour_space_sppsex <- data.frame(species = sapply(strsplit(rownames(colour_space), split = "-"), "[", 1),
                                  sex = sapply(strsplit(rownames(colour_space), split = "-"), "[", 2), 
                                  colour_space)

# Join threat data to colourspace values
threat_vals <- threat_matrix |> 
  inner_join(colour_space_sppsex, by = join_by("jetz_species" == "species"))
# this throws a warning but it's just because we have male and female colourspace values together
# - it's not a problem

# remove duplicates based on second-order code
threat_vals <- threat_vals |> 
  distinct(second_ord_code, jetz_species, sex, .keep_all = TRUE)

# convert to long data so I can boxplot with PC facet
threat_vals_long <- tidyr::pivot_longer(
  threat_vals,
  cols = tidyr::starts_with("PC"),
  names_to = "PC",
  values_to = "val"
)

# make boxplot, faceted by sex and PC
threat_vals_long |> 
  filter(
    !is.na(sex),
    !is.na(ex_driver),
    #   iucn_cat != "LC",
    PC %in% axes[2]
  ) |> 
  ggplot(aes(x = ex_driver, y = val, fill = ex_driver)) + 
  geom_boxplot(outliers = F, notch = T) + 
  facet_grid(rows = vars(PC), cols = vars(sex))

# loop over PC axes of choice and make models
pc_aov_mods <- vector("list", length(axes))
names(pc_boxplots) <- axes
pc_tukey_mods <- vector("list", length(axes))
names(pc_boxplots) <- axes
for(axis in axes){
  
  dat <- threat_vals[, c("jetz_species", "second_ord_code", "ex_driver", "sex", axis)]
  colnames(dat) <- c("jetz_species", "second_ord_code", "ex_driver", "sex", "pc_vals")
  
    # boxplot
  bp <- dat |> 
    filter(
      !is.na(sex),
      !is.na(ex_driver)
      #   iucn_cat != "LC"
    ) |> 
    ggplot(aes(x = ex_driver, y = pc_vals, fill = ex_driver)) + 
    geom_boxplot(outliers = F) + 
    facet_grid(~ sex)
  
  pc_boxplots[[axis]] <- bp 
  
} 
