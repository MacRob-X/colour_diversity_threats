# Compare mean colour pattern diversity for species threatened by different IUCN threat groups
# Robert MacDonald

# Clear environment
rm(list=ls())

# Load libraries ----
library(dplyr)
library(ggplot2)
library(ggridges)

# Load custom functions ----
source(
  here::here(
    "02_scripts", "R", "06_pc_axes_threat_functions.R"
  )
)

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

# load colour pattern space UMAP (created in Chapter 1 - patch-pipeline)
umap_path <- paste0("G:/My Drive/patch-pipeline/2_Patches/3_OutputData/", clade, "/2_PCA_ColourPattern_spaces/2_UMAP/", clade, ".matchedsex.patches.nn.25.mindist.0.1.lab.UMAP.rds")
umap <- readRDS(umap_path)[["layout"]]

# Data preparation ----


# Exclude Data Deficient (DD), Extinct (EX) and Extinct in the Wild (EW) species
threat_matrix <- threat_matrix |> 
  filter(
    !(iucn_cat %in% c("DD", "EX", "EW"))
  )

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


# First, let's make a version with a binary variable for each extinction driver
threat_matrix_bin <- threat_matrix |> 
  mutate(
    acc_mort = ifelse(second_ord_code %in% acc_mort_codes, 1, 0),
    clim_chan = ifelse(second_ord_code %in% clim_chan_codes, 1, 0),
    hab_loss = ifelse(second_ord_code %in% hab_loss_codes, 1, 0),
    hunt_col = ifelse(second_ord_code %in% hunt_col_codes, 1, 0),
    invas_spec = ifelse(second_ord_code %in% invas_spec_codes, 1, 0),
    other = ifelse(second_ord_code %in% other_codes, 1, 0),
    pollut = ifelse(second_ord_code %in% pollut_codes, 1, 0),
  ) |> 
  select(
    jetz_species, iucn_cat, acc_mort, clim_chan, hab_loss, hunt_col, invas_spec, other, pollut 
  ) |> 
  distinct() |> 
  group_by(jetz_species, iucn_cat) %>%
  # Apply the 'max' function across all remaining columns
  # This ensures if a threat is 1 in ANY row, it becomes 1 in the result
  summarise(across(everything(), max), .groups = 'drop')

# can make a long version of this
threat_bin_long <- threat_matrix_bin |> 
  tidyr::pivot_longer(
    c(acc_mort, clim_chan, hab_loss, hunt_col, invas_spec, other, pollut),
    names_to = "threat_type",
    values_to = "threat_present"
  )

# Alternate version with a categorical variable instead of binary
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
# threat_matrix <- threat_matrix |> 
#   mutate(
#     ex_driver = ifelse(
#       ex_driver == "FLAG",
#       NA,
#       ex_driver
#     )
#   )
# NOTE 18/02/2026
# I don't actually want to do the above snippet, as I need to treat these species differently to
# species for which we truly have no threat data (which will be NA)
# Instead I will mark these as "no_sig_driver" so that I can include them in downstream 
# analysis with confidence
threat_matrix <- threat_matrix |>
  mutate(
    ex_driver = ifelse(
      ex_driver == "FLAG",
      "no_sig_driver",
      ex_driver
    )
  )




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


# Plotting ----

# get threat data with UMAP axes
umap_sppsex <- data.frame(species = sapply(strsplit(rownames(umap), split = "-"), "[", 1),
                          sex = sapply(strsplit(rownames(umap), split = "-"), "[", 2), 
                          umap)


threat_umap <- threat_matrix |> 
  inner_join(umap_sppsex, by = join_by("jetz_species" == "species")) |> 
  rename(
    UMAP1 = X1,
    UMAP2 = X2
  )
# this throws a warning but it's just because we have male and female data together
# - it's not a problem

# remove duplicates based on second-order code
threat_umap_clean <- threat_umap |> 
  distinct(second_ord_code, jetz_species, sex, .keep_all = TRUE)

# add PC values
threat_umap_clean <- threat_umap_clean |> 
  inner_join(
    colour_space_sppsex,
    by = c("jetz_species" = "species", "sex" = "sex")
  )


# Proportional 2D density plots (The Juice) ----

# set focal extinction driver
focal_threat <- "hunt_col"
# set axes (PCs or UMAP axes)
ax_1 <- "PC1"
ax_2 <- "PC2"

# Get proportional density of species threatened by focal threat vs ALL other species
prop_dens <- prop_dens_2d(
  threat_umap_clean,
  focal_threat = focal_threat,
  x_axis = ax_1, y_axis = ax_2,
  threatened_spp_only = FALSE,
  n_bins = 200,
  return_params = TRUE
)


# plot as heatmap
plot_prop_dens_2d(
  prop_dens
)

# Plot multiple threats together
ex_drivers <- c("clim_chan", "hab_loss", "pollut", "hunt_col", "acc_mort", "invas_spec")
names(ex_drivers) <- ex_drivers
# First get the proportional densities
prop_densities <- pbapply::pblapply(
  ex_drivers,
  function(driver){
    
    prop_dens <- prop_dens_2d(
      threat_umap_clean,
      focal_threat = driver,
      x_axis = ax_1, y_axis = ax_2,
      threatened_spp_only = FALSE,
      n_bins = 200,
      return_params = TRUE
    )
    
  }
)
# Then plot on same figure
# First get global maximum density value
global_max_dens <- max(
  unlist(
    lapply(
      prop_densities,
      function(dens){
        max <- max(abs(dens$prop_2dkd$z))
      }
      )
  )
)
n_plots <- length(prop_densities)
n_rows <- 3
if((n_plots %% n_rows) == 0){
  n_cols <- n_plots / n_rows
} else {
  n_cols <- (n_plots %/% n_rows) + 1
}

# Save as PNG
png(
  filename = here::here(
    "04_output_plots", "06_pc_axes_threats", "01_proportional_density_plots",
    paste0("extinction_drivers_vs_all_species_", ax_1, ax_2, ".png")
  ), 
  width = 210, height = 297,
  units = "mm",
  res = 300
)
par(
  mfrow = c(n_rows, n_cols),
  mar = c(4.1, 4.1, 3.1, 2.1)
)
lapply(
  prop_densities,
  plot_prop_dens_2d,
  title = TRUE,
  max_density_val = global_max_dens
)
dev.off()


# UMAP density plot by threat
threat_umap_clean |> 
  bind_rows(
    threat_umap_clean |> 
      mutate(
        ex_driver = "all"
      )
  ) |> 
  filter(
    iucn_cat != "LC",
    !(ex_driver %in% c(NA, "other"))
  ) |> 
  mutate(
    ex_driver = factor(ex_driver, levels = c("all", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
    sex = factor(sex, levels = c("M", "F"))
  ) |> 
  # filter(
  #   ex_driver != "all"
  # ) |> 
  ggplot(aes(x = X1, y = X2)) + 
  geom_density_2d_filled() + 
  facet_wrap(~ ex_driver)


threat_colour_clean |> 
  bind_rows(
    threat_colour_clean |> 
      mutate(
        ex_driver = "all"
      )
  ) |> 
  filter(
    iucn_cat != "LC",
    !(ex_driver %in% c(NA, "other"))
  ) |> 
  mutate(
    ex_driver = factor(ex_driver, levels = c("all", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
    sex = factor(sex, levels = c("M", "F"))
  ) |> 
  # filter(
  #   ex_driver != "all"
  # ) |> 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_density_2d_filled() + 
  facet_wrap(~ ex_driver)

# boxplot of PC values by threat type - exclude LC species
# Compare to all threatened species
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
    PC %in% paste0("PC", 2)
  ) |> 
  mutate(
    ex_driver = factor(ex_driver, levels = c("all", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
    sex = factor(sex, levels = c("M", "F"))
  ) |> 
  ggplot(aes(x = ex_driver, y = PC_value, fill = ex_driver)) + 
  geom_boxplot(outliers = F) + 
  facet_grid(~ sex)

# Density plot of all species (excluding LC) vs species threatened by each threat type,
# faceted by PC axis (first 7 PCs only)
threat_colour_long |> 
  bind_rows(
    threat_colour_long |> 
      filter(
        ex_driver != "hunt_col"
      ) |> 
      mutate(
        ex_driver = "all"
      ) |> 
      distinct()
  ) |> 
  filter(
    iucn_cat != "LC",
    !(ex_driver %in% c(NA, "other")),
    PC %in% paste0("PC", 1)
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

# function to do the above
plot_pc_distrib <- function(pc_axis, threat_type, long_data, absolute = FALSE) {
  
  if(absolute == TRUE){
    long_data$PC_value <- abs(long_data$PC_value)
  }
  
  p <- long_data |> 
    bind_rows(
      long_data |> 
        filter(
          ex_driver != threat_type
        ) |> 
        mutate(
          ex_driver = "all"
        ) |> distinct()
    ) |> 
    filter(
      iucn_cat != "LC",
      !(ex_driver %in% c(NA, "other")),
      PC == pc_axis
    ) |> 
    mutate(
      ex_driver = factor(ex_driver, levels = c("all", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
      sex = factor(sex, levels = c("M", "F"))
    ) |> 
    filter(
      ex_driver %in% c("all", threat_type)
    ) |> 
    ggplot(aes(x = PC_value, fill = ex_driver)) + 
    geom_density(alpha = 0.4) + 
    scale_fill_discrete(name = "Threats", labels = c("All threats", "Specific threat")) + 
    labs(title = names(threat_types)[threat_types == threat_type]) + 
    geom_segment(x = 0, xend = 0, y = 0, yend = 0.04, lwd = 0.3, colour = "grey30") + # Vertical line at 0
    facet_grid(cols = vars(sex)) + 
    theme_minimal() + 
    theme(
      #      legend.position = "none",
      plot.title = element_text(size = 10, vjust = -4, hjust = -0.1),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    )
  
  return(p)
}

threat_types <- c(hab_loss = "hab_loss", hunt_col = "hunt_col", invas_spec = "invas_spec", clim_chan = "clim_chan", pollut = "pollut", acc_mort = "acc_mort")
pc_axis <- "PC2"
pc_plots <- lapply(
  threat_types,
  plot_pc_distrib,
  pc_axis = pc_axis,
  long_data = threat_colour_long,
  absolute = FALSE
)
p_pc1_allthreats <- ggpubr::ggarrange(plotlist = pc_plots, common.legend = TRUE, ncol = 2, nrow = 3)
ggpubr::annotate_figure(
  p_pc1_allthreats,
  bottom = ggpubr::text_grob(pc_axis),
  left = ggpubr::text_grob("Density", rot = 90)
)


# Same plot, but with absolute PC values - this is effectively distance to centroid of each PC
abs_pc_plots <- lapply(
  threat_types,
  plot_pc_distrib,
  pc_axis = pc_axis,
  long_data = threat_colour_long,
  absolute = TRUE
)
abs_allthreats <- ggpubr::ggarrange(plotlist = abs_pc_plots, common.legend = TRUE, ncol = 2, nrow = 3)
ggpubr::annotate_figure(
  abs_allthreats,
  bottom = ggpubr::text_grob(pc_axis),
  left = ggpubr::text_grob("Density", rot = 90)
)


# Statistics ----

# Perform a Cramer-von Mises test to check for difference in distribution of PC1 values between 
# species threatened by hunting and collection and species not threatened by hunting and collection
# Exclude LC species

pc_axis <- "PC1"
focal_driver <- "hunt_col"

focal_distrib <- threat_colour_long |> 
  filter(
    iucn_cat != "LC", 
    PC == pc_axis, 
    ex_driver == focal_driver
    ) |> 
  pull(PC_value)
non_focal_distrib <- threat_colour_long |> 
  filter(
    iucn_cat != "LC", 
    PC == pc_axis, 
    ex_driver != focal_driver
    ) |> 
  pull(PC_value)

cvm <- twosamples::cvm_test(
  a = focal_distrib,
  b = non_focal_distrib,
)
cvm
plot(cvm)

# Calculate Standard Effect Size of CVM test statistic by randomly sampling from the
# non-hunting&collection distribution 
# Note that I think this is essentially what twosamples::cvm_test does under the hood anyway
# first calculate CVM test statistic for the observed distribution
cvm_obs <- twosamples::cvm_stat(
  a = focal_distrib,
  b = non_focal_distrib,
)
all_distrib <- threat_colour_long |> 
  filter(
    iucn_cat != "LC",
    PC == pc_axis
  ) |> 
  pull(PC_value)
random_samples <- lapply(1:2000, 
                         function(i){
                           sample(
                             all_distrib,
                             size = length(focal_distrib),
                             replace = FALSE
                           )
                         }
                         ) 
cvm_stats_random <- unlist(
  lapply(
    random_samples,
    function(sample){
      cvm <- twosamples::cvm_stat(
        a = sample,
        b = non_focal_distrib,
      )
      return(cvm)
    }
  )
)

# Now calculate SES by calculating (obs - mean(null)) / sd(null)
# SES > 1.96 indicates deviation from difference between distributions expected by chance
ses <- (cvm_obs - mean(cvm_stats_random)) / sd(cvm_stats_random)
# plot histogram of random sample distribution with abserved value as red dotted line
hist(cvm_stats_random, breaks = 50)
abline(v = cvm_obs, col = "red")
# compare with the bootstrap distribution of twosamples::cvm_test
plot(cvm)
# it looks like my SES code is doing essentially the same thing as twosamples::cvm_test
# So, no need to manually do the random sampling myself - I can let the premade function do it



# Now apply across each extinction driver and PC
# could parallelise this but it only takes ~30s for 7 axes
# define threat types and focal PCs first 
threat_types <- c(hab_loss = "hab_loss", hunt_col = "hunt_col", invas_spec = "invas_spec", clim_chan = "clim_chan", pollut = "pollut", acc_mort = "acc_mort")
pcs <- paste0("PC", 1:3)

# Run across PCs/extinction drivers
pc_ex_drive_res <- pbapply::pblapply(
  
  threat_types,
  
  function(threat_type){
    
    cvm <- pbapply::pblapply(
      pcs, 
      test_pc_threat,
      threat_colour_long = threat_colour_long,
      threat = threat_type,
      n = 1000,
      seed = 42,
      return_distrib = FALSE,
      stat_test = "dts"
    )
    names(cvm) <- pcs
    cvm <- data.table::rbindlist(cvm, idcol = "PC")
    # Bonferroni correction for p-values
    # cvm$p_adjusted <- p.adjust(cvm$p_value, method = "bonferroni")
    return(cvm)
    
  }
  
)
pc_ex_drive_res <- do.call(rbind, pc_ex_drive_res)
# Bonferroni correction for p-values
pc_ex_drive_res$p_adjusted <- p.adjust(pc_ex_drive_res$p_value, method = "bonferroni")

# compare to results just from using twosamples::wass_test (which uses essentially the
# same permutation test to determine a count-based p-value)
res_ts <- pbapply::pblapply(
  
  threat_types,
  
  function(threat_type){
    
    pc_stats <- pbapply::pblapply(
      pcs,
      function(pc_axis){
        
        focal_distrib <- threat_colour_long |> 
          filter(
            iucn_cat != "LC", 
            PC == pc_axis, 
            ex_driver == threat_type
          ) |> 
          pull(PC_value)
        non_focal_distrib <- threat_colour_long |> 
          filter(
            iucn_cat != "LC", 
            PC == pc_axis, 
            ex_driver != threat_type
          ) |> 
          pull(PC_value)
        
        test_stat <- twosamples::wass_test(
          focal_distrib,
          non_focal_distrib,
          nboots = 1000,
          keep.boots = F
        )
        test_stat <- as.data.frame(t(as.data.frame(test_stat)))
        test_stat$PC <- pc_axis
        return(test_stat)
        
      }
        )
    
    threat_res <- do.call(rbind, pc_stats)
  }
  
)
res_ts <- data.table::rbindlist(res_ts, idcol = "extinction_driver")
colnames(res_ts) <- c("extinction_driver", "test_stat", "p_value", "PC")
res_ts <- res_ts |> 
  select(PC, extinction_driver, test_stat, p_value)

ordered_ts_res <- res_ts |> 
  arrange(p_value)
ordered_wass_res <- wass_res_new |> 
  arrange(p_value)
identical(ordered_ts_res[, c("PC", "extinction_driver")], ordered_wass_res[, c("PC", "extinction_driver")])
full_res <- ordered_ts_res |> 
  full_join(
    ordered_wass_res, by = c("PC", "extinction_driver")
  )
plot(p_value.y ~ p_value.x, data = full_res)
mod <- lm(p_value.y ~ p_value.x, data = full_res)
mod
summary(mod)
# essentially the same p-value order - I'm guessing the tiny amount of variation is because of
# stochasticity in the permutation testing
# check that all significant results are significant under both methods, and same for non-significant
# results
identical(full_res[full_res$p_value.x < 0.05, ], full_res[full_res$p_value.y < 0.05, ])
identical(full_res[full_res$p_value.x > 0.05, ], full_res[full_res$p_value.y > 0.05, ])

# same for adjusted p-values
full_res <- full_res |> 
  mutate(
    p_value.x_adj = p.adjust(p_value.x),
    p_value.y_adj = p.adjust(p_value.y)
  )
identical(full_res[full_res$p_value.x_adj < 0.05, ], full_res[full_res$p_value.y_adj < 0.05, ])
identical(full_res[full_res$p_value.x_adj > 0.05, ], full_res[full_res$p_value.y_adj > 0.05, ])


# inspect results
pc_ex_drive_res |> 
  # filter(
  #   abs(ses) > 1.96
  # ) |> 
  arrange(
    p_adjusted
  )


# Keep only the axes for which p-value (from twosamples results) < 0.05
# Use this because SES is unreliable for skewed distributions of null statistics
sig_pc_exdrive <- res_ts |> 
  filter(
    p_value < 0.05
  )
sig_combos <- sig_pc_exdrive |> 
  select(
    PC, extinction_driver
  ) |> 
  mutate(
    combo = paste(PC, extinction_driver, sep = "_")
  )

# Test for SPECIFIC distributional differences in each of these axis/threat combinations
# --> Mean shift (via lm/ANOVA)
# --> Shift in variance (via Levene's test)

# Mean shift and variance inequality
# I want to run a t-test  on each significant PC/extinction driver combination, comparing the 
# distribution of the significant combination with that of the distribution of the same PC values of
# other threatened species

ttest_res <- pbapply::pblapply(
  1:nrow(sig_combos),
  function(combo_number){
    
    # define focal combination
    focal_combo <- sig_combos[combo_number, ]
    
    focal_distrib <- threat_colour_long |> 
      filter(
        iucn_cat != "LC",
      ) |> 
      filter(
        PC == focal_combo$PC,
        ex_driver == focal_combo$extinction_driver
      )
    non_focal_distrib <- threat_colour_long |> 
      filter(
        iucn_cat != "LC"
      ) |> 
      filter(
        PC == focal_combo$PC,
        ex_driver != focal_combo$extinction_driver
      ) |> 
      mutate(
        ex_driver = "non_focal"
      )
    mod_dat <- focal_distrib |> 
      bind_rows(
        non_focal_distrib
      ) |> 
      select(
        PC, ex_driver, PC_value
      )
    
    # calculate observed t statistic (for mean shift) and Levene's test (for variance inequality)
    obs_t_stat <- t.test(PC_value ~ ex_driver, data = mod_dat)$statistic
    obs_f_val <- suppressWarnings(car::leveneTest(PC_value ~ ex_driver, data = mod_dat)["group", "F value"])
    
    null_distribs <- lapply(1:1000, 
                            function(i){
                              sample_rows <- sample(
                                1:nrow(mod_dat),
                                size = length(focal_distrib$PC_value),
                                replace = FALSE
                              )
                              sample_distrib <- mod_dat[, c("ex_driver", "PC_value")]
                              sample_distrib[sample_rows, "ex_driver"] <- focal_combo$extinction_driver
                              sample_distrib[-sample_rows, "ex_driver"] <- "non_focal"
                              return(sample_distrib)
                            }
    )
    
    # Get null distribution of t statistics
    null_t_stats <- unlist(
      lapply(
        null_distribs,
        function(sample){
          null_t_stat <- t.test(PC_value ~ ex_driver, data = sample)$statistic
          return(null_t_stat)
        }
      )
    )
    # Calculate SES for t statistics
    null_t_mean <- mean(null_t_stats)
    null_t_sd <- sd(null_t_stats)
    mean_shift_es <- obs_t_stat - null_t_mean
    mean_shift_ses <- mean_shift_es / null_t_sd
    
    # Get null distribution of F values
    null_f_vals <- unlist(
      lapply(
        null_distribs,
        function(sample){
          null_f_val <- suppressWarnings(car::leveneTest(PC_value ~ ex_driver, data = sample)["group", "F value"])
          return(null_f_val)
        }
      )
    )
    # Calculate SES for F values
    null_f_mean <- mean(null_f_vals)
    null_f_sd <- sd(null_f_vals)
    var_inequal_es <- obs_f_val - null_f_mean
    var_inequal_ses <- var_inequal_es / null_f_sd
    
    res <- data.frame(
      PC = focal_combo$PC, 
      ex_driver = focal_combo$extinction_driver, 
      mean_shift_obs = obs_t_stat, 
      mean_shift_null_mean = null_t_mean, 
      mean_shift_null_sd = null_t_sd, 
      mean_shift_es = mean_shift_es,
      mean_shift_ses = mean_shift_ses, 
      var_inequal_obs = obs_f_val, 
      var_inequal_null_mean = null_f_mean, 
      var_inequal_sd = null_f_sd, 
      var_inequal_es = var_inequal_es,
      var_inequal_ses = var_inequal_ses)
    
    return(res)
    
  }
)
results <- do.call(rbind, ttest_res)

# plot significant meanshift drivers
results |> 
  filter(
    abs(mean_shift_ses) > 1.96
  ) |> 
  ggplot(aes(x = ex_driver, y = mean_shift_es, fill = ex_driver)) + 
  geom_col() + 
  facet_wrap(~ PC)

# and significant variance inequality drivers
results |> 
  filter(
    abs(var_inequal_ses) > 1.96
  ) |> 
  ggplot(aes(x = ex_driver, y = var_inequal_es, fill = ex_driver)) + 
  geom_col() + 
  facet_wrap(~ PC)

# check if inequality of variance SES is associated with mean shift SES
mod <- lm(abs(mean_shift_ses) ~ abs(var_inequal_ses), data = results)
summary(mod)
# no, not associated

# Same, but using Levene's test for equality of variances between two distributions - this will
# tell me if one distribution is more clustered around the mean than another


focal_distrib <- threat_colour_long |> 
  filter(
    iucn_cat != "LC",
  ) |> 
  filter(
    PC == focal_combo$PC,
    ex_driver == focal_combo$extinction_driver
  )
non_focal_distrib <- threat_colour_long |> 
  filter(
    iucn_cat != "LC"
  ) |> 
  filter(
    PC == focal_combo$PC,
    ex_driver != focal_combo$extinction_driver
  ) |> 
  mutate(
    ex_driver = "non_focal"
  )
mod_dat <- focal_distrib |> 
  bind_rows(
    non_focal_distrib
  ) |> 
  select(
    PC, ex_driver, PC_value
  )

obs_t_stat <- t.test(PC_value ~ ex_driver, data = mod_dat)$statistic
null_distribs <- lapply(1:1000, 
                        function(i){
                          sample_rows <- sample(
                            1:nrow(mod_dat),
                            size = length(focal_distrib$PC_value),
                            replace = FALSE
                          )
                          sample_distrib <- mod_dat[, c("ex_driver", "PC_value")]
                          sample_distrib[sample_rows, "ex_driver"] <- focal_combo$extinction_driver
                          sample_distrib[-sample_rows, "ex_driver"] <- "non_focal"
                          return(sample_distrib)
                        }
)
null_t_stats <- unlist(
  lapply(
    null_distribs,
    function(sample){
      null_t_stat <- t.test(PC_value ~ ex_driver, data = sample)$statistic
      return(null_t_stat)
    }
  )
)
null_mean <- mean(null_t_stats)
null_sd <- sd(null_t_stats)
ses <- (obs_t_stat - null_mean) / null_sd

t_mod <- t.test(PC_value ~ ex_driver, data = mod_dat)
t_mod
summary(t_mod)

threat_col_sig_combos <- threat_colour_long |> 
  filter(
    iucn_cat != "LC",
  ) |> 
  mutate(
    combo = paste(PC, ex_driver, sep = "_")
  ) |> 
  filter(
    combo %in% sig_combos$combo
  )

library(lme4)
mean_shift_mod <- lme4::glmer(PC_value ~ 1 + ex_driver + (1 | PC), data = threat_col_sig_combos)

mean_shift_mod
summary(mean_shift_mod)

threat_col_sig_combos |> 
  ggplot(aes(x = ex_driver, y = PC_value, fill = ex_driver)) + 
  geom_boxplot() + 
  stat_smooth(method = "lm", fullrange = T) + 
  facet_wrap(~ PC)

plot(mean_shift_mod)

mean_shift_mods <- pbapply::pblapply(
  sig_combos$combo,
  function(var_combo){
    
    combos <- sig_combos |> 
      combo = var_combo
    
    focal_distrib <- threat_col_sig_combos |> 
      filter(
        combo == var_combo,
      ) |> 
      pull(
        PC_value
      )
    non_focal_distrib <- threat_col_sig_combos |> 
      filter(
        pc == unique(combos$PC),
        ex_driver 
      )
    
    mod <- lm()
    
  }
)





### ALTERNATIVE HIERARCHICAL MODELLING APPROACH ----

# use the clean (non-long data)

# Add an extra column called 'focal' threat - this needs to be specified each time the model is run
# Here we'll use hunting and collection as the focal variable
# Remove duplicate rows and also filter out non-threatened (LC) species
focal_threat <- "hunt_col"
focal_pc <- "PC3"
mod_data <- threat_colour_clean |> 
  mutate(
    focal_threat = as.factor(ifelse(ex_driver == focal_threat, 1, 0))
  ) |> 
  distinct(
    jetz_species, focal_threat,
    .keep_all = TRUE
  ) |> 
  filter(
    iucn_cat != "LC"
  )

# set up model
form <- formula(get(focal_pc) ~ focal_threat)
mod <- lm(form, data = mod_data)
summary(mod)
# I think this is essentailly just testing for a difference in means between the two groups - 
# not holistic difference in distribution
# CONCLUSION: continue with my previous methodology (random resampling combined with CVM test)