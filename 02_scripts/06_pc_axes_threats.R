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
# Colour space to use
space <- "lab"
# Exclude past threats?
exclude_past_threats <- FALSE
# Remove extinct (EX) species?
remove_extinct <- TRUE

# Load data ----

# Load threat data with centroid distances (binary variables version)
if(remove_extinct){
  ext_par <- "remove_extinct_"
} else {
  ext_par <- "keep_extinct_"
}
if(exclude_past_threats) {
  past_threat_par <- "exclude_past_threats_"
} else {
  past_threat_par <- NULL
}
if(latest == TRUE){
  filename <- paste0(ext_par, past_threat_par, "centroid_jetz_threat_matrix_latest_2026-01-07.csv")
} else if(latest == FALSE){
  filename <- paste0(ext_par, past_threat_par, "centroid_jetz_threat_matrix_", cutoff_year, "_cutoff_year.csv")
}
threat_data <- read.csv(
  here::here(
    "03_output_data", "03a_threat_centroid",
    filename
  )
)

# load colour pattern space UMAP (created in Chapter 1 - patch-pipeline)
umap_path <- paste0("G:/My Drive/patch-pipeline/2_Patches/3_OutputData/", clade, "/2_PCA_ColourPattern_spaces/2_UMAP/", clade, ".matchedsex.patches.nn.25.mindist.0.1.", space, ".UMAP.rds")
umap <- readRDS(umap_path)[["layout"]]

# Data preparation ----


# Exclude Data Deficient (DD), Extinct (EX) and Extinct in the Wild (EW) species
# threat_matrix <- threat_matrix %>% 
#   filter(
#     !(iucn_cat %in% c("DD", "EX", "EW"))
#   )


# can make a long version of threat matrix (pivot on threats)
threat_bin_long <- threat_data %>% 
  tidyr::pivot_longer(
    c(acc_mort, clim_chan, hab_loss, hunt_col, invas_spec, other, pollut),
    names_to = "threat_type",
    values_to = "threat_present"
  )



# Pivot longer on PC axes (so we can facet plots by axis)
threat_pc_long <- threat_data %>% 
  tidyr::pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "PC_value"
  ) 

# get wide threat data with UMAP axes
umap_sppsex <- data.frame(species = sapply(strsplit(rownames(umap), split = "-"), "[", 1),
                          sex = sapply(strsplit(rownames(umap), split = "-"), "[", 2), 
                          umap)


threat_umap <- threat_data %>% 
  inner_join(umap_sppsex, by = join_by("jetz_species" == "species", "sex")) %>% 
  rename(
    UMAP1 = X1,
    UMAP2 = X2
  )
# this throws a warning but it's just because we have male and female data together
# - it's not a problem

# Plotting ----


# Proportional 2D density plots (The Juice) ----

# set focal extinction driver
focal_threat <- "hunt_col"
# set axes (PCs or UMAP axes)
ax_1 <- "PC1"
ax_2 <- "PC2"

# Get proportional density of species threatened by focal threat vs ALL other species
prop_dens <- prop_dens_2d(
  threat_umap,
  focal_threat = focal_threat,
  x_axis = ax_1, y_axis = ax_2,
  threatened_spp_only = FALSE,
  vs = "other_species",
  n_bins = 200,
  return_params = TRUE
)


# plot as heatmap
plot_prop_dens_2d(
  prop_dens
)

# plot non-proportional (raw) density of focal threat species, for comparison
pollut_dens <- MASS::kde2d(
  x = threat_umap[which(threat_umap[[focal_threat]] == 1), ax_1], 
  y = threat_umap[which(threat_umap[[focal_threat]] == 1), ax_2],
  n = 200,
  lims = c(range(threat_umap[[ax_1]]), range(threat_umap[[ax_2]]))
)
plot_prop_dens_2d(pollut_dens)
# plot density of all species, for comparison
all_dens <- MASS::kde2d(
  x = threat_umap[, ax_1], 
  y = threat_umap[, ax_2],
  n = 200,
  lims = c(range(threat_umap[[ax_1]]), range(threat_umap[[ax_2]]))
)
plot_prop_dens_2d(all_dens)

# Plot multiple threats together
# Set vs parameter
# Either 'all_species' to plot density of species threatened by focal threat proportional to
# all other species (to see how the denisty of the threat differs from the original density of the 
# space) or 'other_species' to plot density of species threatened by focal threat proportional to
# species NOT threatened by the focal threat (i.e. excluding the focal species, to see how the 
# density of the threat differs from the density of the rest of the species)
# Actually the difference between these methods is very minor, at least when not restricting
# to threatened species only
vs_par <- "other_species"
ex_drivers <- c("clim_chan", "hab_loss", "pollut", "hunt_col", "acc_mort", "invas_spec")
names(ex_drivers) <- ex_drivers
# First get the proportional densities
prop_densities <- pbapply::pblapply(
  ex_drivers,
  function(driver){
    
    prop_dens <- prop_dens_2d(
      threat_umap,
      focal_threat = driver,
      x_axis = ax_1, y_axis = ax_2,
      threatened_spp_only = FALSE,
      vs = vs_par,
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
    paste0("extinction_drivers_vs_", vs_par, "_", space, "_", ax_1, ax_2, ".png")
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
  max_density_val = global_max_dens,
  asp = TRUE
)
dev.off()


# UMAP density plot by threat (non-proportional)
threat_umap %>% 
  bind_rows(
    threat_umap %>% 
      mutate(
        all_spec = 1
      )
  ) %>% 
  mutate(
    all_spec = ifelse(is.na(all_spec), 0, all_spec)
  ) %>% 
  tidyr::pivot_longer(
    c(acc_mort, clim_chan, hab_loss, hunt_col, invas_spec, other, pollut, all_spec),
    names_to = "threat_type",
    values_to = "threat_present"
  ) %>% 
  filter(
    iucn_cat != "LC",
    !(threat_type %in% c(NA, "other")),
    threat_present == 1
  )  %>%  
  mutate(
    threat_type = factor(threat_type, levels = c("all_spec", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
    sex = factor(sex, levels = c("M", "F"))
  ) %>% 
  # filter(
  #   ex_driver != "all"
  # ) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2)) + 
  geom_density_2d_filled() + 
  facet_wrap(~ threat_type)



# Same for PC1/2
threat_data %>% 
  bind_rows(
    threat_data %>% 
      mutate(
        all_spec = 1
      )
  ) %>% 
  mutate(
    all_spec = ifelse(is.na(all_spec), 0, all_spec)
  ) %>% 
  tidyr::pivot_longer(
    c(acc_mort, clim_chan, hab_loss, hunt_col, invas_spec, other, pollut, all_spec),
    names_to = "threat_type",
    values_to = "threat_present"
  ) %>% 
  filter(
    iucn_cat != "LC",
    !(threat_type %in% "other"),
    threat_present == 1
  )  %>%  
  mutate(
    threat_type = factor(threat_type, levels = c("all_spec", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
    sex = factor(sex, levels = c("M", "F"))
  ) %>% 
  # filter(
  #   ex_driver != "all"
  # ) %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_density_2d_filled() + 
  facet_wrap(~ threat_type)


# boxplot of PC values by threat type - exclude LC species
# Compare to all threatened species
threat_pc_long %>% 
  bind_rows(
    threat_pc_long %>% 
      mutate(
        all_spec = 1
      )
  ) %>% 
  mutate(
    all_spec = ifelse(is.na(all_spec), 0, all_spec)
  ) %>% 
  tidyr::pivot_longer(
    c(acc_mort, clim_chan, hab_loss, hunt_col, invas_spec, other, pollut, all_spec),
    names_to = "threat_type",
    values_to = "threat_present"
  ) %>% 
  filter(
    iucn_cat != "LC",
    !(threat_type %in% c(NA, "other")),
    PC %in% paste0("PC", 4),
    threat_present == 1
  ) %>% 
  mutate(
    threat_type = factor(threat_type, levels = c("all_spec", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
    sex = factor(sex, levels = c("M", "F"))
  ) %>% 
  ggplot(aes(x = threat_type, y = PC_value, fill = threat_type)) + 
  geom_boxplot(outliers = F) + 
  facet_grid(~ sex)

# Density plot of all species (excluding LC) vs species threatened by each threat type,
# faceted by PC axis (first 7 PCs only)
threat_pc_long %>% 
  bind_rows(
    threat_pc_long %>% 
      mutate(
        all_spec = 1
      )
  ) %>% 
  mutate(
    all_spec = ifelse(is.na(all_spec), 0, all_spec)
  ) %>% 
  tidyr::pivot_longer(
    c(acc_mort, clim_chan, hab_loss, hunt_col, invas_spec, other, pollut, all_spec),
    names_to = "threat_type",
    values_to = "threat_present"
  ) %>% 
  filter(
    iucn_cat != "LC",
    !(threat_type %in% c(NA, "other")),
    PC %in% paste0("PC", 5:8),
    threat_present == 1
  ) %>% 
  mutate(
    threat_type = factor(threat_type, levels = c("all_spec", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
    sex = factor(sex, levels = c("M", "F"))
  ) %>% 
  filter(
    threat_type %in% c("all_spec", "hunt_col")
  ) %>% 
  ggplot(aes(x = PC_value, fill = threat_type)) + 
  geom_density(alpha = 0.4) + 
  geom_segment(x = 0, xend = 0, y = 0, yend = 0.04, lwd = 0.3, colour = "grey30") + # Vertical line at 0
  facet_grid(rows = vars(PC), cols = vars(sex)) + 
  theme_minimal()


# function to do the above
plot_pc_distrib <- function(pc_axis, threat_type, long_data, absolute = FALSE) {
  
  if(absolute == TRUE){
    long_data$PC_value <- abs(long_data$PC_value)
  }
  
  p <- long_data %>% 
      bind_rows(
        long_data %>% 
          mutate(
            all_spec = 1
          )
      ) %>% 
        mutate(
          all_spec = ifelse(is.na(all_spec), 0, all_spec)
        ) %>% 
        tidyr::pivot_longer(
          c(acc_mort, clim_chan, hab_loss, hunt_col, invas_spec, other, pollut, all_spec),
          names_to = "ex_driver",
          values_to = "threat_present"
        ) %>% 
    filter(
      iucn_cat != "LC",
      !(ex_driver %in% c(NA, "other")),
      PC == pc_axis,
      threat_present == 1
    ) %>% 
    mutate(
      ex_driver = factor(ex_driver, levels = c("all_spec", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
      sex = factor(sex, levels = c("M", "F"))
    ) %>% 
    filter(
      ex_driver %in% c("all_spec", threat_type)
    ) %>% 
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
pc_axis <- "PC1"
pc_plots <- lapply(
  threat_types,
  plot_pc_distrib,
  pc_axis = pc_axis,
  long_data = threat_pc_long,
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
  long_data = threat_pc_long,
  absolute = TRUE
)
abs_allthreats <- ggpubr::ggarrange(plotlist = abs_pc_plots, common.legend = TRUE, ncol = 2, nrow = 3)
ggpubr::annotate_figure(
  abs_allthreats,
  bottom = ggpubr::text_grob(pc_axis),
  left = ggpubr::text_grob("Density", rot = 90)
)


# Statistics ----

# Perform a Cramer-von Mises test to check for difference in distribution of PC values between 
# species threatened by focal extinction driver and species not threatened by focal driver
# excludes threatened species for which there is
# no threat data
# This is the correct thing to do because we cannot say for sure that these threatened species
# are not threatened by the focal threat
# It includes species which are threatened by non-significant drivers of extinction, as 
# determined by Stewart et al 2025 Nat Ecol Evol (see Supplementary_Dataset.xlsx for details)


# apply across each extinction driver and PC
# could parallelise this but it only takes ~30s for 7 axes
# define threat types and focal PCs first 
threat_types <- c(hab_loss = "hab_loss", hunt_col = "hunt_col", invas_spec = "invas_spec", clim_chan = "clim_chan", pollut = "pollut", acc_mort = "acc_mort")
pcs <- paste0("PC", 1:7)

# # Run across PCs/extinction drivers
# pc_ex_drive_res <- pbapply::pblapply(
#   
#   threat_types,
#   
#   function(threat_type){
#     
#     cvm <- pbapply::pblapply(
#       pcs, 
#       test_pc_threat,
#       threat_pc_long = threat_pc_long,
#       threat = threat_type,
#       n = 1000,
#       seed = 42,
#       return_distrib = FALSE,
#       stat_test = "wass"
#     )
#     names(cvm) <- pcs
#     cvm <- data.table::rbindlist(cvm, idcol = "PC")
#     # Bonferroni correction for p-values
#     # cvm$p_adjusted <- p.adjust(cvm$p_value, method = "bonferroni")
#     return(cvm)
#     
#   }
#   
# )
# pc_ex_drive_res <- do.call(rbind, pc_ex_drive_res)
# # Bonferroni correction for p-values
# pc_ex_drive_res$p_adjusted <- p.adjust(pc_ex_drive_res$p_value, method = "bonferroni")

# compare to results just from using twosamples::wass_test (which uses essentially the
# same permutation test to determine a count-based p-value)
# ## THIS IS WHAT I SHOULD ACTUALLY USE - no benefit of using my own code to do the same thing
# and twosamples does the same thing much faster
sexes <- c("both_sexes")
threatened_spp_only_par <- FALSE
res_twosamples <- lapply(
  sexes,
  function(sex){

    sex_res <- lapply(
      
      threat_types,
      
      function(threat_type){

        pc_stats <- lapply(
          pcs,
          function(pc){
            # filter to focal PC here to save overhead
            dat <- threat_pc_long[threat_pc_long$PC == pc, ]
            stat <- test_pc_threat_twosamples(pc, 
                                              threat_colour_long = dat, 
                                              threat_type = threat_type,
                                              focal_sex = sex,
                                              n_boots = 1000,
                                              stat_test = "wass",
                                              threatened_spp_only = threatened_spp_only_par)
            return(stat)
          }
        )
        
        threat_res <- do.call(rbind, pc_stats)
      }
      
    )
    sex_res <- data.table::rbindlist(sex_res, idcol = "extinction_driver")
  }
)
# format results
res_twosamples <- do.call(rbind, res_twosamples)
colnames(res_twosamples) <- c("extinction_driver", "test_stat", "p_value", "PC", "sex")
res_twosamples <- res_twosamples %>% 
  select(PC, extinction_driver, sex, test_stat, p_value)

# Bonferroni correction for p-values (within sex)
res_twosamples <- res_twosamples %>% 
  group_by(
    sex
  ) %>% 
  mutate(
    p_adjusted = p.adjust(p_value, method = "bonferroni")
  )
# NOTE that this is not necessary - I only want to identify which threat types might differentially
# threaten which PC axes - better to keep the net wide here as the mean shift/variance tests
# will add an additional filter

# Save results as CSV
# set filename based on parameters
if(threatened_spp_only_par == TRUE){
  wass_filename <- "PC_ex_driver_wasserstein_distance_threatened_spp_only.csv"
} else if(threatened_spp_only_par == FALSE){
  wass_filename <- "PC_ex_driver_wasserstein_distance_all_spp.csv"
}

write.csv(
  res_twosamples,
  file = here::here(
    "03_output_data", "06_pc_axes_threats",
    wass_filename
  )
)

# Comparison to my own function for calculating SES/ES/p-value
# NOTE that I am no longer using my own function - it does the same thing as the twosamples function
# and is much slower
# This code snippet is no longer necessary
# ordered_ts_res <- res_twosamples %>% 
#   arrange(p_value)
# ordered_wass_res <- wass_res_new %>% 
#   arrange(p_value)
# identical(ordered_ts_res[, c("PC", "extinction_driver")], ordered_wass_res[, c("PC", "extinction_driver")])
# full_res <- ordered_ts_res %>% 
#   full_join(
#     ordered_wass_res, by = c("PC", "extinction_driver")
#   )
# plot(p_value.y ~ p_value.x, data = full_res)
# mod <- lm(p_value.y ~ p_value.x, data = full_res)
# mod
# summary(mod)
# # essentially the same p-value order - I'm guessing the tiny amount of variation is because of
# # stochasticity in the permutation testing
# # check that all significant results are significant under both methods, and same for non-significant
# # results
# identical(full_res[full_res$p_value.x < 0.05, ], full_res[full_res$p_value.y < 0.05, ])
# identical(full_res[full_res$p_value.x > 0.05, ], full_res[full_res$p_value.y > 0.05, ])
# 
# # same for adjusted p-values
# full_res <- full_res %>% 
#   mutate(
#     p_value.x_adj = p.adjust(p_value.x),
#     p_value.y_adj = p.adjust(p_value.y)
#   )
# identical(full_res[full_res$p_value.x_adj < 0.05, ], full_res[full_res$p_value.y_adj < 0.05, ])
# identical(full_res[full_res$p_value.x_adj > 0.05, ], full_res[full_res$p_value.y_adj > 0.05, ])

### NOTE 2026-03-25
### HAVE INCORPORATED SEX-SPECIFIC ANALYSIS UP TO HERE, BUT NEED TO INCORPORATE IT INTO THE BELOW


# Keep only the axes for which p-value (from twosamples results) < 0.05
# Use this because SES is unreliable for skewed distributions of null statistics
sig_pc_exdrive <- res_twosamples %>% 
  filter(
    p_value < 0.05
  )
sig_combos <- sig_pc_exdrive %>% 
  select(
    PC, extinction_driver, sex
  ) %>% 
  mutate(
    combo = paste(PC, extinction_driver, sex, sep = "_")
  )

# Test for SPECIFIC distributional differences in each of these axis/threat combinations
# --> Mean shift (via difference in means, with permutation test)
# --> Shift in variance (via difference in standard deviation, via permutation test)

# 2026-10-07 New method for detecting mean shifts and variance differences
# test function
# first filter data to only the focal PC (essential to save overhead, otherwise function
# times out if you call it in an lapply because the dataset is so big)
focal_combo_num <- 9
focal_pc <- sig_combos[[focal_combo_num, "PC"]]
test_data <- threat_pc_long[threat_pc_long[["PC"]] == focal_pc, ]
debug(test_meanshift_varshift)
test_res <- test_meanshift_varshift(
  combo_number = focal_combo_num,
  sig_combos_df = sig_combos,
  filtered_colour_data = test_data,
  n_boots = 1000,
  threatened_spp_only = FALSE,
  focal_sex = "both_sexes"
)

# now run across all combos
spec_test_res <- lapply(
  1:nrow(sig_combos),
  function(focal_combo_num){
    # get focal PC
    focal_pc <- sig_combos[[focal_combo_num, "PC"]]
    # get focal threat
    focal_threat <- sig_combos[[focal_combo_num, "extinction_driver"]]
    # print which combo working on
    message("Focal combo ", focal_combo_num, " of ", nrow(sig_combos), " (", focal_pc, ", ", focal_threat, ")")
    # filter data to focal PC only
    focal_dat <- threat_pc_long[threat_pc_long[["PC"]] == focal_pc, ]
    # run testing function for mean shift and log variance ratio
    focal_res <- test_meanshift_varshift(
      combo_number = focal_combo_num,
      sig_combos_df = sig_combos,
      filtered_colour_data = focal_dat,
      n_boots = 1000,
      threatened_spp_only = FALSE,
      focal_sex = "both_sexes"
    )
  }
)
spec_test_res <- do.call(rbind, spec_test_res)

## DEPRECATED - not run
# Mean shift and variance inequality
# I want to run a t-test and Levene's test on each significant PC/extinction driver combination,
# comparing the distribution of the significant combination with that of the distribution of the 
# same PC values of other threatened species
# Note that this code will NOT work if you want to test sexes separately - only both together

# spec_test_res <- pbapply::pblapply(
#   1:nrow(sig_combos),
#   ttest_levenetest_pc_threat,
#   sig_combos_df = sig_combos, # this gives us the PC axis, threat type and focal sex
#   threat_pc_long = threat_pc_long,
#   n_boots = 1000,
#   threatened_spp_only = threatened_spp_only_par
# )
# spec_test_res <- do.call(rbind, spec_test_res)

# Save results as CSV
# set filename based on parameters
if(threatened_spp_only_par == TRUE){
  spec_filename <- "PC_ex_driver_meanshift-logvar_threatened_spp_only.csv"
} else if(threatened_spp_only_par == FALSE){
  spec_filename <- "PC_ex_driver_meanshift-logvar_all_spp.csv"
}

write.csv(
  spec_test_res,
  file = here::here(
    "03_output_data", "06_pc_axes_threats",
    spec_filename
  ), 
  row.names = FALSE
)

# reload results for plotting
threatened_spp_only_par <- FALSE
# set filename based on parameters
if(threatened_spp_only_par == TRUE){
  spec_filename <- "PC_ex_driver_ttestlevtest_threatened_spp_only.csv"
} else if(threatened_spp_only_par == FALSE){
  spec_filename <- "PC_ex_driver_ttestlevtest_all_spp.csv"
}
spec_test_res <- read.csv(
  file = here::here(
    "03_output_data", "06_pc_axes_threats",
    spec_filename
  )
)

# plot significant meanshift drivers - vertical bars
mean_shift_plot <- spec_test_res %>% 
  filter(
    mean_shift_p < 0.05
  ) %>% 
  ggplot(aes(x = ex_driver, y = mean_shift_es, fill = ex_driver)) + 
  geom_col() + 
  facet_wrap(~ PC) + 
  labs(x = element_blank(), y = "Effect size (mean shift)", fill = "Extinction driver") +
  theme_bw() + 
  theme( axis.text.x = element_blank(), # Remove x axis tick labels
         axis.text.y = element_blank(), # Remove y axis tick labels
         axis.ticks = element_blank(),  # Remove ticks
         legend.position = "inside", 
         legend.title.position = "top",
         legend.position.inside = c(0.8, 0.15), 
         legend.direction = "horizontal", 
         legend.text.position = "bottom") 
mean_shift_plot

# plot significant meanshift drivers - horizontal bars
mean_shift_plot <- spec_test_res %>% 
  # filter(
  #   mean_shift_p < 0.05
  # ) %>% 
  mutate(
    signif = ifelse(mean_shift_p < 0.05, "y", "n")
  ) %>% 
  filter(PC %in% paste0("PC", 1:6)) %>% 
  ggplot(aes(y = ex_driver, x = mean_shift_ses, fill = signif)) + 
  geom_col() + 
  geom_vline(xintercept = 0) + 
  facet_wrap(~ PC, dir = "v") + 
  labs(y = "Extinction driver", x = paste0("Effect size (mean shift)")) + 
  scale_y_discrete(
    labels = c(
      "pollut" = "Pollution",
      "invas_spec" = "Invasive species",
      "hunt_col" = "Hunting & Collection",
      "hab_loss" =  "Habitat Loss",
      "clim_chan" =  "Climate Change",
      "acc_mort" = "Accidental Mortality"
    )
  ) + 
  scale_fill_discrete(palette = c("grey70", "grey25")) + 
  theme_bw() + 
  theme(legend.position = "none")
mean_shift_plot

# and significant variance inequality drivers - vertical bars
log_var_ratio_plot <- spec_test_res %>% 
  filter(
    log_var_ratio_p < 0.05
  ) %>% 
  filter(PC %in% paste0("PC", 1:6)) %>% 
  ggplot(aes(x = ex_driver, y = log_var_ratio_es, fill = ex_driver)) + 
  geom_col() + 
  facet_wrap(~ PC) + 
  labs(x = element_blank(), y = "Effect size (variance inequality)", fill = "Extinction driver") +
  theme_bw() + 
  theme( axis.text.x = element_blank(), # Remove x axis tick labels
         axis.text.y = element_blank(), # Remove y axis tick labels
         axis.ticks = element_blank(),  # Remove ticks
         legend.position = "inside", 
         legend.title.position = "top",
         legend.position.inside = c(0.8, 0.15), 
         legend.direction = "horizontal", 
         legend.text.position = "bottom") 
log_var_ratio_plot

# horizontal bars
log_var_ratio_plot <- spec_test_res %>% 
  # filter(
  #   log_var_ratio_p < 0.05
  # ) %>% 
  mutate(
    signif = ifelse(log_var_ratio_p < 0.05, "y", "n"),
    sd_diff = ifelse(delta_sd > 0, "pos", "neg")
  ) %>% 
  filter(PC %in% paste0("PC", 1:6)) %>% 
  ggplot(aes(y = ex_driver, x = log_var_ratio_ses, fill = signif)) + 
  geom_col(linewidth = 1) + 
  geom_vline(xintercept = 0) + 
  facet_wrap(~ PC) + 
  labs(y = "Extinction driver", x = "Effect size (variance inequality)") + 
  scale_y_discrete(
    labels = c(
      "pollut" = "Pollution",
      "invas_spec" = "Invasive species",
      "hunt_col" = "Hunting & Collection",
      "hab_loss" =  "Habitat Loss",
      "clim_chan" =  "Climate Change",
      "acc_mort" = "Accidental Mortality"
    )
  ) +
  scale_fill_discrete(palette = c(adjustcolor("lightgrey"), "grey25")) + 
#  scale_fill_discrete(palette = c(adjustcolor("darkblue", alpha = 0.8), adjustcolor("darkred", alpha = 0.8))) +
  theme_bw() + 
  theme(legend.position = "none")
log_var_ratio_plot

# save plots
if(threatened_spp_only_par == TRUE){
  ms_plot_filename <- "PC_ex_driver_ttestmeanshift_threatened_spp_only.svg"
  vi_plot_filename <- "PC_ex_driver_levtestvarinequality_threatened_spp_only.svg"
} else if(threatened_spp_only_par == FALSE){
  ms_plot_filename <- "PC_ex_driver_ttestmeanshift_all_spp.svg"
  vi_plot_filename <- "PC_ex_driver_levtestvarinequality_all_spp.svg"
}
ggsave(
  filename = ms_plot_filename,
  plot = mean_shift_plot, 
  device = "svg", 
  path = here::here(
    "04_output_plots", "06_pc_axes_threats", "02_distributional_differences"
  ), 
  width = 180, height = 140, units = "mm"
)
ggsave(
  filename = vi_plot_filename,
  plot = var_inequal_plot, 
  device = "svg", 
  path = here::here(
    "04_output_plots", "06_pc_axes_threats", "02_distributional_differences"
  ), 
  width = 180, height = 140, units = "mm"
)


# compute correlation matrix to see if suites of threats act together
# in terms of mean shifts
library(RColorBrewer)
col_pal <- colorRampPalette(brewer.pal(10, "PiYG"))(256)

data_cor <- tidyr::pivot_wider(spec_test_res[, c("PC", "ex_driver", "mean_shift_obs")], names_from = PC, values_from = mean_shift_obs)
data_cor_rownames <- data_cor$ex_driver
data_cor <- t(as.matrix(data_cor[, paste0("PC", 1:7)]))
data_cor[which(is.na(data_cor))] <- 0
colnames(data_cor) <- data_cor_rownames
cor_mat_meanshift <- cor(data_cor)
diag(cor_mat_meanshift) <- NA
heatmap(cor_mat_meanshift, symm = TRUE, col = col_pal)


# delve more deeply into threat covariance
# plot correlation between these two traits across PCs
data_cor <- tidyr::pivot_wider(spec_test_res[, c("PC", "ex_driver", "mean_shift_obs")], names_from = PC, values_from = mean_shift_obs)
data_cor_rownames <- data_cor$ex_driver
data_cor <- t(data_cor[, paste0("PC", 1:7)])
colnames(data_cor) <- data_cor_rownames
data_cor <- as.data.frame(data_cor)
data_cor$PC <- rownames(data_cor)

data_cor %>% 
  ggplot(aes(x = hab_loss, y = clim_chan)) + 
  geom_point() + 
  geom_smooth(method = "lm")
GGally::ggpairs(
  data_cor[, c("hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")],
  lower = list(continuous = GGally::wrap("smooth", method = "lm", colour = "blue"))
  )


# check if inequality of variance SES is associated with mean shift SES
mod <- lm(abs(mean_shift_ses) ~ abs(var_inequal_ses), data = spec_test_res)
summary(mod)
# no, not associated

# Same, but using Levene's test for equality of variances between two distributions - this will
# tell me if one distribution is more clustered around the mean than another


focal_distrib <- threat_pc_long %>% 
  filter(
    iucn_cat != "LC",
  ) %>% 
  filter(
    PC == focal_combo$PC,
    ex_driver == focal_combo$extinction_driver
  )
non_focal_distrib <- threat_pc_long %>% 
  filter(
    iucn_cat != "LC"
  ) %>% 
  filter(
    PC == focal_combo$PC,
    ex_driver != focal_combo$extinction_driver
  ) %>% 
  mutate(
    ex_driver = "non_focal"
  )
mod_dat <- focal_distrib %>% 
  bind_rows(
    non_focal_distrib
  ) %>% 
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

threat_col_sig_combos <- threat_pc_long %>% 
  filter(
    iucn_cat != "LC",
  ) %>% 
  mutate(
    combo = paste(PC, ex_driver, sep = "_")
  ) %>% 
  filter(
    combo %in% sig_combos$combo
  )

library(lme4)
mean_shift_mod <- lme4::glmer(PC_value ~ 1 + ex_driver + (1 | PC), data = threat_col_sig_combos)

mean_shift_mod
summary(mean_shift_mod)

threat_col_sig_combos %>% 
  ggplot(aes(x = ex_driver, y = PC_value, fill = ex_driver)) + 
  geom_boxplot() + 
  stat_smooth(method = "lm", fullrange = T) + 
  facet_wrap(~ PC)

plot(mean_shift_mod)

mean_shift_mods <- pbapply::pblapply(
  sig_combos$combo,
  function(var_combo){
    
    combos <- sig_combos %>% 
      combo = var_combo
    
    focal_distrib <- threat_col_sig_combos %>% 
      filter(
        combo == var_combo,
      ) %>% 
      pull(
        PC_value
      )
    non_focal_distrib <- threat_col_sig_combos %>% 
      filter(
        pc == unique(combos$PC),
        ex_driver 
      )
    
    mod <- lm()
    
  }
)




